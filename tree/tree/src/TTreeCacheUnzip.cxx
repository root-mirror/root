// Authors: Rene Brun   04/06/2006
//          Leandro Franco   10/04/2008
//          Fabrizio Furano (CERN) Aug 2009

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** 
\class TTreeCacheUnzip
\ingroup tree

A TTreeCache which exploits parallelized decompression of its own content.

*/

#include "TTreeCacheUnzip.h"
#include "TBranch.h"
#include "TChain.h"
#include "TEnv.h"
#include "TEventList.h"
#include "TFile.h"
#include "TMath.h"
#include "TROOT.h"
#include "TMutex.h"
#include "ROOT/RMakeUnique.hxx"

#ifdef R__USE_IMT
#include "ROOT/TThreadExecutor.hxx"
#include "ROOT/TTaskGroup.hxx"
#endif

extern "C" void R__unzip(Int_t *nin, UChar_t *bufin, Int_t *lout, char *bufout, Int_t *nout);
extern "C" int R__unzip_header(Int_t *nin, UChar_t *bufin, Int_t *lout);

TTreeCacheUnzip::EParUnzipMode TTreeCacheUnzip::fgParallel = TTreeCacheUnzip::kDisable;

// The unzip cache does not consume memory by itself, it just allocates in advance
// mem blocks which are then picked as they are by the baskets.
// Hence there is no good reason to limit it too much
Double_t TTreeCacheUnzip::fgRelBuffSize = .5;

ClassImp(TTreeCacheUnzip);

////////////////////////////////////////////////////////////////////////////////
/// Clear all baskets' state arrays.

void TTreeCacheUnzip::UnzipState::Clear(Int_t size) {
   for (Int_t i = 0; i < size; i++) {
      if (!fUnzipLen.empty()) fUnzipLen[i] = 0;
      if (fUnzipChunks) {
         if (fUnzipChunks[i]) fUnzipChunks[i].reset();
      }
      if (fUnzipStatus) fUnzipStatus[i].store(0);
   }
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TTreeCacheUnzip::UnzipState::IsUntouched(Int_t index) const {
   return fUnzipStatus[index].load() == kUntouched;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TTreeCacheUnzip::UnzipState::IsProgress(Int_t index) const {
   return fUnzipStatus[index].load() == kProgress;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TTreeCacheUnzip::UnzipState::IsFinished(Int_t index) const {
   return fUnzipStatus[index].load() == kFinished;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the basket is unzipped already. We must make sure the length in
/// fUnzipLen is larger than 0.

Bool_t TTreeCacheUnzip::UnzipState::IsUnzipped(Int_t index) const {
   return (fUnzipStatus[index].load() == kFinished) && (fUnzipChunks[index].get()) && (fUnzipLen[index] > 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset all baskets' state arrays. This function is only called by main
/// thread and parallel processing from upper layers should be disabled such
/// as IMT in TTree::GetEntry(). Other threads should not call this function 
/// since it is not thread-safe.

void TTreeCacheUnzip::UnzipState::Reset(Int_t oldSize, Int_t newSize) {
   std::vector<Int_t>       aUnzipLen    = std::vector<Int_t>(newSize, 0);
   std::unique_ptr<char[]> *aUnzipChunks = new std::unique_ptr<char[]>[newSize];
   std::atomic<Byte_t>     *aUnzipStatus = new std::atomic<Byte_t>[newSize];

   for (Int_t i = 0; i < newSize; ++i)
      aUnzipStatus[i].store(0);

   for (Int_t i = 0; i < oldSize; i++) {
      aUnzipLen[i]    = fUnzipLen[i];
      aUnzipChunks[i] = std::move(fUnzipChunks[i]);
      aUnzipStatus[i].store(fUnzipStatus[i].load());
   }

   if (fUnzipChunks) delete [] fUnzipChunks;
   if (fUnzipStatus) delete [] fUnzipStatus;

   fUnzipLen    = aUnzipLen;
   fUnzipChunks = aUnzipChunks;
   fUnzipStatus = aUnzipStatus;
}

////////////////////////////////////////////////////////////////////////////////
/// Set cache as finished. 
/// There are three scenarios that a basket is set as finished:
///    1. The basket has already been unzipped.
///    2. The thread is aborted from unzipping process.
///    3. To avoid other tasks/threads work on this basket,
///       main thread marks the basket as finished and designates itself
///       to unzip this basket.

void TTreeCacheUnzip::UnzipState::SetFinished(Int_t index) {
   fUnzipLen[index] = 0;
   fUnzipChunks[index].reset();
   fUnzipStatus[index].store((Byte_t)kFinished);
}

////////////////////////////////////////////////////////////////////////////////

void TTreeCacheUnzip::UnzipState::SetMissed(Int_t index) {
   fUnzipChunks[index].reset();
   fUnzipStatus[index].store((Byte_t)kFinished);
}

////////////////////////////////////////////////////////////////////////////////

void TTreeCacheUnzip::UnzipState::SetUnzipped(Int_t index, char* buf, Int_t len) {
   // Update status array at the very end because we need to be synchronous with the main thread.
   fUnzipLen[index] = len;
   fUnzipChunks[index].reset(buf);
   fUnzipStatus[index].store((Byte_t)kFinished);
}

////////////////////////////////////////////////////////////////////////////////
/// Start unzipping the basket if it is untouched yet.

Bool_t TTreeCacheUnzip::UnzipState::TryUnzipping(Int_t index) {
   Byte_t oldValue = kUntouched;
   Byte_t newValue = kProgress;
   return fUnzipStatus[index].compare_exchange_weak(oldValue, newValue, std::memory_order_release, std::memory_order_relaxed);
}

////////////////////////////////////////////////////////////////////////////////

TTreeCacheUnzip::TTreeCacheUnzip() : TTreeCache(),
   fAsyncReading(kFALSE),
   fEmpty(kTRUE),
   fCycle(0),
   fNseekMax(0),
   fUnzipGroupSize(0),
   fUnzipBufferSize(0),
   fNFound(0),
   fNMissed(0),
   fNStalls(0),
   fNUnzip(0)
{
   // Default Constructor.
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TTreeCacheUnzip::TTreeCacheUnzip(TTree *tree, Int_t buffersize) : TTreeCache(tree,buffersize),
   fAsyncReading(kFALSE),
   fEmpty(kTRUE),
   fCycle(0),
   fNseekMax(0),
   fUnzipGroupSize(0),
   fUnzipBufferSize(0),
   fNFound(0),
   fNMissed(0),
   fNStalls(0),
   fNUnzip(0)
{
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// Initialization procedure common to all the constructors.

void TTreeCacheUnzip::Init()
{
#ifdef R__USE_IMT
   fUnzipTaskGroup.reset();
#endif
   fIOMutex = std::make_unique<TMutex>(kTRUE);

   fCompBuffer = new char[16384];
   fCompBufferSize = 16384;

   fUnzipGroupSize = 102400; // Each task unzips at least 100 KB

   if (fgParallel == kDisable) {
      fParallel = kFALSE;
   }
   else if(fgParallel == kEnable || fgParallel == kForce) {
      fUnzipBufferSize = Long64_t(fgRelBuffSize * GetBufferSize());

      if(gDebug > 0)
         Info("TTreeCacheUnzip", "Enabling Parallel Unzipping");

      fParallel = kTRUE;

   }
   else {
      Warning("TTreeCacheUnzip", "Parallel Option unknown");
   }

   // Check if asynchronous reading is supported by this TFile specialization
   if (gEnv->GetValue("TFile.AsyncReading", 1)) {
      if (fFile && !(fFile->ReadBufferAsync(0, 0)))
         fAsyncReading = kTRUE;
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Destructor. (in general called by the TFile destructor)

TTreeCacheUnzip::~TTreeCacheUnzip()
{
   ResetCache();
   fUnzipState.Clear(fNseekMax);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a branch to the list of branches to be stored in the cache
/// this function is called by TBranch::GetBasket
/// Returns:
///  - 0 branch added or already included
///  - -1 on error

Int_t TTreeCacheUnzip::AddBranch(TBranch *b, Bool_t subbranches /*= kFALSE*/)
{
   return TTreeCache::AddBranch(b, subbranches);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a branch to the list of branches to be stored in the cache
/// this function is called by TBranch::GetBasket
/// Returns:
///  - 0 branch added or already included
///  - -1 on error

Int_t TTreeCacheUnzip::AddBranch(const char *branch, Bool_t subbranches /*= kFALSE*/)
{
   return TTreeCache::AddBranch(branch, subbranches);
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TTreeCacheUnzip::FillBuffer()
{

   if (fNbranches <= 0) return kFALSE;

   // Fill the cache buffer with the branches in the cache.
   fIsTransferred = kFALSE;

   TTree *tree = ((TBranch*)fBranches->UncheckedAt(0))->GetTree();
   Long64_t entry = tree->GetReadEntry();

   // If the entry is in the range we previously prefetched, there is
   // no point in retrying.   Note that this will also return false
   // during the training phase (fEntryNext is then set intentional to
   // the end of the training phase).
   if (fEntryCurrent <= entry  && entry < fEntryNext) return kFALSE;

   // Triggered by the user, not the learning phase
   if (entry == -1)  entry = 0;

   TTree::TClusterIterator clusterIter = tree->GetClusterIterator(entry);
   fEntryCurrent = clusterIter();
   fEntryNext = clusterIter.GetNextEntry();

   if (fEntryCurrent < fEntryMin) fEntryCurrent = fEntryMin;
   if (fEntryMax <= 0) fEntryMax = tree->GetEntries();
   if (fEntryNext > fEntryMax) fEntryNext = fEntryMax;

   // Check if owner has a TEventList set. If yes we optimize for this
   // Special case reading only the baskets containing entries in the
   // list.
   TEventList *elist = fTree->GetEventList();
   Long64_t chainOffset = 0;
   if (elist) {
      if (fTree->IsA() ==TChain::Class()) {
         TChain *chain = (TChain*)fTree;
         Int_t t = chain->GetTreeNumber();
         chainOffset = chain->GetTreeOffset()[t];
      }
   }

   //clear cache buffer
   TFileCacheRead::Prefetch(0,0);

   //store baskets
   for (Int_t i = 0; i < fNbranches; i++) {
      TBranch *b = (TBranch*)fBranches->UncheckedAt(i);
      if (b->GetDirectory() == 0) continue;
      if (b->GetDirectory()->GetFile() != fFile) continue;
      Int_t nb = b->GetMaxBaskets();
      Int_t *lbaskets   = b->GetBasketBytes();
      Long64_t *entries = b->GetBasketEntry();
      if (!lbaskets || !entries) continue;
      //we have found the branch. We now register all its baskets
      //from the requested offset to the basket below fEntrymax
      Int_t blistsize = b->GetListOfBaskets()->GetSize();
      for (Int_t j=0;j<nb;j++) {
         // This basket has already been read, skip it
         if (j<blistsize && b->GetListOfBaskets()->UncheckedAt(j)) continue;

         Long64_t pos = b->GetBasketSeek(j);
         Int_t len = lbaskets[j];
         if (pos <= 0 || len <= 0) continue;
         //important: do not try to read fEntryNext, otherwise you jump to the next autoflush
         if (entries[j] >= fEntryNext) continue;
         if (entries[j] < entry && (j < nb - 1 && entries[j+1] <= entry)) continue;
         if (elist) {
            Long64_t emax = fEntryMax;
            if (j < nb - 1) emax = entries[j+1] - 1;
            if (!elist->ContainsRange(entries[j] + chainOffset, emax + chainOffset)) continue;
         }
         fNReadPref++;

         TFileCacheRead::Prefetch(pos, len);
      }
      if (gDebug > 0) printf("Entry: %lld, registering baskets branch %s, fEntryNext=%lld, fNseek=%d, fNtot=%d\n", entry, ((TBranch*)fBranches->UncheckedAt(i))->GetName(), fEntryNext, fNseek, fNtot);
   }

   // Now fix the size of the status arrays
   ResetCache();
   fIsLearning = kFALSE;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Change the underlying buffer size of the cache.
/// Returns:
///  - 0 if the buffer content is still available
///  - 1 if some or all of the buffer content has been made unavailable
///  - -1 on error

Int_t TTreeCacheUnzip::SetBufferSize(Int_t buffersize)
{
   Int_t res = TTreeCache::SetBufferSize(buffersize);
   if (res < 0) {
      return res;
   }
   fUnzipBufferSize = Long64_t(fgRelBuffSize * GetBufferSize());
   ResetCache();
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the minimum and maximum entry number to be processed
/// this information helps to optimize the number of baskets to read
/// when prefetching the branch buffers.

void TTreeCacheUnzip::SetEntryRange(Long64_t emin, Long64_t emax)
{
   TTreeCache::SetEntryRange(emin, emax);
}

////////////////////////////////////////////////////////////////////////////////
/// It's the same as TTreeCache::StopLearningPhase but we guarantee that
/// we start the unzipping just after getting the buffers

void TTreeCacheUnzip::StopLearningPhase()
{
   TTreeCache::StopLearningPhase();
}

////////////////////////////////////////////////////////////////////////////////
///update pointer to current Tree and recompute pointers to the branches in the cache

void TTreeCacheUnzip::UpdateBranches(TTree *tree)
{
   TTreeCache::UpdateBranches(tree);
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// From now on we have the methods concerning the threading part of the cache //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Static function that returns the parallel option
/// (to indicate an additional thread)

TTreeCacheUnzip::EParUnzipMode TTreeCacheUnzip::GetParallelUnzip()
{
   return fgParallel;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function that tells wether the multithreading unzipping is activated

Bool_t TTreeCacheUnzip::IsParallelUnzip()
{
   if (fgParallel == kEnable || fgParallel == kForce)
      return kTRUE;

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function that (de)activates multithreading unzipping
///
/// The possible options are:
///  - kEnable _Enable_ it, which causes an automatic detection and launches the
///    additional thread if the number of cores in the machine is greater than
///    one
///  - kDisable _Disable_ will not activate the additional thread.
///  - kForce _Force_ will start the additional thread even if there is only one
///    core. the default will be taken as kEnable.
///
/// Returns 0 if there was an error, 1 otherwise.

Int_t TTreeCacheUnzip::SetParallelUnzip(TTreeCacheUnzip::EParUnzipMode option)
{
   if(fgParallel == kEnable || fgParallel == kForce || fgParallel == kDisable) {
      fgParallel = option;
      return 1;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// From now on we have the methods concerning the unzipping part of the cache //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Read the logical record header from the buffer buf.
/// That must be the pointer tho the header part not the object by itself and
/// must contain data of at least maxbytes
/// Returns nread;
///
/// In output arguments:
///
/// -  nbytes : number of bytes in record
///             if negative, this is a deleted record
///             if 0, cannot read record, wrong value of argument first
/// -  objlen : uncompressed object size
/// -  keylen : length of logical record header
///
/// Note that the arguments objlen and keylen are returned only
/// if maxbytes >=16
/// Note: This was adapted from TFile... so some things dont apply

Int_t TTreeCacheUnzip::GetRecordHeader(char *buf, Int_t maxbytes, Int_t &nbytes, Int_t &objlen, Int_t &keylen)
{
   Version_t versionkey;
   Short_t klen;
   UInt_t datime;
   Int_t nb = 0, olen;
   Int_t nread = maxbytes;
   frombuf(buf, &nb);
   nbytes = nb;
   if (nb < 0) return nread;
   //   const Int_t headerSize = Int_t(sizeof(nb) +sizeof(versionkey) +sizeof(olen) +sizeof(datime) +sizeof(klen));
   const Int_t headerSize = 16;
   if (nread < headerSize) return nread;
   frombuf(buf, &versionkey);
   frombuf(buf, &olen);
   frombuf(buf, &datime);
   frombuf(buf, &klen);
   if (!olen) olen = nbytes - klen;
   objlen = olen;
   keylen = klen;
   return nread;
}

////////////////////////////////////////////////////////////////////////////////
/// This will delete the list of buffers that are in the unzipping cache
/// and will reset certain values in the cache.
/// This name is ambiguos because the method doesn't reset the whole cache,
/// only the part related to the unzipping
/// Note: This method is completely different from TTreeCache::ResetCache(),
/// in that method we were cleaning the prefetching buffer while here we
/// delete the information about the unzipped buffers

void TTreeCacheUnzip::ResetCache()
{
   // Reset all the lists and wipe all the chunks
   fCycle++;
   fUnzipState.Clear(fNseekMax);

   if(fNseekMax < fNseek){
      if (gDebug > 0)
         Info("ResetCache", "Changing fNseekMax from:%d to:%d", fNseekMax, fNseek);

      fUnzipState.Reset(fNseekMax, fNseek);
      fNseekMax = fNseek;
   }
   fEmpty = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// This inflates a basket in the cache.. passing the data to a new
/// buffer that will only wait there to be read...
/// This function is responsible to update corresponding elements in 
/// fUnzipStatus, fUnzipChunks and fUnzipLen. Since we use atomic variables
/// in fUnzipStatus to exclusively unzip the basket, we must update 
/// fUnzipStatus after fUnzipChunks and fUnzipLen and make sure fUnzipChunks
/// and fUnzipLen are ready before main thread fetch the data.

Int_t TTreeCacheUnzip::UnzipCache(Int_t index)
{
   Int_t myCycle;
   const Int_t hlen = 128;
   Int_t objlen = 0, keylen = 0;
   Int_t nbytes = 0;
   Int_t readbuf = 0;

   Long64_t rdoffs = 0;
   Int_t rdlen = 0;

   // To synchronize with the 'paging'
   myCycle = fCycle;
   rdoffs = fSeek[index];
   rdlen = fSeekLen[index];

   Int_t loc = -1;
   if (!fNseek || fIsLearning) {
      return 1;
   }

   if ((myCycle != fCycle) || !fIsTransferred)  {
      fUnzipState.SetFinished(index); // Set it as not done, main thread will take charge
      return 1;
   }

   // Prepare a memory buffer of adequate size
   char* locbuff = 0;
   if (rdlen > 16384) {
      locbuff = new char[rdlen];
   } else if (rdlen * 3 < 16384) {
      locbuff = new char[rdlen * 2];
   } else {
      locbuff = new char[16384];
   }

   readbuf = ReadBufferExt(locbuff, rdoffs, rdlen, loc);

   if (readbuf <= 0) {
      fUnzipState.SetFinished(index); // Set it as not done, main thread will take charge
      if (locbuff) delete [] locbuff;
      return -1;
   }

   GetRecordHeader(locbuff, hlen, nbytes, objlen, keylen);

   Int_t len = (objlen > nbytes - keylen) ? keylen + objlen : nbytes;
   // If the single unzipped chunk is really too big, reset it to not processable
   // I.e. mark it as done but set the pointer to 0
   // This block will be unzipped synchronously in the main thread
   // TODO: ROOT internally breaks zipped buffers into 16MB blocks, we can probably still unzip in parallel.
   if (len > 4 * fUnzipBufferSize) {
           if (gDebug > 0)
                   Info("UnzipCache", "Block %d is too big, skipping.", index);

           fUnzipState.SetFinished(index); // Set it as not done, main thread will take charge
           if (locbuff) delete [] locbuff;
           return 0;
   }

   // Unzip it into a new blk
   char *ptr = 0;
   Int_t loclen = UnzipBuffer(&ptr, locbuff);
   if ((loclen > 0) && (loclen == objlen + keylen)) {
      if ((myCycle != fCycle) || !fIsTransferred) {
         fUnzipState.SetFinished(index); // Set it as not done, main thread will take charge
         if (locbuff) delete [] locbuff;
         return 1;
      }
      fUnzipState.SetUnzipped(index, ptr, loclen); // Set it as done
      fNUnzip++;
   } else {
      fUnzipState.SetFinished(index); // Set it as not done, main thread will take charge
   }

   if (locbuff) delete [] locbuff;
   return 0;
}

#ifdef R__USE_IMT
////////////////////////////////////////////////////////////////////////////////
/// We create a TTaskGroup and asynchronously maps each group of baskets(> 100 kB in total)
/// to a task. In TTaskGroup, we use TThreadExecutor to do the actually work of unzipping 
/// a group of basket. The purpose of creating TTaskGroup is to avoid competing with main thread.

Int_t TTreeCacheUnzip::CreateTasks()
{
   auto mapFunction = [&]() {
      auto unzipFunction = [&](const std::vector<Int_t> &indices) {
         // If cache is invalidated and we should return immediately.
         if (!fIsTransferred) return nullptr;

         for (auto ii : indices) {
            if(fUnzipState.TryUnzipping(ii)) {
               Int_t res = UnzipCache(ii);
               if(res)
                  if (gDebug > 0)
                     Info("UnzipCache", "Unzipping failed or cache is in learning state");
            }
         }
         return nullptr;
      };

      Int_t accusz = 0;
      std::vector<std::vector<Int_t>> basketIndices;
      std::vector<Int_t> indices;
      if (fUnzipGroupSize <= 0) fUnzipGroupSize = 102400;
      for (Int_t i = 0; i < fNseek; i++) {
         while (accusz < fUnzipGroupSize) {
            accusz += fSeekLen[i];
            indices.push_back(i);
            i++;
            if (i >= fNseek) break;
         }
         if (i < fNseek) i--;
         basketIndices.push_back(indices);
         indices.clear();
         accusz = 0;
      }
      ROOT::TThreadExecutor pool;
      pool.Foreach(unzipFunction, basketIndices);
   };

   fUnzipTaskGroup.reset(new ROOT::Experimental::TTaskGroup());
   fUnzipTaskGroup->Run(mapFunction);

   return 0;
}
#endif

////////////////////////////////////////////////////////////////////////////////
/// We try to read a buffer that has already been unzipped
/// Returns -1 in case of read failure, 0 in case it's not in the
/// cache and n>0 in case read from cache (number of bytes copied).
/// pos and len are the original values as were passed to ReadBuffer
/// but instead we will return the inflated buffer.
/// Note!! : If *buf == 0 we will allocate the buffer and it will be the
/// responsability of the caller to free it... it is useful for example
/// to pass it to the creator of TBuffer

Int_t TTreeCacheUnzip::GetUnzipBuffer(char **buf, Long64_t pos, Int_t len, Bool_t *free)
{
   Int_t res = 0;
   Int_t loc = -1;

   // We go straight to TTreeCache/TfileCacheRead, in order to get the info we need
   //  pointer to the original zipped chunk
   //  its index in the original unsorted offsets lists
   //
   // Actually there are situations in which copying the buffer is not
   // useful. But the choice is among doing once more a small memcpy or a binary search in a large array. I prefer the former.
   // Also, here we prefer not to trigger the (re)population of the chunks in the TFileCacheRead. That is
   // better to be done in the main thread.

   Int_t myCycle = fCycle;

   if (fParallel && !fIsLearning) {

      if(fNseekMax < fNseek){
         if (gDebug > 0)
            Info("GetUnzipBuffer", "Changing fNseekMax from:%d to:%d", fNseekMax, fNseek);

         fUnzipState.Reset(fNseekMax, fNseek);
         fNseekMax = fNseek;
      }

      loc = (Int_t)TMath::BinarySearch(fNseek, fSeekSort, pos);
      if ((fCycle == myCycle) && (loc >= 0) && (loc < fNseek) && (pos == fSeekSort[loc])) {

         // The buffer is, at minimum, in the file cache. We must know its index in the requests list
         // In order to get its info
         Int_t seekidx = fSeekIndex[loc];

         do {

            // If the block is ready we get it immediately.
            // And also we don't have to alloc the blks. This is supposed to be
            // the main thread of the app.
            if (fUnzipState.IsUnzipped(seekidx)) {
               if(!(*buf)) {
                  *buf = fUnzipState.fUnzipChunks[seekidx].get();
                  fUnzipState.fUnzipChunks[seekidx].release();
                  *free = kTRUE;
               } else {
                  memcpy(*buf, fUnzipState.fUnzipChunks[seekidx].get(), fUnzipState.fUnzipLen[seekidx]);
                  fUnzipState.fUnzipChunks[seekidx].reset();
                  *free = kFALSE;
               }

               fNFound++;
               return fUnzipState.fUnzipLen[seekidx];
            }

            // If the requested basket is being unzipped by a background task, we try to steal a blk to unzip.
            Int_t reqi = -1;
            
            if (fUnzipState.IsProgress(seekidx)) {
               if (fEmpty) {
                  for (Int_t ii = 0; ii < fNseek; ++ii) {
                     Int_t idx = (seekidx + 1 + ii) % fNseek;
                     if (fUnzipState.IsUntouched(idx)) {
                        if(fUnzipState.TryUnzipping(idx)) {
                           reqi = idx;
                           break;
                        }
                     }
                  }
                  if (reqi < 0) {
                     fEmpty = kFALSE;
                  } else {
                     UnzipCache(reqi);
                  }
               }
 
               if ( myCycle != fCycle ) {
                  if (gDebug > 0)
                     Info("GetUnzipBuffer", "Sudden paging Break!!! fNseek: %d, fIsLearning:%d",
                          fNseek, fIsLearning);

                  seekidx = -1;
                  break;
               }
            }

         } while (fUnzipState.IsProgress(seekidx));

         // Here the block is not pending. It could be done or aborted or not yet being processed.
         if ( (seekidx >= 0) && (fUnzipState.IsUnzipped(seekidx)) ) {
            if(!(*buf)) {
              *buf = fUnzipState.fUnzipChunks[seekidx].get();
               fUnzipState.fUnzipChunks[seekidx].release();
               *free = kTRUE;
            } else {
               memcpy(*buf, fUnzipState.fUnzipChunks[seekidx].get(), fUnzipState.fUnzipLen[seekidx]);
               fUnzipState.fUnzipChunks[seekidx].reset();
               *free = kFALSE;
            }

            fNStalls++;
            return fUnzipState.fUnzipLen[seekidx];
         } else {
            // This is a complete miss. We want to avoid the background tasks
            // to try unzipping this block in the future.
            fUnzipState.SetMissed(seekidx);
         }
      } else {
         loc = -1;
         fIsTransferred = kFALSE;
      }
   }

   if (len > fCompBufferSize) {
      if(fCompBuffer) delete [] fCompBuffer;
      fCompBuffer = new char[len];
      fCompBufferSize = len;
   } else {
      if (fCompBufferSize > len * 4) {
         if(fCompBuffer) delete [] fCompBuffer;
         fCompBuffer = new char[len*2];
         fCompBufferSize = len * 2;
      }
   }

   res = 0;
   if (!ReadBufferExt(fCompBuffer, pos, len, loc)) {
      // Cache is invalidated and we need to wait for all unzipping tasks to befinished before fill new baskets in cache.
#ifdef R__USE_IMT
      if(ROOT::IsImplicitMTEnabled() && fUnzipTaskGroup) {
         fUnzipTaskGroup->Cancel();
         fUnzipTaskGroup.reset();
      }
#endif
      {
         // Fill new baskets into cache.
         R__LOCKGUARD(fIOMutex.get());
	      fFile->Seek(pos);
	      res = fFile->ReadBuffer(fCompBuffer, len);
      } // end of lock scope
#ifdef R__USE_IMT
      if(ROOT::IsImplicitMTEnabled()) {
         CreateTasks();
      }
#endif
   }

   if (res) res = -1;

   if (!res) {
      res = UnzipBuffer(buf, fCompBuffer);
      *free = kTRUE;
   }

   if (!fIsLearning) {
      fNMissed++;
   }
   
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// static function: Sets the unzip relatibe buffer size

void TTreeCacheUnzip::SetUnzipRelBufferSize(Float_t relbufferSize)
{
   fgRelBuffSize = relbufferSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the size for the unzipping cache... by default it should be
/// two times the size of the prefetching cache

void TTreeCacheUnzip::SetUnzipBufferSize(Long64_t bufferSize)
{
   fUnzipBufferSize = bufferSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Unzips a ROOT specific buffer... by reading the header at the beginning.
/// returns the size of the inflated buffer or -1 if error
/// Note!! : If *dest == 0 we will allocate the buffer and it will be the
/// responsability of the caller to free it... it is useful for example
/// to pass it to the creator of TBuffer
/// src is the original buffer with the record (header+compressed data)
/// *dest is the inflated buffer (including the header)

Int_t TTreeCacheUnzip::UnzipBuffer(char **dest, char *src)
{
   Int_t  uzlen = 0;
   Bool_t alloc = kFALSE;

   // Here we read the header of the buffer
   const Int_t hlen = 128;
   Int_t nbytes = 0, objlen = 0, keylen = 0;
   GetRecordHeader(src, hlen, nbytes, objlen, keylen);

   if (!(*dest)) {
      /* early consistency check */
      UChar_t *bufcur = (UChar_t *) (src + keylen);
      Int_t nin, nbuf;
      if(objlen > nbytes - keylen && R__unzip_header(&nin, bufcur, &nbuf) != 0) {
         Error("UnzipBuffer", "Inconsistency found in header (nin=%d, nbuf=%d)", nin, nbuf);
         uzlen = -1;
         return uzlen;
      }
      Int_t l = keylen + objlen;
      *dest = new char[l];
      alloc = kTRUE;
   }
   // Must unzip the buffer
   // fSeekPos[ind]; adress of zipped buffer
   // fSeekLen[ind]; len of the zipped buffer
   // &fBuffer[fSeekPos[ind]]; memory address

   // This is similar to TBasket::ReadBasketBuffers
   Bool_t oldCase = objlen == nbytes - keylen
      && ((TBranch*)fBranches->UncheckedAt(0))->GetCompressionLevel() != 0
      && fFile->GetVersion() <= 30401;

   if (objlen > nbytes-keylen || oldCase) {

      // Copy the key
      memcpy(*dest, src, keylen);
      uzlen += keylen;

      char *objbuf = *dest + keylen;
      UChar_t *bufcur = (UChar_t *) (src + keylen);
      Int_t nin, nbuf;
      Int_t nout = 0;
      Int_t noutot = 0;

      while (1) {
         Int_t hc = R__unzip_header(&nin, bufcur, &nbuf);
         if (hc != 0) break;
         if (gDebug > 2)
            Info("UnzipBuffer", " nin:%d, nbuf:%d, bufcur[3] :%d, bufcur[4] :%d, bufcur[5] :%d ",
                 nin, nbuf, bufcur[3], bufcur[4], bufcur[5]);
         if (oldCase && (nin > objlen || nbuf > objlen)) {
            if (gDebug > 2)
               Info("UnzipBuffer", "oldcase objlen :%d ", objlen);

            //buffer was very likely not compressed in an old version
            memcpy(*dest + keylen, src + keylen, objlen);
            uzlen += objlen;
            return uzlen;
         }

         R__unzip(&nin, bufcur, &nbuf, objbuf, &nout);

         if (gDebug > 2)
            Info("UnzipBuffer", "R__unzip nin:%d, bufcur:%p, nbuf:%d, objbuf:%p, nout:%d",
                 nin, bufcur, nbuf, objbuf, nout);

         if (!nout) break;
         noutot += nout;
         if (noutot >= objlen) break;
         bufcur += nin;
         objbuf += nout;
      }

      if (noutot != objlen) {
         Error("UnzipBuffer", "nbytes = %d, keylen = %d, objlen = %d, noutot = %d, nout=%d, nin=%d, nbuf=%d",
               nbytes,keylen,objlen, noutot,nout,nin,nbuf);
         uzlen = -1;
         if(alloc) delete [] *dest;
         *dest = 0;
         return uzlen;
      }
      uzlen += objlen;
   } else {
      memcpy(*dest, src, keylen);
      uzlen += keylen;
      memcpy(*dest + keylen, src + keylen, objlen);
      uzlen += objlen;
   }
   return uzlen;
}

////////////////////////////////////////////////////////////////////////////////

void  TTreeCacheUnzip::Print(Option_t* option) const {

   printf("******TreeCacheUnzip statistics for file: %s ******\n",fFile->GetName());
   printf("Max allowed mem for pending buffers: %lld\n", fUnzipBufferSize);
   printf("Number of blocks unzipped by threads: %d\n", fNUnzip);
   printf("Number of hits: %d\n", fNFound);
   printf("Number of stalls: %d\n", fNStalls);
   printf("Number of misses: %d\n", fNMissed);

   TTreeCache::Print(option);
}

////////////////////////////////////////////////////////////////////////////////

Int_t TTreeCacheUnzip::ReadBufferExt(char *buf, Long64_t pos, Int_t len, Int_t &loc) {
   R__LOCKGUARD(fIOMutex.get());
   return TTreeCache::ReadBufferExt(buf, pos, len, loc);
}
