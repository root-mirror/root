/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**
 \class TFileCacheRead
 \ingroup IO

 A cache when reading files over the network.

 A caching system to speed up network I/O, i.e. when there is
 no operating system caching support (like the buffer cache for
 local disk I/O). The cache makes sure that every I/O is done with
 a (large) fixed length buffer thereby avoiding many small I/O's.
 Currently the read cache system is used by the classes TNetFile,
 TXNetFile and TWebFile (via TFile::ReadBuffers()).
 When processing TTree, TChain, a specialized class TTreeCache that
 derives from this class is automatically created.
*/

#include "TEnv.h"
#include "TFile.h"
#include "TFileCacheRead.h"
#include "TFileCacheWrite.h"
#include "TFilePrefetch.h"
#include "TMathBase.h"
#include "TFileBufferRead.h"

ClassImp(TFileCacheRead);

////////////////////////////////////////////////////////////////////////////////
/// Default Constructor.

TFileCacheRead::TFileCacheRead() : TObject(), fEnableBuffering(false)
{
   fBufferSizeMin = 0;
   fBufferSize  = 0;
   fBufferLen   = 0;
   fBytesRead   = 0;
   fNoCacheBytesRead = 0;
   fBytesReadExtra = 0;
   fReadCalls   = 0;
   fNoCacheReadCalls = 0;
   fNseek       = 0;
   fNtot        = 0;
   fNb          = 0;
   fSeekSize    = 0;
   fSeek        = 0;
   fSeekIndex   = 0;
   fSeekSort    = 0;
   fPos         = 0;
   fSeekLen     = 0;
   fSeekSortLen = 0;
   fSeekPos     = 0;
   fLen         = 0;
   fFile        = 0;
   fBuffer      = 0;
   fIsSorted    = kFALSE;
   fIsTransferred = kFALSE;

   //values for the second prefetched block
   fBNseek       = 0;
   fBNtot        = 0;
   fBNb          = 0;
   fBSeekSize    = 0;
   fBSeek        = 0;
   fBSeekSort    = 0;
   fBSeekIndex   = 0;
   fBPos         = 0;
   fBSeekLen     = 0;
   fBSeekSortLen = 0;
   fBSeekPos     = 0;
   fBLen         = 0;
   fBIsSorted    = kFALSE;
   fBIsTransferred=kFALSE;

   fAsyncReading = kFALSE;
   fEnablePrefetching = kFALSE;
   fPrefetch        = 0;
   fPrefetchedBlocks= 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a TFileCacheRead data structure.

TFileCacheRead::TFileCacheRead(TFile *file, Int_t buffersize, TObject *tree)
           : TObject(), fEnableBuffering(false)
{
   if (buffersize <=10000) fBufferSize = 100000;
   else fBufferSize = buffersize;

   fBufferSizeMin = fBufferSize;
   fBufferLen   = 0;
   fBytesRead   = 0;
   fNoCacheBytesRead = 0;
   fBytesReadExtra = 0;
   fReadCalls   = 0;
   fNoCacheReadCalls = 0;
   fNseek       = 0;
   fNtot        = 0;
   fNb          = 0;
   fSeekSize    = 10000;
   fSeek        = new Long64_t[fSeekSize];
   fSeekIndex   = new Int_t[fSeekSize];
   fSeekSort    = new Long64_t[fSeekSize];
   fPos         = new Long64_t[fSeekSize];
   fSeekLen     = new Int_t[fSeekSize];
   fSeekSortLen = new Int_t[fSeekSize];
   fSeekPos     = new Int_t[fSeekSize];
   fLen         = new Int_t[fSeekSize];
   fFile        = file;

   //initialisation for the second block
   fBNseek       = 0;
   fBNtot        = 0;
   fBNb          = 0;
   fBSeekSize    = 10000;
   fBSeek        = new Long64_t[fBSeekSize];
   fBSeekIndex   = new Int_t[fBSeekSize];
   fBSeekSort    = new Long64_t[fBSeekSize];
   fBPos         = new Long64_t[fBSeekSize];
   fBSeekLen     = new Int_t[fBSeekSize];
   fBSeekSortLen = new Int_t[fBSeekSize];
   fBSeekPos     = new Int_t[fBSeekSize];
   fBLen         = new Int_t[fBSeekSize];

   fBuffer = 0;
   fPrefetch = 0;
   fPrefetchedBlocks = 0;

   //initialise the prefetch object and set the cache directory
   // start the thread only if the file is not local
   fEnablePrefetching = gEnv->GetValue("TFile.AsyncPrefetching", 0);
   fEnableBuffering   = gEnv->GetValue("TFile.LocalBuffer", 0);
   bool remoteFile = strcmp(file->GetEndpointUrl()->GetProtocol(), "file");

   if (fEnablePrefetching && remoteFile) {
      SetEnablePrefetchingImpl(true);
   } else if (file && fEnableBuffering && remoteFile) {
      SetEnablePrefetchingImpl(false);
      SetEnableBufferRead(true);
   }
   else { //disable the async pref for local files
      SetEnablePrefetchingImpl(false);
   }

   fIsSorted       = kFALSE;
   fIsTransferred  = kFALSE;
   fBIsSorted      = kFALSE;
   fBIsTransferred = kFALSE;

   if (file) file->SetCacheRead(this, tree);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TFileCacheRead::~TFileCacheRead()
{
   SafeDelete(fPrefetch);
   delete [] fSeek;
   delete [] fSeekIndex;
   delete [] fSeekSort;
   delete [] fPos;
   delete [] fSeekLen;
   delete [] fSeekSortLen;
   delete [] fSeekPos;
   delete [] fLen;
   if (fBuffer)
      delete [] fBuffer;
   delete [] fBSeek;
   delete [] fBSeekIndex;
   delete [] fBSeekSort;
   delete [] fBPos;
   delete [] fBSeekLen;
   delete [] fBSeekSortLen;
   delete [] fBSeekPos;
   delete [] fBLen;
}

////////////////////////////////////////////////////////////////////////////////
/// Close out any threads or asynchronous fetches used by the underlying
/// implementation.
/// This is called by TFile::Close to prevent usage of the file handles
/// after the closing of the file.

void TFileCacheRead::Close(Option_t * /* opt = "" */)
{
   if (fPrefetch) {
      delete fPrefetch;
      fPrefetch = 0;
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Add block of length len at position pos in the list of blocks to
/// be prefetched. If pos <= 0 the current blocks (if any) are reset.

void TFileCacheRead::Prefetch(Long64_t pos, Int_t len)
{
   fIsSorted = kFALSE;
   fIsTransferred = kFALSE;
   if (pos <= 0) {
      fNseek = 0;
      fNtot  = 0;
      return;
   }
   if (fNseek >= fSeekSize) {
      //reallocate buffers
      fSeekSize *= 2;
      Long64_t *aSeek        = new Long64_t[fSeekSize];
      Int_t    *aSeekIndex   = new Int_t[fSeekSize];
      Long64_t *aSeekSort    = new Long64_t[fSeekSize];
      Long64_t *aPos         = new Long64_t[fSeekSize];
      Int_t    *aSeekLen     = new Int_t[fSeekSize];
      Int_t    *aSeekSortLen = new Int_t[fSeekSize];
      Int_t    *aSeekPos     = new Int_t[fSeekSize];
      Int_t    *aLen         = new Int_t[fSeekSize];
      for (Int_t i=0;i<fNseek;i++) {
         aSeek[i]        = fSeek[i];
         aSeekIndex[i]   = fSeekIndex[i];
         aSeekSort[i]    = fSeekSort[i];
         aPos[i]         = fPos[i];
         aSeekLen[i]     = fSeekLen[i];
         aSeekSortLen[i] = fSeekSortLen[i];
         aSeekPos[i]     = fSeekPos[i];
         aLen[i]         = fLen[i];
      }
      delete [] fSeek;
      delete [] fSeekIndex;
      delete [] fSeekSort;
      delete [] fPos;
      delete [] fSeekLen;
      delete [] fSeekSortLen;
      delete [] fSeekPos;
      delete [] fLen;
      fSeek        = aSeek;
      fSeekIndex   = aSeekIndex;
      fSeekSort    = aSeekSort;
      fPos         = aPos;
      fSeekLen     = aSeekLen;
      fSeekSortLen = aSeekSortLen;
      fSeekPos     = aSeekPos;
      fLen         = aLen;
   }

   fSeek[fNseek] = pos;
   fSeekLen[fNseek] = len;
   fNseek++;
   fNtot += len;
}


////////////////////////////////////////////////////////////////////////////////

void TFileCacheRead::SecondPrefetch(Long64_t pos, Int_t len){
   //add a new element and increase the size if necessary
   fBIsSorted = kFALSE;
   if (pos <= 0) {
      fBNseek = 0;
      fBNtot  = 0;
      return;
   }
   if (fBNseek >= fBSeekSize) {
      //reallocate buffers
      fBSeekSize *= 2;
      Long64_t *aSeek        = new Long64_t[fBSeekSize];
      Int_t    *aSeekIndex   = new Int_t[fBSeekSize];
      Long64_t *aSeekSort    = new Long64_t[fBSeekSize];
      Long64_t *aPos         = new Long64_t[fBSeekSize];
      Int_t    *aSeekLen     = new Int_t[fBSeekSize];
      Int_t    *aSeekSortLen = new Int_t[fBSeekSize];
      Int_t    *aSeekPos     = new Int_t[fBSeekSize];
      Int_t    *aLen         = new Int_t[fBSeekSize];
      for (Int_t i=0;i<fBNseek;i++) {
         aSeek[i]        = fBSeek[i];
         aSeekIndex[i]   = fBSeekIndex[i];
         aSeekSort[i]    = fBSeekSort[i];
         aPos[i]         = fBPos[i];
         aSeekLen[i]     = fBSeekLen[i];
         aSeekSortLen[i] = fBSeekSortLen[i];
         aSeekPos[i]     = fBSeekPos[i];
         aLen[i]         = fBLen[i];
      }
      delete [] fBSeek;
      delete [] fBSeekIndex;
      delete [] fBSeekSort;
      delete [] fBPos;
      delete [] fBSeekLen;
      delete [] fBSeekSortLen;
      delete [] fBSeekPos;
      delete [] fBLen;
      fBSeek        = aSeek;
      fBSeekIndex   = aSeekIndex;
      fBSeekSort    = aSeekSort;
      fBPos         = aPos;
      fBSeekLen     = aSeekLen;
      fBSeekSortLen = aSeekSortLen;
      fBSeekPos     = aSeekPos;
      fBLen         = aLen;
   }

   fBSeek[fBNseek] = pos;
   fBSeekLen[fBNseek] = len;
   fBNseek++;
   fBNtot += len;
}


////////////////////////////////////////////////////////////////////////////////
/// Print cache statistics.
///
/// The format is:
///     ******TreeCache statistics for file: cms2.root ******
///     Reading............................: 72761843 bytes in 7 transactions
///     Readahead..........................: 256000 bytes with overhead = 0 bytes
///     Average transaction................: 10394.549000 Kbytes
///     Number of blocks in current cache..: 210, total size: 6280352
///
/// If option = "a" the list of blocks in the cache is printed
/// NB: this function is automatically called by TTreeCache::Print

void TFileCacheRead::Print(Option_t *option) const
{
   TString opt = option;
   opt.ToLower();
   printf("Cached Reading.....................: %lld bytes in %d transactions\n",this->GetBytesRead(), this->GetReadCalls());
   printf("Reading............................: %lld bytes in %d uncached transactions\n",this->GetNoCacheBytesRead(), this->GetNoCacheReadCalls());
   printf("Readahead..........................: %d bytes with overhead = %lld bytes\n",TFile::GetReadaheadSize(),this->GetBytesReadExtra());
   if (this->GetReadCalls() > 0)
      printf("Average transaction................: %f Kbytes\n",0.001*Double_t(this->GetBytesRead())/Double_t(this->GetReadCalls()));
   else
      printf("Average transaction................: No read calls yet\n");
   printf("Number of blocks in current cache..: %d, total size: %d\n",fNseek,fNtot);
   if (fPrefetch){
     printf("Prefetching .......................: %lli blocks\n", fPrefetchedBlocks);
     printf("Prefetching Wait Time..............: %f seconds\n", fPrefetch->GetWaitTime() / 1e+6);
   }
   if (fBufferFile) {
     printf("Buffer file chunks ................: %u blocks\n", fBufferFile->GetCount());
   }

   if (!opt.Contains("a")) return;
   for (Int_t i=0;i<fNseek;i++) {
      if (fIsSorted && !opt.Contains("s")) {
         printf("block: %5d, from: %lld to %lld, len = %d bytes\n",i,fSeekSort[i],fSeekSort[i]+fSeekSortLen[i],fSeekSortLen[i]);
      } else {
         printf("block: %5d, from: %lld to %lld, len = %d bytes\n",i,fSeek[i],fSeek[i]+fSeekLen[i],fSeekLen[i]);
      }
   }
   printf ("Number of long buffers = %d\n",fNb);
   for (Int_t j=0;j<fNb;j++) {
      printf("fPos[%d] = %lld, fLen = %d\n",j,fPos[j],fLen[j]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read buffer at position pos.
///
/// If pos is in the list of prefetched blocks read from fBuffer,
/// otherwise need to make a normal read from file. Returns -1 in case of
/// read error, 0 in case not in cache, 1 in case read from cache.

Int_t TFileCacheRead::ReadBuffer(char *buf, Long64_t pos, Int_t len)
{
   Long64_t fileBytesRead0 = fFile->GetBytesRead();
   Long64_t fileBytesReadExtra0 = fFile->GetBytesReadExtra();
   Int_t fileReadCalls0 = fFile->GetReadCalls();

   Int_t loc = -1;
   Int_t rc = ReadBufferExt(buf, pos, len, loc);

   fBytesRead += fFile->GetBytesRead() - fileBytesRead0;
   fBytesReadExtra += fFile->GetBytesReadExtra() - fileBytesReadExtra0;
   fReadCalls += fFile->GetReadCalls() - fileReadCalls0;

   return rc;
}

////////////////////////////////////////////////////////////////////////////////

Int_t TFileCacheRead::ReadBufferExt(char *buf, Long64_t pos, Int_t len, Int_t &loc)
{
   if (fEnablePrefetching) {
      return ReadBufferExtPrefetch(buf, pos, len, loc);
   } else if (fBufferFile) {
      auto retval = fBufferFile->Pread(buf, len, pos);
      return (retval == -1) ? retval : (len == retval);
   } else {
      return ReadBufferExtNormal(buf, pos, len, loc);
   }
}


////////////////////////////////////////////////////////////////////////////////
///prefetch the first block

Int_t TFileCacheRead::ReadBufferExtPrefetch(char *buf, Long64_t pos, Int_t len, Int_t &loc)
{
   if (fNseek > 0 && !fIsSorted) {
      Sort();
      loc = -1;
      fPrefetch->ReadBlock(fPos, fLen, fNb);
      fPrefetchedBlocks++;
      fIsTransferred = kTRUE;
   }

   //try to prefetch the second block
   if (fBNseek > 0 && !fBIsSorted) {
      SecondSort();
      loc = -1;
      fPrefetch->ReadBlock(fBPos, fBLen, fBNb);
      fPrefetchedBlocks++;
   }

   // in case we are writing and reading to/from this file, we must check
   // if this buffer is in the write cache (not yet written to the file)
   if (TFileCacheWrite *cachew = fFile->GetCacheWrite()) {
      if (cachew->ReadBuffer(buf,pos,len) == 0) {
         fFile->SetOffset(pos+len);
         return 1;
      }
   }

   // try to prefetch from the first block
   if (loc < 0) {
      loc = (Int_t)TMath::BinarySearch(fNseek,fSeekSort,pos);
   }

   if (loc >= 0 && loc < fNseek && pos == fSeekSort[loc]) {
      if (buf && fPrefetch){
         // prefetch with the new method
         fPrefetch->ReadBuffer(buf, pos, len);
         return 1;
      }
   }
   else if (buf && fPrefetch){
      // try to preferch from the second block
      loc = (Int_t)TMath::BinarySearch(fBNseek, fBSeekSort, pos);

      if (loc >= 0 && loc < fBNseek && pos == fBSeekSort[loc]){
         if (fPrefetch->ReadBuffer(buf, pos, len)) {
           return 1;
        }
      }
   }

   return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Base function for ReadBuffer.
///
/// Also gives out the position of the block in the internal buffer.
/// This helps TTreeCacheUnzip to avoid doing twice the binary search.

Int_t TFileCacheRead::ReadBufferExtNormal(char *buf, Long64_t pos, Int_t len, Int_t &loc)
{
   if (fNseek > 0 && !fIsSorted) {
      Sort();
      loc = -1;

      // If ReadBufferAsync is not supported by this implementation...
      if (!fAsyncReading) {
         // Then we use the vectored read to read everything now
         if (fFile->ReadBuffers(fBuffer,fPos,fLen,fNb)) {
            return -1;
         }
         fIsTransferred = kTRUE;
      } else {
         // In any case, we'll start to request the chunks.
         // This implementation simply reads all the chunks in advance
         // in the async way.

         // Use the async readv instead of single reads
         fFile->ReadBuffers(0, 0, 0, 0); //Clear the XrdClient cache
         if (fFile->ReadBuffers(0,fPos,fLen,fNb)) {
            return -1;
         }
         fIsTransferred = kTRUE;
      }
   }

   // in case we are writing and reading to/from this file, we much check
   // if this buffer is in the write cache (not yet written to the file)
   if (TFileCacheWrite *cachew = fFile->GetCacheWrite()) {
      if (cachew->ReadBuffer(buf,pos,len) == 0) {
         fFile->SetOffset(pos+len);
         return 1;
      }
   }

   // If asynchronous reading is supported by this implementation...
   if (fAsyncReading) {

         // Now we dont have to look for it in the local buffer
         // if it's async, we expect that the communication library
         // will handle it more efficiently than we can do here

      Int_t retval;
      if (loc < 0)
         loc = (Int_t)TMath::BinarySearch(fNseek,fSeekSort,pos);

      // We use the internal list just to notify if the list is to be reconstructed
      if (loc >= 0 && loc < fNseek && pos == fSeekSort[loc]) {
         // Block found, the caller will get it

         if (buf) {
            // disable cache to avoid infinite recursion
            if (fFile->ReadBuffer(buf, pos, len)) {
               return -1;
            }
            fFile->SetOffset(pos+len);
         }

         retval = 1;
      } else {
         // Block not found in the list, we report it as a miss
         retval = 0;
      }

      if (gDebug > 0)
         Info("ReadBuffer","pos=%lld, len=%d, retval=%d, loc=%d, "
              "fseekSort[loc]=%lld, fSeekLen[loc]=%d",
              pos, len, retval, loc, fSeekSort[loc], fSeekLen[loc]);

      return retval;
   } else {

      if (loc < 0)
         loc = (Int_t)TMath::BinarySearch(fNseek, fSeekSort, pos);

      if (loc >= 0 && loc <fNseek && pos == fSeekSort[loc]) {
         if (buf) {
            memcpy(buf,&fBuffer[fSeekPos[loc]],len);
            fFile->SetOffset(pos+len);
         }
         return 1;
      }
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the file using this cache and reset the current blocks (if any).

void TFileCacheRead::SetFile(TFile *file, TFile::ECacheAction action)
{
   fFile = file;

   if (fAsyncReading) {
      // If asynchronous reading is not supported by this TFile specialization
      // we use sync primitives, hence we need the local buffer
      if (file && file->ReadBufferAsync(0, 0)) {
         fAsyncReading = kFALSE;
         fBuffer       = new char[fBufferSize];
      }
   }

   if (action == TFile::kDisconnect)
      Prefetch(0,0);

   if (fPrefetch) {
      if (action == TFile::kDisconnect)
         SecondPrefetch(0, 0);
      fPrefetch->SetFile(file, action);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Sort buffers to be prefetched in increasing order of positions.
/// Merge consecutive blocks if necessary.

void TFileCacheRead::Sort()
{
   if (!fNseek) return;
   TMath::Sort(fNseek,fSeek,fSeekIndex,kFALSE);
   Int_t i;
   Int_t nb = 0;
   Int_t effectiveNseek = 0;
   for (i=0;i<fNseek;i++) {
      // Skip duplicates
      Int_t ind = fSeekIndex[i];
      if (effectiveNseek!=0 && fSeek[ind]==fSeekSort[effectiveNseek-1])
      {
         if (fSeekSortLen[effectiveNseek-1] < fSeekLen[ind]) {
            fSeekSortLen[effectiveNseek-1] = fSeekLen[ind];
         }
         continue;
      }
      fSeekSort[effectiveNseek] = fSeek[ind];
      fSeekSortLen[effectiveNseek] = fSeekLen[ind];
      ++effectiveNseek;
   }
   fNseek = effectiveNseek;
   if (fNtot > fBufferSizeMin) {
      fBufferSize = fNtot + 100;
      delete [] fBuffer;
      fBuffer = 0;
      // If ReadBufferAsync is not supported by this implementation
      // it means that we are using sync primitives, hence we need the local buffer
      if (!fAsyncReading)
         fBuffer = new char[fBufferSize];
   }
   fPos[0]  = fSeekSort[0];
   fLen[0]  = fSeekSortLen[0];
   fSeekPos[0] = 0;
   for (i=1;i<fNseek;i++) {
      fSeekPos[i] = fSeekPos[i-1] + fSeekSortLen[i-1];
      //in the test below 16 MBytes is pure empirirical and may depend on the file system.
      //increasing this number must be done with care, as it may increase
      //the job real time (mismatch with OS buffers)
      if ((fSeekSort[i] != fSeekSort[i-1]+fSeekSortLen[i-1]) ||
          (fLen[nb] > 16000000)) {
         nb++;
         fPos[nb] = fSeekSort[i];
         fLen[nb] = fSeekSortLen[i];
      } else {
         fLen[nb] += fSeekSortLen[i];
      }
   }
   fNb = nb+1;
   fIsSorted = kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Sort buffers to be prefetched in increasing order of positions.
///
/// Merge consecutive blocks if necessary.

void TFileCacheRead::SecondSort()
{
   if (!fBNseek) return;
   TMath::Sort(fBNseek,fBSeek,fBSeekIndex,kFALSE);
   Int_t i;
   Int_t nb = 0;
   Int_t effectiveNseek = 0;
   for (i=0;i<fBNseek;i++) {
      // Skip duplicates
      Int_t ind = fBSeekIndex[i];
      if (effectiveNseek!=0 && fBSeek[ind]==fBSeekSort[effectiveNseek-1])
      {
         if (fBSeekSortLen[effectiveNseek-1] < fBSeekLen[ind]) {
            fBSeekSortLen[effectiveNseek-1] = fBSeekLen[ind];
         }
         continue;
      }
      fBSeekSort[effectiveNseek] = fBSeek[ind];
      fBSeekSortLen[effectiveNseek] = fBSeekLen[ind];
      ++effectiveNseek;
   }
   fBNseek = effectiveNseek;
   if (fBNtot > fBufferSizeMin) {
      fBufferSize = fBNtot + 100;
      delete [] fBuffer;
      fBuffer = 0;
      // If ReadBufferAsync is not supported by this implementation
      // it means that we are using sync primitives, hence we need the local buffer
      if (!fAsyncReading)
         fBuffer = new char[fBufferSize];
   }
   fBPos[0]  = fBSeekSort[0];
   fBLen[0]  = fBSeekSortLen[0];
   fBSeekPos[0] = 0;
   for (i=1;i<fBNseek;i++) {
      fBSeekPos[i] = fBSeekPos[i-1] + fBSeekSortLen[i-1];
      //in the test below 16 MBytes is pure empirirical and may depend on the file system.
      //increasing this number must be done with care, as it may increase
      //the job real time (mismatch with OS buffers)
      if ((fBSeekSort[i] != fBSeekSort[i-1]+fBSeekSortLen[i-1]) ||
         (fBLen[nb] > 16000000)) {
         nb++;
         fBPos[nb] = fBSeekSort[i];
         fBLen[nb] = fBSeekSortLen[i];
      } else {
         fBLen[nb] += fBSeekSortLen[i];
      }
   }
   fBNb = nb+1;
   fBIsSorted = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////

TFilePrefetch* TFileCacheRead::GetPrefetchObj(){
   return this->fPrefetch;
}


////////////////////////////////////////////////////////////////////////////////

void TFileCacheRead::WaitFinishPrefetch()
{
   if ( fEnablePrefetching && fPrefetch ) {
      fPrefetch->WaitFinishPrefetch();
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Sets the buffer size.
///
/// If the current prefetch list is too large to fit in
/// the new buffer some or all of the prefetch blocks are dropped. The
/// requested buffersize must be greater than zero.
/// Return values:
///   - 0 if the prefetch block lists remain unchanged
///   - 1 if some or all blocks have been removed from the prefetch list
///   - -1 on error

Int_t TFileCacheRead::SetBufferSize(Int_t buffersize)
{
   if (buffersize <= 0) return -1;
   if (buffersize <=10000) buffersize = 100000;

   if (buffersize == fBufferSize) {
      fBufferSizeMin = buffersize;
      return 0;
   }

   Bool_t inval = kFALSE;

   // the cached data is too large to fit in the new buffer size mark data unavailable
   if (fNtot > buffersize) {
      Prefetch(0, 0);
      inval = kTRUE;
   }
   if (fBNtot > buffersize) {
      SecondPrefetch(0, 0);
      inval = kTRUE;
   }

   char *np = 0;
   if (!fEnablePrefetching && !fAsyncReading) {
      char *pres = 0;
      if (fIsTransferred) {
         // will need to preserve buffer data
         pres = fBuffer;
         fBuffer = 0;
      }
      delete [] fBuffer;
      fBuffer = 0;
      np = new char[buffersize];
      if (pres) {
         memcpy(np, pres, fNtot);
      }
      delete [] pres;
   }

   delete [] fBuffer;
   fBuffer = np;
   fBufferSizeMin = buffersize;
   fBufferSize = buffersize;

   if (inval) {
      return 1;
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the buffer mode of this file.
///
/// If setBuffer is `true`, then the TFileCacheRead will store a remote file's
/// contents on a local disk directory (can be controlled by Cache.Directory).
/// When the TFile is closed, the buffer is automatically cleaned up by the
/// operating system.
///
void TFileCacheRead::SetEnableBufferRead(bool setBuffer) {
   if (setBuffer && !fBufferFile) {fBufferFile.reset(new TFileBufferRead(fFile));}
   else if (!setBuffer) {fBufferFile.release();}
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of blocks fetched to the local disk buffer.
///
/// If local disk buffering is disabled, this returns 0.
///
Long_t TFileCacheRead::GetBufferedBlocks() const {
   return fBufferFile ? fBufferFile->GetCount() : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the prefetching mode of this file.
///
/// If 'setPrefetching', enable the asynchronous prefetching
/// (using TFilePrefetch) and if the gEnv and rootrc
/// variable Cache.Directory is set, also enable the local
/// caching of the prefetched blocks.
/// if 'setPrefetching', the old prefetcher is enabled is
/// the gEnv and rootrc variable is TFile.AsyncReading

void TFileCacheRead::SetEnablePrefetching(Bool_t setPrefetching)
{
   SetEnablePrefetchingImpl(setPrefetching);
}

////////////////////////////////////////////////////////////////////////////////
/// TFileCacheRead implementation of SetEnablePrefetching.
///
/// This function is called from the constructor and should not be virtual.

void TFileCacheRead::SetEnablePrefetchingImpl(Bool_t setPrefetching)
{
   fEnablePrefetching = setPrefetching;

   if (!fPrefetch && fEnablePrefetching) {
      fPrefetch = new TFilePrefetch(fFile);
      const char* cacheDir = gEnv->GetValue("Cache.Directory", "");
      if (strcmp(cacheDir, ""))
        if (!fPrefetch->SetCache((char*) cacheDir))
           fprintf(stderr, "Error while trying to set the cache directory: %s.\n", cacheDir);
      if (fPrefetch->ThreadStart()){
         fprintf(stderr,"Error stating prefetching thread. Disabling prefetching.\n");
         fEnablePrefetching = 0;
      }
   } else if (fPrefetch && !fEnablePrefetching) {
       SafeDelete(fPrefetch);
       fPrefetch = NULL;
   }

   //environment variable used to switch to the new method of reading asynchronously
   if (fEnablePrefetching) {
      fAsyncReading = kFALSE;
   }
   else {
      fAsyncReading = gEnv->GetValue("TFile.AsyncReading", 0);
      if (fAsyncReading) {
         // Check if asynchronous reading is supported by this TFile specialization
         fAsyncReading = kFALSE;
         if (fFile && !(fFile->ReadBufferAsync(0, 0)))
            fAsyncReading = kTRUE;
         }
      if (!fAsyncReading && fBuffer == 0) {
         // we use sync primitives, hence we need the local buffer
         fBuffer = new char[fBufferSize];
      }
   }
}

