// @(#)root/tree:$Id: 4e77188fbf1e7fd026a984989de66663c49b12fc $
// Author: Rene Brun   19/01/96
/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <chrono>

#include "TBasket.h"
#include "TBuffer.h"
#include "TBufferFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TMath.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "TVirtualMutex.h"
#include "TVirtualPerfStats.h"
#include "TTimeStamp.h"
#include "ROOT/TIOFeatures.hxx"
#include "RZip.h"

#include <bitset>

const UInt_t kDisplacementMask = 0xFF000000;  // In the streamer the two highest bytes of
                                              // the fEntryOffset are used to stored displacement.

ClassImp(TBasket);

/** \class TBasket
\ingroup tree

Manages buffers for branches of a Tree.

See picture in TTree.
*/

////////////////////////////////////////////////////////////////////////////////
/// Default contructor.

TBasket::TBasket()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor used during reading.

TBasket::TBasket(TDirectory *motherDir) : TKey(motherDir)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Basket normal constructor, used during writing.

TBasket::TBasket(const char *name, const char *title, TBranch *branch)
   : TKey(branch->GetDirectory()), fBufferSize(branch->GetBasketSize()), fNevBufSize(branch->GetEntryOffsetLen()),
     fHeaderOnly(kTRUE), fIOBits(branch->GetIOFeatures().GetFeatures())
{
   SetName(name);
   SetTitle(title);
   fClassName   = "TBasket";
   fBuffer = nullptr;
   fBufferRef   = new TBufferFile(TBuffer::kWrite, fBufferSize);
   fVersion    += 1000;
   if (branch->GetDirectory()) {
      TFile *file = branch->GetFile();
      fBufferRef->SetParent(file);
   }
   if (branch->GetTree()) {
#ifdef R__USE_IMT
      fCompressedBufferRef = branch->GetTransientBuffer(fBufferSize);
#else
      fCompressedBufferRef = branch->GetTree()->GetTransientBuffer(fBufferSize);
#endif
      fOwnsCompressedBuffer = kFALSE;
      if (!fCompressedBufferRef) {
         fCompressedBufferRef = new TBufferFile(TBuffer::kRead, fBufferSize);
         fOwnsCompressedBuffer = kTRUE;
      }
   }
   fBranch = branch;
   Streamer(*fBufferRef);
   fKeylen      = fBufferRef->Length();
   fObjlen      = fBufferSize - fKeylen;
   fLast        = fKeylen;
   fBuffer      = 0;
   fHeaderOnly  = kFALSE;
   if (fNevBufSize) {
      fEntryOffset = new Int_t[fNevBufSize];
      for (Int_t i=0;i<fNevBufSize;i++) fEntryOffset[i] = 0;
   }
   branch->GetTree()->IncrementTotalBuffers(fBufferSize);
}

////////////////////////////////////////////////////////////////////////////////
/// Basket destructor.

TBasket::~TBasket()
{
   if (fDisplacement) delete [] fDisplacement;
   ResetEntryOffset();
   if (fBufferRef) delete fBufferRef;
   fBufferRef = 0;
   fBuffer = 0;
   fDisplacement= 0;
   // Note we only delete the compressed buffer if we own it
   if (fCompressedBufferRef && fOwnsCompressedBuffer) {
      delete fCompressedBufferRef;
      fCompressedBufferRef = 0;
   }
   // TKey::~TKey will use fMotherDir to attempt to remove they key
   // from the directory's list of key.  A basket is never in that list
   // and in some cases (eg. f = new TFile(); TTree t; delete f;) the
   // directory is gone before the TTree.
   fMotherDir = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Increase the size of the current fBuffer up to newsize.

void TBasket::AdjustSize(Int_t newsize)
{
   if (fBuffer == fBufferRef->Buffer()) {
      fBufferRef->Expand(newsize);
      fBuffer = fBufferRef->Buffer();
   } else {
      fBufferRef->Expand(newsize);
   }
   fBranch->GetTree()->IncrementTotalBuffers(newsize-fBufferSize);
   fBufferSize  = newsize;
   fLastWriteBufferSize[0] = newsize;
   fLastWriteBufferSize[1] = 0;
   fLastWriteBufferSize[2] = 0;
   fNextBufferSizeRecord = 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy the basket of this branch onto the file to.

Long64_t TBasket::CopyTo(TFile *to)
{
   fBufferRef->SetWriteMode();
   Int_t nout = fNbytes - fKeylen;
   fBuffer = fBufferRef->Buffer();
   Create(nout, to);
   fBufferRef->SetBufferOffset(0);
   fHeaderOnly = kTRUE;
   Streamer(*fBufferRef);
   fHeaderOnly = kFALSE;
   Int_t nBytes = WriteFileKeepBuffer(to);

   return nBytes>0 ? nBytes : -1;
}

////////////////////////////////////////////////////////////////////////////////
///  Delete fEntryOffset array.

void TBasket::DeleteEntryOffset()
{
   ResetEntryOffset();
   fNevBufSize  = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Drop buffers of this basket if it is not the current basket.

Int_t TBasket::DropBuffers()
{
   if (!fBuffer && !fBufferRef) return 0;

   if (fDisplacement) delete [] fDisplacement;
   ResetEntryOffset();
   if (fBufferRef)    delete fBufferRef;
   if (fCompressedBufferRef && fOwnsCompressedBuffer) delete fCompressedBufferRef;
   fBufferRef   = 0;
   fCompressedBufferRef = 0;
   fBuffer      = 0;
   fDisplacement= 0;
   fEntryOffset = 0;
   fBranch->GetTree()->IncrementTotalBuffers(-fBufferSize);
   return fBufferSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculates the entry offset array, if possible.
///
/// Result is cached, meaning that this should only be invoked once per basket.

Int_t *TBasket::GetCalculatedEntryOffset()
{
   if (fEntryOffset != reinterpret_cast<Int_t *>(-1)) {
      return fEntryOffset;
   }

   if (R__unlikely(!fBranch)) {
      Error("GetCalculatedEntryOffset", "Basket entry offset calculation requested, but no associated TBranch!");
      return nullptr;
   }
   if (R__unlikely(fBranch->GetNleaves() != 1)) {
      Error("GetCalculatedEntryOffset", "Branch contains multiple leaves - unable to calculated entry offsets!");
      return nullptr;
   }
   TLeaf *leaf = static_cast<TLeaf *>((*fBranch->GetListOfLeaves())[0]);
   fEntryOffset = leaf->GenerateOffsetArray(fKeylen, fNevBuf);
   return fEntryOffset;
}

////////////////////////////////////////////////////////////////////////////////
/// Determine whether we can generate the offset array when this branch is read.
///

Bool_t TBasket::CanGenerateOffsetArray()
{
   if (fBranch->GetNleaves() != 1) {
      return kFALSE;
   }
   TLeaf *leaf = static_cast<TLeaf *>((*fBranch->GetListOfLeaves())[0]);
   return leaf->CanGenerateOffsetArray();
}

////////////////////////////////////////////////////////////////////////////////
/// Get pointer to buffer for internal entry.

Int_t TBasket::GetEntryPointer(Int_t entry)
{
   Int_t offset;
   Int_t *entryOffset = GetEntryOffset();
   if (entryOffset)  offset = entryOffset[entry];
   else              offset = fKeylen + entry*fNevBufSize;
   fBufferRef->SetBufferOffset(offset);
   return offset;
}

////////////////////////////////////////////////////////////////////////////////
/// Load basket buffers in memory without unziping.
/// This function is called by TTreeCloner.
/// The function returns 0 in case of success, 1 in case of error.

Int_t TBasket::LoadBasketBuffers(Long64_t pos, Int_t len, TFile *file, TTree *tree)
{
   if (fBufferRef) {
      // Reuse the buffer if it exist.
      fBufferRef->Reset();

      // We use this buffer both for reading and writing, we need to
      // make sure it is properly sized for writing.
      fBufferRef->SetWriteMode();
      if (fBufferRef->BufferSize() < len) {
         fBufferRef->Expand(len);
      }
      fBufferRef->SetReadMode();
   } else {
      fBufferRef = new TBufferFile(TBuffer::kRead, len);
   }
   fBufferRef->SetParent(file);
   char *buffer = fBufferRef->Buffer();
   file->Seek(pos);
   TFileCacheRead *pf = tree->GetReadCache(file);
   if (pf) {
      TVirtualPerfStats* temp = gPerfStats;
      if (tree->GetPerfStats()) gPerfStats = tree->GetPerfStats();
      Int_t st = pf->ReadBuffer(buffer,pos,len);
      if (st < 0) {
         return 1;
      } else if (st == 0) {
         // fOffset might have been changed via TFileCacheRead::ReadBuffer(), reset it
         file->Seek(pos);
         // If we are using a TTreeCache, disable reading from the default cache
         // temporarily, to force reading directly from file
         TTreeCache *fc = dynamic_cast<TTreeCache*>(file->GetCacheRead());
         if (fc) fc->Disable();
         Int_t ret = file->ReadBuffer(buffer,len);
         if (fc) fc->Enable();
         pf->AddNoCacheBytesRead(len);
         pf->AddNoCacheReadCalls(1);
         if (ret) {
            return 1;
         }
      }
      gPerfStats = temp;
      // fOffset might have been changed via TFileCacheRead::ReadBuffer(), reset it
      file->SetOffset(pos + len);
   } else {
      TVirtualPerfStats* temp = gPerfStats;
      if (tree->GetPerfStats() != 0) gPerfStats = tree->GetPerfStats();
      if (file->ReadBuffer(buffer,len)) {
         gPerfStats = temp;
         return 1; //error while reading
      }
      else gPerfStats = temp;
   }

   fBufferRef->SetReadMode();
   fBufferRef->SetBufferOffset(0);
   Streamer(*fBufferRef);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove the first dentries of this basket, moving entries at
/// dentries to the start of the buffer.

void TBasket::MoveEntries(Int_t dentries)
{
   Int_t i;

   if (dentries >= fNevBuf) return;
   Int_t bufbegin;
   Int_t moved;

   Int_t *entryOffset = GetEntryOffset();
   if (entryOffset) {
      bufbegin = entryOffset[dentries];
      moved = bufbegin-GetKeylen();

      // First store the original location in the fDisplacement array
      // and record the new start offset

      if (!fDisplacement) {
         fDisplacement = new Int_t[fNevBufSize];
      }
      for (i = 0; i<(fNevBufSize-dentries); ++i) {
         fDisplacement[i] = entryOffset[i + dentries];
         entryOffset[i] = entryOffset[i + dentries] - moved;
      }
      for (i = fNevBufSize-dentries; i<fNevBufSize; ++i) {
         fDisplacement[i] = 0;
         entryOffset[i] = 0;
      }

   } else {
      // If there is no EntryOffset array, this means
      // that each entry has the same size and that
      // it does not point to other objects (hence there
      // is no need for a displacement array).
      bufbegin = GetKeylen() + dentries*fNevBufSize;
      moved = bufbegin-GetKeylen();
   }
   TBuffer *buf = GetBufferRef();
   char *buffer = buf->Buffer();
   memmove(buffer+GetKeylen(),buffer+bufbegin,buf->Length()-bufbegin);
   buf->SetBufferOffset(buf->Length()-moved);
   fNevBuf -= dentries;
}

#define OLD_CASE_EXPRESSION fObjlen==fNbytes-fKeylen && GetBranch()->GetCompressionLevel()!=0 && file->GetVersion()<=30401
////////////////////////////////////////////////////////////////////////////////
/// By-passing buffer unzipping has been requested and is
/// possible (only 1 entry in this basket).

Int_t TBasket::ReadBasketBuffersUncompressedCase()
{
   fBuffer = fBufferRef->Buffer();

   // Make sure that the buffer is set at the END of the data
   fBufferRef->SetBufferOffset(fNbytes);

   // Indicate that this buffer is weird.
   fBufferRef->SetBit(TBufferFile::kNotDecompressed);

   // Usage of this mode assume the existance of only ONE
   // entry in this basket.
   ResetEntryOffset();
   delete [] fDisplacement; fDisplacement = 0;

   fBranch->GetTree()->IncrementTotalBuffers(fBufferSize);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// We always create the TBuffer for the basket but it hold the buffer from the cache.

Int_t TBasket::ReadBasketBuffersUnzip(char* buffer, Int_t size, Bool_t mustFree, TFile* file)
{
   if (fBufferRef) {
      fBufferRef->SetBuffer(buffer, size, mustFree);
      fBufferRef->SetReadMode();
      fBufferRef->Reset();
   } else {
      fBufferRef = new TBufferFile(TBuffer::kRead, size, buffer, mustFree);
   }
   fBufferRef->SetParent(file);

   Streamer(*fBufferRef);

   if (IsZombie()) {
      return -1;
   }

   Bool_t oldCase = OLD_CASE_EXPRESSION;

   if ((fObjlen > fNbytes-fKeylen || oldCase) && TestBit(TBufferFile::kNotDecompressed) && (fNevBuf==1)) {
      return TBasket::ReadBasketBuffersUncompressedCase();
   }

   fBuffer = fBufferRef->Buffer();
   return fObjlen+fKeylen;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize a buffer for reading if it is not already initialized

static inline TBuffer* R__InitializeReadBasketBuffer(TBuffer* bufferRef, Int_t len, TFile* file)
{
   TBuffer* result;
   if (R__likely(bufferRef)) {
      bufferRef->SetReadMode();
      Int_t curBufferSize = bufferRef->BufferSize();
      if (curBufferSize < len) {
         // Experience shows that giving 5% "wiggle-room" decreases churn.
         bufferRef->Expand(Int_t(len*1.05));
      }
      bufferRef->Reset();
      result = bufferRef;
   } else {
      result = new TBufferFile(TBuffer::kRead, len);
   }
   result->SetParent(file);
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the compressed buffer; either from the TTree or create a local one.

void inline TBasket::InitializeCompressedBuffer(Int_t len, TFile* file)
{
   Bool_t compressedBufferExists = fCompressedBufferRef != NULL;
   fCompressedBufferRef = R__InitializeReadBasketBuffer(fCompressedBufferRef, len, file);
   if (R__unlikely(!compressedBufferExists)) {
      fOwnsCompressedBuffer = kTRUE;
   }
}

void TBasket::ResetEntryOffset()
{
   if (fEntryOffset != reinterpret_cast<Int_t *>(-1)) {
      delete[] fEntryOffset;
   }
   fEntryOffset = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Read basket buffers in memory and cleanup.
///
/// Read a basket buffer. Check if buffers of previous ReadBasket
/// should not be dropped. Remember, we keep buffers in memory up to
/// fMaxVirtualSize.
/// The function returns 0 in case of success, 1 in case of error
/// This function was modified with the addition of the parallel
/// unzipping, it will try to get the unzipped file from the cache
/// receiving only a pointer to that buffer (so we shall not
/// delete that pointer), although we get a new buffer in case
/// it's not found in the cache.
/// There is a lot of code duplication but it was necesary to assure
/// the expected behavior when there is no cache.

Int_t TBasket::ReadBasketBuffers(Long64_t pos, Int_t len, TFile *file)
{
   if(!fBranch->GetDirectory()) {
      return -1;
   }

   Bool_t oldCase;
   char *rawUncompressedBuffer, *rawCompressedBuffer;
   Int_t uncompressedBufferLen;

   // See if the cache has already unzipped the buffer for us.
   TFileCacheRead *pf = nullptr;
   {
      R__LOCKGUARD_IMT(gROOTMutex); // Lock for parallel TTree I/O
      pf = fBranch->GetTree()->GetReadCache(file);
   }
   if (pf) {
      Int_t res = -1;
      Bool_t free = kTRUE;
      char *buffer = nullptr;
      res = pf->GetUnzipBuffer(&buffer, pos, len, &free);
      if (R__unlikely(res >= 0)) {
         len = ReadBasketBuffersUnzip(buffer, res, free, file);
         // Note that in the kNotDecompressed case, the above function will return 0;
         // In such a case, we should stop processing
         if (len <= 0) return -len;
         goto AfterBuffer;
      }
   }

   // Determine which buffer to use, so that we can avoid a memcpy in case of
   // the basket was not compressed.
   TBuffer* readBufferRef;
   if (R__unlikely(fBranch->GetCompressionLevel()==0)) {
      // Initialize the buffer to hold the uncompressed data.
      fBufferRef = R__InitializeReadBasketBuffer(fBufferRef, len, file);
      readBufferRef = fBufferRef;
   } else {
      // Initialize the buffer to hold the compressed data.
      fCompressedBufferRef = R__InitializeReadBasketBuffer(fCompressedBufferRef, len, file);
      readBufferRef = fCompressedBufferRef;
   }

   // fBufferSize is likely to be change in the Streamer call (below)
   // and we will re-add the new size later on.
   fBranch->GetTree()->IncrementTotalBuffers(-fBufferSize);

   if (!readBufferRef) {
      Error("ReadBasketBuffers", "Unable to allocate buffer.");
      return 1;
   }

   if (pf) {
      TVirtualPerfStats* temp = gPerfStats;
      if (fBranch->GetTree()->GetPerfStats() != 0) gPerfStats = fBranch->GetTree()->GetPerfStats();
      Int_t st = 0;
      {
         R__LOCKGUARD_IMT(gROOTMutex); // Lock for parallel TTree I/O
         st = pf->ReadBuffer(readBufferRef->Buffer(),pos,len);
      }
      if (st < 0) {
         return 1;
      } else if (st == 0) {
         // Read directly from file, not from the cache
         // If we are using a TTreeCache, disable reading from the default cache
         // temporarily, to force reading directly from file
         R__LOCKGUARD_IMT(gROOTMutex);  // Lock for parallel TTree I/O
         TTreeCache *fc = dynamic_cast<TTreeCache*>(file->GetCacheRead());
         if (fc) fc->Disable();
         Int_t ret = file->ReadBuffer(readBufferRef->Buffer(),pos,len);
         if (fc) fc->Enable();
         pf->AddNoCacheBytesRead(len);
         pf->AddNoCacheReadCalls(1);
         if (ret) {
            return 1;
         }
      }
      gPerfStats = temp;
   } else {
      // Read from the file and unstream the header information.
      TVirtualPerfStats* temp = gPerfStats;
      if (fBranch->GetTree()->GetPerfStats() != 0) gPerfStats = fBranch->GetTree()->GetPerfStats();
      R__LOCKGUARD_IMT(gROOTMutex);  // Lock for parallel TTree I/O
      if (file->ReadBuffer(readBufferRef->Buffer(),pos,len)) {
         gPerfStats = temp;
         return 1;
      }
      else gPerfStats = temp;
   }
   Streamer(*readBufferRef);
   if (IsZombie()) {
      return 1;
   }

   rawCompressedBuffer = readBufferRef->Buffer();

   // Are we done?
   if (R__unlikely(readBufferRef == fBufferRef)) // We expect most basket to be compressed.
   {
      if (R__likely(fObjlen+fKeylen == fNbytes)) {
         // The basket was really not compressed as expected.
         goto AfterBuffer;
      } else {
         // Well, somehow the buffer was compressed anyway, we have the compressed data in the uncompressed buffer
         // Make sure the compressed buffer is initialized, and memcpy.
         InitializeCompressedBuffer(len, file);
         if (!fCompressedBufferRef) {
            Error("ReadBasketBuffers", "Unable to allocate buffer.");
            return 1;
         }
         fBufferRef->Reset();
         rawCompressedBuffer = fCompressedBufferRef->Buffer();
         memcpy(rawCompressedBuffer, fBufferRef->Buffer(), len);
      }
   }

   // Initialize buffer to hold the uncompressed data
   // Note that in previous versions we didn't allocate buffers until we verified
   // the zip headers; this is no longer beforehand as the buffer lifetime is scoped
   // to the TBranch.
   uncompressedBufferLen = len > fObjlen+fKeylen ? len : fObjlen+fKeylen;
   fBufferRef = R__InitializeReadBasketBuffer(fBufferRef, uncompressedBufferLen, file);
   rawUncompressedBuffer = fBufferRef->Buffer();
   fBuffer = rawUncompressedBuffer;

   oldCase = OLD_CASE_EXPRESSION;
   // Case where ROOT thinks the buffer is compressed.  Copy over the key and uncompress the object
   if (fObjlen > fNbytes-fKeylen || oldCase) {
      if (R__unlikely(TestBit(TBufferFile::kNotDecompressed) && (fNevBuf==1))) {
         return ReadBasketBuffersUncompressedCase();
      }

      // Optional monitor for zip time profiling.
      Double_t start = 0;
      if (R__unlikely(gPerfStats)) {
         start = TTimeStamp();
      }

      memcpy(rawUncompressedBuffer, rawCompressedBuffer, fKeylen);
      char *rawUncompressedObjectBuffer = rawUncompressedBuffer+fKeylen;
      UChar_t *rawCompressedObjectBuffer = (UChar_t*)rawCompressedBuffer+fKeylen;
      Int_t nin, nbuf;
      Int_t nout = 0, noutot = 0, nintot = 0;

      // Unzip all the compressed objects in the compressed object buffer.
      while (1) {
         // Check the header for errors.
         if (R__unlikely(R__unzip_header(&nin, rawCompressedObjectBuffer, &nbuf) != 0)) {
            Error("ReadBasketBuffers", "Inconsistency found in header (nin=%d, nbuf=%d)", nin, nbuf);
            break;
         }
         if (R__unlikely(oldCase && (nin > fObjlen || nbuf > fObjlen))) {
            //buffer was very likely not compressed in an old version
            memcpy(rawUncompressedBuffer+fKeylen, rawCompressedObjectBuffer+fKeylen, fObjlen);
            goto AfterBuffer;
         }

         R__unzip(&nin, rawCompressedObjectBuffer, &nbuf, (unsigned char*) rawUncompressedObjectBuffer, &nout);
         if (!nout) break;
         noutot += nout;
         nintot += nin;
         if (noutot >= fObjlen) break;
         rawCompressedObjectBuffer += nin;
         rawUncompressedObjectBuffer += nout;
      }

      // Make sure the uncompressed numbers are consistent with header.
      if (R__unlikely(noutot != fObjlen)) {
         Error("ReadBasketBuffers", "fNbytes = %d, fKeylen = %d, fObjlen = %d, noutot = %d, nout=%d, nin=%d, nbuf=%d", fNbytes,fKeylen,fObjlen, noutot,nout,nin,nbuf);
         fBranch->GetTree()->IncrementTotalBuffers(fBufferSize);
         return 1;
      }
      len = fObjlen+fKeylen;
      TVirtualPerfStats* temp = gPerfStats;
      if (fBranch->GetTree()->GetPerfStats() != 0) gPerfStats = fBranch->GetTree()->GetPerfStats();
      if (R__unlikely(gPerfStats)) {
         gPerfStats->UnzipEvent(fBranch->GetTree(),pos,start,nintot,fObjlen);
      }
      gPerfStats = temp;
   } else {
      // Nothing is compressed - copy over wholesale.
      memcpy(rawUncompressedBuffer, rawCompressedBuffer, len);
   }

AfterBuffer:

   fBranch->GetTree()->IncrementTotalBuffers(fBufferSize);

   // Read offsets table if needed.
   // If there's no EntryOffsetLen in the branch -- or the fEntryOffset is marked to be calculated-on-demand --
   // then we skip reading out.
   if (!fBranch->GetEntryOffsetLen() || (fEntryOffset == reinterpret_cast<Int_t *>(-1))) {
      return 0;
   }
   // At this point, we're required to read out an offset array.
   ResetEntryOffset(); // TODO: every basket, we reset the offset array.  Is this necessary?
                       // Could we instead switch to std::vector?
   fBufferRef->SetBufferOffset(fLast);
   fBufferRef->ReadArray(fEntryOffset);
   if (R__unlikely(!fEntryOffset)) {
      fEntryOffset = new Int_t[fNevBuf+1];
      fEntryOffset[0] = fKeylen;
      Warning("ReadBasketBuffers","basket:%s has fNevBuf=%d but fEntryOffset=0, pos=%lld, len=%d, fNbytes=%d, fObjlen=%d, trying to repair",GetName(),fNevBuf,pos,len,fNbytes,fObjlen);
      return 0;
   }
   if (fIOBits & static_cast<UChar_t>(TBasket::EIOBits::kGenerateOffsetMap)) {
      // In this case, we cannot regenerate the offset array at runtime -- but we wrote out an array of
      // sizes instead of offsets (as sizes compress much better).
      fEntryOffset[0] = fKeylen;
      for (Int_t idx = 1; idx < fNevBuf + 1; idx++) {
         fEntryOffset[idx] += fEntryOffset[idx - 1];
      }
   }
   fReadEntryOffset = kTRUE;
   // Read the array of diplacement if any.
   delete [] fDisplacement;
   fDisplacement = 0;
   if (fBufferRef->Length() != len) {
      // There is more data in the buffer!  It is the displacement
      // array.  If len is less than TBuffer::kMinimalSize the actual
      // size of the buffer is too large, so we can not use the
      // fBufferRef->BufferSize()
      fBufferRef->ReadArray(fDisplacement);
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read basket buffers in memory and cleanup
///
/// Read first bytes of a logical record starting at position pos
/// return record length (first 4 bytes of record).

Int_t TBasket::ReadBasketBytes(Long64_t pos, TFile *file)
{
   const Int_t len = 128;
   char buffer[len];
   Int_t keylen;
   file->GetRecordHeader(buffer, pos,len, fNbytes, fObjlen, keylen);
   fKeylen = keylen;
   return fNbytes;
}

////////////////////////////////////////////////////////////////////////////////
/// Disown all references to the internal buffer - some other object likely now
/// owns it.
///
/// This TBasket is now useless and invalid until it is told to adopt a buffer.
void TBasket::DisownBuffer()
{
   fBufferRef = NULL;
}


////////////////////////////////////////////////////////////////////////////////
/// Adopt a buffer from an external entity
void TBasket::AdoptBuffer(TBuffer *user_buffer)
{
   delete fBufferRef;
   fBufferRef = user_buffer;
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the read basket TBuffer memory allocation if needed.
///
/// This allows to reduce the number of memory allocation while avoiding to
/// always use the maximum size.

void TBasket::ReadResetBuffer(Int_t basketnumber)
{
   // By default, we don't reallocate.
   fResetAllocation = false;
#ifdef R__TRACK_BASKET_ALLOC_TIME
   fResetAllocationTime = 0;
#endif

   // Downsize the buffer if needed.

   const auto maxbaskets = fBranch->GetMaxBaskets();
   if (!fBufferRef || basketnumber >= maxbaskets)
      return;

   Int_t curSize = fBufferRef->BufferSize();

   const Float_t target_mem_ratio = fBranch->GetTree()->GetTargetMemoryRatio();
   const auto basketbytes = fBranch->GetBasketBytes();

   Int_t max_size = basketbytes[basketnumber];
   for(Int_t b = basketnumber + 1; (b < maxbaskets) && (b < (basketnumber+10)); ++b) {
      max_size = std::max(max_size, basketbytes[b]);
   }

   Float_t cx = 1;
   if (fBranch->GetZipBytes())
      cx = (Float_t)fBranch->GetTotBytes()/fBranch->GetZipBytes();

   Int_t target_size = static_cast<Int_t>(cx * target_mem_ratio * Float_t(max_size));

   if (target_size && (curSize > target_size)) {
      /// Only reduce the size if significant enough?
      Int_t newSize = max_size + 512 - max_size % 512; // Wiggle room and alignment, as above.
      // We only bother with a resize if it saves 8KB (two normal memory pages).
      if ((newSize <= curSize - 8 * 1024) &&
          (static_cast<Float_t>(curSize) / static_cast<Float_t>(newSize) > target_mem_ratio))
      {
         if (gDebug > 0) {
            Info("ReadResetBuffer",
                 "Resizing %d to %d bytes (was %d); next 10 sizes are [%d, %d, %d, %d, %d, %d, %d, %d, %d, %d]. cx=%f ratio=%f max_size = %d ",
                 basketnumber, newSize, curSize,
                 basketbytes[basketnumber],
                 (basketnumber + 1) < maxbaskets ? basketbytes[basketnumber + 1] : 0,
                 (basketnumber + 2) < maxbaskets ? basketbytes[basketnumber + 2] : 0,
                 (basketnumber + 3) < maxbaskets ? basketbytes[basketnumber + 3] : 0,
                 (basketnumber + 4) < maxbaskets ? basketbytes[basketnumber + 4] : 0,
                 (basketnumber + 5) < maxbaskets ? basketbytes[basketnumber + 5] : 0,
                 (basketnumber + 6) < maxbaskets ? basketbytes[basketnumber + 6] : 0,
                 (basketnumber + 7) < maxbaskets ? basketbytes[basketnumber + 7] : 0,
                 (basketnumber + 8) < maxbaskets ? basketbytes[basketnumber + 8] : 0,
                 (basketnumber + 9) < maxbaskets ? basketbytes[basketnumber + 9] : 0,
                 cx, target_mem_ratio, max_size);
         }
         fResetAllocation = true;
#ifdef R__TRACK_BASKET_ALLOC_TIME
         std::chrono::time_point<std::chrono::system_clock> start, end;
         start = std::chrono::high_resolution_clock::now();
#endif
         fBufferRef->Expand(newSize, kFALSE); // Expand without copying the existing data.
#ifdef R__TRACK_BASKET_ALLOC_TIME
         end = std::chrono::high_resolution_clock::now();
         auto us = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
         fResetAllocationTime = us.count();
#endif
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the write basket to the starting state. i.e. as it was after calling
/// the constructor (and potentially attaching a TBuffer.)
/// Reduce memory used by fEntryOffset and the TBuffer if needed ..

void TBasket::WriteReset()
{
   // By default, we don't reallocate.
   fResetAllocation = false;
#ifdef R__TRACK_BASKET_ALLOC_TIME
   fResetAllocationTime = 0;
#endif

   // Name, Title, fClassName, fBranch
   // stay the same.

   // Downsize the buffer if needed.
   // See if our current buffer size is significantly larger (>2x) than the historical average.
   // If so, try decreasing it at this flush boundary to closer to the size from OptimizeBaskets
   // (or this historical average).
   Int_t curSize = fBufferRef->BufferSize();
   // fBufferLen at this point is already reset, so use indirect measurements
   Int_t curLen = (GetObjlen() + GetKeylen());
   Long_t newSize = -1;
   if (curSize > 2*curLen)
   {
      Long_t curBsize = fBranch->GetBasketSize();
      if (curSize > 2*curBsize ) {
         Long_t avgSize = (Long_t)(fBranch->GetTotBytes() / (1+fBranch->GetWriteBasket())); // Average number of bytes per basket so far
         if (curSize > 2*avgSize) {
            newSize = curBsize;
            if (curLen > newSize) {
               newSize = curLen;
            }
            if (avgSize > newSize) {
               newSize = avgSize;
            }
            newSize = newSize + 512 - newSize%512;  // Wiggle room and alignment (512 is same as in OptimizeBaskets)
         }
      }
   }
   // If fBufferRef grew since we last saw it, shrink it to "target memory ratio" of the occupied size
   // This discourages us from having poorly-occupied buffers on branches with little variability.
   //
   // Does not help protect against a burst in event sizes, but does help in the cases where the basket
   // size jumps from 4MB to 8MB while filling the basket, but we only end up utilizing 4.1MB.
   //
   // The above code block is meant to protect against extremely large events.

   Float_t target_mem_ratio = fBranch->GetTree()->GetTargetMemoryRatio();
   Int_t max_size = TMath::Max(fLastWriteBufferSize[0], std::max(fLastWriteBufferSize[1], fLastWriteBufferSize[2]));
   Int_t target_size = static_cast<Int_t>(target_mem_ratio * Float_t(max_size));
   if (max_size && (curSize > target_size) && (newSize == -1)) {
      newSize = target_size;
      newSize = newSize + 512 - newSize % 512; // Wiggle room and alignment, as above.
      // We only bother with a resize if it saves 8KB (two normal memory pages).
      if ((newSize > curSize - 8 * 1024) ||
          (static_cast<Float_t>(curSize) / static_cast<Float_t>(newSize) < target_mem_ratio)) {
         newSize = -1;
      } else if (gDebug > 0) {
         Info("Reset", "Resizing to %ld bytes (was %d); last three sizes were [%d, %d, %d].", newSize, curSize,
              fLastWriteBufferSize[0], fLastWriteBufferSize[1], fLastWriteBufferSize[2]);
      }
   }

   if (newSize != -1) {
      fResetAllocation = true;
#ifdef R__TRACK_BASKET_ALLOC_TIME
      std::chrono::time_point<std::chrono::system_clock> start, end;
      start = std::chrono::high_resolution_clock::now();
#endif
      fBufferRef->Expand(newSize,kFALSE);     // Expand without copying the existing data.
#ifdef R__TRACK_BASKET_ALLOC_TIME
      end = std::chrono::high_resolution_clock::now();
      auto us = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
      fResetAllocationTime = us.count();
#endif
   }

   // Record the actual occupied size of the buffer.
   fLastWriteBufferSize[fNextBufferSizeRecord] = curLen;
   fNextBufferSizeRecord = (fNextBufferSizeRecord + 1) % 3;

   TKey::Reset();

   Int_t newNevBufSize = fBranch->GetEntryOffsetLen();
   if (newNevBufSize==0) {
      ResetEntryOffset();
   } else if ((newNevBufSize != fNevBufSize) || (fEntryOffset <= reinterpret_cast<Int_t *>(-1))) {
      ResetEntryOffset();
      fEntryOffset = new Int_t[newNevBufSize];
   }
   fNevBufSize = newNevBufSize;

   fNevBuf      = 0;
   Int_t *storeEntryOffset = fEntryOffset;
   fEntryOffset = 0;
   Int_t *storeDisplacement = fDisplacement;
   fDisplacement= 0;
   fBuffer      = 0;

   fBufferRef->Reset();
   fBufferRef->SetWriteMode();

   fHeaderOnly  = kTRUE;
   fLast        = 0;  //Must initialize before calling Streamer()

   Streamer(*fBufferRef);

   fKeylen      = fBufferRef->Length();
   fObjlen      = fBufferSize - fKeylen;
   fLast        = fKeylen;
   fBuffer      = 0;
   fHeaderOnly  = kFALSE;
   fDisplacement= storeDisplacement;
   fEntryOffset = storeEntryOffset;
   if (fNevBufSize) {
      for (Int_t i=0;i<fNevBufSize;i++) fEntryOffset[i] = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set read mode of basket.

void TBasket::SetReadMode()
{
   fLast = fBufferRef->Length();
   fBufferRef->SetReadMode();
}

////////////////////////////////////////////////////////////////////////////////
/// Set write mode of basket.

void TBasket::SetWriteMode()
{
   fBufferRef->SetWriteMode();
   fBufferRef->SetBufferOffset(fLast);
}

////////////////////////////////////////////////////////////////////////////////
/// Stream a class object.

void TBasket::Streamer(TBuffer &b)
{
   // As in TBranch::GetBasket, this is used as a half-hearted measure to suppress
   // the error reporting when many failures occur.
   static std::atomic<Int_t> nerrors(0);

   char flag;
   if (b.IsReading()) {
      TKey::Streamer(b); //this must be first
      Version_t v = b.ReadVersion();
      b >> fBufferSize;
      // NOTE: we now use the upper-bit of the fNevBufSize to see if we have serialized any of the
      // optional IOBits.  If that bit is set, we immediately read out the IOBits; to replace this
      // (minimal) safeguard against corruption, we will set aside the upper-bit of fIOBits to do
      // the same thing (the fact this bit is reserved is tested in the unit tests).  If there is
      // someday a need for more than 7 IOBits, we'll widen the field using the same trick.
      //
      // We like to keep this safeguard because we immediately will allocate a buffer based on
      // the value of fNevBufSize -- and would like to avoid wildly inappropriate allocations.
      b >> fNevBufSize;
      if (fNevBufSize < 0) {
         fNevBufSize = -fNevBufSize;
         b >> fIOBits;
         if (!fIOBits || (fIOBits & (1 << 7))) {
            Error("TBasket::Streamer",
                  "The value of fNevBufSize (%d) or fIOBits (%d) is incorrect ; setting the buffer to a zombie.",
                  -fNevBufSize, fIOBits);
            MakeZombie();
            fNevBufSize = 0;
         } else if (fIOBits && (fIOBits & ~static_cast<Int_t>(EIOBits::kSupported))) {
            nerrors++;
            if (nerrors < 10) {
               Error("Streamer", "The value of fIOBits (%s) contains unknown flags (supported flags "
                                 "are %s), indicating this was written with a newer version of ROOT "
                                 "utilizing critical IO features this version of ROOT does not support."
                                 "  Refusing to deserialize.",
                     std::bitset<32>(static_cast<Int_t>(fIOBits)).to_string().c_str(),
                     std::bitset<32>(static_cast<Int_t>(EIOBits::kSupported)).to_string().c_str());
            } else if (nerrors == 10) {
               Error("Streamer", "Maximum number of errors has been reported; disabling further messages"
                                 "from this location until the process exits.");
            }
            fNevBufSize = 0;
            MakeZombie();
         }
      }
      b >> fNevBuf;
      b >> fLast;
      b >> flag;
      if (fLast > fBufferSize) fBufferSize = fLast;
      Bool_t mustGenerateOffsets = false;
      if (flag >= 80) {
         mustGenerateOffsets = true;
         flag -= 80;
      }
      if (!mustGenerateOffsets && flag && (flag % 10 != 2)) {
         ResetEntryOffset();
         fEntryOffset = new Int_t[fNevBufSize];
         if (fNevBuf) b.ReadArray(fEntryOffset);
         if (20<flag && flag<40) {
            for(int i=0; i<fNevBuf; i++){
               fEntryOffset[i] &= ~kDisplacementMask;
            }
         }
         if (flag>40) {
            fDisplacement = new Int_t[fNevBufSize];
            b.ReadArray(fDisplacement);
         }
      } else if (mustGenerateOffsets) {
         // We currently believe that in all cases when offsets can be generated, then the
         // displacement array must be zero.
         assert(flag <= 40);
         fEntryOffset = reinterpret_cast<Int_t *>(-1);
      }
      if (flag == 1 || flag > 10) {
         fBufferRef = new TBufferFile(TBuffer::kRead,fBufferSize);
         fBufferRef->SetParent(b.GetParent());
         char *buf  = fBufferRef->Buffer();
         if (v > 1) b.ReadFastArray(buf,fLast);
         else       b.ReadArray(buf);
         fBufferRef->SetBufferOffset(fLast);
         // This is now done in the TBranch streamer since fBranch might not
         // yet be set correctly.
         //   fBranch->GetTree()->IncrementTotalBuffers(fBufferSize);
      }
   } else {

      TKey::Streamer(b);   //this must be first
      b.WriteVersion(TBasket::IsA());
      if (fBufferRef) {
         Int_t curLast = fBufferRef->Length();
         if (!fHeaderOnly && !fSeekKey && curLast > fLast) fLast = curLast;
      }
      if (fLast > fBufferSize) fBufferSize = fLast;

      b << fBufferSize;
      if (fIOBits) {
         b << -fNevBufSize;
         b << fIOBits;
      } else {
         b << fNevBufSize;
      }
      b << fNevBuf;
      b << fLast;
      Bool_t mustGenerateOffsets = fEntryOffset && fNevBuf &&
                                   (fIOBits & static_cast<UChar_t>(TBasket::EIOBits::kGenerateOffsetMap)) &&
                                   CanGenerateOffsetArray();
      // We currently believe that in all cases when offsets can be generated, then the
      // displacement array must be zero.
      assert(!mustGenerateOffsets || fDisplacement == nullptr);
      if (fHeaderOnly) {
         flag = mustGenerateOffsets ? 80 : 0;
         b << flag;
      } else {
         // On return from this function, we are guaranteed that fEntryOffset
         // is either a valid pointer or nullptr.
         if (fNevBuf) {
            GetEntryOffset();
         }
         flag = 1;
         if (!fNevBuf || !fEntryOffset)
            flag = 2;
         if (fBufferRef)     flag += 10;
         if (fDisplacement)  flag += 40;
         // Test if we can skip writing out the offset map.
         if (mustGenerateOffsets) {
            flag += 80;
         }
         b << flag;

         if (!mustGenerateOffsets && fEntryOffset && fNevBuf) {
            b.WriteArray(fEntryOffset, fNevBuf);
            if (fDisplacement) b.WriteArray(fDisplacement, fNevBuf);
         }
         if (fBufferRef) {
            char *buf  = fBufferRef->Buffer();
            b.WriteFastArray(buf, fLast);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Update basket header and EntryOffset table.

void TBasket::Update(Int_t offset, Int_t skipped)
{
   Int_t *entryOffset = GetEntryOffset();
   if (entryOffset) {
      if (fNevBuf+1 >= fNevBufSize) {
         Int_t newsize = TMath::Max(10,2*fNevBufSize);
         Int_t *newoff = TStorage::ReAllocInt(fEntryOffset, newsize,
                                              fNevBufSize);
         if (fDisplacement) {
            Int_t *newdisp = TStorage::ReAllocInt(fDisplacement, newsize,
                                                  fNevBufSize);
            fDisplacement = newdisp;
         }
         fEntryOffset  = newoff;
         fNevBufSize   = newsize;

         //Update branch only for the first 10 baskets
         if (fBranch->GetWriteBasket() < 10) {
            fBranch->SetEntryOffsetLen(newsize);
         }
      }
      fEntryOffset[fNevBuf] = offset;

      if (skipped!=offset && !fDisplacement){
         fDisplacement = new Int_t[fNevBufSize];
         for (Int_t i = 0; i<fNevBufSize; i++) fDisplacement[i] = fEntryOffset[i];
      }
      if (fDisplacement) {
         fDisplacement[fNevBuf] = skipped;
         fBufferRef->SetBufferDisplacement(skipped);
      }
   }

   fNevBuf++;
}

////////////////////////////////////////////////////////////////////////////////
/// Write buffer of this basket on the current file.
///
/// The function returns the number of bytes committed to the memory.
/// If a write error occurs, the number of bytes returned is -1.
/// If no data are written, the number of bytes returned is 0.

Int_t TBasket::WriteBuffer()
{
   const Int_t kWrite = 1;

   TFile *file = fBranch->GetFile(kWrite);
   if (!file) return 0;
   if (!file->IsWritable()) {
      return -1;
   }
   fMotherDir = file; // fBranch->GetDirectory();

   // This mutex prevents multiple TBasket::WriteBuffer invocations from interacting
   // with the underlying TFile at once - TFile is assumed to *not* be thread-safe.
   //
   // The only parallelism we'd like to exploit (right now!) is the compression
   // step - everything else should be serialized at the TFile level.
#ifdef R__USE_IMT
   std::unique_lock<std::mutex> sentry(file->fWriteMutex);
#endif  // R__USE_IMT

   if (R__unlikely(fBufferRef->TestBit(TBufferFile::kNotDecompressed))) {
      // Read the basket information that was saved inside the buffer.
      Bool_t writing = fBufferRef->IsWriting();
      fBufferRef->SetReadMode();
      fBufferRef->SetBufferOffset(0);

      Streamer(*fBufferRef);
      if (writing) fBufferRef->SetWriteMode();
      Int_t nout = fNbytes - fKeylen;

      fBuffer = fBufferRef->Buffer();

      Create(nout,file);
      fBufferRef->SetBufferOffset(0);
      fHeaderOnly = kTRUE;

      Streamer(*fBufferRef);         //write key itself again
      int nBytes = WriteFileKeepBuffer();
      fHeaderOnly = kFALSE;
      return nBytes>0 ? fKeylen+nout : -1;
   }

   // Transfer fEntryOffset table at the end of fBuffer.
   fLast = fBufferRef->Length();
   Int_t *entryOffset = GetEntryOffset();
   if (entryOffset) {
      Bool_t hasOffsetBit = fIOBits & static_cast<UChar_t>(TBasket::EIOBits::kGenerateOffsetMap);
      if (!CanGenerateOffsetArray()) {
         // If we have set the offset map flag, but cannot dynamically generate the map, then
         // we should at least convert the offset array to a size array.  Note that we always
         // write out (fNevBuf+1) entries to match the original case.
         if (hasOffsetBit) {
            for (Int_t idx = fNevBuf; idx > 0; idx--) {
               entryOffset[idx] -= entryOffset[idx - 1];
            }
            entryOffset[0] = 0;
         }
         fBufferRef->WriteArray(entryOffset, fNevBuf + 1);
         // Convert back to offset format: keeping both sizes and offsets in-memory were considered,
         // but it seems better to use CPU than memory.
         if (hasOffsetBit) {
            entryOffset[0] = fKeylen;
            for (Int_t idx = 1; idx < fNevBuf + 1; idx++) {
               entryOffset[idx] += entryOffset[idx - 1];
            }
         }
      } else if (!hasOffsetBit) { // In this case, write out as normal
         fBufferRef->WriteArray(entryOffset, fNevBuf + 1);
      }
      if (fDisplacement) {
         fBufferRef->WriteArray(fDisplacement, fNevBuf + 1);
         delete[] fDisplacement;
         fDisplacement = 0;
      }
   }

   Int_t lbuf, nout, noutot, bufmax, nzip;
   lbuf       = fBufferRef->Length();
   fObjlen    = lbuf - fKeylen;

   fHeaderOnly = kTRUE;
   fCycle = fBranch->GetWriteBasket();
   Int_t cxlevel = fBranch->GetCompressionLevel();
   if (cxlevel == ROOT::RCompressionSetting::ELevel::kInherit)
      cxlevel = file->GetCompressionLevel();
   ROOT::RCompressionSetting::EAlgorithm::EValues cxAlgorithm = static_cast<ROOT::RCompressionSetting::EAlgorithm::EValues>(fBranch->GetCompressionAlgorithm());
   if (cxAlgorithm == ROOT::RCompressionSetting::EAlgorithm::kInherit)
      cxAlgorithm = static_cast<ROOT::RCompressionSetting::EAlgorithm::EValues>(file->GetCompressionAlgorithm());
   if (cxlevel > 0) {
      Int_t nbuffers = 1 + (fObjlen - 1) / kMAXZIPBUF;
      Int_t buflen = fKeylen + fObjlen + 9 * nbuffers + 28; //add 28 bytes in case object is placed in a deleted gap
      InitializeCompressedBuffer(buflen, file);
      if (!fCompressedBufferRef) {
         Warning("WriteBuffer", "Unable to allocate the compressed buffer");
         return -1;
      }
      fCompressedBufferRef->SetWriteMode();
      fBuffer = fCompressedBufferRef->Buffer();
      char *objbuf = fBufferRef->Buffer() + fKeylen;
      char *bufcur = &fBuffer[fKeylen];
      noutot = 0;
      nzip   = 0;
      for (Int_t i = 0; i < nbuffers; ++i) {
         if (i == nbuffers - 1) bufmax = fObjlen - nzip;
         else bufmax = kMAXZIPBUF;
         // Compress the buffer.  Note that we allow multiple TBasket compressions to occur at once
         // for a given TFile: that's because the compression buffer when we use IMT is no longer
         // shared amongst several threads.
#ifdef R__USE_IMT
         sentry.unlock();
#endif  // R__USE_IMT
         // NOTE this is declared with C linkage, so it shouldn't except.  Also, when
         // USE_IMT is defined, we are guaranteed that the compression buffer is unique per-branch.
         // (see fCompressedBufferRef in constructor).
         R__zipMultipleAlgorithm(cxlevel, &bufmax, objbuf, &bufmax, bufcur, &nout, cxAlgorithm);
#ifdef R__USE_IMT
         sentry.lock();
#endif  // R__USE_IMT

         // test if buffer has really been compressed. In case of small buffers
         // when the buffer contains random data, it may happen that the compressed
         // buffer is larger than the input. In this case, we write the original uncompressed buffer
         if (nout == 0 || nout >= fObjlen) {
            nout = fObjlen;
            // We used to delete fBuffer here, we no longer want to since
            // the buffer (held by fCompressedBufferRef) might be re-used later.
            fBuffer = fBufferRef->Buffer();
            Create(fObjlen,file);
            fBufferRef->SetBufferOffset(0);

            Streamer(*fBufferRef);         //write key itself again
            if ((nout+fKeylen)>buflen) {
               Warning("WriteBuffer","Possible memory corruption due to compression algorithm, wrote %d bytes past the end of a block of %d bytes. fNbytes=%d, fObjLen=%d, fKeylen=%d",
                  (nout+fKeylen-buflen),buflen,fNbytes,fObjlen,fKeylen);
            }
            goto WriteFile;
         }
         bufcur += nout;
         noutot += nout;
         objbuf += kMAXZIPBUF;
         nzip   += kMAXZIPBUF;
      }
      nout = noutot;
      Create(noutot,file);
      fBufferRef->SetBufferOffset(0);

      Streamer(*fBufferRef);         //write key itself again
      memcpy(fBuffer,fBufferRef->Buffer(),fKeylen);
   } else {
      fBuffer = fBufferRef->Buffer();
      Create(fObjlen,file);
      fBufferRef->SetBufferOffset(0);

      Streamer(*fBufferRef);         //write key itself again
      nout = fObjlen;
   }

WriteFile:
   Int_t nBytes = WriteFileKeepBuffer();
   fHeaderOnly = kFALSE;
   return nBytes>0 ? fKeylen+nout : -1;
}
