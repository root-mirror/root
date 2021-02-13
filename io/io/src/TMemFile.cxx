// @(#)root/io:$Id$
// Author: Philippe Canal, May 2011

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**
\class TMemFile TMemFile.cxx
\ingroup IO

A TMemFile is like a normal TFile except that it reads and writes
only from memory.
*/

#include "TBufferFile.h"
#include "TMemFile.h"
#include "TError.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TArrayC.h"
#include "TKey.h"
#include "TClass.h"
#include "TVirtualMutex.h"
#include <errno.h>
#include <stdio.h>
#include <sys/stat.h>

// The following snippet is used for developer-level debugging
#define TMemFile_TRACE
#ifndef TMemFile_TRACE
#define TRACE(x) \
  Debug("TMemFile", "%s", x);
#else
#define TRACE(x);
#endif

ClassImp(TMemFile);

////////////////////////////////////////////////////////////////////////////////
/// Constructor allocating the memory buffer.
///
/// \param size: size of the buffer to be allocated. A value of -1 means that
///              no allocation should happen, leaving fBuffer and fSize at 0.
///
/// \param previous: previous TMemBlock, used to set up the linked list.

TMemFile::TMemBlock::TMemBlock(Long64_t size, TMemBlock *previous) : fPrevious(previous)
{
   // size will be -1 when copying an existing buffer into fBuffer.
   if (size != -1) {
      fBuffer = new UChar_t[size];
      fSize = size;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor not allocating the memory buffer, for external ownership.

TMemFile::TMemBlock::TMemBlock(UChar_t *data, Long64_t size) : fBuffer(data), fSize(size)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Usual destructors.  Delete the block memory.

TMemFile::TMemBlock::~TMemBlock()
{
   delete fNext;
   delete [] fBuffer;
}

////////////////////////////////////////////////////////////////////////////////

void TMemFile::TMemBlock::CreateNext(Long64_t size)
{
   R__ASSERT(fNext == nullptr);
   fNext = new TMemBlock(size,this);
}

////////////////////////////////////////////////////////////////////////////////
/// Parse option strings and set fOption.
TMemFile::EMode TMemFile::ParseOption(Option_t *option)
{
   fOption = option;
   fOption.ToUpper();
   if (fOption == "NEW")  fOption = "CREATE";

   EMode mode = EMode::kRead;
   if (fOption == "CREATE")
      mode = EMode::kCreate;
   else if (fOption == "RECREATE")
      mode = EMode::kRecreate;
   else if (fOption == "UPDATE")
      mode = EMode::kUpdate;
   else {
      fOption = "READ";
   }

   return mode;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor to create a TMemFile re-using external C-Style storage.

TMemFile::TMemFile(const char *path, const ZeroCopyView_t &datarange)
   : TFile(path, "WEB", "read-only TMemFile", 0 /*compress*/),
     fBlockList(reinterpret_cast<UChar_t *>(const_cast<char *>(datarange.fStart)), datarange.fSize),
     fSize(datarange.fSize), fBlockSeek(&(fBlockList))
{
   fD = 0;
   fOption = "READ";
   fWritable = kFALSE;

   // This is read-only, so become a zombie if created with an empty buffer
   if (!fBlockList.fBuffer) {
      MakeZombie();
      gDirectory = gROOT;
      return;
   }

   Init(/* create */ false);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor to create a TMemFile re-using external storage.

TMemFile::TMemFile(const char *path, ExternalDataPtr_t data)
   : TMemFile(path, ZeroCopyView_t(data->data(), data->size()))
{
   fExternalData = data;
}

////////////////////////////////////////////////////////////////////////////////////
/// Constructor to create a read-only TMemFile using an std::unique_ptr<TBufferFile>

TMemFile::TMemFile(const char *name, std::unique_ptr<TBufferFile> buffer)
   : TMemFile(name, ZeroCopyView_t(buffer->Buffer(), (size_t)buffer->BufferSize()))
{
   assert(!fD && !fWritable);

   fIsOwnedByROOT = true;

   // Note: We need to release the buffer here to avoid double delete.
   // The memory of a TBufferFile is allocated with new[], so we can let
   // TMemBlock delete it, as its destructor calls "delete [] fBuffer;"
   buffer.release();
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Usual Constructor.
/// The defBlockSize parameter defines the size of the blocks of memory allocated
/// when expanding the underlying TMemFileBuffer. If the value 0 is passed, the
/// default block size, fgDefaultBlockSize, is adopted.
/// See the TFile constructor for details.

TMemFile::TMemFile(const char *path, Option_t *option, const char *ftitle, Int_t compress, Long64_t defBlockSize)
   : TMemFile(path, nullptr, -1, option, ftitle, compress, defBlockSize)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Usual Constructor.  See the TFile constructor for details. Copy data from buffer.

TMemFile::TMemFile(const char *path, char *buffer, Long64_t size, Option_t *option, const char *ftitle, Int_t compress,
                   Long64_t defBlockSize)
   : TFile(path, "WEB", ftitle, compress), fBlockList(size), fIsOwnedByROOT(kTRUE), fSize(size),
     fBlockSeek(&(fBlockList))
{
   fDefaultBlockSize = defBlockSize == 0LL ? fgDefaultBlockSize : defBlockSize;

   EMode optmode = ParseOption(option);

   if (NeedsToWrite(optmode)) {
      Int_t mode = O_RDWR | O_CREAT;
      if (optmode == EMode::kRecreate) mode |= O_TRUNC;

      fD = TMemFile::SysOpen(path, O_RDWR | O_CREAT, 0644);
      if (fD == -1) {
         SysError("TMemFile", "file %s can not be opened", path);
         goto zombie;
      }
      fWritable = kTRUE;

   } else {
      fD = TMemFile::SysOpen(path, O_RDONLY, 0644);
      if (fD == -1) {
         SysError("TMemFile", "file %s can not be opened for reading", path);
         goto zombie;
      }
      fWritable = kFALSE;
   }

   if (buffer)
      SysWriteImpl(fD,buffer,size);

   Init(!NeedsExistingFile(optmode));
   return;

zombie:
   // Error in opening file; make this a zombie
   MakeZombie();
   gDirectory = gROOT;
}

////////////////////////////////////////////////////////////////////////////////
/// Copying the content of the TMemFile into another TMemFile.

TMemFile::TMemFile(const TMemFile &orig)
   : TFile(orig.GetEndpointUrl()->GetUrl(), "WEB", orig.GetTitle(), orig.GetCompressionSettings()),
     fBlockList(orig.GetEND()), fExternalData(orig.fExternalData), fIsOwnedByROOT(orig.fIsOwnedByROOT),
     fSize(orig.GetEND()), fBlockSeek(&(fBlockList))
{
   EMode optmode = ParseOption(orig.fOption);

   fD = orig.fD; // not really used, so it is okay to have the same value.
   fWritable = orig.fWritable;

   if (!IsExternalData()) {
      // We intentionally allocated just one big buffer for this object.
      orig.CopyTo(fBlockList.fBuffer,fSize);
   }

   Init(!NeedsExistingFile(optmode)); // A copy is
}


////////////////////////////////////////////////////////////////////////////////
/// Close and clean-up file.

TMemFile::~TMemFile()
{
   // Need to call close, now as it will need both our virtual table
   // and the content of the list of blocks
   Close();
   if (IsExternalData()) {
      // Do not delete external buffer, we don't own it.
      fBlockList.fBuffer = nullptr;
      // We must not get extra blocks, as writing is disabled for external data!
      R__ASSERT(!fBlockList.fNext && "External block is not the only one!");
   }
   TRACE("destroy")
}

////////////////////////////////////////////////////////////////////////////////
/// Copy the binary representation of the TMemFile into
/// the memory area starting at 'to' and of length at most 'maxsize'
/// returns the number of bytes actually copied.

Long64_t TMemFile::CopyTo(void *to, Long64_t maxsize) const
{
   Long64_t len = GetSize();
   if (len > maxsize) {
      len = maxsize;
   }
   Long64_t storedSysOffset   = fSysOffset;
   Long64_t storedBlockOffset = fBlockOffset;
   TMemBlock *storedBlockSeek = fBlockSeek;

   const_cast<TMemFile*>(this)->SysSeek(fD, 0, SEEK_SET);
   len = const_cast<TMemFile*>(this)->SysReadImpl(fD, to, len);

   const_cast<TMemFile*>(this)->fBlockSeek   = storedBlockSeek;
   const_cast<TMemFile*>(this)->fBlockOffset = storedBlockOffset;
   const_cast<TMemFile*>(this)->fSysOffset   = storedSysOffset;
   return len;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy the binary representation of the TMemFile into
/// the TBuffer tobuf

void TMemFile::CopyTo(TBuffer &tobuf) const
{
   const TMemBlock *current = &fBlockList;
   while(current) {
      tobuf.WriteFastArray(current->fBuffer,current->fSize);
      current = current->fNext;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the current size of the memory file

Long64_t TMemFile::GetSize() const
{
   // We could also attempt to read it from the beginning of the buffer
   return fSize;
}

////////////////////////////////////////////////////////////////////////////////

void TMemFile::Print(Option_t *option /* = "" */) const
{
   Printf("TMemFile: name=%s, title=%s, option=%s", GetName(), GetTitle(), GetOption());
   if (strcmp(option,"blocks")==0) {
      const TMemBlock *current = &fBlockList;
      Int_t counter = 0;
      while(current) {
         Printf("TMemBlock: %d size=%lld addr=%p curr=%p prev=%p next=%p",
                counter,current->fSize,current->fBuffer,
                current,current->fPrevious,current->fNext);
         current = current->fNext;
         ++counter;
      }
   } else {
      GetList()->R__FOR_EACH(TObject,Print)(option);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Wipe all the data from the permanent buffer but keep, the in-memory object
/// alive.

void TMemFile::ResetAfterMerge(TFileMergeInfo *info)
{
   ResetObjects(this,info);

   fNbytesKeys = 0;
   fSeekKeys = 0;

   fMustFlush = kTRUE;
   fInitDone = kFALSE;

   if (fFree) {
      fFree->Delete();
      delete fFree;
      fFree      = nullptr;
   }
   fWritten      = 0;
   fSumBuffer    = 0;
   fSum2Buffer   = 0;
   fBytesRead    = 0;
   fBytesReadExtra = 0;
   fBytesWrite   = 0;
   delete fClassIndex;
   fClassIndex   = nullptr;
   fSeekInfo     = 0;
   fNbytesInfo   = 0;
   delete fProcessIDs;
   fProcessIDs   = nullptr;
   fNProcessIDs  = 0;
   fOffset       = 0;
   fCacheRead    = 0;
   fCacheWrite   = 0;
   fReadCalls    = 0;
   if (fFree) {
      fFree->Delete();
      delete fFree;
      fFree = nullptr;
   }

   fSysOffset   = 0;
   fBlockSeek   = &fBlockList;
   fBlockOffset = 0;
   {
      R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfFiles()->Remove(this);
   }

   {
      TDirectory::TContext ctxt(this);
      Init(kTRUE);

      // And now we need re-initilize the directories ...

      TIter   next(this->GetList());
      TObject *idcur;
      while ((idcur = next())) {
         if (idcur->IsA() == TDirectoryFile::Class()) {
            ((TDirectoryFile*)idcur)->ResetAfterMerge(info);
         }
      }

   }
}

////////////////////////////////////////////////////////////////////////////////
/// Wipe all the data from the permanent buffer but keep, the in-memory object
/// alive.

void TMemFile::ResetObjects(TDirectoryFile *directory, TFileMergeInfo *info) const
{
   if (directory->GetListOfKeys()) {
      TIter next(directory->GetListOfKeys());
      TKey *key;
      while( (key = (TKey*)next()) ) {
         if (nullptr == directory->GetList()->FindObject(key->GetName())) {
            Warning("ResetObjects","Key/Object %s is not attached to the directory %s and can not be ResetAfterMerge correctly",
                    key->GetName(),directory->GetName());
         }
      }
      directory->GetListOfKeys()->Delete("slow");
   }

   TString listHargs;
   listHargs.Form("(TFileMergeInfo*)0x%lx",(ULong_t)info);

   TIter   next(directory->GetList());
   TObject *idcur;
   while ((idcur = next())) {
      TClass *objcl = idcur->IsA();
      if (objcl == TDirectoryFile::Class()) {
         ResetObjects((TDirectoryFile*)idcur,info);
      } else if (objcl->GetResetAfterMerge()) {
         (objcl->GetResetAfterMerge())(idcur,info);
      } else if (idcur->IsA()->GetMethodWithPrototype("ResetAfterMerge", "TFileMergeInfo*") ) {
         Int_t error = 0;
         idcur->Execute("ResetAfterMerge", listHargs.Data(), &error);
         if (error) {
            Error("ResetObjects", "calling ResetAfterMerge() on '%s' failed.",
                  idcur->GetName());
         }
      } else {
//         Error("ResetObjects","In %s, we can not reset %s (not ResetAfterMerge function)",
//               directory->GetName(),idcur->GetName());
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read specified number of bytes from current offset into the buffer.
/// See documentation for TFile::SysRead().

Int_t TMemFile::SysReadImpl(Int_t, void *buf, Long64_t len)
{
   TRACE("READ")

   if (fBlockSeek == nullptr || fBlockSeek->fBuffer == nullptr) {
      errno = EBADF;
      gSystem->SetErrorStr("The memory file is not open.");
      return 0;
   } else {
      // Don't read past the end.
      if (fSysOffset + len > fSize) {
         len = fSize - fSysOffset;
      }

      if (fBlockOffset+len <= fBlockSeek->fSize) {
         // 'len' does not go past the end of the current block,
         // so let's make a simple copy.
         memcpy(buf,fBlockSeek->fBuffer+fBlockOffset,len);
         fBlockOffset += len;
      } else {
         // We are going to have to copy data from more than one
         // block.

         // First copy the end of the first block.
         Int_t sublen = fBlockSeek->fSize - fBlockOffset;
         memcpy(buf,fBlockSeek->fBuffer+fBlockOffset,sublen);

         // Move to the next.
         buf = (char*)buf + sublen;
         Int_t len_left = len - sublen;
         fBlockSeek = fBlockSeek->fNext;

         // Copy all the full blocks that are covered by the request.
         while (len_left > fBlockSeek->fSize) {
            R__ASSERT(fBlockSeek);

            memcpy(buf, fBlockSeek->fBuffer, fBlockSeek->fSize);
            buf = (char*)buf + fBlockSeek->fSize;
            len_left -= fBlockSeek->fSize;
            fBlockSeek = fBlockSeek->fNext;
         }

         // Copy the data from the last block.
         R__ASSERT(fBlockSeek);
         memcpy(buf,fBlockSeek->fBuffer, len_left);
         fBlockOffset = len_left;

      }
      fSysOffset += len;
      return len;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Read specified number of bytes from current offset into the buffer.
/// See documentation for TFile::SysRead().

Int_t TMemFile::SysRead(Int_t fd, void *buf, Int_t len)
{
   return SysReadImpl(fd, buf, len);
}

////////////////////////////////////////////////////////////////////////////////
/// Seek to a specified position in the file.  See TFile::SysSeek().
/// Note that TMemFile does not support seeks when the file is open for write.

Long64_t TMemFile::SysSeek(Int_t, Long64_t offset, Int_t whence)
{
   TRACE("SEEK")
   if (whence == SEEK_SET) {
      fSysOffset = offset;
      fBlockSeek = &fBlockList;
      Long64_t counter = 0;
      while(fBlockSeek->fNext && (counter+fBlockSeek->fSize) < fSysOffset)
      {
         counter += fBlockSeek->fSize;
         fBlockSeek = fBlockSeek->fNext;
      }
      fBlockOffset = fSysOffset - counter;  // If we seek past the 'end' of the file, we now have fBlockOffset > fBlockSeek->fSize
   } else if (whence == SEEK_CUR) {

      if (offset == 0) {
         // nothing to do, really
      } else if (offset > 0) {
         // Move forward.
         if ( (fBlockOffset+offset) < fBlockSeek->fSize) {
            fSysOffset += offset;
            fBlockOffset += offset;
         } else {
            Long64_t counter = fSysOffset;
            fSysOffset += offset;
            while(fBlockSeek->fNext && counter < fSysOffset)
            {
               counter += fBlockSeek->fSize;
               fBlockSeek = fBlockSeek->fNext;
            }
            fBlockOffset = fSysOffset - counter; // If we seek past the 'end' of the file, we now have fBlockOffset > fBlockSeek->fSize
         }
      } else {
         // Move backward in the file (offset < 0).
         Long64_t counter = fSysOffset;
         fSysOffset += offset;
         if (fSysOffset < 0) {
            SysError("TMemFile", "Unable to seek past the beginning of file");
            fSysOffset   = 0;
            fBlockSeek   = &fBlockList;
            fBlockOffset = 0;
            return -1;
         } else {
            if (offset+fBlockOffset >= 0) {
               // We are just moving in the current block.
               fBlockOffset += offset;
            } else {
               while(fBlockSeek->fPrevious && counter > fSysOffset)
               {
                  counter -= fBlockSeek->fSize;
                  fBlockSeek = fBlockSeek->fPrevious;
               }
               fBlockOffset = fSysOffset - counter;
            }
         }
      }
   } else if (whence == SEEK_END) {
      if (offset > 0) {
         SysError("TMemFile", "Unable to seek past end of file");
         return -1;
      }
      if (fSize == -1) {
         SysError("TMemFile", "Unable to seek to end of file");
         return -1;
      }
      fSysOffset = fSize;
   } else {
      SysError("TMemFile", "Unknown whence!");
      return -1;
   }
   return fSysOffset;
}

////////////////////////////////////////////////////////////////////////////////
/// Open a file in 'MemFile'.

Int_t TMemFile::SysOpen(const char * /* pathname */, Int_t /* flags */, UInt_t /* mode */)
{
   if (!fBlockList.fBuffer) {
      fBlockList.fBuffer = new UChar_t[fDefaultBlockSize];
      fBlockList.fSize = fDefaultBlockSize;
      fSize = fDefaultBlockSize;
   }
   if (fBlockList.fBuffer) {
      return 0;
   } else {
      return -1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Close the mem file.

Int_t TMemFile::SysClose(Int_t /* fd */)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Write a buffer into the file.

Int_t TMemFile::SysWriteImpl(Int_t /* fd */, const void *buf, Long64_t len)
{
   TRACE("WRITE")

   if (IsExternalData()) {
      gSystem->SetErrorStr("A memory file with shared data is read-only.");
      return 0;
   }

   if (fBlockList.fBuffer == 0) {
      errno = EBADF;
      gSystem->SetErrorStr("The memory file is not open.");
      return 0;
   } else {
      if (fBlockOffset+len <= fBlockSeek->fSize) {
         // 'len' does not go past the end of the current block,
         // so let's make a simple copy.
         memcpy(fBlockSeek->fBuffer+fBlockOffset,buf,len);
         fBlockOffset += len;
      } else {
         // We are going to have to copy data into more than one
         // block.

         // First copy to the end of the first block.
         Int_t sublen = fBlockSeek->fSize - fBlockOffset;
         memcpy(fBlockSeek->fBuffer+fBlockOffset,buf,sublen);

         // Move to the next.
         buf = (char*)buf + sublen;
         Int_t len_left = len - sublen;
         if (!fBlockSeek->fNext) {
            fBlockSeek->CreateNext(fDefaultBlockSize);
            fSize += fDefaultBlockSize;
         }
         fBlockSeek = fBlockSeek->fNext;

         // Copy all the full blocks that are covered by the request.
         while (len_left > fBlockSeek->fSize) {
            R__ASSERT(fBlockSeek);

            memcpy(fBlockSeek->fBuffer, buf, fBlockSeek->fSize);
            buf = (char*)buf + fBlockSeek->fSize;
            len_left -= fBlockSeek->fSize;
            if (!fBlockSeek->fNext) {
               fBlockSeek->CreateNext(fDefaultBlockSize);
               fSize += fDefaultBlockSize;
            }
            fBlockSeek = fBlockSeek->fNext;
         }

         // Copy the data from the last block.
         R__ASSERT(fBlockSeek);
         memcpy(fBlockSeek->fBuffer, buf, len_left);
         fBlockOffset = len_left;

      }
      fSysOffset += len;
      return len;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write a buffer into the file.

Int_t TMemFile::SysWrite(Int_t fd, const void *buf, Int_t len)
{
   return SysWriteImpl(fd,buf,len);
}

////////////////////////////////////////////////////////////////////////////////
/// Perform a stat on the file; see TFile::SysStat().

Int_t TMemFile::SysStat(Int_t, Long_t* /* id */, Long64_t* /* size */, Long_t* /* flags */, Long_t* /* modtime */)
{
   MayNotUse("SysStat");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Sync remaining data to disk.
/// Nothing to do here.

Int_t TMemFile::SysSync(Int_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Simply calls TSystem::ResetErrno().

void TMemFile::ResetErrno() const
{
   TSystem::ResetErrno();
}
