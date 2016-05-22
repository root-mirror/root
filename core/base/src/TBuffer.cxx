// @(#)root/base:$Id: 6da0b5b613bbcfaa3a5cd4074e7b2be2448dfb31 $
// Author: Fons Rademakers   04/05/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TBuffer
Buffer base class used for serializing objects.
*/

#include "TBuffer.h"
#include "TClass.h"
#include "TProcessID.h"
#include "TROOT.h"

const Int_t  kExtraSpace        = 8;   // extra space at end of buffer (used for free block count)

ClassImp(TBuffer)

////////////////////////////////////////////////////////////////////////////////
/// The user has provided memory than we don't own, thus we can not extent it
/// either.

static char *R__NoReAllocChar(char *, size_t, size_t)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create an I/O buffer object. Mode should be either TBuffer::kRead or
/// TBuffer::kWrite. By default the I/O buffer has a size of
/// TBuffer::kInitialSize (1024) bytes.
///
/// \param[in] mode Read mode or write mode of this TBuffer
/// \param[in] def  Determining whether this TKey is created with endianness defined by global gROOT or defined by buffBigEndian
/// \param[in] buffBigEndian Determining the endianness of this TKey's buffer only if def is kFALSE
///
/// TBuffer might store StreamerInfo which is always stored as big endian on TFile. Therefore, we should get rid of the effect
/// of gROOT->IsBufBigEndian() and always define its buffer as big endian. For more explanation, take a look at the function:
///
///     TKey::TKey(TDirectory* motherDir, const TKey &orig, UShort_t pidOffset, Bool_t def, Bool_t buffBigEndian);
///
/// It provides details of why we need 'def' and 'buffBigEndian' and how to define them.

TBuffer::TBuffer(EMode mode, Bool_t def, Bool_t buffBigEndian)
{
   fBufSize      = kInitialSize;
   fMode         = mode;
   fVersion      = 0;
   fParent       = 0;
   if (def) {
      fBufBigEndian = gROOT->IsBufBigEndian();
   } else {
      if (buffBigEndian)
         fBufBigEndian = kTRUE;
      else
         fBufBigEndian = kFALSE;
   }

   SetBit(kIsOwner);

   fBuffer = new char[fBufSize+kExtraSpace];

   fBufCur = fBuffer;
   fBufMax = fBuffer + fBufSize;

   SetReAllocFunc( 0 );
}

////////////////////////////////////////////////////////////////////////////////
/// Create an I/O buffer object. Mode should be either TBuffer::kRead or
/// TBuffer::kWrite.

TBuffer::TBuffer(EMode mode, Int_t bufsiz, Bool_t def, Bool_t buffBigEndian)
{
   if (bufsiz < kMinimalSize) bufsiz = kMinimalSize;
   fBufSize      = bufsiz;
   fMode         = mode;
   fVersion      = 0;
   fParent       = 0;
   if (def) {
      fBufBigEndian = gROOT->IsBufBigEndian();
   } else {
      if (buffBigEndian) 
         fBufBigEndian = kTRUE;
      else
         fBufBigEndian = kFALSE;
   }

   SetBit(kIsOwner);

   fBuffer = new char[fBufSize+kExtraSpace];

   fBufCur = fBuffer;
   fBufMax = fBuffer + fBufSize;

   SetReAllocFunc( 0 );
}

////////////////////////////////////////////////////////////////////////////////
/// Create an I/O buffer object. Mode should be either TBuffer::kRead or
/// TBuffer::kWrite. By default the I/O buffer has a size of
/// TBuffer::kInitialSize (1024) bytes. An external buffer can be passed
/// to TBuffer via the buf argument. By default this buffer will be adopted
/// unless adopt is false.
///
/// If the new buffer is _not_ adopted and no memory allocation routine
/// is provided, a Fatal error will be issued if the Buffer attempts to
/// expand.

TBuffer::TBuffer(EMode mode, Int_t bufsiz, void *buf, Bool_t adopt, ReAllocCharFun_t reallocfunc, Bool_t def, Bool_t buffBigEndian)
{
   fBufSize      = bufsiz;
   fMode         = mode;
   fVersion      = 0;
   fParent       = 0;
   if (def) {
      fBufBigEndian = gROOT->IsBufBigEndian();
   } else {
      if (buffBigEndian) 
         fBufBigEndian = kTRUE;
      else
         fBufBigEndian = kFALSE;
   }

   SetBit(kIsOwner);

   if (buf) {
      fBuffer = (char *)buf;
      if ( (fMode&kWrite)!=0 ) {
         fBufSize -= kExtraSpace;
      }
      if (!adopt) ResetBit(kIsOwner);
   } else {
      if (fBufSize < kMinimalSize) {
         fBufSize = kMinimalSize;
      }
      fBuffer = new char[fBufSize+kExtraSpace];
   }
   fBufCur = fBuffer;
   fBufMax = fBuffer + fBufSize;

   SetReAllocFunc( reallocfunc );

   if (buf && ( (fMode&kWrite)!=0 ) && fBufSize < 0) {
      Expand( kMinimalSize );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Delete an I/O buffer object.

TBuffer::~TBuffer()
{
   if (TestBit(kIsOwner)) {
      //printf("Deleting fBuffer=%lx\n", fBuffer);
      delete [] fBuffer;
   }
   fBuffer = 0;
   fParent = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Automatically calculate a new size and expand the buffer to fit at least size_needed.
/// The goals is to minimize the number of memory allocation and the memory allocation
/// which avoiding too much memory wastage.
///
/// If the size_needed is larger than the current size, the policy
/// is to expand to double the current size or the size_needed which ever is largest.

void TBuffer::AutoExpand(Int_t size_needed)
{
   if (size_needed > fBufSize) {
      if (size_needed > 2*fBufSize) {
         Expand(size_needed);
      } else {
         Expand(2*fBufSize);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Sets a new buffer in an existing TBuffer object. If newsiz=0 then the
/// new buffer is expected to have the same size as the previous buffer.
/// The current buffer position is reset to the start of the buffer.
/// If the TBuffer owned the previous buffer, it will be deleted prior
/// to accepting the new buffer. By default the new buffer will be
/// adopted unless adopt is false.
///
/// If the new buffer is _not_ adopted and no memory allocation routine
/// is provided, a Fatal error will be issued if the Buffer attempts to
/// expand.

void TBuffer::SetBuffer(void *buf, UInt_t newsiz, Bool_t adopt, ReAllocCharFun_t reallocfunc)
{
   if (fBuffer && TestBit(kIsOwner))
      delete [] fBuffer;

   if (adopt)
      SetBit(kIsOwner);
   else
      ResetBit(kIsOwner);

   fBuffer = (char *)buf;
   fBufCur = fBuffer;
   if (newsiz > 0) {
      if ( (fMode&kWrite)!=0 ) {
         fBufSize = newsiz - kExtraSpace;
      } else {
         fBufSize = newsiz;
      }
   }
   fBufMax = fBuffer + fBufSize;

   SetReAllocFunc( reallocfunc );

   if (buf && ( (fMode&kWrite)!=0 ) && fBufSize < 0) {
      Expand( kMinimalSize );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Expand (or shrink) the I/O buffer to newsize bytes.
/// If copy is true (the default), the existing content of the
/// buffer is preserved, otherwise the buffer is returned zero-ed out.
///
/// In order to avoid losing data, if the current length is greater than
/// the requested size, we only shrink down to the current length.

void TBuffer::Expand(Int_t newsize, Bool_t copy)
{
   Int_t l  = Length();
   if ( (l > newsize) && copy ) {
      newsize = l;
   }
   if ( (fMode&kWrite)!=0 ) {
      fBuffer  = fReAllocFunc(fBuffer, newsize+kExtraSpace,
                              copy ? fBufSize+kExtraSpace : 0);
   } else {
      fBuffer  = fReAllocFunc(fBuffer, newsize,
                              copy ? fBufSize : 0);
   }
   if (fBuffer == 0) {
      if (fReAllocFunc == TStorage::ReAllocChar) {
         Fatal("Expand","Failed to expand the data buffer using TStorage::ReAllocChar.");
      } if (fReAllocFunc == R__NoReAllocChar) {
         Fatal("Expand","Failed to expand the data buffer because TBuffer does not own it and no custom memory reallocator was provided.");
      } else {
         Fatal("Expand","Failed to expand the data buffer using custom memory reallocator 0x%lx.", (Long_t)fReAllocFunc);
      }
   }
   fBufSize = newsize;
   fBufCur  = fBuffer + l;
   fBufMax  = fBuffer + fBufSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to parent of this buffer.

TObject *TBuffer::GetParent() const
{
   return fParent;
}

////////////////////////////////////////////////////////////////////////////////
/// Set parent owning this buffer.

void TBuffer::SetParent(TObject *parent)
{
   fParent = parent;
}
////////////////////////////////////////////////////////////////////////////////
/// Return the reallocation method currently used.

ReAllocCharFun_t TBuffer::GetReAllocFunc() const
{
   return fReAllocFunc;
}

////////////////////////////////////////////////////////////////////////////////
/// Set which memory reallocation method to use.  If reallocafunc is null,
/// reset it to the default value (TStorage::ReAlloc)

void  TBuffer::SetReAllocFunc(ReAllocCharFun_t reallocfunc )
{
   if (reallocfunc) {
      fReAllocFunc = reallocfunc;
   } else {
      if (TestBit(kIsOwner)) {
         fReAllocFunc = TStorage::ReAllocChar;
      } else {
         fReAllocFunc = R__NoReAllocChar;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set buffer in read mode.

void TBuffer::SetReadMode()
{
   if ( (fMode&kWrite)!=0 ) {
      // We had reserved space for the free block count,
      // release it,
      fBufSize += kExtraSpace;
   }
   fMode = kRead;
}

////////////////////////////////////////////////////////////////////////////////
/// Set buffer in write mode.

void TBuffer::SetWriteMode()
{
   if ( (fMode&kWrite)==0 ) {
      // We had not yet reserved space for the free block count,
      // reserve it now.
      fBufSize -= kExtraSpace;
   }
   fMode = kWrite;
}

////////////////////////////////////////////////////////////////////////////////
/// Forward to TROOT::GetClass().

TClass *TBuffer::GetClass(const std::type_info &typeinfo)
{
   return TClass::GetClass(typeinfo);
}

////////////////////////////////////////////////////////////////////////////////
/// Forward to TROOT::GetClass().

TClass *TBuffer::GetClass(const char *className)
{
   return TClass::GetClass(className);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the current Process-ID.

TProcessID *TBuffer::ReadProcessID(UShort_t pidf)
{
   if (!pidf) return TProcessID::GetPID(); //may happen when cloning an object
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Always return 0 (current processID).

UShort_t TBuffer::WriteProcessID(TProcessID *)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Push a new data cache area onto the list of area to be used for
/// temporarily store 'missing' data members.

void TBuffer::PushDataCache(TVirtualArray *obj)
{
   fCacheStack.push_back(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the 'current' data cache area from the list of area to be used for
/// temporarily store 'missing' data members.

TVirtualArray *TBuffer::PeekDataCache() const
{
   if (fCacheStack.empty()) return 0;
   return fCacheStack.back();
}

////////////////////////////////////////////////////////////////////////////////
/// Pop and Return the 'current' data cache area from the list of area to be used for
/// temporarily store 'missing' data members.

TVirtualArray *TBuffer::PopDataCache()
{
   TVirtualArray *val = PeekDataCache();
   fCacheStack.pop_back();
   return val;
}

