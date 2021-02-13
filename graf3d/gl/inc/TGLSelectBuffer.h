// @(#)root/gl:$Id$
// Author:  Matevz Tadel, Feb 2007

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGLSelectBuffer
#define ROOT_TGLSelectBuffer

#include <Rtypes.h>

class TGLSelectRecordBase;

#include <vector>
#include <utility>

/**************************************************************************/
// TGLSelectBuffer
/**************************************************************************/

class TGLSelectBuffer
{
protected:
   Int_t   fBufSize;  // Size of buffer.
   UInt_t* fBuf;      // Actual buffer.

   Int_t   fNRecords; // Number of records as returned by glRenderMode.

   typedef std::pair<UInt_t, UInt_t*>  RawRecord_t;
   typedef std::vector<RawRecord_t>   vRawRecord_t;

   vRawRecord_t fSortedRecords;

   static Int_t fgMaxBufSize;

public:
   TGLSelectBuffer();
   virtual ~TGLSelectBuffer();

   Int_t   GetBufSize()  const { return fBufSize; }
   UInt_t* GetBuf()      const { return fBuf; }
   Int_t   GetNRecords() const { return fNRecords; }

   Bool_t CanGrow();
   void   Grow();

   void ProcessResult(Int_t glResult);

   UInt_t* RawRecord(Int_t i) { return fSortedRecords[i].second; }

   Int_t   SelectRecord(TGLSelectRecordBase& rec, Int_t i);

   ClassDef(TGLSelectBuffer, 0) // OpenGL select buffer with depth sorting.
};

#endif
