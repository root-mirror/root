// @(#)root/io:$Id$
// Author: Rene Brun   28/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TFree.h"
#include "TList.h"
#include "TFile.h"
#include "Bytes.h"
#include <iostream>

ClassImp(TFree);

/**
\class TFree
\ingroup IO
Service class for TFile.

Each file has a linked list of free segments. Each free segment is described
by its firts and last address.
When an object is written to a file, a new Key (see TKey)
is created. The first free segment big enough to accomodate the object
is used.
If the object size has a length corresponding to the size of the free segment,
the free segment is deleted from the list of free segments.
When an object is deleted from a file, a new TFree object is generated.
If the deleted object is contiguous to an already deleted object, the free
segments are merged in one single segment.
*/

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TFree::TFree()
{
   fFirst = fLast = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor for a free segment.

TFree::TFree(TList *lfree, Long64_t first, Long64_t last)
{
   fFirst = first;
   fLast  = last;
   lfree->Add(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a new free segment to the list of free segments.
///
///   - if last just precedes an existing free segment, then first becomes
///     the new starting location of the free segment.
///   - if first just follows an existing free segment, then last becomes
///     the new ending location of the free segment.
///   - if first just follows an existing free segment AND last just precedes
///     an existing free segment, these two segments are merged into
///     one single segment.
///

TFree *TFree::AddFree(TList *lfree, Long64_t first, Long64_t last)
{
   TFree *idcur = this;
   while (idcur) {
      Long64_t curfirst = idcur->GetFirst();
      Long64_t curlast  = idcur->GetLast();
      if (curlast == first-1) {
         idcur->SetLast(last);
         TFree *idnext = (TFree*)lfree->After(idcur);
         if (idnext == 0) return idcur;
         if (idnext->GetFirst() > last+1) return idcur;
         idcur->SetLast( idnext->GetLast() );
         lfree->Remove(idnext);
         delete idnext;
         return idcur;
      }
      if (curfirst == last+1) {
         idcur->SetFirst(first);
         return idcur;
      }
      if (first < curfirst) {
         TFree * newfree = new TFree();
         newfree->SetFirst(first);
         newfree->SetLast(last);
         lfree->AddBefore(idcur, newfree);
         return newfree;
      }
      idcur = (TFree*)lfree->After(idcur);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TFree::~TFree()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Encode fre structure into output buffer.

void TFree::FillBuffer(char *&buffer)
{
   Version_t version = TFree::Class_Version();
   if (fLast > TFile::kStartBigFile) version += 1000;
   tobuf(buffer, version);
   // printf("TFree::fillBuffer, fFirst=%lld, fLast=%lld, version=%d\n",fFirst,fLast,version);
   if (version > 1000) {
      tobuf(buffer, fFirst);
      tobuf(buffer, fLast);
   } else {
      tobuf(buffer, (Int_t)fFirst);
      tobuf(buffer, (Int_t)fLast);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the best free segment where to store nbytes.

TFree *TFree::GetBestFree(TList *lfree, Int_t nbytes)
{
   TFree *idcur = this;
   if (idcur == 0) return 0;
   TFree *idcur1 = 0;
   do {
      Long64_t nleft = Long64_t(idcur->fLast - idcur->fFirst +1);
      if (nleft == nbytes) {
         // Found an exact match
         return idcur;
      }
      if(nleft > (Long64_t)(nbytes+3)) {
         if (idcur1 == 0) {
            idcur1=idcur;
         }
      }
      idcur = (TFree*)lfree->After(idcur);
   } while (idcur !=0);

   // return first segment >nbytes
   if (idcur1) return idcur1;

   // try big file
   idcur = (TFree*)lfree->Last();
   Long64_t last = idcur->fLast+1000000000LL;
   idcur->SetLast(last);
   return idcur;
}

////////////////////////////////////////////////////////////////////////////////
/// List free segment contents.

void TFree::ls(Option_t *) const
{
   std::cout <<"Free Segment: "<<fFirst<<"\t"<<fLast<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Decode one free structure from input buffer

void TFree::ReadBuffer(char *&buffer)
{
   Version_t version;
   frombuf(buffer, &version);
   if (version > 1000) {
      frombuf(buffer, &fFirst);
      frombuf(buffer, &fLast);
   } else {
      Int_t first,last;
      frombuf(buffer, &first);  fFirst = (Long64_t)first;
      frombuf(buffer, &last);   fLast  = (Long64_t)last;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// return number of bytes occupied by this TFree on permanent storage

Int_t TFree::Sizeof() const
{
   // printf("TFree::Sizeof, fFirst=%lld, fLast=%lld, version=%d\n",fFirst,fLast, (fLast > TFile::kStartBigFile));
   if (fLast > TFile::kStartBigFile) return 18;
   else                              return 10;
}

