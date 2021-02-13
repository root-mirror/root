// @(#)root/tree:$Id$
// Author: Rene Brun   12/01/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TLeafD
\ingroup tree

A TLeaf for a 64 bit floating point data type.
*/

#include "TLeafD.h"
#include "TBranch.h"
#include "TBuffer.h"
#include "TClonesArray.h"
#include <iostream>

ClassImp(TLeafD);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor for LeafD.

TLeafD::TLeafD(): TLeaf()
{
   fLenType = 8;
   fMinimum = 0;
   fMaximum = 0;
   fValue   = 0;
   fPointer = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a LeafD.

TLeafD::TLeafD(TBranch *parent, const char *name, const char *type)
   :TLeaf(parent, name,type)
{
   fLenType = 8;
   fMinimum = 0;
   fMaximum = 0;
   fValue   = 0;
   fPointer = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Default destructor for a LeafD.

TLeafD::~TLeafD()
{
   if (ResetAddress(0,kTRUE)) delete [] fValue;
}

////////////////////////////////////////////////////////////////////////////////
/// Export element from local leaf buffer to ClonesArray.

void TLeafD::Export(TClonesArray *list, Int_t n)
{
   Int_t j = 0;
   for (Int_t i=0;i<n;i++) {
      memcpy((char*)list->UncheckedAt(i) + fOffset,&fValue[j], 8*fLen);
      j += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Pack leaf elements in Basket output buffer.

void TLeafD::FillBasket(TBuffer &b)
{
   Int_t len = GetLen();
   if (fPointer) fValue = *fPointer;
   b.WriteFastArray(fValue,len);
}

////////////////////////////////////////////////////////////////////////////////
/// Import element from ClonesArray into local leaf buffer.

void TLeafD::Import(TClonesArray *list, Int_t n)
{
   const Double_t kDoubleUndefined = -9999.;
   Int_t j = 0;
   char *clone;
   for (Int_t i=0;i<n;i++) {
      clone = (char*)list->UncheckedAt(i);
      if (clone) memcpy(&fValue[j],clone + fOffset, 8*fLen);
      else       memcpy(&fValue[j],&kDoubleUndefined,  8*fLen);
      j += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Prints leaf value.

void TLeafD::PrintValue(Int_t l) const
{
   Double_t *value = (Double_t *)GetValuePointer();
   printf("%g",value[l]);
}

////////////////////////////////////////////////////////////////////////////////
/// Read leaf elements from Basket input buffer.

void TLeafD::ReadBasket(TBuffer &b)
{
   if (!fLeafCount && fNdata == 1) {
      b.ReadDouble(fValue[0]);
   }else {
      if (fLeafCount) {
         Long64_t entry = fBranch->GetReadEntry();
         if (fLeafCount->GetBranch()->GetReadEntry() != entry) {
            fLeafCount->GetBranch()->GetEntry(entry);
         }
         Int_t len = Int_t(fLeafCount->GetValue());
         if (len > fLeafCount->GetMaximum()) {
            printf("ERROR leaf:%s, len=%d and max=%d\n",GetName(),len,fLeafCount->GetMaximum());
            len = fLeafCount->GetMaximum();
         }
         fNdata = len*fLen;
         b.ReadFastArray(fValue,len*fLen);
      } else {
         b.ReadFastArray(fValue,fLen);
      }
   }
}

// Deserialize N events from an input buffer.
bool TLeafD::ReadBasketFast(TBuffer &input_buf, Long64_t N) {
   if (R__unlikely(fLeafCount)) { return false; }
   return input_buf.ByteSwapBuffer(fLen*N, kDouble_t);
}

////////////////////////////////////////////////////////////////////////////////
/// Read leaf elements from Basket input buffer and export buffer to
/// TClonesArray objects.

void TLeafD::ReadBasketExport(TBuffer &b, TClonesArray *list, Int_t n)
{
   b.ReadFastArray(fValue,n*fLen);

   Int_t j = 0;
   for (Int_t i=0;i<n;i++) {
      memcpy((char*)list->UncheckedAt(i) + fOffset,&fValue[j], 8*fLen);
      j += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read a double from std::istream s and store it into the branch buffer.

void TLeafD::ReadValue(std::istream &s, Char_t /*delim = ' '*/)
{
   Double_t *value = (Double_t*)GetValuePointer();
   for (Int_t i=0;i<fLen;i++) s >> value[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Set leaf buffer data address.

void TLeafD::SetAddress(void *add)
{
   if (ResetAddress(add) && (add!= fValue)) {
      delete [] fValue;
   }
   if (add) {
      if (TestBit(kIndirectAddress)) {
         fPointer = (Double_t**) add;
         Int_t ncountmax = fLen;
         if (fLeafCount) ncountmax = fLen*(fLeafCount->GetMaximum() + 1);
         if ((fLeafCount && ncountmax > Int_t(fLeafCount->GetValue())) ||
             ncountmax > fNdata || *fPointer == 0) {
            if (*fPointer) delete [] *fPointer;
            if (ncountmax > fNdata) fNdata = ncountmax;
            *fPointer = new Double_t[fNdata];
         }
         fValue = *fPointer;
      } else {
         fValue = (Double_t*)add;
      }
   } else {
      fValue = new Double_t[fNdata];
      fValue[0] = 0;
   }
}
