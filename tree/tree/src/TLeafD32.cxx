// @(#)root/tree:$Id$
// Author: Simon Spies 23/02/19

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TLeafD32
\ingroup tree

A TLeaf for a 24 bit truncated floating point data type.
*/

#include "TLeafD32.h"
#include "TBranch.h"
#include "TBuffer.h"
#include "TClonesArray.h"
#include "TStreamerElement.h"
#include <iostream>

ClassImp(TLeafD32);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor for LeafD32.

TLeafD32::TLeafD32() : TLeaf()
{
   fLenType = 8;
   fMinimum = 0;
   fMaximum = 0;
   fValue = nullptr;
   fPointer = nullptr;
   fElement = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a LeafD32.

TLeafD32::TLeafD32(TBranch *parent, const char *name, const char *type) : TLeaf(parent, name, type)
{
   fLenType = 8;
   fMinimum = 0;
   fMaximum = 0;
   fValue = nullptr;
   fPointer = nullptr;
   fElement = nullptr;
   fTitle = type;

   if (strchr(type, '['))
      fElement = new TStreamerElement(Form("%s_Element", name), type, 0, 0, "Double32_t");
}

////////////////////////////////////////////////////////////////////////////////
/// Default destructor for a LeafD32.

TLeafD32::~TLeafD32()
{
   if (ResetAddress(nullptr, kTRUE))
      delete[] fValue;

   if (fElement)
      delete fElement;
}

////////////////////////////////////////////////////////////////////////////////
/// Export element from local leaf buffer to ClonesArray.

void TLeafD32::Export(TClonesArray *list, Int_t n)
{
   Int_t j = 0;
   for (Int_t i = 0; i < n; i++) {
      memcpy((char *)list->UncheckedAt(i) + fOffset, &fValue[j], 8 * fLen);
      j += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Pack leaf elements in Basket output buffer.

void TLeafD32::FillBasket(TBuffer &b)
{
   Int_t len = GetLen();
   if (fPointer)
      fValue = *fPointer;
   b.WriteFastArrayDouble32(fValue, len, fElement);
}

////////////////////////////////////////////////////////////////////////////////
/// Import element from ClonesArray into local leaf buffer.

void TLeafD32::Import(TClonesArray *list, Int_t n)
{
   const Double32_t kDoubleUndefined = -9999.;
   Int_t j = 0;
   for (Int_t i = 0; i < n; i++) {
      auto clone = (char *)list->UncheckedAt(i);
      if (clone)
         memcpy(&fValue[j], clone + fOffset, 8 * fLen);
      else
         memcpy(&fValue[j], &kDoubleUndefined, 8 * fLen);
      j += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Prints leaf value.

void TLeafD32::PrintValue(Int_t l) const
{
   auto value = (Double32_t *)GetValuePointer();
   printf("%g", value[l]);
}

////////////////////////////////////////////////////////////////////////////////
/// Read leaf elements from Basket input buffer.

void TLeafD32::ReadBasket(TBuffer &b)
{
   if (!fLeafCount && fNdata == 1) {
      b.ReadDouble32(fValue, fElement);
   } else {
      if (fLeafCount) {
         Long64_t entry = fBranch->GetReadEntry();
         if (fLeafCount->GetBranch()->GetReadEntry() != entry) {
            fLeafCount->GetBranch()->GetEntry(entry);
         }
         auto len = Int_t(fLeafCount->GetValue());
         if (len > fLeafCount->GetMaximum()) {
            printf("ERROR leaf:%s, len=%d and max=%d\n", GetName(), len, fLeafCount->GetMaximum());
            len = fLeafCount->GetMaximum();
         }
         fNdata = len * fLen;
         b.ReadFastArrayDouble32(fValue, len * fLen, fElement);
      } else {
         b.ReadFastArrayDouble32(fValue, fLen, fElement);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read leaf elements from Basket input buffer and export buffer to
/// TClonesArray objects.

void TLeafD32::ReadBasketExport(TBuffer &b, TClonesArray *list, Int_t n)
{
   b.ReadFastArrayDouble32(fValue, n * fLen, fElement);

   Int_t j = 0;
   for (Int_t i = 0; i < n; i++) {
      memcpy((char *)list->UncheckedAt(i) + fOffset, &fValue[j], 8 * fLen);
      j += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read a double from std::istream s and store it into the branch buffer.

void TLeafD32::ReadValue(std::istream &s, Char_t /*delim = ' '*/)
{
   auto value = (Double32_t *)GetValuePointer();
   for (Int_t i = 0; i < fLen; i++)
      s >> value[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Set leaf buffer data address.

void TLeafD32::SetAddress(void *add)
{
   if (ResetAddress(add) && (add != fValue)) {
      delete[] fValue;
   }
   if (add) {
      if (TestBit(kIndirectAddress)) {
         fPointer = (Double32_t **)add;
         Int_t ncountmax = fLen;
         if (fLeafCount)
            ncountmax = fLen * (fLeafCount->GetMaximum() + 1);
         if ((fLeafCount && ncountmax > Int_t(fLeafCount->GetValue())) || ncountmax > fNdata || *fPointer == nullptr) {
            if (*fPointer)
               delete[] * fPointer;
            if (ncountmax > fNdata)
               fNdata = ncountmax;
            *fPointer = new Double32_t[fNdata];
         }
         fValue = *fPointer;
      } else {
         fValue = (Double32_t *)add;
      }
   } else {
      fValue = new Double32_t[fNdata];
      fValue[0] = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TLeafD32.

void TLeafD32::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TLeafD32::Class(), this);

      if (fTitle.Contains("["))
	 fElement = new TStreamerElement(Form("%s_Element", fName.Data()), fTitle.Data(), 0, 0, "Double32_t");
   } else {
      R__b.WriteClassBuffer(TLeafD32::Class(), this);
   }
}
