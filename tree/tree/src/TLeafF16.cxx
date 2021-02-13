// @(#)root/tree:$Id$
// Author: Simon Spies 23/02/19

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TLeafF16
\ingroup tree

A TLeaf for a 24 bit truncated floating point data type.
*/

#include "TLeafF16.h"
#include "TBranch.h"
#include "TBuffer.h"
#include "TClonesArray.h"
#include "TStreamerElement.h"
#include <iostream>

ClassImp(TLeafF16);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor for LeafF16.

TLeafF16::TLeafF16() : TLeaf()
{
   fLenType = 4;
   fMinimum = 0;
   fMaximum = 0;
   fValue = nullptr;
   fPointer = nullptr;
   fElement = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a LeafF16.

TLeafF16::TLeafF16(TBranch *parent, const char *name, const char *type) : TLeaf(parent, name, type)
{
   fLenType = 4;
   fMinimum = 0;
   fMaximum = 0;
   fValue = nullptr;
   fPointer = nullptr;
   fElement = nullptr;
   fTitle = type;

   if (strchr(type, '['))
      fElement = new TStreamerElement(Form("%s_Element", name), type, 0, 0, "Float16_t");
}

////////////////////////////////////////////////////////////////////////////////
/// Default destructor for a LeafF16.

TLeafF16::~TLeafF16()
{
   if (ResetAddress(nullptr, kTRUE))
      delete[] fValue;

   if (fElement)
      delete fElement;
}

////////////////////////////////////////////////////////////////////////////////
/// Export element from local leaf buffer to ClonesArray.

void TLeafF16::Export(TClonesArray *list, Int_t n)
{
   Float16_t *value = fValue;
   for (Int_t i = 0; i < n; i++) {
      auto first = (char *)list->UncheckedAt(i);
      auto ff = (Float16_t *)&first[fOffset];
      for (Int_t j = 0; j < fLen; j++) {
         ff[j] = value[j];
      }
      value += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Pack leaf elements in Basket output buffer.

void TLeafF16::FillBasket(TBuffer &b)
{
   Int_t len = GetLen();
   if (fPointer)
      fValue = *fPointer;
   b.WriteFastArrayFloat16(fValue, len, fElement);
}

////////////////////////////////////////////////////////////////////////////////
/// Import element from ClonesArray into local leaf buffer.

void TLeafF16::Import(TClonesArray *list, Int_t n)
{
   const Float16_t kFloatUndefined = -9999.;
   Int_t j = 0;
   for (Int_t i = 0; i < n; i++) {
      auto clone = (char *)list->UncheckedAt(i);
      if (clone)
         memcpy(&fValue[j], clone + fOffset, 4 * fLen);
      else
         memcpy(&fValue[j], &kFloatUndefined, 4 * fLen);
      j += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Prints leaf value.

void TLeafF16::PrintValue(Int_t l) const
{
   auto value = (Float16_t *)GetValuePointer();
   printf("%g", value[l]);
}

////////////////////////////////////////////////////////////////////////////////
/// Read leaf elements from Basket input buffer.

void TLeafF16::ReadBasket(TBuffer &b)
{
   if (!fLeafCount && fNdata == 1) {
      b.ReadFloat16(fValue, fElement);
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
         b.ReadFastArrayFloat16(fValue, len * fLen, fElement);
      } else {
         b.ReadFastArrayFloat16(fValue, fLen, fElement);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read leaf elements from Basket input buffer and export buffer to
/// TClonesArray objects.

void TLeafF16::ReadBasketExport(TBuffer &b, TClonesArray *list, Int_t n)
{
   if (n * fLen == 1) {
      b.ReadFloat16(fValue, fElement);
   } else {
      b.ReadFastArrayFloat16(fValue, n * fLen, fElement);
   }

   Float16_t *value = fValue;
   for (Int_t i = 0; i < n; i++) {
      auto first = (char *)list->UncheckedAt(i);
      auto ff = (Float16_t *)&first[fOffset];
      for (Int_t j = 0; j < fLen; j++) {
         ff[j] = value[j];
      }
      value += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read a float from std::istream s and store it into the branch buffer.

void TLeafF16::ReadValue(std::istream &s, Char_t /*delim = ' '*/)
{
   auto value = (Float16_t *)GetValuePointer();
   for (Int_t i = 0; i < fLen; i++)
      s >> value[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Set leaf buffer data address.

void TLeafF16::SetAddress(void *add)
{
   if (ResetAddress(add) && (add != fValue)) {
      delete[] fValue;
   }

   if (add) {
      if (TestBit(kIndirectAddress)) {
         fPointer = (Float16_t **)add;
         Int_t ncountmax = fLen;
         if (fLeafCount)
            ncountmax = fLen * (fLeafCount->GetMaximum() + 1);
         if ((fLeafCount && ncountmax > Int_t(fLeafCount->GetValue())) || ncountmax > fNdata || *fPointer == nullptr) {
            if (*fPointer)
               delete[] * fPointer;
            if (ncountmax > fNdata)
               fNdata = ncountmax;
            *fPointer = new Float16_t[fNdata];
         }
         fValue = *fPointer;
      } else {
         fValue = (Float16_t *)add;
      }
   } else {
      fValue = new Float16_t[fNdata];
      fValue[0] = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TLeafF16.

void TLeafF16::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TLeafF16::Class(), this);

      if (fTitle.Contains("["))
	 fElement = new TStreamerElement(Form("%s_Element", fName.Data()), fTitle.Data(), 0, 0, "Float16_t");
   } else {
      R__b.WriteClassBuffer(TLeafF16::Class(), this);
   }
}
