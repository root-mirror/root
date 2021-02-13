// @(#)root/tree:$Id$
// Author: Enrico Guiraud

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TLeafG
\ingroup tree

A TLeaf for a 64 bit Integer data type.
*/

#include "TLeafG.h"
#include "TBranch.h"
#include "TBuffer.h"
#include "TClonesArray.h"
#include <iostream>

ClassImp(TLeafG);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor for LeafL.

TLeafG::TLeafG(): TLeaf()
{
   fLenType = sizeof(Long_t);
   fMinimum = 0;
   fMaximum = 0;
   fValue   = 0;
   fPointer = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a LeafL.

TLeafG::TLeafG(TBranch *parent, const char *name, const char *type)
   :TLeaf(parent, name,type)
{
   fLenType = sizeof(Long_t);
   fMinimum = 0;
   fMaximum = 0;
   fValue   = 0;
   fPointer = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Default destructor for a LeafL.

TLeafG::~TLeafG()
{
   if (ResetAddress(0,kTRUE)) delete [] fValue;
}

////////////////////////////////////////////////////////////////////////////////
/// Export element from local leaf buffer to ClonesArray.

void TLeafG::Export(TClonesArray *list, Int_t n)
{
   Long_t *value = fValue;
   for (Int_t i=0;i<n;i++) {
      char *first = (char*)list->UncheckedAt(i);
      Long_t *ii = (Long_t*)&first[fOffset];
      for (Int_t j=0;j<fLen;j++) {
         ii[j] = value[j];
      }
      value += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Pack leaf elements in Basket output buffer.

void TLeafG::FillBasket(TBuffer &b)
{
   Int_t i;
   Int_t len = GetLen();
   if (fPointer) fValue = *fPointer;
   if (IsRange()) {
      if (fValue[0] > fMaximum) fMaximum = fValue[0];
   }
   if (IsUnsigned()) {
      for (i=0;i<len;i++) b << (ULong_t)fValue[i];
   } else {
      b.WriteFastArray(fValue,len);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns name of leaf type.

const char *TLeafG::GetTypeName() const
{
   if (fIsUnsigned) return "ULong_t";
   return "Long_t";
}

////////////////////////////////////////////////////////////////////////////////
/// Returns current value of leaf.
/// - if leaf is a simple type, i must be set to 0
/// - if leaf is an array, i is the array element number to be returned

Double_t TLeafG::GetValue(Int_t i) const
{
   if (fIsUnsigned) return (Double_t)((ULong_t)fValue[i]);
   return fValue[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Returns current value of leaf
/// - if leaf is a simple type, i must be set to 0
/// - if leaf is an array, i is the array element number to be returned

LongDouble_t TLeafG::GetValueLongDouble(Int_t i) const
{
   if (fIsUnsigned) return (LongDouble_t)((ULong_t)fValue[i]);
   return fValue[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Copy/set fMinimum and fMaximum to include/be wide than those of the parameter

Bool_t TLeafG::IncludeRange(TLeaf *input)
{
    if (input) {
        if (input->GetMaximum() > this->GetMaximum())
            this->SetMaximum( input->GetMaximum() );
        if (input->GetMinimum() < this->GetMinimum())
            this->SetMinimum( input->GetMinimum() );
        return kTRUE;
    } else {
        return kFALSE;
    }
}

////////////////////////////////////////////////////////////////////////////////
/// Import element from ClonesArray into local leaf buffer.

void TLeafG::Import(TClonesArray *list, Int_t n)
{
   const Int_t kIntUndefined = -9999;
   Int_t j = 0;
   char *clone;
   for (Int_t i=0;i<n;i++) {
      clone = (char*)list->UncheckedAt(i);
      if (clone) memcpy(&fValue[j],clone + fOffset, 8*fLen);
      else       memcpy(&fValue[j],&kIntUndefined,  8*fLen);
      j += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Prints leaf value.

void TLeafG::PrintValue(Int_t l) const
{
   if (fIsUnsigned) {
      ULong_t *uvalue = (ULong_t*)GetValuePointer();
      printf("%lu",uvalue[l]);
   } else {
      Long_t *value = (Long_t*)GetValuePointer();
      printf("%ld",value[l]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read leaf elements from Basket input buffer.

void TLeafG::ReadBasket(TBuffer &b)
{
   if (!fLeafCount && fNdata == 1) {
      b.ReadLong(fValue[0]);
   } else {
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

////////////////////////////////////////////////////////////////////////////////
/// Deserialize input by performing byteswap as needed.
bool TLeafG::ReadBasketFast(TBuffer& input_buf, Long64_t N)
{
   if (R__unlikely(fLeafCount)) {return false;}
   return input_buf.ByteSwapBuffer(fLen*N, kLong_t);
}

////////////////////////////////////////////////////////////////////////////////
/// Read leaf elements from Basket input buffer and export buffer to
/// TClonesArray objects.

void TLeafG::ReadBasketExport(TBuffer &b, TClonesArray *list, Int_t n)
{
   if (n*fLen == 1) {
      b >> fValue[0];
   } else {
      b.ReadFastArray(fValue,n*fLen);
   }
   Long_t *value = fValue;
   for (Int_t i=0;i<n;i++) {
      char *first = (char*)list->UncheckedAt(i);
      Long_t *ii = (Long_t*)&first[fOffset];
      for (Int_t j=0;j<fLen;j++) {
         ii[j] = value[j];
      }
      value += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read a long integer from std::istream s and store it into the branch buffer.

void TLeafG::ReadValue(std::istream &s, Char_t /*delim = ' '*/)
{
#if defined(_MSC_VER) && (_MSC_VER<1300)
   printf("Due to a bug in VC++6, the function TLeafG::ReadValue is dummy\n");
#else
   if (fIsUnsigned) {
      ULong_t *uvalue = (ULong_t*)GetValuePointer();
      for (Int_t i=0;i<fLen;i++) s >> uvalue[i];
   } else {
      Long_t *value = (Long_t*)GetValuePointer();
      for (Int_t i=0;i<fLen;i++) s >> value[i];
   }
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Set leaf buffer data address.

void TLeafG::SetAddress(void *add)
{
   if (ResetAddress(add) && (add!=fValue)) {
      delete [] fValue;
   }
   if (add) {
      if (TestBit(kIndirectAddress)) {
         fPointer = (Long_t**) add;
         Int_t ncountmax = fLen;
         if (fLeafCount) ncountmax = fLen*(fLeafCount->GetMaximum() + 1);
         if ((fLeafCount && ncountmax > Int_t(fLeafCount->GetValue())) ||
             ncountmax > fNdata || *fPointer == 0) {
            if (*fPointer) delete [] *fPointer;
            if (ncountmax > fNdata) fNdata = ncountmax;
            *fPointer = new Long_t[fNdata];
         }
         fValue = *fPointer;
      } else {
         fValue = (Long_t*)add;
      }
   } else {
      fValue = new Long_t[fNdata];
      fValue[0] = 0;
   }
}

