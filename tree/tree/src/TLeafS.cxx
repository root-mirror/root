// @(#)root/tree:$Id$
// Author: Rene Brun   12/01/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TLeafS
\ingroup tree

A TLeaf for a 16 bit Integer data type.
*/

#include "TLeafS.h"
#include "TBranch.h"
#include "TBuffer.h"
#include "TClonesArray.h"
#include <iostream>

ClassImp(TLeafS);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor for LeafS.

TLeafS::TLeafS(): TLeaf()
{
   fValue = 0;
   fPointer = 0;
   fMinimum = 0;
   fMaximum = 0;
   fLenType = 2;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a LeafS.

TLeafS::TLeafS(TBranch *parent, const char *name, const char *type)
   :TLeaf(parent,name,type)
{
   fLenType = 2;
   fMinimum = 0;
   fMaximum = 0;
   fValue   = 0;
   fPointer = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Default destructor for a LeafS.

TLeafS::~TLeafS()
{
   if (ResetAddress(0,kTRUE)) delete [] fValue;
}

////////////////////////////////////////////////////////////////////////////////
/// Export element from local leaf buffer to ClonesArray.

void TLeafS::Export(TClonesArray *list, Int_t n)
{
   Int_t j = 0;
   for (Int_t i=0;i<n;i++) {
      memcpy((char*)list->UncheckedAt(i) + fOffset,&fValue[j], 2*fLen);
      j += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Pack leaf elements in Basket output buffer.

void TLeafS::FillBasket(TBuffer &b)
{
   Int_t i;
   Int_t len = GetLen();
   if (fPointer) fValue = *fPointer;
   if (IsRange()) {
      if (fValue[0] > fMaximum) fMaximum = fValue[0];
   }
   if (IsUnsigned()) {
      for (i=0;i<len;i++) b << (UShort_t)fValue[i];
   } else {
      b.WriteFastArray(fValue,len);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns name of leaf type.

const char *TLeafS::GetTypeName() const
{
   if (fIsUnsigned) return "UShort_t";
   return "Short_t";
}

////////////////////////////////////////////////////////////////////////////////
/// Returns current value of leaf.
/// - if leaf is a simple type, i must be set to 0
/// - if leaf is an array, i is the array element number to be returned

Double_t TLeafS::GetValue(Int_t i) const
{
   if (fIsUnsigned) return (UShort_t)fValue[i];
   return fValue[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Copy/set fMinimum and fMaximum to include/be wide than those of the parameter

Bool_t TLeafS::IncludeRange(TLeaf *input)
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

void TLeafS::Import(TClonesArray *list, Int_t n)
{
   const Short_t kShortUndefined = -9999;
   Int_t j = 0;
   char *clone;
   for (Int_t i=0;i<n;i++) {
      clone = (char*)list->UncheckedAt(i);
      if (clone) memcpy(&fValue[j],clone + fOffset, 2*fLen);
      else       memcpy(&fValue[j],&kShortUndefined,  2*fLen);
      j += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Prints leaf value.

void TLeafS::PrintValue(Int_t l) const
{
   if (fIsUnsigned) {
      UShort_t *uvalue = (UShort_t*)GetValuePointer();
      printf("%u",uvalue[l]);
   } else {
      Short_t *value = (Short_t*)GetValuePointer();
      printf("%d",value[l]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read leaf elements from Basket input buffer.

void TLeafS::ReadBasket(TBuffer &b)
{
   if (!fLeafCount && fNdata == 1) {
      b.ReadShort(fValue[0]);
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

////////////////////////////////////////////////////////////////////////////////
/// Deserialize input by performing byteswap as needed.
bool TLeafS::ReadBasketFast(TBuffer& input_buf, Long64_t N)
{
   if (R__unlikely(fLeafCount)) {return false;}
   return input_buf.ByteSwapBuffer(fLen*N, kShort_t);
}

////////////////////////////////////////////////////////////////////////////////
/// Read leaf elements from Basket input buffer and export buffer to
/// TClonesArray objects.

void TLeafS::ReadBasketExport(TBuffer &b, TClonesArray *list, Int_t n)
{
   if (n*fLen == 1) {
      b >> fValue[0];
   } else {
      b.ReadFastArray(fValue,n*fLen);
   }

   Int_t j = 0;
   for (Int_t i=0;i<n;i++) {
      memcpy((char*)list->UncheckedAt(i) + fOffset,&fValue[j], 2*fLen);
      j += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read a integer integer from std::istream s and store it into the branch buffer.

void TLeafS::ReadValue(std::istream &s, Char_t /*delim = ' '*/)
{
   if (fIsUnsigned) {
      UShort_t *uvalue = (UShort_t*)GetValuePointer();
      for (Int_t i=0;i<fLen;i++) s >> uvalue[i];
   } else {
      Short_t *value = (Short_t*)GetValuePointer();
      for (Int_t i=0;i<fLen;i++) s >> value[i];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set leaf buffer data address.

void TLeafS::SetAddress(void *add)
{
   if (ResetAddress(add) && (add!=fValue)) {
      delete [] fValue;
   }
   if (add) {
      if (TestBit(kIndirectAddress)) {
         fPointer = (Short_t**) add;
         Int_t ncountmax = fLen;
         if (fLeafCount) ncountmax = fLen*(fLeafCount->GetMaximum() + 1);
         if ((fLeafCount && ncountmax > Int_t(fLeafCount->GetValue())) ||
             ncountmax > fNdata || *fPointer == 0) {
            if (*fPointer) delete [] *fPointer;
            if (ncountmax > fNdata) fNdata = ncountmax;
            *fPointer = new Short_t[fNdata];
         }
         fValue = *fPointer;
      } else {
         fValue = (Short_t*)add;
      }
   } else {
      fValue = new Short_t[fNdata];
      fValue[0] = 0;
   }
}
