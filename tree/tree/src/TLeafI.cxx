// @(#)root/tree:$Id$
// Author: Rene Brun   12/01/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TLeafI
\ingroup tree

A TLeaf for an Integer data type.
*/

#include "TLeafI.h"
#include "TBranch.h"
#include "TBuffer.h"
#include "TClonesArray.h"
#include <iostream>

ClassImp(TLeafI);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor for LeafI.

TLeafI::TLeafI(): TLeaf()
{
   fLenType = 4;
   fMinimum = 0;
   fMaximum = 0;
   fValue   = 0;
   fPointer = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a LeafI.

TLeafI::TLeafI(TBranch *parent, const char *name, const char *type)
   :TLeaf(parent, name,type)
{
   fLenType = 4;
   fMinimum = 0;
   fMaximum = 0;
   fValue   = 0;
   fPointer = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Default destructor for a LeafI.

TLeafI::~TLeafI()
{
   if (ResetAddress(0,kTRUE)) delete [] fValue;
}

////////////////////////////////////////////////////////////////////////////////
/// Export element from local leaf buffer to ClonesArray.

void TLeafI::Export(TClonesArray *list, Int_t n)
{
   Int_t *value = fValue;
   for (Int_t i=0;i<n;i++) {
      char *first = (char*)list->UncheckedAt(i);
      Int_t *ii = (Int_t*)&first[fOffset];
      for (Int_t j=0;j<fLen;j++) {
         ii[j] = value[j];
      }
      value += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Pack leaf elements in Basket output buffer.

void TLeafI::FillBasket(TBuffer &b)
{
   Int_t i;
   Int_t len = GetLen();
   if (fPointer) fValue = *fPointer;
   if (IsRange()) {
      if (fValue[0] > fMaximum) fMaximum = fValue[0];
   }
   if (IsUnsigned()) {
      for (i=0;i<len;i++) b << (UInt_t)fValue[i];
   } else {
      b.WriteFastArray(fValue,len);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns name of leaf type.

const char *TLeafI::GetTypeName() const
{
   if (fIsUnsigned) return "UInt_t";
   return "Int_t";
}

////////////////////////////////////////////////////////////////////////////////
/// Returns current value of leaf
/// - if leaf is a simple type, i must be set to 0
/// - if leaf is an array, i is the array element number to be returned

Double_t TLeafI::GetValue(Int_t i) const
{
   if (fIsUnsigned) return (UInt_t)fValue[i];
   return fValue[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Copy/set fMinimum and fMaximum to include/be wide than those of the parameter

Bool_t TLeafI::IncludeRange(TLeaf *input)
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

void TLeafI::Import(TClonesArray *list, Int_t n)
{
   const Int_t kIntUndefined = -9999;
   Int_t j = 0;
   char *clone;
   for (Int_t i=0;i<n;i++) {
      clone = (char*)list->UncheckedAt(i);
      if (clone) memcpy(&fValue[j],clone + fOffset, 4*fLen);
      else       memcpy(&fValue[j],&kIntUndefined,  4*fLen);
      j += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Prints leaf value.

void TLeafI::PrintValue(Int_t l) const
{
   if (fIsUnsigned) {
      UInt_t *uvalue = (UInt_t*)GetValuePointer();
      printf("%u",uvalue[l]);
   } else {
      Int_t *value = (Int_t*)GetValuePointer();
      printf("%d",value[l]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read leaf elements from Basket input buffer.

void TLeafI::ReadBasket(TBuffer &b)
{
   if (!fLeafCount && fNdata == 1) {
      b.ReadInt(fValue[0]);
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
bool TLeafI::ReadBasketFast(TBuffer& input_buf, Long64_t N)
{
   if (R__unlikely(fLeafCount)) {return false;}
   return input_buf.ByteSwapBuffer(fLen*N, kInt_t);
}

////////////////////////////////////////////////////////////////////////////////
/// Read leaf elements from Basket input buffer and export buffer to
/// TClonesArray objects.

void TLeafI::ReadBasketExport(TBuffer &b, TClonesArray *list, Int_t n)
{
   if (n*fLen == 1) {
      b >> fValue[0];
   } else {
      b.ReadFastArray(fValue,n*fLen);
   }
   Int_t *value = fValue;
   for (Int_t i=0;i<n;i++) {
      char *first = (char*)list->UncheckedAt(i);
      Int_t *ii = (Int_t*)&first[fOffset];
      for (Int_t j=0;j<fLen;j++) {
         ii[j] = value[j];
      }
      value += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read an integer from std::istream s and store it into the branch buffer.

void TLeafI::ReadValue(std::istream &s, Char_t /*delim = ' '*/)
{
   if (fIsUnsigned) {
      UInt_t *uvalue = (UInt_t*)GetValuePointer();
      for (Int_t i=0;i<fLen;i++) s >> uvalue[i];
   } else {
      Int_t *value = (Int_t*)GetValuePointer();
      for (Int_t i=0;i<fLen;i++) s >> value[i];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set leaf buffer data address.

void TLeafI::SetAddress(void *add)
{
   if (ResetAddress(add) && (add!= fValue)) {
      delete [] fValue;
   }
   if (add) {
      if (TestBit(kIndirectAddress)) {
         fPointer = (Int_t**) add;
         Int_t ncountmax = fLen;
         if (fLeafCount) ncountmax = fLen*(fLeafCount->GetMaximum() + 1);
         if ((fLeafCount && ncountmax > Int_t(fLeafCount->GetValue())) ||
             ncountmax > fNdata || *fPointer == 0) {
            if (*fPointer) delete [] *fPointer;
            if (ncountmax > fNdata) fNdata = ncountmax;
            *fPointer = new Int_t[fNdata];
         }
         fValue = *fPointer;
      } else {
         fValue = (Int_t*)add;
      }
   } else {
      fValue = new Int_t[fNdata];
      fValue[0] = 0;
   }
}

