// @(#)root/tree:$Id$
// Author: Rene Brun   17/03/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TLeafC
\ingroup tree

A TLeaf for a variable length string.
*/

#include "TLeafC.h"
#include "TBuffer.h"
#include "TBranch.h"
#include "TBasket.h"
#include "TClonesArray.h"
#include "strlcpy.h"

#include <string>
#include <iostream>

ClassImp(TLeafC);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor for LeafC.

TLeafC::TLeafC(): TLeaf()
{
   fLenType = 1;
   fMinimum = 0;
   fMaximum = 0;
   fValue   = 0;
   fPointer = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a LeafC.

TLeafC::TLeafC(TBranch *parent, const char *name, const char *type)
   :TLeaf(parent, name,type)
{
   fLenType = 1;
   fMinimum = 0;
   fMaximum = 0;
   fValue   = 0;
   fPointer = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Default destructor for a LeafC.

TLeafC::~TLeafC()
{
   if (ResetAddress(0,kTRUE)) delete [] fValue;
}

////////////////////////////////////////////////////////////////////////////////
/// Export element from local leaf buffer to ClonesArray.

void TLeafC::Export(TClonesArray *list, Int_t n)
{
   Int_t j = 0;
   for (Int_t i=0;i<n;i++) {
      memcpy((char*)list->UncheckedAt(i) + fOffset,&fValue[j], 1);
      j += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Pack leaf elements in Basket output buffer.

void TLeafC::FillBasket(TBuffer &b)
{
   if (fPointer) fValue = *fPointer;
   Int_t len = strlen(fValue);
   if (len >= fMaximum) fMaximum = len+1;
   if (len >= fLen)     fLen = len+1;
   b.WriteFastArrayString(fValue,len);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns name of leaf type.

const char *TLeafC::GetTypeName() const
{
   if (fIsUnsigned) return "UChar_t";
   return "Char_t";
}

////////////////////////////////////////////////////////////////////////////////
/// Copy/set fMinimum and fMaximum to include/be wide than those of the parameter

Bool_t TLeafC::IncludeRange(TLeaf *input)
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

void TLeafC::Import(TClonesArray *list, Int_t n)
{
   Int_t j = 0;
   for (Int_t i=0;i<n;i++) {
      memcpy(&fValue[j],(char*)list->UncheckedAt(i) + fOffset, 1);
      j += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Prints leaf value.

void TLeafC::PrintValue(Int_t) const
{
   char *value = (char*)GetValuePointer();
   printf("%s",value);
}

////////////////////////////////////////////////////////////////////////////////
/// Read leaf elements from Basket input buffer.

void TLeafC::ReadBasket(TBuffer &b)
{
   // Try to deal with the file written during the time where len was not
   // written to disk when len was == 0.
   Int_t readbasket = GetBranch()->GetReadBasket();
   TBasket *basket = GetBranch()->GetBasket(readbasket);
   if (!basket) {
      fValue[0] = '\0';
      return;
   }
   Int_t* entryOffset = basket->GetEntryOffset();
   if (entryOffset) {
      Long64_t first = GetBranch()->GetBasketEntry()[readbasket];
      Long64_t entry = GetBranch()->GetReadEntry();
      if ( (readbasket == GetBranch()->GetWriteBasket() && (entry+1) == GetBranch()->GetEntries()) /* Very last entry */
               ||
               (readbasket <  GetBranch()->GetWriteBasket() && (entry+1) == GetBranch()->GetBasketEntry()[readbasket+1] ) /* Last entry of the basket */
           )
         {
            if ( entryOffset[entry-first] == basket->GetLast() ) /* The 'read' point is at the end of the basket */
               {
                  // Empty string
                  fValue[0] = '\0';
                  return;
               }
         }
      else if ( entryOffset[entry-first] == entryOffset[entry-first+1] ) /* This string did not use up any space in the buffer */
         {
            // Empty string
            fValue[0] = '\0';
            return;
         }
   }
   b.ReadFastArrayString(fValue,fLen);
}

////////////////////////////////////////////////////////////////////////////////
/// Read leaf elements from Basket input buffer and export buffer to
/// TClonesArray objects.

void TLeafC::ReadBasketExport(TBuffer &b, TClonesArray *list, Int_t n)
{
   UChar_t len;
   b >> len;
   if (len) {
      if (len >= fLen) len = fLen-1;
      b.ReadFastArray(fValue,len);
      fValue[len] = 0;
   } else {
      fValue[0] = 0;
   }

   Int_t j = 0;
   for (Int_t i=0;i<n;i++) {
      memcpy((char*)list->UncheckedAt(i) + fOffset,&fValue[j], 1);
      j += fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read a string from std::istream s up to delimiter and store it into the branch
/// buffer.

void TLeafC::ReadValue(std::istream &s, Char_t delim /*= ' '*/)
{
   std::string temp;
   std::getline(s, temp, delim);
   if (TestBit(kNewValue) &&
        (temp.length()+1 > ((UInt_t)fNdata))) {
      // Grow buffer if needed and we created the buffer.
      fNdata = ((UInt_t)temp.size()) + 1;
      if (TestBit(kIndirectAddress) && fPointer) {
         delete [] *fPointer;
         *fPointer = new char[fNdata];
      } else {
         fValue = new char[fNdata];
      }
   }
   strlcpy(fValue,temp.c_str(),fNdata);
}

////////////////////////////////////////////////////////////////////////////////
/// Set leaf buffer data address.

void TLeafC::SetAddress(void *add)
{
   if (ResetAddress(add)) {
      delete [] fValue;
   }
   if (add) {
      if (TestBit(kIndirectAddress)) {
         fPointer = (char**)add;
         Int_t ncountmax = fLen;
         if (fLeafCount) ncountmax = fLen*(fLeafCount->GetMaximum() + 1);
         if ((fLeafCount && ncountmax > Int_t(fLeafCount->GetValue())) ||
             ncountmax > fNdata || *fPointer == 0) {
            if (*fPointer) delete [] *fPointer;
            if (ncountmax > fNdata) fNdata = ncountmax;
            *fPointer = new char[fNdata];
         }
         fValue = *fPointer;
      } else {
         fValue = (char*)add;
      }
   }
   else {
      fValue = new char[fNdata];
      fValue[0] = 0;
   }
}
