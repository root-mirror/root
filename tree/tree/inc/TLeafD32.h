// @(#)root/tree:$Id$
// Author: Simon Spies 23/02/19

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TLeafD32
#define ROOT_TLeafD32

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TLeafD32                                                             //
//                                                                      //
// A TLeaf for a 24 bit truncated floating point data type.             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TLeaf.h"

class TStreamerElement;


class TLeafD32 : public TLeaf {

protected:
   Double32_t fMinimum;        ///<  Minimum value if leaf range is specified
   Double32_t fMaximum;        ///<  Maximum value if leaf range is specified
   Double32_t *fValue;         ///<! Pointer to data buffer
   Double32_t **fPointer;      ///<! Address of pointer to data buffer
   TStreamerElement *fElement; ///<! StreamerElement used for TBuffer read / write

public:
   TLeafD32();
   TLeafD32(TBranch *parent, const char *name, const char *type);
   virtual ~TLeafD32();

   virtual void Export(TClonesArray *list, Int_t n);
   virtual void FillBasket(TBuffer &b);
   const char *GetTypeName() const { return "Double32_t"; }
   Double_t GetValue(Int_t i = 0) const;
   virtual void *GetValuePointer() const { return fValue; }
   virtual void Import(TClonesArray *list, Int_t n);
   virtual void PrintValue(Int_t i = 0) const;
   virtual void ReadBasket(TBuffer &b);
   virtual void ReadBasketExport(TBuffer &b, TClonesArray *list, Int_t n);
   virtual void ReadValue(std::istream &s, Char_t delim = ' ');
   virtual void SetAddress(void *add = 0);

   ClassDef(TLeafD32, 1); // A TLeaf for a 24 bit truncated floating point data type.
};

// if leaf is a simple type, i must be set to 0
// if leaf is an array, i is the array element number to be returned
inline Double_t TLeafD32::GetValue(Int_t i) const
{
   return fValue[i];
}

#endif
