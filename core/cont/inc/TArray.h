// @(#)root/cont:$Id$
// Author: Fons Rademakers   21/10/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TArray
#define ROOT_TArray


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TArray                                                               //
//                                                                      //
// Abstract array base class. Used by TArrayC, TArrayS, TArrayI,        //
// TArrayL, TArrayF and TArrayD.                                        //
// Data member is public for historical reasons.                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"
#include <cstring>   // used in all TArray derived classes for memset()

class TBuffer;

class TArray {

protected:
   Bool_t        BoundsOk(const char *where, Int_t at) const;
   Bool_t        OutOfBoundsError(const char *where, Int_t i) const;

public:
   Int_t     fN;            //Number of array elements

   TArray(): fN(0) { }
   TArray(Int_t n): fN(n) { }
   TArray(const TArray &a): fN(a.fN) { }
   TArray         &operator=(const TArray &rhs)
     {if(this!=&rhs) fN = rhs.fN; return *this; }
   virtual        ~TArray() { fN = 0; }

   Int_t          GetSize() const { return fN; }
   virtual void   Set(Int_t n) = 0;

   virtual Double_t GetAt(Int_t i) const = 0;
   virtual void   SetAt(Double_t v, Int_t i) = 0;

   static TArray *ReadArray(TBuffer &b, const TClass *clReq);
   static void    WriteArray(TBuffer &b, const TArray *a);

   friend TBuffer &operator<<(TBuffer &b, const TArray *obj);

   ClassDef(TArray,1)  //Abstract array base class
};

#if defined R__TEMPLATE_OVERLOAD_BUG
template <>
#endif
inline TBuffer &operator>>(TBuffer &buf, TArray *&obj)
{
   // Read TArray object from buffer.

   obj = (TArray *) TArray::ReadArray(buf, TArray::Class());
   return buf;
}

#if defined R__TEMPLATE_OVERLOAD_BUG
template <>
#endif
TBuffer &operator<<(TBuffer &b, const TArray *obj);

inline Bool_t TArray::BoundsOk(const char *where, Int_t at) const
{
   return (at < 0 || at >= fN)
                  ? OutOfBoundsError(where, at)
                  : kTRUE;
}

#endif
