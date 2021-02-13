// @(#)root/cont:$Id$
// Author: Rene Brun   11/02/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TClonesArray
#define ROOT_TClonesArray


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TClonesArray                                                         //
//                                                                      //
// An array of clone TObjects. The array expands automatically when     //
// adding elements (shrinking can be done explicitly).                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObjArray.h"

class TClass;

class TClonesArray : public TObjArray {

protected:
   TClass       *fClass;       //!Pointer to the class of the elements
   TObjArray    *fKeep;        //!Saved copies of pointers to objects

public:
   enum EStatusBits {
      kBypassStreamer = BIT(12),  // Class Streamer not called (default)
      kForgetBits     = BIT(15)   // Do not create branches for fBits, fUniqueID
   };

   TClonesArray();
   TClonesArray(const char *classname, Int_t size = 1000, Bool_t call_dtor = kFALSE);
   TClonesArray(const TClass *cl, Int_t size = 1000, Bool_t call_dtor = kFALSE);
   TClonesArray(const TClonesArray& tc);
   virtual         ~TClonesArray();
   TClonesArray& operator=(const TClonesArray& tc);
   virtual void     Compress();
   virtual void     Clear(Option_t *option="");
   virtual void     Delete(Option_t *option="");
   virtual void     Expand(Int_t newSize);
   virtual void     ExpandCreate(Int_t n);
   virtual void     ExpandCreateFast(Int_t n);
   TClass          *GetClass() const { return fClass; }
   virtual void     SetOwner(Bool_t enable = kTRUE);

   void             AddFirst(TObject *) { MayNotUse("AddFirst"); }
   void             AddLast(TObject *) { MayNotUse("AddLast"); }
   void             AddAt(TObject *, Int_t) { MayNotUse("AddAt"); }
   void             AddAtAndExpand(TObject *, Int_t) { MayNotUse("AddAtAndExpand"); }
   Int_t            AddAtFree(TObject *) { MayNotUse("AddAtFree"); return 0; }
   void             AddAfter(const TObject *, TObject *) { MayNotUse("AddAfter"); }
   void             AddBefore(const TObject *, TObject *) { MayNotUse("AddBefore"); }
   void             BypassStreamer(Bool_t bypass=kTRUE);
   Bool_t           CanBypassStreamer() const { return TestBit(kBypassStreamer); }
   TObject         *ConstructedAt(Int_t idx);
   TObject         *ConstructedAt(Int_t idx, Option_t *clear_options);
   void             SetClass(const char *classname,Int_t size=1000);
   void             SetClass(const TClass *cl,Int_t size=1000);

   void             AbsorbObjects(TClonesArray *tc);
   void             AbsorbObjects(TClonesArray *tc, Int_t idx1, Int_t idx2);
   void             MultiSort(Int_t nTCs, TClonesArray** tcs, Int_t upto = kMaxInt);
   virtual TObject *RemoveAt(Int_t idx);
   virtual TObject *Remove(TObject *obj);
   virtual void     RemoveRange(Int_t idx1, Int_t idx2);
   virtual void     Sort(Int_t upto = kMaxInt);

   TObject         *New(Int_t idx);
   TObject         *AddrAt(Int_t idx);
   TObject         *&operator[](Int_t idx);
   TObject         *operator[](Int_t idx) const;

   ClassDef(TClonesArray,4)  //An array of clone objects
};

inline TObject *TClonesArray::AddrAt(Int_t idx)
{
   return operator[](idx);
}

#endif
