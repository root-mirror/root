// @(#)root/meta:$Id$
// Author: Fons Rademakers   08/02/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TBaseClass
#define ROOT_TBaseClass


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TBaseClass                                                           //
//                                                                      //
// Description of a base class.                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TDictionary.h"
#include "TClassRef.h"

#include <atomic>

class TBrowser;
class TClass;

class TBaseClass : public TDictionary {
#ifndef __CLING__
   using AtomicInt_t = std::atomic<Int_t>;
   static_assert(sizeof(std::atomic<Int_t>) == sizeof(Int_t),
                 "We requiqre atomic<int> and <int> to have the same size but they are not");
#else
   // std::atomic is not yet supported in the I/O, so
   // we hide them from Cling
   using AtomicInt_t = Int_t;
#endif

private:
   TBaseClass(const TBaseClass &) = delete;
   TBaseClass&operator=(const TBaseClass &) = delete;

private:
   BaseClassInfo_t    *fInfo;      //!pointer to CINT base class info
   TClassRef           fClassPtr;  // pointer to the base class TClass
   TClass             *fClass;     //!pointer to parent class
   AtomicInt_t         fDelta;     // BaseClassInfo_t offset (INT_MAX if unset)
   mutable AtomicInt_t fProperty;  // BaseClassInfo_t's properties
   Int_t               fSTLType;   // cache of IsSTLContainer()

public:
   TBaseClass(BaseClassInfo_t *info = nullptr, TClass *cl = nullptr);
   virtual     ~TBaseClass();
   virtual void   Browse(TBrowser *b);
   const char    *GetTitle() const;
   TClass        *GetClassPointer(Bool_t load=kTRUE);
   Int_t          GetDelta();
   Bool_t         IsFolder() const {return kTRUE;}
   ROOT::ESTLType IsSTLContainer();
   Long_t         Property() const;
   void           SetClass(TClass* cl) { fClass = cl; }

   ClassDef(TBaseClass,2)  //Description of a base class
};

#endif
