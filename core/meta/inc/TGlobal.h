// @(#)root/meta:$Id$
// Author: Rene Brun   13/11/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGlobal
#define ROOT_TGlobal


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGlobal                                                              //
//                                                                      //
// Global variables class (global variables are obtained from CINT).    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TDictionary.h"

#include <functional>


class TGlobal : public TDictionary {

private:
   DataMemberInfo_t  *fInfo;      //!pointer to CINT data member info

public:
   TGlobal(DataMemberInfo_t *info = 0);
   TGlobal (const TGlobal &);
   TGlobal &operator=(const TGlobal &);

   virtual       ~TGlobal();
   virtual Int_t  GetArrayDim() const;
   virtual DeclId_t GetDeclId() const;
   virtual Int_t  GetMaxIndex(Int_t dim) const;
   virtual void  *GetAddress() const;
   virtual const char *GetTypeName() const;
   virtual const char *GetFullTypeName() const;
   virtual Bool_t IsValid();
   virtual Long_t Property() const;
   virtual bool   Update(DataMemberInfo_t *info);

   ClassDef(TGlobal,2)  //Global variable class
};

   //Class to map the "funcky" globals and be able to add them to the list of globals.
   class TGlobalMappedFunction: public TGlobal {
   public:
      typedef void* (*GlobalFunc_t)();
      typedef std::function<void*()> GlobalFunctor_t;
      TGlobalMappedFunction(const char* name, const char* type,
                            GlobalFunc_t funcPtr):fFuncPtr(funcPtr)
      {
         SetNameTitle(name, type);
      }
      TGlobalMappedFunction(const char* name, const char* type,
                            void *dummyPtr, GlobalFunctor_t functor): fFuncPtr((GlobalFunc_t)dummyPtr), fFunctor(functor)
      {
         SetNameTitle(name, type);
      }
      virtual       ~TGlobalMappedFunction() {}
      Int_t          GetArrayDim() const { return 0;}
      DeclId_t       GetDeclId() const { return (DeclId_t)(fFuncPtr); } // Used as DeclId because of uniqueness
      Int_t          GetMaxIndex(Int_t /*dim*/) const { return -1; }
      void          *GetAddress() const { return !fFunctor ? (*fFuncPtr)() : fFunctor(); }
      const char    *GetTypeName() const { return fTitle; }
      const char    *GetFullTypeName() const { return fTitle; }
      Long_t         Property() const { return 0; }
      virtual bool   Update(DataMemberInfo_t * /*info*/) { return false; }
      static  void   Add(TGlobalMappedFunction* gmf);

   private:

      GlobalFunc_t fFuncPtr; // Function to call to get the address

      GlobalFunctor_t fFunctor; // functor which correctly returns pointer

      TGlobalMappedFunction &operator=(const TGlobal &); // not implemented.
      // Some of the special ones are created before the list is create e.g gFile
      // We need to buffer them.
      static TList&  GetEarlyRegisteredGlobals();

      friend class TROOT;
   };

#endif
