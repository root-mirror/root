// @(#)root/io:$Id$
// Author: Markus Frank 28/10/04

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGenCollectionProxy.h"
#include "TVirtualStreamerInfo.h"
#include "TStreamerElement.h"
#include "TClassEdit.h"
#include "TClass.h"
#include "TError.h"
#include "TEnum.h"
#include "TROOT.h"
#include "TInterpreter.h" // For gInterpreterMutex
#include "TVirtualMutex.h"
#include "TStreamerInfoActions.h"
#include "THashTable.h"
#include "THashList.h"
#include <cstdlib>

#define MESSAGE(which,text)

/**
\class TGenVectorProxy
\ingroup IO
Local optimization class.

Collection proxies get copied. On copy we switch the type of the
proxy to the concrete STL type. The concrete types are optimized
for element access.
*/

class TGenVectorProxy : public TGenCollectionProxy {
public:
   // Standard Destructor
   TGenVectorProxy(const TGenCollectionProxy& c) : TGenCollectionProxy(c)
   {
   }
   // Standard Destructor
   virtual ~TGenVectorProxy()
{
   }
   // Return the address of the value at index 'idx'
   virtual void* At(UInt_t idx)
{
      if ( fEnv && fEnv->fObject ) {
         fEnv->fIdx = idx;
         switch( idx ) {
         case 0:
            return fEnv->fStart = fFirst.invoke(fEnv);
         default:
            if (! fEnv->fStart ) fEnv->fStart = fFirst.invoke(fEnv);
            return ((char*)fEnv->fStart) + fValDiff*idx;
         }
      }
      Fatal("TGenVectorProxy","At> Logic error - no proxy object set.");
      return 0;
   }
   // Call to delete/destruct individual item
   virtual void DeleteItem(Bool_t force, void* ptr) const
   {
      if ( force && ptr ) {
         if ( fVal->fProperties&kNeedDelete) {
            TVirtualCollectionProxy *proxy = fVal->fType->GetCollectionProxy();
            TPushPop helper(proxy,ptr);
            proxy->Clear("force");
         }
         fVal->DeleteItem(ptr);
      }
   }
};

/**
\class TGenVectorBoolProxy
\ingroup IO
Local optimization class.

Collection proxies get copied. On copy we switch the type of the
proxy to the concrete STL type. The concrete types are optimized
for element access.
*/
class TGenVectorBoolProxy : public TGenCollectionProxy {
   Bool_t fLastValue;

public:
   TGenVectorBoolProxy(const TGenCollectionProxy& c) : TGenCollectionProxy(c), fLastValue(false)
   {
      // Standard Constructor.
   }
   virtual ~TGenVectorBoolProxy()
   {
      // Standard Destructor.
   }
   virtual void* At(UInt_t idx)
   {
      // Return the address of the value at index 'idx'

      // However we can 'take' the address of the content of std::vector<bool>.
      if ( fEnv && fEnv->fObject ) {
         auto vec = (std::vector<bool> *)(fEnv->fObject);
         fLastValue = (*vec)[idx];
         fEnv->fIdx = idx;
         return &fLastValue;
      }
      Fatal("TGenVectorProxy","At> Logic error - no proxy object set.");
      return 0;
   }

   virtual void DeleteItem(Bool_t force, void* ptr) const
   {
      // Call to delete/destruct individual item
      if ( force && ptr ) {
         fVal->DeleteItem(ptr);
      }
   }
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  class TGenBitsetProxy
//
//   Local optimization class.
//
//   Collection proxies get copied. On copy we switch the type of the
//   proxy to the concrete STL type. The concrete types are optimized
//   for element access.
//
//////////////////////////////////////////////////////////////////////////
class TGenBitsetProxy : public TGenCollectionProxy {

public:
   TGenBitsetProxy(const TGenCollectionProxy& c) : TGenCollectionProxy(c)
   {
      // Standard Constructor.
   }
   virtual ~TGenBitsetProxy()
   {
      // Standard Destructor.
   }
   virtual void* At(UInt_t idx)
   {
      // Return the address of the value at index 'idx'

      if ( fEnv && fEnv->fObject ) {
         switch( idx ) {
            case 0:
               fEnv->fStart = fFirst.invoke(fEnv);
               fEnv->fIdx = idx;
               break;
            default:
               fEnv->fIdx = idx - fEnv->fIdx;
               if (! fEnv->fStart ) fEnv->fStart = fFirst.invoke(fEnv);
               fNext.invoke(fEnv);
               fEnv->fIdx = idx;
               break;
         }
         typedef ROOT::TCollectionProxyInfo::Environ<std::pair<size_t,Bool_t> > EnvType_t;
         EnvType_t *e = (EnvType_t*)fEnv;
         return &(e->fIterator.second);
      }
      Fatal("TGenVectorProxy","At> Logic error - no proxy object set.");
      return 0;
   }

   virtual void DeleteItem(Bool_t force, void* ptr) const
   {
      // Call to delete/destruct individual item
      if ( force && ptr ) {
         fVal->DeleteItem(ptr);
      }
   }
};

/*
\class TGenListProxy
\ingroup IO
Local optimization class.

Collection proxies get copied. On copy we switch the type of the
proxy to the concrete STL type. The concrete types are optimized
for element access.
**/

class TGenListProxy : public TGenVectorProxy {
public:
   // Standard Destructor
   TGenListProxy(const TGenCollectionProxy& c) : TGenVectorProxy(c)
{
   }
   // Standard Destructor
   virtual ~TGenListProxy()
{
   }
   // Return the address of the value at index 'idx'
   void* At(UInt_t idx)
{
      if ( fEnv && fEnv->fObject ) {
         switch( idx ) {
         case 0:
            fEnv->fIdx = idx;
            return fEnv->fStart = fFirst.invoke(fEnv);
         default:  {
            fEnv->fIdx = idx - fEnv->fIdx;
            if (! fEnv->fStart ) fEnv->fStart = fFirst.invoke(fEnv);
            void* result = fNext.invoke(fEnv);
            fEnv->fIdx = idx;
            return result;
         }
         }
      }
      Fatal("TGenListProxy","At> Logic error - no proxy object set.");
      return 0;
   }
};

/**
\class TGenSetProxy
\ingroup IO
Localoptimization class.

Collection proxies get copied. On copy we switch the type of the
proxy to the concrete STL type. The concrete types are optimized
for element access.
*/

class TGenSetProxy : public TGenVectorProxy {
public:
   // Standard Destructor
   TGenSetProxy(const TGenCollectionProxy& c) : TGenVectorProxy(c)
{
   }
   // Standard Destructor
   virtual ~TGenSetProxy()
{
   }
   // Return the address of the value at index 'idx'
   void* At(UInt_t idx)
{
      if ( fEnv && fEnv->fObject ) {
         if ( fEnv->fUseTemp ) {
            return (((char*)fEnv->fTemp)+idx*fValDiff);
         }
         switch( idx ) {
         case 0:
            fEnv->fIdx = idx;
            return fEnv->fStart = fFirst.invoke(fEnv);
         default:  {
            fEnv->fIdx = idx - fEnv->fIdx;
            if (! fEnv->fStart ) fEnv->fStart = fFirst.invoke(fEnv);
            void* result = fNext.invoke(fEnv);
            fEnv->fIdx = idx;
            return result;
         }
         }
      }
      Fatal("TGenSetProxy","At> Logic error - no proxy object set.");
      return 0;
   }
};

/**
\class TGenMapProxy
\ingroup IO
Localoptimization class.

Collection proxies get copied. On copy we switch the type of the
proxy to the concrete STL type. The concrete types are optimized
for element access.
*/

class TGenMapProxy : public TGenSetProxy {
public:
   // Standard Destructor
   TGenMapProxy(const TGenCollectionProxy& c) : TGenSetProxy(c)
{
   }
   // Standard Destructor
   virtual ~TGenMapProxy()
{
   }
   // Call to delete/destruct individual item
   virtual void DeleteItem(Bool_t force, void* ptr) const
   {
      if (force) {
         if ( fKey->fProperties&kNeedDelete) {
            TVirtualCollectionProxy *proxy = fKey->fType->GetCollectionProxy();
            TPushPop helper(proxy,fKey->fCase&kIsPointer ? *(void**)ptr : ptr);
            proxy->Clear("force");
         }
         if ( fVal->fProperties&kNeedDelete) {
            TVirtualCollectionProxy *proxy = fVal->fType->GetCollectionProxy();
            char *addr = ((char*)ptr)+fValOffset;
            TPushPop helper(proxy,fVal->fCase&kIsPointer ? *(void**)addr : addr);
            proxy->Clear("force");
         }
      }
      if ( fKey->fCase&kIsPointer ) {
         fKey->DeleteItem(*(void**)ptr);
      }
      if ( fVal->fCase&kIsPointer ) {
         char *addr = ((char*)ptr)+fValOffset;
         fVal->DeleteItem(*(void**)addr);
      }
   }
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGenCollectionProxy::Value::Value(const std::string& inside_type, Bool_t silent)
{
   std::string inside = (inside_type.find("const ")==0) ? inside_type.substr(6) : inside_type;
   fCase = 0;
   fProperties = 0;
   fCtor = 0;
   fDtor = 0;
   fDelete = 0;
   fSize = std::string::npos;
   fKind = kNoType_t;

   // Let's treat the unique_ptr case
   bool nameChanged = false;
   std::string intype = TClassEdit::GetNameForIO(inside.c_str(), TClassEdit::EModType::kNone, &nameChanged);

   bool isPointer = nameChanged; // unique_ptr is considered a pointer
   // The incoming name is normalized (it comes from splitting the name of a TClass),
   // so all we need to do is drop the last trailing star (if any) and record that information.
   if (!nameChanged && intype[intype.length()-1] == '*') {
      isPointer = true;
      intype.pop_back();
      if (intype[intype.length()-1] == '*') {
         // The value is a pointer to a pointer
         if (!silent)
            Warning("TGenCollectionProxy::Value::Value", "I/O not supported for collection of pointer to pointer: %s", inside_type.c_str());
         fSize = sizeof(void*);
         fKind = kVoid_t;
         return;
      }
   }

   if ( intype.substr(0,6) == "string" || intype.substr(0,11) == "std::string" ) {
      fCase = kBIT_ISSTRING;
      fType = TClass::GetClass("string");
      fCtor = fType->GetNew();
      fDtor = fType->GetDestructor();
      fDelete = fType->GetDelete();
      if (isPointer) {
         fCase |= kIsPointer;
         fSize = sizeof(void*);
      } else {
         fSize = sizeof(std::string);
      }
   }
   else {
      // In the case where we have an emulated class,
      // if the class is nested (in a class or a namespace),
      // calling G__TypeInfo ti(inside.c_str());
      // might fail because CINT does not known the nesting
      // scope, so let's first look for an emulated class:

      fType = TClass::GetClass(intype.c_str(),kTRUE,silent);

      if (fType) {
         if (isPointer) {
            fCase |= kIsPointer;
            fSize = sizeof(void*);
            if (fType == TString::Class()) {
               fCase |= kBIT_ISTSTRING;
            }
         }
         fCase  |= kIsClass;
         fCtor   = fType->GetNew();
         fDtor   = fType->GetDestructor();
         fDelete = fType->GetDelete();
      } else {
         R__LOCKGUARD(gInterpreterMutex);

         // Try to avoid autoparsing.

         THashTable *typeTable = dynamic_cast<THashTable*>( gROOT->GetListOfTypes() );
         assert(typeTable && "The type of the list of type has changed");

         TDataType *fundType = (TDataType *)typeTable->THashTable::FindObject( intype.c_str() );
         if (fundType && fundType->GetType() < 0x17 && fundType->GetType() > 0) {
            fKind = (EDataType)fundType->GetType();
            // R__ASSERT((fKind>0 && fKind<0x17) || (fKind==-1&&(prop&kIsPointer)) );

            fCase |= kIsFundamental;
            if (isPointer) {
               fCase |= kIsPointer;
               fSize = sizeof(void*);
            } else {
               fSize = fundType->Size();
            }
         } else if (TEnum::GetEnum( intype.c_str(), TEnum::kNone) ) {
            // This is a known enum.
            fCase = kIsEnum;
            fSize = sizeof(Int_t);
            fKind = kInt_t;
            if (isPointer) {
               fCase |= kIsPointer;
               fSize = sizeof(void*);
            }
         } else {
            // This fallback solution should be hardly used ...
            // One of the common use case is to 'discover' that this is a
            // collection for the content of which we do not have a dictionary
            // which can happen at least in the following cases:
            //    - empty emulated collection
            //    - emulated collection of enums
            // In those two cases there is no StreamerInfo stored in the file
            // for the content.

            // R__ASSERT("FallBack, should be hardly used.");

            TypeInfo_t *ti = gCling->TypeInfo_Factory();
            gCling->TypeInfo_Init(ti,inside.c_str());
            if ( !gCling->TypeInfo_IsValid(ti) ) {
               if (isPointer) {
                  fCase |= kIsPointer;
                  fSize = sizeof(void*);
               }
               // Since we already search for GetClass earlier, this should
               // never be true.
//               fType = TClass::GetClass(intype.c_str(),kTRUE,silent);
//               if (fType) {
//                  fCase  |= kIsClass;
//                  fCtor   = fType->GetNew();
//                  fDtor   = fType->GetDestructor();
//                  fDelete = fType->GetDelete();
//               }
//               else {
                  // either we have an Emulated enum or a really unknown class!
                  // let's just claim its an enum :(
                  fCase = kIsEnum;
                  fSize = sizeof(Int_t);
                  fKind = kInt_t;
//               }
            }
            else {
               Long_t prop = gCling->TypeInfo_Property(ti);
               if ( prop&kIsPointer ) {
                  fSize = sizeof(void*);
               }
               if ( prop&kIsStruct ) {
                  prop |= kIsClass;
               }

               if ( prop&kIsClass ) {
                  // We can get here in the case where the value if forward declared or
                  // is an std::pair that can not be (yet) emulated (eg. "std::pair<int,void*>")
                  fSize = std::string::npos;
                  if (!silent)
                     Error("TGenCollectionProxy", "Could not retrieve the TClass for %s", intype.c_str());
//                fType = TClass::GetClass(intype.c_str(),kTRUE,silent);
//                R__ASSERT(fType);
//                fCtor   = fType->GetNew();
//                fDtor   = fType->GetDestructor();
//                fDelete = fType->GetDelete();
               }
               else if ( prop&kIsFundamental ) {
                  fundType = gROOT->GetType( intype.c_str() );
                  if (fundType==0) {
                     if (intype != "long double" && !silent) {
                        Error("TGenCollectionProxy","Unknown fundamental type %s",intype.c_str());
                     }
                     fSize = sizeof(int);
                     fKind = kInt_t;
                  } else {
                     fKind = (EDataType)fundType->GetType();
                     fSize = gCling->TypeInfo_Size(ti);
                     R__ASSERT((fKind>0 && fKind<0x17) || (fKind==-1&&(prop&kIsPointer)) );
                  }
               }
               else if ( prop&kIsEnum ) {
                  fSize = sizeof(int);
                  fKind = kInt_t;
               }
               fCase = prop & (kIsPointer|kIsFundamental|kIsEnum|kIsClass);
               if (fType == TString::Class() && (fCase&kIsPointer)) {
                  fCase |= kBIT_ISTSTRING;
               }
            }
            gCling->TypeInfo_Delete(ti);
         }
      }
      if (fType) {
         TVirtualCollectionProxy *proxy = fType->GetCollectionProxy();
         if (proxy && (proxy->GetProperties() & kNeedDelete)) {
            fProperties |= kNeedDelete;
         }
      }
   }
   if ( fSize == std::string::npos ) {
      if ( fType == 0 ) {
         // The caller should check the validity by calling IsValid()
      } else {
         fSize = fType->Size();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the Value has been properly initialized.

Bool_t TGenCollectionProxy::Value::IsValid()
{


   return fSize != std::string::npos;
}

void TGenCollectionProxy::Value::DeleteItem(void* ptr)
{
   // Delete an item.

   if ( ptr && fCase&kIsPointer ) {
      if ( fDelete ) {
         (*fDelete)(ptr);
      }
      else if ( fType ) {
         fType->Destructor(ptr);
      }
      else {
         ::operator delete(ptr);
      }
   }
}

/**
 \class TGenCollectionProxy TGenCollectionProxy.cxx
 \ingroup IO

 Proxy around an arbitrary container, which implements basic
 functionality and iteration.

 The purpose of this implementation
 is to shield any generated dictionary implementation from the
 underlying streamer/proxy implementation and only expose
 the creation functions.

 In particular this is used to implement splitting and abstract
 element access of any container. Access to compiled code is necessary
 to implement the abstract iteration sequence and functionality like
 size(), clear(), resize(). resize() may be a void operation.
*/

////////////////////////////////////////////////////////////////////////////////
/// Build a proxy for an emulated container.

TGenCollectionProxy::TGenCollectionProxy(const TGenCollectionProxy& copy)
   : TVirtualCollectionProxy(copy.fClass),
     fTypeinfo(copy.fTypeinfo)
{
   fEnv            = 0;
   fName           = copy.fName;
   fPointers       = copy.fPointers;
   fSTL_type       = copy.fSTL_type;
   fSize.call      = copy.fSize.call;
   fNext.call      = copy.fNext.call;
   fFirst.call     = copy.fFirst.call;
   fClear.call     = copy.fClear.call;
   fResize         = copy.fResize;
   fDestruct       = copy.fDestruct;
   fConstruct      = copy.fConstruct;
   fFeed           = copy.fFeed;
   fCollect        = copy.fCollect;
   fCreateEnv.call = copy.fCreateEnv.call;
   fValOffset      = copy.fValOffset;
   fValDiff        = copy.fValDiff;
   fValue          = copy.fValue.load(std::memory_order_relaxed) ? new Value(*copy.fValue) : 0;
   fVal            = copy.fVal   ? new Value(*copy.fVal)   : 0;
   fKey            = copy.fKey   ? new Value(*copy.fKey)   : 0;
   fOnFileClass    = copy.fOnFileClass;
   fReadMemberWise = new TObjArray(TCollection::kInitCapacity,-1);
   fConversionReadMemberWise = 0;
   fWriteMemberWise = 0;
   fProperties     = copy.fProperties;
   fFunctionCreateIterators    = copy.fFunctionCreateIterators;
   fFunctionCopyIterator       = copy.fFunctionCopyIterator;
   fFunctionNextIterator       = copy.fFunctionNextIterator;
   fFunctionDeleteIterator     = copy.fFunctionDeleteIterator;
   fFunctionDeleteTwoIterators = copy.fFunctionDeleteTwoIterators;
}

////////////////////////////////////////////////////////////////////////////////
/// Build a proxy for a collection whose type is described by 'collectionClass'.

TGenCollectionProxy::TGenCollectionProxy(Info_t info, size_t iter_size)
   : TVirtualCollectionProxy(0),
     fTypeinfo(info)
{
   fEnv             = 0;
   fSize.call       = 0;
   fFirst.call      = 0;
   fNext.call       = 0;
   fClear.call      = 0;
   fResize          = 0;
   fDestruct        = 0;
   fConstruct       = 0;
   fCollect         = 0;
   fCreateEnv.call  = 0;
   fFeed            = 0;
   fValue           = 0;
   fKey             = 0;
   fVal             = 0;
   fValOffset       = 0;
   fValDiff         = 0;
   fPointers        = false;
   fOnFileClass     = 0;
   fSTL_type        = ROOT::kNotSTL;
   Env_t e;
   if ( iter_size > sizeof(e.fIterator) ) {
      Fatal("TGenCollectionProxy",
            "%s %s are too large:%ld bytes. Maximum is:%ld bytes",
            "Iterators for collection",
            fClass->GetName(),
            (Long_t)iter_size,
            (Long_t)sizeof(e.fIterator));
   }
   fReadMemberWise = new TObjArray(TCollection::kInitCapacity,-1);
   fConversionReadMemberWise   = 0;
   fWriteMemberWise            = 0;
   fFunctionCreateIterators    = 0;
   fFunctionCopyIterator       = 0;
   fFunctionNextIterator       = 0;
   fFunctionDeleteIterator     = 0;
   fFunctionDeleteTwoIterators = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Build a proxy for a collection whose type is described by 'collectionClass'.

TGenCollectionProxy::TGenCollectionProxy(const ROOT::TCollectionProxyInfo &info, TClass *cl)
   : TVirtualCollectionProxy(cl),
     fTypeinfo(info.fInfo), fOnFileClass(0)
{
   fEnv            = 0;
   fValDiff        = info.fValueDiff;
   fValOffset      = info.fValueOffset;
   fSize.call      = info.fSizeFunc;
   fResize         = info.fResizeFunc;
   fNext.call      = info.fNextFunc;
   fFirst.call     = info.fFirstFunc;
   fClear.call     = info.fClearFunc;
   fConstruct      = info.fConstructFunc;
   fDestruct       = info.fDestructFunc;
   fFeed           = info.fFeedFunc;
   fCollect        = info.fCollectFunc;
   fCreateEnv.call = info.fCreateEnv;

   if (cl) {
      fName = cl->GetName();
   }
   CheckFunctions();

   fValue           = 0;
   fKey             = 0;
   fVal             = 0;
   fPointers        = false;
   fSTL_type        = ROOT::kNotSTL;

   Env_t e;
   if ( info.fIterSize > sizeof(e.fIterator) ) {
      Fatal("TGenCollectionProxy",
            "%s %s are too large:%ld bytes. Maximum is:%ld bytes",
            "Iterators for collection",
            fClass->GetName(),
            (Long_t)info.fIterSize,
            (Long_t)sizeof(e.fIterator));
   }
   fReadMemberWise = new TObjArray(TCollection::kInitCapacity,-1);
   fConversionReadMemberWise   = 0;
   fWriteMemberWise            = 0;
   fFunctionCreateIterators    = info.fCreateIterators;
   fFunctionCopyIterator       = info.fCopyIterator;
   fFunctionNextIterator       = info.fNext;
   fFunctionDeleteIterator     = info.fDeleteSingleIterator;
   fFunctionDeleteTwoIterators = info.fDeleteTwoIterators;
}

namespace {
   template <class vec>
   void clearVector(vec& v)
   {
      // Clear out the proxies.

      for(typename vec::iterator i=v.begin(); i != v.end(); ++i) {
         typename vec::value_type e = *i;
         if ( e ) {
            delete e;
         }
      }
      v.clear();
   }
}
////////////////////////////////////////////////////////////////////////////////
/// Standard destructor

TGenCollectionProxy::~TGenCollectionProxy()
{
   clearVector(fProxyList);
   clearVector(fProxyKept);
   clearVector(fStaged);

   if ( fValue.load() ) delete fValue.load();
   if ( fVal   ) delete fVal;
   if ( fKey   ) delete fKey;

   delete fReadMemberWise;
   if (fConversionReadMemberWise) {
      std::map<std::string, TObjArray*>::iterator it;
      std::map<std::string, TObjArray*>::iterator end = fConversionReadMemberWise->end();
      for( it = fConversionReadMemberWise->begin(); it != end; ++it ) {
         delete it->second;
      }
      delete fConversionReadMemberWise;
      fConversionReadMemberWise = 0;
   }
   delete fWriteMemberWise;
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual copy constructor

TVirtualCollectionProxy* TGenCollectionProxy::Generate() const
{
   if ( !fValue.load(std::memory_order_relaxed) ) Initialize(kFALSE);

   if( fPointers )
      return new TGenCollectionProxy(*this);

   switch(fSTL_type) {
      case ROOT::kSTLbitset: {
         return new TGenBitsetProxy(*this);
      }
      case ROOT::kSTLvector: {
         if ((*fValue).fKind == kBool_t) {
            return new TGenVectorBoolProxy(*this);
         } else {
            return new TGenVectorProxy(*this);
         }
      }
      case ROOT::kSTLlist:
      case ROOT::kSTLforwardlist:
         return new TGenListProxy(*this);
      case ROOT::kSTLmap:
      case ROOT::kSTLunorderedmap:
      case ROOT::kSTLmultimap:
      case ROOT::kSTLunorderedmultimap:
         return new TGenMapProxy(*this);
      case ROOT::kSTLset:
      case ROOT::kSTLunorderedset:
      case ROOT::kSTLmultiset:
      case ROOT::kSTLunorderedmultiset:
         return new TGenSetProxy(*this);
      default:
         return new TGenCollectionProxy(*this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Proxy initializer

TGenCollectionProxy *TGenCollectionProxy::Initialize(Bool_t silent) const
{
   TGenCollectionProxy* p = const_cast<TGenCollectionProxy*>(this);
   if ( fValue.load() ) return p;
   return p->InitializeEx(silent);
}

////////////////////////////////////////////////////////////////////////////////
/// Check existence of function pointers

void TGenCollectionProxy::CheckFunctions() const
{
   if ( 0 == fSize.call ) {
      Fatal("TGenCollectionProxy","No 'size' function pointer for class %s present.",fName.c_str());
   }
   if ( 0 == fResize ) {
      Fatal("TGenCollectionProxy","No 'resize' function for class %s present.",fName.c_str());
   }
   if ( 0 == fNext.call  ) {
      Fatal("TGenCollectionProxy","No 'next' function for class %s present.",fName.c_str());
   }
   if ( 0 == fFirst.call ) {
      Fatal("TGenCollectionProxy","No 'begin' function for class %s present.",fName.c_str());
   }
   if ( 0 == fClear.call ) {
      Fatal("TGenCollectionProxy","No 'clear' function for class %s present.",fName.c_str());
   }
   if ( 0 == fConstruct ) {
      Fatal("TGenCollectionProxy","No 'block constructor' function for class %s present.",fName.c_str());
   }
   if ( 0 == fDestruct ) {
      Fatal("TGenCollectionProxy","No 'block destructor' function for class %s present.",fName.c_str());
   }
   if ( 0 == fFeed ) {
      Fatal("TGenCollectionProxy","No 'data feed' function for class %s present.",fName.c_str());
   }
   if ( 0 == fCollect ) {
      Fatal("TGenCollectionProxy","No 'data collect' function for class %s present.",fName.c_str());
   }
   if (0 == fCreateEnv.call ) {
      Fatal("TGenCollectionProxy","No 'environment creation' function for class %s present.",fName.c_str());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Utility routine to issue a Fatal error is the Value object is not valid

static TGenCollectionProxy::Value *R__CreateValue(const std::string &name, Bool_t silent)
{
   TGenCollectionProxy::Value *val = new TGenCollectionProxy::Value( name, silent );
   if ( !val->IsValid() ) {
      Fatal("TGenCollectionProxy","Could not find %s!",name.c_str());
   }
   return val;
}

////////////////////////////////////////////////////////////////////////////////
/// Proxy initializer

TGenCollectionProxy *TGenCollectionProxy::InitializeEx(Bool_t silent)
{
   R__LOCKGUARD(gInterpreterMutex);
   if (fValue.load()) return this;

   TClass *cl = fClass ? fClass.GetClass() : TClass::GetClass(fTypeinfo,kTRUE,silent);
   if ( cl ) {
      fEnv    = 0;
      fName   = cl->GetName();
      fPointers   = false;
      int nested = 0;
      std::vector<std::string> inside;
      int num = TClassEdit::GetSplit(cl->GetName(),inside,nested);
      if ( num > 1 ) {
         std::string nam;
         Value* newfValue = nullptr;
         if ( inside[0].find("stdext::hash_") != std::string::npos )
            inside[0].replace(3,10,"::");
         if ( inside[0].find("__gnu_cxx::hash_") != std::string::npos )
            inside[0].replace(0,16,"std::");
         fSTL_type = TClassEdit::STLKind(inside[0]);
         switch ( fSTL_type ) {
            case ROOT::kSTLmap:
            case ROOT::kSTLunorderedmap:
            case ROOT::kSTLmultimap:
            case ROOT::kSTLunorderedmultimap:
            case ROOT::kSTLset:
            case ROOT::kSTLunorderedset:
            case ROOT::kSTLmultiset:
            case ROOT::kSTLunorderedmultiset:
            case ROOT::kSTLbitset: // not really an associate container but it has no real iterator.
               fProperties |= kIsAssociative;
               if (num > 3 && !inside[3].empty()) {
                  if (! TClassEdit::IsDefAlloc(inside[3].c_str(),inside[0].c_str())) {
                     fProperties |= kCustomAlloc;
                  }
               }
               break;
         };

         int slong = sizeof(void*);
         switch ( fSTL_type ) {
            case ROOT::kSTLmap:
            case ROOT::kSTLunorderedmap:
            case ROOT::kSTLmultimap:
            case ROOT::kSTLunorderedmultimap:
               nam = "pair<"+inside[1]+","+inside[2];
               nam += (nam[nam.length()-1]=='>') ? " >" : ">";

               fVal   = R__CreateValue(inside[2], silent);
               fKey   = R__CreateValue(inside[1], silent);

               {
                  TInterpreter::SuspendAutoParsing autoParseRaii(gCling);
                  if (0==TClass::GetClass(nam.c_str(), true, false, fValOffset, fValDiff)) {
                     // We need to emulate the pair
                     TVirtualStreamerInfo::Factory()->GenerateInfoForPair(inside[1], inside[2], silent, fValOffset, fValDiff);
                  } else {
                     TClass *paircl = TClass::GetClass(nam.c_str());
                     if (paircl->GetClassSize() != fValDiff) {
                        if (paircl->GetState() >= TClass::kInterpreted)
                           Fatal("InitializeEx",
                                 "The %s for %s reports a class size that is inconsistent with the one registered "
                                 "through the CollectionProxy for %s:  %d vs %d\n",
                                 paircl->IsLoaded() ? "dictionary" : "interpreter information for", nam.c_str(),
                                 cl->GetName(), (int)paircl->GetClassSize(), (int)fValDiff);
                        else {
                           gROOT->GetListOfClasses()->Remove(paircl);
                           TClass *newpaircl = TClass::GetClass(nam.c_str(), true, false, fValOffset, fValDiff);
                           if (newpaircl == paircl || newpaircl->GetClassSize() != fValDiff)
                              Fatal("InitializeEx",
                                    "The TClass creation for %s did not get the right size: %d instead of%d\n",
                                    nam.c_str(), (int)paircl->GetClassSize(), (int)fValDiff);
                           paircl->ReplaceWith(newpaircl);
                           delete paircl;
                        }
                     }
                  }
               }
               newfValue = R__CreateValue(nam, silent);

               fPointers = (0 != (fKey->fCase&kIsPointer));
               if (fPointers || (0 != (fKey->fProperties&kNeedDelete))) {
                  fProperties |= kNeedDelete;
               }
               if ( 0 == fValDiff ) {
                  fValDiff = fKey->fSize + fVal->fSize;
                  fValDiff += (slong - fKey->fSize%slong)%slong;
                  fValDiff += (slong - fValDiff%slong)%slong;
               }
               if ( 0 == fValOffset ) {
                  fValOffset = fKey->fSize;
                  fValOffset += (slong - fKey->fSize%slong)%slong;
               }
               break;
            case ROOT::kSTLbitset:
               inside[1] = "bool";
               // Intentional fall through
            default:
               newfValue = R__CreateValue(inside[1], silent);

               fVal   = new Value(*newfValue);
               if ( 0 == fValDiff ) {
                  fValDiff = fVal->fSize;
                  fValDiff += (slong - fValDiff%slong)%slong;
               }
               if (num > 2 && !inside[2].empty()) {
                  if (! TClassEdit::IsDefAlloc(inside[2].c_str(),inside[0].c_str())) {
                     fProperties |= kCustomAlloc;
                  }
               }
               break;
         }

         fPointers = fPointers || (0 != (fVal->fCase&kIsPointer));
         if (fPointers || (0 != (fVal->fProperties&kNeedDelete))) {
            fProperties |= kNeedDelete;
         }
         fClass = cl;
         //fValue must be set last since we use it to indicate that we are initialized
         fValue = newfValue;
         return this;
      }
      Fatal("TGenCollectionProxy","Components of %s not analysed!",cl->GetName());
   }
   Fatal("TGenCollectionProxy","Collection class %s not found!",fTypeinfo.name());
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the TClass representing the container

TClass *TGenCollectionProxy::GetCollectionClass() const
{
   return fClass ? fClass : Initialize(kFALSE)->fClass;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the type of collection see TClassEdit::ESTLType

Int_t TGenCollectionProxy::GetCollectionType() const
{
   if (!fValue.load(std::memory_order_relaxed)) {
      Initialize(kFALSE);
   }
   return fSTL_type;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the offset between two consecutive value_types (memory layout).

ULong_t TGenCollectionProxy::GetIncrement() const {
   if (!fValue.load(std::memory_order_relaxed)) {
      Initialize(kFALSE);
   }
   return fValDiff;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the sizeof the collection object.

UInt_t TGenCollectionProxy::Sizeof() const
{
   return fClass->Size();
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the content is of type 'pointer to'

Bool_t TGenCollectionProxy::HasPointers() const
{
   // Initialize proxy in case it hasn't been initialized yet
   if( !fValue.load(std::memory_order_relaxed) )
      Initialize(kFALSE);

   // The content of a map and multimap is always a 'pair' and hence
   // fPointers means "Flag to indicate if containee has pointers (key or value)"
   // so we need to ignore its value for map and multimap;
   return fPointers && !(fSTL_type == ROOT::kSTLmap || fSTL_type == ROOT::kSTLmultimap ||
                         fSTL_type == ROOT::kSTLunorderedmap || fSTL_type == ROOT::kSTLunorderedmultimap);
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the TClass representing the content.

TClass *TGenCollectionProxy::GetValueClass() const
{
   auto value = fValue.load(std::memory_order_relaxed);
   if (!value) {
      Initialize(kFALSE);
      value = fValue.load(std::memory_order_relaxed);
   }
   return value ? (*value).fType.GetClass() : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// If the content is a simple numerical value, return its type (see TDataType)

EDataType TGenCollectionProxy::GetType() const
{
   auto value = fValue.load(std::memory_order_relaxed);
   if (!value) {
      Initialize(kFALSE);
      value = fValue.load(std::memory_order_relaxed);
   }
   return value ? (*value).fKind : kNoType_t;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the address of the value at index 'idx'

void* TGenCollectionProxy::At(UInt_t idx)
{
   if ( fEnv && fEnv->fObject ) {
      switch (fSTL_type) {
      case ROOT::kSTLvector:
         if ((*fValue).fKind == kBool_t) {
            auto vec = (std::vector<bool> *)(fEnv->fObject);
            fEnv->fLastValueVecBool = (*vec)[idx];
            fEnv->fIdx = idx;
            return &(fEnv->fLastValueVecBool);
         }
      // intentional fall through
      case ROOT::kROOTRVec:
         fEnv->fIdx = idx;
         switch( idx ) {
         case 0:
            return fEnv->fStart = fFirst.invoke(fEnv);
         default:
            if (! fEnv->fStart ) fEnv->fStart = fFirst.invoke(fEnv);
            return ((char*)fEnv->fStart) + fValDiff*idx;
         }
      case ROOT::kSTLbitset: {
         switch (idx) {
         case 0:
            fEnv->fStart = fFirst.invoke(fEnv);
            fEnv->fIdx = idx;
            break;
         default:
            fEnv->fIdx = idx - fEnv->fIdx;
            if (!fEnv->fStart) fEnv->fStart = fFirst.invoke(fEnv);
            fNext.invoke(fEnv);
            fEnv->fIdx = idx;
            break;
         }
         typedef ROOT::TCollectionProxyInfo::Environ <std::pair<size_t, Bool_t>> EnvType_t;
         EnvType_t *e = (EnvType_t *) fEnv;
         return &(e->fIterator.second);
      }
      case ROOT::kSTLset:
      case ROOT::kSTLunorderedset:
      case ROOT::kSTLmultiset:
      case ROOT::kSTLunorderedmultiset:
      case ROOT::kSTLmap:
      case ROOT::kSTLunorderedmap:
      case ROOT::kSTLmultimap:
      case ROOT::kSTLunorderedmultimap:
         if ( fEnv->fUseTemp ) {
            return (((char*)fEnv->fTemp)+idx*fValDiff);
         }
         // Intentional fall through.
      default:
         switch( idx ) {
         case 0:
            fEnv->fIdx = idx;
            return fEnv->fStart = fFirst.invoke(fEnv);
         default:  {
            fEnv->fIdx = idx - fEnv->fIdx;
            if (! fEnv->fStart ) fEnv->fStart = fFirst.invoke(fEnv);
            void* result = fNext.invoke(fEnv);
            fEnv->fIdx = idx;
            return result;
         }
         }
      }
   }
   Fatal("TGenCollectionProxy","At> Logic error - no proxy object set.");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Clear the emulated collection.

void TGenCollectionProxy::Clear(const char* opt)
{
   if ( fEnv && fEnv->fObject ) {
      if ( (fProperties & kNeedDelete) && opt && *opt=='f' ) {
         size_t i, n = *(size_t*)fSize.invoke(fEnv);
         if ( n > 0 ) {
            for (i=0; i<n; ++i)
               DeleteItem(true, TGenCollectionProxy::At(i));
         }
      }
      fClear.invoke(fEnv);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the current size of the container

UInt_t TGenCollectionProxy::Size() const
{
   if ( fEnv && fEnv->fObject ) {
      if (fEnv->fUseTemp) {
         return fEnv->fSize;
      } else {
         return *(size_t*)fSize.invoke(fEnv);
      }
   }
   Fatal("TGenCollectionProxy","Size> Logic error - no proxy object set.");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Resize the container

void TGenCollectionProxy::Resize(UInt_t n, Bool_t force)
{
   if ( fEnv && fEnv->fObject ) {
      if ( force && fPointers ) {
         size_t i, nold = *(size_t*)fSize.invoke(fEnv);
         if ( n != nold ) {
            for (i=n; i<nold; ++i)
               DeleteItem(true, *(void**)TGenCollectionProxy::At(i));
         }
      }
      MESSAGE(3, "Resize(n)" );
      fEnv->fSize = n;
      fResize(fEnv->fObject,fEnv->fSize);
      return;
   }
   Fatal("TGenCollectionProxy","Resize> Logic error - no proxy object set.");
}

////////////////////////////////////////////////////////////////////////////////
/// Allocate the needed space.
/// For associative collection, this returns a TStaging object that
/// need to be deleted manually __or__ returned by calling Commit(TStaging*)

void* TGenCollectionProxy::Allocate(UInt_t n, Bool_t /* forceDelete */ )
{
   if ( fEnv && fEnv->fObject ) {
      switch ( fSTL_type ) {
         case ROOT::kSTLset:
         case ROOT::kSTLunorderedset:
         case ROOT::kSTLmultiset:
         case ROOT::kSTLunorderedmultiset:
         case ROOT::kSTLmap:
         case ROOT::kSTLunorderedmap:
         case ROOT::kSTLmultimap:
         case ROOT::kSTLunorderedmultimap:{
            if ( (fProperties & kNeedDelete) )
               Clear("force");
            else
               fClear.invoke(fEnv);
            // Commit no longer use the environment and thus no longer decrease
            // the count.  Consequently we no longer should increase it here.
            // ++fEnv->fRefCount;
            fEnv->fSize  = n;

            TStaging *s;
            if (fStaged.empty()) {
               s = new TStaging(n,fValDiff);
            } else {
               s = fStaged.back();
               fStaged.pop_back();
               s->Resize(n);
            }
            fConstruct(s->GetContent(),s->GetSize());

            s->SetTarget(fEnv->fObject);

            fEnv->fTemp = s->GetContent();
            fEnv->fUseTemp = kTRUE;
            fEnv->fStart = fEnv->fTemp;

            return s;
         }
         case ROOT::kSTLvector:
         case ROOT::kSTLlist:
         case ROOT::kSTLforwardlist:
         case ROOT::kSTLdeque:
         case ROOT::kROOTRVec:
            if( (fProperties & kNeedDelete) ) {
               Clear("force");
            }
            fEnv->fSize = n;
            fResize(fEnv->fObject,n);
            return fEnv->fObject;

        case ROOT::kSTLbitset: {
            TStaging *s;
            if (fStaged.empty()) {
               s = new TStaging(n,fValDiff);
            } else {
               s = fStaged.back();
               fStaged.pop_back();
               s->Resize(n);
            }
            s->SetTarget(fEnv->fObject);

            fEnv->fTemp = s->GetContent();
            fEnv->fUseTemp = kTRUE;
            fEnv->fStart = fEnv->fTemp;

            return s;
        }
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Insert data into the container where data is a C-style array of the actual type contained in the collection
/// of the given size.   For associative container (map, etc.), the data type is the pair<key,value>.

void TGenCollectionProxy::Insert(const void *data, void *container, size_t size)
{
   fFeed((void*)data,container,size);
}

////////////////////////////////////////////////////////////////////////////////
/// Commit the change.

void TGenCollectionProxy::Commit(void* from)
{
   if (fProperties & kIsAssociative) {
//      case ROOT::kSTLmap:
//      case ROOT::kSTLmultimap:
//      case ROOT::kSTLset:
//      case ROOT::kSTLmultiset:
      if ( from ) {
         TStaging *s = (TStaging*) from;
         if ( s->GetTarget() ) {
            fFeed(s->GetContent(),s->GetTarget(),s->GetSize());
         }
         fDestruct(s->GetContent(),s->GetSize());
         s->SetTarget(0);
         fStaged.push_back(s);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add an object.

void TGenCollectionProxy::PushProxy(void *objstart)
{
   if ( !fValue.load(std::memory_order_relaxed) ) Initialize(kFALSE);
   if ( !fProxyList.empty() ) {
      EnvironBase_t* back = fProxyList.back();
      if ( back->fObject == objstart ) {
         ++back->fRefCount;
         fProxyList.push_back(back);
         fEnv = back;
         return;
      }
   }
   EnvironBase_t* e    = 0;
   if ( fProxyKept.empty() ) {
      e = (EnvironBase_t*)fCreateEnv.invoke();
      e->fTemp  = 0;
      e->fUseTemp = kFALSE;
   }
   else {
      e = fProxyKept.back();
      fProxyKept.pop_back();
   }
   e->fSize     = 0;
   e->fRefCount = 1;
   e->fObject   = objstart;
   e->fStart    = 0;
   e->fIdx      = 0;
   // ::memset(e->buff,0,sizeof(e->buff));
   fProxyList.push_back(e);
   fEnv = e;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove the last object.

void TGenCollectionProxy::PopProxy()
{
   if ( !fProxyList.empty() ) {
      EnvironBase_t* e = fProxyList.back();
      if ( --e->fRefCount <= 0 ) {
         fProxyKept.push_back(e);
         e->fUseTemp = kFALSE;
      }
      fProxyList.pop_back();
   }
   fEnv = fProxyList.empty() ? 0 : fProxyList.back();
}

////////////////////////////////////////////////////////////////////////////////
/// Call to delete/destruct individual item.

void TGenCollectionProxy::DeleteItem(Bool_t force, void* ptr) const
{
   if ( force && ptr ) {
      switch (fSTL_type) {
         case ROOT::kSTLmap:
         case ROOT::kSTLunorderedmap:
         case ROOT::kSTLmultimap:
         case ROOT::kSTLunorderedmultimap:{
            if ( fKey->fCase&kIsPointer ) {
               if (fKey->fProperties&kNeedDelete) {
                  TVirtualCollectionProxy *proxy = fKey->fType->GetCollectionProxy();
                  TPushPop helper(proxy,*(void**)ptr);
                  proxy->Clear("force");
               }
               fKey->DeleteItem(*(void**)ptr);
             } else {
               if (fKey->fProperties&kNeedDelete) {
                  TVirtualCollectionProxy *proxy = fKey->fType->GetCollectionProxy();
                  TPushPop helper(proxy,ptr);
                  proxy->Clear("force");
               }
            }
            char *addr = ((char*)ptr)+fValOffset;
            if ( fVal->fCase&kIsPointer ) {
               if ( fVal->fProperties&kNeedDelete) {
                  TVirtualCollectionProxy *proxy = fVal->fType->GetCollectionProxy();
                  TPushPop helper(proxy,*(void**)addr);
                  proxy->Clear("force");
               }
               fVal->DeleteItem(*(void**)addr);
           } else {
               if ( fVal->fProperties&kNeedDelete) {
                  TVirtualCollectionProxy *proxy = fVal->fType->GetCollectionProxy();
                  TPushPop helper(proxy,addr);
                  proxy->Clear("force");
               }
            }
            break;
         }
         default: {
            if ( fVal->fCase&kIsPointer ) {
               if (fVal->fProperties&kNeedDelete) {
                  TVirtualCollectionProxy *proxy = fVal->fType->GetCollectionProxy();
                  TPushPop helper(proxy,*(void**)ptr);
                  proxy->Clear("force");
               }
               fVal->DeleteItem(*(void**)ptr);
            } else {
               if (fVal->fProperties&kNeedDelete) {
                  TVirtualCollectionProxy *proxy = fVal->fType->GetCollectionProxy();
                  TPushPop helper(proxy,ptr);
                  proxy->Clear("force");
               }
            }
            break;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

void TGenCollectionProxy::ReadBuffer(TBuffer & /* b */, void * /* obj */, const TClass * /* onfileClass */)
{
   MayNotUse("TGenCollectionProxy::ReadBuffer(TBuffer &, void *, const TClass *)");
}

////////////////////////////////////////////////////////////////////////////////

void TGenCollectionProxy::ReadBuffer(TBuffer & /* b */, void * /* obj */)
{
   MayNotUse("TGenCollectionProxy::ReadBuffer(TBuffer &, void *)");
}

////////////////////////////////////////////////////////////////////////////////
/// Streamer Function.

void TGenCollectionProxy::Streamer(TBuffer &buff)
{
   if ( fEnv ) {
      GetCollectionClass()->Streamer( fEnv->fObject, buff );
      return;
   }
   Fatal("TGenCollectionProxy","Streamer> Logic error - no proxy object set.");
}

////////////////////////////////////////////////////////////////////////////////
/// Streamer I/O overload

void TGenCollectionProxy::Streamer(TBuffer &buff, void *objp, int /* siz */ )
{
   TPushPop env(this, objp);
   Streamer(buff);
}

////////////////////////////////////////////////////////////////////////////////
/// TClassStreamer IO overload

void TGenCollectionProxy::operator()(TBuffer &b, void *objp)
{
   Streamer(b, objp, 0);
}


struct TGenCollectionProxy__SlowIterator {
   TVirtualCollectionProxy *fProxy;
   UInt_t fIndex;
   TGenCollectionProxy__SlowIterator(TVirtualCollectionProxy *proxy) : fProxy(proxy), fIndex(0) {}
};

////////////////////////////////////////////////////////////////////////////////

void TGenCollectionProxy__SlowCreateIterators(void * /* collection */, void **begin_arena, void **end_arena, TVirtualCollectionProxy *proxy)
{
   new (*begin_arena) TGenCollectionProxy__SlowIterator(proxy);
   *(UInt_t*)*end_arena = proxy->Size();
}

////////////////////////////////////////////////////////////////////////////////

void *TGenCollectionProxy__SlowNext(void *iter, const void *end)
{
   TGenCollectionProxy__SlowIterator *iterator = (TGenCollectionProxy__SlowIterator*)iter;
   if (iterator->fIndex != *(UInt_t*)end) {
      void *result = iterator->fProxy->At(iterator->fIndex);
      ++(iterator->fIndex);
      return result;
   } else {
      return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////

void * TGenCollectionProxy__SlowCopyIterator(void *dest, const void *source)
{
   *(TGenCollectionProxy__SlowIterator*)dest = *(TGenCollectionProxy__SlowIterator*)source;
   return dest;
}

////////////////////////////////////////////////////////////////////////////////
/// Nothing to do

void TGenCollectionProxy__SlowDeleteSingleIterators(void *)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Nothing to do

void TGenCollectionProxy__SlowDeleteTwoIterators(void *, void *)
{
}


////////////////////////////////////////////////////////////////////////////////
/// We can safely assume that the std::vector layout does not really depend on
/// the content!

void TGenCollectionProxy__VectorCreateIterators(void *obj, void **begin_arena, void **end_arena, TVirtualCollectionProxy*)
{
   std::vector<char> *vec = (std::vector<char>*)obj;
   if (vec->empty()) {
      *begin_arena = 0;
      *end_arena = 0;
      return;
   }
   *begin_arena = &(*vec->begin());
#ifdef R__VISUAL_CPLUSPLUS
   *end_arena = &(*(vec->end()-1)) + 1; // On windows we can not dererence the end iterator at all.
#else
   // coverity[past_the_end] Safe on other platforms
   *end_arena = &(*vec->end());
#endif

}

////////////////////////////////////////////////////////////////////////////////
/// Should not be used.

void *TGenCollectionProxy__VectorNext(void *, const void *)
{
   R__ASSERT(0);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////

void *TGenCollectionProxy__VectorCopyIterator(void *dest, const void *source)
{
   *(void**)dest = *(void**)source;
   return dest;
}

////////////////////////////////////////////////////////////////////////////////
/// Nothing to do

void TGenCollectionProxy__VectorDeleteSingleIterators(void *)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Nothing to do

void TGenCollectionProxy__VectorDeleteTwoIterators(void *, void *)
{
}



////////////////////////////////////////////////////////////////////////////////

void TGenCollectionProxy__StagingCreateIterators(void *obj, void **begin_arena, void **end_arena, TVirtualCollectionProxy *)
{
   TGenCollectionProxy::TStaging * s = (TGenCollectionProxy::TStaging*)obj;
   *begin_arena = s->GetContent();
   *end_arena = s->GetEnd();
}

////////////////////////////////////////////////////////////////////////////////
/// Should not be used.

void *TGenCollectionProxy__StagingNext(void *, const void *)
{
   R__ASSERT(0);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////

void *TGenCollectionProxy__StagingCopyIterator(void *dest, const void *source)
{
   *(void**)dest = *(void**)source;
   return dest;
}

////////////////////////////////////////////////////////////////////////////////
/// Nothing to do

void TGenCollectionProxy__StagingDeleteSingleIterators(void *)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Nothing to do

void TGenCollectionProxy__StagingDeleteTwoIterators(void *, void *)
{
}


////////////////////////////////////////////////////////////////////////////////
/// See typedef void (*CreateIterators_t)(void *collection, void *&begin_arena, void *&end_arena);
/// begin_arena and end_arena should contain the location of memory arena  of size fgIteratorSize.
/// If the collection iterator are of that size or less, the iterators will be constructed in place in those location (new with placement)
/// Otherwise the iterators will be allocated via a regular new and their address returned by modifying the value of begin_arena and end_arena.

TVirtualCollectionProxy::CreateIterators_t TGenCollectionProxy::GetFunctionCreateIterators(Bool_t read)
{
   if (read) {
      if ( !fValue.load(std::memory_order_relaxed) ) InitializeEx(kFALSE);
      if ( (fProperties & kIsAssociative) && read)
         return TGenCollectionProxy__StagingCreateIterators;
   }

   if ( fFunctionCreateIterators ) return fFunctionCreateIterators;

   if ( !fValue.load(std::memory_order_relaxed) ) InitializeEx(kFALSE);

//   fprintf(stderr,"GetFunctinCreateIterator for %s will give: ",fClass.GetClassName());
//   if (fSTL_type==ROOT::kSTLvector || (fProperties & kIsEmulated))
//      fprintf(stderr,"vector/emulated iterator\n");
//   else if ( (fProperties & kIsAssociative) && read)
//      fprintf(stderr,"an associative read iterator\n");
//   else
//      fprintf(stderr,"a generic iterator\n");

   // TODO could we do better than SlowCreateIterators for RVec?
   if (fSTL_type==ROOT::kSTLvector || (fProperties & kIsEmulated))
      return fFunctionCreateIterators = TGenCollectionProxy__VectorCreateIterators;
   else if ( (fProperties & kIsAssociative) && read)
      return TGenCollectionProxy__StagingCreateIterators;
   else
      return fFunctionCreateIterators = TGenCollectionProxy__SlowCreateIterators;
}

////////////////////////////////////////////////////////////////////////////////
/// See typedef void (*CopyIterator_t)(void *&dest, const void *source);
/// Copy the iterator source, into dest.   dest should contain should contain the location of memory arena  of size fgIteratorSize.
/// If the collection iterator are of that size or less, the iterator will be constructed in place in this location (new with placement)
/// Otherwise the iterator will be allocated via a regular new and its address returned by modifying the value of dest.

TVirtualCollectionProxy::CopyIterator_t TGenCollectionProxy::GetFunctionCopyIterator(Bool_t read)
{
   if (read) {
      if ( !fValue.load(std::memory_order_relaxed) ) InitializeEx(kFALSE);
      if ( (fProperties & kIsAssociative) && read)
         return TGenCollectionProxy__StagingCopyIterator;
   }

   if ( fFunctionCopyIterator ) return fFunctionCopyIterator;

   if ( !fValue.load(std::memory_order_relaxed) ) InitializeEx(kFALSE);

   // TODO can we do better than the default for RVec?
   if (fSTL_type==ROOT::kSTLvector || (fProperties & kIsEmulated))
      return fFunctionCopyIterator = TGenCollectionProxy__VectorCopyIterator;
   else if ( (fProperties & kIsAssociative) && read)
      return TGenCollectionProxy__StagingCopyIterator;
   else
      return fFunctionCopyIterator = TGenCollectionProxy__SlowCopyIterator;
}

////////////////////////////////////////////////////////////////////////////////
/// See typedef void* (*Next_t)(void *iter, void *end);
/// iter and end should be pointer to respectively an iterator to be incremented and the result of colleciton.end()
/// 'Next' will increment the iterator 'iter' and return 0 if the iterator reached the end.
/// If the end is not reached, 'Next' will return the address of the content unless the collection contains pointers in
/// which case 'Next' will return the value of the pointer.

TVirtualCollectionProxy::Next_t TGenCollectionProxy::GetFunctionNext(Bool_t read)
{
   if (read) {
      if ( !fValue.load(std::memory_order_relaxed) ) InitializeEx(kFALSE);
      if ( (fProperties & kIsAssociative) && read)
         return TGenCollectionProxy__StagingNext;
   }

   if ( fFunctionNextIterator ) return fFunctionNextIterator;

   if ( !fValue.load(std::memory_order_relaxed) ) InitializeEx(kFALSE);

   // TODO can we do better than the default for RVec?
   if (fSTL_type==ROOT::kSTLvector || (fProperties & kIsEmulated))
      return fFunctionNextIterator = TGenCollectionProxy__VectorNext;
   else if ( (fProperties & kIsAssociative) && read)
      return TGenCollectionProxy__StagingNext;
   else
      return fFunctionNextIterator = TGenCollectionProxy__SlowNext;
}

////////////////////////////////////////////////////////////////////////////////
/// See typedef void (*DeleteIterator_t)(void *iter);
/// If the sizeof iterator is greater than fgIteratorArenaSize, call delete on the addresses,
/// Otherwise just call the iterator's destructor.

TVirtualCollectionProxy::DeleteIterator_t TGenCollectionProxy::GetFunctionDeleteIterator(Bool_t read)
{
   if (read) {
      if ( !fValue.load(std::memory_order_relaxed) ) InitializeEx(kFALSE);
      if ( (fProperties & kIsAssociative) && read)
         return TGenCollectionProxy__StagingDeleteSingleIterators;
   }

   if ( fFunctionDeleteIterator ) return fFunctionDeleteIterator;

   if ( !fValue.load(std::memory_order_relaxed) ) InitializeEx(kFALSE);

   // TODO can we do better than the default for RVec?
   if (fSTL_type==ROOT::kSTLvector || (fProperties & kIsEmulated))
      return fFunctionDeleteIterator = TGenCollectionProxy__VectorDeleteSingleIterators;
   else if ( (fProperties & kIsAssociative) && read)
      return TGenCollectionProxy__StagingDeleteSingleIterators;
   else
      return fFunctionDeleteIterator = TGenCollectionProxy__SlowDeleteSingleIterators;
}

////////////////////////////////////////////////////////////////////////////////
/// See typedef void (*DeleteTwoIterators_t)(void *begin, void *end);
/// If the sizeof iterator is greater than fgIteratorArenaSize, call delete on the addresses,
/// Otherwise just call the iterator's destructor.

TVirtualCollectionProxy::DeleteTwoIterators_t TGenCollectionProxy::GetFunctionDeleteTwoIterators(Bool_t read)
{
   if (read) {
      if ( !fValue.load(std::memory_order_relaxed) ) InitializeEx(kFALSE);
      if ( (fProperties & kIsAssociative) && read)
         return TGenCollectionProxy__StagingDeleteTwoIterators;
   }

   if ( fFunctionDeleteTwoIterators ) return fFunctionDeleteTwoIterators;

   if ( !fValue.load(std::memory_order_relaxed) ) InitializeEx(kFALSE);

   // TODO could RVec use something faster than SlowCopyIterator?
   if (fSTL_type==ROOT::kSTLvector || (fProperties & kIsEmulated))
      return fFunctionDeleteTwoIterators = TGenCollectionProxy__VectorDeleteTwoIterators;
   else if ( (fProperties & kIsAssociative) && read)
      return TGenCollectionProxy__StagingDeleteTwoIterators;
   else
      return fFunctionDeleteTwoIterators = TGenCollectionProxy__SlowDeleteTwoIterators;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the set of action necessary to stream in this collection member-wise coming from
/// the old value class layout refered to by 'version'.

TStreamerInfoActions::TActionSequence *TGenCollectionProxy::GetConversionReadMemberWiseActions(TClass *oldClass, Int_t version)
{
   if (oldClass == 0) {
      return 0;
   }
   TObjArray* arr = 0;
   TStreamerInfoActions::TActionSequence *result = 0;
   if (fConversionReadMemberWise) {
      std::map<std::string, TObjArray*>::iterator it;

      it = fConversionReadMemberWise->find( oldClass->GetName() );

      if( it != fConversionReadMemberWise->end() ) {
         arr = it->second;
      }

      if (arr) {
         result = (TStreamerInfoActions::TActionSequence *)arr->At(version);
         if (result) {
            return result;
         }
      }
   }

   // Need to create it.
   TClass *valueClass = GetValueClass();
   if (valueClass == 0) {
      return 0;
   }
   TVirtualStreamerInfo *info = valueClass->GetConversionStreamerInfo(oldClass,version);
   if (info == 0) {
      return 0;
   }
   result = TStreamerInfoActions::TActionSequence::CreateReadMemberWiseActions(info,*this);

   if (!arr) {
      arr = new TObjArray(version+10, -1);
      if (!fConversionReadMemberWise) {
         fConversionReadMemberWise = new std::map<std::string, TObjArray*>();
      }
      (*fConversionReadMemberWise)[oldClass->GetName()] = arr;
   }
   arr->AddAtAndExpand( result, version );

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the set of action necessary to stream in this collection member-wise coming from
/// the old value class layout refered to by 'version'.

TStreamerInfoActions::TActionSequence *TGenCollectionProxy::GetReadMemberWiseActions(Int_t version)
{
   TStreamerInfoActions::TActionSequence *result = 0;
   if (version < (fReadMemberWise->GetSize()-1)) { // -1 because the 'index' starts at -1
      result = (TStreamerInfoActions::TActionSequence *)fReadMemberWise->At(version);
   }
   if (result == 0) {
      // Need to create it.
      TClass *valueClass = GetValueClass();
      TVirtualStreamerInfo *info = 0;
      if (valueClass) {
         info = valueClass->GetStreamerInfo(version);
      }
      result = TStreamerInfoActions::TActionSequence::CreateReadMemberWiseActions(info,*this);
      fReadMemberWise->AddAtAndExpand(result,version);
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the set of action necessary to stream out this collection member-wise.

TStreamerInfoActions::TActionSequence *TGenCollectionProxy::GetWriteMemberWiseActions()
{
  TStreamerInfoActions::TActionSequence *result = fWriteMemberWise;
  if (result == 0) {
     // Need to create it.
     TClass *valueClass = GetValueClass();
     TVirtualStreamerInfo *info = 0;
     if (valueClass) {
        info = valueClass->GetStreamerInfo();
     }
     result = TStreamerInfoActions::TActionSequence::CreateWriteMemberWiseActions(info,*this);
     fWriteMemberWise=result;
  }
  return result;
}
