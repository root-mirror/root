// @(#)root/meta:$Id$
// Author: Markus Frank 20/05/2005

/*************************************************************************
* Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
* All rights reserved.                                                  *
*                                                                       *
* For the licensing terms see $ROOTSYS/LICENSE.                         *
* For the list of contributors see $ROOTSYS/README/CREDITS.             *
*************************************************************************/

#include "TClass.h"
#include "TError.h"
#include "TInterpreter.h"
#include "TIsAProxy.h"

#include <map>
#include <type_traits>

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TClass                                                               //
//                                                                      //
// TIsAProxy implementation class.                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

namespace {
   struct DynamicType {
      // Helper class to enable typeid on any address
      // Used in code similar to:
      //    typeid( * (DynamicType*) void_ptr );
      virtual ~DynamicType() {}
   };

   typedef std::map<const void*, TClass*> ClassMap_t; // Internal type map
   inline ClassMap_t *GetMap(const void* p)
   {
      return (ClassMap_t*)p;
   }

   inline ClassMap_t::value_type* ToPair(void*p) 
   {
     return (ClassMap_t::value_type*)p;
   }
}

//______________________________________________________________________________
TIsAProxy::TIsAProxy(const std::type_info& typ, void* ctxt)
   : fType(&typ), fClass(nullptr), fLast(nullptr), fContext(ctxt),
     fSubTypesReaders(0), fSubTypesWriteLockTaken(false),
     fVirtual(false), fInit(false)
{
   // Standard initializing constructor

   ::new(fSubTypes) ClassMap_t();
   static_assert(sizeof(ClassMap_t)<=sizeof(fSubTypes), "ClassMap size is to large for array");
}

//______________________________________________________________________________
TIsAProxy::~TIsAProxy()
{
   // Standard destructor

   ClassMap_t* m = GetMap(fSubTypes);
   m->clear();
   m->~ClassMap_t();
}

//______________________________________________________________________________
void TIsAProxy::SetClass(TClass *cl)
{
   // Set class pointer
   //   This method is not thread safe
   GetMap(fSubTypes)->clear();
   fClass = cl;
   fLast = nullptr;
}

//______________________________________________________________________________
TClass* TIsAProxy::operator()(const void *obj)
{
   // IsA callback

   if ( !fInit )  {
      if ( !fClass && fType ) {
         auto cls = TClass::GetClass(*fType);
         TClass* expected = nullptr;
         fClass.compare_exchange_strong(expected,cls);
      }
      if ( !fClass) return nullptr;
      fVirtual = (*fClass).ClassProperty() & kClassHasVirtual;
      fInit = kTRUE;
   }
   if ( !obj || !fVirtual )  {
      return fClass.load();
   }
   // Avoid the case that the first word is a virtual_base_offset_table instead of
   // a virtual_function_table
   Long_t offset = **(Long_t**)obj;
   if ( offset == 0 ) return fClass.load();
   
   DynamicType* ptr = (DynamicType*)obj;
   const std::type_info* typ = &typeid(*ptr);
   
   if ( typ == fType )  {
     return fClass.load();
   }
   auto last = ToPair(fLast.load());
   if ( last && typ == last->first )  {
     return last->second;
   }
   // Check if type is already in sub-class cache
   if ( nullptr == (last = ToPair(FindSubType(typ)) ) )  {
     
     // Last resort: lookup root class
     auto cls = TClass::GetClass(*typ);
     last = ToPair(CacheSubType(typ,cls));
   }
   fLast.store(last);

   return last == nullptr? nullptr: last->second;
}
//______________________________________________________________________________
inline void* TIsAProxy::FindSubType(const type_info* type) const
{
   bool needToWait = true;
   do {
     ++fSubTypesReaders;

     //See if there is a writer, if there is we need to release
     // our reader count so that the writer can proceed
     if(fSubTypesWriteLockTaken) {
       --fSubTypesReaders;
       while(fSubTypesWriteLockTaken) {}
     } else {
       needToWait = false;
     }
   } while(needToWait);

   void* returnValue =nullptr;
   auto const map = GetMap(fSubTypes);

   auto found = map->find(type);
   if(found != map->end()) {
      returnValue = &(*found);
   }
   --fSubTypesReaders;
   return returnValue;
}

//______________________________________________________________________________
void* TIsAProxy::CacheSubType(const type_info* type, TClass* cls)
{
   //See if another thread has the write lock, wait if it does
   bool expected = false;
   while(not fSubTypesWriteLockTaken.compare_exchange_strong(expected,true) ) {expected = false;};

   //See if there are any readers
   while(fSubTypesReaders > 0);

   auto map = GetMap(fSubTypes);
   auto ret = map->emplace(type,cls);

   fSubTypesWriteLockTaken = false;
   return &(*(ret.first));
}
