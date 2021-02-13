// @(#)root/cont:$Id$
// Author: Philippe Canal 20/08/2003

/*************************************************************************
 * Copyright (C) 1995-2003, Rene Brun, Fons Rademakers and al.           *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TVirtualCollectionProxy
#define ROOT_TVirtualCollectionProxy

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TVirtualCollectionProxy                                              //
//                                                                      //
// Virtual interface of a proxy object for a collection class           //
// In particular this is used to implement splitting, emulation,        //
// and TTreeFormula access to STL containers.                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TClassRef.h"
#include "TDataType.h"

// Macro indicating the version of the Collection Proxy interface followed
// by this ROOT build (See also Reflex/Builder/CollectionProxy.h).

#define ROOT_COLLECTIONPROXY_VERSION 3

class TClass;
namespace TStreamerInfoActions {
   class TActionSequence;
}

class TVirtualCollectionProxy {
private:
   TVirtualCollectionProxy(const TVirtualCollectionProxy&) = delete;
   TVirtualCollectionProxy& operator=(const TVirtualCollectionProxy&) = delete;

protected:
   TClassRef fClass;
   UInt_t    fProperties;
   friend class TClass;

public:
   enum EProperty {
      // No longer used
      // kIsInitialized = BIT(1),
      kIsAssociative = BIT(2),
      kIsEmulated    = BIT(3),
      kNeedDelete    = BIT(4),  // Flag to indicate that this collection that contains directly or indirectly (only via other collection) some pointers that will need explicit deletions.
      kCustomAlloc   = BIT(5)   // The collection has a custom allocator.
   };

   class TPushPop {
      // Helper class that insures that push and pop are done when entering
      // and leaving a C++ context (even in the presence of exceptions)
   public:
      TVirtualCollectionProxy *fProxy;
      inline TPushPop(TVirtualCollectionProxy *proxy,
         void *objectstart) : fProxy(proxy) { fProxy->PushProxy(objectstart); }
      inline ~TPushPop() { fProxy->PopProxy(); }
   private:
      TPushPop(const TPushPop&) = delete;
      TPushPop& operator=(const TPushPop&) = delete;
   };

   TVirtualCollectionProxy() : fClass(), fProperties(0) {};
   TVirtualCollectionProxy(TClass *cl) : fClass(cl), fProperties(0) {};

   virtual TVirtualCollectionProxy* Generate() const = 0; // Returns an object of the actual CollectionProxy class
   virtual ~TVirtualCollectionProxy() {};

   virtual TClass   *GetCollectionClass() const { return fClass; }
   // Return a pointer to the TClass representing the container

   virtual Int_t     GetCollectionType() const = 0;
   // Return the type of collection see TClassEdit::ESTLType

   virtual ULong_t   GetIncrement() const = 0;
   // Return the offset between two consecutive value_types (memory layout).

   virtual Int_t     GetProperties() const { return fProperties; }
   // Return miscallenous properties of the proxy see TVirtualCollectionProxy::EProperty

   virtual void     *New() const {
      // Return a new container object
      return fClass.GetClass()==0 ? 0 : fClass->New();
   }
   virtual void     *New(void *arena) const {
      // Execute the container constructor
      return fClass.GetClass()==0 ? 0 : fClass->New(arena);
   }
   virtual TClass::ObjectPtr NewObject() const {
      // Return a new container object
      return fClass.GetClass()==0 ? TClass::ObjectPtr{} : fClass->NewObject();
   }
   virtual TClass::ObjectPtr NewObject(void *arena) const {
      // Execute the container constructor
      return fClass.GetClass()==0 ? TClass::ObjectPtr{} : fClass->NewObject(arena);
   }

   virtual void     *NewArray(Int_t nElements) const {
      // Return a new container object
      return fClass.GetClass()==0 ? 0 : fClass->NewArray(nElements);
   }
   virtual void     *NewArray(Int_t nElements, void *arena) const {
      // Execute the container constructor
      return fClass.GetClass()==0 ? 0 : fClass->NewArray(nElements, arena);
   }
   virtual TClass::ObjectPtr NewObjectArray(Int_t nElements) const {
      // Return a new container object
      return fClass.GetClass()==0 ? TClass::ObjectPtr{} : fClass->NewObjectArray(nElements);
   }
   virtual TClass::ObjectPtr NewObjectArray(Int_t nElements, void *arena) const {
      // Execute the container constructor
      return fClass.GetClass()==0 ? TClass::ObjectPtr{} : fClass->NewObjectArray(nElements, arena);
   }

   virtual void      Destructor(void *p, Bool_t dtorOnly = kFALSE) const {
      // Execute the container destructor
      TClass* cl = fClass.GetClass();
      if (cl) cl->Destructor(p, dtorOnly);
   }

   virtual void      DeleteArray(void *p, Bool_t dtorOnly = kFALSE) const {
      // Execute the container array destructor
      TClass* cl = fClass.GetClass();
      if (cl) cl->DeleteArray(p, dtorOnly);
   }

   virtual UInt_t    Sizeof() const = 0;
   // Return the sizeof the collection object.

   virtual void      PushProxy(void *objectstart) = 0;
   // Set the address of the container being proxied and keep track of the previous one.

   virtual void      PopProxy() = 0;
   // Reset the address of the container being proxied to the previous container

   virtual Bool_t    HasPointers() const = 0;
   // Return true if the content is of type 'pointer to'

   virtual TClass   *GetValueClass() const = 0;
   // Return a pointer to the TClass representing the content.

   virtual EDataType GetType() const = 0;
   // If the content is a simple numerical value, return its type (see TDataType)

   virtual void     *At(UInt_t idx) = 0;
   // Return the address of the value at index 'idx'

   virtual void      Clear(const char *opt = "") = 0;
   // Clear the container

   virtual UInt_t    Size() const = 0;
   // Return the current size of the container

   virtual void*     Allocate(UInt_t n, Bool_t forceDelete) = 0;

   virtual void      Commit(void*) = 0;

   virtual void      Insert(const void *data, void *container, size_t size) = 0;
   // Insert data into the container where data is a C-style array of the actual type contained in the collection
   // of the given size.   For associative container (map, etc.), the data type is the pair<key,value>.

           char     *operator[](UInt_t idx) const { return (char*)(const_cast<TVirtualCollectionProxy*>(this))->At(idx); }

   // MemberWise actions
   virtual TStreamerInfoActions::TActionSequence *GetConversionReadMemberWiseActions(TClass *oldClass, Int_t version) = 0;
   virtual TStreamerInfoActions::TActionSequence *GetReadMemberWiseActions(Int_t version) = 0;
   virtual TStreamerInfoActions::TActionSequence *GetWriteMemberWiseActions() = 0;

   // Set of functions to iterate easily throught the collection
   static const Int_t fgIteratorArenaSize = 16; // greater than sizeof(void*) + sizeof(UInt_t)

   typedef void (*CreateIterators_t)(void *collection, void **begin_arena, void **end_arena, TVirtualCollectionProxy *proxy);
   virtual CreateIterators_t GetFunctionCreateIterators(Bool_t read = kTRUE) = 0;
   // begin_arena and end_arena should contain the location of a memory arena of size fgIteratorSize.
   // If the collection iterator are of that size or less, the iterators will be constructed in place in those location (new with placement)
   // Otherwise the iterators will be allocated via a regular new and their address returned by modifying the value of begin_arena and end_arena.

   typedef void* (*CopyIterator_t)(void *dest, const void *source);
   virtual CopyIterator_t GetFunctionCopyIterator(Bool_t read = kTRUE) = 0;
   // Copy the iterator source, into dest.   dest should contain the location of a memory arena of size fgIteratorSize.
   // If the collection iterator is of that size or less, the iterator will be constructed in place in this location (new with placement)
   // Otherwise the iterator will be allocated via a regular new.
   // The actual address of the iterator is returned in both case.

   typedef void* (*Next_t)(void *iter, const void *end);
   virtual Next_t GetFunctionNext(Bool_t read = kTRUE) = 0;
   // iter and end should be pointers to respectively an iterator to be incremented and the result of collection.end()
   // If the iterator has not reached the end of the collection, 'Next' increment the iterator 'iter' and return 0 if
   // the iterator reached the end.
   // If the end was not reached, 'Next' returns the address of the content pointed to by the iterator before the
   // incrementation ; if the collection contains pointers, 'Next' will return the value of the pointer.

   typedef void (*DeleteIterator_t)(void *iter);
   typedef void (*DeleteTwoIterators_t)(void *begin, void *end);

   virtual DeleteIterator_t GetFunctionDeleteIterator(Bool_t read = kTRUE) = 0;
   virtual DeleteTwoIterators_t GetFunctionDeleteTwoIterators(Bool_t read = kTRUE) = 0;
   // If the size of the iterator is greater than fgIteratorArenaSize, call delete on the addresses,
   // Otherwise just call the iterator's destructor.

};

#endif
