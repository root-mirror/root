// @(#)root/io:$Id$
// Author: Markus Frank  28/10/04

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifndef ROOT_TEmulatedCollectionProxy
#define ROOT_TEmulatedCollectionProxy

#include "TGenCollectionProxy.h"

#include <vector>

class TEmulatedCollectionProxy : public TGenCollectionProxy  {

   // Friend declaration
   friend class TCollectionProxy;

public:
   // Container type definition
   typedef std::vector<char>  Cont_t;
   // Pointer to container type
   typedef Cont_t            *PCont_t;
protected:

   // Some hack to avoid const-ness
   virtual TGenCollectionProxy* InitializeEx(Bool_t silent);

   // Object input streamer
   void ReadItems(int nElements, TBuffer &b);

   // Object output streamer
   void WriteItems(int nElements, TBuffer &b);

   // Shrink the container
   void Shrink(UInt_t nCurr, UInt_t left, Bool_t force);

   // Expand the container
   void Expand(UInt_t nCurr, UInt_t left);

private:
   TEmulatedCollectionProxy &operator=(const TEmulatedCollectionProxy &); // Not implemented.

public:
   // Virtual copy constructor
   virtual TVirtualCollectionProxy* Generate() const;

   // Copy constructor
   TEmulatedCollectionProxy(const TEmulatedCollectionProxy& copy);

   // Initializing constructor
   TEmulatedCollectionProxy(const char* cl_name, Bool_t silent);

   // Standard destructor
   virtual ~TEmulatedCollectionProxy();

   // Virtual constructor
   virtual void* New()   const             {  return new Cont_t;         }

   // Virtual in-place constructor
   virtual void* New(void* memory)   const {  return new(memory) Cont_t; }

   // Virtual constructor
   virtual TClass::ObjectPtr NewObject()   const             {  return {new Cont_t, nullptr};         }

   // Virtual in-place constructor
   virtual TClass::ObjectPtr NewObject(void* memory)   const {  return {new(memory) Cont_t, nullptr}; }

   // Virtual array constructor
   virtual void* NewArray(Int_t nElements)   const             {  return new Cont_t[nElements];         }

   // Virtual in-place constructor
   virtual void* NewArray(Int_t nElements, void* memory)   const {  return new(memory) Cont_t[nElements]; }

   // Virtual array constructor
   virtual TClass::ObjectPtr NewObjectArray(Int_t nElements)   const             {  return {new Cont_t[nElements], nullptr};         }

   // Virtual in-place constructor
   virtual TClass::ObjectPtr NewObjectArray(Int_t nElements, void* memory)   const {  return {new(memory) Cont_t[nElements], nullptr}; }

   // Virtual destructor
   virtual void  Destructor(void* p, Bool_t dtorOnly = kFALSE) const;

   // Virtual array destructor
   virtual void  DeleteArray(void* p, Bool_t dtorOnly = kFALSE) const;

   // TVirtualCollectionProxy overload: Return the sizeof the collection object.
   virtual UInt_t Sizeof() const           {  return sizeof(Cont_t);     }

   // Return the address of the value at index 'idx'
   virtual void *At(UInt_t idx);

   // Clear the container
   virtual void Clear(const char *opt = "");

   // Resize the container
   virtual void Resize(UInt_t n, Bool_t force_delete);

   // Return the current size of the container
   virtual UInt_t Size() const;

   // Block allocation of containees
   virtual void* Allocate(UInt_t n, Bool_t forceDelete);

   // Block commit of containees
   virtual void Commit(void* env);

   // Insert data into the container where data is a C-style array of the actual type contained in the collection
   // of the given size.   For associative container (map, etc.), the data type is the pair<key,value>.
   virtual void  Insert(const void *data, void *container, size_t size);

   // Read portion of the streamer
   virtual void ReadBuffer(TBuffer &buff, void *pObj);
   virtual void ReadBuffer(TBuffer &buff, void *pObj, const TClass *onfile);

   // Streamer for I/O handling
   virtual void Streamer(TBuffer &refBuffer);

   // Streamer I/O overload
   virtual void Streamer(TBuffer &buff, void *pObj, int siz) {
      TGenCollectionProxy::Streamer(buff,pObj,siz);
   }

   // Check validity of the proxy itself
   Bool_t IsValid() const;
};

#endif
