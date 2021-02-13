// @(#)root/cont:$Id$
// Author: Philippe Canal 20/08/2010

/*************************************************************************
 * Copyright (C) 1995-2003, Rene Brun, Fons Rademakers and al.           *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TVirtualCollectionIterators
#define ROOT_TVirtualCollectionIterators

/**
\class TVirtualCollectionIterators
\ingroup IO
Small helper class to generically acquire and release iterators.
*/

#include "TVirtualCollectionProxy.h"
#include "TError.h"
#include <vector>

class TVirtualCollectionIterators
{
private:
   TVirtualCollectionIterators() = delete; // Intentionally unimplemented.
   TVirtualCollectionIterators(const TVirtualCollectionIterators&) = delete; // Intentionally unimplemented.

public:
   // Note when the collection is a vector, fBegin and fEnd points to
   // the start and end of the memory content rather than to the address
   // of iterators (saving one dereference when being used).

   typedef TVirtualCollectionProxy::CreateIterators_t CreateIterators_t;
   typedef TVirtualCollectionProxy::DeleteTwoIterators_t DeleteTwoIterators_t;

   char  fBeginBuffer[TVirtualCollectionProxy::fgIteratorArenaSize];
   char  fEndBuffer[TVirtualCollectionProxy::fgIteratorArenaSize];
   void *fBegin; // Pointer to the starting iterator (collection->begin())
   void *fEnd;   // Pointer to the ending iterator (collection->end())
   CreateIterators_t    fCreateIterators;
   DeleteTwoIterators_t fDeleteTwoIterators;

   TVirtualCollectionIterators(TVirtualCollectionProxy *proxy, Bool_t read_from_file = kTRUE) : fBegin( &(fBeginBuffer[0]) ), fEnd(&(fEndBuffer[0])), fCreateIterators(0), fDeleteTwoIterators(0)
   {
      // Constructor given a collection proxy.

      //         memset(fBeginBuffer,0,TVirtualCollectionProxy::fgIteratorArenaSize);
      //         memset(fEndBuffer,0,TVirtualCollectionProxy::fgIteratorArenaSize);
      if (proxy) {
         fCreateIterators = proxy->GetFunctionCreateIterators(read_from_file);
         fDeleteTwoIterators = proxy->GetFunctionDeleteTwoIterators(read_from_file);
      } else {
         ::Fatal("TIterators::TIterators","Created with out a collection proxy!\n");
      }
   }

   TVirtualCollectionIterators(CreateIterators_t creator, DeleteTwoIterators_t destruct) : fBegin( &(fBeginBuffer[0]) ), fEnd(&(fEndBuffer[0])), fCreateIterators(creator), fDeleteTwoIterators(destruct)
   {
      // Constructor given the creation and delete routines.
   }

   inline void CreateIterators(void *collection, TVirtualCollectionProxy *proxy)
   {
      // Initialize the fBegin and fEnd iterators.

      fCreateIterators(collection, &fBegin, &fEnd, proxy);
   }

   inline ~TVirtualCollectionIterators()
   {
      // Destructor.

      if (fBegin != &(fBeginBuffer[0])) {
         // assert(end != endbuf);
         fDeleteTwoIterators(fBegin,fEnd);
      }
   }
};


class TGenericCollectionIterator
{
protected:
   TVirtualCollectionIterators fIterators;

   // The actual implementation.
   class RegularIterator;
   class VectorIterator;

   TGenericCollectionIterator() = delete;
   TGenericCollectionIterator(const TGenericCollectionIterator&) = delete;

   TGenericCollectionIterator(void *collection, TVirtualCollectionProxy *proxy, Bool_t read_from_file = kTRUE) :
      fIterators(proxy,read_from_file)
   {
      // Regular constructor.

      fIterators.CreateIterators(collection,proxy);
   }

   virtual ~TGenericCollectionIterator()
   {
      // Regular destructor.
   }

public:

   virtual void *Next() = 0;

   virtual void* operator*() const = 0;

   virtual operator bool() const = 0;

   TGenericCollectionIterator& operator++() { Next(); return *this; }

   static TGenericCollectionIterator *New(void *collection, TVirtualCollectionProxy *proxy);
};

class TGenericCollectionIterator::RegularIterator : public TGenericCollectionIterator {
   typedef TVirtualCollectionProxy::Next_t Next_t;

   Next_t    fNext;
   void     *fCurrent;
   bool      fStarted : 1;

public:
   RegularIterator(void *collection, TVirtualCollectionProxy *proxy, Bool_t read_from_file) :
      TGenericCollectionIterator(collection,proxy,read_from_file),
      fNext( proxy->GetFunctionNext(read_from_file) ),
      fCurrent(0),
      fStarted(kFALSE)
   {
   }

   void *Next() {
      fStarted = kTRUE;
      fCurrent = fNext(fIterators.fBegin,fIterators.fEnd);
      return fCurrent;
   }

   virtual void* operator*() const { return fCurrent; }

   operator bool() const { return fStarted ? fCurrent != 0 : kTRUE; }

};

class TGenericCollectionIterator::VectorIterator : public TGenericCollectionIterator {

   ULong_t fIncrement;
   Bool_t  fHasPointer;

   inline void *GetValue() const {
      if ((bool)*this) return fHasPointer ? *(void**)fIterators.fBegin : fIterators.fBegin;
      else return 0;
   }

public:
   VectorIterator(void *collection, TVirtualCollectionProxy *proxy, Bool_t read_from_file) :
      TGenericCollectionIterator(collection,proxy,read_from_file),
      fIncrement(proxy->GetIncrement()),
      fHasPointer(proxy->HasPointers())
   {
   }

   void *Next() {
      if ( ! (bool)*this ) return 0;
      void *result = GetValue();
      fIterators.fBegin = ((char*)fIterators.fBegin) + fIncrement;
      return result;
   }

   virtual void* operator*() const { return GetValue(); }

   operator bool() const { return fIterators.fBegin != fIterators.fEnd; }

};

inline TGenericCollectionIterator *TGenericCollectionIterator::New(void *collection, TVirtualCollectionProxy *proxy)
{
   if (proxy->GetCollectionType() == ROOT::kSTLvector) {
      return new VectorIterator(collection, proxy, kFALSE);
   } else {
      return new RegularIterator(collection, proxy, kFALSE);
   }
}

/**
\class TVirtualCollectionPtrIterators
\ingroup IO
*/
class TVirtualCollectionPtrIterators
{
public:
   typedef TVirtualCollectionProxy::Next_t Next_t;
   typedef TVirtualCollectionProxy::CopyIterator_t Copy_t;
   typedef TVirtualCollectionProxy::CreateIterators_t CreateIterators_t;
   typedef TVirtualCollectionProxy::DeleteIterator_t Delete_t;
   typedef TVirtualCollectionProxy::DeleteTwoIterators_t DeleteTwoIterators_t;

private:
   TVirtualCollectionPtrIterators(); // Intentionally unimplemented.
   TVirtualCollectionPtrIterators(const TVirtualCollectionPtrIterators&); // Intentionally unimplemented.

   CreateIterators_t    fCreateIterators;
   DeleteTwoIterators_t fDeleteTwoIterators;

   Bool_t fAllocated;

   char  fRawBeginBuffer[TVirtualCollectionProxy::fgIteratorArenaSize];
   char  fRawEndBuffer[TVirtualCollectionProxy::fgIteratorArenaSize];

   struct TInternalIterator {
   private:
      TInternalIterator &operator=(const TInternalIterator&); // intentionally not implemented
   public:
      TInternalIterator() : fCopy(0),fDelete(0),fNext(0),fIter(0) {}
      TInternalIterator(const TInternalIterator &source) : fCopy(source.fCopy),fDelete(source.fDelete),fNext(source.fNext),fIter(0) {}

      Copy_t    fCopy;
      Delete_t  fDelete;
      Next_t    fNext;

      void     *fIter;
   };

   TInternalIterator fBeginBuffer;
   TInternalIterator fEndBuffer;

public:
   // Note when the collection is a vector, fBegin and fEnd points to
   // the start and end of the memory content rather than to the address
   // of iterators (saving one dereference when being used).

   void *fBegin; // Pointer to the starting iterator (collection->begin())
   void *fEnd;   // Pointer to the ending iterator (collection->end())

   TVirtualCollectionPtrIterators(TVirtualCollectionProxy *proxy) : fCreateIterators(0), fDeleteTwoIterators(0), fAllocated(kFALSE),
                                                                    fBegin( &(fRawBeginBuffer[0]) ),
                                                                    fEnd( &(fRawEndBuffer[0]) )
   {
      //         memset(fBeginBuffer,0,TVirtualCollectionProxy::fgIteratorArenaSize);
      //         memset(fEndBuffer,0,TVirtualCollectionProxy::fgIteratorArenaSize);
      if (proxy) {
         fCreateIterators = proxy->GetFunctionCreateIterators();
         fDeleteTwoIterators = proxy->GetFunctionDeleteTwoIterators();

         fEndBuffer.fCopy = fBeginBuffer.fCopy = proxy->GetFunctionCopyIterator();
         fEndBuffer.fNext = fBeginBuffer.fNext = proxy->GetFunctionNext();
         fEndBuffer.fDelete = fBeginBuffer.fDelete = proxy->GetFunctionDeleteIterator();
      } else {
         ::Fatal("TIterators::TIterators","Created with out a collection proxy!\n");
      }
   }

   inline void CreateIterators(void *collection, TVirtualCollectionProxy *proxy)
   {
      // Initialize the fBegin and fEnd iterators.

      fBegin = &(fRawBeginBuffer[0]);
      fEnd = &(fRawEndBuffer[0]);
      fCreateIterators(collection, &fBegin, &fEnd, proxy);
      if (fBegin != &(fRawBeginBuffer[0])) {
         // The iterator where too large to buffer in the  buffer
         fAllocated = kTRUE;
      }
      fBeginBuffer.fIter = fBegin;
      fEndBuffer.fIter = fEnd;
      fBegin = &fBeginBuffer;
      fEnd = &fEndBuffer;
   }

   inline ~TVirtualCollectionPtrIterators()
   {
      if (fAllocated) {
         // assert(end != endbuf);
         fDeleteTwoIterators(fBeginBuffer.fIter,fEndBuffer.fIter);
      }
   }

   static void *Next(void *iter, const void *end)
   {
      TInternalIterator *internal_iter = (TInternalIterator*) iter;
      TInternalIterator *internal_end = (TInternalIterator*) end;

      void **ptr = (void**)internal_iter->fNext(internal_iter->fIter,internal_end->fIter);
      if(ptr) return *ptr;
      else return 0;
   }

   static void DeleteIterator(void *iter)
   {
      TInternalIterator *internal_iter = (TInternalIterator*) iter;
      if (internal_iter->fDelete) {
         internal_iter->fDelete(internal_iter->fIter);
      }
   }

   static void *CopyIterator(void *dest, const void *source)
   {
      TInternalIterator *internal_source = (TInternalIterator*)source;
      TInternalIterator *internal_dest = new TInternalIterator(*internal_source);

      void *newiter = internal_source->fCopy(dest,internal_source->fIter);
      if (newiter == dest) {
         internal_dest->fDelete = 0;
      }
      internal_dest->fIter = newiter;
      return internal_dest;
   }
};

// Specialization of TVirtualCollectionIterators when we know the collection
// to be a vector (hence there is nothing to delete at the end).
struct TVirtualVectorIterators
{
private:
   TVirtualVectorIterators(const TVirtualVectorIterators&); // Intentionally unimplemented.

public:
   // Note when the collection is a vector, fBegin and fEnd points to
   // the start and end of the memory content rather than to the address
   // of iterators (saving one dereference when being used).

   typedef TVirtualCollectionProxy::CreateIterators_t CreateIterators_t;

   void *fBegin; // Pointer to the starting iterator (collection->begin())
   void *fEnd;   // Pointer to the ending iterator (collection->end())

   TVirtualVectorIterators(TVirtualCollectionProxy * /* proxy */) : fBegin(0), fEnd(0)
   {
      // fCreateIterators = proxy->GetFunctionCreateIterators();
   }

   TVirtualVectorIterators(CreateIterators_t /* creator */) : fBegin(0), fEnd(0)
   {
      // fCreateIterators = creator;
   }

   TVirtualVectorIterators() : fBegin(0), fEnd(0)
   {
      // Default constructor.
   }

   inline void CreateIterators(void *collection)
   {
      // Initialize the fBegin and fEnd iterators.

      // We can safely assume that the std::vector layout does not really depend on
      // the content!
      std::vector<char> *vec = (std::vector<char>*)collection;
      if (vec->empty()) {
         fBegin = 0;
         fEnd = 0;
         return;
      }
      fBegin= &(*vec->begin());
#ifdef R__VISUAL_CPLUSPLUS
      fEnd = &(*(vec->end()-1)) + 1; // On windows we can not dererence the end iterator at all.
#else
      // coverity[past_the_end] Safe on other platforms
      fEnd = &(*vec->end());
#endif
      //fCreateIterators(collection, &fBegin, &fEnd);
   }
};

#endif // ROOT_TVirtualCollectionIterators

