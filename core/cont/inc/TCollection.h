// @(#)root/cont:$Id$
// Author: Fons Rademakers   13/08/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TCollection
#define ROOT_TCollection


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCollection                                                          //
//                                                                      //
// Collection abstract base class. This class inherits from TObject     //
// because we want to be able to have collections of collections.       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"

#include "TIterator.h"

#include "TString.h"

#include "TVirtualRWMutex.h"

#include <cassert>

class TClass;
class TObjectTable;
class TVirtualMutex;
class TIter;

const Bool_t kIterForward  = kTRUE;
const Bool_t kIterBackward = !kIterForward;

R__EXTERN TVirtualMutex *gCollectionMutex;

// #define R__CHECK_COLLECTION_MULTI_ACCESS

// When R__CHECK_COLLECTION_MULTI_ACCESS is turned on (defined),
// the normal (not locked) ROOT TCollections are instrumented with a
// pseudo read-write lock which does not halt the execution but detects
// and report concurrent access to the same collections.
// Multiple readers are allowed.
// Multiple concurrent writer is reported as a Conflict
// Readers access while a write is running is reported as Conflict
// Re-entrant writing call by the same Writer thread are allowed.
// Entering a writing section by a single Reader thread is allowed.

#ifdef R__CHECK_COLLECTION_MULTI_ACCESS
#include <atomic>
#include <thread>
#include <unordered_multiset>
#endif

class TCollection : public TObject {

#ifdef R__CHECK_COLLECTION_MULTI_ACCESS
public:
   class TErrorLock {
      // Warn when multiple thread try to acquire the same 'lock'
      std::atomic<std::thread::id> fWriteCurrent;
      std::atomic<size_t> fWriteCurrentRecurse;
      std::atomic<size_t> fReadCurrentRecurse;
      std::unordered_multiset<std::thread::id> fReadSet;
      std::atomic_flag fSpinLockFlag;

      void Lock(const TCollection *collection, const char *function);

      void Unlock();

      void ReadLock(const TCollection *collection, const char *function);

      void ReadUnlock();

      void ConflictReport(std::thread::id holder, const char *accesstype, const TCollection *collection,
                          const char *function);

   public:
      TErrorLock() : fWriteCurrent(), fWriteCurrentRecurse(0), fReadCurrentRecurse(0)
      {
         std::atomic_flag_clear(&fSpinLockFlag);
      }

      class WriteGuard {
         TErrorLock *fLock;

      public:
         WriteGuard(TErrorLock &lock, const TCollection *collection, const char *function) : fLock(&lock)
         {
            fLock->Lock(collection, function);
         }
         ~WriteGuard() { fLock->Unlock(); }
      };

      class ReadGuard {
         TErrorLock *fLock;

      public:
         ReadGuard(TErrorLock &lock, const TCollection *collection, const char *function) : fLock(&lock)
         {
            fLock->ReadLock(collection, function);
         }
         ~ReadGuard() { fLock->ReadUnlock(); }
      };
   };

   mutable TErrorLock fLock; //! Special 'lock' to detect multiple access to a collection.

#define R__COLLECTION_WRITE_GUARD() TCollection::TErrorLock::WriteGuard wg(fLock, this, __PRETTY_FUNCTION__)
#define R__COLLECTION_READ_GUARD() TCollection::TErrorLock::ReadGuard rg(fLock, this, __PRETTY_FUNCTION__)

#define R__COLLECTION_ITER_GUARD(collection) \
   TCollection::TErrorLock::ReadGuard rg(collection->fLock, collection, __PRETTY_FUNCTION__)

#else

#define R__COLLECTION_WRITE_GUARD()
#define R__COLLECTION_READ_GUARD()
#define R__COLLECTION_ITER_GUARD(collection)

#endif

private:
   static TCollection  *fgCurrentCollection;  //used by macro R__FOR_EACH
   static TObjectTable *fgGarbageCollection;  //used by garbage collector
   static Bool_t        fgEmptyingGarbage;    //used by garbage collector
   static Int_t         fgGarbageStack;       //used by garbage collector

   TCollection(const TCollection &);    //private and not-implemented, collections
   void operator=(const TCollection &); //are too complex to be automatically copied

protected:
   enum EStatusBits {
      kIsOwner   = BIT(14),
      // BIT(15) is used by TClonesArray and TMap
      kUseRWLock = BIT(16)
   };

   TString   fName;               //name of the collection
   Int_t     fSize;               //number of elements in collection

   TCollection() : fName(), fSize(0) { }

   virtual void        PrintCollectionHeader(Option_t* option) const;
   virtual const char* GetCollectionEntryName(TObject* entry) const;
   virtual void        PrintCollectionEntry(TObject* entry, Option_t* option, Int_t recurse) const;

public:
   enum { kInitCapacity = 16, kInitHashTableCapacity = 17 };

   virtual            ~TCollection();
   virtual void       Add(TObject *obj) = 0;
   void               AddVector(TObject *obj1, ...);
   virtual void       AddAll(const TCollection *col);
   Bool_t             AssertClass(TClass *cl) const;
   void               Browse(TBrowser *b);
   Int_t              Capacity() const { return fSize; }
   virtual void       Clear(Option_t *option="") = 0;
   virtual TObject   *Clone(const char *newname="") const;
   Int_t              Compare(const TObject *obj) const;
   Bool_t             Contains(const char *name) const { return FindObject(name) != 0; }
   Bool_t             Contains(const TObject *obj) const { return FindObject(obj) != 0; }
   virtual void       Delete(Option_t *option="") = 0;
   virtual void       Draw(Option_t *option="");
   virtual void       Dump() const ;
   virtual TObject   *FindObject(const char *name) const;
   TObject           *operator()(const char *name) const;
   virtual TObject   *FindObject(const TObject *obj) const;
   virtual Int_t      GetEntries() const { return GetSize(); }
   virtual const char *GetName() const;
   virtual TObject  **GetObjectRef(const TObject *obj) const = 0;
   /// Return the *capacity* of the collection, i.e. the current total amount of space that has been allocated so far.
   /// Same as `Capacity`. Use `GetEntries` to get the number of elements currently in the collection.
   virtual Int_t      GetSize() const { return fSize; }
   virtual Int_t      GrowBy(Int_t delta) const;
   ULong_t            Hash() const { return fName.Hash(); }
   Bool_t             IsArgNull(const char *where, const TObject *obj) const;
   virtual Bool_t     IsEmpty() const { return GetSize() <= 0; }
   virtual Bool_t     IsFolder() const { return kTRUE; }
   Bool_t             IsOwner() const { return TestBit(kIsOwner); }
   Bool_t             IsSortable() const { return kTRUE; }
   virtual void       ls(Option_t *option="") const ;
   virtual Bool_t     Notify();
   virtual TIterator *MakeIterator(Bool_t dir = kIterForward) const = 0;
   virtual TIterator *MakeReverseIterator() const { return MakeIterator(kIterBackward); }
   virtual void       Paint(Option_t *option="");
   virtual void       Print(Option_t *option="") const;
   virtual void       Print(Option_t *option, Int_t recurse) const;
   virtual void       Print(Option_t *option, const char* wildcard, Int_t recurse=1) const;
   virtual void       Print(Option_t *option, TPRegexp& regexp, Int_t recurse=1) const;
   virtual void       RecursiveRemove(TObject *obj);
   virtual TObject   *Remove(TObject *obj) = 0;
   virtual void       RemoveAll(TCollection *col);
   void               RemoveAll() { Clear(); }
   void               SetCurrentCollection();
   void               SetName(const char *name) { fName = name; }
   virtual void       SetOwner(Bool_t enable = kTRUE);
   virtual bool       UseRWLock();
   virtual Int_t      Write(const char *name=0, Int_t option=0, Int_t bufsize=0);
   virtual Int_t      Write(const char *name=0, Int_t option=0, Int_t bufsize=0) const;

   R__ALWAYS_INLINE Bool_t IsUsingRWLock() const { return TestBit(TCollection::kUseRWLock); }

   static TCollection  *GetCurrentCollection();
   static void          StartGarbageCollection();
   static void          GarbageCollect(TObject *obj);
   static void          EmptyGarbageCollection();

   TIter begin() const;
   TIter end() const;

   ClassDef(TCollection,3)  //Collection abstract base class
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TIter                                                                //
//                                                                      //
// Iterator wrapper. Type of iterator used depends on type of           //
// collection.                                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TIter {

private:
   TIterator    *fIterator;         //collection iterator

protected:
   TIter() : fIterator(nullptr) { }

public:
   TIter(const TCollection *col, Bool_t dir = kIterForward)
         : fIterator(col ? col->MakeIterator(dir) : 0) { }
   TIter(TIterator *it) : fIterator(it) { }
   TIter(const TIter &iter);
   TIter &operator=(const TIter &rhs);
   virtual ~TIter() { SafeDelete(fIterator); }
   TObject           *operator()() { return Next(); }
   TObject           *Next() { return fIterator ? fIterator->Next() : nullptr; }
   const TCollection *GetCollection() const { return fIterator ? fIterator->GetCollection() : nullptr; }
   Option_t          *GetOption() const { return fIterator ? fIterator->GetOption() : ""; }
   void               Reset() { if (fIterator) fIterator->Reset(); }
   TIter             &operator++() { Next(); return *this; }
   Bool_t             operator==(const TIter &aIter) const {
      if (fIterator == nullptr)
         return aIter.fIterator == nullptr || **aIter.fIterator == nullptr;
      if (aIter.fIterator == nullptr)
         return fIterator == nullptr || **fIterator == nullptr;
      return *fIterator == *aIter.fIterator;
   }
   Bool_t             operator!=(const TIter &aIter) const {
      return !(*this == aIter);
   }
   TObject           *operator*() const { return fIterator ? *(*fIterator): nullptr; }
   TIter             &Begin();
   static TIter       End();

   ClassDef(TIter,0)  //Iterator wrapper
};

template <class T>
class TIterCategory: public TIter, public std::iterator_traits<typename T::Iterator_t> {

public:
   TIterCategory(const TCollection *col, Bool_t dir = kIterForward) : TIter(col, dir) { }
   TIterCategory(TIterator *it) : TIter(it) { }
   virtual ~TIterCategory() { }
   TIterCategory &Begin() { TIter::Begin(); return *this; }
   static TIterCategory End() { return TIterCategory(static_cast<TIterator*>(nullptr)); }
};


inline TIter TCollection::begin() const { return ++(TIter(this)); }
inline TIter TCollection::end() const { return TIter::End(); }

namespace ROOT {
namespace Internal {

const TCollection &EmptyCollection();
bool ContaineeInheritsFrom(TClass *cl, TClass *base);

/// @brief Internal help class implmenting an iterator for TRangeDynCast.
template <class Containee> // Containee must derive from TObject.
class TRangeDynCastIterator : public TIter {
   static_assert(std::is_base_of<TObject, Containee>::value, "Containee type must inherit from TObject");

   /// This is a workaround against ClassDefInline not supporting classes
   /// missing their default constructor or having them private.
   template <class T>
   friend class ROOT::Internal::ClassDefGenerateInitInstanceLocalInjector;

   TRangeDynCastIterator() = default;

public:
   using TIter::TIter;
   TRangeDynCastIterator(const TIter &iter) : TIter(iter) {}

   Containee *operator()() = delete;

   Containee *Next() { return dynamic_cast<Containee *>(TIter::Next()); }
   Containee *operator*() const { return dynamic_cast<Containee *>(TIter::operator*()); }

   ClassDefInline(TRangeDynCastIterator, 0);
};

} // namespace Internal

namespace Detail {

/// @brief TTypedIter is a typed version of TIter.
///
/// This requires the collection to contains elements of the type requested
/// (or a derived class).  Any deviation from this expectation
/// will only be caught/reported by an assert in debug builds.
///
/// This is best used with a TClonesArray, for other cases prefered TRangeDynCast.
///
/// The typical use is:
/// ```{.cpp}
///    TTypedIter<TBaseClass> next(tbaseClassClonesArrayPtr);
///    while(auto bcl = next()) {
///       ... use bcl as a TBaseClass*
///    }
/// ```
template <class Containee> // Containee must derive from TObject.
class TTypedIter : public TIter {
   static_assert(std::is_base_of<TObject, Containee>::value, "Containee type must inherit from TObject");

   /// This is a workaround against ClassDefInline not supporting classes
   /// missing their default constructor or having them private.
   template <class T>
   friend class ROOT::Internal::ClassDefGenerateInitInstanceLocalInjector;

   TTypedIter() = default;

   static Containee *StaticCast(TObject *obj)
   {
      assert(!obj || ROOT::Internal::ContaineeInheritsFrom(obj->IsA(), Containee::Class()));
      return static_cast<Containee *>(obj);
   }

public:
   using TIter::TIter;
   TTypedIter(const TIter &iter) : TIter(iter) {}

   Containee *operator()() { return StaticCast(TIter::Next()); }
   Containee *Next() { return StaticCast(TIter::Next()); }
   Containee *operator*() const { return StaticCast(TIter::operator*()); }

   ClassDefInline(TTypedIter, 0);
};

/// @brief TRangeStaticCast is an adaptater class that allows the typed iteration
/// through a TCollection.  This requires the collection to contains element
/// of the type requested (or a derived class).  Any deviation from this expectation
/// will only be caught/reported by an assert in debug builds.
///
/// This is best used with a TClonesArray, for other cases prefered TRangeDynCast.
///
/// The typical use is:
/// ```{.cpp}
///    for(auto bcl : TRangeStaticCast<TBaseClass>( *tbaseClassClonesArrayPtr )) {
///        assert(bcl && bcl->IsA()->InheritsFrom(TBaseClass::Class()));
///        ... use bcl as a TBaseClass*
///    }
///    for(auto bcl : TRangeStaticCast<TBaseClass>( tbaseClassClonesArrayPtr )) {
///        assert(bcl && bcl->IsA()->InheritsFrom(TBaseClass::Class()));
///        ... use bcl as a TBaseClass*
///    }
/// ```
template <class T>
class TRangeStaticCast {
   const TCollection &fCollection;

public:
   TRangeStaticCast(const TCollection &col) : fCollection(col) {}
   TRangeStaticCast(const TCollection *col) : fCollection(col != nullptr ? *col : ROOT::Internal::EmptyCollection()) {}

   TTypedIter<T> begin() const { return fCollection.begin(); }
   TTypedIter<T> end() const { return fCollection.end(); }
};

} // namespace Detail
} // namespace ROOT

/// @brief TRangeDynCast is an adaptater class that allows the typed iteration
/// through a TCollection.
///
/// The typical use is:
/// ```{.cpp}
///    for(auto bcl : TRangeDynCast<TBaseClass>( *cl->GetListOfBases() )) {
///        if (!bcl) continue;
///        ... use bcl as a TBaseClass*
///    }
///    for(auto bcl : TRangeDynCast<TBaseClass>( cl->GetListOfBases() )) {
///        if (!bcl) continue;
///        ... use bcl as a TBaseClass*
///    }
/// ```
template <class T>
class TRangeDynCast {
   const TCollection &fCollection;

public:
   TRangeDynCast(const TCollection &col) : fCollection(col) {}
   TRangeDynCast(const TCollection *col) : fCollection(col != nullptr ? *col : ROOT::Internal::EmptyCollection()) {}

   ROOT::Internal::TRangeDynCastIterator<T> begin() const { return fCollection.begin(); }
   ROOT::Internal::TRangeDynCastIterator<T> end() const { return fCollection.end(); }
};

// Zero overhead macros in case not compiled with thread support
#if defined (_REENTRANT) || defined (WIN32)

#define R__COLL_COND_MUTEX(mutex) this->IsUsingRWLock() ? mutex : nullptr

#define R__COLLECTION_READ_LOCKGUARD(mutex) ::ROOT::TReadLockGuard _R__UNIQUE_(R__readguard)(R__COLL_COND_MUTEX(mutex))
#define R__COLLECTION_READ_LOCKGUARD_NAMED(name,mutex) ::ROOT::TReadLockGuard _NAME2_(R__readguard,name)(R__COLL_COND_MUTEX(mutex))

#define R__COLLECTION_WRITE_LOCKGUARD(mutex) ::ROOT::TWriteLockGuard _R__UNIQUE_(R__readguard)(R__COLL_COND_MUTEX(mutex))
#define R__COLLECTION_WRITE_LOCKGUARD_NAMED(name,mutex) ::ROOT::TWriteLockGuard _NAME2_(R__readguard,name)(R__COLL_COND_MUTEX(mutex))

#else

#define R__COLLECTION_READ_LOCKGUARD(mutex) (void)mutex
#define R__COLLECTION_COLLECTION_READ_LOCKGUARD_NAMED(name,mutex) (void)mutex

#define R__COLLECTION_WRITE_LOCKGUARD(mutex) (void)mutex
#define R__COLLECTION_WRITE_LOCKGUARD_NAMED(name,mutex) (void)mutex

#endif

//---- R__FOR_EACH macro -------------------------------------------------------

// Macro to loop over all elements of a list of type "type" while executing
// procedure "proc" on each element

#define R__FOR_EACH(type,proc) \
    SetCurrentCollection(); \
    TIter _NAME3_(nxt_,type,proc)(TCollection::GetCurrentCollection()); \
    type *_NAME3_(obj_,type,proc); \
    while ((_NAME3_(obj_,type,proc) = (type*) _NAME3_(nxt_,type,proc)())) \
       _NAME3_(obj_,type,proc)->proc

#endif
