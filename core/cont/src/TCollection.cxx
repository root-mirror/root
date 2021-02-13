// @(#)root/cont:$Id$
// Author: Fons Rademakers   13/08/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TCollection
\ingroup Containers
Collection abstract base class. This class describes the base
protocol all collection classes have to implement. The ROOT
collection classes always store pointers to objects that inherit
from TObject. They never adopt the objects. Therefore, it is the
user's responsibility to take care of deleting the actual objects
once they are not needed anymore. In exceptional cases, when the
user is 100% sure nothing else is referencing the objects in the
collection, one can delete all objects and the collection at the
same time using the Delete() function.

Collections can be iterated using an iterator object (see
TIterator). Depending on the concrete collection class there may be
some additional methods of iterating. See the respective classes.

TCollection inherits from TObject since we want to be able to have
collections of collections.

In a later release the collections may become templatized.
*/

#include "TCollection.h"
#include "Varargs.h"
#include "TBuffer.h"
#include "TClass.h"
#include "TROOT.h"
#include "TBrowser.h"
#include "TObjectTable.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "TVirtualMutex.h"
#include "TError.h"
#include "TSystem.h"
#include "TObjArray.h"
#include <iostream>
#include <sstream>

#include "TSpinLockGuard.h"

TVirtualMutex *gCollectionMutex = nullptr;

TCollection   *TCollection::fgCurrentCollection = nullptr;
TObjectTable  *TCollection::fgGarbageCollection = nullptr;
Bool_t         TCollection::fgEmptyingGarbage   = kFALSE;
Int_t          TCollection::fgGarbageStack      = 0;

ClassImp(TCollection);
ClassImp(TIter);

#ifdef R__CHECK_COLLECTION_MULTI_ACCESS

void TCollection::TErrorLock::ConflictReport(std::thread::id holder, const char *accesstype,
                                             const TCollection *collection, const char *function)
{

   auto local = std::this_thread::get_id();
   std::stringstream cur, loc;
   if (holder == std::thread::id())
      cur << "None";
   else
      cur << "0x" << std::hex << holder;
   loc << "0x" << std::hex << local;

   //   std::cerr << "Error in " << function << ": Access (" << accesstype << ") to a collection (" <<
   //   collection->IsA()->GetName() << ":" << collection <<
   //   ") from multiple threads at a time. holder=" << "0x" << std::hex << holder << " readers=" << fReadSet.size() <<
   //   "0x" << std::hex << local << std::endl;

   ::Error(function,
           "Access (%s) to a collection (%s:%p) from multiple threads at a time. holder=%s readers=%lu intruder=%s",
           accesstype, collection->IsA()->GetName(), collection, cur.str().c_str(), fReadSet.size(), loc.str().c_str());

   std::set<std::thread::id> tmp;
   for (auto r : fReadSet) tmp.insert(r);
   for (auto r : tmp) {
      std::stringstream reader;
      reader << "0x" << std::hex << r;
      ::Error(function, " Readers includes %s", reader.str().c_str());
   }
   gSystem->StackTrace();
}

void TCollection::TErrorLock::Lock(const TCollection *collection, const char *function)
{
   auto local = std::this_thread::get_id();

   std::thread::id holder;

   if (fWriteCurrent.compare_exchange_strong(holder, local)) {
      // fWriteCurrent was the default id and is now local.
      ++fWriteCurrentRecurse;
      // std::cerr << "#" << "0x" << std::hex << local << " acquired first " << collection << " lock:" << this <<
      // std::endl;

      // Now check if there is any readers lingering
      if (fReadCurrentRecurse) {
         if (fReadSet.size() > 1 || fReadSet.find(local) != fReadSet.end()) {
            ConflictReport(std::thread::id(), "WriteLock while ReadLock taken", collection, function);
         }
      }
   } else {
      // fWriteCurrent was not the default id and is still the 'holder' thread id
      // this id is now also in the holder variable
      if (holder == local) {
         // The holder was actually this thread, no problem there, we
         // allow re-entrancy.
         // std::cerr << "#" << "0x" << std::hex << local << " re-entered " << fWriteCurrentRecurse << " " << collection
         // << " lock:" << this << std::endl;
      } else {
         ConflictReport(holder, "WriteLock", collection, function);
      }
      ++fWriteCurrentRecurse;
   }
}

void TCollection::TErrorLock::Unlock()
{
   auto local = std::this_thread::get_id();
   auto none = std::thread::id();

   --fWriteCurrentRecurse;
   if (fWriteCurrentRecurse == 0) {
      if (fWriteCurrent.compare_exchange_strong(local, none)) {
         // fWriteCurrent was local and is now none.

         // std::cerr << "#" << "0x" << std::hex << local << " zero and cleaned : " << std::dec << fWriteCurrentRecurse
         // << " 0x" << std::hex << fWriteCurrent.load() << " lock:" << this << std::endl;
      } else {
         // fWriteCurrent was not local, just live it as is.

         // std::cerr << "#" << "0x" << std::hex << local << " zero but somebody else : " << "0x" << std::hex <<
         // fWriteCurrent.load() << " lock:" << this << std::endl;
      }
   } else {
      // std::cerr << "#" << "0x" << std::hex << local << " still holding " << "0x" << std::hex << fWriteCurrentRecurse
      // << " lock:" << this << std::endl;
   }

   // std::cerr << "#" << "0x" << std::hex << local << " ended with : " << std::dec << fWriteCurrentRecurse << " 0x" <<
   // std::hex << fWriteCurrent.load() << " lock:" << this << std::endl;
}

void TCollection::TErrorLock::ReadLock(const TCollection *collection, const char *function)
{
   auto local = std::this_thread::get_id();

   {
      ROOT::Internal::TSpinLockGuard guard(fSpinLockFlag);
      fReadSet.insert(local); // this is not thread safe ...
   }
   ++fReadCurrentRecurse;

   if (fWriteCurrentRecurse) {
      auto holder = fWriteCurrent.load();
      if (holder != local) ConflictReport(holder, "ReadLock with WriteLock taken", collection, function);
   }
}

void TCollection::TErrorLock::ReadUnlock()
{
   auto local = std::this_thread::get_id();
   {
      ROOT::Internal::TSpinLockGuard guard(fSpinLockFlag);
      fReadSet.erase(local); // this is not thread safe ...
   }
   --fReadCurrentRecurse;
}

#endif // R__CHECK_COLLECTION_MULTI_ACCESS

////////////////////////////////////////////////////////////////////////////////
/// TNamed destructor.

TCollection::~TCollection()
{
   // Required since we overload TObject::Hash.
   ROOT::CallRecursiveRemoveIfNeeded(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Add all objects from collection col to this collection.

void TCollection::AddAll(const TCollection *col)
{
   TIter next(col);
   TObject *obj;

   while ((obj = next()))
      Add(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Add all arguments to the collection. The list of objects must be
/// terminated by 0, e.g.: l.AddVector(o1, o2, o3, o4, 0);

void TCollection::AddVector(TObject *va_(obj1), ...)
{
   va_list ap;
   va_start(ap, va_(obj1));
   TObject *obj;

   Add(va_(obj1));
   while ((obj = va_arg(ap, TObject *)))
      Add(obj);
   va_end(ap);
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure all objects in this collection inherit from class cl.

Bool_t TCollection::AssertClass(TClass *cl) const
{
   TObject *obj;
   TIter    next(this);
   Bool_t   error = kFALSE;

   if (!cl) {
      Error("AssertClass", "class == 0");
      return kTRUE;
   }

   for (int i = 0; (obj = next()); i++)
      if (!obj->InheritsFrom(cl)) {
         Error("AssertClass", "element %d is not an instance of class %s (%s)",
               i, cl->GetName(), obj->ClassName());
         error = kTRUE;
      }
   return error;
}

////////////////////////////////////////////////////////////////////////////////
/// Browse this collection (called by TBrowser).
/// If b=0, there is no Browse call TObject::Browse(0) instead.
///         This means TObject::Inspect() will be invoked indirectly

void TCollection::Browse(TBrowser *b)
{
   TIter next(this);
   TObject *obj;

   if (b)
      while ((obj = next())) b->Add(obj);
   else
      TObject::Browse(b);
}

////////////////////////////////////////////////////////////////////////////////
/// Make a clone of an collection using the Streamer facility.
/// If newname is specified, this will be the name of the new collection.

TObject *TCollection::Clone(const char *newname) const
{
   TCollection *new_collection = (TCollection*)TObject::Clone(newname);
   if (newname && strlen(newname)) new_collection->SetName(newname);
   return new_collection;
}


////////////////////////////////////////////////////////////////////////////////
/// Compare two TCollection objects. Returns 0 when equal, -1 when this is
/// smaller and +1 when bigger (like strcmp()).

Int_t TCollection::Compare(const TObject *obj) const
{
   if (this == obj) return 0;
   return fName.CompareTo(obj->GetName());
}

////////////////////////////////////////////////////////////////////////////////
/// Draw all objects in this collection.

void TCollection::Draw(Option_t *option)
{
   TIter next(this);
   TObject *object;

   while ((object = next())) {
      object->Draw(option);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Dump all objects in this collection.

void TCollection::Dump() const
{
   TIter next(this);
   TObject *object;

   while ((object = next())) {
      object->Dump();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Find an object in this collection using its name. Requires a sequential
/// scan till the object has been found. Returns 0 if object with specified
/// name is not found.

TObject *TCollection::FindObject(const char *name) const
{
   TIter next(this);
   TObject *obj;

   while ((obj = next()))
      if (!strcmp(name, obj->GetName())) return obj;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Find an object in this collection by name.

TObject *TCollection::operator()(const char *name) const
{
   return FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Find an object in this collection using the object's IsEqual()
/// member function. Requires a sequential scan till the object has
/// been found. Returns 0 if object is not found.
/// Typically this function is overridden by a more efficient version
/// in concrete collection classes (e.g. THashTable).

TObject *TCollection::FindObject(const TObject *obj) const
{
   TIter next(this);
   TObject *ob;

   while ((ob = next()))
      if (ob->IsEqual(obj)) return ob;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return name of this collection.
/// if no name, return the collection class name.

const char *TCollection::GetName() const
{
   if (fName.Length() > 0) return fName.Data();
   return ClassName();
}

////////////////////////////////////////////////////////////////////////////////
/// Increase the collection's capacity by delta slots.

Int_t TCollection::GrowBy(Int_t delta) const
{
   if (delta < 0) {
      Error("GrowBy", "delta < 0");
      delta = Capacity();
   }
   return Capacity() + TMath::Range(2, kMaxInt - Capacity(), delta);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true if object is a null pointer.

Bool_t  TCollection::IsArgNull(const char *where, const TObject *obj) const
{
   return obj ? kFALSE : (Error(where, "argument is a null pointer"), kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// List (ls) all objects in this collection.
/// Wildcarding supported, eg option="xxx*" lists only objects
/// with names xxx*.

void TCollection::ls(Option_t *option) const
{
   TROOT::IndentLevel();
   std::cout <<"OBJ: " << IsA()->GetName() << "\t" << GetName() << "\t" << GetTitle() << " : "
        << Int_t(TestBit(kCanDelete)) << std::endl;

   TRegexp re(option,kTRUE);
   TIter next(this);
   TObject *object;
   char *star = 0;
   if (option) star = (char*)strchr(option,'*');

   TROOT::IncreaseDirLevel();
   while ((object = next())) {
      if (star) {
         TString s = object->GetName();
         if (s != option && s.Index(re) == kNPOS) continue;
      }
      object->ls(option);
   }
   TROOT::DecreaseDirLevel();
}

////////////////////////////////////////////////////////////////////////////////
/// 'Notify' all objects in this collection.
Bool_t TCollection::Notify()
{
   Bool_t success = true;
   for (auto obj : *this) success &= obj->Notify();
   return success;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint all objects in this collection.

void TCollection::Paint(Option_t *option)
{
   this->R__FOR_EACH(TObject,Paint)(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Print the collection header.

void TCollection::PrintCollectionHeader(Option_t*) const
{
   TROOT::IndentLevel();
   printf("Collection name='%s', class='%s', size=%d\n",
          GetName(), ClassName(), GetSize());
}

////////////////////////////////////////////////////////////////////////////////
/// For given collection entry return the string that is used to
/// identify the object and, potentially, perform wildcard/regexp
/// filtering on.

const char* TCollection::GetCollectionEntryName(TObject* entry) const
{
   return entry->GetName();
}

////////////////////////////////////////////////////////////////////////////////
/// Print the collection entry.

void TCollection::PrintCollectionEntry(TObject* entry, Option_t* option, Int_t recurse) const
{
   TCollection* coll = dynamic_cast<TCollection*>(entry);
   if (coll) {
      coll->Print(option, recurse);
   } else {
      TROOT::IndentLevel();
      entry->Print(option);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Default print for collections, calls Print(option, 1).
/// This will print the collection header and Print() methods of
/// all the collection entries.
///
/// If you want to override Print() for a collection class, first
/// see if you can accomplish it by overriding the following protected
/// methods:
/// ~~~ {.cpp}
///   void        PrintCollectionHeader(Option_t* option) const;
///   const char* GetCollectionEntryName(TObject* entry) const;
///   void        PrintCollectionEntry(TObject* entry, Option_t* option, Int_t recurse) const;
/// ~~~
/// Otherwise override the `Print(Option_t *option, Int_t)`
/// variant. Remember to declare:
/// ~~~ {.cpp}
///   using TCollection::Print;
/// ~~~
/// somewhere close to the method declaration.

void TCollection::Print(Option_t *option) const
{
   Print(option, 1);
}

////////////////////////////////////////////////////////////////////////////////
/// Print the collection header and its elements.
///
/// If recurse is non-zero, descend into printing of
/// collection-entries with recurse - 1.
/// This means, if recurse is negative, the recursion is infinite.
///
/// Option is passed recursively.

void TCollection::Print(Option_t *option, Int_t recurse) const
{
   PrintCollectionHeader(option);

   if (recurse != 0)
   {
      TIter next(this);
      TObject *object;

      TROOT::IncreaseDirLevel();
      while ((object = next())) {
         PrintCollectionEntry(object, option, recurse - 1);
      }
      TROOT::DecreaseDirLevel();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print the collection header and its elements that match the wildcard.
///
/// If recurse is non-zero, descend into printing of
/// collection-entries with recurse - 1.
/// This means, if recurse is negative, the recursion is infinite.
///
/// Option is passed recursively, but wildcard is only used on the
/// first level.

void TCollection::Print(Option_t *option, const char* wildcard, Int_t recurse) const
{
   PrintCollectionHeader(option);

   if (recurse != 0)
   {
      if (!wildcard) wildcard = "";
      TRegexp re(wildcard, kTRUE);
      Int_t nch = strlen(wildcard);
      TIter next(this);
      TObject *object;

      TROOT::IncreaseDirLevel();
      while ((object = next())) {
         TString s = GetCollectionEntryName(object);
         if (nch == 0 || s == wildcard || s.Index(re) != kNPOS) {
            PrintCollectionEntry(object, option, recurse - 1);
         }
      }
      TROOT::DecreaseDirLevel();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print the collection header and its elements that match the regexp.
///
/// If recurse is non-zero, descend into printing of
/// collection-entries with recurse - 1.
/// This means, if recurse is negative, the recursion is infinite.
///
/// Option is passed recursively, but regexp is only used on the
/// first level.

void TCollection::Print(Option_t *option, TPRegexp& regexp, Int_t recurse) const
{
   PrintCollectionHeader(option);

   if (recurse != 0)
   {
      TIter next(this);
      TObject *object;

      TROOT::IncreaseDirLevel();
      while ((object = next())) {
         TString s = GetCollectionEntryName(object);
         if (regexp.MatchB(s)) {
            PrintCollectionEntry(object, option, recurse - 1);
         }
      }
      TROOT::DecreaseDirLevel();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Remove object from this collection and recursively remove the object
/// from all other objects (and collections).

void TCollection::RecursiveRemove(TObject *obj)
{
   if (!obj) return;

   // Scan list and remove obj in the list itself
   while (Remove(obj))
      ;

   // Scan again the list and invoke RecursiveRemove for all objects
   TIter next(this);
   TObject *object;

   while ((object = next())) {
      if (object->TestBit(kNotDeleted)) object->RecursiveRemove(obj);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Remove all objects in collection col from this collection.

void TCollection::RemoveAll(TCollection *col)
{
   TIter next(col);
   TObject *obj;

   while ((obj = next()))
      Remove(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Stream all objects in the collection to or from the I/O buffer.

void TCollection::Streamer(TBuffer &b)
{
   Int_t nobjects;
   TObject *obj;
   UInt_t R__s, R__c;

   if (b.IsReading()) {
      Version_t v = b.ReadVersion(&R__s, &R__c);
      if (v > 2)
         TObject::Streamer(b);
      if (v > 1)
         fName.Streamer(b);
      b >> nobjects;
      for (Int_t i = 0; i < nobjects; i++) {
         b >> obj;
         Add(obj);
      }
      b.CheckByteCount(R__s, R__c,TCollection::IsA());
   } else {
      R__c = b.WriteVersion(TCollection::IsA(), kTRUE);
      TObject::Streamer(b);
      fName.Streamer(b);
      nobjects = GetSize();
      b << nobjects;

      TIter next(this);

      while ((obj = next())) {
         b << obj;
      }
      b.SetByteCount(R__c, kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write all objects in this collection. By default all objects in
/// the collection are written individually (each object gets its
/// own key). Note, this is recursive, i.e. objects in collections
/// in the collection are also written individually. To write all
/// objects using a single key specify a name and set option to
/// TObject::kSingleKey (i.e. 1).

Int_t TCollection::Write(const char *name, Int_t option, Int_t bsize) const
{
   if ((option & kSingleKey)) {
      return TObject::Write(name, option, bsize);
   } else {
      option &= ~kSingleKey;
      Int_t nbytes = 0;
      TIter next(this);
      TObject *obj;
      while ((obj = next())) {
         nbytes += obj->Write(name, option, bsize);
      }
      return nbytes;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write all objects in this collection. By default all objects in
/// the collection are written individually (each object gets its
/// own key). Note, this is recursive, i.e. objects in collections
/// in the collection are also written individually. To write all
/// objects using a single key specify a name and set option to
/// TObject::kSingleKey (i.e. 1).

Int_t TCollection::Write(const char *name, Int_t option, Int_t bsize)
{
   return ((const TCollection*)this)->Write(name,option,bsize);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the globally accessible collection.

TCollection *TCollection::GetCurrentCollection()
{
   return fgCurrentCollection;
}

////////////////////////////////////////////////////////////////////////////////
/// Set this collection to be the globally accessible collection.

void TCollection::SetCurrentCollection()
{
   fgCurrentCollection = this;
}

////////////////////////////////////////////////////////////////////////////////
/// Set up for garbage collection.

void TCollection::StartGarbageCollection()
{
   R__LOCKGUARD2(gCollectionMutex);
   if (!fgGarbageCollection) {
      fgGarbageCollection = new TObjectTable;
      fgEmptyingGarbage   = kFALSE;
      fgGarbageStack      = 0;
   }
   fgGarbageStack++;
}

////////////////////////////////////////////////////////////////////////////////
/// Do the garbage collection.

void TCollection::EmptyGarbageCollection()
{
   R__LOCKGUARD2(gCollectionMutex);
   if (fgGarbageStack > 0) fgGarbageStack--;
   if (fgGarbageCollection && fgGarbageStack == 0 && fgEmptyingGarbage == kFALSE) {
      fgEmptyingGarbage = kTRUE;
      fgGarbageCollection->Delete();
      fgEmptyingGarbage = kFALSE;
      SafeDelete(fgGarbageCollection);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add to the list of things to be cleaned up.

void TCollection::GarbageCollect(TObject *obj)
{
   {
      R__LOCKGUARD2(gCollectionMutex);
      if (fgGarbageCollection) {
         if (!fgEmptyingGarbage) {
            fgGarbageCollection->Add(obj);
            return;
         }
      }
   }
   delete obj;
}

////////////////////////////////////////////////////////////////////////////////
/// Set whether this collection is the owner (enable==true)
/// of its content.  If it is the owner of its contents,
/// these objects will be deleted whenever the collection itself
/// is deleted. The objects might also be deleted or destructed when Clear
/// is called (depending on the collection).

void TCollection::SetOwner(Bool_t enable)
{
   if (enable)
      SetBit(kIsOwner);
   else
      ResetBit(kIsOwner);
}

////////////////////////////////////////////////////////////////////////////////
/// Set this collection to use a RW lock upon access, making it thread safe.
/// Return the previous state.
///
/// Note: To test whether the usage is enabled do:
///    collection->TestBit(TCollection::kUseRWLock);

bool TCollection::UseRWLock()
{
   bool prev = TestBit(TCollection::kUseRWLock);
   SetBit(TCollection::kUseRWLock);
   return prev;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy a TIter. This involves allocating a new TIterator of the right
/// sub class and assigning it with the original.

TIter::TIter(const TIter &iter)
{
   if (iter.fIterator) {
      fIterator = iter.GetCollection()->MakeIterator();
      fIterator->operator=(*iter.fIterator);
   } else
      fIterator = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Assigning an TIter to another. This involves allocating a new TIterator
/// of the right sub class and assigning it with the original.

TIter &TIter::operator=(const TIter &rhs)
{
   if (this != &rhs) {
      if (rhs.fIterator) {
         delete fIterator;
         fIterator = rhs.GetCollection()->MakeIterator();
         fIterator->operator=(*rhs.fIterator);
      }
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Pointing to the first element of the container.

TIter &TIter::Begin()
{
   fIterator->Reset();
   fIterator->Next();
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Pointing to the element after the last - to a nullptr value in our case.

TIter TIter::End()
{
   return TIter(static_cast<TIterator*>(nullptr));
}

////////////////////////////////////////////////////////////////////////////////
/// Return an empty collection for use with nullptr TRangeCast

const TCollection &ROOT::Internal::EmptyCollection()
{
   static TObjArray sEmpty;
   return sEmpty;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if 'cl' inherits from 'base'.

bool ROOT::Internal::ContaineeInheritsFrom(TClass *cl, TClass *base)
{
   return cl->InheritsFrom(base);
}
