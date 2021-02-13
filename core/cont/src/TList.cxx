// @(#)root/cont:$Id$
// Author: Fons Rademakers   10/08/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TList
\ingroup Containers
A doubly linked list.

All classes inheriting from TObject can be
inserted in a TList. Before being inserted into the list the object
pointer is wrapped in a TObjLink object which contains, besides
the object pointer also a previous and next pointer.

There are several ways to iterate over a TList; in order of preference, if
not forced by other constraints:

  1. (Preferred way) Using the C++ range-based `for` or `begin()` / `end()`:
~~~ {.cpp}
         for(const auto&& obj: *GetListOfPrimitives())
            obj->Write();
~~~
  2. Using the R__FOR_EACH macro:
~~~ {.cpp}
         GetListOfPrimitives()->R__FOR_EACH(TObject,Paint)(option);
~~~
  3. Using the TList iterator TListIter (via the wrapper class TIter):
~~~ {.cpp}
         TIter next(GetListOfPrimitives());
         while (TObject *obj = next())
            obj->Draw(next.GetOption());
~~~
  4. Using the TList iterator TListIter and std::for_each algorithm:
~~~ {.cpp}
         // A function object, which will be applied to each element
         // of the given range.
         struct STestFunctor {
            bool operator()(TObject *aObj) {
               ...
               return true;
            }
         }
         ...
         ...
         TIter iter(mylist);
         for_each( iter.Begin(), TIter::End(), STestFunctor() );
~~~
  5. Using the TObjLink list entries (that wrap the TObject*):
~~~ {.cpp}
         TObjLink *lnk = GetListOfPrimitives()->FirstLink();
         while (lnk) {
            lnk->GetObject()->Draw(lnk->GetOption());
            lnk = lnk->Next();
         }
~~~
  6. Using the TList's After() and Before() member functions:
~~~ {.cpp}
         TFree *idcur = this;
         while (idcur) {
            ...
            ...
            idcur = (TFree*)GetListOfFree()->After(idcur);
         }
~~~
Methods 3, 4 and 5 can also easily iterate backwards using either
a backward TIter (using argument kIterBackward) or by using
LastLink() and lnk->Prev() or by using the Before() member.
*/

#include "TList.h"
#include "TClass.h"
#include "TROOT.h"
#include "TVirtualMutex.h"
#include "TBuffer.h"

#include <string>

using namespace std;

ClassImp(TList);

////////////////////////////////////////////////////////////////////////////////
/// Delete the list. Objects are not deleted unless the TList is the
/// owner (set via SetOwner()).

TList::~TList()
{
   Clear();
}

////////////////////////////////////////////////////////////////////////////////
/// Add object at the beginning of the list.

void TList::AddFirst(TObject *obj)
{
   R__COLLECTION_WRITE_GUARD();

   if (IsArgNull("AddFirst", obj)) return;

   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   if (!fFirst) {
      fFirst = NewLink(obj);
      fLast = fFirst;
   } else {
      auto t = NewLink(obj);
      t->fNext = fFirst;
      fFirst->fPrev = t;
      fFirst = t;
   }
   fSize++;
   Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Add object at the beginning of the list and also store option.
/// Storing an option is useful when one wants to change the behaviour
/// of an object a little without having to create a complete new
/// copy of the object. This feature is used, for example, by the Draw()
/// method. It allows the same object to be drawn in different ways.

void TList::AddFirst(TObject *obj, Option_t *opt)
{
   R__COLLECTION_WRITE_GUARD();

   if (IsArgNull("AddFirst", obj)) return;

   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   if (!fFirst) {
      fFirst = NewOptLink(obj, opt);
      fLast = fFirst;
   } else {
      auto t = NewOptLink(obj, opt);
      t->fNext = fFirst;
      fFirst->fPrev = t;
      fFirst = t;
   }
   fSize++;
   Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Add object at the end of the list.

void TList::AddLast(TObject *obj)
{
   R__COLLECTION_WRITE_GUARD();

   if (IsArgNull("AddLast", obj)) return;

   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   if (!fFirst) {
      fFirst = NewLink(obj);
      fLast  = fFirst;
   } else
      fLast = NewLink(obj, fLast);
   fSize++;
   Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Add object at the end of the list and also store option.
/// Storing an option is useful when one wants to change the behaviour
/// of an object a little without having to create a complete new
/// copy of the object. This feature is used, for example, by the Draw()
/// method. It allows the same object to be drawn in different ways.

void TList::AddLast(TObject *obj, Option_t *opt)
{
   R__COLLECTION_WRITE_GUARD();

   if (IsArgNull("AddLast", obj)) return;

   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   if (!fFirst) {
      fFirst = NewOptLink(obj, opt);
      fLast  = fFirst;
   } else
      fLast = NewOptLink(obj, opt, fLast);
   fSize++;
   Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Insert object before object before in the list.

void TList::AddBefore(const TObject *before, TObject *obj)
{
   R__COLLECTION_WRITE_GUARD();

   if (IsArgNull("AddBefore", obj)) return;

   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   if (!before)
      TList::AddFirst(obj);
   else {
      Int_t    idx;
      TObjLink *t = FindLink(before, idx);
      if (!t) {
         Error("AddBefore", "before not found, object not added");
         return;
      }
      if (t == fFirst.get())
         TList::AddFirst(obj);
      else {
         NewLink(obj, t->fPrev.lock());
         fSize++;
         Changed();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Insert object before the specified ObjLink object. If before = 0 then add
/// to the head of the list. An ObjLink can be obtained by looping over a list
/// using the above describe iterator method 3.

void TList::AddBefore(TObjLink *before, TObject *obj)
{
   R__COLLECTION_WRITE_GUARD();

   if (IsArgNull("AddBefore", obj)) return;

   if (!before)
      TList::AddFirst(obj);
   else {
      if (before == fFirst.get())
         TList::AddFirst(obj);
      else {
         NewLink(obj, before->fPrev.lock());
         fSize++;
         Changed();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Insert object after object after in the list.

void TList::AddAfter(const TObject *after, TObject *obj)
{
   R__COLLECTION_WRITE_GUARD();

   if (IsArgNull("AddAfter", obj)) return;

   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   if (!after)
      TList::AddLast(obj);
   else {
      Int_t    idx;
      TObjLink *t = FindLink(after, idx);
      if (!t) {
         Error("AddAfter", "after not found, object not added");
         return;
      }
      if (t == fLast.get())
         TList::AddLast(obj);
      else {
         NewLink(obj, t->shared_from_this());
         fSize++;
         Changed();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Insert object after the specified ObjLink object. If after = 0 then add
/// to the tail of the list. An ObjLink can be obtained by looping over a list
/// using the above describe iterator method 3.

void TList::AddAfter(TObjLink *after, TObject *obj)
{
   R__COLLECTION_WRITE_GUARD();

   if (IsArgNull("AddAfter", obj)) return;

   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   if (!after)
      TList::AddLast(obj);
   else {
      if (after == fLast.get())
         TList::AddLast(obj);
      else {
         NewLink(obj, after->shared_from_this());
         fSize++;
         Changed();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Insert object at position idx in the list.

void TList::AddAt(TObject *obj, Int_t idx)
{
   R__COLLECTION_WRITE_GUARD();

   if (IsArgNull("AddAt", obj)) return;

   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   TObjLink *lnk = LinkAt(idx);
   if (!lnk)
      TList::AddLast(obj);
   else if (lnk == fFirst.get())
      TList::AddFirst(obj);
   else {
      NewLink(obj, lnk->fPrev.lock());
      fSize++;
      Changed();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the object after object obj. Obj is found using the
/// object's IsEqual() method.  Returns 0 if obj is last in list.

TObject *TList::After(const TObject *obj) const
{
   R__COLLECTION_WRITE_GUARD();

   TObjLink *t;

   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   auto cached = fCache.lock();
   if (cached.get() && cached->GetObject() && cached->GetObject()->IsEqual(obj)) {
      t = cached.get();
      ((TList*)this)->fCache = cached->fNext;  //cast const away, fCache should be mutable
   } else {
      Int_t idx;
      t = FindLink(obj, idx);
      if (t) ((TList*)this)->fCache = t->fNext;
   }

   if (t && t->Next())
      return t->Next()->GetObject();
   else
      return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the object at position idx. Returns 0 if idx is out of range.

TObject *TList::At(Int_t idx) const
{
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);
   R__COLLECTION_WRITE_GUARD();

   TObjLink *lnk = LinkAt(idx);
   if (lnk) return lnk->GetObject();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the object before object obj. Obj is found using the
/// object's IsEqual() method.  Returns 0 if obj is first in list.

TObject *TList::Before(const TObject *obj) const
{
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);
   R__COLLECTION_WRITE_GUARD();

   TObjLink *t;

   auto cached = fCache.lock();
   if (cached.get() && cached->GetObject() && cached->GetObject()->IsEqual(obj)) {
      t = cached.get();
      ((TList*)this)->fCache = cached->fPrev;  //cast const away, fCache should be mutable
   } else {
      Int_t idx;
      t = FindLink(obj, idx);
      if (t) ((TList*)this)->fCache = t->fPrev;
   }

   if (t && t->Prev())
      return t->Prev()->GetObject();
   else
      return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove all objects from the list. Does not delete the objects
/// unless the TList is the owner (set via SetOwner()) and option
/// "nodelete" is not set.
/// If option="nodelete" then don't delete any heap objects that were
/// marked with the kCanDelete bit, otherwise these objects will be
/// deleted (this option is used by THashTable::Clear()).

void TList::Clear(Option_t *option)
{
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);
   R__COLLECTION_WRITE_GUARD();

   Bool_t nodel = option ? (!strcmp(option, "nodelete") ? kTRUE : kFALSE) : kFALSE;

   if (!nodel && IsOwner()) {
      Delete(option);
      return;
   }

   // In some case, for example TParallelCoord, a list (the pad's list of
   // primitives) will contain both the container and the containees
   // (the TParallelCoordVar) but if the Clear is being called from
   // the destructor of the container of this list, one of the first
   // thing done will be the remove the container (the pad) for the
   // list (of Primitives of the canvas) that was connecting it
   // (indirectly) to the list of cleanups.
   // Note: The Code in TParallelCoordVar was changed (circa June 2017),
   // to no longer have this behavior and thus rely on this code (by moving
   // from using Draw to Paint) but the structure might still exist elsewhere
   // so we keep this comment here.

   // To preserve this connection (without introducing one when there was none),
   // we re-use fCache to inform RecursiveRemove of the node currently
   // being cleared/deleted.
   while (fFirst) {
      auto tlk = fFirst;
      fFirst = fFirst->fNext;
      fSize--;


      // Make node available to RecursiveRemove
      tlk->fNext.reset();
      tlk->fPrev.reset();
      fCache = tlk;

      // delete only heap objects marked OK to clear
      auto obj = tlk->GetObject();
      if (!nodel && obj) {
         if (!obj->TestBit(kNotDeleted)) {
            Error("Clear", "A list is accessing an object (%p) already deleted (list name = %s)",
                  obj, GetName());
         } else if (obj->IsOnHeap()) {
            if (obj->TestBit(kCanDelete)) {
               if (obj->TestBit(kNotDeleted)) {
                  TCollection::GarbageCollect(obj);
               }
            }
         }
      }
      // delete tlk;
   }
   fFirst.reset();
   fLast.reset();
   fCache.reset();
   fSize = 0;
   Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Remove all objects from the list AND delete all heap based objects.
/// If option="slow" then keep list consistent during delete. This allows
/// recursive list operations during the delete (e.g. during the dtor
/// of an object in this list one can still access the list to search for
/// other not yet deleted objects).

void TList::Delete(Option_t *option)
{
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);
   R__COLLECTION_WRITE_GUARD();

   Bool_t slow = option ? (!strcmp(option, "slow") ? kTRUE : kFALSE) : kFALSE;

   TList removeDirectory; // need to deregister these from their directory

   if (slow) {

      // In some case, for example TParallelCoord, a list (the pad's list of
      // primitives) will contain both the container and the containees
      // (the TParallelCoorVar) but if the Clear is being called from
      // the destructor of the container of this list, one of the first
      // thing done will be the remove the container (the pad) for the
      // list (of Primitives of the canvas) that was connecting it
      // (indirectly) to the list of cleanups.

      // To preserve this connection (without introducing one when there was none),
      // we re-use fCache to inform RecursiveRemove of the node currently
      // being cleared/deleted.
      while (fFirst) {
         auto tlk = fFirst;
         fFirst = fFirst->fNext;
         fSize--;

         // Make node available to RecursiveRemove
         tlk->fNext.reset();
         tlk->fPrev.reset();
         fCache = tlk;

         // delete only heap objects
         auto obj = tlk->GetObject();
         if (obj && !obj->TestBit(kNotDeleted))
            Error("Delete", "A list is accessing an object (%p) already deleted (list name = %s)",
                  obj, GetName());
         else if (obj && obj->IsOnHeap())
            TCollection::GarbageCollect(obj);
         else if (obj && obj->IsA()->GetDirectoryAutoAdd())
            removeDirectory.Add(obj);

         // delete tlk;
      }

      fFirst.reset();
      fLast.reset();
      fCache.reset();
      fSize  = 0;

   } else {

      auto first = fFirst;    //pointer to first entry in linked list
      fFirst.reset();
      fLast.reset();
      fCache.reset();
      fSize  = 0;
      while (first) {
         auto tlk = first;
         first = first->fNext;
         // delete only heap objects
         auto obj = tlk->GetObject();
         tlk->SetObject(nullptr);
         if (obj && !obj->TestBit(kNotDeleted))
            Error("Delete", "A list is accessing an object (%p) already deleted (list name = %s)",
                  obj, GetName());
         else if (obj && obj->IsOnHeap())
            TCollection::GarbageCollect(obj);
         else if (obj && obj->IsA()->GetDirectoryAutoAdd())
            removeDirectory.Add(obj);

         // The formerly first token, when tlk goes out of scope has no more references
         // because of the fFirst.reset()
      }
   }

   // These objects cannot expect to have a valid TDirectory anymore;
   // e.g. because *this is the TDirectory's list of objects. Even if
   // not, they are supposed to be deleted, so we can as well unregister
   // them from their directory, even if they are stack-based:
   TIter iRemDir(&removeDirectory);
   TObject* dirRem = 0;
   while ((dirRem = iRemDir())) {
      (*dirRem->IsA()->GetDirectoryAutoAdd())(dirRem, 0);
   }
   Changed();
}

#if 0
////////////////////////////////////////////////////////////////////////////////
/// Delete a TObjLink object.

void TList::DeleteLink(TObjLink *lnk)
{
   R__COLLECTION_WRITE_GUARD();

   lnk->fNext = lnk->fPrev = 0;
   lnk->fObject = 0;
   delete lnk;
}
#endif

////////////////////////////////////////////////////////////////////////////////
/// Find an object in this list using its name. Requires a sequential
/// scan till the object has been found. Returns 0 if object with specified
/// name is not found. This method overrides the generic FindObject()
/// of TCollection for efficiency reasons.

TObject *TList::FindObject(const char *name) const
{
   R__COLLECTION_READ_GUARD();

   if (!name)
      return nullptr;

   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

   for (TObjLink *lnk = FirstLink(); lnk != nullptr; lnk = lnk->Next()) {
      if (TObject *obj = lnk->GetObject()) {
         const char *objname = obj->GetName();
         if (objname && strcmp(name, objname) == 0)
            return obj;
      }
   }

   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Find an object in this list using the object's IsEqual()
/// member function. Requires a sequential scan till the object has
/// been found. Returns 0 if object is not found.
/// This method overrides the generic FindObject() of TCollection for
/// efficiency reasons.

TObject *TList::FindObject(const TObject *obj) const
{
   R__COLLECTION_READ_GUARD();

   if (!obj)
      return nullptr;

   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

   TObjLink *lnk = FirstLink();

   while (lnk) {
      TObject *ob = lnk->GetObject();
      if (ob->IsEqual(obj)) return ob;
      lnk = lnk->Next();
   }
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the TObjLink object that contains object obj. In idx it returns
/// the position of the object in the list.

TObjLink *TList::FindLink(const TObject *obj, Int_t &idx) const
{
   R__COLLECTION_READ_GUARD();

   if (!obj)
      return nullptr;

   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

   if (!fFirst) return 0;

   TObject *object;
   TObjLink *lnk = fFirst.get();
   idx = 0;

   while (lnk) {
      object = lnk->GetObject();
      if (object) {
         if (object->TestBit(kNotDeleted)) {
            if (object->IsEqual(obj)) return lnk;
         }
      }
      lnk = lnk->Next();
      idx++;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the first object in the list. Returns 0 when list is empty.

TObject *TList::First() const
{
   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);
   R__COLLECTION_READ_GUARD();

   if (fFirst) return fFirst->GetObject();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return address of pointer to obj

TObject **TList::GetObjectRef(const TObject *obj) const
{
   R__COLLECTION_READ_GUARD();

   if (!obj)
   return nullptr;

   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

   TObjLink *lnk = FirstLink();

   while (lnk) {
      TObject *ob = lnk->GetObject();
      if (ob->IsEqual(obj)) return lnk->GetObjectRef();
      lnk = lnk->Next();
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the last object in the list. Returns 0 when list is empty.

TObject *TList::Last() const
{
   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);
   R__COLLECTION_READ_GUARD();

   if (fLast) return fLast->GetObject();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the TObjLink object at index idx.

TObjLink *TList::LinkAt(Int_t idx) const
{
   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);
   R__COLLECTION_READ_GUARD();

   Int_t    i = 0;
   TObjLink *lnk = fFirst.get();
   while (i < idx && lnk) {
      i++;
      lnk = lnk->Next();
   }
   return lnk;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a list iterator.

TIterator *TList::MakeIterator(Bool_t dir) const
{
   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);
   R__COLLECTION_READ_GUARD();

   return new TListIter(this, dir);
}

////////////////////////////////////////////////////////////////////////////////
/// Return a new TObjLink.

TList::TObjLinkPtr_t TList::NewLink(TObject *obj, const TObjLinkPtr_t &prev)
{
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);
   R__COLLECTION_WRITE_GUARD();

   auto newlink = std::make_shared<TObjLink>(obj);
   if (prev) {
      InsertAfter(newlink, prev);
   }
   return newlink;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a new TObjOptLink (a TObjLink that also stores the option).

TList::TObjLinkPtr_t TList::NewOptLink(TObject *obj, Option_t *opt, const TObjLinkPtr_t &prev)
{
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);
   R__COLLECTION_WRITE_GUARD();

   auto newlink = std::make_shared<TObjOptLink>(obj, opt);
   if (prev) {
      InsertAfter(newlink, prev);
   }
   return newlink;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove object from this collection and recursively remove the object
/// from all other objects (and collections).

void TList::RecursiveRemove(TObject *obj)
{
   R__COLLECTION_WRITE_GUARD();

   if (!obj) return;

   // When fCache is set and has no previous and next node, it represents
   // the node being cleared and/or deleted.
   {
      auto cached = fCache.lock();
      if (cached && cached->fNext.get() == nullptr && cached->fPrev.lock().get() == nullptr) {
         TObject *ob = cached->GetObject();
         if (ob && ob->TestBit(kNotDeleted)) {
            ob->RecursiveRemove(obj);
         }
      }
   }

   if (!fFirst.get())
      return;

   auto lnk  = fFirst;
   decltype(lnk) next;
   while (lnk.get()) {
      next = lnk->fNext;
      TObject *ob = lnk->GetObject();
      if (ob && ob->TestBit(kNotDeleted)) {
         if (ob->IsEqual(obj)) {
            lnk->SetObject(nullptr);
            if (lnk == fFirst) {
               fFirst = next;
               if (lnk == fLast)
                  fLast = fFirst;
               else
                  fFirst->fPrev.reset();
               // DeleteLink(lnk);
            } else if (lnk == fLast) {
               fLast = lnk->fPrev.lock();
               fLast->fNext.reset();
               // DeleteLink(lnk);
            } else {
               lnk->Prev()->fNext = next;
               lnk->Next()->fPrev = lnk->fPrev;
               // DeleteLink(lnk);
            }
            fSize--;
            fCache.reset();
            Changed();
         } else
            ob->RecursiveRemove(obj);
      }
      lnk = next;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Remove object from the list.

TObject *TList::Remove(TObject *obj)
{
   R__COLLECTION_WRITE_GUARD();

   if (!obj) return 0;

   Int_t    idx;
   TObjLink *lnk = FindLink(obj, idx);

   if (!lnk) return 0;

   // return object found, which may be (pointer wise) different than the
   // input object (depending on what IsEqual() is doing)

   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   TObject *ob = lnk->GetObject();
   lnk->SetObject(nullptr);
   if (lnk == fFirst.get()) {
      fFirst = lnk->fNext;
      // lnk is still alive as we have either fLast
      // or the 'new' fFirst->fPrev pointing to it.
      if (lnk == fLast.get()) {
         fLast.reset();
         fFirst.reset();
      } else
         fFirst->fPrev.reset();
      //DeleteLink(lnk);
   } else if (lnk == fLast.get()) {
      fLast = lnk->fPrev.lock();
      fLast->fNext.reset();
      //DeleteLink(lnk);
   } else {
      lnk->Next()->fPrev = lnk->fPrev;
      lnk->Prev()->fNext = lnk->fNext;
      //DeleteLink(lnk);
   }
   fSize--;
   fCache.reset();
   Changed();

   return ob;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove object link (and therefore the object it contains)
/// from the list.

TObject *TList::Remove(TObjLink *lnk)
{
   R__COLLECTION_WRITE_GUARD();

   if (!lnk) return 0;

   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   TObject *obj = lnk->GetObject();
   lnk->SetObject(nullptr);
   if (lnk == fFirst.get()) {
      fFirst = lnk->fNext;
      // lnk is still alive as we have either fLast
      // or the 'new' fFirst->fPrev pointing to it.
      if (lnk == fLast.get()) {
         fLast.reset();
         fFirst.reset();
      } else
         fFirst->fPrev.reset();
      // DeleteLink(lnk);
   } else if (lnk == fLast.get()) {
      fLast = lnk->fPrev.lock();
      fLast->fNext.reset();
      // DeleteLink(lnk);
   } else {
      lnk->Next()->fPrev = lnk->fPrev;
      lnk->Prev()->fNext = lnk->fNext;
      // DeleteLink(lnk);
   }
   fSize--;
   fCache.reset();
   Changed();

   return obj;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove the last object of the list.

void TList::RemoveLast()
{
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);
   R__COLLECTION_WRITE_GUARD();

   TObjLink *lnk = fLast.get();
   if (!lnk) return;

   lnk->SetObject(nullptr);
   if (lnk == fFirst.get()) {
      fFirst.reset();
      fLast.reset();
   } else {
      fLast = lnk->fPrev.lock();
      fLast->fNext.reset();
   }
   // DeleteLink(lnk);

   fSize--;
   fCache.reset();
   Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Sort linked list. Real sorting is done in private function DoSort().
/// The list can only be sorted when is contains objects of a sortable
/// class.

void TList::Sort(Bool_t order)
{
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);
   R__COLLECTION_WRITE_GUARD();

   if (!fFirst) return;

   fAscending = order;

   if (!fFirst->GetObject()->IsSortable()) {
      Error("Sort", "objects in list are not sortable");
      return;
   }

   DoSort(&fFirst, fSize);

   // correct back links
   std::shared_ptr<TObjLink> ol, lnk = fFirst;

   if (lnk.get()) lnk->fPrev.reset();
   while ((ol = lnk)) {
      lnk = lnk->fNext;
      if (lnk)
         lnk->fPrev = ol;
      else
         fLast = ol;
   }
   fSorted = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Compares the objects stored in the TObjLink objects.
/// Depending on the flag IsAscending() the function returns
/// true if the object in l1 <= l2 (ascending) or l2 <= l1 (descending).

Bool_t TList::LnkCompare(const TObjLinkPtr_t &l1, const TObjLinkPtr_t &l2)
{
   Int_t cmp = l1->GetObject()->Compare(l2->GetObject());

   if ((IsAscending() && cmp <=0) || (!IsAscending() && cmp > 0))
      return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Sort linked list.

std::shared_ptr<TObjLink> *TList::DoSort(std::shared_ptr<TObjLink> *head, Int_t n)
{
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);
   R__COLLECTION_WRITE_GUARD();

   std::shared_ptr<TObjLink> p1, p2, *h2, *t2;

   switch (n) {
      case 0:
         return head;

      case 1:
         return &((*head)->fNext);

      case 2:
         p2 = (p1 = *head)->fNext;
         if (LnkCompare(p1, p2)) return &(p2->fNext);
         p1->fNext = (*head = p2)->fNext;
         return &((p2->fNext = p1)->fNext);
   }

   int m;
   n -= m = n / 2;

   t2 = DoSort(h2 = DoSort(head, n), m);

   if (LnkCompare((p1 = *head), (p2 = *h2))) {
      do {
         if (!--n) return *h2 = p2, t2;
      } while (LnkCompare((p1 = *(head = &(p1->fNext))), p2));
   }

   while (1) {
      *head = p2;
      do {
         if (!--m) return *h2 = *t2, *t2 = p1, h2;
      } while (!LnkCompare(p1, (p2 = *(head = &(p2->fNext)))));
      *head = p1;
      do {
         if (!--n) return *h2 = p2, t2;
      } while (LnkCompare((p1 = *(head = &(p1->fNext))), p2));
   }
}

/** \class TObjLink
Wrapper around a TObject so it can be stored in a TList.
*/

////////////////////////////////////////////////////////////////////////////////
/// Insert a new link in the chain.

void TList::InsertAfter(const TObjLinkPtr_t &newlink, const TObjLinkPtr_t &prev)
{
   newlink->fNext = prev->fNext;
   newlink->fPrev = prev;
   prev->fNext = newlink;
   if (newlink->fNext)
      newlink->fNext->fPrev = newlink;
}

/** \class TListIter
Iterator of linked list.
*/

ClassImp(TListIter);

////////////////////////////////////////////////////////////////////////////////
/// Create a new list iterator. By default the iteration direction
/// is kIterForward. To go backward use kIterBackward.

TListIter::TListIter(const TList *l, Bool_t dir)
        : fList(l), fCurCursor(0), fCursor(0), fDirection(dir), fStarted(kFALSE)
{
   R__COLLECTION_ITER_GUARD(fList);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy ctor.

TListIter::TListIter(const TListIter &iter) : TIterator(iter)
{
   R__COLLECTION_ITER_GUARD(iter.fList);

   fList      = iter.fList;
   fCurCursor = iter.fCurCursor;
   fCursor    = iter.fCursor;
   fDirection = iter.fDirection;
   fStarted   = iter.fStarted;
}

////////////////////////////////////////////////////////////////////////////////
/// Overridden assignment operator.

TIterator &TListIter::operator=(const TIterator &rhs)
{

   const TListIter *rhs1 = dynamic_cast<const TListIter *>(&rhs);
   if (this != &rhs && rhs1) {
      R__COLLECTION_ITER_GUARD(rhs1->fList);
      TIterator::operator=(rhs);
      fList      = rhs1->fList;
      fCurCursor = rhs1->fCurCursor;
      fCursor    = rhs1->fCursor;
      fDirection = rhs1->fDirection;
      fStarted   = rhs1->fStarted;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Overloaded assignment operator.

TListIter &TListIter::operator=(const TListIter &rhs)
{
   if (this != &rhs) {
      R__COLLECTION_ITER_GUARD(rhs.fList);
      TIterator::operator=(rhs);
      fList      = rhs.fList;
      fCurCursor = rhs.fCurCursor;
      fCursor    = rhs.fCursor;
      fDirection = rhs.fDirection;
      fStarted   = rhs.fStarted;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Return next object in the list. Returns 0 when no more objects in list.

TObject *TListIter::Next()
{
   if (!fList) return 0;

   R__COLLECTION_ITER_GUARD(fList);

   if (fDirection == kIterForward) {
      if (!fStarted) {
         fCursor = fList->fFirst;
         fStarted = kTRUE;
      }
      fCurCursor = fCursor;
      if (fCursor) {
         auto next = fCursor = fCursor->NextSP();
      }
   } else {
      if (!fStarted) {
         fCursor = fList->fLast;
         fStarted = kTRUE;
      }
      fCurCursor = fCursor;
      if (fCursor) fCursor = fCursor->PrevSP();
   }

   if (fCurCursor) return fCurCursor->GetObject();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the object option stored in the list.

Option_t *TListIter::GetOption() const
{
   if (fCurCursor) return fCurCursor->GetOption();
   return "";
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the object option stored in the list.

void TListIter::SetOption(Option_t *option)
{
   if (fCurCursor) fCurCursor->SetOption(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset list iterator.

void TListIter::Reset()
{
   R__COLLECTION_ITER_GUARD(fList);
   fStarted = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// This operator compares two TIterator objects.

Bool_t TListIter::operator!=(const TIterator &aIter) const
{
   if (IsA() == aIter.IsA()) {
      // We compared equal only two iterator of the same type.
      // Since this is a function of TListIter, we consequently know that
      // both this and aIter are of type inheriting from TListIter.
      const TListIter &iter(dynamic_cast<const TListIter &>(aIter));
      return (fCurCursor != iter.fCurCursor);
   }
   return false; // for base class we don't implement a comparison
}

////////////////////////////////////////////////////////////////////////////////
/// This operator compares two TListIter objects.

Bool_t TListIter::operator!=(const TListIter &aIter) const
{
   return (fCurCursor != aIter.fCurCursor);
}

////////////////////////////////////////////////////////////////////////////////
/// Stream all objects in the collection to or from the I/O buffer.

void TList::Streamer(TBuffer &b)
{
   R__COLLECTION_WRITE_GUARD();

   Int_t nobjects;
   UChar_t nch;
   Int_t nbig;
   TObject *obj;
   UInt_t R__s, R__c;

   if (b.IsReading()) {
      Clear(); // Get rid of old data if any.
      Version_t v = b.ReadVersion(&R__s, &R__c);
      if (v > 3) {
         TObject::Streamer(b);
         fName.Streamer(b);
         b >> nobjects;
         string readOption;
         for (Int_t i = 0; i < nobjects; i++) {
            b >> obj;
            b >> nch;
            if (v > 4 && nch == 255)  {
               b >> nbig;
            } else {
               nbig = nch;
            }
            readOption.resize(nbig,'\0');
            b.ReadFastArray((char*) readOption.data(),nbig);
            if (obj) { // obj can be null if the class had a custom streamer and we do not have the shared library nor a streamerInfo.
               if (nch) {
                  Add(obj,readOption.c_str());
               } else {
                  Add(obj);
               }
            }
         }
         b.CheckByteCount(R__s, R__c,TList::IsA());
         return;
      }

      //  process old versions when TList::Streamer was in TCollection::Streamer
      if (v > 2)
         TObject::Streamer(b);
      if (v > 1)
         fName.Streamer(b);
      b >> nobjects;
      for (Int_t i = 0; i < nobjects; i++) {
         b >> obj;
         Add(obj);
      }
      b.CheckByteCount(R__s, R__c,TList::IsA());

   } else {
      R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

      R__c = b.WriteVersion(TList::IsA(), kTRUE);
      TObject::Streamer(b);
      fName.Streamer(b);
      nobjects = GetSize();
      b << nobjects;

      TObjLink *lnk = fFirst.get();
      while (lnk) {
         obj = lnk->GetObject();
         b << obj;

         nbig = strlen(lnk->GetAddOption());
         if (nbig > 254) {
            nch = 255;
            b << nch;
            b << nbig;
         } else {
            nch = UChar_t(nbig);
            b << nch;
         }
         b.WriteFastArray(lnk->GetAddOption(),nbig);

         lnk = lnk->Next();
      }
      b.SetByteCount(R__c, kTRUE);
   }
}
