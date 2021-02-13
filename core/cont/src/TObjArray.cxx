// @(#)root/cont:$Id$
// Author: Fons Rademakers   11/09/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TObjArray
\ingroup Containers
An array of TObjects. The array expands automatically when
objects are added (shrinking can be done by hand using Expand(),
how nice to have meaningful names -:)).
Use operator[] to have "real" array behaviour.

Note on ownership and copy:
By default the TObjArray does not own the objects it points to and
will not delete them unless explicitly asked (via a call to the
Delete member function).   To assign ownership of the content to
the array, call:
~~~ {.cpp}
     myarr->SetOwner(kTRUE);
~~~
When the array owns its content a call to Clear or the deletion of
the array itself will lead to the deletion of its contents.

You can either make a shallow copy of the array:
~~~ {.cpp}
     otherarr = new TObjArray(*myarr);
    *otherarr = *myarr;
~~~
in which case ownership (if any) is not transfered but the other
array points to the same object as the original array.  Note that
if the content of either array is deleted the other array is not
notified in any way (i.e. still points to the now deleted objects).

You can also make a deep copy of the array:
~~~ {.cpp}
     otherarr = (TObjArray*)myarr->Clone();
~~~
in which case the array and the content are both duplicated (i.e.
otherarr and myarr do not point to the same objects).  If myarr
is set to the be the owner of its content, otherarr will also be
set to the owner of its own content.
*/

#include "TObjArray.h"
#include "TError.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TVirtualMutex.h"
#include <stdlib.h>

ClassImp(TObjArray);

////////////////////////////////////////////////////////////////////////////////
/// Create an object array. Using s one can set the array size (default is
/// kInitCapacity=16) and lowerBound can be used to set the array lowerbound
/// index (default is 0).

TObjArray::TObjArray(Int_t s, Int_t lowerBound)
{
   if (s < 0) {
      Warning("TObjArray", "size (%d) < 0", s);
      s = TCollection::kInitCapacity;
   } else if (s == 0)
      s = TCollection::kInitCapacity;
   fCont = 0;
   Init(s, lowerBound);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a copy of TObjArray a. Note, does not copy the kIsOwner flag.

TObjArray::TObjArray(const TObjArray &a) : TSeqCollection()
{
   fCont = 0;
   Init(a.fSize, a.fLowerBound);

   for (Int_t i = 0; i < fSize; i++)
      fCont[i] = a.fCont[i];

   fLast = a.fLast;
   fName = a.fName;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete an array. Objects are not deleted unless the TObjArray is the
/// owner (set via SetOwner()).

TObjArray::~TObjArray()
{
   if (IsOwner())
      Delete();

   TStorage::Dealloc(fCont);
   fCont = 0;
   fSize = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator. Note, unsets the kIsOwner flag.

TObjArray& TObjArray::operator=(const TObjArray &a)
{
   if (this != &a) {
      R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

      if (IsOwner())
         Delete();
      SetOwner(kFALSE);

      Init(a.fSize, a.fLowerBound);

      for (Int_t i = 0; i < fSize; i++)
         fCont[i] = a.fCont[i];

      fLast = a.fLast;
      fName = a.fName;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the object at position i. Returns address at position 0
/// if i is out of bounds. Result may be used as an lvalue.

TObject *&TObjArray::operator[](Int_t i)
{
   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

   int j = i-fLowerBound;
   if (j >= 0 && j < fSize) {
      fLast = TMath::Max(j, GetAbsLast());
      Changed();
      return fCont[j];
   }
   BoundsOk("operator[]", i);
   fLast = -2; // invalidate fLast since the result may be used as an lvalue
   return fCont[0];
}

////////////////////////////////////////////////////////////////////////////////
/// Return the object at position at. Returns 0 if i is out of bounds.

TObject *TObjArray::operator[](Int_t i) const
{
   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

   int j = i-fLowerBound;
   if (j >= 0 && j < fSize) return fCont[j];
   BoundsOk("operator[] const", i);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Add object in the first slot of the array. This will overwrite the
/// first element that might have been there. To have insertion semantics
/// use either a TList or a TOrdCollection.

void TObjArray::AddFirst(TObject *obj)
{
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   fCont[0] = obj;
   if (fLast == -1)
      fLast = 0;
   Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Add object in the next empty slot in the array. Expand the array
/// if necessary.

void TObjArray::AddLast(TObject *obj)
{
   AddAtAndExpand(obj, GetAbsLast()+1+fLowerBound);
}

////////////////////////////////////////////////////////////////////////////////
/// Add object in the slot before object before. If before=0 add object
/// in the first slot. Note that this will overwrite any object that
/// might have already been in this slot. For insertion semantics use
/// either a TList or a TOrdCollection.

void TObjArray::AddBefore(const TObject *before, TObject *obj)
{
   if (!before)
      AddFirst(obj);
   else {
      R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

      Int_t idx = IndexOf(before) - fLowerBound;
      if (idx == -1) {
         Error("AddBefore", "before not found, object not added");
         return;
      }
      if (idx == 0) {
         Error("AddBefore", "cannot add before lowerbound (%d)", fLowerBound);
         return;
      }
      AddAt(obj, idx+fLowerBound-1);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add object in the slot after object after. If after=0 add object in
/// the last empty slot. Note that this will overwrite any object that
/// might have already been in this slot. For insertion semantics use
/// either a TList or a TOrdCollection.

void TObjArray::AddAfter(const TObject *after, TObject *obj)
{
   if (!after)
      AddLast(obj);
   else {
      R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

      Int_t idx = IndexOf(after) - fLowerBound;
      if (idx == -1) {
         Error("AddAfter", "after not found, object not added");
         return;
      }
      AddAtAndExpand(obj, idx+fLowerBound+1);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add object at position idx. If idx is larger than the current size
/// of the array, expand the array (double its size).

void TObjArray::AddAtAndExpand(TObject *obj, Int_t idx)
{
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   if (idx < fLowerBound) {
      Error("AddAt", "out of bounds at %d in %lx", idx, (Long_t)this);
      return;
   }
   if (idx-fLowerBound >= fSize)
      Expand(TMath::Max(idx-fLowerBound+1, GrowBy(fSize)));
   fCont[idx-fLowerBound] = obj;
   fLast = TMath::Max(idx-fLowerBound, GetAbsLast());
   Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Add object at position ids. Give an error when idx is out of bounds
/// (i.e. the array is not expanded).

void TObjArray::AddAt(TObject *obj, Int_t idx)
{
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   if (!BoundsOk("AddAt", idx)) return;

   fCont[idx-fLowerBound] = obj;
   fLast = TMath::Max(idx-fLowerBound, GetAbsLast());
   Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the position of the new object.
/// Find the first empty cell or AddLast if there is no empty cell

Int_t  TObjArray::AddAtFree(TObject *obj)
{
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   if (Last()) {    // <---------- This is to take in account "empty" TObjArray's
      Int_t i;
      for (i = 0; i < fSize; i++)
         if (!fCont[i]) {         // Add object at position i
            fCont[i] = obj;
            fLast = TMath::Max(i, GetAbsLast());
            Changed();
            return i+fLowerBound;
         }
   }
   AddLast(obj);
   return GetLast();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the object after obj. Returns 0 if obj is last object.

TObject *TObjArray::After(const TObject *obj) const
{
   if (!obj) return 0;

   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

   Int_t idx = IndexOf(obj) - fLowerBound;
   if (idx == -1 || idx == fSize-1) return 0;

   return fCont[idx+1];
}

////////////////////////////////////////////////////////////////////////////////
/// Return the object before obj. Returns 0 if obj is first object.

TObject *TObjArray::Before(const TObject *obj) const
{
   if (!obj) return 0;

   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

   Int_t idx = IndexOf(obj) - fLowerBound;
   if (idx == -1 || idx == 0) return 0;

   return fCont[idx-1];
}

////////////////////////////////////////////////////////////////////////////////
/// Remove all objects from the array. Does not delete the objects
/// unless the TObjArray is the owner (set via SetOwner()).

void TObjArray::Clear(Option_t *)
{
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   if (IsOwner())
      Delete();
   else
      Init(fSize, fLowerBound);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove empty slots from array.

void TObjArray::Compress()
{
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   Int_t j = 0;

   for (Int_t i = 0; i < fSize; i++) {
      if (fCont[i]) {
         fCont[j] = fCont[i];
         j++;
      }
   }

   fLast = j - 1;

   for ( ; j < fSize; j++)
      fCont[j] = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove all objects from the array AND delete all heap based objects.

void TObjArray::Delete(Option_t * /* opt */)
{
   // In some case, for example TParallelCoord, a list (the pad's list of
   // primitives) will contain both the container and the containees
   // (the TParallelCoordVar) but if the Clear is being called from
   // the destructor of the container of this list, one of the first
   // thing done will be the remove the container (the pad) from the
   // list (of Primitives of the canvas) that was connecting it
   // (indirectly) to the list of cleanups.
   // Note: The Code in TParallelCoordVar was changed (circa June 2017),
   // to no longer have this behavior and thus rely on this code (by moving
   // from using Draw to Paint) but the structure might still exist elsewhere
   // so we keep this comment here.

   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   // Since we set fCont[i] only after the deletion is completed, we do not
   // lose the connection and thus do not need to take any special action.
   for (Int_t i = 0; i < fSize; i++) {
      if (fCont[i] && fCont[i]->IsOnHeap()) {
         TCollection::GarbageCollect(fCont[i]);
         fCont[i] = 0;
      }
   }

   Init(fSize, fLowerBound);
}

////////////////////////////////////////////////////////////////////////////////
/// Expand or shrink the array to newSize elements.

void TObjArray::Expand(Int_t newSize)
{
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   if (newSize < 0) {
      Error ("Expand", "newSize must be positive (%d)", newSize);
      return;
   }
   if (newSize == fSize)
      return;
   if (newSize < fSize) {
      // if the array is shrunk check whether there are nonempty entries
      for (Int_t j = newSize; j < fSize; j++)
         if (fCont[j]) {
            Error ("Expand", "expand would cut off nonempty entries at %d", j);
            return;
         }
   }
   fCont = (TObject**) TStorage::ReAlloc(fCont, newSize * sizeof(TObject*),
                                         fSize * sizeof(TObject*));
   fSize = newSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Find an object in this collection using its name. Requires a sequential
/// scan till the object has been found. Returns 0 if object with specified
/// name is not found.

TObject *TObjArray::FindObject(const char *name) const
{
   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

   Int_t nobjects = GetAbsLast()+1;
   for (Int_t i = 0; i < nobjects; ++i) {
      TObject *obj = fCont[i];
      if (obj && 0==strcmp(name, obj->GetName())) return obj;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Find an object in this collection using the object's IsEqual()
/// member function. Requires a sequential scan till the object has
/// been found. Returns 0 if object is not found.
/// Typically this function is overridden by a more efficient version
/// in concrete collection classes (e.g. THashTable).

TObject *TObjArray::FindObject(const TObject *iobj) const
{
   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

   Int_t nobjects = GetAbsLast()+1;
   for (Int_t i = 0; i < nobjects; ++i) {
      TObject *obj = fCont[i];
      if (obj && obj->IsEqual(iobj)) return obj;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream all objects in the array to or from the I/O buffer.

void TObjArray::Streamer(TBuffer &b)
{
   UInt_t R__s, R__c;
   Int_t nobjects;
   if (b.IsReading()) {
      Version_t v = b.ReadVersion(&R__s, &R__c);
      if (v > 2)
         TObject::Streamer(b);
      if (v > 1)
         fName.Streamer(b);

      if (GetEntriesFast() > 0) Clear();

      b >> nobjects;
      b >> fLowerBound;
      if (nobjects >= fSize) Expand(nobjects);
      fLast = -1;
      TObject *obj;
      for (Int_t i = 0; i < nobjects; i++) {
         obj = (TObject*) b.ReadObjectAny(TObject::Class());
         if (obj) {
            fCont[i] = obj;
            fLast = i;
         }
      }
      Changed();
      b.CheckByteCount(R__s, R__c,TObjArray::IsA());
   } else {
      R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

      R__c = b.WriteVersion(TObjArray::IsA(), kTRUE);
      TObject::Streamer(b);
      fName.Streamer(b);
      nobjects = GetAbsLast()+1;
      b << nobjects;
      b << fLowerBound;

      for (Int_t i = 0; i < nobjects; i++) {
         b << fCont[i];
      }
      b.SetByteCount(R__c, kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the object in the first slot.

TObject *TObjArray::First() const
{
   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

   return fCont[0];
}

////////////////////////////////////////////////////////////////////////////////
/// Return the object in the last filled slot. Returns 0 if no entries.

TObject *TObjArray::Last() const
{
   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

   if (fLast == -1)
      return 0;
   else
      return fCont[GetAbsLast()];
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of objects in array (i.e. number of non-empty slots).
/// Attention: use this method ONLY if you want to know the number of
/// non-empty slots. This function loops over the complete array and
/// is therefore very slow when applied in a loop. Most of the time you
/// better use GetEntriesFast() (only in case when there are no empty slots).

Int_t TObjArray::GetEntries() const
{
   Int_t cnt = 0;

   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

   for (Int_t i = 0; i < fSize; i++)
      if (fCont[i]) cnt++;

   return cnt;
}

////////////////////////////////////////////////////////////////////////////////
/// Return absolute index to last object in array. Returns -1 in case
/// array is empty.

Int_t TObjArray::GetAbsLast() const
{
   // For efficiency we need sometimes to update fLast so we have
   // to cast const away. Ugly, but making GetAbsLast() not const breaks
   // many other const functions.

   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

   if (fLast == -2) {
      for (Int_t i = fSize-1; i >= 0; i--)
         if (fCont[i]) {
            R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);
            ((TObjArray*)this)->fLast = i;
            return fLast;
         }
      R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);
      ((TObjArray*)this)->fLast = -1;
   }
   return fLast;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of objects in array (i.e. number of non-empty slots).
/// This is a thread-unsafe version of GetEntriesFast. Use it only if sure
/// it will not be invoked concurrently.

Int_t TObjArray::GetEntriesUnsafe() const
{
   if (R__unlikely(fLast == -2))
      return GetEntriesFast();
   else
      return fLast + 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return index of last object in array. Returns lowerBound-1 in case
/// array is empty.

Int_t TObjArray::GetLast() const
{
   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

   return fLowerBound+GetAbsLast();
}

////////////////////////////////////////////////////////////////////////////////
/// Return address of pointer obj. If obj is 0 returns address of container.

TObject **TObjArray::GetObjectRef(const TObject *obj) const
{
   if (!obj)
      return fCont;

   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

   Int_t index = IndexOf(obj);
   return &fCont[index];
}

////////////////////////////////////////////////////////////////////////////////
///  - obj != 0 Return index of object in array.
///             Returns lowerBound-1 in case array doesn't contain the obj.
///
///  - obj == 0 Return the index of the first empty slot.
///             Returns lowerBound-1 in case array doesn't contain any empty slot.

Int_t TObjArray::IndexOf(const TObject *obj) const
{
   Int_t i;

   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

   if (obj) {
      for (i = 0; i < fSize; i++)
         if (fCont[i] && fCont[i]->IsEqual(obj))
            return i+fLowerBound;
   } else {    // Look for the first empty slot
      for (i = 0; i < fSize; i++)
         if (!fCont[i])
            return i+fLowerBound;
   }

   return fLowerBound-1;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize a TObjArray.

void TObjArray::Init(Int_t s, Int_t lowerBound)
{
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   if (fCont && fSize != s) {
      TStorage::Dealloc(fCont);
      fCont = 0;
   }

   fSize = s;

   if (!fCont)
      fCont = (TObject**) TStorage::Alloc(fSize*sizeof(TObject*)); //new TObject* [fSize];
   memset(fCont, 0, fSize*sizeof(TObject*));
   fLowerBound = lowerBound;
   fLast = -1;
   Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns an array iterator.

TIterator *TObjArray::MakeIterator(Bool_t dir) const
{
   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);
   return new TObjArrayIter(this, dir);
}

////////////////////////////////////////////////////////////////////////////////
/// Generate an out-of-bounds error. Always returns false.

Bool_t TObjArray::OutOfBoundsError(const char *where, Int_t i) const
{
   Error(where, "index %d out of bounds (size: %d, this: 0x%lx)", i, fSize, (Long_t)this);
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove object from this collection and recursively remove the object
/// from all other objects (and collections).

void TObjArray::RecursiveRemove(TObject *obj)
{
   if (!obj) return;

   // We need to have the write lock even-though we are 'just'
   // reading as any insert or remove during the iteration will
   // invalidate fatally the cursor (e.g. might skip some items)
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   for (int i = 0; i < fSize; i++) {
      if (fCont[i] && fCont[i]->TestBit(kNotDeleted) && fCont[i]->IsEqual(obj)) {
         fCont[i] = 0;
         // recalculate array size
         if (i == fLast)
            do {
               fLast--;
            } while (fLast >= 0 && fCont[fLast] == 0);
         Changed();
      } else if (fCont[i] && fCont[i]->TestBit(kNotDeleted))
         fCont[i]->RecursiveRemove(obj);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Remove object at index idx.

TObject *TObjArray::RemoveAt(Int_t idx)
{
   if (!BoundsOk("RemoveAt", idx)) return 0;

   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   int i = idx-fLowerBound;

   TObject *obj = 0;
   if (fCont[i]) {
      obj = fCont[i];
      fCont[i] = 0;
      // recalculate array size
      if (i == fLast)
         do {
            fLast--;
         } while (fLast >= 0 && fCont[fLast] == 0);
      Changed();
   }
   return obj;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove object from array.

TObject *TObjArray::Remove(TObject *obj)
{
   if (!obj) return 0;

   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   Int_t idx = IndexOf(obj) - fLowerBound;

   if (idx == -1) return 0;

   TObject *ob = fCont[idx];
   fCont[idx] = 0;
   // recalculate array size
   if (idx == fLast)
      do {
         fLast--;
      } while (fLast >= 0 && fCont[fLast] == 0);
   Changed();
   return ob;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove objects from index idx1 to idx2 included.

void TObjArray::RemoveRange(Int_t idx1, Int_t idx2)
{
   if (!BoundsOk("RemoveRange", idx1)) return;
   if (!BoundsOk("RemoveRange", idx2)) return;

   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   idx1 -= fLowerBound;
   idx2 -= fLowerBound;

   Bool_t change = kFALSE;
   for (TObject **obj = fCont+idx1; obj <= fCont+idx2; obj++) {
      if (*obj) {
         *obj = 0;
         change = kTRUE;
      }
   }

   // recalculate array size
   if (change) Changed();
   if (idx1 < fLast || fLast > idx2) return;
   do { fLast--; } while (fLast >= 0 && fCont[fLast] == 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Set index of last object in array, effectively truncating the
/// array. Use carefully since whenever last position has to be
/// recalculated, e.g. after a Remove() or Sort() it will be reset
/// to the last non-empty slot. If last is -2 this will force the
/// recalculation of the last used slot.
/// If last is -1, this effectively truncate the array completely.

void TObjArray::SetLast(Int_t last)
{
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   if (last == -2 || last == -1)
      fLast = last;
   else if (BoundsOk("SetLast", last))
      fLast = last - fLowerBound;
}

////////////////////////////////////////////////////////////////////////////////
/// Randomize objects inside the array, i.e. permute randomly objects.
/// With fLast being the index of the last entry in the array, the following
/// algorithm is applied to the array:
///
///  - for each entry j between 0 and fLast, another entry k is chosen
///    randomly between 0 and fLast.
///  - the objects at j and k are swapped.
///  - this process is repeated ntimes (ntimes = 1 by default).

void TObjArray::Randomize(Int_t ntimes)
{
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   for (Int_t i = 0; i < ntimes; i++) {
      for (Int_t j = 0; j < fLast; j++) {
#ifdef R__WIN32
         Int_t k = (Int_t)(0.5+fLast*Double_t(rand())/Double_t((RAND_MAX+1.0)));
#else
         Int_t k = (Int_t)(0.5+fLast*Double_t(random())/Double_t((RAND_MAX+1.0)));
#endif
         if (k == j) continue;
         TObject *obj = fCont[j];
         fCont[j] = fCont[k];
         fCont[k] = obj;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// If objects in array are sortable (i.e. IsSortable() returns true
/// for all objects) then sort array.

void TObjArray::Sort(Int_t upto)
{
   R__COLLECTION_WRITE_LOCKGUARD(ROOT::gCoreMutex);

   if (GetAbsLast() == -1 || fSorted) return;
   for (Int_t i = 0; i < fSize; i++)
      if (fCont[i]) {
         if (!fCont[i]->IsSortable()) {
            Error("Sort", "objects in array are not sortable");
            return;
         }
      }

   QSort(fCont, 0, TMath::Min(fSize, upto-fLowerBound));

   fLast   = -2;
   fSorted = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Find object using a binary search. Array must first have been sorted.
/// Search can be limited by setting upto to desired index.

Int_t TObjArray::BinarySearch(TObject *op, Int_t upto)
{
   R__COLLECTION_READ_LOCKGUARD(ROOT::gCoreMutex);

   Int_t   base, position, last, result = 0;
   TObject *op2;

   if (!op) return -1;

   if (!fSorted) {
      Error("BinarySearch", "array must first be sorted");
      return -1;
   }

   base = 0;
   last = TMath::Min(fSize, upto-fLowerBound) - 1;

   while (last >= base) {
      position = (base+last) / 2;
      op2 = fCont[position];
      if (op2 && (result = op->Compare(op2)) == 0)
         return position + fLowerBound;
      if (!op2 || result < 0)
         last = position-1;
      else
         base = position+1;
   }
   return -1;
}

/** \class TObjArrayIter
Iterator of object array.
*/

ClassImp(TObjArrayIter);

////////////////////////////////////////////////////////////////////////////////
/// Create array iterator. By default the iteration direction
/// is kIterForward. To go backward use kIterBackward.

TObjArrayIter::TObjArrayIter(const TObjArray *arr, Bool_t dir)
{
   fArray     = arr;
   fDirection = dir;
   Reset();
}

////////////////////////////////////////////////////////////////////////////////
/// Copy ctor.

TObjArrayIter::TObjArrayIter(const TObjArrayIter &iter) : TIterator(iter)
{
   fArray     = iter.fArray;
   fDirection = iter.fDirection;
   fCursor    = iter.fCursor;
   fCurCursor = iter.fCurCursor;
}

////////////////////////////////////////////////////////////////////////////////
/// Overridden assignment operator.

TIterator &TObjArrayIter::operator=(const TIterator &rhs)
{
   if (this != &rhs && rhs.IsA() == TObjArrayIter::Class()) {
      const TObjArrayIter &rhs1 = (const TObjArrayIter &)rhs;
      fArray     = rhs1.fArray;
      fDirection = rhs1.fDirection;
      fCursor    = rhs1.fCursor;
      fCurCursor = rhs1.fCurCursor;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Overloaded assignment operator.

TObjArrayIter &TObjArrayIter::operator=(const TObjArrayIter &rhs)
{
   if (this != &rhs) {
      fArray     = rhs.fArray;
      fDirection = rhs.fDirection;
      fCursor    = rhs.fCursor;
      fCurCursor = rhs.fCurCursor;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Return next object in array. Returns 0 when no more objects in array.

TObject *TObjArrayIter::Next()
{
   if (fDirection == kIterForward) {
      for ( ; fCursor < fArray->Capacity() && fArray->fCont[fCursor] == 0;
              fCursor++) { }

      fCurCursor = fCursor;
      if (fCursor < fArray->Capacity()) {
         return fArray->fCont[fCursor++];
      }
   } else {
      for ( ; fCursor >= 0 && fArray->fCont[fCursor] == 0;
              fCursor--) { }

      fCurCursor = fCursor;
      if (fCursor >= 0) {
         return fArray->fCont[fCursor--];
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Reset array iterator.

void TObjArrayIter::Reset()
{
   if (fDirection == kIterForward)
      fCursor = 0;
   else
      fCursor = fArray->Capacity() - 1;

   fCurCursor = fCursor;
}

////////////////////////////////////////////////////////////////////////////////
/// This operator compares two TIterator objects.

Bool_t TObjArrayIter::operator!=(const TIterator &aIter) const
{
   if (aIter.IsA() == TObjArrayIter::Class()) {
      const TObjArrayIter &iter(dynamic_cast<const TObjArrayIter &>(aIter));
      return (fCurCursor != iter.fCurCursor);
   }
   return false; // for base class we don't implement a comparison
}

////////////////////////////////////////////////////////////////////////////////
/// This operator compares two TObjArrayIter objects.

Bool_t TObjArrayIter::operator!=(const TObjArrayIter &aIter) const
{
   return (fCurCursor != aIter.fCurCursor);
}

////////////////////////////////////////////////////////////////////////////////
/// Return current object or nullptr.

TObject *TObjArrayIter::operator*() const
{
   return (((fCurCursor >= 0) && (fCurCursor < fArray->Capacity())) ?
           fArray->fCont[fCurCursor] : nullptr);
}
