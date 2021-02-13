// @(#)root/cont:$Id$
// Author: Rene Brun  02/10/2001

/*************************************************************************
 * Copyright (C) 1995-2001, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TRefArray
\ingroup Containers
An array of references to TObjects. The array expands automatically
when  objects are added (shrinking can be done by hand using Expand() )

The TRefArray can be filled with:
~~~ {.cpp}
    array.Add(obj)
    array.AddAt(obj,i)
    but not array[i] = obj  !!!
~~~
The array elements can be retrieved with:
~~~ {.cpp}
    TObject *obj = array.At(i);
~~~
By default the TRefArray 'points' to the current process and can only
receive object that have been created in this process.
To point the TRefArray to a different process do:
~~~ {.cpp}
    TRefArray array( processId );
~~~
For example, if 'obj' is an instance that was created in the different
process and you do:
~~~ {.cpp}
    TRefArray array( TProcessID::GetProcessWithUID( obj ) );
~~~
Then
~~~ {.cpp}
    array.Add(obj);
~~~
is correct (obj comes from the process the array is pointed to
while
~~~ {.cpp}
    TObject *nobj = new TObject;
    array.Add(nobj);
~~~
is incorrect since 'nobj' was created in a different process than the
one the array is pointed to. In this case you will see error message:
~~~ {.cpp}
    Error in <TRefArray::AddAtAndExpand>: The object at 0x... is not
    registered in the process the TRefArray point to
    (pid = ProcessID../....)
~~~
When a TRefArray is Streamed, only the pointer unique id is written,
not the referenced object. TRefArray may be assigned to different
branches of one Tree or several Trees.
The branch containing the TRefArray can be read before or after the
array (eg TClonesArray, STL vector,..) of the referenced objects.

See an example in $ROOTSYS/test/Event.h

### RESTRICTIONS when using TRefArray

  - Elements in a TRefArray cannot point to a TFile or TDirectory.
  - All elements of a TRefArray must be set in the same process,
    In particular, one cannot modify some elements of the array in
    a different process.

Use an array of TRef when one of the above restrictions is met.

The number of TRef handled by a single process id is limited to
16777215 (see TRef for more detail).   When the TProcessID is full
(has seen 16777215 objects), we switch to new one TProcessID
maximum 65535 including the TProcessIDs read from file).
However TRefArray can not switch to new TProcessID if they already
contain objects.

When the TProcessID has been switched due to overflow and an new
object is added to an existing, empty TRefArray, you will see:

~~~ {.cpp}
Warning in <TRefArray::AddAtAndExpand>: The ProcessID for the 0x5f83819e8 has been switched to ProcessID4/6c89f37e-8259-11e2-9717-166ee183beef:4
~~~
If the TRefArray was not empty, you will see:

~~~ {.cpp}
Error in <TRefArray::AddAtAndExpand>: The object at %p can not be registered in the process the TRefArray points to (pid = ProcessID4/6c89f37e-8259-11e2-9717-166ee183beef) because the ProcessID has too many objects and the TRefArray already contains other objects.
~~~
When running out of TProcessIds, you will see:

~~~ {.cpp}
Warning in <TProcessID::AddProcessID>: Maximum number of TProcessID (65535) is almost reached (one left).  TRef will stop being functional when the limit is reached.
Fatal in <TProcessID::AddProcessID>: Maximum number of TProcessID (65535) has been reached.  TRef are not longer functional.
~~~
*/

#include "TRefArray.h"
#include "TRefTable.h"
#include "TBuffer.h"
#include "TError.h"
#include "TBits.h"

ClassImp(TRefArray);

////////////////////////////////////////////////////////////////////////////////
/// default constructor

TRefArray::TRefArray(TProcessID *pid)
{
   fPID  = pid ? pid : TProcessID::GetSessionProcessID();
   fUIDs = 0;
   fSize = 0;
   fLast = -1;
   fLowerBound = 0;
   Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Create an object array. Using s one can set the array size
/// and lowerBound can be used to set the array lowerbound
/// index (default is 0).

TRefArray::TRefArray(Int_t s, TProcessID *pid)
{
   if (s < 0) {
      Warning("TRefArray", "size (%d) < 0", s);
      s = TCollection::kInitCapacity;
   }

   fPID  = pid ? pid : TProcessID::GetSessionProcessID();
   fUIDs = 0;
   Init(s, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Create an object array. Using s one can set the array size
/// and lowerBound can be used to set the array lowerbound
/// index (default is 0).

TRefArray::TRefArray(Int_t s, Int_t lowerBound, TProcessID *pid)
{
   if (s < 0) {
      Warning("TRefArray", "size (%d) < 0", s);
      s = TCollection::kInitCapacity;
   }

   fPID  = pid ? pid : TProcessID::GetSessionProcessID();
   fUIDs = 0;
   Init(s, lowerBound);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a copy of TRefArray a.

TRefArray::TRefArray(const TRefArray &a) : TSeqCollection()
{
   fPID  = a.fPID;
   fUIDs = 0;
   Init(a.fSize, a.fLowerBound);

   for (Int_t i = 0; i < fSize; i++)
      fUIDs[i] = a.fUIDs[i];

   fLast = a.fLast;
   fName = a.fName;
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TRefArray& TRefArray::operator=(const TRefArray &a)
{
   if (this != &a) {
      // Copy this by hand because the assignment operator
      // of TCollection is private
      fName   = a.fName;
      fSize   = a.fSize;
      fSorted = a.fSorted;

      fPID = a.fPID;
      Init(a.fSize, a.fLowerBound);

      for (Int_t i = 0; i < fSize; i++)
         fUIDs[i] = a.fUIDs[i];

      fLast = a.fLast;
      fName = a.fName;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Usual destructor (The object pointed to by the array are never deleted).

TRefArray::~TRefArray()
{
   if (fUIDs) delete [] fUIDs;
   fPID  = 0;
   fUIDs = 0;
   fSize = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Private/static function, check for validity of pid.

Bool_t TRefArray::GetObjectUID(Int_t &uid, TObject *obj, const char *methodname)
{
   // Check if the object can belong here.
   Bool_t valid = kTRUE;
   if (obj->TestBit(kHasUUID)) {
      valid = kFALSE;
   } else if (obj->TestBit(kIsReferenced)) {
      valid = (fPID == TProcessID::GetProcessWithUID(obj));
      if (valid) {
         uid = obj->GetUniqueID();
      } else {
         if (GetAbsLast() < 0) {
            // The container is empty, we can switch the ProcessID.
            fPID = TProcessID::GetProcessWithUID(obj);
            valid = kTRUE;
            if (gDebug > 3)
               Info(TString::Format("TRefArray::%s",methodname),"The ProcessID for the %p has been switched to %s/%s:%d.",
                    this,fPID->GetName(),fPID->GetTitle(),fPID->GetUniqueID());
        }
      }
   } else {
      // If we could, we would just add the object to the
      // TRefArray's ProcessID.  For now, just check the
      // ProcessID it would be added to, i.e the current one,
      // is not full.

      if (!(TProcessID::GetObjectCount() >= 16777215)) {
         valid = (fPID == TProcessID::GetSessionProcessID());
         if (valid) {
            uid = TProcessID::AssignID(obj);
         }
      } else {
         // The AssignID will create a new TProcessID.
         if (GetAbsLast() < 0) {
            // If we are empty, we can handle it.
            uid = TProcessID::AssignID(obj);
            fPID = TProcessID::GetProcessWithUID(obj);
            Warning(TString::Format("TRefArray::%s",methodname),"The ProcessID for the %p has been switched to %s/%s:%d. There are too many referenced objects.",
                    this,fPID->GetName(),fPID->GetTitle(),fPID->GetUniqueID());
            return kTRUE;
        } else {
            Error(TString::Format("TRefArray::%s",methodname),"The object at %p can not be registered in the process the TRefArray points to (pid = %s/%s) because the ProcessID has too many objects and the TRefArray already contains other objects.",obj,fPID->GetName(),fPID->GetTitle());
            return kFALSE;
         }
      }
   }

   if (!valid) {
      ::Error(TString::Format("TRefArray::%s",methodname),
              "The object at %p is not registered in the process the TRefArray points to (pid = %s/%s)",obj,fPID->GetName(),fPID->GetTitle());
   }
   return valid;
}

////////////////////////////////////////////////////////////////////////////////
/// Add object in the first slot of the array. This will overwrite the
/// first element that might have been there. To have insertion semantics
/// use either a TList or a TOrdCollection.

void TRefArray::AddFirst(TObject *obj)
{
   if (!obj) return;

   // Check if the object can belong here
   Int_t uid;
   if (GetObjectUID(uid, obj, "AddFirst")) {
      fUIDs[0] = uid;
      Changed();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add object in the next empty slot in the array. Expand the array
/// if necessary.

void TRefArray::AddLast(TObject *obj)
{
   AddAtAndExpand(obj, GetAbsLast()+1+fLowerBound);
}

////////////////////////////////////////////////////////////////////////////////
/// Add object in the slot before object before. If before=0 add object
/// in the first slot. Note that this will overwrite any object that
/// might have already been in this slot. For insertion semantics use
/// either a TList or a TOrdCollection.

void TRefArray::AddBefore(const TObject *before, TObject *obj)
{
   if (!before)
      AddFirst(obj);
   else {
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

void TRefArray::AddAfter(const TObject *after, TObject *obj)
{
   if (!after)
      AddLast(obj);
   else {
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

void TRefArray::AddAtAndExpand(TObject *obj, Int_t idx)
{
   if (!obj) return;
   if (idx < fLowerBound) {
      Error("AddAt", "out of bounds at %d in %lx", idx, (Long_t)this);
      return;
   }
   if (idx-fLowerBound >= fSize)
      Expand(TMath::Max(idx-fLowerBound+1, GrowBy(fSize)));

   // Check if the object can belong here
   Int_t uid;
   if (GetObjectUID(uid, obj, "AddAtAndExpand")) {
      fUIDs[idx-fLowerBound] = uid;
      fLast = TMath::Max(idx-fLowerBound, GetAbsLast());
      Changed();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add object at position ids. Give an error when idx is out of bounds
/// (i.e. the array is not expanded).

void TRefArray::AddAt(TObject *obj, Int_t idx)
{
   if (!obj) return;
   if (!BoundsOk("AddAt", idx)) return;

   // Check if the object can belong here
   Int_t uid;
   if (GetObjectUID(uid, obj, "AddAt")) {
      fUIDs[idx-fLowerBound] = uid;;
      fLast = TMath::Max(idx-fLowerBound, GetAbsLast());
      Changed();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the position of the new object.
/// Find the first empty cell or AddLast if there is no empty cell

Int_t  TRefArray::AddAtFree(TObject *obj)
{
   if (!obj) return 0;
   if (Last()) {    // <---------- This is to take in account "empty" TRefArray's
      Int_t i;
      for (i = 0; i < fSize; i++)
         if (!fUIDs[i]) {         // Add object at position i
            // Check if the object can belong here
            Int_t uid;
            if (GetObjectUID(uid, obj, "AddAtFree")) {
               fUIDs[i] = uid;
               fLast = TMath::Max(i, GetAbsLast());
               Changed();
               return i+fLowerBound;
            }
         }
   }
   AddLast(obj);
   return GetLast();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the object after obj. Returns 0 if obj is last object.

TObject *TRefArray::After(const TObject *obj) const
{
   if (!obj || !fPID) return 0;

   Int_t idx = IndexOf(obj) - fLowerBound;
   if (idx == -1 || idx == fSize-1) return 0;

   return fPID->GetObjectWithID(fUIDs[idx+1]);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the object before obj. Returns 0 if obj is first object.

TObject *TRefArray::Before(const TObject *obj) const
{
   if (!obj || !fPID) return 0;

   Int_t idx = IndexOf(obj) - fLowerBound;
   if (idx == -1 || idx == 0) return 0;

   return fPID->GetObjectWithID(fUIDs[idx-1]);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove all objects from the array.

void TRefArray::Clear(Option_t *)
{
   fLast = - 1;

   for (Int_t j=0 ; j < fSize; j++) fUIDs[j] = 0;

   Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Remove empty slots from array.

void TRefArray::Compress()
{
   Int_t j = 0;

   for (Int_t i = 0; i < fSize; i++) {
      if (fUIDs[i]) {
         fUIDs[j] = fUIDs[i];
         j++;
      }
   }

   fLast = j - 1;

   for ( ; j < fSize; j++) fUIDs[j] = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove all objects from the array and free the internal memory.

void TRefArray::Delete(Option_t *)
{
   fLast = -1;

   fSize = 0;
   if (fUIDs) {
      delete [] fUIDs;
      fUIDs = 0;
   }

   Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Expand or shrink the array to newSize elements.

void TRefArray::Expand(Int_t newSize)
{
   if (newSize < 0) {
      Error ("Expand", "newSize must be positive (%d)", newSize);
      return;
   }
   if (newSize == fSize) return;
   UInt_t *temp = fUIDs;
   if (newSize != 0) {
      fUIDs = new UInt_t[newSize];
      if (newSize < fSize) memcpy(fUIDs,temp, newSize*sizeof(UInt_t));
      else {
         memcpy(fUIDs,temp,fSize*sizeof(UInt_t));
         memset(&fUIDs[fSize],0,(newSize-fSize)*sizeof(UInt_t));
      }
   } else {
      fUIDs = 0;
   }
   if (temp) delete [] temp;
   fSize = newSize;
}

////////////////////////////////////////////////////////////////////////////////
///the reference may be in the TRefTable

TObject *TRefArray::GetFromTable(Int_t idx) const
{
   TRefTable *table = TRefTable::GetRefTable();
   if (table) {
      table->SetUID(fUIDs[idx], fPID);
      table->Notify();
      return fPID->GetObjectWithID(fUIDs[idx]);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream all objects in the array to or from the I/O buffer.

void TRefArray::Streamer(TBuffer &R__b)
{
   UInt_t R__s, R__c;
   Int_t nobjects;
   UShort_t pidf;
   if (R__b.IsReading()) {
      R__b.ReadVersion(&R__s, &R__c);
      TObject::Streamer(R__b);
      fName.Streamer(R__b);
      R__b >> nobjects;
      R__b >> fLowerBound;
      if (nobjects >= fSize) Expand(nobjects);
      fLast = -1;
      R__b >> pidf;
      pidf += R__b.GetPidOffset();
      fPID = R__b.ReadProcessID(pidf);
      if (gDebug > 1) printf("Reading TRefArray, pidf=%d, fPID=%lx, nobjects=%d\n",pidf,(Long_t)fPID,nobjects);
      for (Int_t i = 0; i < nobjects; i++) {
         R__b >> fUIDs[i];
         if (fUIDs[i] != 0) fLast = i;
         if (gDebug > 1) {
            printf(" %d",fUIDs[i]);
            if ((i > 0 && i%10 == 0) || (i == nobjects-1)) printf("\n");
         }
      }
      memset(&fUIDs[fLast+1], 0, (fSize - fLast - 1)*sizeof(fUIDs[0]));
      Changed();
      R__b.CheckByteCount(R__s, R__c,TRefArray::IsA());
   } else {
      R__c = R__b.WriteVersion(TRefArray::IsA(), kTRUE);
      TObject::Streamer(R__b);
      fName.Streamer(R__b);
      nobjects = GetAbsLast()+1;
      R__b << nobjects;
      R__b << fLowerBound;
      pidf = R__b.WriteProcessID(fPID);
      R__b << pidf;
      if (gDebug > 1) printf("Writing TRefArray, pidf=%d, fPID=%lx, nobjects=%d\n",pidf,(Long_t)fPID,nobjects);

      for (Int_t i = 0; i < nobjects; i++) {
         R__b << fUIDs[i];
         if (gDebug > 1) {
            printf(" %d",fUIDs[i]);
            if ((i > 0 && i%10 == 0) || (i == nobjects-1)) printf("\n");
         }
      }
      R__b.SetByteCount(R__c, kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the object in the first slot.

TObject *TRefArray::First() const
{
   return fPID->GetObjectWithID(fUIDs[0]);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the object in the last filled slot. Returns 0 if no entries.

TObject *TRefArray::Last() const
{
   if (fLast == -1)
      return 0;
   else
      return fPID->GetObjectWithID(fUIDs[GetAbsLast()]);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of objects in array (i.e. number of non-empty slots).
/// Attention: use this method ONLY if you want to know the number of
/// non-empty slots. This function loops over the complete array and
/// is therefore very slow when applied in a loop. Most of the time you
/// better use GetLast()+1.

Int_t TRefArray::GetEntries() const
{
   Int_t cnt = 0;

   for (Int_t i = 0; i < fSize; i++)
      if (fUIDs[i]) cnt++;

   return cnt;
}

////////////////////////////////////////////////////////////////////////////////
/// Return absolute index to last object in array. Returns -1 in case
/// array is empty.

Int_t TRefArray::GetAbsLast() const
{
   // For efficiency we need sometimes to update fLast so we have
   // to cast const away. Ugly, but making GetAbsLast() not const breaks
   // many other const functions.
   if (fLast == -2) {
      for (Int_t i = fSize-1; i >= 0; i--)
         if (fUIDs[i]) {
            ((TRefArray*)this)->fLast = i;
            return fLast;
         }
      ((TRefArray*)this)->fLast = -1;
   }
   return fLast;
}

////////////////////////////////////////////////////////////////////////////////
/// Return index of last object in array. Returns lowerBound-1 in case
/// array is empty.

Int_t TRefArray::GetLast() const
{
   return fLowerBound+GetAbsLast();
}

////////////////////////////////////////////////////////////////////////////////
/// Return address of pointer obj.

TObject **TRefArray::GetObjectRef(const TObject *) const
{
   //Int_t index = IndexOf(obj);
   //return &fCont[index];
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return UID of element at.

UInt_t TRefArray::GetUID(Int_t at) const
{
   int j = at-fLowerBound;
   if (j >= 0 && j < fSize) {
      if (!fPID) return 0;
      return fUIDs[j];
   }
   BoundsOk("At", at);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
///  - obj != 0 Return index of object in array.
///             Returns lowerBound-1 in case array doesn't contain the obj.
///
///  - obj == 0 Return the index of the first empty slot.
///             Returns lowerBound-1 in case array doesn't contain any empty slot.

Int_t TRefArray::IndexOf(const TObject *obj) const
{
   Int_t i;
   if (obj) {
      if (!TProcessID::IsValid(fPID)) {
         return fLowerBound-1;
      }
      for (i = 0; i < fSize; i++)
         if (fUIDs[i] && fPID->GetObjectWithID(fUIDs[i]) == obj)
            return i+fLowerBound;
   } else {    // Look for the first empty slot
      for (i = 0; i < fSize; i++)
         if (!fUIDs[i])
            return i+fLowerBound;
   }

   return fLowerBound-1;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize a TRefArray.

void TRefArray::Init(Int_t s, Int_t lowerBound)
{
   if (fUIDs && fSize != s) {
      delete [] fUIDs;
      fUIDs = 0;
   }

   fSize = s;

   if (fSize) {
      fUIDs = new UInt_t[fSize];
      for (Int_t i=0;i<s;i++) fUIDs[i] = 0;
   } else {
      fUIDs = 0;
   }
   fLowerBound = lowerBound;
   fLast = -1;
   Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns an array iterator.

TIterator *TRefArray::MakeIterator(Bool_t dir) const
{
   return new TRefArrayIter(this, dir);
}

////////////////////////////////////////////////////////////////////////////////
/// Generate an out-of-bounds error. Always returns false.

Bool_t TRefArray::OutOfBoundsError(const char *where, Int_t i) const
{
   Error(where, "index %d out of bounds (size: %d, this: 0x%lx)", i, fSize, (Long_t)this);
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove object at index idx.

TObject *TRefArray::RemoveAt(Int_t idx)
{
   if (!BoundsOk("RemoveAt", idx)) return 0;

   int i = idx-fLowerBound;

   TObject *obj = 0;
   if (fUIDs[i]) {
      if (!TProcessID::IsValid(fPID)) return 0;
      obj = fPID->GetObjectWithID(fUIDs[i]);
      fUIDs[i] = 0;
      // recalculate array size
      if (i == fLast)
         do {
            fLast--;
         } while (fLast >= 0 && fUIDs[fLast] == 0);
      Changed();
   }

   return obj;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove object from array.

TObject *TRefArray::Remove(TObject *obj)
{
   if (!obj) return 0;

   Int_t idx = IndexOf(obj) - fLowerBound;

   if (idx == -1) return 0;

   TObject *ob = fPID->GetObjectWithID(fUIDs[idx]);
   fUIDs[idx] = 0;
   // recalculate array size
   if (idx == fLast)
      do {
         fLast--;
      } while (fLast >= 0 && fUIDs[fLast] == 0);
   Changed();
   return ob;
}

////////////////////////////////////////////////////////////////////////////////
/// Set index of last object in array, effectively truncating the
/// array. Use carefully since whenever last position has to be
/// recalculated, e.g. after a Remove() or Sort() it will be reset
/// to the last non-empty slot. If last is -2 this will force the
/// recalculation of the last used slot.

void TRefArray::SetLast(Int_t last)
{
   if (last == -2)
      fLast = -2;
   else if (BoundsOk("SetLast", last))
      fLast = last - fLowerBound;
}

////////////////////////////////////////////////////////////////////////////////
/// If objects in array are sortable (i.e. IsSortable() returns true
/// for all objects) then sort array.

void TRefArray::Sort(Int_t)
{
   Error("Sort","Function not yet implemented");
/*
   if (GetAbsLast() == -1 || fSorted) return;
   for (Int_t i = 0; i < fSize; i++)
      if (fUIDs[i]) {
         if (!fUIDs[i]->IsSortable()) {
            Error("Sort", "objects in array are not sortable");
            return;
         }
      }

   QSort(fUIDs, 0, TMath::Min(fSize, upto-fLowerBound));

   fLast   = -2;
   fSorted = kTRUE;
*/
}

////////////////////////////////////////////////////////////////////////////////
/// Find object using a binary search. Array must first have been sorted.
/// Search can be limited by setting upto to desired index.

Int_t TRefArray::BinarySearch(TObject *, Int_t)
{
   Error("BinarySearch","Function not yet implemented");
/*
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
      //position = (base+last) / 2;
      //op2 = fCont[position];
      //if (op2 && (result = op->Compare(op2)) == 0)
      //   return position + fLowerBound;
      //if (!op2 || result < 0)
      //   last = position-1;
      //else
      //   base = position+1;
   }
*/
   return -1;
}

/** \class TRefArrayIter
Iterator of object array.
*/

ClassImp(TRefArrayIter);

////////////////////////////////////////////////////////////////////////////////
/// Create array iterator. By default the iteration direction
/// is kIterForward. To go backward use kIterBackward.

TRefArrayIter::TRefArrayIter(const TRefArray *arr, Bool_t dir)
{
   fArray     = arr;
   fDirection = dir;
   Reset();
}

////////////////////////////////////////////////////////////////////////////////
/// Copy ctor.

TRefArrayIter::TRefArrayIter(const TRefArrayIter &iter) : TIterator(iter)
{
   fArray     = iter.fArray;
   fDirection = iter.fDirection;
   fCursor    = iter.fCursor;
   fCurCursor = iter.fCurCursor;
}

////////////////////////////////////////////////////////////////////////////////
/// Overridden assignment operator.

TIterator &TRefArrayIter::operator=(const TIterator &rhs)
{
   if (this != &rhs && rhs.IsA() == TRefArrayIter::Class()) {
      const TRefArrayIter &rhs1 = (const TRefArrayIter &)rhs;
      fArray     = rhs1.fArray;
      fDirection = rhs1.fDirection;
      fCursor    = rhs1.fCursor;
      fCurCursor = rhs1.fCurCursor;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Overloaded assignment operator.

TRefArrayIter &TRefArrayIter::operator=(const TRefArrayIter &rhs)
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

TObject *TRefArrayIter::Next()
{
   if (fDirection == kIterForward) {
      for ( ; fCursor < fArray->Capacity() && fArray->At(fCursor+fArray->LowerBound()) == 0;
              fCursor++) { }

      fCurCursor = fCursor;
      if (fCursor < fArray->Capacity()) {
         fCursor++;
         return fArray->At(fCurCursor+fArray->LowerBound());
      }
   } else {
      for ( ; fCursor >= 0 && fArray->At(fCursor) == 0;
              fCursor--) { }

      fCurCursor = fCursor;
      if (fCursor >= 0) {
         fCursor--;
         return fArray->At(fCurCursor+fArray->LowerBound());
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Reset array iterator.

void TRefArrayIter::Reset()
{
   if (fDirection == kIterForward)
      fCursor = 0;
   else
      fCursor = fArray->Capacity() - 1;

   fCurCursor = fCursor;
}

////////////////////////////////////////////////////////////////////////////////
/// This operator compares two TIterator objects.

Bool_t TRefArrayIter::operator!=(const TIterator &aIter) const
{
   if (aIter.IsA() == TRefArrayIter::Class()) {
      const TRefArrayIter &iter(dynamic_cast<const TRefArrayIter &>(aIter));
      return (fCurCursor != iter.fCurCursor);
   }
   return false; // for base class we don't implement a comparison
}

////////////////////////////////////////////////////////////////////////////////
/// This operator compares two TRefArrayIter objects.

Bool_t TRefArrayIter::operator!=(const TRefArrayIter &aIter) const
{
   return (fCurCursor != aIter.fCurCursor);
}

////////////////////////////////////////////////////////////////////////////////
/// Return current object or nullptr.

TObject *TRefArrayIter::operator*() const
{
   return (((fCurCursor >= 0) && (fCurCursor < fArray->Capacity())) ?
           fArray->At(fCurCursor) : nullptr);
}
