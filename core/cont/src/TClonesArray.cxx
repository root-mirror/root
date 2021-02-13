// @(#)root/cont:$Id$
// Author: Rene Brun   11/02/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TClonesArray
\ingroup Containers
An array of clone (identical) objects. Memory for the objects
stored in the array is allocated only once in the lifetime of the
clones array. All objects must be of the same class. For the rest
this class has the same properties as TObjArray.

To reduce the very large number of new and delete calls in large
loops like this (O(100000) x O(10000) times new/delete):
~~~ {.cpp}
  TObjArray a(10000);
  while (TEvent *ev = (TEvent *)next()) {      // O(100000) events
     for (int i = 0; i < ev->Ntracks; i++) {   // O(10000) tracks
        a[i] = new TTrack(x,y,z,...);
        ...
        ...
     }
     ...
     a.Delete();
  }
~~~
One better uses a TClonesArray which reduces the number of
new/delete calls to only O(10000):
~~~ {.cpp}
  TClonesArray a("TTrack", 10000);
  while (TEvent *ev = (TEvent *)next()) {      // O(100000) events
     for (int i = 0; i < ev->Ntracks; i++) {   // O(10000) tracks
        new(a[i]) TTrack(x,y,z,...);
        ...
        ...
     }
     ...
     a.Delete(); // or a.Clear() or a.Clear("C")
  }
~~~
To reduce the number of call to the constructor (especially useful
if the user class requires memory allocation), the object can be
added (and constructed when needed) using ConstructedAt which only
calls the constructor once per slot.
~~~ {.cpp}
  TClonesArray a("TTrack", 10000);
  while (TEvent *ev = (TEvent *)next()) {      // O(100000) events
     for (int i = 0; i < ev->Ntracks; i++) {   // O(10000) tracks
        TTrack *track = (TTrack*)a.ConstructedAt(i);
        track->Set(x,y,z,....);
        ...
        ...
     }
     ...
     a.Clear(); // or a.Clear("C");
  }
~~~
Note: the only supported way to add objects to a TClonesArray is
via the new with placement method or the ConstructedAt method.
The other Add() methods ofTObjArray and its base classes are not
allowed.

Considering that a new/delete costs about 70 mus on a 300 MHz HP,
O(10^9) new/deletes will save about 19 hours.

### NOTE 1

C/C++ offers the possibility of allocating and deleting memory.
Forgetting to delete allocated memory is a programming error called a
"memory leak", i.e. the memory of your process grows and eventually
your program crashes. Even if you *always* delete the allocated
memory, the recovered space may not be efficiently reused. The
process knows that there are portions of free memory, but when you
allocate it again, a fresh piece of memory is grabbed. Your program
is free from semantic errors, but the total memory of your process
still grows, because your program's memory is full of "holes" which
reduce the efficiency of memory access; this is called "memory
fragmentation". Moreover new / delete are expensive operations in
terms of CPU time.

Without entering into technical details, TClonesArray allows you to
"reuse" the same portion of memory for new/delete avoiding memory
fragmentation and memory growth and improving the performance by
orders of magnitude. Every time the memory of the TClonesArray has
to be reused, the Clear() method is used. To provide its benefits,
each TClonesArray must be allocated *once* per process and disposed
of (deleted) *only when not needed any more*.

So a job should see *only one* deletion for each TClonesArray,
which should be Clear()ed during the job several times. Deleting a
TClonesArray is a double waste. Not only you do not avoid memory
fragmentation, but you worsen it because the TClonesArray itself
is a rather heavy structure, and there is quite some code in the
destructor, so you have more memory fragmentation and slower code.

### NOTE 2

When investigating misuse of TClonesArray, please make sure of the following:

   - Use Clear() or Clear("C") instead of Delete(). This will improve
     program execution time.
   - TClonesArray object classes containing pointers allocate memory.
     To avoid causing memory leaks, special Clear("C") must be used
     for clearing TClonesArray. When option "C" is specified, ROOT
     automatically executes the Clear() method (by default it is
     empty contained in TObject). This method must be overridden in
     the relevant TClonesArray object class, implementing the reset
     procedure for pointer objects.
   - If the objects are added using the placement new then the Clear must
     deallocate the memory.
   - If the objects are added using TClonesArray::ConstructedAt then the
     heap-based memory can stay allocated and reused as the constructor is
     not called for already constructed/added object.
   - To reduce memory fragmentation, please make sure that the
     TClonesArrays are not destroyed and created on every event. They
     must only be constructed/destructed at the beginning/end of the
     run.
*/

#include "TClonesArray.h"

#include "TError.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TClass.h"
#include "TObject.h"
#include "TObjectTable.h"
#include "snprintf.h"

#include <cstdlib>

ClassImp(TClonesArray);

// To allow backward compatibility of TClonesArray of v5 TF1 objects
// that were stored member-wise.
using Updater_t = void (*)(Int_t nobjects, TObject **from, TObject **to);
Updater_t gClonesArrayTF1Updater = nullptr;
Updater_t gClonesArrayTFormulaUpdater = nullptr;

bool R__SetClonesArrayTF1Updater(Updater_t func) {
   gClonesArrayTF1Updater = func;
   return true;
}

bool R__SetClonesArrayTFormulaUpdater(Updater_t func) {
   gClonesArrayTFormulaUpdater = func;
   return true;
}

/// Internal Utility routine to correctly release the memory for an object
static inline void R__ReleaseMemory(TClass *cl, TObject *obj)
{
   if (obj && obj->TestBit(TObject::kNotDeleted)) {
      // -- The TObject destructor has not been called.
      cl->Destructor(obj);
   } else {
      // -- The TObject destructor was called, just free memory.
      //
      // remove any possible entries from the ObjectTable
      if (TObject::GetObjectStat() && gObjectTable) {
         gObjectTable->RemoveQuietly(obj);
      }
      ::operator delete(obj);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Default Constructor.

TClonesArray::TClonesArray() : TObjArray()
{
   fClass      = 0;
   fKeep       = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create an array of clone objects of classname. The class must inherit from
/// TObject.
/// The second argument s indicates an approximate number of objects
/// that will be entered in the array. If more than s objects are entered,
/// the array will be automatically expanded.
///
/// The third argument is not used anymore and only there for backward
/// compatibility reasons.

TClonesArray::TClonesArray(const char *classname, Int_t s, Bool_t) : TObjArray(s)
{
   fKeep = 0;
   SetClass(classname,s);
}

////////////////////////////////////////////////////////////////////////////////
/// Create an array of clone objects of class cl. The class must inherit from
/// TObject.
/// The second argument, s, indicates an approximate number of objects
/// that will be entered in the array. If more than s objects are entered,
/// the array will be automatically expanded.
///
/// The third argument is not used anymore and only there for backward
/// compatibility reasons.

TClonesArray::TClonesArray(const TClass *cl, Int_t s, Bool_t) : TObjArray(s)
{
   fKeep = 0;
   SetClass(cl,s);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy ctor.

TClonesArray::TClonesArray(const TClonesArray& tc): TObjArray(tc)
{
   fKeep = new TObjArray(tc.fSize);
   fClass = tc.fClass;

   BypassStreamer(kTRUE);

   for (Int_t i = 0; i < fSize; i++) {
      if (tc.fCont[i]) fCont[i] = tc.fCont[i]->Clone();
      fKeep->fCont[i] = fCont[i];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TClonesArray& TClonesArray::operator=(const TClonesArray& tc)
{
   if (this == &tc) return *this;

   if (fClass != tc.fClass) {
      Error("operator=", "cannot copy TClonesArray's when classes are different");
      return *this;
   }

   if (tc.fSize > fSize)
      Expand(TMath::Max(tc.fSize, GrowBy(fSize)));

   Int_t i;

   for (i = 0; i < fSize; i++)
      if (fKeep->fCont[i]) {
         R__ReleaseMemory(fClass,fKeep->fCont[i]);
         fKeep->fCont[i] = nullptr;
         fCont[i] = nullptr;
      }

   BypassStreamer(kTRUE);

   for (i = 0; i < tc.fSize; i++) {
      if (tc.fCont[i]) fKeep->fCont[i] = tc.fCont[i]->Clone();
      fCont[i] = fKeep->fCont[i];
   }

   fLast = tc.fLast;
   Changed();
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete a clones array.

TClonesArray::~TClonesArray()
{
   if (fKeep) {
      for (Int_t i = 0; i < fKeep->fSize; i++) {
         R__ReleaseMemory(fClass,fKeep->fCont[i]);
         fKeep->fCont[i] = nullptr;
      }
   }
   SafeDelete(fKeep);

   // Protect against erroneously setting of owner bit
   SetOwner(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// When the kBypassStreamer bit is set, the automatically
/// generated Streamer can call directly TClass::WriteBuffer.
/// Bypassing the Streamer improves the performance when writing/reading
/// the objects in the TClonesArray. However there is a drawback:
/// When a TClonesArray is written with split=0 bypassing the Streamer,
/// the StreamerInfo of the class in the array being optimized,
/// one cannot use later the TClonesArray with split>0. For example,
/// there is a problem with the following scenario:
///  1. A class Foo has a TClonesArray of Bar objects
///  2. The Foo object is written with split=0 to Tree T1.
///     In this case the StreamerInfo for the class Bar is created
///     in optimized mode in such a way that data members of the same type
///     are written as an array improving the I/O performance.
///  3. In a new program, T1 is read and a new Tree T2 is created
///     with the object Foo in split>1
///  4. When the T2 branch is created, the StreamerInfo for the class Bar
///     is created with no optimization (mandatory for the split mode).
///     The optimized Bar StreamerInfo is going to be used to read
///     the TClonesArray in T1. The result will be Bar objects with
///     data member values not in the right sequence.
/// The solution to this problem is to call BypassStreamer(kFALSE)
/// for the TClonesArray. In this case, the normal Bar::Streamer function
/// will be called. The Bar::Streamer function works OK independently
/// if the Bar StreamerInfo had been generated in optimized mode or not.

void TClonesArray::BypassStreamer(Bool_t bypass)
{
   if (bypass)
      SetBit(kBypassStreamer);
   else
      ResetBit(kBypassStreamer);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove empty slots from array.

void TClonesArray::Compress()
{
   Int_t j = 0, je = 0;

   TObject **tmp = new TObject* [fSize];

   for (Int_t i = 0; i < fSize; i++) {
      if (fCont[i]) {
         fCont[j] = fCont[i];
         fKeep->fCont[j] = fKeep->fCont[i];
         j++;
      } else {
         tmp[je] = fKeep->fCont[i];
         je++;
      }
   }

   fLast = j - 1;

   Int_t jf = 0;
   for ( ; j < fSize; j++) {
      fCont[j] = 0;
      fKeep->fCont[j] = tmp[jf];
      jf++;
   }

   delete [] tmp;

   R__ASSERT(je == jf);
}

////////////////////////////////////////////////////////////////////////////////
/// Get an object at index 'idx' that is guaranteed to have been constructed.
/// It might be either a freshly allocated object or one that had already been
/// allocated (and assumingly used).  In the later case, it is the callers
/// responsibility to insure that the object is returned to a known state,
/// usually by calling the Clear method on the TClonesArray.
///
/// Tests to see if the destructor has been called on the object.
/// If so, or if the object has never been constructed the class constructor is called using
/// New().  If not, return a pointer to the correct memory location.
/// This explicitly to deal with TObject classes that allocate memory
/// which will be reset (but not deallocated) in their Clear()
/// functions.

TObject *TClonesArray::ConstructedAt(Int_t idx)
{
   TObject *obj = (*this)[idx];
   if ( obj && obj->TestBit(TObject::kNotDeleted) ) {
      return obj;
   }
   return (fClass) ? static_cast<TObject*>(fClass->New(obj)) : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get an object at index 'idx' that is guaranteed to have been constructed.
/// It might be either a freshly allocated object or one that had already been
/// allocated (and assumingly used).  In the later case, the function Clear
/// will be called and passed the value of 'clear_options'
///
/// Tests to see if the destructor has been called on the object.
/// If so, or if the object has never been constructed the class constructor is called using
/// New().  If not, return a pointer to the correct memory location.
/// This explicitly to deal with TObject classes that allocate memory
/// which will be reset (but not deallocated) in their Clear()
/// functions.

TObject *TClonesArray::ConstructedAt(Int_t idx, Option_t *clear_options)
{
   TObject *obj = (*this)[idx];
   if ( obj && obj->TestBit(TObject::kNotDeleted) ) {
      obj->Clear(clear_options);
      return obj;
   }
   return (fClass) ? static_cast<TObject*>(fClass->New(obj)) : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Clear the clones array. Only use this routine when your objects don't
/// allocate memory since it will not call the object dtors.
/// However, if the class in the TClonesArray implements the function
/// Clear(Option_t *option) and if option = "C" the function Clear()
/// is called for all objects in the array. In the function Clear(), one
/// can delete objects or dynamic arrays allocated in the class.
/// This procedure is much faster than calling TClonesArray::Delete().
/// When the option starts with "C+", eg "C+xyz" the objects in the array
/// are in turn cleared with the option "xyz"

void TClonesArray::Clear(Option_t *option)
{
   if (option && option[0] == 'C') {
      const char *cplus = strstr(option,"+");
      if (cplus) {
         cplus = cplus + 1;
      } else {
         cplus = "";
      }
      Int_t n = GetEntriesFast();
      for (Int_t i = 0; i < n; i++) {
         TObject *obj = UncheckedAt(i);
         if (obj) {
            obj->Clear(cplus);
            obj->ResetBit( kHasUUID );
            obj->ResetBit( kIsReferenced );
            obj->SetUniqueID( 0 );
         }
      }
   }

   // Protect against erroneously setting of owner bit
   SetOwner(kFALSE);

   TObjArray::Clear();
}

////////////////////////////////////////////////////////////////////////////////
/// Clear the clones array. Use this routine when your objects allocate
/// memory (e.g. objects inheriting from TNamed or containing TStrings
/// allocate memory). If not you better use Clear() since if is faster.

void TClonesArray::Delete(Option_t *)
{
   if ( fClass->TestBit(TClass::kIsEmulation) ) {
      // In case of emulated class, we can not use the delete operator
      // directly, it would use the wrong destructor.
      for (Int_t i = 0; i < fSize; i++) {
         if (fCont[i] && fCont[i]->TestBit(kNotDeleted)) {
            fClass->Destructor(fCont[i],kTRUE);
         }
      }
   } else {
      for (Int_t i = 0; i < fSize; i++) {
         if (fCont[i] && fCont[i]->TestBit(kNotDeleted)) {
            fCont[i]->~TObject();
         }
      }
   }

   // Protect against erroneously setting of owner bit.
   SetOwner(kFALSE);

   TObjArray::Clear();
}

////////////////////////////////////////////////////////////////////////////////
/// Expand or shrink the array to newSize elements.

void TClonesArray::Expand(Int_t newSize)
{
   if (newSize < 0) {
      Error ("Expand", "newSize must be positive (%d)", newSize);
      return;
   }
   if (!fKeep) {
      Error("ExpandCreate", "Not initialized properly, fKeep is still a nullptr");
      return;
   }
   if (newSize == fSize)
      return;
   if (newSize < fSize) {
      // release allocated space in fKeep and set to 0 so
      // Expand() will shrink correctly
      for (int i = newSize; i < fSize; i++)
         if (fKeep->fCont[i]) {
            R__ReleaseMemory(fClass,fKeep->fCont[i]);
            fKeep->fCont[i] = nullptr;
         }
   }

   TObjArray::Expand(newSize);
   fKeep->Expand(newSize);
}

////////////////////////////////////////////////////////////////////////////////
/// Expand or shrink the array to n elements and create the clone
/// objects by calling their default ctor. If n is less than the current size
/// the array is shrunk and the allocated space is freed.
/// This routine is typically used to create a clonesarray into which
/// one can directly copy object data without going via the
/// "new (arr[i]) MyObj()" (i.e. the vtbl is already set correctly).

void TClonesArray::ExpandCreate(Int_t n)
{
   if (n < 0) {
      Error("ExpandCreate", "n must be positive (%d)", n);
      return;
   }
   if (!fKeep) {
      Error("ExpandCreate", "Not initialized properly, fKeep is still a nullptr");
      return;
   }
   if (n > fSize)
      Expand(TMath::Max(n, GrowBy(fSize)));

   Int_t i;
   for (i = 0; i < n; i++) {
      if (!fKeep->fCont[i]) {
         fKeep->fCont[i] = (TObject*)fClass->New();
      } else if (!fKeep->fCont[i]->TestBit(kNotDeleted)) {
         // The object has been deleted (or never initialized)
         fClass->New(fKeep->fCont[i]);
      }
      fCont[i] = fKeep->fCont[i];
   }

   for (i = n; i < fSize; i++)
      if (fKeep->fCont[i]) {
         R__ReleaseMemory(fClass,fKeep->fCont[i]);
         fKeep->fCont[i] = nullptr;
         fCont[i] = nullptr;
      }

   fLast = n - 1;
   Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Expand or shrink the array to n elements and create the clone
/// objects by calling their default ctor. If n is less than the current size
/// the array is shrunk but the allocated space is _not_ freed.
/// This routine is typically used to create a clonesarray into which
/// one can directly copy object data without going via the
/// "new (arr[i]) MyObj()" (i.e. the vtbl is already set correctly).
/// This is a simplified version of ExpandCreate used in the TTree mechanism.

void TClonesArray::ExpandCreateFast(Int_t n)
{
   Int_t oldSize = fKeep->GetSize();
   if (n > fSize)
      Expand(TMath::Max(n, GrowBy(fSize)));

   Int_t i;
   for (i = 0; i < n; i++) {
      if (i >= oldSize || !fKeep->fCont[i]) {
         fKeep->fCont[i] = (TObject*)fClass->New();
      } else if (!fKeep->fCont[i]->TestBit(kNotDeleted)) {
         // The object has been deleted (or never initialized)
         fClass->New(fKeep->fCont[i]);
      }
      fCont[i] = fKeep->fCont[i];
   }
   if (fLast >= n) {
      memset(fCont + n, 0, (fLast - n + 1) * sizeof(TObject*));
   }
   fLast = n - 1;
   Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Remove object at index idx.

TObject *TClonesArray::RemoveAt(Int_t idx)
{
   if (!BoundsOk("RemoveAt", idx)) return 0;

   int i = idx-fLowerBound;

   if (fCont[i] && fCont[i]->TestBit(kNotDeleted)) {
      fCont[i]->~TObject();
   }

   if (fCont[i]) {
      fCont[i] = 0;
      // recalculate array size
      if (i == fLast)
         do { fLast--; } while (fLast >= 0 && fCont[fLast] == 0);
      Changed();
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove object from array.

TObject *TClonesArray::Remove(TObject *obj)
{
   if (!obj) return 0;

   Int_t i = IndexOf(obj) - fLowerBound;

   if (i == -1) return 0;

   if (fCont[i] && fCont[i]->TestBit(kNotDeleted)) {
      fCont[i]->~TObject();
   }

   fCont[i] = 0;
   // recalculate array size
   if (i == fLast)
      do { fLast--; } while (fLast >= 0 && fCont[fLast] == 0);
   Changed();
   return obj;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove objects from index idx1 to idx2 included.

void TClonesArray::RemoveRange(Int_t idx1, Int_t idx2)
{
   if (!BoundsOk("RemoveRange", idx1)) return;
   if (!BoundsOk("RemoveRange", idx2)) return;

   idx1 -= fLowerBound;
   idx2 -= fLowerBound;

   Bool_t change = kFALSE;
   for (TObject **obj=fCont+idx1; obj<=fCont+idx2; obj++) {
      if (!*obj) continue;
      if ((*obj)->TestBit(kNotDeleted)) {
         (*obj)->~TObject();
      }
      *obj = 0;
      change = kTRUE;
   }

   // recalculate array size
   if (change) Changed();
   if (idx1 < fLast || fLast > idx2) return;
   do { fLast--; } while (fLast >= 0 && fCont[fLast] == 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Create an array of clone objects of class cl. The class must inherit from
/// TObject.
/// The second argument s indicates an approximate number of objects
/// that will be entered in the array. If more than s objects are entered,
/// the array will be automatically expanded.
///
/// NB: This function should not be called in the TClonesArray is already
///     initialized with a class.

void TClonesArray::SetClass(const TClass *cl, Int_t s)
{
   if (fKeep) {
      Error("SetClass", "TClonesArray already initialized with another class");
      return;
   }
   fClass = (TClass*)cl;
   if (!fClass) {
      MakeZombie();
      Error("SetClass", "called with a null pointer");
      return;
   }
   const char *classname = fClass->GetName();
   if (!fClass->IsTObject()) {
      MakeZombie();
      Error("SetClass", "%s does not inherit from TObject", classname);
      return;
   }
   if (fClass->GetBaseClassOffset(TObject::Class())!=0) {
      MakeZombie();
      Error("SetClass", "%s must inherit from TObject as the left most base class.", classname);
      return;
   }
   Int_t nch = strlen(classname)+2;
   char *name = new char[nch];
   snprintf(name,nch, "%ss", classname);
   SetName(name);
   delete [] name;

   fKeep = new TObjArray(s);

   BypassStreamer(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
///see TClonesArray::SetClass(const TClass*)

void TClonesArray::SetClass(const char *classname, Int_t s)
{
   SetClass(TClass::GetClass(classname),s);
}


////////////////////////////////////////////////////////////////////////////////
/// A TClonesArray is always the owner of the object it contains.
/// However the collection its inherits from (TObjArray) does not.
/// Hence this member function needs to be a nop for TClonesArray.

void TClonesArray::SetOwner(Bool_t /* enable */)
{
   // Nothing to be done.
}

////////////////////////////////////////////////////////////////////////////////
/// If objects in array are sortable (i.e. IsSortable() returns true
/// for all objects) then sort array.

void TClonesArray::Sort(Int_t upto)
{
   Int_t nentries = GetAbsLast()+1;
   if (nentries <= 0 || fSorted) return;
   for (Int_t i = 0; i < fSize; i++)
      if (fCont[i]) {
         if (!fCont[i]->IsSortable()) {
            Error("Sort", "objects in array are not sortable");
            return;
         }
      }

   QSort(fCont, fKeep->fCont, 0, TMath::Min(nentries, upto-fLowerBound));

   fLast   = -2;
   fSorted = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Write all objects in array to the I/O buffer. ATTENTION: empty slots
/// are also stored (using one byte per slot). If you don't want this
/// use a TOrdCollection or TList.

void TClonesArray::Streamer(TBuffer &b)
{
   // Important Note: if you modify this function, remember to also modify
   // TConvertClonesArrayToProxy accordingly

   Int_t   nobjects;
   char    nch;
   TString s, classv;
   UInt_t R__s, R__c;

   if (b.IsReading()) {
      Version_t v = b.ReadVersion(&R__s, &R__c);
      if (v == 3) {
         const Int_t kOldBypassStreamer = BIT(14);
         if (TestBit(kOldBypassStreamer)) BypassStreamer();
      }
      if (v > 2)
         TObject::Streamer(b);
      if (v > 1)
         fName.Streamer(b);
      s.Streamer(b);
      classv = s;
      Int_t clv = 0;
      Ssiz_t pos = s.Index(";");
      if (pos != kNPOS) {
         classv = s(0, pos);
         s = s(pos+1, s.Length()-pos-1);
         clv = s.Atoi();
      }
      TClass *cl = TClass::GetClass(classv);
      if (!cl) {
         Error("Streamer", "expecting class %s but it was not found by TClass::GetClass\n",
               classv.Data());
         b.CheckByteCount(R__s, R__c,TClonesArray::IsA());
         return;
      }

      b >> nobjects;
      if (nobjects < 0)
         nobjects = -nobjects;  // still there for backward compatibility
      b >> fLowerBound;
      if (fClass == 0) {
         fClass = cl;
         if (fKeep == 0) {
            fKeep  = new TObjArray(fSize);
            Expand(nobjects);
         }
      } else if (cl != fClass && classv == fClass->GetName()) {
         // If fClass' name is different from classv, the user has intentionally changed
         // the target class, so we must not override it.
         fClass = cl;
         //this case may happen when switching from an emulated class to the real class
         //may not be an error. fClass may point to a deleted object
         //Error("Streamer", "expecting objects of type %s, finding objects"
         //   " of type %s", fClass->GetName(), cl->GetName());
         //return;
      }

      // make sure there are enough slots in the fKeep array
      if (fKeep->GetSize() < nobjects)
         Expand(nobjects);

      //reset fLast. nobjects may be 0
      Int_t oldLast = fLast;
      fLast = nobjects-1;

      //TStreamerInfo *sinfo = fClass->GetStreamerInfo(clv);
      if (CanBypassStreamer() && !b.TestBit(TBuffer::kCannotHandleMemberWiseStreaming)) {
         for (Int_t i = 0; i < nobjects; i++) {
            if (!fKeep->fCont[i]) {
               fKeep->fCont[i] = (TObject*)fClass->New();
            } else if (!fKeep->fCont[i]->TestBit(kNotDeleted)) {
               // The object has been deleted (or never initialized)
               fClass->New(fKeep->fCont[i]);
            }

            fCont[i] = fKeep->fCont[i];
         }
         if (clv < 8 && classv == "TF1") {
            // To allow backward compatibility of TClonesArray of v5 TF1 objects
            // that were stored member-wise.
            TClonesArray temp("ROOT::v5::TF1Data");
            temp.ExpandCreate(nobjects);
            b.ReadClones(&temp, nobjects, clv);
            // And now covert the v5 into the current
            if (gClonesArrayTF1Updater)
               gClonesArrayTF1Updater(nobjects, temp.GetObjectRef(nullptr), this->GetObjectRef(nullptr));
         } else if (clv <= 8 && clv > 3 && clv != 6 && classv == "TFormula") {
            // To allow backwar compatibility of TClonesArray of v5 TF1 objects
            // that were stored member-wise.
            TClonesArray temp("ROOT::v5::TFormula");
            temp.ExpandCreate(nobjects);
            b.ReadClones(&temp, nobjects, clv);
            // And now covert the v5 into the current
            if (gClonesArrayTFormulaUpdater)
               gClonesArrayTFormulaUpdater(nobjects, temp.GetObjectRef(nullptr), this->GetObjectRef(nullptr));
         } else {
            // sinfo->ReadBufferClones(b,this,nobjects,-1,0);
            b.ReadClones(this, nobjects, clv);
         }
      } else {
         for (Int_t i = 0; i < nobjects; i++) {
            b >> nch;
            if (nch) {
               if (!fKeep->fCont[i])
                  fKeep->fCont[i] = (TObject*)fClass->New();
               else if (!fKeep->fCont[i]->TestBit(kNotDeleted)) {
                  // The object has been deleted (or never initialized)
                  fClass->New(fKeep->fCont[i]);
               }

               fCont[i] = fKeep->fCont[i];
               b.StreamObject(fKeep->fCont[i]);
            }
         }
      }
      for (Int_t i = TMath::Max(nobjects,0); i < oldLast+1; ++i) fCont[i] = 0;
      Changed();
      b.CheckByteCount(R__s, R__c,TClonesArray::IsA());
   } else {
      //Make sure TStreamerInfo is not optimized, otherwise it will not be
      //possible to support schema evolution in read mode.
      //In case the StreamerInfo has already been computed and optimized,
      //one must disable the option BypassStreamer
      b.ForceWriteInfoClones(this);

      // make sure the status of bypass streamer is part of the buffer
      // (bits in TObject), so that when reading the object the right
      // mode is used, independent of the method (e.g. written via
      // TMessage, received and stored to a file and then later read via
      // TBufferFile)
      Bool_t bypass = kFALSE;
      if (b.TestBit(TBuffer::kCannotHandleMemberWiseStreaming)) {
         bypass = CanBypassStreamer();
         BypassStreamer(kFALSE);
      }

      R__c = b.WriteVersion(TClonesArray::IsA(), kTRUE);
      TObject::Streamer(b);
      fName.Streamer(b);
      s.Form("%s;%d", fClass->GetName(), fClass->GetClassVersion());
      s.Streamer(b);
      nobjects = GetEntriesFast();
      b << nobjects;
      b << fLowerBound;
      if (CanBypassStreamer()) {
         b.WriteClones(this,nobjects);
      } else {
         for (Int_t i = 0; i < nobjects; i++) {
            if (!fCont[i]) {
               nch = 0;
               b << nch;
            } else {
               nch = 1;
               b << nch;
               b.StreamObject(fCont[i]);
            }
         }
      }
      b.SetByteCount(R__c, kTRUE);

      if (bypass)
         BypassStreamer();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to reserved area in which a new object of clones
/// class can be constructed. This operator should not be used for
/// lefthand side assignments, like a[2] = xxx. Only like,
/// new (a[2]) myClass, or xxx = a[2]. Of course right hand side usage
/// is only legal after the object has been constructed via the
/// new operator or via the New() method. To remove elements from
/// the clones array use Remove() or RemoveAt().

TObject *&TClonesArray::operator[](Int_t idx)
{
   if (idx < 0) {
      Error("operator[]", "out of bounds at %d in %lx", idx, (Long_t)this);
      return fCont[0];
   }
   if (!fClass) {
      Error("operator[]", "invalid class specified in TClonesArray ctor");
      return fCont[0];
   }
   if (idx >= fSize)
      Expand(TMath::Max(idx+1, GrowBy(fSize)));

   if (!fKeep->fCont[idx]) {
      fKeep->fCont[idx] = (TObject*) TStorage::ObjectAlloc(fClass->Size());
      // Reset the bit so that:
      //    obj = myClonesArray[i];
      //    obj->TestBit(TObject::kNotDeleted)
      // will behave correctly.
      // TObject::kNotDeleted is one of the higher bit that is not settable via the public
      // interface. But luckily we are its friend.
      fKeep->fCont[idx]->fBits &= ~kNotDeleted;
   }
   fCont[idx] = fKeep->fCont[idx];

   fLast = TMath::Max(idx, GetAbsLast());
   Changed();

   return fCont[idx];
}

////////////////////////////////////////////////////////////////////////////////
/// Return the object at position idx. Returns 0 if idx is out of bounds.

TObject *TClonesArray::operator[](Int_t idx) const
{
   if (idx < 0 || idx >= fSize) {
      Error("operator[]", "out of bounds at %d in %lx", idx, (Long_t)this);
      return 0;
   }

   return fCont[idx];
}

////////////////////////////////////////////////////////////////////////////////
/// Create an object of type fClass with the default ctor at the specified
/// index. Returns 0 in case of error.

TObject *TClonesArray::New(Int_t idx)
{
   if (idx < 0) {
      Error("New", "out of bounds at %d in %lx", idx, (Long_t)this);
      return 0;
   }
   if (!fClass) {
      Error("New", "invalid class specified in TClonesArray ctor");
      return 0;
   }

   return (TObject *)fClass->New(operator[](idx));
}

//______________________________________________________________________________
//
// The following functions are utilities implemented by Jason Detwiler
// (jadetwiler@lbl.gov)
//
////////////////////////////////////////////////////////////////////////////////
/// Directly move the object pointers from tc without cloning (copying).
/// This TClonesArray takes over ownership of all of tc's object
/// pointers. The tc array is left empty upon return.

void TClonesArray::AbsorbObjects(TClonesArray *tc)
{
   // tests
   if (tc == 0 || tc == this || tc->GetEntriesFast() == 0) return;
   AbsorbObjects(tc, 0, tc->GetEntriesFast() - 1);
}

////////////////////////////////////////////////////////////////////////////////
/// Directly move the range of object pointers from tc without cloning
/// (copying).
/// This TClonesArray takes over ownership of all of tc's object pointers
/// from idx1 to idx2. The tc array is re-arranged by return.

void TClonesArray::AbsorbObjects(TClonesArray *tc, Int_t idx1, Int_t idx2)
{
   // tests
   if (tc == 0 || tc == this || tc->GetEntriesFast() == 0) return;
   if (fClass != tc->fClass) {
      Error("AbsorbObjects", "cannot absorb objects when classes are different");
      return;
   }

   if (idx1 > idx2) {
      Error("AbsorbObjects", "range is not valid: idx1>idx2");
      return;
   }
   if (idx2 >= tc->GetEntriesFast()) {
      Error("AbsorbObjects", "range is not valid: idx2 out of bounds");
      return;
   }

   // cache the sorted status
   Bool_t wasSorted = IsSorted() && tc->IsSorted() &&
                      (Last() == 0 || Last()->Compare(tc->First()) == -1);

   // expand this
   Int_t oldSize = GetEntriesFast();
   Int_t newSize = oldSize + (idx2-idx1+1);
   if(newSize > fSize)
      Expand(newSize);

   // move
   for (Int_t i = idx1; i <= idx2; i++) {
      Int_t newindex = oldSize+i -idx1;
      fCont[newindex] = tc->fCont[i];
      R__ReleaseMemory(fClass,fKeep->fCont[newindex]);
      (*fKeep)[newindex] = (*(tc->fKeep))[i];
      tc->fCont[i] = 0;
      (*(tc->fKeep))[i] = 0;
   }

   // cleanup
   for (Int_t i = idx2+1; i < tc->GetEntriesFast(); i++) {
      tc->fCont[i-(idx2-idx1+1)] = tc->fCont[i];
      (*(tc->fKeep))[i-(idx2-idx1+1)] = (*(tc->fKeep))[i];
      tc->fCont[i] = 0;
      (*(tc->fKeep))[i] = 0;
   }
   tc->fLast = tc->GetEntriesFast() - 2 - (idx2 - idx1);
   fLast = newSize-1;
   if (!wasSorted)
      Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Sort multiple TClonesArrays simultaneously with this array.
/// If objects in array are sortable (i.e. IsSortable() returns true
/// for all objects) then sort array.

void TClonesArray::MultiSort(Int_t nTCs, TClonesArray** tcs, Int_t upto)
{
   Int_t nentries = GetAbsLast()+1;
   if (nentries <= 1 || fSorted) return;
   Bool_t sortedCheck = kTRUE;
   for (Int_t i = 0; i < fSize; i++) {
      if (fCont[i]) {
         if (!fCont[i]->IsSortable()) {
            Error("MultiSort", "objects in array are not sortable");
            return;
         }
      }
      if (sortedCheck && i > 1) {
         if (ObjCompare(fCont[i], fCont[i-1]) < 0) sortedCheck = kFALSE;
      }
   }
   if (sortedCheck) {
      fSorted = kTRUE;
      return;
   }

   for (int i = 0; i < nTCs; i++) {
      if (tcs[i] == this) {
         Error("MultiSort", "tcs[%d] = \"this\"", i);
         return;
      }
      if (tcs[i]->GetEntriesFast() != GetEntriesFast()) {
         Error("MultiSort", "tcs[%d] has length %d != length of this (%d)",
               i, tcs[i]->GetEntriesFast(), this->GetEntriesFast());
         return;
      }
   }

   int nBs = nTCs*2+1;
   TObject*** b = new TObject**[nBs];
   for (int i = 0; i < nTCs; i++) {
      b[2*i]   = tcs[i]->fCont;
      b[2*i+1] = tcs[i]->fKeep->fCont;
   }
   b[nBs-1] = fKeep->fCont;
   QSort(fCont, nBs, b, 0, TMath::Min(nentries, upto-fLowerBound));
   delete [] b;

   fLast = -2;
   fSorted = kTRUE;
}
