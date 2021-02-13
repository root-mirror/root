// @(#)root/cont:$Id$
// Author: Rene Brun   28/09/2001

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TProcessID
\ingroup Base

A TProcessID identifies a ROOT job in a unique way in time and space.
The TProcessID title consists of a TUUID object which provides a globally
unique identifier (for more see TUUID.h).

A TProcessID is automatically created by the TROOT constructor.
When a TFile contains referenced objects (see TRef), the TProcessID
object is written to the file.
If a file has been written in multiple sessions (same machine or not),
a TProcessID is written for each session.
These objects are used by the class TRef to uniquely identified
any TObject pointed by a TRef.

When a referenced object is read from a file (its bit kIsReferenced is set),
this object is entered into the objects table of the corresponding TProcessID.
Each TFile has a list of TProcessIDs (see TFile::fProcessIDs) also
accessible via TProcessID::fgPIDs (for all files).
When this object is deleted, it is removed from the table via the cleanup
mechanism invoked by the TObject destructor.

Each TProcessID has a table (TObjArray *fObjects) that keeps track
of all referenced objects. If a referenced object has a fUniqueID set,
a pointer to this unique object may be found via fObjects->At(fUniqueID).
In the same way, when a TRef::GetObject is called, GetObject uses
its own fUniqueID to find the pointer to the referenced object.
See TProcessID::GetObjectWithID and PutObjectWithID.

When a referenced object is deleted, its slot in fObjects is set to null.
//
See also TProcessUUID: a specialized TProcessID to manage the single list
of TUUIDs.
*/

#include "TProcessID.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TExMap.h"
#include "TVirtualMutex.h"
#include "TError.h"
#include "snprintf.h"

TObjArray  *TProcessID::fgPIDs   = nullptr; //pointer to the list of TProcessID
TProcessID *TProcessID::fgPID    = nullptr; //pointer to the TProcessID of the current session
std::atomic_uint TProcessID::fgNumber(0); //Current referenced object instance count
TExMap     *TProcessID::fgObjPIDs = nullptr; //Table (pointer,pids)

ClassImp(TProcessID);

static std::atomic<TProcessID *> gIsValidCache;
using PIDCacheContent_t = std::pair<Int_t, TProcessID*>;
static std::atomic<PIDCacheContent_t *> gGetProcessWithUIDCache;

////////////////////////////////////////////////////////////////////////////////
/// Return hash value for this object.

static inline ULong_t Void_Hash(const void *ptr)
{
   return TString::Hash(&ptr, sizeof(void*));
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TProcessID::TProcessID()
{
   // MSVC doesn't support fSpinLock=ATOMIC_FLAG_INIT; in the class definition
   // and Apple LLVM version 7.3.0 (clang-703.0.31) warns about:
   // fLock(ATOMIC_FLAG_INIT)
   // ^~~~~~~~~~~~~~~~
   //    c++/v1/atomic:1779:26: note: expanded from macro 'ATOMIC_FLAG_INIT'
   //  #define ATOMIC_FLAG_INIT {false}
   // So reset the flag instead.
   std::atomic_flag_clear( &fLock );

   fCount = 0;
   fObjects = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TProcessID::~TProcessID()
{
   delete fObjects;
   fObjects = 0;

   TProcessID *This = this; // We need a referencable value for the 1st argument
   gIsValidCache.compare_exchange_strong(This, nullptr);

   auto current = gGetProcessWithUIDCache.load();
   if (current && current->second == this) {
      gGetProcessWithUIDCache.compare_exchange_strong(current, nullptr);
      delete current;
   }

   R__WRITE_LOCKGUARD(ROOT::gCoreMutex);
   fgPIDs->Remove(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Static function to add a new TProcessID to the list of PIDs.

TProcessID *TProcessID::AddProcessID()
{
   R__WRITE_LOCKGUARD(ROOT::gCoreMutex);

   if (fgPIDs && fgPIDs->GetEntriesFast() >= 65534) {
      if (fgPIDs->GetEntriesFast() == 65534) {
         ::Warning("TProcessID::AddProcessID","Maximum number of TProcessID (65535) is almost reached (one left).  TRef will stop being functional when the limit is reached.");
      } else {
         ::Fatal("TProcessID::AddProcessID","Maximum number of TProcessID (65535) has been reached.  TRef are not longer functional.");
      }
   }

   TProcessID *pid = new TProcessID();

   if (!fgPIDs) {
      fgPID  = pid;
      fgPIDs = new TObjArray(10);
      gROOT->GetListOfCleanups()->Add(fgPIDs);
   }
   UShort_t apid = fgPIDs->GetEntriesFast();
   pid->IncrementCount();

   fgPIDs->Add(pid);
   // if (apid == 0) for(int incr=0; incr < 65533; ++incr) fgPIDs->Add(0); // NOTE: DEBUGGING ONLY MUST BE REMOVED!
   char name[20];
   snprintf(name,20,"ProcessID%d",apid);
   pid->SetName(name);
   pid->SetUniqueID((UInt_t)apid);
   TUUID u;
   //apid = fgPIDs->GetEntriesFast();
   pid->SetTitle(u.AsString());
   return pid;
}

////////////////////////////////////////////////////////////////////////////////
/// static function returning the ID assigned to obj
/// If the object is not yet referenced, its kIsReferenced bit is set
/// and its fUniqueID set to the current number of referenced objects so far.

UInt_t TProcessID::AssignID(TObject *obj)
{
   R__WRITE_LOCKGUARD(ROOT::gCoreMutex);

   UInt_t uid = obj->GetUniqueID() & 0xffffff;
   if (obj == fgPID->GetObjectWithID(uid)) return uid;
   if (obj->TestBit(kIsReferenced)) {
      fgPID->PutObjectWithID(obj,uid);
      return uid;
   }
   if (fgNumber >= 16777215) {
      // This process id is 'full', we need to use a new one.
      fgPID = AddProcessID();
      fgNumber = 0;
      for(Int_t i = 0; i < fgPIDs->GetLast()+1; ++i) {
         TProcessID *pid = (TProcessID*)fgPIDs->At(i);
         if (pid && pid->fObjects && pid->fObjects->GetEntries() == 0) {
            pid->Clear();
         }
      }
   }
   fgNumber++;
   obj->SetBit(kIsReferenced);
   uid = fgNumber;
   // if (fgNumber<10) fgNumber = 16777213; // NOTE: DEBUGGING ONLY MUST BE REMOVED!
   if ( fgPID->GetUniqueID() < 255 ) {
      obj->SetUniqueID( (uid & 0xffffff) + (fgPID->GetUniqueID()<<24) );
   } else {
      obj->SetUniqueID( (uid & 0xffffff) + 0xff000000 /* 255 << 24 */ );
   }
   fgPID->PutObjectWithID(obj,uid);
   return uid;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize fObjects.

void TProcessID::CheckInit()
{
   if (!fObjects) {
       while (fLock.test_and_set(std::memory_order_acquire));  // acquire lock
       if (!fObjects) fObjects = new TObjArray(100);
       fLock.clear(std::memory_order_release);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// static function (called by TROOT destructor) to delete all TProcessIDs

void TProcessID::Cleanup()
{
   R__WRITE_LOCKGUARD(ROOT::gCoreMutex);

   fgPIDs->Delete();
   gROOT->GetListOfCleanups()->Remove(fgPIDs);
   delete fgPIDs;
   fgPIDs = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// delete the TObjArray pointing to referenced objects
/// this function is called by TFile::Close("R")

void TProcessID::Clear(Option_t *)
{
   if (GetUniqueID()>254 && fObjects && fgObjPIDs) {
      // We might have many references registered in the map
      for(Int_t i = 0; i < fObjects->GetSize(); ++i) {
         TObject *obj = fObjects->UncheckedAt(i);
         if (obj) {
            ULong64_t hash = Void_Hash(obj);
            fgObjPIDs->Remove(hash,(Long64_t)obj);
            (*fObjects)[i] = 0;
         }
      }
   }
   delete fObjects; fObjects = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// The reference fCount is used to delete the TProcessID
/// in the TFile destructor when fCount = 0

Int_t TProcessID::DecrementCount()
{
   fCount--;
   if (fCount < 0) fCount = 0;
   return fCount;
}

////////////////////////////////////////////////////////////////////////////////
/// static function returning a pointer to TProcessID number pid in fgPIDs

TProcessID *TProcessID::GetProcessID(UShort_t pid)
{
   return (TProcessID*)fgPIDs->At(pid);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the (static) number of process IDs.

UInt_t TProcessID::GetNProcessIDs()
{
   return fgPIDs ? fgPIDs->GetLast()+1 : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// static function returning a pointer to TProcessID with its pid
/// encoded in the highest byte of uid

TProcessID *TProcessID::GetProcessWithUID(UInt_t uid, const void *obj)
{

   Int_t pid = (uid>>24)&0xff;
   if (pid==0xff) {
      // Look up the pid in the table (pointer,pid)
      if (fgObjPIDs==0) return 0;
      ULong_t hash = Void_Hash(obj);

      R__READ_LOCKGUARD(ROOT::gCoreMutex);
      pid = fgObjPIDs->GetValue(hash,(Long_t)obj);
      return (TProcessID*)fgPIDs->At(pid);
   } else {
      auto current = gGetProcessWithUIDCache.load();
      if (current && current->first == pid)
         return current->second;

      R__READ_LOCKGUARD(ROOT::gCoreMutex);
      auto res = (TProcessID*)fgPIDs->At(pid);

      auto next = new PIDCacheContent_t(pid, res);
      auto old = gGetProcessWithUIDCache.exchange(next);
      delete old;

      return res;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// static function returning a pointer to TProcessID with its pid
/// encoded in the highest byte of obj->GetUniqueID()

TProcessID *TProcessID::GetProcessWithUID(const TObject *obj)
{
   return GetProcessWithUID(obj->GetUniqueID(),obj);
}

////////////////////////////////////////////////////////////////////////////////
/// static function returning the pointer to the session TProcessID

TProcessID *TProcessID::GetSessionProcessID()
{
   return fgPID;
}

////////////////////////////////////////////////////////////////////////////////
/// Increase the reference count to this object.

Int_t TProcessID::IncrementCount()
{
   CheckInit();
   ++fCount;
   return fCount;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the current referenced object count
/// fgNumber is incremented every time a new object is referenced

UInt_t TProcessID::GetObjectCount()
{
   return fgNumber;
}

////////////////////////////////////////////////////////////////////////////////
/// returns the TObject with unique identifier uid in the table of objects

TObject *TProcessID::GetObjectWithID(UInt_t uidd)
{
   Int_t uid = uidd & 0xffffff;  //take only the 24 lower bits

   if (fObjects==0 || uid >= fObjects->GetSize()) return 0;
   return fObjects->UncheckedAt(uid);
}

////////////////////////////////////////////////////////////////////////////////
/// static: returns pointer to current TProcessID

TProcessID *TProcessID::GetPID()
{
   return fgPID;
}

////////////////////////////////////////////////////////////////////////////////
/// static: returns array of TProcessIDs

TObjArray *TProcessID::GetPIDs()
{
   return fgPIDs;
}


////////////////////////////////////////////////////////////////////////////////
/// static function. return kTRUE if pid is a valid TProcessID

Bool_t TProcessID::IsValid(TProcessID *pid)
{
   if (gIsValidCache == pid)
      return kTRUE;

   R__READ_LOCKGUARD(ROOT::gCoreMutex);

   if (fgPIDs==0) return kFALSE;
   if (fgPIDs->IndexOf(pid) >= 0) {
      gIsValidCache = pid;
      return kTRUE;
   }
    if (pid == (TProcessID*)gROOT->GetUUIDs()) {
      gIsValidCache = pid;
      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// stores the object at the uid th slot in the table of objects
/// The object uniqued is set as well as its kMustCleanup bit

void TProcessID::PutObjectWithID(TObject *obj, UInt_t uid)
{
   R__LOCKGUARD_IMT(gROOTMutex); // Lock for parallel TTree I/O

   if (uid == 0) uid = obj->GetUniqueID() & 0xffffff;

   if (!fObjects) fObjects = new TObjArray(100);
   fObjects->AddAtAndExpand(obj,uid);

   obj->SetBit(kMustCleanup);
   if ( (obj->GetUniqueID()&0xff000000)==0xff000000 ) {
      // We have more than 255 pids we need to store this
      // pointer in the table(pointer,pid) since there is no
      // more space in fUniqueID
      if (fgObjPIDs==0) fgObjPIDs = new TExMap;
      ULong_t hash = Void_Hash(obj);

      // We use operator() rather than Add() because
      // if the address has already been registered, we want to
      // update it's uniqueID (this can easily happen when the
      // referenced object have been stored in a TClonesArray.
      (*fgObjPIDs)(hash, (Long_t)obj) = GetUniqueID();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// called by the object destructor
/// remove reference to obj from the current table if it is referenced

void TProcessID::RecursiveRemove(TObject *obj)
{
   if (!fObjects) return;
   if (!obj->TestBit(kIsReferenced)) return;
   UInt_t uid = obj->GetUniqueID() & 0xffffff;
   if (obj == GetObjectWithID(uid)) {
      R__WRITE_LOCKGUARD(ROOT::gCoreMutex);
      // Only attempt to remove from the map the items that are already
      // registered (because they are associated with a TProcessID with index 
      // greater than 255.  Attempting to remove an item that is not in the map
      // issues a Warning message.  
      if (fgObjPIDs && ((obj->GetUniqueID()&0xff000000)==0xff000000)) {
         ULong64_t hash = Void_Hash(obj);
         fgObjPIDs->Remove(hash,(Long64_t)obj);
      }
      (*fObjects)[uid] = 0; // Avoid recalculation of fLast (compared to ->RemoveAt(uid))
   }
}


////////////////////////////////////////////////////////////////////////////////
/// static function to set the current referenced object count
/// fgNumber is incremented every time a new object is referenced

void TProcessID::SetObjectCount(UInt_t number)
{
   fgNumber = number;
}
