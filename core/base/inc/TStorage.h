// @(#)root/base:$Id$
// Author: Fons Rademakers   29/07/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TStorage
#define ROOT_TStorage


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TStorage                                                             //
//                                                                      //
// Storage manager.                                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

// #include "RConfigure.h" // included via Rtypes.h
#include "Rtypes.h"

typedef void (*FreeHookFun_t)(void*, void *addr, size_t);
typedef void *(*ReAllocFun_t)(void*, size_t);
typedef void *(*ReAllocCFun_t)(void*, size_t, size_t);
typedef char *(*ReAllocCharFun_t)(char*, size_t, size_t);


class TStorage {

private:
   static size_t         fgMaxBlockSize;       // largest block allocated
   static FreeHookFun_t  fgFreeHook;           // function called on free
   static void          *fgFreeHookData;       // data used by this function
   static ReAllocFun_t   fgReAllocHook;        // custom ReAlloc
   static ReAllocCFun_t  fgReAllocCHook;       // custom ReAlloc with length check
   static Bool_t         fgHasCustomNewDelete; // true if using ROOT's new/delete

   //----- Private bits, clients can only test but not change them
   enum {
      kIsOnHeap      = 0x01000000,    ///< object is on heap
   };

public:
   static const UInt_t   kObjectAllocMemValue = 0x99999999;
                                               // magic number for ObjectAlloc

public:
   virtual ~TStorage() { }

   static ULong_t        GetHeapBegin();
   static ULong_t        GetHeapEnd();
   static FreeHookFun_t  GetFreeHook();
   static void          *GetFreeHookData();
   static size_t         GetMaxBlockSize();
   static void          *Alloc(size_t size);
   static void           Dealloc(void *ptr);
   static void          *ReAlloc(void *vp, size_t size);
   static void          *ReAlloc(void *vp, size_t size, size_t oldsize);
   static char          *ReAllocChar(char *vp, size_t size, size_t oldsize);
   static Int_t         *ReAllocInt(Int_t *vp, size_t size, size_t oldsize);
   static void          *ObjectAlloc(size_t size);
   static void          *ObjectAllocArray(size_t size);
   static void          *ObjectAlloc(size_t size, void *vp);
   static void           ObjectDealloc(void *vp);
#ifdef R__SIZEDDELETE
   static void           ObjectDealloc(void *vp, size_t size);
#endif
   static void           ObjectDealloc(void *vp, void *ptr);

   static void EnterStat(size_t size, void *p);
   static void RemoveStat(void *p);
   static void PrintStatistics();
   static void SetMaxBlockSize(size_t size);
   static void SetFreeHook(FreeHookFun_t func, void *data);
   static void SetReAllocHooks(ReAllocFun_t func1, ReAllocCFun_t func2);
   static void SetCustomNewDelete();
   static void EnableStatistics(int size= -1, int ix= -1);

   static Bool_t HasCustomNewDelete();

   // only valid after call to a TStorage allocating method
   static void   AddToHeap(ULong_t begin, ULong_t end);
   static Bool_t IsOnHeap(void *p);

   static Bool_t FilledByObjectAlloc(volatile const UInt_t* const member);
   static void UpdateIsOnHeap(volatile const UInt_t &uniqueID, volatile UInt_t &bits);

   ClassDef(TStorage,0)  //Storage manager class
};

inline Bool_t TStorage::FilledByObjectAlloc(volatile const UInt_t *const member) {
   //called by TObject's constructor to determine if object was created by call to new

   // This technique is necessary as there is one stack per thread
   // and we can not rely on comparison with the current stack memory position.
   // Note that a false positive (this routine returning true for an object
   // created on the stack) requires the previous stack value to have been
   // set to exactly kObjectAllocMemValue at exactly the right position (i.e.
   // where this object's fUniqueID is located.
   // The consequence of a false positive will be visible if and only if
   //   the object is auto-added to a TDirectory (i.e. TTree, TH*, TGraph,
   //      TEventList) or explicitly added to the directory by the user
   // and
   //   the TDirectory (or TFile) object is created on the stack *before*
   //      the object.
   // The consequence would be that those objects would be deleted twice, once
   // by the TDirectory and once automatically when going out of scope
   // (and thus quite visible).  A false negative (which is not posible with
   // this implementation) would have been a silent memory leak.

   // This will be reported by valgrind as uninitialized memory reads for
   // object created on the stack, use $ROOTSYS/etc/valgrind-root.supp
R__INTENTIONALLY_UNINIT_BEGIN
   return *member == kObjectAllocMemValue;
R__INTENTIONALLY_UNINIT_END
}

// Assign the kIsOnHeap bit in 'bits' based on the pattern seen in uniqueID.
// See Storage::FilledByObjectAlloc for details.
// This routine is marked as inline with attribute noinline so that it never
// inlined and thus can be used in a valgrind suppression file to suppress
// the known/intentional uninitialized memory read but still be a 'quick'
// function call to avoid losing performance at object creation.
// Moving the function into the source file, results in doubling of the
// overhead (compared to inlining)
R__NEVER_INLINE void TStorage::UpdateIsOnHeap(volatile const UInt_t &uniqueID, volatile UInt_t &bits) {
   if (TStorage::FilledByObjectAlloc(&uniqueID))
      bits |= kIsOnHeap;
   else
      bits &= ~kIsOnHeap;
}


inline size_t TStorage::GetMaxBlockSize() { return fgMaxBlockSize; }

inline void TStorage::SetMaxBlockSize(size_t size) { fgMaxBlockSize = size; }

inline FreeHookFun_t TStorage::GetFreeHook() { return fgFreeHook; }

namespace ROOT {
namespace Internal {
using FreeIfTMapFile_t = bool(void*);
R__EXTERN FreeIfTMapFile_t *gFreeIfTMapFile;
R__EXTERN void *gMmallocDesc;
}
}

#endif
