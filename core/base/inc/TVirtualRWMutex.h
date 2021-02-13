// @(#)root/base:$Id$
// Author: Philippe Canal, 2017

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TVirtualRWMutex
#define ROOT_TVirtualRWMutex


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TVirtualRWMutex                                                      //
//                                                                      //
// This class implements a read-write mutex interface. The actual work  //
// is done via TRWSpinLock which is available as soon as the thread     //
// library is loaded.                                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TVirtualMutex.h"

#include <memory>

namespace ROOT {

class TVirtualRWMutex;

// Global mutex set in TThread::Init
// Use either R__READ_LOCKGUARD(ROOT::gCoreMutex);
//         or R__WRITE_LOCKGUARD(ROOT::gCoreMutex);
R__EXTERN TVirtualRWMutex *gCoreMutex;

class TVirtualRWMutex : public TVirtualMutex  {

public:
   // The following are opaque type and are never really declared
   // The specific implementation of TInterpreter will cast the
   // value of pointer to this types to the correct type (possibly
   // distinct from these)
   class Hint_t;

   /// \class State
   /// Earlier lock state as returned by `GetState()` that can be passed to
   /// `Restore()`
   struct State {
      virtual ~State(); // implemented in TVirtualMutex.cxx
   };

   struct StateAndRecurseCount {
      /// State of gCoreMutex when the first interpreter-related function was invoked.
      std::unique_ptr<ROOT::TVirtualRWMutex::State> fState;

      /// Interpreter-related functions will push the "entry" lock state to *this.
      /// Recursive calls will do that, too - but we must only forget about the lock
      /// state once this recursion count went to 0.
      Int_t fRecurseCount = 0;

      operator bool() const { return (bool)fState; }
   };

   /// \class StateDelta
   /// State as returned by `GetStateDelta()` that can be passed to
   /// `Restore()`
   struct StateDelta {
      virtual ~StateDelta(); // implemented in TVirtualMutex.cxx
   };

   virtual Hint_t *ReadLock() = 0;
   virtual void ReadUnLock(Hint_t *) = 0;
   virtual Hint_t *WriteLock() = 0;
   virtual void WriteUnLock(Hint_t *) = 0;

   Int_t Lock() override { WriteLock(); return 1; }
   Int_t TryLock() override { WriteLock(); return 1; }
   Int_t UnLock() override { WriteUnLock(nullptr); return 1; }
   Int_t CleanUp() override { WriteUnLock(nullptr); return 1; }

   virtual std::unique_ptr<State> GetStateBefore() = 0;
   virtual std::unique_ptr<StateDelta> Rewind(const State& earlierState) = 0;
   virtual void Apply(std::unique_ptr<StateDelta> &&delta) = 0;

   TVirtualRWMutex *Factory(Bool_t /*recursive*/ = kFALSE) override = 0;

   ClassDefOverride(TVirtualRWMutex, 0)  // Virtual mutex lock class
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TReadLockGuard                                                       //
//                                                                      //
// This class provides RW mutex resource management in a guaranteed and //
// exception safe way. Use like this:                                   //
// {                                                                    //
//    TReadLockGuard guard(mutex);                                      //
//    ... // read something                                             //
// }                                                                    //
// when guard goes out of scope the mutex is unlocked in the TLockGuard //
// destructor. The exception mechanism takes care of calling the dtors  //
// of local objects so it is exception safe.                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TReadLockGuard {

private:
   TVirtualRWMutex *const fMutex;
   TVirtualRWMutex::Hint_t *fHint;

   TReadLockGuard(const TReadLockGuard&) = delete;
   TReadLockGuard& operator=(const TReadLockGuard&) = delete;

public:
   TReadLockGuard(TVirtualRWMutex *mutex) : fMutex(mutex), fHint(nullptr) {
      if (fMutex) fHint = fMutex->ReadLock();
   }

   ~TReadLockGuard() { if (fMutex) fMutex->ReadUnLock(fHint); }

   ClassDefNV(TReadLockGuard,0)  // Exception safe read locking/unlocking of mutex
};

class TWriteLockGuard {

private:
   TVirtualRWMutex *const fMutex;
   TVirtualRWMutex::Hint_t *fHint;

   TWriteLockGuard(const TWriteLockGuard&) = delete;
   TWriteLockGuard& operator=(const TWriteLockGuard&) = delete;

public:
   TWriteLockGuard(TVirtualRWMutex *mutex) : fMutex(mutex), fHint(nullptr) {
      if (fMutex) fHint = fMutex->WriteLock();
   }

   ~TWriteLockGuard() { if (fMutex) fMutex->WriteUnLock(fHint); }

   ClassDefNV(TWriteLockGuard,0)  // Exception safe read locking/unlocking of mutex
};

} // namespace ROOT.

// Zero overhead macros in case not compiled with thread support
#if defined (_REENTRANT) || defined (WIN32)

#define R__READ_LOCKGUARD(mutex) ::ROOT::TReadLockGuard _R__UNIQUE_(R__readguard)(mutex)
#define R__READ_LOCKGUARD_NAMED(name,mutex) ::ROOT::TReadLockGuard _NAME2_(R__readguard,name)(mutex)

#define R__WRITE_LOCKGUARD(mutex) ::ROOT::TWriteLockGuard _R__UNIQUE_(R__readguard)(mutex)
#define R__WRITE_LOCKGUARD_NAMED(name,mutex) ::ROOT::TWriteLockGuard _NAME2_(R__readguard,name)(mutex)

#else

#define R__READ_LOCKGUARD(mutex) (void)mutex
#define R__READ_LOCKGUARD_NAMED(name,mutex) (void)mutex

#define R__WRITE_LOCKGUARD(mutex) (void)mutex
#define R__WRITE_LOCKGUARD_NAMED(name,mutex) (void)mutex

#endif


#endif
