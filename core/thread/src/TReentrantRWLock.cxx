// @(#)root/thread:$Id$
// Authors: Enric Tejedor CERN  12/09/2016
//          Philippe Canal FNAL 12/09/2016

/*************************************************************************
 * Copyright (C) 1995-2016, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TReentrantRWLock
    \brief An implementation of a reentrant read-write lock with a
           configurable internal mutex/lock (default Spin Lock).

This class provides an implementation of a reentrant read-write lock
that uses an internal lock and a condition variable to synchronize
readers and writers when necessary.

The implementation allows a single reader to take the write lock without
releasing the reader lock.  It also allows the writer to take a read lock.
In other word, the lock is re-entrant for both reading and writing.

The implementation tries to make faster the scenario when readers come
and go but there is no writer. In that case, readers will not pay the
price of taking the internal spin lock.

Moreover, this RW lock tries to be fair with writers, giving them the
possibility to claim the lock and wait for only the remaining readers,
thus preventing starvation.
*/

#include "ROOT/TReentrantRWLock.hxx"
#include "ROOT/TSpinMutex.hxx"
#include "TMutex.h"
#include "TError.h"
#include <assert.h>

using namespace ROOT;

#ifdef NDEBUG
# define R__MAYBE_AssertReadCountLocIsFromCurrentThread(READERSCOUNTLOC)
# define R__MAYBE_ASSERT_WITH_LOCAL_LOCK(where,msg,what)
#else
# define R__MAYBE_AssertReadCountLocIsFromCurrentThread(READERSCOUNTLOC) \
   AssertReadCountLocIsFromCurrentThread(READERSCOUNTLOC)
#define R__MAYBE_ASSERT_WITH_LOCAL_LOCK(where, msg, what) \
   {                                                      \
      std::unique_lock<MutexT> lock(fMutex);              \
      auto local = fRecurseCounts.GetLocal();             \
      if (!(what))                                        \
         Error(where, "%s", msg);                         \
   }
#endif

Internal::UniqueLockRecurseCount::UniqueLockRecurseCount()
{
   static bool singleton = false;
   if (singleton) {
      ::Fatal("UniqueLockRecurseCount Ctor", "Only one TReentrantRWLock using a UniqueLockRecurseCount is allowed.");
   }
   singleton = true;
}


////////////////////////////////////////////////////////////////////////////
/// Acquire the lock in read mode.
template <typename MutexT, typename RecurseCountsT>
TVirtualRWMutex::Hint_t *TReentrantRWLock<MutexT, RecurseCountsT>::ReadLock()
{
   ++fReaderReservation;

   // if (fReaders == std::numeric_limits<decltype(fReaders)>::max()) {
   //    ::Fatal("TRWSpinLock::WriteLock", "Too many recursions in TRWSpinLock!");
   // }

   auto local = fRecurseCounts.GetLocal();

   TVirtualRWMutex::Hint_t *hint = nullptr;

   if (!fWriter) {
      // There is no writer, go freely to the critical section
      ++fReaders;
      --fReaderReservation;

      hint = fRecurseCounts.IncrementReadCount(local, fMutex);

   } else if (fRecurseCounts.IsCurrentWriter(local)) {

      --fReaderReservation;
      // This can run concurrently with another thread trying to get
      // the read lock and ending up in the next section ("Wait for writers, if any")
      // which need to also get the local readers count and thus can
      // modify the map.
      hint = fRecurseCounts.IncrementReadCount(local, fMutex);
      ++fReaders;

   } else {
      // A writer claimed the RW lock, we will need to wait on the
      // internal lock
      --fReaderReservation;

      std::unique_lock<MutexT> lock(fMutex);

      // Wait for writers, if any
      if (fWriter && fRecurseCounts.IsNotCurrentWriter(local)) {
         auto readerCount = fRecurseCounts.GetLocalReadersCount(local);
         if (readerCount == 0)
            fCond.wait(lock, [this] { return !fWriter; });
         // else
         //   There is a writer **but** we have outstanding readers
         //   locks, this must mean that the writer is actually
         //   waiting on this thread to release its read locks.
         //   This can be done in only two ways:
         //     * request the writer lock
         //     * release the reader lock
         //   Either way, this thread needs to proceed to
         //   be able to reach a point whether it does one
         //   of the two.
      }

      hint = fRecurseCounts.IncrementReadCount(local);

      // This RW lock now belongs to the readers
      ++fReaders;

      lock.unlock();
   }

   return hint;
}

//////////////////////////////////////////////////////////////////////////
/// Release the lock in read mode.
template <typename MutexT, typename RecurseCountsT>
void TReentrantRWLock<MutexT, RecurseCountsT>::ReadUnLock(TVirtualRWMutex::Hint_t *hint)
{
   size_t *localReaderCount;
   if (!hint) {
      // This should be very rare.
      auto local = fRecurseCounts.GetLocal();
      std::lock_guard<MutexT> lock(fMutex);
      localReaderCount = &(fRecurseCounts.GetLocalReadersCount(local));
   } else {
      localReaderCount = reinterpret_cast<size_t*>(hint);
   }

   --fReaders;
   if (fWriterReservation && fReaders == 0) {
      // We still need to lock here to prevent interleaving with a writer
      std::lock_guard<MutexT> lock(fMutex);

      --(*localReaderCount);

      // Make sure you wake up a writer, if any
      // Note: spurrious wakeups are okay, fReaders
      // will be checked again in WriteLock
      fCond.notify_all();
   } else {

      --(*localReaderCount);
   }
}

//////////////////////////////////////////////////////////////////////////
/// Acquire the lock in write mode.
template <typename MutexT, typename RecurseCountsT>
TVirtualRWMutex::Hint_t *TReentrantRWLock<MutexT, RecurseCountsT>::WriteLock()
{
   ++fWriterReservation;

   std::unique_lock<MutexT> lock(fMutex);

   auto local = fRecurseCounts.GetLocal();

   // Release this thread's reader lock(s)
   auto &readerCount = fRecurseCounts.GetLocalReadersCount(local);
   TVirtualRWMutex::Hint_t *hint = reinterpret_cast<TVirtualRWMutex::Hint_t *>(&readerCount);

   fReaders -= readerCount;

   // Wait for other writers, if any
   if (fWriter && fRecurseCounts.IsNotCurrentWriter(local)) {
      if (readerCount && fReaders == 0) {
         // we decrease fReaders to zero, let's wake up the
         // other writer.
         fCond.notify_all();
      }
      fCond.wait(lock, [this] { return !fWriter; });
   }

   // Claim the lock for this writer
   fWriter = true;
   fRecurseCounts.SetIsWriter(local);

   // Wait until all reader reservations finish
   while (fReaderReservation) {
   };

   // Wait for remaining readers
   fCond.wait(lock, [this] { return fReaders == 0; });

   // Restore this thread's reader lock(s)
   fReaders += readerCount;

   --fWriterReservation;

   lock.unlock();

   return hint;
}

//////////////////////////////////////////////////////////////////////////
/// Release the lock in write mode.
template <typename MutexT, typename RecurseCountsT>
void TReentrantRWLock<MutexT, RecurseCountsT>::WriteUnLock(TVirtualRWMutex::Hint_t *)
{
   // We need to lock here to prevent interleaving with a reader
   std::lock_guard<MutexT> lock(fMutex);

   if (!fWriter || fRecurseCounts.fWriteRecurse == 0) {
      Error("TReentrantRWLock::WriteUnLock", "Write lock already released for %p", this);
      return;
   }

   fRecurseCounts.DecrementWriteCount();

   if (!fRecurseCounts.fWriteRecurse) {
      fWriter = false;

      auto local = fRecurseCounts.GetLocal();
      fRecurseCounts.ResetIsWriter(local);

      // Notify all potential readers/writers that are waiting
      fCond.notify_all();
   }
}
namespace {
template <typename MutexT, typename RecurseCountsT>
struct TReentrantRWLockState: public TVirtualRWMutex::State {
   size_t *fReadersCountLoc = nullptr;
   int fReadersCount = 0;
   size_t fWriteRecurse = 0;
};

template <typename MutexT, typename RecurseCountsT>
struct TReentrantRWLockStateDelta: public TVirtualRWMutex::StateDelta {
   size_t *fReadersCountLoc = nullptr;
   int fDeltaReadersCount = 0;
   int fDeltaWriteRecurse = 0;
};
}

//////////////////////////////////////////////////////////////////////////
/// Get the lock state before the most recent write lock was taken.

template <typename MutexT, typename RecurseCountsT>
std::unique_ptr<TVirtualRWMutex::State>
TReentrantRWLock<MutexT, RecurseCountsT>::GetStateBefore()
{
   using State_t = TReentrantRWLockState<MutexT, RecurseCountsT>;

   if (!fWriter) {
      Error("TReentrantRWLock::GetStateBefore()", "Must be write locked!");
      return nullptr;
   }

   auto local = fRecurseCounts.GetLocal();
   if (fRecurseCounts.IsNotCurrentWriter(local)) {
      Error("TReentrantRWLock::GetStateBefore()", "Not holding the write lock!");
      return nullptr;
   }

   std::unique_ptr<State_t> pState(new State_t);
   {
      std::lock_guard<MutexT> lock(fMutex);
      pState->fReadersCountLoc = &(fRecurseCounts.GetLocalReadersCount(local));
   }
   pState->fReadersCount = *(pState->fReadersCountLoc);
   // *Before* the most recent write lock (that is required by GetStateBefore())
   // was taken, the write recursion level was `fWriteRecurse - 1`
   pState->fWriteRecurse = fRecurseCounts.fWriteRecurse - 1;

#if __GNUC__ < 7
   // older version of gcc can not convert implicitly from
   // unique_ptr of derived to unique_ptr of base
   using BaseState_t = TVirtualRWMutex::State;
   return std::unique_ptr<BaseState_t>(pState.release());
#else
   return pState;
#endif
}

//////////////////////////////////////////////////////////////////////////
/// Rewind to an earlier mutex state, returning the delta.

template <typename MutexT, typename RecurseCountsT>
std::unique_ptr<TVirtualRWMutex::StateDelta>
TReentrantRWLock<MutexT, RecurseCountsT>::Rewind(const State &earlierState) {
   using State_t = TReentrantRWLockState<MutexT, RecurseCountsT>;
   using StateDelta_t = TReentrantRWLockStateDelta<MutexT, RecurseCountsT>;
   auto& typedState = static_cast<const State_t&>(earlierState);

   R__MAYBE_AssertReadCountLocIsFromCurrentThread(typedState.fReadersCountLoc);

   std::unique_ptr<StateDelta_t> pStateDelta(new StateDelta_t);
   pStateDelta->fReadersCountLoc = typedState.fReadersCountLoc;
   pStateDelta->fDeltaReadersCount = *typedState.fReadersCountLoc - typedState.fReadersCount;
   pStateDelta->fDeltaWriteRecurse = fRecurseCounts.fWriteRecurse - typedState.fWriteRecurse;

   if (pStateDelta->fDeltaReadersCount < 0) {
      Error("TReentrantRWLock::Rewind", "Inconsistent read lock count!");
      return nullptr;
   }

   if (pStateDelta->fDeltaWriteRecurse < 0) {
      Error("TReentrantRWLock::Rewind", "Inconsistent write lock count!");
      return nullptr;
   }

   auto hint = reinterpret_cast<TVirtualRWMutex::Hint_t *>(typedState.fReadersCountLoc);
   if (pStateDelta->fDeltaWriteRecurse != 0) {
      R__MAYBE_ASSERT_WITH_LOCAL_LOCK("TReentrantRWLock::Rewind",
                                      "Lock rewinded from a thread that does not own the Write lock",
                                      fWriter && !fRecurseCounts.IsNotCurrentWriter(local));  // local defined in macro

      // Claim a recurse-state +1 to be able to call Unlock() below.
      fRecurseCounts.fWriteRecurse = typedState.fWriteRecurse + 1;
      // Release this thread's write lock
      WriteUnLock(hint);
   } else {
      Error("TReentrantRWLock::Rewind", "has been called with no write lock held.");
   }

   if (pStateDelta->fDeltaReadersCount != 0) {
      // Claim a recurse-state +1 to be able to call Unlock() below.
      *typedState.fReadersCountLoc = typedState.fReadersCount + 1;
      // ReadersCount and Readers can be read/updated before actually holding the R/W lock.
      // (See WriteLock and ReadLock methods) but fReaders is an atomic, so that should be fine
      // but that's weird, that does not account to other change in fReaders during between
      // the snapshot and the rewind ... humm unless the lock held is a WriteLock
      // (the actual use case) in which case there is no other thread that can update fReaders
      // and we also assume that the "user code" is balanced and release all read locks it takes.
      fReaders = typedState.fReadersCount + 1;
      // Release this thread's reader lock(s)
      ReadUnLock(hint);
   }
   // else earlierState and *this are identical!

   return std::unique_ptr<TVirtualRWMutex::StateDelta>(std::move(pStateDelta));
}

//////////////////////////////////////////////////////////////////////////
/// Re-apply a delta.

template <typename MutexT, typename RecurseCountsT>
void TReentrantRWLock<MutexT, RecurseCountsT>::Apply(std::unique_ptr<StateDelta> &&state) {
   if (!state) {
      Error("TReentrantRWLock::Apply", "Cannot apply empty delta!");
      return;
   }

   using StateDelta_t = TReentrantRWLockStateDelta<MutexT, RecurseCountsT>;
   const StateDelta_t* typedDelta = static_cast<const StateDelta_t*>(state.get());

   if (typedDelta->fDeltaWriteRecurse < 0) {
      Error("TReentrantRWLock::Apply", "Negative write recurse count delta!");
      return;
   }
   if (typedDelta->fDeltaReadersCount < 0) {
      Error("TReentrantRWLock::Apply", "Negative read count delta!");
      return;
   }
   R__MAYBE_AssertReadCountLocIsFromCurrentThread(typedDelta->fReadersCountLoc);

   if (typedDelta->fDeltaWriteRecurse != 0) {
      WriteLock();
      fRecurseCounts.fWriteRecurse += typedDelta->fDeltaWriteRecurse - 1;
   }
   if (typedDelta->fDeltaReadersCount != 0) {
      ReadLock();
      // "- 1" due to ReadLock() above.
      fReaders += typedDelta->fDeltaReadersCount - 1;
      *typedDelta->fReadersCountLoc += typedDelta->fDeltaReadersCount - 1;
   }
}

//////////////////////////////////////////////////////////////////////////
/// Assert that presumedLocalReadersCount really matches the local read count.
/// Print an error message if not.

template <typename MutexT, typename RecurseCountsT>
void TReentrantRWLock<MutexT, RecurseCountsT>::AssertReadCountLocIsFromCurrentThread(const size_t* presumedLocalReadersCount)
{
   auto local = fRecurseCounts.GetLocal();
   size_t* localReadersCount;
   {
      std::lock_guard<MutexT> lock(fMutex);
      localReadersCount = &(fRecurseCounts.GetLocalReadersCount(local));
   }
   if (localReadersCount != presumedLocalReadersCount) {
      Error("TReentrantRWLock::AssertReadCountLocIsFromCurrentThread", "ReadersCount is from different thread!");
   }
}


namespace ROOT {
template class TReentrantRWLock<ROOT::TSpinMutex, ROOT::Internal::RecurseCounts>;
template class TReentrantRWLock<TMutex, ROOT::Internal::RecurseCounts>;
template class TReentrantRWLock<std::mutex, ROOT::Internal::RecurseCounts>;

template class TReentrantRWLock<ROOT::TSpinMutex, ROOT::Internal::UniqueLockRecurseCount>;
template class TReentrantRWLock<TMutex, ROOT::Internal::UniqueLockRecurseCount>;
template class TReentrantRWLock<std::mutex, ROOT::Internal::UniqueLockRecurseCount>;
}
