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

#ifndef ROOT_TReentrantRWLock
#define ROOT_TReentrantRWLock

#include "ThreadLocalStorage.h"
#include "ROOT/TSpinMutex.hxx"
#include "TVirtualRWMutex.h"

#include <atomic>
#include <condition_variable>
#include <thread>
#include <unordered_map>

namespace ROOT {
namespace Internal {
struct UniqueLockRecurseCount {
   using Hint_t = TVirtualRWMutex::Hint_t;

   struct LocalCounts {
      size_t fReadersCount = 0;
      bool fIsWriter = false;
   };
   size_t fWriteRecurse = 0; ///<! Number of re-entry in the lock by the same thread.

   UniqueLockRecurseCount();

   using local_t = LocalCounts*;

   local_t GetLocal(){
      TTHREAD_TLS_DECL(LocalCounts, gLocal);
      return &gLocal;
   }

   Hint_t *IncrementReadCount(local_t &local) {
      ++(local->fReadersCount);
      return reinterpret_cast<TVirtualRWMutex::Hint_t *>(&(local->fReadersCount));
   }

   template <typename MutexT>
   Hint_t *IncrementReadCount(local_t &local, MutexT &) {
      return IncrementReadCount(local);
   }

   Hint_t *DecrementReadCount(local_t &local) {
      --(local->fReadersCount);
      return reinterpret_cast<TVirtualRWMutex::Hint_t *>(&(local->fReadersCount));
   }

   template <typename MutexT>
   Hint_t *DecrementReadCount(local_t &local, MutexT &) {
      return DecrementReadCount(local);
   }

   void ResetReadCount(local_t &local, int newvalue) {
      local->fReadersCount = newvalue;
   }

   bool IsCurrentWriter(local_t &local) { return local->fIsWriter; }
   bool IsNotCurrentWriter(local_t &local) { return !local->fIsWriter; }

   void SetIsWriter(local_t &local)
   {
      // if (fWriteRecurse == std::numeric_limits<decltype(fWriteRecurse)>::max()) {
      //    ::Fatal("TRWSpinLock::WriteLock", "Too many recursions in TRWSpinLock!");
      // }
      ++fWriteRecurse;
      local->fIsWriter = true;
   }

   void DecrementWriteCount() { --fWriteRecurse; }

   void ResetIsWriter(local_t &local) { local->fIsWriter = false; }

   size_t &GetLocalReadersCount(local_t &local) { return local->fReadersCount; }
};

struct RecurseCounts {
   using Hint_t = TVirtualRWMutex::Hint_t;
   using ReaderColl_t = std::unordered_map<std::thread::id, size_t>;
   size_t fWriteRecurse =  0; ///<! Number of re-entry in the lock by the same thread.

   std::thread::id fWriterThread; ///<! Holder of the write lock
   ReaderColl_t fReadersCount;    ///<! Set of reader thread ids

   using local_t = std::thread::id;

   local_t GetLocal() const { return std::this_thread::get_id(); }

   Hint_t *IncrementReadCount(local_t &local) {
      auto &count = fReadersCount[local];
      ++(count);
      return reinterpret_cast<TVirtualRWMutex::Hint_t *>(&count);
   }

   template <typename MutexT>
   Hint_t *IncrementReadCount(local_t &local, MutexT &mutex)
   {
      std::unique_lock<MutexT> lock(mutex);
      return IncrementReadCount(local);
   }

   Hint_t *DecrementReadCount(local_t &local) {
      auto &count = fReadersCount[local];
      --count;
      return reinterpret_cast<TVirtualRWMutex::Hint_t *>(&count);
   }

   template <typename MutexT>
   Hint_t *DecrementReadCount(local_t &local, MutexT &mutex)
   {
      std::unique_lock<MutexT> lock(mutex);
      return DecrementReadCount(local);
   }

   void ResetReadCount(local_t &local, int newvalue) {
      fReadersCount[local] = newvalue;
   }

   bool IsCurrentWriter(local_t &local) const { return fWriterThread == local; }
   bool IsNotCurrentWriter(local_t &local) const { return fWriterThread != local; }

   void SetIsWriter(local_t &local)
   {
      // if (fWriteRecurse == std::numeric_limits<decltype(fWriteRecurse)>::max()) {
      //    ::Fatal("TRWSpinLock::WriteLock", "Too many recursions in TRWSpinLock!");
      // }
      ++fWriteRecurse;
      fWriterThread = local;
   }

   void DecrementWriteCount() { --fWriteRecurse; }

   void ResetIsWriter(local_t & /* local */) { fWriterThread = std::thread::id(); }

   size_t &GetLocalReadersCount(local_t &local) { return fReadersCount[local]; }


};
} // Internal

template <typename MutexT = ROOT::TSpinMutex, typename RecurseCountsT = Internal::RecurseCounts>
class TReentrantRWLock {
private:

   std::atomic<int> fReaders;           ///<! Number of readers
   std::atomic<int> fReaderReservation; ///<! A reader wants access
   std::atomic<int> fWriterReservation; ///<! A writer wants access
   std::atomic<bool> fWriter;           ///<! Is there a writer?
   MutexT fMutex;                       ///<! RWlock internal mutex
   std::condition_variable_any fCond;   ///<! RWlock internal condition variable

   RecurseCountsT fRecurseCounts;        ///<! Trackers for re-entry in the lock by the same thread.

   // size_t fWriteRecurse;                ///<! Number of re-entry in the lock by the same thread.

   // std::thread::id fWriterThread; ///<! Holder of the write lock
   // ReaderColl_t fReadersCount;    ///<! Set of reader thread ids

   void AssertReadCountLocIsFromCurrentThread(const size_t* presumedLocalReadersCount);

public:
   using State = TVirtualRWMutex::State;
   using StateDelta = TVirtualRWMutex::StateDelta;

   ////////////////////////////////////////////////////////////////////////
   /// Regular constructor.
   TReentrantRWLock() : fReaders(0), fReaderReservation(0), fWriterReservation(0), fWriter(false) {}

   TVirtualRWMutex::Hint_t *ReadLock();
   void ReadUnLock(TVirtualRWMutex::Hint_t *);
   TVirtualRWMutex::Hint_t *WriteLock();
   void WriteUnLock(TVirtualRWMutex::Hint_t *);

   std::unique_ptr<State> GetStateBefore();
   std::unique_ptr<StateDelta> Rewind(const State &earlierState);
   void Apply(std::unique_ptr<StateDelta> &&delta);
   };
} // end of namespace ROOT

#endif
