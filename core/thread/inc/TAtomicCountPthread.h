// @(#)root/thread:$Id$
// Author: Fons Rademakers   14/11/06

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TAtomicCountPthread
#define ROOT_TAtomicCountPthread

#ifndef ROOT_TAtomicCount
# error This should be #included only by TAtomicCount.h. Please #include TAtomicCount.h.
#endif //ROOT_TAtomicCount

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TAtomicCountPthread                                                  //
//                                                                      //
// Class providing atomic operations on a long. Setting, getting,       //
// incrementing and decrementing are atomic, thread safe, operations.   //
//                                                                      //
// This implementation uses pthread mutexes for locking. This clearly   //
// is less efficient than the version using asm locking instructions    //
// as in TAtomicCountGcc.h, but better than nothing.                    //
//                                                                      //
// ATTENTION: Don't use this file directly, it is included by           //
//            TAtomicCount.h.                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <pthread.h>

class TAtomicCount {
private:
   Long_t                  fCnt;     // counter
   mutable pthread_mutex_t fMutex;   // mutex used to lock counter

   TAtomicCount(const TAtomicCount &) = delete;
   TAtomicCount &operator=(const TAtomicCount &) = delete;

   class LockGuard {
   private:
      pthread_mutex_t &fM;  // mutex to be guarded
   public:
      LockGuard(pthread_mutex_t &m): fM(m) { pthread_mutex_lock(&fM); }
      ~LockGuard() { pthread_mutex_unlock(&fM); }
   };

public:
   explicit TAtomicCount(Long_t v): fCnt(v) {
      pthread_mutex_init(&fMutex, 0);
   }

   ~TAtomicCount() { pthread_mutex_destroy(&fMutex); }

   void operator++() {
      LockGuard lock(fMutex);
      ++fCnt;
   }

   Long_t operator--() {
      LockGuard lock(fMutex);
      return --fCnt;
   }

   operator long() const {
      LockGuard lock(fMutex);
      return fCnt;
   }

   void Set(Long_t v) {
      LockGuard lock(fMutex);
      fCnt = v;
   }

   Long_t Get() const {
      LockGuard lock(fMutex);
      return fCnt;
   }
 };

#endif
