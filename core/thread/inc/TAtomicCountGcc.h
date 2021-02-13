// @(#)root/thread:$Id$
// Author: Fons Rademakers   14/11/06

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TAtomicCountGcc
#define ROOT_TAtomicCountGcc


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TAtomicCountGcc                                                      //
//                                                                      //
// Class providing atomic operations on a long. Setting, getting,       //
// incrementing and decrementing are atomic, thread safe, operations.   //
//                                                                      //
// This implementation uses GNU libstdc++ v3 atomic primitives, see     //
// http://gcc.gnu.org/onlinedocs/porting/Thread-safety.html.            //
//                                                                      //
// ATTENTION: Don't use this file directly, it is included by           //
//            TAtomicCount.h.                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TAtomicCount
# error This should be #included only by TAtomicCount.h. Please #include TAtomicCount.h.
#endif //ROOT_TAtomicCount

#if __GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 2) || \
    (defined(__APPLE_CC__) && __APPLE_CC__ > 5000  && !defined(MAC_OS_X_VERSION_10_6))
#include <bits/atomicity.h>
#else
#include <ext/atomicity.h>
#endif

#if defined(__GLIBCXX__) // g++ 3.4+

using __gnu_cxx::__atomic_add;
using __gnu_cxx::__exchange_and_add;

#endif

class TAtomicCount {
private:
   mutable _Atomic_word fCnt;   // counter

   TAtomicCount(const TAtomicCount &) = delete;
   TAtomicCount &operator=(const TAtomicCount &) = delete;

public:
   explicit TAtomicCount(Long_t v) : fCnt(v) { }
   void operator++() { __atomic_add(&fCnt, 1); }
   Long_t operator--() { return __exchange_and_add(&fCnt, -1) - 1; }
   operator long() const { return __exchange_and_add(&fCnt, 0); }
   void Set(Long_t v) {
      fCnt = v;
#ifdef _GLIBCXX_WRITE_MEM_BARRIER
#if !(defined(__INTEL_COMPILER) && defined(__ia64__)) //ICC doesn't support inline asm on IA-64
      _GLIBCXX_WRITE_MEM_BARRIER;
#endif
#endif
   }
   Long_t Get() const { return __exchange_and_add(&fCnt, 0); }
};

#endif
