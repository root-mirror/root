// @(#)root/thread:$Id$
// Author: Fons Rademakers   14/11/06

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TWin32AtomicCount
#define ROOT_TWin32AtomicCount

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TWin32AtomicCount                                                    //
//                                                                      //
// Class providing atomic operations on a long. Setting, getting,       //
// incrementing and decrementing are atomic, thread safe, operations.   //
//                                                                      //
// This implementation uses the Win32 InterLocked API for locking.      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TAtomicCount
#error "Do not use TWin32AtomicCount.h directly. #include \"TAtomicCount.h\" instead."
#endif // ROOT_TAtomicCount

# include "Windows4Root.h"

class TAtomicCount {
private:
   Long_t fCnt;   // counter

   TAtomicCount(const TAtomicCount &) = delete;
   TAtomicCount &operator=(const TAtomicCount &) = delete;

public:
   explicit TAtomicCount(Long_t v) : fCnt(v) { }
   void operator++() { _InterlockedIncrement(&fCnt); }
   Long_t operator--() { return _InterlockedDecrement(&fCnt); }
   operator long() const { return static_cast<long const volatile &>(fCnt); }
   void Set(Long_t v) { fCnt = v; }
   Long_t Get() const { return static_cast<long const volatile &>(fCnt); }
};

#endif
