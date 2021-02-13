// @(#)root/thread:$Id$
// Author: Fons Rademakers   14/11/06

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TAtomicCount
#define ROOT_TAtomicCount


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TAtomicCount                                                         //
//                                                                      //
// Class providing atomic operations on a long. Setting, getting,       //
// incrementing and decrementing are atomic, thread safe, operations.   //
//                                                                      //
//  TAtomicCount a(n);                                                  //
//                                                                      //
//    (n is convertible to long)                                        //
//                                                                      //
//    Effects: Constructs an TAtomicCount with an initial value of n.   //
//                                                                      //
//  long(a);                                                            //
//                                                                      //
//    Returns: (long) the current value of a.                           //
//                                                                      //
//  ++a;                                                                //
//                                                                      //
//    Effects: Atomically increments the value of a.                    //
//    Returns: nothing.                                                 //
//                                                                      //
//  --a;                                                                //
//                                                                      //
//    Effects: Atomically decrements the value of a.                    //
//    Returns: (long) zero if the new value of a is zero,               //
//             unspecified non-zero value otherwise                     //
//             (usually the new value).                                 //
//                                                                      //
//  a.Set(n);                                                           //
//                                                                      //
//    Effects: Set a to the value n.                                    //
//    Returns: nothing.                                                 //
//                                                                      //
//  a.Get();                                                            //
//                                                                      //
//    Returns: (long) the current value of a.                           //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "RtypesCore.h"
#include "RConfigure.h"

#if (defined(__GLIBCPP__) || defined(__GLIBCXX__)) && !defined(__CINT__)
#include "TAtomicCountGcc.h"
#elif defined(_WIN32) && !defined(__CINT__)
#include "TWin32AtomicCount.h"
#elif defined(R__HAS_PTHREAD) && !defined(__CINT__)
#include "TAtomicCountPthread.h"
#else
class TAtomicCount {
private:
   Long_t  fCnt;   // counter

   TAtomicCount(const TAtomicCount &) = delete;
   TAtomicCount &operator=(const TAtomicCount &) = delete;

public:
   explicit TAtomicCount(Long_t v) : fCnt(v) { }
   void operator++() { ++fCnt; }
   Long_t operator--() { return --fCnt; }
   operator long() const { return fCnt; }
   void Set(Long_t v) { fCnt = v; }
   Long_t Get() const { return fCnt; }
};
#endif

#endif
