// @(#)root/thread:$Id$
// Author: Fons Rademakers   02/07/97 (Revised: G Ganis, Nov 2015)

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TSemaphore
#define ROOT_TSemaphore

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TSemaphore                                                           //
//                                                                      //
// This class implements a counting semaphore. Use a semaphore          //
// to synchronize threads.                                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <mutex>
#include <condition_variable>

#include "TObject.h"

class TSemaphore : public TObject {

private:
   std::mutex              fMutex;   // semaphore mutex
   std::condition_variable fCond;    // semaphore condition variable
   Int_t                   fValue;   // semaphore value
   UInt_t                  fWakeups; // wakeups

   TSemaphore(const TSemaphore &s) = delete;
   TSemaphore& operator=(const TSemaphore &s) = delete;

public:
   TSemaphore(Int_t initial = 1);
   virtual ~TSemaphore() { }

   Int_t  Wait();
   Int_t  Wait(Int_t millisec);
   Int_t  TryWait();
   Int_t  Post();

   ClassDef(TSemaphore, 0)  // Counting semaphore
};

#endif
