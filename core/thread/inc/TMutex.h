// @(#)root/thread:$Id$
// Author: Fons Rademakers   26/06/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMutex
#define ROOT_TMutex


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMutex                                                               //
//                                                                      //
// This class implements mutex locks. A mutex is a mutual exclusive     //
// lock. The actual work is done via the TMutexImp class (either        //
// TPosixMutex or TWin32Mutex).                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TVirtualMutex.h"
#include "TMutexImp.h"


class TMutex : public TVirtualMutex {

friend class TCondition;
friend class TThread;

private:
   TMutexImp  *fMutexImp;   // pointer to mutex implementation

   TMutex(const TMutex&) = delete;
   TMutex& operator=(const TMutex&) = delete;

public:
   TMutex(Bool_t recursive = kFALSE);
   virtual ~TMutex() { delete fMutexImp; }

   Int_t  Lock();
   Int_t  TryLock();
   Int_t  UnLock();
   Int_t  CleanUp();

   // Compatibility with standard library
   void lock() { TMutex::Lock(); }
   void unlock() { TMutex::UnLock(); }

   TVirtualMutex *Factory(Bool_t recursive = kFALSE);

   ClassDef(TMutex,0)  // Mutex lock class
};

#endif
