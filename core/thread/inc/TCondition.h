// @(#)root/thread:$Id$
// Author: Fons Rademakers   01/07/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TCondition
#define ROOT_TCondition


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCondition                                                           //
//                                                                      //
// This class implements a condition variable. Use a condition variable //
// to signal threads. The actual work is done via the TConditionImp     //
// class (either TPosixCondition or TWin32Condition).                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TConditionImp.h"

class TMutex;


class TCondition : public TObject {

friend class TThread;

private:
   TConditionImp  *fConditionImp;  // pointer to condition variable implementation
   TMutex         *fMutex;         // mutex used around Wait() and TimedWait()
   Bool_t          fPrivateMutex;  // is fMutex our private mutex

   TCondition(const TCondition&) = delete;
   TCondition& operator=(const TCondition&) = delete;

public:
   TCondition(TMutex *m = nullptr);
   virtual ~TCondition();

   TMutex *GetMutex() const;

   Int_t   Wait();
   Int_t   TimedWait(ULong_t secs, ULong_t nanoSecs);
   Int_t   TimedWaitRelative(ULong_t ms);
   Int_t   Signal() { if (fConditionImp) return fConditionImp->Signal(); return -1; }
   Int_t   Broadcast() { if (fConditionImp) return fConditionImp->Broadcast(); return -1; }

   ClassDef(TCondition,0)  // Condition variable class
};

#endif
