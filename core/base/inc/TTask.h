// @(#)root/base:$Id$
// Author: Rene Brun   02/09/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TTask
#define ROOT_TTask


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TTask                                                                //
//                                                                      //
// Base class for recursive execution of tasks.                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

#ifdef R__LESS_INCLUDES
class TList;
#else
#include "TList.h"
#endif

class TBrowser;


class TTask : public TNamed {

protected:
   TList        *fTasks;        //List of Tasks
   TString       fOption;       //Option specified in ExecuteTask
   Int_t         fBreakin;      //=1 if a break point set at task extry
   Int_t         fBreakout;     //=1 if a break point set at task exit
   Bool_t        fHasExecuted;  //True if task has executed
   Bool_t        fActive;       //true if task is active

   static TTask *fgBeginTask;    //pointer to task initiator
   static TTask *fgBreakPoint;   //pointer to current break point

private:

public:
   TTask();
   TTask(const char* name, const char *title);
   virtual ~TTask();
   TTask(const TTask &task);
   TTask& operator=(const TTask& tt);

   virtual void  Abort();  // *MENU*
   virtual void  Add(TTask *task);
   virtual void  Browse(TBrowser *b);
   virtual void  CleanTasks();
   virtual void  Clear(Option_t *option="");
   virtual void  Continue(); // *MENU*
   virtual void  Exec(Option_t *option);
   virtual void  ExecuteTask(Option_t *option="0");  // *MENU*
   virtual void  ExecuteTasks(Option_t *option);
   Int_t         GetBreakin() const { return fBreakin; }
   Int_t         GetBreakout() const { return fBreakout; }
         Bool_t  IsActive() const { return fActive; }
         Bool_t  IsFolder() const { return kTRUE; }
   virtual void  ls(Option_t *option="*") const;  // *MENU*
           void  SetActive(Bool_t active=kTRUE) { fActive = active; } // *TOGGLE*
           void  SetBreakin(Int_t breakin=1) { fBreakin = breakin; } // *TOGGLE*
           void  SetBreakout(Int_t breakout=1) { fBreakout = breakout; } // *TOGGLE*
   TList        *GetListOfTasks() const { return fTasks; }

   ClassDef(TTask,1)  //Base class for tasks
};

#endif
