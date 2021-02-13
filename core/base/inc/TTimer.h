// @(#)root/base:$Id$
// Author: Fons Rademakers   28/11/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TTimer
#define ROOT_TTimer


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TTimer                                                               //
//                                                                      //
// Handles synchronous and a-synchronous timer events. You can use      //
// this class in one of the following ways:                             //
//    - Sub-class TTimer and override the Notify() method.              //
//    - Re-implement the TObject::HandleTimer() method in your class    //
//      and pass a pointer to this object to timer, see the SetObject() //
//      method.                                                         //
//    - Pass an interpreter command to timer, see SetCommand() method.  //
//    - Create a TTimer, connect its Timeout() signal to the            //
//      appropriate methods. Then when the time is up it will emit a    //
//      Timeout() signal and call connected slots.                      //
//                                                                      //
//  Minimum timeout interval is defined in TSystem::ESysConstants as    //
//  kItimerResolution (currently 10 ms).                                //
//                                                                      //
//  Signal/slots example:                                               //
//       TTimer *timer = new TTimer();                                  //
//       timer->Connect("Timeout()", "myObjectClassName",               //
//                      myObject, "TimerDone()");                       //
//       timer->Start(2000, kTRUE);   // 2 seconds single-shot          //
//                                                                      //
//    // Timeout signal is emitted repeadetly with minimum timeout      //
//    // timer->Start(0, kFALSE);                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TSysEvtHandler.h"
#include "TTime.h"
#include "TString.h"



class TTimer : public TSysEvtHandler {

protected:
   TTime     fTime;        // time out time in ms
   TTime     fAbsTime;     // absolute time out time in ms
   Bool_t    fTimeout;     // true if timer has timed out
   Bool_t    fSync;        // true if synchrounous timer
   Bool_t    fIntSyscalls; // true is a-synchronous timer is to interrupt system calls
   UInt_t    fTimeID;      // the system ID of this timer (for WIN32)
   TObject  *fObject;      // object to be notified (if any)
   TString   fCommand;     // interpreter command to be executed

private:
   TTimer(const TTimer&) = delete;
   TTimer& operator=(const TTimer&) = delete;

public:
   TTimer(Long_t milliSec = 0, Bool_t mode = kTRUE);
   TTimer(TObject *obj, Long_t milliSec, Bool_t mode = kTRUE);
   TTimer(const char *command, Long_t milliSec, Bool_t mode = kTRUE);
   virtual ~TTimer() { Remove(); }

   Bool_t         CheckTimer(const TTime &now);
   const char    *GetCommand() const { return fCommand.Data(); }
   TObject       *GetObject() { return fObject; }
   TTime          GetTime() const { return fTime; }
   UInt_t         GetTimerID() { return fTimeID;}
   TTime          GetAbsTime() const { return fAbsTime; }
   Bool_t         HasTimedOut() const { return fTimeout; }
   Bool_t         IsSync() const { return fSync; }
   Bool_t         IsAsync() const { return !fSync; }
   Bool_t         IsInterruptingSyscalls() const { return fIntSyscalls; }
   virtual Bool_t Notify();
   void           Add() { TurnOn(); }
   void           Remove() { TurnOff(); }
   void           Reset();
   void           SetCommand(const char *command);
   void           SetObject(TObject *object);
   void           SetInterruptSyscalls(Bool_t set = kTRUE);
   void           SetTime(Long_t milliSec) { fTime = milliSec; }
   void           SetTimerID(UInt_t id = 0) { fTimeID = id; }
   virtual void   Start(Long_t milliSec = -1, Bool_t singleShot = kFALSE);
   virtual void   Stop() { TurnOff(); }
   virtual void   TurnOn();                         //*SIGNAL*
   virtual void   TurnOff();                        //*SIGNAL*
   virtual void   Timeout() { Emit("Timeout()"); }  //*SIGNAL*

   static void    SingleShot(Int_t milliSec, const char *receiver_class,
                             void *receiver, const char *method);

   ClassDef(TTimer,0)  //Handle timer event
};

#endif
