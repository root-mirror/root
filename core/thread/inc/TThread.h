// @(#)root/thread:$Id$
// Author: Fons Rademakers   02/07/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TThread
#define ROOT_TThread


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TThread                                                              //
//                                                                      //
// This class implements threads. A thread is an execution environment  //
// much lighter than a process. A single process can have multiple      //
// threads. The actual work is done via the TThreadImp class (either    //
// TPosixThread or TWin32Thread).                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TTimer.h"
#include <stdarg.h>

#ifdef R__LESS_INCLUDES
class TCondition;
#else
#include "TCondition.h"
#endif

class TMutex;
class TThreadImp;

class TThread : public TNamed {

friend class TThreadImp;
friend class TPosixThread;
friend class TThreadTimer;
friend class TThreadCleaner;
friend class TWin32Thread;
friend class TThreadTearDownGuard;
friend class TJoinHelper;

public:

   typedef void *(*VoidRtnFunc_t)(void *);
   typedef void  (*VoidFunc_t)(void *);

   enum EPriority {
      kLowPriority,
      kNormalPriority,
      kHighPriority
   };

   enum EState {
      kInvalidState,            // thread was not created properly
      kNewState,                // thread object exists but hasn't started
      kRunningState,            // thread is running
      kTerminatedState,         // thread has terminated but storage has not
                                // yet been reclaimed (i.e. waiting to be joined)
      kFinishedState,           // thread has finished
      kCancelingState,          // thread in process of canceling
      kCanceledState,           // thread has been canceled
      kDeletingState            // thread in process of deleting
   };

private:
   TThread       *fNext;                  // pointer to next thread
   TThread       *fPrev;                  // pointer to prev thread
   TThread      **fHolder;                // pointer to holder of this (delete only)
   EPriority      fPriority;              // thread priority
   EState         fState;                 // thread state
   EState         fStateComing;           // coming thread state
   Long_t         fId;                    // thread id
   Long_t         fHandle;                // Win32 thread handle
   Bool_t         fDetached;              // kTRUE if thread is Detached
   Bool_t         fNamed;                 // kTRUE if thread is Named
   VoidRtnFunc_t  fFcnRetn;               // void* start function of thread
   VoidFunc_t     fFcnVoid;               // void  start function of thread
   void          *fThreadArg;             // thread start function arguments
   void          *fClean;                 // support of cleanup structure
   char           fComment[100];          // thread specific state comment

   static TThreadImp      *fgThreadImp;   // static pointer to thread implementation
   static char  * volatile fgXAct;        // Action name to do by main thread
   static void ** volatile fgXArr;        // pointer to control array of void pointers for action
   static volatile Int_t   fgXAnb;        // size of array above
   static volatile Int_t   fgXArt;        // return XA flag
   static Long_t           fgMainId;      // thread id of main thread
   static TThread         *fgMain;        // pointer to chain of TThread's
   static TMutex          *fgMainMutex;   // mutex to protect chain of threads
   static TMutex          *fgXActMutex;   // mutex to protect XAction
   static TCondition      *fgXActCondi;   // condition for XAction

   // Private Member functions
   void           Constructor();
   void           SetComment(const char *txt = nullptr)
                     { fComment[0] = 0; if (txt) { strncpy(fComment, txt, 99); fComment[99] = 0; } }
   void           DoError(Int_t level, const char *location, const char *fmt, va_list va) const;
   void           ErrorHandler(int level, const char *location, const char *fmt, va_list ap) const;
   static void    Init();
   static void   *Function(void *ptr);
   static Int_t   XARequest(const char *xact, Int_t nb, void **ar, Int_t *iret);
   static void    AfterCancel(TThread *th);
   static void    **GetTls(Int_t k);

   TThread(const TThread&) = delete;
   TThread& operator=(const TThread&) = delete;

public:
   TThread(VoidRtnFunc_t fn, void *arg = nullptr, EPriority pri = kNormalPriority);
   TThread(VoidFunc_t fn, void *arg = nullptr, EPriority pri = kNormalPriority);
   TThread(const char *thname, VoidRtnFunc_t fn, void *arg = nullptr, EPriority pri = kNormalPriority);
   TThread(const char *thname, VoidFunc_t fn, void *arg = nullptr, EPriority pri = kNormalPriority);
   TThread(Long_t id = 0);
   virtual ~TThread();

   Int_t            Kill();
   Int_t            Run(void *arg = nullptr);
   void             SetPriority(EPriority pri);
   void             Delete(Option_t *option="") { TObject::Delete(option); }
   EPriority        GetPriority() const { return fPriority; }
   EState           GetState() const { return fState; }
   Long_t           GetId() const { return fId; }
   static void      Ps();
   static void      ps() { Ps(); }

   static void      Initialize();
   static Bool_t    IsInitialized();

   Long_t           Join(void **ret = nullptr);
   static Long_t    Join(Long_t id, void **ret = nullptr);

   static Int_t     Exit(void *ret = nullptr);
   static Int_t     Exists();
   static TThread  *GetThread(Long_t id);
   static TThread  *GetThread(const char *name);

   static Int_t     Lock();                  //User's lock of main mutex
   static Int_t     TryLock();               //User's try lock of main mutex
   static Int_t     UnLock();                //User's unlock of main mutex
   static TThread  *Self();
   static Long_t    SelfId();
   static Int_t     Sleep(ULong_t secs, ULong_t nanos = 0);
   static Int_t     GetTime(ULong_t *absSec, ULong_t *absNanoSec);

   static Int_t     Delete(TThread *&th);
   static void    **Tsd(void *dflt, Int_t k);

   // Cancellation
   // there are two types of TThread cancellation:
   //    DEFERRED     - Cancellation only in user provided cancel-points
   //    ASYNCHRONOUS - In any point
   //    DEFERRED is more safe, it is DEFAULT.
   static Int_t     SetCancelOn();
   static Int_t     SetCancelOff();
   static Int_t     SetCancelAsynchronous();
   static Int_t     SetCancelDeferred();
   static Int_t     CancelPoint();
   static Int_t     Kill(Long_t id);
   static Int_t     Kill(const char *name);
   static Int_t     CleanUpPush(void *free, void *arg = nullptr);
   static Int_t     CleanUpPop(Int_t exe = 0);
   static Int_t     CleanUp();

   // XActions
   static void      Printf(const char *fmt, ...)   // format and print
#if defined(__GNUC__) && !defined(__CINT__)
   __attribute__((format(printf, 1, 2)))
#endif
   ;
   static void      XAction();

   ClassDef(TThread,0)  // Thread class
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TThreadCleaner                                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TThreadCleaner {
public:
   TThreadCleaner() { }
   ~TThreadCleaner();
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TThreadTimer                                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TThreadTimer : public TTimer {
public:
   // if this time is less or equal to kItimerResolution, TUnixSystem::DispatchOneEvent i
   // can not exit and have its caller react to the other TTimer's actions (like the request
   // to stop the event loop) until there is another type of event.
   TThreadTimer(Long_t ms = kItimerResolution + 10);
   Bool_t Notify();
};

#endif
