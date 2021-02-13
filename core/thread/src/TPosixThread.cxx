// @(#)root/thread:$Id$
// Author: Fons Rademakers   02/07/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TPosixThread                                                         //
//                                                                      //
// This class provides an interface to the posix thread routines.       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TPosixThread.h"

#include "TThread.h"

ClassImp(TPosixThread);


////////////////////////////////////////////////////////////////////////////////
/// Create a pthread. Returns 0 on success, otherwise an error number will
/// be returned.

Int_t TPosixThread::Run(TThread *th)
{
   int det;
   pthread_t id;
   pthread_attr_t *attr = new pthread_attr_t;

   pthread_attr_init(attr);

   // Set detach state
   det = (th->fDetached) ? PTHREAD_CREATE_DETACHED : PTHREAD_CREATE_JOINABLE;

   pthread_attr_setdetachstate(attr, det);

   // See e.g. https://developer.apple.com/library/mac/#qa/qa1419/_index.html
   // MacOS has only 512k of stack per thread; Linux has 2MB.
   const size_t requiredStackSize = 1024*1024*2;
   size_t stackSize = 0;
   if (!pthread_attr_getstacksize(attr, &stackSize)
       && stackSize < requiredStackSize) {
      pthread_attr_setstacksize(attr, requiredStackSize);
   }
   int ierr = pthread_create(&id, attr, &TThread::Function, th);
   if (!ierr) th->fId = (Long_t) id;

   pthread_attr_destroy(attr);
   delete attr;

   return ierr;
}

////////////////////////////////////////////////////////////////////////////////
/// Join  suspends  the  execution  of the calling thread until the
/// thread identified by th terminates, either by  calling  pthread_exit
/// or by being cancelled. Returns 0 on success, otherwise an error number will
/// be returned.

Int_t TPosixThread::Join(TThread *th, void **ret)
{
   return pthread_join((pthread_t) th->fId, ret);
}

////////////////////////////////////////////////////////////////////////////////
/// Terminates the execution of the calling thread. Return 0.

Int_t TPosixThread::Exit(void *ret)
{
   pthread_exit(ret);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Cancellation is the mechanism by which a thread can terminate the
/// execution of another thread. Returns 0 on success, otherwise an error
/// number will be returned.

Int_t TPosixThread::Kill(TThread *th)
{
   return pthread_cancel((pthread_t) th->fId);
}

////////////////////////////////////////////////////////////////////////////////
/// Turn off the cancellation state of the calling thread. Returns 0 on
/// success, otherwise an error number will be returned.

Int_t TPosixThread::SetCancelOff()
{
   return pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Turn on the cancellation state of the calling thread. Returns 0 on
/// success, otherwise an error number will be returned.

Int_t TPosixThread::SetCancelOn()
{
   return pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the cancellation response type of the calling thread to
/// asynchronous, i.e. cancel as soon as the cancellation request
/// is received.

Int_t TPosixThread::SetCancelAsynchronous()
{
   return pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the cancellation response type of the calling thread to
/// deferred, i.e. cancel only at next cancellation point.
/// Returns 0 on success, otherwise an error number will be returned.

Int_t TPosixThread::SetCancelDeferred()
{
   return pthread_setcanceltype (PTHREAD_CANCEL_DEFERRED, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Introduce an explicit cancellation point. Returns 0.

Int_t TPosixThread::CancelPoint()
{
   int istate;
   pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, &istate);
   pthread_testcancel();
   pthread_setcancelstate(istate, 0);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Add thread cleanup function.

Int_t TPosixThread::CleanUpPush(void **main, void *free, void *arg)
{
   // pthread_cleanup_push(free, arg);
   if (!free) Error("CleanUpPush", "cleanup rountine = 0");
   new TPosixThreadCleanUp(main, free, arg);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Pop thread cleanup function from stack.

Int_t TPosixThread::CleanUpPop(void **main,Int_t exe)
{
   //  pthread_cleanup_pop(exe); // happy pthread future

   if (!main || !*main) return 1;
   TPosixThreadCleanUp *l = (TPosixThreadCleanUp*)(*main);
   if (!l->fRoutine) Error("CleanUpPop", "cleanup routine = 0");
   if (exe && l->fRoutine) ((void (*)(void*))(l->fRoutine))(l->fArgument);
   *main = l->fNext;  delete l;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Default thread cleanup routine.

Int_t TPosixThread::CleanUp(void **main)
{
   if (gDebug > 0)
      Info("Cleanup", "cleanup 0x%lx", (Long_t)*main);
   while (!CleanUpPop(main, 1)) { }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the thread identifier for the calling thread.

Long_t TPosixThread::SelfId()
{
   return (Long_t) pthread_self();
}

//   Clean Up section. PTHREAD implementations of cleanup after cancel are
//   too different and often too bad. Temporary I invent my own bicycle.
//                                                              V.Perev.

////////////////////////////////////////////////////////////////////////////////
///cleanup function

TPosixThreadCleanUp::TPosixThreadCleanUp(void **main, void *routine, void *arg)
{
   fNext = (TPosixThreadCleanUp*)*main;
   fRoutine = routine; fArgument = arg;
   *main  = this;
}
