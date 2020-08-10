// @(#)root/base:$Id$
// Author: Fons Rademakers   14/02/2002

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TVirtualMutex
\ingroup Base

This class implements a mutex interface. The actual work is done via
TMutex which is available as soon as the thread library is loaded.

and

TLockGuard

This class provides mutex resource management in a guaranteed and
exception safe way. Use like this:
~~~ {.cpp}
{
   TLockGuard guard(mutex);
   ... // do something
}
~~~
when guard goes out of scope the mutex is unlocked in the TLockGuard
destructor. The exception mechanism takes care of calling the dtors
of local objects so it is exception safe.
*/

#include "TVirtualMutex.h"
#include "TVirtualRWMutex.h"

ClassImp(TVirtualMutex);
ClassImp(TLockGuard);

// Global mutex set in TThread::Init protecting creation
// of other (preferably local) mutexes. Note that in this
// concept gGlobalMutex must be used in TStorage to prevent
// lockup of the system (see TMutex::Factory)
TVirtualMutex *gGlobalMutex = 0;

// From TVirtualRWMutex.h:
ROOT::TVirtualRWMutex::State::~State() = default;
ROOT::TVirtualRWMutex::StateDelta::~StateDelta() = default;
