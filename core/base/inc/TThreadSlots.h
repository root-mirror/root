// @(#)root/base:$Id$
// Author: Philippe Canal   09/30/2011

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TThreadSlots
#define ROOT_TThreadSlots

#include "RtypesCore.h"

namespace ROOT {

   enum EThreadSlotReservation {
      // Describe the system wide slot pre-allocation in the TThread
      // 'special data' storage array ; meant to be used as thread local
      // storage.  (See TThread::Tsd)
      //
      // Slot 0 through 19 can be used for user application
      // Slot 20 and above are reserved for the global system
      kMaxUserThreadSlot   = 20,

      // Slot reserved by ROOT's packages.
      kPadThreadSlot       = 20,
      kClassThreadSlot     = 21,
      kDirectoryThreadSlot = 22,
      kFileThreadSlot      = 23,
      kPerfStatsThreadSlot = 24,

      kMaxThreadSlot       = 25  // Size of the array of thread local slots in TThread
   };
}

#ifndef __CINT__
R__EXTERN void **(*gThreadTsd)(void*,Int_t);
#endif

#endif // ROOT_TThreadSlots
