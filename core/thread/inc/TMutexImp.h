// @(#)root/thread:$Id$
// Author: Fons Rademakers   01/07/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMutexImp
#define ROOT_TMutexImp


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMutexImp                                                            //
//                                                                      //
// This class provides an abstract interface to the OS dependent mutex  //
// classes (TPosixMutex and TWin32Mutex).                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"
#include "TObject.h"

#include "TVirtualMutex.h"


class TMutexImp : public TObject {

public:
   TMutexImp() { }
   virtual ~TMutexImp() { }

   virtual Int_t  Lock() = 0;
   virtual Int_t  TryLock() = 0;
   virtual Int_t  UnLock() = 0;

   virtual std::unique_ptr<TVirtualMutex::State> Reset() = 0;
   virtual void Restore(std::unique_ptr<TVirtualMutex::State>&&) = 0;

   ClassDef(TMutexImp,0)  // Mutex lock implementation ABC
};

#endif
