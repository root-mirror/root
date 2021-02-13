// @(#)root/base:$Id$
// Author: Matevz Tadel   16/08/2006

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TObjectSpy.h"
#include "TSeqCollection.h"
#include "TROOT.h"
#include "TVirtualMutex.h"

/** \class TObjectRefSpy
    \class TObjectSpy
\ingroup Base

Monitors objects for deletion and reflects the deletion by reverting
the internal pointer to zero. When this pointer is zero we know the
object has been deleted. This avoids the unsafe TestBit(kNotDeleted)
hack. The spied object must have the kMustCleanup bit set otherwise
you will get an error.
*/

ClassImp(TObjectSpy);
ClassImp(TObjectRefSpy);

////////////////////////////////////////////////////////////////////////////////
/// Register the object that must be spied. The object must have the
/// kMustCleanup bit set. If the object has been deleted during a
/// RecusiveRemove() operation, GetObject() will return 0.

TObjectSpy::TObjectSpy(TObject *obj, Bool_t fixMustCleanupBit) :
   TObject(), fObj(obj), fResetMustCleanupBit(kFALSE)
{
   {
      R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfCleanups()->Add(this);
   }
   if (fObj && !fObj->TestBit(kMustCleanup)) {
      if (fixMustCleanupBit) {
         fResetMustCleanupBit = kTRUE;
         fObj->SetBit(kMustCleanup, kTRUE);
      } else {
         Error("TObjectSpy", "spied object must have the kMustCleanup bit set");
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup.

TObjectSpy::~TObjectSpy()
{
   if (fObj && fResetMustCleanupBit)
      fObj->SetBit(kMustCleanup, kFALSE);
   R__LOCKGUARD(gROOTMutex);
   gROOT->GetListOfCleanups()->Remove(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the object pointer to zero if the object is deleted in the
/// RecursiveRemove() operation.

void TObjectSpy::RecursiveRemove(TObject *obj)
{
   if (obj == fObj) {
      fObj = nullptr;
      fResetMustCleanupBit = kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set obj as the spy target.

void TObjectSpy::SetObject(TObject *obj, Bool_t fixMustCleanupBit)
{
   if (fObj && fResetMustCleanupBit)
      fObj->SetBit(kMustCleanup, kFALSE);
   fResetMustCleanupBit = kFALSE;

   fObj = obj;

   if (fObj && !fObj->TestBit(kMustCleanup)) {
      if (fixMustCleanupBit) {
         fResetMustCleanupBit = kTRUE;
         fObj->SetBit(kMustCleanup, kTRUE);
      } else {
         Error("TObjectSpy", "spied object must have the kMustCleanup bit set");
      }
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Register the object that must be spied. The object must have the
/// kMustCleanup bit set. If the object has been deleted during a
/// RecusiveRemove() operation, GetObject() will return 0.

TObjectRefSpy::TObjectRefSpy(TObject *&obj, Bool_t fixMustCleanupBit) :
   fObj(obj), fResetMustCleanupBit(kFALSE)
{
   {
      R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfCleanups()->Add(this);
   }
   if (fObj && !fObj->TestBit(kMustCleanup)) {
      if (fixMustCleanupBit) {
         fResetMustCleanupBit = kTRUE;
         fObj->SetBit(kMustCleanup, kTRUE);
      } else {
         Error("TObjectSpy", "spied object must have the kMustCleanup bit set");
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup.

TObjectRefSpy::~TObjectRefSpy()
{
   if (fObj && fResetMustCleanupBit)
      fObj->SetBit(kMustCleanup, kFALSE);
   R__LOCKGUARD(gROOTMutex);
   gROOT->GetListOfCleanups()->Remove(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the object pointer to zero if the object is deleted in the
/// RecursiveRemove() operation.

void TObjectRefSpy::RecursiveRemove(TObject *obj)
{
   if (obj == fObj) {
      fObj = nullptr;
      fResetMustCleanupBit = kFALSE;
   }
}
