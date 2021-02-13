// @(#)root/base:$Id$
// Author: Rene Brun   05/12/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TVirtualPad.h"
#include "TBuffer.h"
#include "TThreadSlots.h"

/** \class TVirtualPad
\ingroup Base

TVirtualPad is an abstract base class for the Pad and Canvas classes.
*/

Int_t (*gThreadXAR)(const char *xact, Int_t nb, void **ar, Int_t *iret) = 0;

////////////////////////////////////////////////////////////////////////////////
/// Return the current pad for the current thread.

TVirtualPad *&TVirtualPad::Pad()
{
   static TVirtualPad *currentPad = 0;
   if (!gThreadTsd)
      return currentPad;
   else
      return *(TVirtualPad**)(*gThreadTsd)(&currentPad,ROOT::kPadThreadSlot);
}

ClassImp(TVirtualPad);

////////////////////////////////////////////////////////////////////////////////
/// VirtualPad default constructor

TVirtualPad::TVirtualPad() : TAttPad()
{
   fResizing = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// VirtualPad constructor

TVirtualPad::TVirtualPad(const char *, const char *, Double_t,
           Double_t, Double_t, Double_t, Color_t color, Short_t , Short_t)
          : TAttPad()
{
   fResizing = kFALSE;

   SetFillColor(color);
   SetFillStyle(1001);
}

////////////////////////////////////////////////////////////////////////////////
/// VirtualPad destructor

TVirtualPad::~TVirtualPad()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TVirtualPad.

void TVirtualPad::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         R__b.ReadClassBuffer(TVirtualPad::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TObject::Streamer(R__b);
      TAttLine::Streamer(R__b);
      TAttFill::Streamer(R__b);
      TAttPad::Streamer(R__b);
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TVirtualPad::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Should always return false unless you have non-standard picking.

Bool_t TVirtualPad::PadInSelectionMode() const
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Should always return false, unless you can highlight selected object in pad.

Bool_t TVirtualPad::PadInHighlightMode() const
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Does nothing, unless you implement your own picking.
/// When complex object containing sub-objects (which can be picked)
/// is painted in a pad, this "top-level" object is pushed into
/// the selectables stack.

void TVirtualPad::PushTopLevelSelectable(TObject * /*object*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Does nothing, unless you implement your own picking.
/// "Complete" object, or part of complex object, which
/// can be picked.

void TVirtualPad::PushSelectableObject(TObject * /*object*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Does nothing, unless you implement your own picking.
/// Remove top level selectable and all its' children.

void TVirtualPad::PopTopLevelSelectable()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Scope-guards ctor, pushe the object on stack.

TPickerStackGuard::TPickerStackGuard(TObject *obj)
{
   gPad->PushTopLevelSelectable(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Guard does out of scope, pop object from stack.

TPickerStackGuard::~TPickerStackGuard()
{
   gPad->PopTopLevelSelectable();
}
