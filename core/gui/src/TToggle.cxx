// @(#)root/meta:$Id$
// Author: Piotr Golonka   30/07/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TToggle

This class defines toggling facility for both - object's method or
variables.
Assume that user provides an object with a two-state field , and
methods to Get/Set value of this field. This object enables to switch
values via this method when the only thing you know about the field
is the name of the method (or method itself) which sets the field.
This facility is required in context Pop-Up menu, when the only
information about how to toggle a field is a name of methhod which
sets it.
This class may be also used for toggling an integer variable,
which may be important while building universal objects...
When user provides a "set-method" of name SetXXX this object tries
automaticaly find a matching "get-method" by lookin for a method
with name GetXXX, IsXXX or HasXXX for given object.
*/

#include "TMethod.h"
#include "TMethodCall.h"
#include "TToggle.h"
#include "TDataMember.h"
#include "snprintf.h"

ClassImp(TToggle);

////////////////////////////////////////////////////////////////////////////////
/// TToggle default constructor. You have to initialize it before using
/// by making a call to SetToggledVariable() or SetToggledObject().

TToggle::TToggle()
{
   fState       =  kFALSE;
   fValue       = -1;
   fOnValue     =  1;
   fOffValue    =  0;
   fInitialized =  0;
   fObject      =  0;
   fGetter      =  0;
   fSetter      =  0;
   fTglVariable =  0;
}

////////////////////////////////////////////////////////////////////////////////
/// Initializes object for use with a variable - you pass it via reference
/// so it will be modified by Toggle.

void TToggle::SetToggledVariable(Int_t &var)
{
   fTglVariable=&var;
   fValue=var;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the state of Toggle according to its current value and
/// fOnValue, returns true if they match.

Bool_t TToggle::GetState()
{
   if (fInitialized)
      if (fGetter) fGetter->Execute(fObject, fValue);
   return (fState = (fValue == fOnValue));
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the value of toggle to fOnValue or fOffValue according to passed
/// argument.

void TToggle::SetState(Bool_t state)
{
   if (fInitialized) {
      char stringon[20];
      char stringoff[20];
      snprintf(stringon,sizeof(stringon),"%li",fOnValue);
      snprintf(stringoff,sizeof(stringoff),"%li",fOffValue);

      fSetter->Execute(fObject, state ? stringon:stringoff);
      fState=state;
      fValue= ( state ? fOnValue : fOffValue);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the value of toggle and modifies its state according to whether
/// the value is equal to fOnValue.

void TToggle::SetValue(Long_t val)
{
   if (fInitialized) {
      char stringval[20];
      snprintf(stringval,sizeof(stringval),"%li",val);
      fSetter->Execute(fObject, stringval);
      fState=(val==fOnValue);
      fValue= val;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Toggles the Values and State of this object and connected data!

void TToggle::Toggle()
{
   if (fInitialized){
      if (fTglVariable){
         *fTglVariable = !(*fTglVariable);
         fValue=(*fTglVariable);
         fState=*fTglVariable;
      }
      if (fGetter && fSetter){
         fGetter->Execute(fObject,fValue);
         fValue=( (fValue==fOnValue) ? fOffValue:fOnValue);
         fState=(!(fValue!=fOnValue));
         char stringon[20];
         snprintf(stringon,sizeof(stringon),"%li",fValue);
         fSetter->Execute(fObject, stringon);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Initializes it to toggle an object's datamember using this object's
/// method.

void TToggle::SetToggledObject(TObject *obj, TMethod *anymethod)
{
   fObject = obj;
   TDataMember *m = anymethod->FindDataMember();
   if (!m) {
      // try to see if the TMethod has a getter associated via the *GETTER=
      // comment string
      if (anymethod->GetterMethod()) {
         fGetter = anymethod->GetterMethod();
         fSetter = anymethod->SetterMethod();
         fInitialized = 1;
      } else
         Error("SetToggledObject", "cannot determine getter method for %s", anymethod->GetName());
   } else {
      fGetter = m->GetterMethod(obj->IsA());
      fSetter = m->SetterMethod(obj->IsA());
      fInitialized = 1;
   }
}
