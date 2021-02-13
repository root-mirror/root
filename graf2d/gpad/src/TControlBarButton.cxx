// @(#)root/gpad:$Id$
// Author: Nenad Buncic   20/02/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TControlBarButton
\ingroup gpad
This class defines the control bar buttons

Created by the TControlBar. Not for general consumption.
*/

#include "TControlBarButton.h"
#include "TVirtualPad.h"
#include "TError.h"
#include "TApplication.h"

#include <cctype>

const char *kBStr = "BUTTON";
const char *kDStr = "DRAWNBUTTON";
const char *kSStr = "SEPARATOR";


ClassImp(TControlBarButton);

////////////////////////////////////////////////////////////////////////////////
/// Default control bar button ctor.

TControlBarButton::TControlBarButton() : TNamed()
{
   fType   = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create control bar button.

TControlBarButton::TControlBarButton(const char *label, const char *action,
                                     const char *hint, const char *type)
   : TNamed(label, hint)
{
   SetType(type);
   SetAction(action);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute control bar button command.

void TControlBarButton::Action()
{
   if (!fAction.IsNull()) {

      gApplication->ProcessLine(fAction.Data());

      if (gPad) gPad->Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set action to be executed by this button.

void TControlBarButton::SetAction(const char *action)
{
   if (action) {
      char *s = Strip(action);
      fAction = s;
      delete [] s;
   } else
      Error("SetAction", "action missing");
}


////////////////////////////////////////////////////////////////////////////////
/// Set button type. Type can be either "button", "drawnbutton" or
/// "separator". String is case insensitive. Default is "button".

void TControlBarButton::SetType(const char *type)
{
   fType = kButton;

   if (type && *type) {
      if (!strcasecmp(type, kBStr))
         fType = kButton;
      else if (!strcasecmp(type, kDStr))
         fType = kDrawnButton;
      else if (!strcasecmp(type, kSStr))
         fType = kSeparator;
      else
         Error("SetType", "unknown type '%s' !\n\t(choice of: %s, %s, %s)",
               type, kBStr, kDStr, kSStr);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set button type. Type can be either kButton, kDrawnButton or kSeparator.
/// Default is kButton.

void TControlBarButton::SetType(Int_t type)
{
   switch (type) {

      case kButton:
      case kDrawnButton:
      case kSeparator:
         fType = type;
         break;

      default:
         fType = kButton;
         Error("SetType", "unknown type: %d !\n\t(choice of: %d, %d, %d)",
               type, kButton, kDrawnButton, kSeparator);
   }
}
