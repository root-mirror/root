// @(#)root/graf:$Id$
// Author: Rene Brun   05/03/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <cstdio>

#include "TVirtualPad.h"
#include "TClass.h"
#include "TLink.h"
#include "snprintf.h"

ClassImp(TLink);

/** \class TLink
\ingroup BasicGraphics

Special TText object used to show hyperlinks.

In the example below created by TObject::Inspect, TLinks are used
to show pointers to other objects.
Clicking on one link, inspect the corresponding object.
*/

////////////////////////////////////////////////////////////////////////////////
/// Link default constructor.

TLink::TLink() : TText()
{
   fLink  = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor to define a link object.
///
/// pointer points to any kind of object.

TLink::TLink(Double_t x, Double_t y, void *pointer)
           : TText(x, y, "")
{
   fLink  = pointer;
   static char line[16];
   snprintf(line,16,"->%lx ", (Long_t)pointer);
   SetTitle(line);
}

////////////////////////////////////////////////////////////////////////////////
/// Link default destructor.

TLink::~TLink()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Execute action corresponding to one event.
///
///  This member function is called when a link is clicked with the locator
///
///  If mouse is clicked on a link text, the object pointed by the link
///  is Inspected

void TLink::ExecuteEvent(Int_t event, Int_t, Int_t)
{
   if (!gPad) return;

   if (event == kMouseMotion)
      gPad->SetCursor(kHand);

   if (event != kButton1Up) return;

   if (TestBit(kIsStarStar)) return;
   TObject *idcur = (TObject*)fLink;
   if (!idcur) return;
   TClass *cl = TClass::GetClass(GetName());
   if (!cl) return;

   // check if link points to a TObject
   TClass *c1 = (TClass*)cl->GetBaseClass("TObject");
   if (!c1) return;

   idcur->Inspect();
}
