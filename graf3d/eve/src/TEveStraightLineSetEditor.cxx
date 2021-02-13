// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveStraightLineSetEditor.h"
#include "TEveStraightLineSet.h"

#include "TGButton.h"

/** \class TEveStraightLineSetEditor
\ingroup TEve
Editor for TEveStraightLineSet class.
*/

ClassImp(TEveStraightLineSetEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveStraightLineSetEditor::TEveStraightLineSetEditor(const TGWindow *p, Int_t width, Int_t height,
                                                     UInt_t options, Pixel_t back) :
   TGedFrame(p, width, height, options | kVerticalFrame, back),
   fM(0)
                                                    // Initialize widget pointers to 0
{
   MakeTitle("TEveStraightLineSet");

   TGHorizontalFrame* frame = new TGHorizontalFrame(this);

   fRnrMarkers = new TGCheckButton(frame, "RnrMarkers");
   frame->AddFrame(fRnrMarkers, new TGLayoutHints(kLHintsLeft, 1, 2, 1, 1));
   fRnrMarkers->Connect
      ("Toggled(Bool_t)",
       "TEveStraightLineSetEditor", this, "DoRnrMarkers()");

   fRnrLines = new TGCheckButton(frame, "RnrLines");
   frame->AddFrame(fRnrLines, new TGLayoutHints(kLHintsLeft, 2, 1, 1, 1));
   fRnrLines->Connect
      ("Toggled(Bool_t)",
       "TEveStraightLineSetEditor", this, "DoRnrLines()");

   AddFrame(frame, new TGLayoutHints(kLHintsTop, 0, 0, 0, 0));
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

void TEveStraightLineSetEditor::SetModel(TObject* obj)
{
   fM = dynamic_cast<TEveStraightLineSet*>(obj);

   // Set values of widgets
   fRnrMarkers->SetState(fM->GetRnrMarkers() ? kButtonDown : kButtonUp);
   fRnrLines->SetState(fM->GetRnrLines() ? kButtonDown : kButtonUp);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for RnrMarkers.

void TEveStraightLineSetEditor::DoRnrMarkers()
{
   fM->SetRnrMarkers(fRnrMarkers->IsOn());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for RnrLines.

void TEveStraightLineSetEditor::DoRnrLines()
{
   fM->SetRnrLines(fRnrLines->IsOn());
   Update();
}
