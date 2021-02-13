// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveLineEditor.h"
#include "TEveLine.h"

#include "TGButton.h"

/** \class TEveLineEditor
\ingroup TEve
Editor for TEveLine class.
*/

ClassImp(TEveLineEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveLineEditor::TEveLineEditor(const TGWindow *p, Int_t width, Int_t height,
                               UInt_t options, Pixel_t back) :
   TGedFrame  (p, width, height, options | kVerticalFrame, back),
   fM         (0),
   fRnrLine   (0),
   fRnrPoints (0),
   fSmooth    (0)
{
   fPriority = 20;
   {
      TGHorizontalFrame* f = new TGHorizontalFrame(this);

      fRnrPoints = new TGCheckButton(f, "Draw Marker");
      f->AddFrame(fRnrPoints, new TGLayoutHints(kLHintsLeft, 2,1,0,0));
      fRnrPoints->Connect("Toggled(Bool_t)", "TEveLineEditor", this, "DoRnrPoints()");

      fRnrLine  = new TGCheckButton(f, "Draw Line");
      f->AddFrame(fRnrLine, new TGLayoutHints(kLHintsLeft, 1,2,0,0));
      fRnrLine->Connect("Toggled(Bool_t)", "TEveLineEditor", this, "DoRnrLine()");

      AddFrame(f, new TGLayoutHints(kLHintsTop, 0,0,2,1));
   }
   fSmooth  = new TGCheckButton(this, "Smooth line");
   AddFrame(fSmooth, new TGLayoutHints(kLHintsNormal, 2,1,2,1));
   fSmooth->Connect("Toggled(Bool_t)", "TEveLineEditor", this, "DoSmooth()");
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

void TEveLineEditor::SetModel(TObject* obj)
{
   fM = dynamic_cast<TEveLine*>(obj);

   fRnrLine  ->SetState(fM->fRnrLine  ? kButtonDown : kButtonUp);
   fRnrPoints->SetState(fM->fRnrPoints ? kButtonDown : kButtonUp);
   fSmooth->SetState(fM->fSmooth ? kButtonDown : kButtonUp);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for RnrLine.

void TEveLineEditor::DoRnrLine()
{
   fM->SetRnrLine(fRnrLine->IsOn());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for RnrPoints.

void TEveLineEditor::DoRnrPoints()
{
   fM->SetRnrPoints(fRnrPoints->IsOn());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for anti-alias.

void TEveLineEditor::DoSmooth()
{
   fM->SetSmooth(fSmooth->IsOn());
   Update();
}
