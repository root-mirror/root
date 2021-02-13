// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveShapeEditor.h"
#include "TEveShape.h"

#include "TColor.h"

#include "TGLabel.h"
#include "TGButton.h"
#include "TGNumberEntry.h"
#include "TGColorSelect.h"


/** \class TEveShapeEditor
\ingroup TEve
GUI editor for TEveShape.
*/

ClassImp(TEveShapeEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveShapeEditor::TEveShapeEditor(const TGWindow *p, Int_t width, Int_t height,
             UInt_t options, Pixel_t back) :
   TGedFrame(p, width, height, options | kVerticalFrame, back),
   fM(0),
   fLineWidth(0),
   fLineColor(0),
   fDrawFrame(0),
   fHighlightFrame(0)
{
   MakeTitle("TEveShape");

   {
      TGCompositeFrame *f = new TGCompositeFrame(this, 80, 20, kHorizontalFrame);

      TGLabel *l = new TGLabel(f, "LineColor:");
      f->AddFrame(l, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 0));
      fLineColor = new TGColorSelect(f, 0, -1);
      fLineColor->Connect("ColorSelected(Pixel_t)", "TEveShapeEditor", this, "DoLineColor(Pixel_t)");
      f->AddFrame(fLineColor, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));

      fLineWidth = new TGNumberEntry(f, 0., 4, -1,
                                     TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAPositive,
                                     TGNumberFormat::kNELLimitMinMax, 0.1, 20.0);
      fLineWidth->GetNumberEntry()->SetToolTipText("Line width of outline.");
      fLineWidth->Connect("ValueSet(Long_t)", "TEveShapeEditor", this, "DoLineWidth()");
      f->AddFrame(fLineWidth, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));

      AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 0, 0));
   }
   {
      TGHorizontalFrame* f = new TGHorizontalFrame(this);

      fDrawFrame  = new TGCheckButton(f, "Draw Frame");
      f->AddFrame(fDrawFrame, new TGLayoutHints(kLHintsLeft, 1,2,0,0));
      fDrawFrame->Connect("Toggled(Bool_t)", "TEveShapeEditor", this, "DoDrawFrame()");

      fHighlightFrame = new TGCheckButton(f, "Highlight Frame");
      f->AddFrame(fHighlightFrame, new TGLayoutHints(kLHintsLeft, 2,1,0,0));
      fHighlightFrame->Connect("Toggled(Bool_t)", "TEveShapeEditor", this, "DoHighlightFrame()");

      AddFrame(f, new TGLayoutHints(kLHintsTop, 0,0,2,1));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

void TEveShapeEditor::SetModel(TObject* obj)
{
   fM = dynamic_cast<TEveShape*>(obj);

   fLineWidth->SetNumber(fM->fLineWidth);
   fLineColor->SetColor(TColor::Number2Pixel(fM->fLineColor), kFALSE);
   fDrawFrame     ->SetState(fM->fDrawFrame      ? kButtonDown : kButtonUp);
   fHighlightFrame->SetState(fM->fHighlightFrame ? kButtonDown : kButtonUp);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting line with of polygon outline.

void TEveShapeEditor::DoLineWidth()
{
   fM->SetLineWidth(fLineWidth->GetNumber());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting line color of polygon outline.

void TEveShapeEditor::DoLineColor(Pixel_t pixel)
{
   fM->SetLineColor(TColor::GetColor(pixel));
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for DrawFrame.

void TEveShapeEditor::DoDrawFrame()
{
   fM->SetDrawFrame(fDrawFrame->IsOn());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for HighlightFrame.

void TEveShapeEditor::DoHighlightFrame()
{
   fM->SetHighlightFrame(fHighlightFrame->IsOn());
   Update();
}
