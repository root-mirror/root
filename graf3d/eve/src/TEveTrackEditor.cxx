// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveTrackEditor.h"
#include "TEveTrack.h"
#include "TEveTrackPropagator.h"
#include "TEveTrackPropagatorEditor.h"
#include "TEveManager.h"

#include "TEveGValuators.h"

#include "TGedEditor.h"
#include "TAttMarkerEditor.h"
#include "TGLabel.h"
#include "TG3DLine.h"
#include "TGButton.h"
#include "TGDoubleSlider.h"

/** \class TEveTrackEditor
\ingroup TEve
Editor for TEveTrack class.
*/

ClassImp(TEveTrackEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveTrackEditor::TEveTrackEditor(const TGWindow *p, Int_t width, Int_t height,
                                 UInt_t options, Pixel_t back) :
   TGedFrame(p, width, height, options | kVerticalFrame, back),
   fM(0),
   fRSEditor(0)
{
   MakeTitle("TEveTrack");

   TGHorizontalFrame* f = new TGHorizontalFrame(this);

   fRSEditor =  new TGTextButton(f, "Edit Propagator");
   fRSEditor->Connect("Clicked()", "TEveTrackEditor", this, "DoEditPropagator()");
   f->AddFrame(fRSEditor, new TGLayoutHints(kLHintsLeft, 2, 1, 4, 4));

   AddFrame(f, new TGLayoutHints(kLHintsTop, 0,0,2,1));
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

void TEveTrackEditor::SetModel(TObject* obj)
{
   fM = dynamic_cast<TEveTrack*>(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for EditPropagator.

void TEveTrackEditor::DoEditPropagator()
{
   fGedEditor->SetModel(fGedEditor->GetPad(), fM->GetPropagator(), kButton1Down);
}

/** \class TEveTrackListEditor
\ingroup TEve
Editor for TEveTrackList class.
*/

ClassImp(TEveTrackListEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveTrackListEditor::TEveTrackListEditor(const TGWindow *p,
                                         Int_t width, Int_t height,
                                         UInt_t options, Pixel_t back) :
   TGedFrame(p, width, height, options | kVerticalFrame, back),

   fTC         (0),
   fPtRange    (0),
   fPRange     (0),
   fRSSubEditor(0)
{
   {
      TGHorizontalFrame* f = new TGHorizontalFrame(this);

      fRnrPoints = new TGCheckButton(f, "Draw Marker");
      f->AddFrame(fRnrPoints, new TGLayoutHints(kLHintsLeft, 2,1,0,0));
      fRnrPoints->Connect("Toggled(Bool_t)", "TEveTrackListEditor", this, "DoRnrPoints()");

      fRnrLine  = new TGCheckButton(f, "Draw Line");
      f->AddFrame(fRnrLine, new TGLayoutHints(kLHintsLeft, 1,2,0,0));
      fRnrLine->Connect("Toggled(Bool_t)", "TEveTrackListEditor", this, "DoRnrLine()");

      AddFrame(f, new TGLayoutHints(kLHintsTop, 0,0,2,1));
   }
   {  // --- Selectors
      Int_t labelW = 51;
      Int_t dbW    = 210;

      fPtRange = new TEveGDoubleValuator(this,"Pt rng:", 40, 0);
      fPtRange->SetNELength(6);
      fPtRange->SetLabelWidth(labelW);
      fPtRange->Build();
      fPtRange->GetSlider()->SetWidth(dbW);
      fPtRange->SetLimits(0, 10, TGNumberFormat::kNESRealTwo);
      fPtRange->Connect("ValueSet()",
                        "TEveTrackListEditor", this, "DoPtRange()");
      AddFrame(fPtRange, new TGLayoutHints(kLHintsTop, 1, 1, 4, 1));

      fPRange = new TEveGDoubleValuator(this,"P rng:", 40, 0);
      fPRange->SetNELength(6);
      fPRange->SetLabelWidth(labelW);
      fPRange->Build();
      fPRange->GetSlider()->SetWidth(dbW);
      fPRange->SetLimits(0, 100, TGNumberFormat::kNESRealTwo);
      fPRange->Connect("ValueSet()",
                       "TEveTrackListEditor", this, "DoPRange()");
      AddFrame(fPRange, new TGLayoutHints(kLHintsTop, 1, 1, 4, 1));
   }

   MakeTitle("RenderStyle");
   fRSSubEditor = new TEveTrackPropagatorSubEditor(this);
   fRSSubEditor->Connect("Changed()", "TEveTrackListEditor", this, "Update()");
   AddFrame(fRSSubEditor, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0,0,0,0));
   CreateRefsTab();
}

////////////////////////////////////////////////////////////////////////////////
/// Create tab for control of path-mark display.

void TEveTrackListEditor::CreateRefsTab()
{
   fRefs = CreateEditorTabSubFrame("Refs");

   TGCompositeFrame *title1 = new TGCompositeFrame(fRefs, 145, 10,
                                                   kHorizontalFrame |
                                                   kLHintsExpandX   |
                                                   kFixedWidth      |
                                                   kOwnBackground);
   title1->AddFrame(new TGLabel(title1, "PathMarks"),
                    new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   title1->AddFrame(new TGHorizontal3DLine(title1),
                    new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
   fRefs->AddFrame(title1, new TGLayoutHints(kLHintsTop, 0, 0, 2, 0));

   // path marks
   fRSSubEditor->CreateRefsContainer(fRefs);
   fRSSubEditor->fPMAtt->SetGedEditor((TGedEditor*)gEve->GetEditor());
   fRSSubEditor->fFVAtt->SetGedEditor((TGedEditor*)gEve->GetEditor());
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

void TEveTrackListEditor::SetModel(TObject* obj)
{
   fTC = dynamic_cast<TEveTrackList*>(obj);
   fRnrLine  ->SetState(fTC->GetRnrLine()   ? kButtonDown : kButtonUp);
   fRnrPoints->SetState(fTC->GetRnrPoints() ? kButtonDown : kButtonUp);

   Float_t llim;
   fPtRange->SetValues(fTC->fMinPt, fTC->fMaxPt);
   llim = fTC->fLimPt > 1 ? TMath::Log10(fTC->fLimPt) : 0;
   fPtRange->SetLimits(0, fTC->fLimPt, llim < 2 ? TGNumberFormat::kNESRealTwo : (llim < 3 ? TGNumberFormat::kNESRealOne : TGNumberFormat::kNESInteger));
   fPRange ->SetValues(fTC->fMinP, fTC->fMaxP);
   llim = fTC->fLimP > 1 ? TMath::Log10(fTC->fLimP) : 0;
   fPRange ->SetLimits(0, fTC->fLimP, llim < 2 ? TGNumberFormat::kNESRealTwo : (llim < 3 ? TGNumberFormat::kNESRealOne : TGNumberFormat::kNESInteger));

   fRSSubEditor->SetModel(fTC->GetPropagator());
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for RnrLine.

void TEveTrackListEditor::DoRnrLine()
{
   fTC->SetRnrLine(fRnrLine->IsOn());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for RnrPoints.

void TEveTrackListEditor::DoRnrPoints()
{
   fTC->SetRnrPoints(fRnrPoints->IsOn());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for PtRange.

void TEveTrackListEditor::DoPtRange()
{
   fTC->SelectByPt(fPtRange->GetMin(), fPtRange->GetMax());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for PRange.

void TEveTrackListEditor::DoPRange()
{
   fTC->SelectByP(fPRange->GetMin(), fPRange->GetMax());
   Update();
}

