// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveGridStepperEditor.h"
#include "TEveGridStepper.h"
#include "TEveGValuators.h"

/** \class TEveGridStepperSubEditor
\ingroup TEve
Sub-editor for TEveGridStepper class.
*/

ClassImp(TEveGridStepperSubEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveGridStepperSubEditor::TEveGridStepperSubEditor(const TGWindow *p) :
   TGVerticalFrame(p),
   fM (0),
   fNx(0), fNy(0), fNz(0),
   fDx(0), fDy(0), fDz(0)
{
   Int_t labelW = 15;

   TGHorizontalFrame* hf = new TGHorizontalFrame(this);

   {
      TGGroupFrame* f = new TGGroupFrame(hf, "NumRows", kVerticalFrame);
      f->SetWidth(30);
      hf->AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 0));

      fNx = new TEveGValuator(f,"X:", 200, 0);
      fNx->SetNELength(3);
      fNx->SetLabelWidth(labelW);
      fNx->SetShowSlider(kFALSE);
      fNx->Build();
      fNx->SetLimits(1, 15);
      fNx->Connect("ValueSet(Double_t)",
                   "TEveGridStepperSubEditor", this, "DoNs()");
      f->AddFrame(fNx, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

      fNy = new TEveGValuator(f,"Y:", 200, 0);
      fNy->SetNELength(3);
      fNy->SetLabelWidth(labelW);
      fNy->SetShowSlider(kFALSE);
      fNy->Build();
      fNy->SetLimits(1, 15);
      fNy->Connect("ValueSet(Double_t)",
                   "TEveGridStepperSubEditor", this, "DoNs()");
      f->AddFrame(fNy, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

      fNz = new TEveGValuator(f,"Z:", 200, 0);
      fNz->SetNELength(3);
      fNz->SetLabelWidth(labelW);
      fNz->SetShowSlider(kFALSE);
      fNz->Build();
      fNz->SetLimits(1, 15);
      fNz->Connect("ValueSet(Double_t)",
                   "TEveGridStepperSubEditor", this, "DoNs()");
      f->AddFrame(fNz, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

      //AddFrame(f, new TGLayoutHints(kLHintsExpandX, 2, 0, 0, 0));
   }
   {
      TGGroupFrame* f = new TGGroupFrame(hf, "Step", kVerticalFrame);
      f->SetWidth(130);
      hf->AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 0));

      fDx = new TEveGValuator(f,"X:", 200, 0);
      fDx->SetNELength(5);
      fDx->SetLabelWidth(labelW);
      fDx->SetShowSlider(kFALSE);
      fDx->Build();
      fDx->SetLimits(0.1, 100, 101, TGNumberFormat::kNESRealOne);
      fDx->Connect("ValueSet(Double_t)",
                   "TEveGridStepperSubEditor", this, "DoDs()");
      f->AddFrame(fDx, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

      fDy = new TEveGValuator(f,"Y:", 200, 0);
      fDy->SetNELength(5);
      fDy->SetLabelWidth(labelW);
      fDy->SetShowSlider(kFALSE);
      fDy->Build();
      fDy->SetLimits(0.1, 100, 101, TGNumberFormat::kNESRealOne);
      fDy->Connect("ValueSet(Double_t)",
                   "TEveGridStepperSubEditor", this, "DoDs()");
      f->AddFrame(fDy, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

      fDz = new TEveGValuator(f,"Z:", 200, 0);
      fDz->SetNELength(5);
      fDz->SetLabelWidth(labelW);
      fDz->SetShowSlider(kFALSE);
      fDz->Build();
      fDz->SetLimits(0.1, 100, 101, TGNumberFormat::kNESRealOne);
      fDz->Connect("ValueSet(Double_t)",
                   "TEveGridStepperSubEditor", this, "DoDs()");
      f->AddFrame(fDz, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

      //AddFrame(f, new TGLayoutHints(kLHintsExpandX, 2, 0, 0, 0));
   }
   AddFrame(hf, new TGLayoutHints(kLHintsExpandX, 2, 0, 0, 0));
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

void TEveGridStepperSubEditor::SetModel(TEveGridStepper* m)
{
   fM = m;

   fNx->SetValue(fM->fNx);
   fNy->SetValue(fM->fNy);
   fNz->SetValue(fM->fNz);

   fDx->SetValue(fM->fDx);
   fDy->SetValue(fM->fDy);
   fDz->SetValue(fM->fDz);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit Changed signal.

void TEveGridStepperSubEditor::Changed()
{
   Emit("Changed()");
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for changing fN's.

void TEveGridStepperSubEditor::DoNs()
{
   fM->SetNs((Int_t)fNx->GetValue(), (Int_t)fNy->GetValue(), (Int_t)fNz->GetValue());
   Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for changing fD's.

void TEveGridStepperSubEditor::DoDs()
{
   fM->SetDs(fDx->GetValue(), fDy->GetValue(), fDz->GetValue());
   Changed();
}

/** \class TEveGridStepperEditor
\ingroup TEve
Editor for TEveGridStepper class.
*/

ClassImp(TEveGridStepperEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveGridStepperEditor::TEveGridStepperEditor(const TGWindow *p, Int_t width, Int_t height,
                                             UInt_t options, Pixel_t back) :
   TGedFrame(p, width, height, options | kVerticalFrame, back),
   fM  (0),
   fSE (0)
{
   MakeTitle("TEveGridStepper");

   fSE = new TEveGridStepperSubEditor(this);
   AddFrame(fSE, new TGLayoutHints(kLHintsTop, 2, 0, 2, 2));
   fSE->Connect("Changed()", "TEveGridStepperEditor", this, "Update()");
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

void TEveGridStepperEditor::SetModel(TObject* obj)
{
   fM = dynamic_cast<TEveGridStepper*>(obj);
   fSE->SetModel(fM);
}
