// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveArrowEditor.h"
#include "TEveArrow.h"
#include "TGLabel.h"
#include "TEveGValuators.h"

/** \class TEveArrowEditor
\ingroup TEve
GUI editor for TEveArrow.
*/

ClassImp(TEveArrowEditor);

////////////////////////////////////////////////////////////////////////////////

TEveArrowEditor::TEveArrowEditor(const TGWindow *p, Int_t width, Int_t height,
                                 UInt_t options, Pixel_t back) :
   TGedFrame(p, width, height, options | kVerticalFrame, back),
   fM(0),

   fTubeR(0),
   fConeR(0),
   fConeL(0),

   fOrigin(0),
   fVector(0)
{
   // Constructor.

   MakeTitle("TEveArrow");
   Int_t labelW = 51;

   fTubeR = new TEveGValuator(this, "TubeR:", 90, 0);
   fTubeR->SetLabelWidth(labelW);
   fTubeR->SetNELength(6);
   fTubeR->Build();
   fTubeR->SetLimits(0.001, 1, 101, TGNumberFormat::kNESRealThree);
   fTubeR->Connect("ValueSet(Double_t)", "TEveArrowEditor", this, "DoTubeR()");
   AddFrame(fTubeR, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

   fConeR = new TEveGValuator(this, "ConeR:", 90, 0);
   fConeR->SetLabelWidth(labelW);
   fConeR->SetNELength(6);
   fConeR->Build();
   fConeR->SetLimits(0.001, 1, 101, TGNumberFormat::kNESRealThree);
   fConeR->Connect("ValueSet(Double_t)", "TEveArrowEditor", this, "DoConeR()");
   AddFrame(fConeR, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

   fConeL = new TEveGValuator(this, "ConeL:", 90, 0);
   fConeL->SetLabelWidth(labelW);
   fConeL->SetNELength(6);
   fConeL->Build();
   fConeL->SetLimits(0.001, 1, 101, TGNumberFormat::kNESRealThree);
   fConeL->Connect("ValueSet(Double_t)", "TEveArrowEditor", this, "DoConeL()");
   AddFrame(fConeL, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

   TGFont *font = gClient->GetFont("-adobe-helvetica-bold-r-*-*-12-*-*-*-*-*-iso8859-1");
   {
      TGHorizontalFrame* hfp = new TGHorizontalFrame(this);
      TGLabel* labp = new TGLabel(hfp, "Origin");
      labp->SetTextFont(font);
      hfp->AddFrame(labp);
      AddFrame(hfp,  new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0,0,2,0));
      fOrigin = new TEveGTriVecValuator(this, "Origin", 160, 20);
      fOrigin->SetNELength(6);
      fOrigin->Build(kFALSE, "", "", "");
      fOrigin->SetLimits(-1e5, 1e5, TGNumberFormat::kNESRealThree);
      AddFrame(fOrigin, new TGLayoutHints(kLHintsTop , 0,0,0,0));
   }

   {
      TGHorizontalFrame* hfp = new TGHorizontalFrame(this);
      TGLabel* labp = new TGLabel(hfp, "Vector");
      labp->SetTextFont(font);
      hfp->AddFrame(labp);
      AddFrame(hfp,  new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0,0,2,0));
      fVector = new TEveGTriVecValuator(this, "Vector", 160, 20);
      fVector->SetNELength(6);
      fVector->Build(kFALSE, "", "", "");
      fVector->SetLimits(-1e5, 1e5, TGNumberFormat::kNESRealThree);
      AddFrame(fVector, new TGLayoutHints(kLHintsTop , 0,0,0,0));
   }

   fVector->Connect("ValueSet()", "TEveArrowEditor", this, "DoVertex()");
   fOrigin->Connect("ValueSet()", "TEveArrowEditor", this, "DoVertex()");
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

void TEveArrowEditor::SetModel(TObject* obj)
{
   fM = dynamic_cast<TEveArrow*>(obj);

   fTubeR->SetValue(fM->fTubeR);
   fTubeR->SetValue(fM->fTubeR);
   fConeL->SetValue(fM->fConeL);

   fOrigin->SetValues(&fM->fOrigin[0]);
   fVector->SetValues(&fM->fVector[0]);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for TubeR.

void TEveArrowEditor::DoTubeR()
{
   fM->SetTubeR(fTubeR->GetValue());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for ConeR.

void TEveArrowEditor::DoConeR()
{
   fM->SetConeR(fConeR->GetValue());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for ConeL.

void TEveArrowEditor::DoConeL()
{
   fM->SetConeL(fConeL->GetValue());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for origin and vector.

void TEveArrowEditor::DoVertex()
{
   Double_t v[3];
   fOrigin->GetValues(v);
   fM->fOrigin.Set(v);
   fVector->GetValues(v);
   fM->fVector.Set(v);
   fM->StampGeom();
   Update();
}
