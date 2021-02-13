// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveProjectionManagerEditor.h"
#include "TEveProjectionManager.h"
#include "TEveGValuators.h"

#include "TGComboBox.h"
#include "TGLabel.h"

/** \class TEveProjectionManagerEditor
\ingroup TEve
GUI editor for TEveProjectionManager class.
*/

ClassImp(TEveProjectionManagerEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveProjectionManagerEditor::TEveProjectionManagerEditor(const TGWindow *p,
                                                         Int_t width, Int_t height,
                                                         UInt_t options, Pixel_t back) :
   TGedFrame(p, width, height, options | kVerticalFrame, back),
   fM(0),

   fType(0),
   fDistortion(0),
   fFixR(0), fFixZ(0),
   fPastFixRFac(0), fPastFixZFac(0),
   fCurrentDepth(0),
   fMaxTrackStep(0),

   fCenterX(0),
   fCenterY(0),
   fCenterZ(0)
{
   MakeTitle("TEveProjection");
   {
      TGHorizontalFrame* f = new TGHorizontalFrame(this);
      TGLabel* lab = new TGLabel(f, "Type");
      f->AddFrame(lab, new TGLayoutHints(kLHintsLeft|kLHintsBottom, 1, 31, 1, 2));
      fType = new TGComboBox(f);
      fType->AddEntry("RPhi", TEveProjection::kPT_RPhi);
      fType->AddEntry("RhoZ", TEveProjection::kPT_RhoZ);
      fType->AddEntry("3D",   TEveProjection::kPT_3D);
      TGListBox* lb = fType->GetListBox();
      lb->Resize(lb->GetWidth(), 2*18);
      fType->Resize(80, 20);
      fType->Connect("Selected(Int_t)", "TEveProjectionManagerEditor",
                     this, "DoType(Int_t)");
      f->AddFrame(fType, new TGLayoutHints(kLHintsTop, 1, 1, 2, 4));
      AddFrame(f);
   }

   Int_t nel = 6;
   Int_t labelW = 60;
   fDistortion = new TEveGValuator(this, "Distortion:", 90, 0);
   fDistortion->SetNELength(nel);
   fDistortion->SetLabelWidth(labelW);
   fDistortion->Build();
   fDistortion->SetLimits(0, 50, 101, TGNumberFormat::kNESRealTwo);
   fDistortion->Connect("ValueSet(Double_t)", "TEveProjectionManagerEditor",
                        this, "DoDistortion()");
   AddFrame(fDistortion, new TGLayoutHints(kLHintsTop, 1, 1, 1, 0));


   fFixR = new TEveGValuator(this, "FixedR:", 90, 0);
   fFixR->SetNELength(nel);
   fFixR->SetLabelWidth(labelW);
   fFixR->Build();
   fFixR->SetLimits(0, 1000, 101, TGNumberFormat::kNESRealOne);
   fFixR->SetToolTip("Radius after which scale is kept constant.");
   fFixR->Connect("ValueSet(Double_t)", "TEveProjectionManagerEditor",
                         this, "DoFixR()");
   AddFrame(fFixR, new TGLayoutHints(kLHintsTop, 1, 1, 1, 0));

   fFixZ = new TEveGValuator(this, "FixedZ:", 90, 0);
   fFixZ->SetNELength(nel);
   fFixZ->SetLabelWidth(labelW);
   fFixZ->Build();
   fFixZ->SetLimits(0, 1000, 101, TGNumberFormat::kNESRealOne);
   fFixZ->SetToolTip("Z-coordinate after which scale is kept constant.");
   fFixZ->Connect("ValueSet(Double_t)", "TEveProjectionManagerEditor",
                         this, "DoFixZ()");
   AddFrame(fFixZ, new TGLayoutHints(kLHintsTop, 1, 1, 1, 0));

   fPastFixRFac = new TEveGValuator(this, "ScaleR:", 90, 0);
   fPastFixRFac->SetNELength(nel);
   fPastFixRFac->SetLabelWidth(labelW);
   fPastFixRFac->Build();
   fPastFixRFac->SetLimits(-2, 2, 101, TGNumberFormat::kNESRealTwo);
   fPastFixRFac->SetToolTip("Relative R-scale beyond FixedR.\nExpressed as 10^x.");
   fPastFixRFac->Connect("ValueSet(Double_t)", "TEveProjectionManagerEditor",
                         this, "DoPastFixRFac()");
   AddFrame(fPastFixRFac, new TGLayoutHints(kLHintsTop, 1, 1, 1, 0));

   fPastFixZFac = new TEveGValuator(this, "ScaleZ:", 90, 0);
   fPastFixZFac->SetNELength(nel);
   fPastFixZFac->SetLabelWidth(labelW);
   fPastFixZFac->Build();
   fPastFixZFac->SetLimits(-2, 2, 101, TGNumberFormat::kNESRealTwo);
   fPastFixZFac->SetToolTip("Relative Z-scale beyond FixedZ.\nExpressed as 10^x.");
   fPastFixZFac->Connect("ValueSet(Double_t)", "TEveProjectionManagerEditor",
                         this, "DoPastFixZFac()");
   AddFrame(fPastFixZFac, new TGLayoutHints(kLHintsTop, 1, 1, 1, 0));

   fCurrentDepth = new TEveGValuator(this, "CurrentZ:", 90, 0);
   fCurrentDepth->SetNELength(nel);
   fCurrentDepth->SetLabelWidth(labelW);
   fCurrentDepth->Build();
   fCurrentDepth->SetLimits(-300, 300, 601, TGNumberFormat::kNESRealTwo);
   fCurrentDepth->SetToolTip("Z coordinate of incoming projected object.");
   fCurrentDepth->Connect("ValueSet(Double_t)", "TEveProjectionManagerEditor",
                          this, "DoCurrentDepth()");
   AddFrame(fCurrentDepth, new TGLayoutHints(kLHintsTop, 1, 1, 1, 0));

   fMaxTrackStep = new TEveGValuator(this, "TrackStep:", 90, 0);
   fMaxTrackStep->SetNELength(nel);
   fMaxTrackStep->SetLabelWidth(labelW);
   fMaxTrackStep->Build();
   fMaxTrackStep->SetLimits(1, 100, 100, TGNumberFormat::kNESRealOne);
   fMaxTrackStep->SetToolTip("Maximum step between two consecutive track-points to avoid artefacts due to projective distortions.\nTaken into account automatically during projection procedure.");
   fMaxTrackStep->Connect("ValueSet(Double_t)", "TEveProjectionManagerEditor",
                         this, "DoMaxTrackStep()");
   AddFrame(fMaxTrackStep, new TGLayoutHints(kLHintsTop, 1, 1, 1, 3));

   // --------------------------------

   MakeTitle("Distortion centre");
   fCenterFrame = new TGVerticalFrame(this);

   fCenterX = new TEveGValuator(fCenterFrame, "CenterX:", 90, 0);
   fCenterX->SetNELength(nel);
   fCenterX->SetLabelWidth(labelW);
   fCenterX->Build();
   fCenterX->SetLimits(-5, 5, 501, TGNumberFormat::kNESRealThree);
   fCenterX->SetToolTip("Origin of the projection.");
   fCenterX->Connect("ValueSet(Double_t)", "TEveProjectionManagerEditor",
                     this, "DoCenter()");
   fCenterFrame->AddFrame(fCenterX, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

   fCenterY = new TEveGValuator(fCenterFrame, "CenterY:", 90, 0);
   fCenterY->SetNELength(nel);
   fCenterY->SetLabelWidth(labelW);
   fCenterY->Build();
   fCenterY->SetLimits(-5, 5, 501, TGNumberFormat::kNESRealThree);
   fCenterY->SetToolTip("Origin of the projection.");
   fCenterY->Connect("ValueSet(Double_t)", "TEveProjectionManagerEditor",
                     this, "DoCenter()");
   fCenterFrame->AddFrame(fCenterY, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

   fCenterZ = new TEveGValuator(fCenterFrame, "CenterZ:", 90, 0);
   fCenterZ->SetNELength(nel);
   fCenterZ->SetLabelWidth(labelW);
   fCenterZ->Build();
   fCenterZ->SetLimits(-25, 25, 501, TGNumberFormat::kNESRealThree);
   fCenterZ->SetToolTip("Origin of the projection.");
   fCenterZ->Connect("ValueSet(Double_t)", "TEveProjectionManagerEditor",
                     this, "DoCenter()");
   fCenterFrame->AddFrame(fCenterZ, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

   AddFrame(fCenterFrame, new TGLayoutHints(kLHintsTop, 1, 1, 1, 0));
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

void TEveProjectionManagerEditor::SetModel(TObject* obj)
{
   fM = dynamic_cast<TEveProjectionManager*>(obj);

   fType->Select(fM->GetProjection()->GetType(), kFALSE);
   fDistortion->SetValue(1000.0f * fM->GetProjection()->GetDistortion());
   fFixR->SetValue(fM->GetProjection()->GetFixR());
   fFixZ->SetValue(fM->GetProjection()->GetFixZ());
   fPastFixRFac->SetValue(fM->GetProjection()->GetPastFixRFac());
   fPastFixZFac->SetValue(fM->GetProjection()->GetPastFixZFac());
   fCurrentDepth->SetValue(fM->GetCurrentDepth());
   fMaxTrackStep->SetValue(fM->GetProjection()->GetMaxTrackStep());

   fCenterX->SetValue(fM->GetCenter().fX);
   fCenterY->SetValue(fM->GetCenter().fY);
   fCenterZ->SetValue(fM->GetCenter().fZ);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting of projection type.

void TEveProjectionManagerEditor::DoType(Int_t type)
{
   try
   {
      fM->SetProjection((TEveProjection::EPType_e)type);
      fM->ProjectChildren();
      Update();
   }
   catch (...)
   {
      SetModel(fM);
      throw;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting distortion.

void TEveProjectionManagerEditor::DoDistortion()
{
   fM->GetProjection()->SetDistortion(0.001f * fDistortion->GetValue());
   fM->UpdateName();
   fM->ProjectChildren();
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting fixed radius.

void TEveProjectionManagerEditor::DoFixR()
{
   fM->GetProjection()->SetFixR(fFixR->GetValue());
   fM->ProjectChildren();
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting fixed z-coordinate.

void TEveProjectionManagerEditor::DoFixZ()
{
   fM->GetProjection()->SetFixZ(fFixZ->GetValue());
   fM->ProjectChildren();
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting fixed radius.

void TEveProjectionManagerEditor::DoPastFixRFac()
{
   fM->GetProjection()->SetPastFixRFac(fPastFixRFac->GetValue());
   fM->ProjectChildren();
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting fixed z-coordinate.

void TEveProjectionManagerEditor::DoPastFixZFac()
{
   fM->GetProjection()->SetPastFixZFac(fPastFixZFac->GetValue());
   fM->ProjectChildren();
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting current depth.

void TEveProjectionManagerEditor::DoCurrentDepth()
{
   fM->SetCurrentDepth(fCurrentDepth->GetValue());
   fM->ProjectChildren();
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting fixed z-coordinate.

void TEveProjectionManagerEditor::DoMaxTrackStep()
{
   fM->GetProjection()->SetMaxTrackStep(fMaxTrackStep->GetValue());
   fM->ProjectChildren();
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting center of distortion.

void TEveProjectionManagerEditor::DoCenter()
{
   fM->SetCenter(fCenterX->GetValue(), fCenterY->GetValue(), fCenterZ->GetValue());
   Update();
}

