// @(#):$Id$
// Author: M.Gheata

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoSphereEditor
\ingroup Geometry_builder

Editor for a TGeoSphere.

\image html geom_sphe_pic.png

\image html geom_sphe_ed.png

*/

#include "TGeoSphereEditor.h"
#include "TGeoTabManager.h"
#include "TGeoSphere.h"
#include "TGeoManager.h"
#include "TVirtualGeoPainter.h"
#include "TVirtualPad.h"
#include "TView.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"
#include "TGDoubleSlider.h"

ClassImp(TGeoSphereEditor);

enum ETGeoSphereWid {
   kSPHERE_NAME, kSPHERE_RMIN, kSPHERE_RMAX, kSPHERE_THETA1,
   kSPHERE_THETA2, kSPHERE_PHI1, kSPHERE_PHI2, kSPHERE_PHI, kSPHERE_THETA,
   kSPHERE_APPLY, kSPHERE_UNDO
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor for sphere editor

TGeoSphereEditor::TGeoSphereEditor(const TGWindow *p, Int_t width,
                                   Int_t height, UInt_t options, Pixel_t back)
   : TGeoGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fShape   = 0;
   fRmini = fRmaxi = fTheta1i = fTheta2i = fPhi1i = fPhi2i = 0.0;
   fNamei = "";
   fIsModified = kFALSE;
   fIsShapeEditable = kTRUE;
   fLock = kFALSE;

   // TextEntry for shape name
   MakeTitle("Name");
   fShapeName = new TGTextEntry(this, new TGTextBuffer(50), kSPHERE_NAME);
   fShapeName->Resize(135, fShapeName->GetDefaultHeight());
   fShapeName->SetToolTipText("Enter the sphere name");
   fShapeName->Associate(this);
   AddFrame(fShapeName, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));

   TGTextEntry *nef;
   MakeTitle("Sphere dimensions");
   TGCompositeFrame *compxyz = new TGCompositeFrame(this, 118, 30, kVerticalFrame | kRaisedFrame);
   // Number entry for rmin
   TGCompositeFrame *f1 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kOwnBackground);
   f1->AddFrame(new TGLabel(f1, "Rmin"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fERmin = new TGNumberEntry(f1, 0., 5, kSPHERE_RMIN);
   fERmin->SetNumAttr(TGNumberFormat::kNEANonNegative);
   nef = (TGTextEntry*)fERmin->GetNumberEntry();
   nef->SetToolTipText("Enter the inner radius");
   fERmin->Associate(this);
   fERmin->Resize(100, fERmin->GetDefaultHeight());
   f1->AddFrame(fERmin, new TGLayoutHints(kLHintsRight , 2, 2, 4, 4));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));

   // Number entry for Rmax
   f1 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kOwnBackground);
   f1->AddFrame(new TGLabel(f1, "Rmax"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fERmax = new TGNumberEntry(f1, 0., 5, kSPHERE_RMAX);
   fERmax->SetNumAttr(TGNumberFormat::kNEANonNegative);
   nef = (TGTextEntry*)fERmax->GetNumberEntry();
   nef->SetToolTipText("Enter the outer radius");
   fERmax->Associate(this);
   fERmax->Resize(100, fERmax->GetDefaultHeight());
   f1->AddFrame(fERmax, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));
   AddFrame(compxyz, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));

   MakeTitle("Phi/theta range");
   TGCompositeFrame *f11 = new TGCompositeFrame(this, 150,200, kHorizontalFrame);
   compxyz = new TGCompositeFrame(f11, 75, 200, kHorizontalFrame | kRaisedFrame);
   // Vertical slider
   fSPhi = new TGDoubleVSlider(compxyz,140);
   fSPhi->SetRange(0.,720.);
   compxyz->AddFrame(fSPhi, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 2, 4, 4));
   f1 = new TGCompositeFrame(compxyz, 50, 200, kVerticalFrame);
   f1->AddFrame(new TGLabel(f1, "Phi min."), new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 2));
   fEPhi1 = new TGNumberEntry(f1, 0., 5, kSPHERE_PHI1);
   fEPhi1->Resize(30, fEPhi1->GetDefaultHeight());
   fEPhi1->SetNumAttr(TGNumberFormat::kNEANonNegative);
   nef = (TGTextEntry*)fEPhi1->GetNumberEntry();
   nef->SetToolTipText("Enter the phi1 value");
   fEPhi1->Associate(this);
   f1->AddFrame(fEPhi1, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 2, 2, 2, 2));

   fEPhi2 = new TGNumberEntry(f1, 0., 5, kSPHERE_PHI2);
   fEPhi2->Resize(30, fEPhi2->GetDefaultHeight());
   fEPhi2->SetNumAttr(TGNumberFormat::kNEANonNegative);
   nef = (TGTextEntry*)fEPhi2->GetNumberEntry();
   nef->SetToolTipText("Enter the phi2 value");
   fEPhi2->Associate(this);
   fEPhi2->Resize(30, fEPhi2->GetDefaultHeight());
   f1->AddFrame(fEPhi2, new TGLayoutHints(kLHintsBottom | kLHintsLeft | kLHintsExpandX, 2, 2, 2, 2));
   f1->AddFrame(new TGLabel(f1, "Phi max."), new TGLayoutHints(kLHintsBottom | kLHintsLeft, 2, 2, 2, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

   compxyz->Resize(75,150);
   f11->AddFrame(compxyz, new TGLayoutHints(kLHintsLeft, 0,0,0,0));

   compxyz = new TGCompositeFrame(f11, 75, 200, kHorizontalFrame | kRaisedFrame);
   // Vertical slider
   fSTheta = new TGDoubleVSlider(compxyz,140);
   fSTheta->SetRange(0.,180.);
   compxyz->AddFrame(fSTheta, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 2, 4, 4));
   f1 = new TGCompositeFrame(compxyz, 50, 200, kVerticalFrame);
   f1->AddFrame(new TGLabel(f1, "Theta min."), new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 2));
   fETheta1 = new TGNumberEntry(f1, 0., 5, kSPHERE_THETA1);
   fETheta1->Resize(30, fETheta1->GetDefaultHeight());
   fETheta1->SetNumAttr(TGNumberFormat::kNEANonNegative);
   nef = (TGTextEntry*)fETheta1->GetNumberEntry();
   nef->SetToolTipText("Enter the theta1 value");
   fETheta1->Associate(this);
   f1->AddFrame(fETheta1, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 2, 2, 2, 2));

   fETheta2 = new TGNumberEntry(f1, 0., 5, kSPHERE_THETA2);
   fETheta2->Resize(30, fETheta2->GetDefaultHeight());
   fETheta2->SetNumAttr(TGNumberFormat::kNEANonNegative);
   nef = (TGTextEntry*)fETheta2->GetNumberEntry();
   nef->SetToolTipText("Enter the theta2 value");
   fETheta2->Associate(this);
   f1->AddFrame(fETheta2, new TGLayoutHints(kLHintsBottom | kLHintsLeft | kLHintsExpandX, 2, 2, 2, 2));
   f1->AddFrame(new TGLabel(f1, "Theta max."), new TGLayoutHints(kLHintsBottom | kLHintsLeft, 2, 2, 2, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

   compxyz->Resize(75,150);
   f11->AddFrame(compxyz, new TGLayoutHints(kLHintsRight, 0, 0, 0, 0));

   AddFrame(f11, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));

   // Delayed draw
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth | kSunkenFrame);
   fDelayed = new TGCheckButton(f1, "Delayed draw");
   f1->AddFrame(fDelayed, new TGLayoutHints(kLHintsLeft , 2, 2, 4, 4));
   AddFrame(f1,  new TGLayoutHints(kLHintsLeft, 6, 6, 4, 4));

   // Buttons
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   fApply = new TGTextButton(f1, "Apply");
   f1->AddFrame(fApply, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));
   fApply->Associate(this);
   fUndo = new TGTextButton(f1, "Undo");
   f1->AddFrame(fUndo, new TGLayoutHints(kLHintsRight , 2, 2, 4, 4));
   fUndo->Associate(this);
   AddFrame(f1,  new TGLayoutHints(kLHintsLeft, 6, 6, 4, 4));
   fUndo->SetSize(fApply->GetSize());
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGeoSphereEditor::~TGeoSphereEditor()
{
   TGFrameElement *el;
   TIter next(GetList());
   while ((el = (TGFrameElement *)next())) {
      if (el->fFrame->IsComposite())
         TGeoTabManager::Cleanup((TGCompositeFrame*)el->fFrame);
   }
   Cleanup();
}

////////////////////////////////////////////////////////////////////////////////
/// Connect signals to slots.

void TGeoSphereEditor::ConnectSignals2Slots()
{
   fApply->Connect("Clicked()", "TGeoSphereEditor", this, "DoApply()");
   fUndo->Connect("Clicked()", "TGeoSphereEditor", this, "DoUndo()");
   fShapeName->Connect("TextChanged(const char *)", "TGeoSphereEditor", this, "DoModified()");
   fERmin->Connect("ValueSet(Long_t)", "TGeoSphereEditor", this, "DoRmin()");
   fERmax->Connect("ValueSet(Long_t)", "TGeoSphereEditor", this, "DoRmax()");
   fEPhi1->Connect("ValueSet(Long_t)", "TGeoSphereEditor", this, "DoPhi1()");
   fEPhi2->Connect("ValueSet(Long_t)", "TGeoSphereEditor", this, "DoPhi2()");
   fETheta1->Connect("ValueSet(Long_t)", "TGeoSphereEditor", this, "DoTheta1()");
   fETheta2->Connect("ValueSet(Long_t)", "TGeoSphereEditor", this, "DoTheta2()");
   fSPhi->Connect("PositionChanged()", "TGeoSphereEditor", this, "DoPhi()");
   fSTheta->Connect("PositionChanged()", "TGeoSphereEditor", this, "DoTheta()");
   fInit = kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Connect to a given sphere.

void TGeoSphereEditor::SetModel(TObject* obj)
{
   if (obj == 0 || (obj->IsA()!=TGeoSphere::Class())) {
      SetActive(kFALSE);
      return;
   }
   fShape = (TGeoSphere*)obj;
   fRmini = fShape->GetRmin();
   fRmaxi = fShape->GetRmax();
   fPhi1i = fShape->GetPhi1();
   fPhi2i = fShape->GetPhi2();
   fTheta1i = fShape->GetTheta1();
   fTheta2i = fShape->GetTheta2();
   fNamei = fShape->GetName();
   fShapeName->SetText(fShape->GetName());
   fERmin->SetNumber(fRmini);
   fERmax->SetNumber(fRmaxi);
   fEPhi1->SetNumber(fPhi1i);
   fEPhi2->SetNumber(fPhi2i);
   fETheta1->SetNumber(fTheta1i);
   fETheta2->SetNumber(fTheta2i);
   fSPhi->SetPosition(fPhi1i, fPhi2i);
   fSTheta->SetPosition(fTheta1i, fTheta2i);

   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);

   if (fInit) ConnectSignals2Slots();
   SetActive();
}

////////////////////////////////////////////////////////////////////////////////
/// Check if shape drawing is delayed.

Bool_t TGeoSphereEditor::IsDelayed() const
{
   return (fDelayed->GetState() == kButtonDown);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for name.

void TGeoSphereEditor::DoName()
{
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for applying modifications.

void TGeoSphereEditor::DoApply()
{
   const char *name = fShapeName->GetText();
   if (strcmp(name,fShape->GetName())) fShape->SetName(name);
   Double_t rmin = fERmin->GetNumber();
   Double_t rmax = fERmax->GetNumber();
   Double_t phi1 = fEPhi1->GetNumber();
   Double_t phi2 = fEPhi2->GetNumber();
   if ((phi2-phi1) > 360.001) {
      phi1 = 0.;
      phi2 = 360.;
      fEPhi1->SetNumber(phi1);
      fEPhi2->SetNumber(phi2);
      fLock = kTRUE;
      fSPhi->SetPosition(phi1,phi2);
      fLock = kFALSE;
   }
   Double_t theta1 = fETheta1->GetNumber();
   Double_t theta2 = fETheta2->GetNumber();
   fShape->SetSphDimensions(rmin, rmax, theta1,theta2,phi1,phi2);
   fShape->ComputeBBox();
   fUndo->SetEnabled();
   fApply->SetEnabled(kFALSE);
   if (fPad) {
      if (gGeoManager && gGeoManager->GetPainter() && gGeoManager->GetPainter()->IsPaintingShape()) {
         fShape->Draw();
         fPad->GetView()->ShowAxis();
      } else Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for signaling modifications.

void TGeoSphereEditor::DoModified()
{
   fApply->SetEnabled();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for undoing last operation.

void TGeoSphereEditor::DoUndo()
{
   fERmin->SetNumber(fRmini);
   fERmax->SetNumber(fRmaxi);
   fEPhi1->SetNumber(fPhi1i);
   fEPhi2->SetNumber(fPhi2i);
   fSPhi->SetPosition(fPhi1i,fPhi2i);
   fETheta1->SetNumber(fTheta1i);
   fETheta2->SetNumber(fTheta2i);
   fSTheta->SetPosition(fTheta1i,fTheta2i);
   DoApply();
   fUndo->SetEnabled(kFALSE);
   fApply->SetEnabled(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for Rmin.

void TGeoSphereEditor::DoRmin()
{
   Double_t rmin = fERmin->GetNumber();
   Double_t rmax = fERmax->GetNumber();
   if (rmin <= 0.) {
      rmin = 0.;
      fERmin->SetNumber(rmin);
   }
   if (rmin >= rmax) {
      rmin = rmax - 0.1;
      fERmin->SetNumber(rmin);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for Rmax.

void TGeoSphereEditor::DoRmax()
{
   Double_t rmin = fERmin->GetNumber();
   Double_t rmax = fERmax->GetNumber();
   if (rmax <= 0.) {
      rmax = 0.1;
      fERmax->SetNumber(rmax);
   }
   if (rmax < rmin+1.e-10) {
      rmax = rmin + 0.1;
      if (rmin < 0.) rmin = 0.;
      fERmax->SetNumber(rmax);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for phi1.

void TGeoSphereEditor::DoPhi1()
{
   Double_t phi1 = fEPhi1->GetNumber();
   Double_t phi2 = fEPhi2->GetNumber();
   if (phi1 > 360-1.e-10) {
      phi1 = 0.;
      fEPhi1->SetNumber(phi1);
   }
   if (phi2<phi1+1.e-10) {
      phi1 = phi2 - 0.1;
      fEPhi1->SetNumber(phi1);
   }
   if (!fLock) {
      DoModified();
      fLock = kTRUE;
      fSPhi->SetPosition(phi1,phi2);
   } else fLock = kFALSE;
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for phi2.

void TGeoSphereEditor::DoPhi2()
{
   Double_t phi1 = fEPhi1->GetNumber();
   Double_t phi2 = fEPhi2->GetNumber();
   if (phi2-phi1 > 360.) {
      phi2 -= 360.;
      fEPhi2->SetNumber(phi2);
   }
   if (phi2<phi1+1.e-10) {
      phi2 = phi1 + 0.1;
      fEPhi2->SetNumber(phi2);
   }
   if (!fLock) {
      DoModified();
      fLock = kTRUE;
      fSPhi->SetPosition(phi1,phi2);
   } else fLock = kFALSE;
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for phi slider.

void TGeoSphereEditor::DoPhi()
{
   if (!fLock) {
      DoModified();
      fLock = kTRUE;
      fEPhi1->SetNumber(fSPhi->GetMinPosition());
      fLock = kTRUE;
      fEPhi2->SetNumber(fSPhi->GetMaxPosition());
   } else fLock = kFALSE;
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for theta1.

void TGeoSphereEditor::DoTheta1()
{
   Double_t theta1 = fETheta1->GetNumber();
   Double_t theta2 = fETheta2->GetNumber();
   if (theta2<theta1+1.e-10) {
      theta2 = theta1 + 0.1;
      fETheta2->SetNumber(theta2);
   }
   if (!fLock) {
      DoModified();
      fLock = kTRUE;
      fSTheta->SetPosition(theta1,theta2);
   } else fLock = kFALSE;
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for theta2.

void TGeoSphereEditor::DoTheta2()
{
   Double_t theta1 = fETheta1->GetNumber();
   Double_t theta2 = fETheta2->GetNumber();
   if (theta2<theta1+1.e-10) {
      theta1 = theta2 - 0.1;
      fETheta1->SetNumber(theta1);
   }
   if (!fLock) {
      DoModified();
      fLock = kTRUE;
      fSTheta->SetPosition(theta1,theta2);
   } else fLock = kFALSE;
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for theta slider.

void TGeoSphereEditor::DoTheta()
{
   if (!fLock) {
      DoModified();
      fLock = kTRUE;
      fETheta1->SetNumber(fSTheta->GetMinPosition());
      fLock = kTRUE;
      fETheta2->SetNumber(fSTheta->GetMaxPosition());
   } else fLock = kFALSE;
   if (!IsDelayed()) DoApply();
}
