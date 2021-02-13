// @(#):$Id$
// Author: M.Gheata

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoTrapEditor
\ingroup Geometry_builder

Editor for a TGeoTrap.

\image html geom_trap_pic.png

\image html geom_trap_ed.png

*/

#include "TGeoTrapEditor.h"
#include "TGeoTabManager.h"
#include "TGeoArb8.h"
#include "TGeoManager.h"
#include "TVirtualGeoPainter.h"
#include "TVirtualPad.h"
#include "TView.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"

ClassImp(TGeoTrapEditor);

enum ETGeoTrapWid {
   kTRAP_NAME, kTRAP_H1, kTRAP_BL1, kTRAP_TL1, kTRAP_DZ, kTRAP_ALPHA1,
   kTRAP_SC1, kTRAP_SC2, kTRAP_THETA, kTRAP_PHI, kTRAP_APPLY, kTRAP_UNDO
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor for para editor

TGeoTrapEditor::TGeoTrapEditor(const TGWindow *p, Int_t width,
                                   Int_t height, UInt_t options, Pixel_t back)
   : TGeoGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fShape   = 0;
   fH1i = fBl1i = fTl1i = fDzi = fAlpha1i = fThetai = fPhii = fSci = 0.0;
   fNamei = "";
   fIsModified = kFALSE;
   fIsShapeEditable = kTRUE;

   // TextEntry for shape name
   MakeTitle("Name");
   fShapeName = new TGTextEntry(this, new TGTextBuffer(50), kTRAP_NAME);
   fShapeName->Resize(135, fShapeName->GetDefaultHeight());
   fShapeName->SetToolTipText("Enter the parallelepiped name");
   fShapeName->Associate(this);
   AddFrame(fShapeName, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));

   TGTextEntry *nef;
   MakeTitle("Dimensions");
   // Number entry for H1
   TGCompositeFrame *f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "DY"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEH1 = new TGNumberEntry(f1, 0., 5, kTRAP_H1);
   fEH1->SetNumAttr(TGNumberFormat::kNEAPositive);
   fEH1->Resize(100, fEH1->GetDefaultHeight());
   nef = (TGTextEntry*)fEH1->GetNumberEntry();
   nef->SetToolTipText("Enter the half length in y at low z");
   fEH1->Associate(this);
   f1->AddFrame(fEH1, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // Number entry for Bl1
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "DX1"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEBl1 = new TGNumberEntry(f1, 0., 5, kTRAP_BL1);
   fEBl1->SetNumAttr(TGNumberFormat::kNEAPositive);
   fEBl1->Resize(100, fEBl1->GetDefaultHeight());
   nef = (TGTextEntry*)fEBl1->GetNumberEntry();
   nef->SetToolTipText("Enter the half length in x at low z and y low edge");
   fEBl1->Associate(this);
   f1->AddFrame(fEBl1, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // Number entry for Tl1
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "DX2"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fETl1 = new TGNumberEntry(f1, 0., 5, kTRAP_TL1);
   fETl1->SetNumAttr(TGNumberFormat::kNEAPositive);
   fETl1->Resize(100, fETl1->GetDefaultHeight());
   nef = (TGTextEntry*)fETl1->GetNumberEntry();
   nef->SetToolTipText("Enter the half length in x at low z and y high edge");
   fETl1->Associate(this);
   f1->AddFrame(fETl1, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

    // Number entry for scale factor
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "SC1"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fESc1 = new TGNumberEntry(f1, 0., 5, kTRAP_SC1);
   fESc1->SetNumAttr(TGNumberFormat::kNEAPositive);
   fESc1->Resize(100, fESc1->GetDefaultHeight());
   nef = (TGTextEntry*)fESc1->GetNumberEntry();
   nef->SetToolTipText("Enter the scale factor for lower Z face");
   fESc1->Associate(this);
   f1->AddFrame(fESc1, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

    // Number entry for scale factor
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "SC2"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fESc2 = new TGNumberEntry(f1, 0., 5, kTRAP_SC2);
   fESc2->SetNumAttr(TGNumberFormat::kNEAPositive);
   fESc2->Resize(100, fESc2->GetDefaultHeight());
   nef = (TGTextEntry*)fESc2->GetNumberEntry();
   nef->SetToolTipText("Enter the scale factor for upper Z face");
   fESc2->Associate(this);
   f1->AddFrame(fESc2, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

  // Number entry for dz
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "DZ"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEDz = new TGNumberEntry(f1, 0., 5, kTRAP_DZ);
   fEDz->SetNumAttr(TGNumberFormat::kNEAPositive);
   fEDz->Resize(100, fEDz->GetDefaultHeight());
   nef = (TGTextEntry*)fEDz->GetNumberEntry();
   nef->SetToolTipText("Enter the half-length in Z");
   fEDz->Associate(this);
   f1->AddFrame(fEDz, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // Number entry for Alpha1
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "ALPHA"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEAlpha1 = new TGNumberEntry(f1, 0., 5, kTRAP_ALPHA1);
   fEAlpha1->Resize(100, fEAlpha1->GetDefaultHeight());
   nef = (TGTextEntry*)fEAlpha1->GetNumberEntry();
   nef->SetToolTipText("Enter  angle between centers of x edges an y axis at low z");
   fEAlpha1->Associate(this);
   f1->AddFrame(fEAlpha1, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // Number entry for Theta
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Theta"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fETheta = new TGNumberEntry(f1, 0., 5, kTRAP_THETA);
   fETheta->SetNumAttr(TGNumberFormat::kNEAPositive);
   fETheta->Resize(100, fETheta->GetDefaultHeight());
   nef = (TGTextEntry*)fETheta->GetNumberEntry();
   nef->SetToolTipText("Enter initial  theta");
   fETheta->Associate(this);
   f1->AddFrame(fETheta, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

    // Number entry for Phi
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Phi"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEPhi = new TGNumberEntry(f1, 0., 5, kTRAP_PHI);
   fEPhi->SetNumAttr(TGNumberFormat::kNEAPositive);
   fEPhi->Resize(100, fEPhi->GetDefaultHeight());
   nef = (TGTextEntry*)fEPhi->GetNumberEntry();
   nef->SetToolTipText("Enter initial  phi");
   fEPhi->Associate(this);
   f1->AddFrame(fEPhi, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // Delayed draw
   fDFrame = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth | kSunkenFrame);
   fDelayed = new TGCheckButton(fDFrame, "Delayed draw");
   fDFrame->AddFrame(fDelayed, new TGLayoutHints(kLHintsLeft , 2, 2, 4, 4));
   AddFrame(fDFrame,  new TGLayoutHints(kLHintsLeft, 6, 6, 4, 4));

   // Buttons
   fBFrame = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   fApply = new TGTextButton(fBFrame, "Apply");
   fBFrame->AddFrame(fApply, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));
   fApply->Associate(this);
   fUndo = new TGTextButton(fBFrame, "Undo");
   fBFrame->AddFrame(fUndo, new TGLayoutHints(kLHintsRight , 2, 2, 4, 4));
   fUndo->Associate(this);
   AddFrame(fBFrame,  new TGLayoutHints(kLHintsLeft, 6, 6, 4, 4));
   fUndo->SetSize(fApply->GetSize());
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoTrapEditor::~TGeoTrapEditor()
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

void TGeoTrapEditor::ConnectSignals2Slots()
{
   fApply->Connect("Clicked()", "TGeoTrapEditor", this, "DoApply()");
   fUndo->Connect("Clicked()", "TGeoTrapEditor", this, "DoUndo()");
   fShapeName->Connect("TextChanged(const char *)", "TGeoTrapEditor", this, "DoModified()");
   fEH1->Connect("ValueSet(Long_t)", "TGeoTrapEditor", this, "DoH1()");
   fEBl1->Connect("ValueSet(Long_t)", "TGeoTrapEditor", this, "DoBl1()");
   fETl1->Connect("ValueSet(Long_t)", "TGeoTrapEditor", this, "DoTl1()");
   fEDz->Connect("ValueSet(Long_t)", "TGeoTrapEditor", this, "DoDz()");
   fESc1->Connect("ValueSet(Long_t)", "TGeoTrapEditor", this, "DoSc1()");
   fESc2->Connect("ValueSet(Long_t)", "TGeoTrapEditor", this, "DoSc2()");
   fEAlpha1->Connect("ValueSet(Long_t)", "TGeoTrapEditor", this, "DoAlpha1()");
   fETheta->Connect("ValueSet(Long_t)", "TGeoTrapEditor", this, "DoTheta()");
   fEPhi->Connect("ValueSet(Long_t)", "TGeoTrapEditor", this, "DoPhi()");
   fEH1->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTrapEditor", this, "DoModified()");
   fEBl1->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTrapEditor", this, "DoModified()");
   fETl1->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTrapEditor", this, "DoModified()");
   fEDz->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTrapEditor", this, "DoModified()");
   fESc1->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTrapEditor", this, "DoModified()");
   fESc2->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTrapEditor", this, "DoModified()");
   fEAlpha1->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTrapEditor", this, "DoModified()");
   fETheta->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTrapEditor", this, "DoModified()");
   fEPhi->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTrapEditor", this, "DoModified()");
   fInit = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Connect to the selected object.

void TGeoTrapEditor::SetModel(TObject* obj)
{
   if (obj == 0 || (obj->IsA()!=TGeoTrap::Class())) {
      SetActive(kFALSE);
      return;
   }
   fShape = (TGeoTrap*)obj;
   fH1i = fShape->GetH1();
   fBl1i = fShape->GetBl1();
   fTl1i = fShape->GetTl1();
   fDzi = fShape->GetDz();
   Double_t h2i = fShape->GetH2();
//   Double_t bl2i = fShape->GetBl2();
//   Double_t tl2i = fShape->GetTl2();
   fSci = h2i/fH1i;
   fAlpha1i = fShape->GetAlpha1();
   fThetai = fShape->GetTheta();
   fPhii = fShape->GetPhi();
   const char *sname = fShape->GetName();
   if (!strcmp(sname, fShape->ClassName())) fShapeName->SetText("-no_name");
   else {
      fShapeName->SetText(sname);
      fNamei = sname;
   }
   fEH1->SetNumber(fH1i);
   fEBl1->SetNumber(fBl1i);
   fETl1->SetNumber(fTl1i);
   fEDz->SetNumber(fDzi);
   fESc1->SetNumber(1.);
   fESc2->SetNumber(fSci);
   fEAlpha1->SetNumber(fAlpha1i);
   fETheta->SetNumber(fThetai);
   fEPhi->SetNumber(fPhii);
   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);

   if (fInit) ConnectSignals2Slots();
   SetActive();
}

////////////////////////////////////////////////////////////////////////////////
/// Check if shape drawing is delayed.

Bool_t TGeoTrapEditor::IsDelayed() const
{
   return (fDelayed->GetState() == kButtonDown);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for name.

void TGeoTrapEditor::DoName()
{
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for applying current settings.

void TGeoTrapEditor::DoApply()
{
   const char *name = fShapeName->GetText();
   if (strcmp(name,fShape->GetName())) fShape->SetName(name);
   Double_t sc1 = fESc1->GetNumber();
   Double_t sc2 = fESc2->GetNumber();
   Double_t h1 = sc1*fEH1->GetNumber();
   Double_t bl1 = sc1*fEBl1->GetNumber();
   Double_t tl1 = sc1*fETl1->GetNumber();
   Double_t h2 = sc2*fEH1->GetNumber();
   Double_t bl2 = sc2*fEBl1->GetNumber();
   Double_t tl2 = sc2*fETl1->GetNumber();
   Double_t dz = fEDz->GetNumber();
   Double_t alpha1 = fEAlpha1->GetNumber();
   Double_t theta = fETheta->GetNumber();
   Double_t phi = fEPhi->GetNumber();
   Double_t param[11];
   param[0] = dz;
   param[1] = theta;
   param[2] = phi;
   param[3] = h1;
   param[7] = h2;
   param[4] = bl1;
   param[8] = bl2;
   param[5] = tl1;
   param[9] = tl2;
   param[6] = alpha1;
   param[10] = alpha1;
   fShape->SetDimensions(param);
   fShape->ComputeBBox();
   fUndo->SetEnabled();
   fApply->SetEnabled(kFALSE);
   if (fPad) {
      if (gGeoManager && gGeoManager->GetPainter() && gGeoManager->GetPainter()->IsPaintingShape()) {
         TView *view = fPad->GetView();
         if (!view) {
            fShape->Draw();
            fPad->GetView()->ShowAxis();
         } else {
            view->SetRange(-fShape->GetDX(), -fShape->GetDY(), -fShape->GetDZ(),
                           fShape->GetDX(), fShape->GetDY(), fShape->GetDZ());
            Update();
         }
      } else Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for notifying modifications.

void TGeoTrapEditor::DoModified()
{
   fApply->SetEnabled();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for undoing last operation.

void TGeoTrapEditor::DoUndo()
{
   fEH1->SetNumber(fH1i);
   fEBl1->SetNumber(fBl1i);
   fETl1->SetNumber(fTl1i);
   fESc1->SetNumber(1.);
   fESc2->SetNumber(fSci);
   fEDz->SetNumber(fDzi);
   fEAlpha1->SetNumber(fAlpha1i);
   fETheta->SetNumber(fThetai);
   fEPhi->SetNumber(fPhii);
   DoApply();
   fUndo->SetEnabled(kFALSE);
   fApply->SetEnabled(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for H1.

void TGeoTrapEditor::DoH1()
{
   Double_t h1 = fEH1->GetNumber();
   if (h1<=0) {
      h1 = 0.1;
      fEH1->SetNumber(h1);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for Bl1.

void TGeoTrapEditor::DoBl1()
{
   Double_t bl1 = fEBl1->GetNumber();
   if (bl1<=0) {
      bl1 = 0.1;
      fEBl1->SetNumber(bl1);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for Tl1.

void TGeoTrapEditor::DoTl1()
{
   Double_t tl1 = fETl1->GetNumber();
   if (tl1<=0) {
      tl1 = 0.1;
      fETl1->SetNumber(tl1);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for Z.

void TGeoTrapEditor::DoDz()
{
   Double_t dz = fEDz->GetNumber();
   if (dz<=0) {
      dz = 0.1;
      fEDz->SetNumber(dz);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for H2.

void TGeoTrapEditor::DoSc1()
{
   Double_t sc1 = fESc1->GetNumber();
   if (sc1<=0) {
      sc1 = 0.1;
      fESc1->SetNumber(sc1);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for H2.

void TGeoTrapEditor::DoSc2()
{
   Double_t sc2 = fESc2->GetNumber();
   if (sc2<=0) {
      sc2 = 0.1;
      fESc2->SetNumber(sc2);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for alpha1.

void TGeoTrapEditor::DoAlpha1()
{
   Double_t alpha1 = fEAlpha1->GetNumber();
   if (TMath::Abs(alpha1)>=90) {
      alpha1 = 89.9*TMath::Sign(1.,alpha1);
      fEAlpha1->SetNumber(alpha1);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for theta.

void TGeoTrapEditor::DoTheta()
{
   Double_t theta = fETheta->GetNumber();
   if (theta<0) {
      theta = 0;
      fETheta->SetNumber(theta);
   }
   if (theta>180) {
      theta = 180;
      fETheta->SetNumber(theta);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for phi.

void TGeoTrapEditor::DoPhi()
{
   Double_t phi = fEPhi->GetNumber();
   if (phi<0 || phi>360) {
      phi = 0;
      fEPhi->SetNumber(phi);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

ClassImp(TGeoGtraEditor);

enum ETGeoGtraWid {
   kGTRA_TWIST
};

/** \class TGeoGtraEditor
\ingroup Geometry_builder

Editor for a TGeoGtra.

\image html geom_gtra_pic.png

\image html geom_gtra_ed.png

*/

////////////////////////////////////////////////////////////////////////////////
/// Constructor for gtra editor

TGeoGtraEditor::TGeoGtraEditor(const TGWindow *p, Int_t width,
                               Int_t height, UInt_t options, Pixel_t back)
   : TGeoTrapEditor(p, width, height, options, back)
{
   fTwisti = 0;
   TGTextEntry *nef;
   // Number entry for Twist angle
   TGCompositeFrame *f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "TWIST"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fETwist = new TGNumberEntry(f1, 0., 5, kGTRA_TWIST);
   fETwist->Resize(100, fETwist->GetDefaultHeight());
   nef = (TGTextEntry*)fETwist->GetNumberEntry();
   nef->SetToolTipText("Enter twist angle");
   fETwist->Associate(this);
   f1->AddFrame(fETwist, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));
   TGeoTabManager::MoveFrame(fDFrame, this);
   TGeoTabManager::MoveFrame(fBFrame, this);
   fETwist->Connect("ValueSet(Long_t)", "TGeoGtraEditor", this, "DoTwist()");
   nef->Connect("TextChanged(const char *)", "TGeoGtraEditor", this, "DoModified()");
}
////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoGtraEditor::~TGeoGtraEditor()
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
/// Connect to a given twisted trapezoid.

void TGeoGtraEditor::SetModel(TObject* obj)
{
   if (obj == 0 || (obj->IsA()!=TGeoGtra::Class())) {
      SetActive(kFALSE);
      return;
   }
   fShape = (TGeoTrap*)obj;
   fH1i = fShape->GetH1();
   fBl1i = fShape->GetBl1();
   fTl1i = fShape->GetTl1();
   fDzi = fShape->GetDz();
   Double_t h2i = fShape->GetH2();
//   Double_t bl2i = fShape->GetBl2();
//   Double_t tl2i = fShape->GetTl2();
   fSci = h2i/fH1i;
   fAlpha1i = fShape->GetAlpha1();
   fThetai = fShape->GetTheta();
   fPhii = fShape->GetPhi();
   fTwisti = ((TGeoGtra*)fShape)->GetTwistAngle();
   const char *sname = fShape->GetName();
   if (!strcmp(sname, fShape->ClassName())) fShapeName->SetText("-no_name");
   else {
      fShapeName->SetText(sname);
      fNamei = sname;
   }
   fEH1->SetNumber(fH1i);
   fEBl1->SetNumber(fBl1i);
   fETl1->SetNumber(fTl1i);
   fEDz->SetNumber(fDzi);
   fESc1->SetNumber(1.);
   fESc2->SetNumber(fSci);
   fEAlpha1->SetNumber(fAlpha1i);
   fETheta->SetNumber(fThetai);
   fEPhi->SetNumber(fPhii);
   fETwist->SetNumber(fTwisti);
   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);

   if (fInit) ConnectSignals2Slots();
   SetActive();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for applying current settings.

void TGeoGtraEditor::DoApply()
{
   const char *name = fShapeName->GetText();
   if (strcmp(name,fShape->GetName())) fShape->SetName(name);
   Double_t sc1 = fESc1->GetNumber();
   Double_t sc2 = fESc2->GetNumber();
   Double_t h1 = sc1*fEH1->GetNumber();
   Double_t bl1 = sc1*fEBl1->GetNumber();
   Double_t tl1 = sc1*fETl1->GetNumber();
   Double_t h2 = sc2*fEH1->GetNumber();
   Double_t bl2 = sc2*fEBl1->GetNumber();
   Double_t tl2 = sc2*fETl1->GetNumber();
   Double_t dz = fEDz->GetNumber();
   Double_t alpha1 = fEAlpha1->GetNumber();
   Double_t theta = fETheta->GetNumber();
   Double_t phi = fEPhi->GetNumber();
   Double_t twist = fETwist->GetNumber();
   Double_t param[12];
   param[0] = dz;
   param[1] = theta;
   param[2] = phi;
   param[3] = h1;
   param[7] = h2;
   param[4] = bl1;
   param[8] = bl2;
   param[5] = tl1;
   param[9] = tl2;
   param[6] = alpha1;
   param[10] = alpha1;
   param[11] = twist;
   TGeoGtra *shape = (TGeoGtra*)fShape;
   shape->SetDimensions(param);
   shape->ComputeBBox();
   fUndo->SetEnabled();
   fApply->SetEnabled(kFALSE);
   if (fPad) {
      if (gGeoManager && gGeoManager->GetPainter() && gGeoManager->GetPainter()->IsPaintingShape()) {
         TView *view = fPad->GetView();
         if (!view) {
            fShape->Draw();
            fPad->GetView()->ShowAxis();
         } else {
            view->SetRange(-fShape->GetDX(), -fShape->GetDY(), -fShape->GetDZ(),
                           fShape->GetDX(), fShape->GetDY(), fShape->GetDZ());
            Update();
         }
      } else Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for undoing last operation.

void TGeoGtraEditor::DoUndo()
{
   fEH1->SetNumber(fH1i);
   fEBl1->SetNumber(fBl1i);
   fETl1->SetNumber(fTl1i);
   fESc1->SetNumber(1.);
   fESc2->SetNumber(fSci);
   fEDz->SetNumber(fDzi);
   fEAlpha1->SetNumber(fAlpha1i);
   fETheta->SetNumber(fThetai);
   fEPhi->SetNumber(fPhii);
   fETwist->SetNumber(fTwisti);
   DoApply();
   fUndo->SetEnabled(kFALSE);
   fApply->SetEnabled(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Change the twist angle.

void TGeoGtraEditor::DoTwist()
{
   Double_t twist = fETwist->GetNumber();
   if (twist<=-180 || twist>=180) {
      twist = 0.;
      fETwist->SetNumber(twist);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}
