// @(#):$Id$
// Author: M.Gheata

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoParaEditor
\ingroup Geometry_builder

Editor for a TGeoPara.

\image html geom_para_pic.png

\image html geom_para_ed.png

*/

#include "TGeoParaEditor.h"
#include "TGeoTabManager.h"
#include "TGeoPara.h"
#include "TGeoManager.h"
#include "TVirtualGeoPainter.h"
#include "TVirtualPad.h"
#include "TView.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"

ClassImp(TGeoParaEditor);

enum ETGeoParaWid {
   kPARA_NAME, kPARA_X, kPARA_Y,  kPARA_Z, kPARA_ALPHA,
   kPARA_THETA, kPARA_PHI, kPARA_APPLY, kPARA_UNDO
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor for para editor

TGeoParaEditor::TGeoParaEditor(const TGWindow *p, Int_t width,
                                   Int_t height, UInt_t options, Pixel_t back)
   : TGeoGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fShape   = 0;
   fXi = fYi = fZi = fAlphai = fThetai = fPhii = 0.0;
   fNamei = "";
   fIsModified = kFALSE;
   fIsShapeEditable = kTRUE;

   // TextEntry for shape name
   MakeTitle("Name");
   fShapeName = new TGTextEntry(this, new TGTextBuffer(50), kPARA_NAME);
   fShapeName->Resize(135, fShapeName->GetDefaultHeight());
   fShapeName->SetToolTipText("Enter the parallelepiped name");
   fShapeName->Associate(this);
   AddFrame(fShapeName, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));

   TGTextEntry *nef;
   MakeTitle("Dimensions");
   // Number entry for dx
   TGCompositeFrame *f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "DX"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEDx = new TGNumberEntry(f1, 0., 5, kPARA_X);
   fEDx->SetNumAttr(TGNumberFormat::kNEAPositive);
   fEDx->Resize(100, fEDx->GetDefaultHeight());
   nef = (TGTextEntry*)fEDx->GetNumberEntry();
   nef->SetToolTipText("Enter the half-length in X");
   fEDx->Associate(this);
   f1->AddFrame(fEDx, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // Number entry for dy
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "DY"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEDy = new TGNumberEntry(f1, 0., 5, kPARA_Y);
   fEDy->SetNumAttr(TGNumberFormat::kNEAPositive);
   fEDy->Resize(100, fEDy->GetDefaultHeight());
   nef = (TGTextEntry*)fEDy->GetNumberEntry();
   nef->SetToolTipText("Enter the half-length in Y");
   fEDy->Associate(this);
   f1->AddFrame(fEDy, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // Number entry for dz
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Dz"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEDz = new TGNumberEntry(f1, 0., 5, kPARA_Z);
   fEDz->SetNumAttr(TGNumberFormat::kNEAPositive);
   fEDz->Resize(100, fEDz->GetDefaultHeight());
   nef = (TGTextEntry*)fEDz->GetNumberEntry();
   nef->SetToolTipText("Enter the half-length in Z");
   fEDz->Associate(this);
   f1->AddFrame(fEDz, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // Number entry for Alpha
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Alpha"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEAlpha = new TGNumberEntry(f1, 0., 5, kPARA_ALPHA);
   fEAlpha->Resize(100, fEAlpha->GetDefaultHeight());
   nef = (TGTextEntry*)fEAlpha->GetNumberEntry();
   nef->SetToolTipText("Enter the angle with respect to Y axis [deg]");
   fEAlpha->Associate(this);
   f1->AddFrame(fEAlpha, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // Number entry for Theta
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Theta"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fETheta = new TGNumberEntry(f1, 0., 5, kPARA_THETA);
   fETheta->SetNumAttr(TGNumberFormat::kNEAPositive);
   fETheta->Resize(100, fETheta->GetDefaultHeight());
   nef = (TGTextEntry*)fETheta->GetNumberEntry();
   nef->SetToolTipText("Enter the theta angle of the para axis [deg]");
   fETheta->Associate(this);
   f1->AddFrame(fETheta, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

    // Number entry for Phi
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Phi"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEPhi = new TGNumberEntry(f1, 0., 5, kPARA_PHI);
   fEPhi->SetNumAttr(TGNumberFormat::kNEAPositive);
   fEPhi->Resize(100, fEPhi->GetDefaultHeight());
   nef = (TGTextEntry*)fEPhi->GetNumberEntry();
   nef->SetToolTipText("Enter the phi angle of the para axis [deg]");
   fEPhi->Associate(this);
   f1->AddFrame(fEPhi, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

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
/// Destructor

TGeoParaEditor::~TGeoParaEditor()
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

void TGeoParaEditor::ConnectSignals2Slots()
{
   fApply->Connect("Clicked()", "TGeoParaEditor", this, "DoApply()");
   fUndo->Connect("Clicked()", "TGeoParaEditor", this, "DoUndo()");
   fShapeName->Connect("TextChanged(const char *)", "TGeoParaEditor", this, "DoModified()");
   fEDx->Connect("ValueSet(Long_t)", "TGeoParaEditor", this, "DoX()");
   fEDy->Connect("ValueSet(Long_t)", "TGeoParaEditor", this, "DoY()");
   fEDz->Connect("ValueSet(Long_t)", "TGeoParaEditor", this, "DoZ()");
   fEAlpha->Connect("ValueSet(Long_t)", "TGeoParaEditor", this, "DoAlpha()");
   fETheta->Connect("ValueSet(Long_t)", "TGeoParaEditor", this, "DoTheta()");
   fEPhi->Connect("ValueSet(Long_t)", "TGeoParaEditor", this, "DoPhi()");
   fEDx->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoParaEditor", this, "DoModified()");
   fEDy->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoParaEditor", this, "DoModified()");
   fEDz->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoParaEditor", this, "DoModified()");
   fEAlpha->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoParaEditor", this, "DoModified()");
   fETheta->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoParaEditor", this, "DoModified()");
   fEPhi->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoParaEditor", this, "DoModified()");
   fInit = kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Connect to the selected object.

void TGeoParaEditor::SetModel(TObject* obj)
{
   if (obj == 0 || (obj->IsA()!=TGeoPara::Class())) {
      SetActive(kFALSE);
      return;
   }
   fShape = (TGeoPara*)obj;
   fXi = fShape->GetX();
   fYi = fShape->GetY();
   fZi = fShape->GetZ();
   fAlphai = fShape->GetAlpha();
   fThetai = fShape->GetTheta();
   fPhii = fShape->GetPhi();
   const char *sname = fShape->GetName();
   if (!strcmp(sname, fShape->ClassName())) fShapeName->SetText("-no_name");
   else {
      fShapeName->SetText(sname);
      fNamei = sname;
   }
   fEDx->SetNumber(fXi);
   fEDy->SetNumber(fYi);
   fEDz->SetNumber(fZi);
   fEAlpha->SetNumber(fAlphai);
   fETheta->SetNumber(fThetai);
   fEPhi->SetNumber(fPhii);
   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);

   if (fInit) ConnectSignals2Slots();
   SetActive();
}

////////////////////////////////////////////////////////////////////////////////
/// Check if shape drawing is delayed.

Bool_t TGeoParaEditor::IsDelayed() const
{
   return (fDelayed->GetState() == kButtonDown);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for name.

void TGeoParaEditor::DoName()
{
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for applying current settings.

void TGeoParaEditor::DoApply()
{
   const char *name = fShapeName->GetText();
   if (strcmp(name,fShape->GetName())) fShape->SetName(name);
   Double_t dx = fEDx->GetNumber();
   Double_t dy = fEDy->GetNumber();
   Double_t dz = fEDz->GetNumber();
   Double_t alpha = fEAlpha->GetNumber();
   Double_t theta = fETheta->GetNumber();
   Double_t phi = fEPhi->GetNumber();
   Double_t param[6];
   param[0] = dx;
   param[1] = dy;
   param[2] = dz;
   param[3] = alpha;
   param[4] = theta;
   param[5] = phi;
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

void TGeoParaEditor::DoModified()
{
   fApply->SetEnabled();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for undoing last operation.

void TGeoParaEditor::DoUndo()
{
   fEDx->SetNumber(fXi);
   fEDy->SetNumber(fYi);
   fEDz->SetNumber(fZi);
   fEAlpha->SetNumber(fAlphai);
   fETheta->SetNumber(fThetai);
   fEPhi->SetNumber(fPhii);
   DoApply();
   fUndo->SetEnabled(kFALSE);
   fApply->SetEnabled(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for X.

void TGeoParaEditor::DoX()
{
   Double_t dx = fEDx->GetNumber();
   if (dx<=0) {
      dx = 0.1;
      fEDx->SetNumber(dx);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for Y.

void TGeoParaEditor::DoY()
{
   Double_t dy = fEDy->GetNumber();
   if (dy<=0) {
      dy = 0.1;
      fEDy->SetNumber(dy);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for Z.

void TGeoParaEditor::DoZ()
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
/// Slot for alpha.

void TGeoParaEditor::DoAlpha()
{
   Double_t alpha = fEAlpha->GetNumber();
   if (TMath::Abs(alpha)>=90) {
      alpha = 89.9*TMath::Sign(1.,alpha);
      fEAlpha->SetNumber(alpha);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for theta.

void TGeoParaEditor::DoTheta()
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

void TGeoParaEditor::DoPhi()
{
   Double_t phi = fEPhi->GetNumber();
   if (phi<0 || phi>360) {
      phi = 0;
      fEPhi->SetNumber(phi);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

