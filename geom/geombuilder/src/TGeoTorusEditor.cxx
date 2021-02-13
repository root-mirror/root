// @(#):$Id$
// Author: M.Gheata

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoTorusEditor
\ingroup Geometry_builder

Editor for a TGeoTorus.

\image html geom_torus_pic.png

\image html geom_torus_ed.png

*/

#include "TGeoTorusEditor.h"
#include "TGeoTabManager.h"
#include "TGeoTorus.h"
#include "TGeoManager.h"
#include "TVirtualGeoPainter.h"
#include "TVirtualPad.h"
#include "TView.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"

ClassImp(TGeoTorusEditor);

enum ETGeoTorusWid {
   kTORUS_NAME, kTORUS_R, kTORUS_RMIN,  kTORUS_RMAX, kTORUS_PHI1,
   kTORUS_DPHI, kTORUS_APPLY, kTORUS_UNDO
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor for torus editor

TGeoTorusEditor::TGeoTorusEditor(const TGWindow *p, Int_t width,
                                   Int_t height, UInt_t options, Pixel_t back)
   : TGeoGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fShape   = 0;
   fRi = fRmini = fRmaxi = fPhi1i = fDphii = 0.0;
   fNamei = "";
   fIsModified = kFALSE;
   fIsShapeEditable = kTRUE;

   // TextEntry for shape name
   MakeTitle("Name");
   fShapeName = new TGTextEntry(this, new TGTextBuffer(50), kTORUS_NAME);
   fShapeName->Resize(135, fShapeName->GetDefaultHeight());
   fShapeName->SetToolTipText("Enter the torus name");
   fShapeName->Associate(this);
   AddFrame(fShapeName, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));

   TGTextEntry *nef;
   MakeTitle("Dimensions");
   // Number entry for R.
   TGCompositeFrame *f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "R"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fER = new TGNumberEntry(f1, 0., 5, kTORUS_R);
   fER->SetNumAttr(TGNumberFormat::kNEAPositive);
   fER->Resize(100, fER->GetDefaultHeight());
   nef = (TGTextEntry*)fER->GetNumberEntry();
   nef->SetToolTipText("Enter the axial radius R");
   fER->Associate(this);
   f1->AddFrame(fER, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // Number entry for rmin.
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Rmin"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fERmin = new TGNumberEntry(f1, 0., 5, kTORUS_RMIN);
   fERmin->SetNumAttr(TGNumberFormat::kNEAPositive);
   fERmin->Resize(100, fERmin->GetDefaultHeight());
   nef = (TGTextEntry*)fERmin->GetNumberEntry();
   nef->SetToolTipText("Enter the inner radius Rmin");
   fERmin->Associate(this);
   f1->AddFrame(fERmin, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // Number entry for rmax
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Rmax"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fERmax = new TGNumberEntry(f1, 0., 5, kTORUS_RMAX);
   fERmax->SetNumAttr(TGNumberFormat::kNEAPositive);
   fERmax->Resize(100, fERmax->GetDefaultHeight());
   nef = (TGTextEntry*)fERmax->GetNumberEntry();
   nef->SetToolTipText("Enter the outer radius Rmax");
   fERmax->Associate(this);
   f1->AddFrame(fERmax, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // Number entry for Phi1
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Phi1"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEPhi1 = new TGNumberEntry(f1, 0., 5, kTORUS_PHI1);
   fEPhi1->SetNumAttr(TGNumberFormat::kNEANonNegative);
   fEPhi1->Resize(100, fEPhi1->GetDefaultHeight());
   nef = (TGTextEntry*)fEPhi1->GetNumberEntry();
   nef->SetToolTipText("Enter the starting phi angle[deg]");
   fEPhi1->Associate(this);
   f1->AddFrame(fEPhi1, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // Number entry for Dphi
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Dphi"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEDphi = new TGNumberEntry(f1, 0., 5, kTORUS_DPHI);
   fEDphi->SetNumAttr(TGNumberFormat::kNEAPositive);
   fEDphi->Resize(100, fEDphi->GetDefaultHeight());
   nef = (TGTextEntry*)fEDphi->GetNumberEntry();
   nef->SetToolTipText("Enter the extent phi Dphi [deg]");
   fEDphi->Associate(this);
   f1->AddFrame(fEDphi, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
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

TGeoTorusEditor::~TGeoTorusEditor()
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

void TGeoTorusEditor::ConnectSignals2Slots()
{
   fApply->Connect("Clicked()", "TGeoTorusEditor", this, "DoApply()");
   fUndo->Connect("Clicked()", "TGeoTorusEditor", this, "DoUndo()");
   fShapeName->Connect("TextChanged(const char *)", "TGeoTorusEditor", this, "DoModified()");
   fER->Connect("ValueSet(Long_t)", "TGeoTorusEditor", this, "DoR()");
   fERmin->Connect("ValueSet(Long_t)", "TGeoTorusEditor", this, "DoRmin()");
   fERmax->Connect("ValueSet(Long_t)", "TGeoTorusEditor", this, "DoRmax()");
   fEPhi1->Connect("ValueSet(Long_t)", "TGeoTorusEditor", this, "DoPhi1()");
   fEDphi->Connect("ValueSet(Long_t)", "TGeoTorusEditor", this, "DoDphi()");
   fER->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTorusEditor", this, "DoModified()");
   fERmin->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTorusEditor", this, "DoModified()");
   fERmax->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTorusEditor", this, "DoModified()");
   fEPhi1->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTorusEditor", this, "DoModified()");
   fEDphi->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTorusEditor", this, "DoModified()");
   fInit = kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Connect to the selected object.

void TGeoTorusEditor::SetModel(TObject* obj)
{
   if (obj == 0 || (obj->IsA()!=TGeoTorus::Class())) {
      SetActive(kFALSE);
      return;
   }
   fShape = (TGeoTorus*)obj;
   fRi = fShape->GetR();
   fRmini = fShape->GetRmin();
   fRmaxi = fShape->GetRmax();
   fPhi1i = fShape->GetPhi1();
   fDphii = fShape->GetDphi();
   const char *sname = fShape->GetName();
   if (!strcmp(sname, fShape->ClassName())) fShapeName->SetText("-no_name");
   else {
      fShapeName->SetText(sname);
      fNamei = sname;
   }
   fER->SetNumber(fRi);
   fERmin->SetNumber(fRmini);
   fERmax->SetNumber(fRmaxi);
   fEPhi1->SetNumber(fPhi1i);
   fEDphi->SetNumber(fDphii);
   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);

   if (fInit) ConnectSignals2Slots();
   SetActive();
}

////////////////////////////////////////////////////////////////////////////////
/// Check if shape drawing is delayed.

Bool_t TGeoTorusEditor::IsDelayed() const
{
   return (fDelayed->GetState() == kButtonDown);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for name.

void TGeoTorusEditor::DoName()
{
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for applying current settings.

void TGeoTorusEditor::DoApply()
{
   const char *name = fShapeName->GetText();
   if (strcmp(name,fShape->GetName())) fShape->SetName(name);
   Double_t r = fER->GetNumber();
   Double_t rmax = fERmax->GetNumber();
   Double_t rmin = fERmin->GetNumber();
   Double_t phi = fEPhi1->GetNumber();
   Double_t dphi = fEDphi->GetNumber();
   Double_t param[5];
   param[0] = r;
   param[1] = rmin;
   param[2] = rmax;
   param[3] = phi;
   param[4] = dphi;
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

void TGeoTorusEditor::DoModified()
{
   fApply->SetEnabled();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for undoing last operation.

void TGeoTorusEditor::DoUndo()
{
   fER->SetNumber(fRi);
   fERmin->SetNumber(fRmini);
   fERmax->SetNumber(fRmaxi);
   fEPhi1->SetNumber(fPhi1i);
   fEDphi->SetNumber(fDphii);
   DoApply();
   fUndo->SetEnabled(kFALSE);
   fApply->SetEnabled(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for R.

void TGeoTorusEditor::DoR()
{
   Double_t r = fER->GetNumber();
   Double_t rmax = fERmax->GetNumber();
   if (r<rmax) {
      r = rmax;
      fER->SetNumber(r);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for Rmin.

void TGeoTorusEditor::DoRmin()
{
   Double_t rmin = fERmin->GetNumber();
   Double_t rmax = fERmax->GetNumber();
   if (rmin>rmax) {
      rmin = rmax-0.1;
      fERmin->SetNumber(rmin);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for Rmax.

void TGeoTorusEditor::DoRmax()
{
   Double_t r = fER->GetNumber();
   Double_t rmin = fERmin->GetNumber();
   Double_t rmax = fERmax->GetNumber();
   if (rmax<=rmin) {
      rmax = rmin+0.1;
      fERmax->SetNumber(rmax);
   }
   if (rmax>r) {
      rmax = r;
      fERmax->SetNumber(rmax);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for phi.

void TGeoTorusEditor::DoPhi1()
{
   Double_t phi = fEPhi1->GetNumber();
   if (phi<0 || phi>360) {
      phi = 0;
      fEPhi1->SetNumber(phi);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for Dphi.

void TGeoTorusEditor::DoDphi()
{
   Double_t dphi = fEDphi->GetNumber();
   if (dphi<=0 || dphi>360) {
      dphi = 1;
      fEDphi->SetNumber(dphi);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

