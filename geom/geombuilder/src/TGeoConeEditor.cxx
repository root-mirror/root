// @(#):$Id$
// Author: M.Gheata

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoConeEditor
\ingroup Geometry_builder

Editor for a TGeoCone.

\image html geom_cone_pic.png

\image html geom_cone_ed.png

*/


#include "TGeoConeEditor.h"
#include "TGeoTabManager.h"
#include "TGeoCone.h"
#include "TGeoManager.h"
#include "TVirtualGeoPainter.h"
#include "TVirtualPad.h"
#include "TView.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"
#include "TGDoubleSlider.h"

ClassImp(TGeoConeEditor);

enum ETGeoConeWid {
   kCONE_NAME, kCONE_RMIN1, kCONE_RMIN2, kCONE_RMAX1, kCONE_RMAX2, kCONE_Z,
   kCONE_APPLY, kCONE_UNDO
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor for volume editor

TGeoConeEditor::TGeoConeEditor(const TGWindow *p, Int_t width,
                                   Int_t height, UInt_t options, Pixel_t back)
   : TGeoGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fShape   = 0;
   fRmini1 = fRmaxi1 = fRmini2 = fRmaxi2 = fDzi = 0.0;
   fNamei = "";
   fIsModified = kFALSE;
   fIsShapeEditable = kTRUE;

   // TextEntry for shape name
   MakeTitle("Name");
   fShapeName = new TGTextEntry(this, new TGTextBuffer(50), kCONE_NAME);
   fShapeName->Resize(135, fShapeName->GetDefaultHeight());
   fShapeName->SetToolTipText("Enter the cone name");
   fShapeName->Associate(this);
   AddFrame(fShapeName, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));

   TGTextEntry *nef;
   MakeTitle("Cone dimensions");
   TGCompositeFrame *compxyz = new TGCompositeFrame(this, 118, 30, kVerticalFrame | kRaisedFrame);

   // Number entry for Rmin1
   TGCompositeFrame *f1 = new TGCompositeFrame(compxyz, 155, 30, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Rmin1"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fERmin1 = new TGNumberEntry(f1, 0., 5, kCONE_RMIN1);
   fERmin1->SetNumAttr(TGNumberFormat::kNEANonNegative);
   nef = (TGTextEntry*)fERmin1->GetNumberEntry();
   nef->SetToolTipText("Enter the inner radius");
   fERmin1->Associate(this);
   fERmin1->Resize(100, fERmin1->GetDefaultHeight());
   f1->AddFrame(fERmin1, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));

  // Number entry for Rmax1
   f1 = new TGCompositeFrame(compxyz, 155, 30, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Rmax1"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fERmax1 = new TGNumberEntry(f1, 0., 5, kCONE_RMAX1);
   fERmax1->SetNumAttr(TGNumberFormat::kNEANonNegative);
   nef = (TGTextEntry*)fERmax1->GetNumberEntry();
   nef->SetToolTipText("Enter the outer radius");
   fERmax1->Associate(this);
   fERmax1->Resize(100, fERmax1->GetDefaultHeight());
   f1->AddFrame(fERmax1, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));

   // Number entry for Rmin2
   f1 = new TGCompositeFrame(compxyz, 155, 30, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Rmin2"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fERmin2 = new TGNumberEntry(f1, 0., 5, kCONE_RMIN2);
   fERmin2->SetNumAttr(TGNumberFormat::kNEANonNegative);
   nef = (TGTextEntry*)fERmin2->GetNumberEntry();
   nef->SetToolTipText("Enter the inner radius");
   fERmin2->Associate(this);
   fERmin2->Resize(100, fERmin2->GetDefaultHeight());
   f1->AddFrame(fERmin2, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));

   // Number entry for Rmax2
   f1 = new TGCompositeFrame(compxyz, 155, 30, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Rmax2"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fERmax2 = new TGNumberEntry(f1, 0., 5, kCONE_RMAX2);
   fERmax2->SetNumAttr(TGNumberFormat::kNEANonNegative);
   nef = (TGTextEntry*)fERmax1->GetNumberEntry();
   nef->SetToolTipText("Enter the outer radius");
   fERmax2->Associate(this);
   fERmax2->Resize(100, fERmax2->GetDefaultHeight());
   f1->AddFrame(fERmax2, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));

   // Number entry for dz
   f1 = new TGCompositeFrame(compxyz, 155, 30, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "DZ"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEDz = new TGNumberEntry(f1, 0., 5, kCONE_Z);
   fEDz->SetNumAttr(TGNumberFormat::kNEAPositive);
   nef = (TGTextEntry*)fEDz->GetNumberEntry();
   nef->SetToolTipText("Enter the cone half-lenth in Z");
   fEDz->Associate(this);
   fEDz->Resize(100, fEDz->GetDefaultHeight());
   f1->AddFrame(fEDz, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));

   compxyz->Resize(150,30);
   AddFrame(compxyz, new TGLayoutHints(kLHintsLeft, 6, 6, 4, 4));

   // Delayed draw
   fDFrame = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth | kSunkenFrame);
   fDelayed = new TGCheckButton(fDFrame, "Delayed draw");
   fDFrame->AddFrame(fDelayed, new TGLayoutHints(kLHintsLeft , 2, 2, 2, 2));
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

TGeoConeEditor::~TGeoConeEditor()
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

void TGeoConeEditor::ConnectSignals2Slots()
{
   fApply->Connect("Clicked()", "TGeoConeEditor", this, "DoApply()");
   fUndo->Connect("Clicked()", "TGeoConeEditor", this, "DoUndo()");
   fShapeName->Connect("TextChanged(const char *)", "TGeoConeEditor", this, "DoModified()");
   fERmin1->Connect("ValueSet(Long_t)", "TGeoConeEditor", this, "DoRmin1()");
   fERmin2->Connect("ValueSet(Long_t)", "TGeoConeEditor", this, "DoRmin2()");
   fERmax1->Connect("ValueSet(Long_t)", "TGeoConeEditor", this, "DoRmax1()");
   fERmax2->Connect("ValueSet(Long_t)", "TGeoConeEditor", this, "DoRmax2()");
   fEDz->Connect("ValueSet(Long_t)", "TGeoConeEditor", this, "DoDz()");
   fERmin1->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoConeEditor", this, "DoRmin1()");
   fERmin2->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoConeEditor", this, "DoRmin2()");
   fERmax1->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoConeEditor", this, "DoRmax1()");
   fERmax2->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoConeEditor", this, "DoRmax2()");
   fEDz->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoConeEditor", this, "DoDz()");
   fInit = kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Connect to the selected object.

void TGeoConeEditor::SetModel(TObject* obj)
{
   if (obj == 0 || (obj->IsA()!=TGeoCone::Class())) {
      SetActive(kFALSE);
      return;
   }
   fShape = (TGeoCone*)obj;
   fRmini1 = fShape->GetRmin1();
   fRmini2 = fShape->GetRmin2();
   fRmaxi1 = fShape->GetRmax1();
   fRmaxi2 = fShape->GetRmax2();
   fDzi = fShape->GetDz();
   fNamei = fShape->GetName();
   fShapeName->SetText(fShape->GetName());
   fERmin1->SetNumber(fRmini1);
   fERmin2->SetNumber(fRmini2);
   fERmax1->SetNumber(fRmaxi1);
   fERmax2->SetNumber(fRmaxi2);
   fEDz->SetNumber(fDzi);
   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);

   if (fInit) ConnectSignals2Slots();
   SetActive();
}

////////////////////////////////////////////////////////////////////////////////
/// Check if shape drawing is delayed.

Bool_t TGeoConeEditor::IsDelayed() const
{
   return (fDelayed->GetState() == kButtonDown);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for name.

void TGeoConeEditor::DoName()
{
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
///Slot for applying current parameters.

void TGeoConeEditor::DoApply()
{
   fApply->SetEnabled(kFALSE);
   const char *name = fShapeName->GetText();
   if (strcmp(name,fShape->GetName())) fShape->SetName(name);
   Double_t rmin1 = fERmin1->GetNumber();
   Double_t rmin2 = fERmin2->GetNumber();
   Double_t rmax1 = fERmax1->GetNumber();
   Double_t rmax2 = fERmax2->GetNumber();
   Double_t dz = fEDz->GetNumber();
   if (rmin1<0 || rmin1>rmax1) return;
   if (rmin2<0 || rmin2>rmax2) return;
   if (dz<=0) return;
   if (rmin1==rmax1 && rmin2==rmax2) return;
   fShape->SetConeDimensions(dz, rmin1, rmax1, rmin2, rmax2);
   fShape->ComputeBBox();
   fUndo->SetEnabled();
   if (fPad) {
      if (gGeoManager && gGeoManager->GetPainter() && gGeoManager->GetPainter()->IsPaintingShape()) {
         fShape->Draw();
         fPad->GetView()->ShowAxis();
      } else Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
///Slot for modifying current parameters.

void TGeoConeEditor::DoModified()
{
   fApply->SetEnabled();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for undoing current operation.

void TGeoConeEditor::DoUndo()
{
   fERmin1->SetNumber(fRmini1);
   fERmin2->SetNumber(fRmini2);
   fERmax1->SetNumber(fRmaxi1);
   fERmax2->SetNumber(fRmaxi2);
   fEDz->SetNumber(fDzi);
   DoApply();
   fUndo->SetEnabled(kFALSE);
   fApply->SetEnabled(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for Rmin1

void TGeoConeEditor::DoRmin1()
{
   Double_t rmin1 = fERmin1->GetNumber();
   Double_t rmax1 = fERmax1->GetNumber();
   if (rmin1<0) {
      rmin1 = 0;
      fERmin1->SetNumber(rmin1);
   }
   if (rmin1>rmax1) {
      rmin1 = rmax1;
      fERmin1->SetNumber(rmin1);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for Rmax1

void TGeoConeEditor::DoRmax1()
{
   Double_t rmin1 = fERmin1->GetNumber();
   Double_t rmax1 = fERmax1->GetNumber();
   if (rmax1<rmin1) {
      rmax1 = rmin1;
      fERmax1->SetNumber(rmax1);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for Rmin2

void TGeoConeEditor::DoRmin2()
{
   Double_t rmin2 = fERmin2->GetNumber();
   Double_t rmax2 = fERmax2->GetNumber();
   if (rmin2<0) {
      rmin2 = 0;
      fERmin2->SetNumber(rmin2);
   }
   if (rmin2>rmax2) {
      rmin2 = rmax2;
      fERmin2->SetNumber(rmin2);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for  Rmax2

void TGeoConeEditor::DoRmax2()
{
   Double_t rmin2 = fERmin2->GetNumber();
   Double_t rmax2 = fERmax2->GetNumber();
   if (rmax2<rmin2) {
      rmax2 = rmin2;
      fERmax2->SetNumber(rmax2);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for Dz

void TGeoConeEditor::DoDz()
{
   Double_t dz = fEDz->GetNumber();
   if (dz<=0) {
      dz = 0.1;
      fEDz->SetNumber(dz);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

/** \class TGeoConeSegEditor
\ingroup Geometry_builder

Editor for a cone segment.

\image html geom_cons_pic.png

\image html geom_cons_ed.png

*/

ClassImp(TGeoConeSegEditor);

enum ETGeoConeSegWid {
   kCONESEG_PHI1, kCONESEG_PHI2, kCONESEG_PHI
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor for cone segment editor

TGeoConeSegEditor::TGeoConeSegEditor(const TGWindow *p, Int_t width,
                                   Int_t height, UInt_t options, Pixel_t back)
                 : TGeoConeEditor(p, width, height, options | kVerticalFrame, back)
{
   fLock = kFALSE;
   MakeTitle("Phi range");
   TGTextEntry *nef;
   TGCompositeFrame *compxyz = new TGCompositeFrame(this, 155, 110, kHorizontalFrame | kFixedWidth | kFixedHeight | kRaisedFrame);
   // Vertical slider
   fSPhi = new TGDoubleVSlider(compxyz,100);
   fSPhi->SetRange(0.,720.);
   fSPhi->Resize(fSPhi->GetDefaultWidth(), 100);
   compxyz->AddFrame(fSPhi, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));
   TGCompositeFrame *f1 = new TGCompositeFrame(compxyz, 135, 100, kVerticalFrame | kFixedHeight);
   f1->AddFrame(new TGLabel(f1, "Phi min."), new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 0, 6, 0));
   fEPhi1 = new TGNumberEntry(f1, 0., 5, kCONESEG_PHI1);
   fEPhi1->Resize(100, fEPhi1->GetDefaultHeight());
   fEPhi1->SetNumAttr(TGNumberFormat::kNEANonNegative);
   nef = (TGTextEntry*)fEPhi1->GetNumberEntry();
   nef->SetToolTipText("Enter the phi1 value");
   fEPhi1->Associate(this);
   f1->AddFrame(fEPhi1, new TGLayoutHints(kLHintsTop | kLHintsRight, 2, 2, 2, 2));

   fEPhi2 = new TGNumberEntry(f1, 0., 5, kCONESEG_PHI2);
   fEPhi2->Resize(100, fEPhi2->GetDefaultHeight());
   fEPhi2->SetNumAttr(TGNumberFormat::kNEANonNegative);
   nef = (TGTextEntry*)fEPhi2->GetNumberEntry();
   nef->SetToolTipText("Enter the phi2 value");
   fEPhi2->Associate(this);
   f1->AddFrame(fEPhi2, new TGLayoutHints(kLHintsBottom | kLHintsRight, 2, 2, 2, 2));
   f1->AddFrame(new TGLabel(f1, "Phi max."), new TGLayoutHints(kLHintsBottom, 0, 0, 6, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

//   compxyz->Resize(150,150);
   AddFrame(compxyz, new TGLayoutHints(kLHintsLeft, 6, 6, 4, 4));
   TGeoTabManager::MoveFrame(fDFrame, this);
   TGeoTabManager::MoveFrame(fBFrame, this);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoConeSegEditor::~TGeoConeSegEditor()
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

void TGeoConeSegEditor::ConnectSignals2Slots()
{
   TGeoConeEditor::ConnectSignals2Slots();
   Disconnect(fApply, "Clicked()",(TGeoConeEditor*)this, "DoApply()");
   Disconnect(fUndo, "Clicked()",(TGeoConeEditor*)this, "DoUndo()");
   fApply->Connect("Clicked()", "TGeoConeSegEditor", this, "DoApply()");
   fUndo->Connect("Clicked()", "TGeoConeSegEditor", this, "DoUndo()");
   fEPhi1->Connect("ValueSet(Long_t)", "TGeoConeSegEditor", this, "DoPhi1()");
   fEPhi2->Connect("ValueSet(Long_t)", "TGeoConeSegEditor", this, "DoPhi2()");
//   fEPhi1->GetNumberEntry()->Connect("TextChanged(const char *)","TGeoConeSegEditor", this, "DoPhi1()");
//   fEPhi2->GetNumberEntry()->Connect("TextChanged(const char *)","TGeoConeSegEditor", this, "DoPhi2()");
   fSPhi->Connect("PositionChanged()","TGeoConeSegEditor", this, "DoPhi()");
}

////////////////////////////////////////////////////////////////////////////////
/// Connect to the selected object.

void TGeoConeSegEditor::SetModel(TObject* obj)
{
   if (obj == 0 || (obj->IsA()!=TGeoConeSeg::Class())) {
      SetActive(kFALSE);
      return;
   }
   fShape = (TGeoCone*)obj;
   fRmini1 = fShape->GetRmin1();
   fRmaxi1 = fShape->GetRmax1();
   fRmini2 = fShape->GetRmin2();
   fRmaxi2 = fShape->GetRmax2();
   fDzi = fShape->GetDz();
   fNamei = fShape->GetName();
   fPmini = ((TGeoConeSeg*)fShape)->GetPhi1();
   fPmaxi = ((TGeoConeSeg*)fShape)->GetPhi2();
   fShapeName->SetText(fShape->GetName());
   fEPhi1->SetNumber(fPmini);
   fEPhi2->SetNumber(fPmaxi);
   fSPhi->SetPosition(fPmini,fPmaxi);
   fERmin1->SetNumber(fRmini1);
   fERmax1->SetNumber(fRmaxi1);
   fERmin2->SetNumber(fRmini2);
   fERmax2->SetNumber(fRmaxi2);
   fEDz->SetNumber(fDzi);
   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);

   if (fInit) ConnectSignals2Slots();
   SetActive();
}

////////////////////////////////////////////////////////////////////////////////
///Slot for Phi1

void TGeoConeSegEditor::DoPhi1()
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
/// Slot for Phi2

void TGeoConeSegEditor::DoPhi2()
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
/// Slot for Phi

void TGeoConeSegEditor::DoPhi()
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
/// Slot for applying current parameters.

void TGeoConeSegEditor::DoApply()
{
   fApply->SetEnabled(kFALSE);
   const char *name = fShapeName->GetText();
   if (strcmp(name,fShape->GetName())) fShape->SetName(name);
   Double_t rmin1 = fERmin1->GetNumber();
   Double_t rmax1 = fERmax1->GetNumber();
   if (rmin1<0 || rmax1<rmin1) return;
   Double_t rmin2 = fERmin2->GetNumber();
   Double_t rmax2 = fERmax2->GetNumber();
   if (rmin2<0 || rmax2<rmin2) return;
   Double_t dz = fEDz->GetNumber();
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
   ((TGeoConeSeg*)fShape)->SetConsDimensions(dz, rmin1, rmax1, rmin2,rmax2, phi1, phi2);
   fShape->ComputeBBox();
   fUndo->SetEnabled();
   if (fPad) {
      if (gGeoManager && gGeoManager->GetPainter() && gGeoManager->GetPainter()->IsPaintingShape()) {
         fShape->Draw();
         fPad->GetView()->ShowAxis();
      } else Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for undoing last operation.

void TGeoConeSegEditor::DoUndo()
{
   fERmin1->SetNumber(fRmini1);
   fERmin2->SetNumber(fRmini2);
   fERmax1->SetNumber(fRmaxi1);
   fERmax2->SetNumber(fRmaxi2);
   fEDz->SetNumber(fDzi);
   fEPhi1->SetNumber(fPmini);
   fEPhi2->SetNumber(fPmaxi);
   fSPhi->SetPosition(fPmini,fPmaxi);
   DoApply();
   fUndo->SetEnabled(kFALSE);
   fApply->SetEnabled(kFALSE);
}




