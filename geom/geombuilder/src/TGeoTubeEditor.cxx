// @(#):$Id$
// Author: M.Gheata

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoTubeEditor
\ingroup Geometry_builder

Editor for a TGeoTube.

\image html geom_tube_pic.png

\image html geom_tube_ed.png

*/

#include "TGeoTubeEditor.h"
#include "TGeoTabManager.h"
#include "TGeoTube.h"
#include "TGeoManager.h"
#include "TVirtualGeoPainter.h"
#include "TVirtualPad.h"
#include "TView.h"
#include "TMath.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"
#include "TGDoubleSlider.h"

ClassImp(TGeoTubeEditor);

enum ETGeoTubeWid {
   kTUBE_NAME, kTUBE_RMIN, kTUBE_RMAX, kTUBE_Z,
   kTUBE_APPLY, kTUBE_UNDO
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor for tube editor

TGeoTubeEditor::TGeoTubeEditor(const TGWindow *p, Int_t width,
                                   Int_t height, UInt_t options, Pixel_t back)
   : TGeoGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fShape   = 0;
   fRmini = fRmaxi = fDzi = 0.0;
   fNamei = "";
   fIsModified = kFALSE;
   fIsShapeEditable = kTRUE;

   // TextEntry for shape name
   MakeTitle("Name");
   fShapeName = new TGTextEntry(this, new TGTextBuffer(50), kTUBE_NAME);
   fShapeName->Resize(135, fShapeName->GetDefaultHeight());
   fShapeName->SetToolTipText("Enter the box name");
   fShapeName->Associate(this);
   AddFrame(fShapeName, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));

   TGTextEntry *nef;
   MakeTitle("Tube dimensions");
   TGCompositeFrame *compxyz = new TGCompositeFrame(this, 118, 30, kVerticalFrame | kRaisedFrame);
   // Number entry for rmin
   TGCompositeFrame *f1 = new TGCompositeFrame(compxyz, 155, 30, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Rmin"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fERmin = new TGNumberEntry(f1, 0., 5, kTUBE_RMIN);
   fERmin->SetNumAttr(TGNumberFormat::kNEANonNegative);
   nef = (TGTextEntry*)fERmin->GetNumberEntry();
   nef->SetToolTipText("Enter the inner radius");
   fERmin->Associate(this);
   fERmin->Resize(100,fERmin->GetDefaultHeight());
   f1->AddFrame(fERmin, new TGLayoutHints(kLHintsRight , 2, 2, 2, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));

   // Number entry for Rmax
   f1 = new TGCompositeFrame(compxyz, 155, 30, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Rmax"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fERmax = new TGNumberEntry(f1, 0., 5, kTUBE_RMAX);
   fERmax->SetNumAttr(TGNumberFormat::kNEANonNegative);
   nef = (TGTextEntry*)fERmax->GetNumberEntry();
   nef->SetToolTipText("Enter the outer radius");
   fERmax->Associate(this);
   fERmax->Resize(100,fERmax->GetDefaultHeight());
   f1->AddFrame(fERmax, new TGLayoutHints(kLHintsRight , 2, 2, 2, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));

   // Number entry for dz
   f1 = new TGCompositeFrame(compxyz, 155, 30, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "DZ"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEDz = new TGNumberEntry(f1, 0., 5, kTUBE_Z);
   fEDz->SetNumAttr(TGNumberFormat::kNEAPositive);
   nef = (TGTextEntry*)fEDz->GetNumberEntry();
   nef->SetToolTipText("Enter the tube half-lenth in Z");
   fEDz->Associate(this);
   fEDz->Resize(100,fEDz->GetDefaultHeight());
   f1->AddFrame(fEDz, new TGLayoutHints(kLHintsRight , 2, 2, 2, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));

//   compxyz->Resize(150,30);
   AddFrame(compxyz, new TGLayoutHints(kLHintsLeft, 6, 6, 4, 4));

   // Delayed draw
   fDFrame = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth | kSunkenFrame);
   fDelayed = new TGCheckButton(fDFrame, "Delayed draw");
   fDFrame->AddFrame(fDelayed, new TGLayoutHints(kLHintsLeft , 2, 2, 4, 4));
   AddFrame(fDFrame,  new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

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

TGeoTubeEditor::~TGeoTubeEditor()
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

void TGeoTubeEditor::ConnectSignals2Slots()
{
   fApply->Connect("Clicked()", "TGeoTubeEditor", this, "DoApply()");
   fUndo->Connect("Clicked()", "TGeoTubeEditor", this, "DoUndo()");
   fShapeName->Connect("TextChanged(const char *)", "TGeoTubeEditor", this, "DoModified()");
   fERmin->Connect("ValueSet(Long_t)", "TGeoTubeEditor", this, "DoRmin()");
   fERmax->Connect("ValueSet(Long_t)", "TGeoTubeEditor", this, "DoRmax()");
   fEDz->Connect("ValueSet(Long_t)", "TGeoTubeEditor", this, "DoDz()");
   fERmin->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTubeEditor", this, "DoRmin()");
   fERmax->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTubeEditor", this, "DoRmax()");
   fEDz->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTubeEditor", this, "DoDz()");
   fInit = kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Connect to the selected object.

void TGeoTubeEditor::SetModel(TObject* obj)
{
   if (obj == 0 || (obj->IsA()!=TGeoTube::Class())) {
      SetActive(kFALSE);
      return;
   }
   fShape = (TGeoTube*)obj;
   fRmini = fShape->GetRmin();
   fRmaxi = fShape->GetRmax();
   fDzi = fShape->GetDz();
   fNamei = fShape->GetName();
   fShapeName->SetText(fShape->GetName());
   fERmin->SetNumber(fRmini);
   fERmax->SetNumber(fRmaxi);
   fEDz->SetNumber(fDzi);
   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);

   if (fInit) ConnectSignals2Slots();
   SetActive();
}

////////////////////////////////////////////////////////////////////////////////
/// Check if shape drawing is delayed.

Bool_t TGeoTubeEditor::IsDelayed() const
{
   return (fDelayed->GetState() == kButtonDown);
}

////////////////////////////////////////////////////////////////////////////////
/// Perform name change.

void TGeoTubeEditor::DoName()
{
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for applying modifications.

void TGeoTubeEditor::DoApply()
{
   const char *name = fShapeName->GetText();
   if (strcmp(name,fShape->GetName())) fShape->SetName(name);
   Double_t rmin = fERmin->GetNumber();
   Double_t rmax = fERmax->GetNumber();
   Double_t dz = fEDz->GetNumber();
   fShape->SetTubeDimensions(rmin, rmax, dz);
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

void TGeoTubeEditor::DoModified()
{
   fApply->SetEnabled();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for undoing last operation.

void TGeoTubeEditor::DoUndo()
{
   fERmin->SetNumber(fRmini);
   fERmax->SetNumber(fRmaxi);
   fEDz->SetNumber(fDzi);
   DoApply();
   fUndo->SetEnabled(kFALSE);
   fApply->SetEnabled(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for rmin.

void TGeoTubeEditor::DoRmin()
{
   Double_t rmin = fERmin->GetNumber();
   Double_t rmax = fERmax->GetNumber();
   if (rmax<rmin+1.e-10) {
      rmin = rmax - 0.1;
      fERmin->SetNumber(rmin);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for rmax.

void TGeoTubeEditor::DoRmax()
{
   Double_t rmin = fERmin->GetNumber();
   Double_t rmax = fERmax->GetNumber();
   if (rmax <= 0.) {
      rmax = 0.1;
      fERmax->SetNumber(rmax);
   }
   if (rmax<rmin+1.e-10) {
      rmax = rmin + 0.1;
      fERmax->SetNumber(rmax);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for dz.

void TGeoTubeEditor::DoDz()
{
   Double_t dz = fEDz->GetNumber();
   if (dz<=0) {
      dz = 0.1;
      fEDz->SetNumber(dz);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

/** \class TGeoTubeSegEditor
\ingroup Geometry_builder

Editor for a TGeoTubeSeg.

\image html geom_tubs_pic.png

\image html geom_tubs_ed.png

*/

ClassImp(TGeoTubeSegEditor);

enum ETGeoTubeSegWid {
   kTUBESEG_PHI1, kTUBESEG_PHI2, kTUBESEG_PHI
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor for tube segment editor

TGeoTubeSegEditor::TGeoTubeSegEditor(const TGWindow *p, Int_t width,
                                     Int_t height, UInt_t options, Pixel_t back)
  : TGeoTubeEditor(p, width, height, options | kVerticalFrame, back)
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
   fEPhi1 = new TGNumberEntry(f1, 0., 5, kTUBESEG_PHI1);
   fEPhi1->Resize(100, fEPhi1->GetDefaultHeight());
   fEPhi1->SetNumAttr(TGNumberFormat::kNEANonNegative);
   nef = (TGTextEntry*)fEPhi1->GetNumberEntry();
   nef->SetToolTipText("Enter the phi1 value");
   fEPhi1->Associate(this);
   f1->AddFrame(fEPhi1, new TGLayoutHints(kLHintsTop | kLHintsRight, 2, 2, 2, 2));

   fEPhi2 = new TGNumberEntry(f1, 0., 5, kTUBESEG_PHI2);
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

TGeoTubeSegEditor::~TGeoTubeSegEditor()
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

void TGeoTubeSegEditor::ConnectSignals2Slots()
{
   TGeoTubeEditor::ConnectSignals2Slots();
   Disconnect(fApply, "Clicked()",(TGeoTubeEditor*)this, "DoApply()");
   Disconnect(fUndo, "Clicked()",(TGeoTubeEditor*)this, "DoUndo()");
   fApply->Connect("Clicked()", "TGeoTubeSegEditor", this, "DoApply()");
   fUndo->Connect("Clicked()", "TGeoTubeSegEditor", this, "DoUndo()");
   fEPhi1->Connect("ValueSet(Long_t)", "TGeoTubeSegEditor", this, "DoPhi1()");
   fEPhi2->Connect("ValueSet(Long_t)", "TGeoTubeSegEditor", this, "DoPhi2()");
//   fEPhi1->GetNumberEntry()->Connect("TextChanged(const char *)","TGeoTubeSegEditor", this, "DoPhi1()");
//   fEPhi2->GetNumberEntry()->Connect("TextChanged(const char *)","TGeoTubeSegEditor", this, "DoPhi2()");
   fSPhi->Connect("PositionChanged()","TGeoTubeSegEditor", this, "DoPhi()");
}

////////////////////////////////////////////////////////////////////////////////
/// Connect to the selected object.

void TGeoTubeSegEditor::SetModel(TObject* obj)
{
   if (obj == 0 || (obj->IsA()!=TGeoTubeSeg::Class())) {
      SetActive(kFALSE);
      return;
   }
   fShape = (TGeoTube*)obj;
   fRmini = fShape->GetRmin();
   fRmaxi = fShape->GetRmax();
   fDzi = fShape->GetDz();
   fNamei = fShape->GetName();
   fPmini = ((TGeoTubeSeg*)fShape)->GetPhi1();
   fPmaxi = ((TGeoTubeSeg*)fShape)->GetPhi2();
   fShapeName->SetText(fShape->GetName());
   fEPhi1->SetNumber(fPmini);
   fEPhi2->SetNumber(fPmaxi);
   fSPhi->SetPosition(fPmini,fPmaxi);
   fERmin->SetNumber(fRmini);
   fERmax->SetNumber(fRmaxi);
   fEDz->SetNumber(fDzi);
   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);

   if (fInit) ConnectSignals2Slots();
   SetActive();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for phi1.

void TGeoTubeSegEditor::DoPhi1()
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

void TGeoTubeSegEditor::DoPhi2()
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

void TGeoTubeSegEditor::DoPhi()
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
/// Slot for applying modifications.

void TGeoTubeSegEditor::DoApply()
{
   fApply->SetEnabled(kFALSE);
   const char *name = fShapeName->GetText();
   if (strcmp(name,fShape->GetName())) fShape->SetName(name);
   Double_t rmin = fERmin->GetNumber();
   Double_t rmax = fERmax->GetNumber();
   if (rmin<0 || rmax<rmin) return;
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
   ((TGeoTubeSeg*)fShape)->SetTubsDimensions(rmin, rmax, dz, phi1, phi2);
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

void TGeoTubeSegEditor::DoUndo()
{
   fERmin->SetNumber(fRmini);
   fERmax->SetNumber(fRmaxi);
   fEDz->SetNumber(fDzi);
   fEPhi1->SetNumber(fPmini);
   fEPhi2->SetNumber(fPmaxi);
   fSPhi->SetPosition(fPmini,fPmaxi);
   DoApply();
   fUndo->SetEnabled(kFALSE);
   fApply->SetEnabled(kFALSE);
}

/** \class TGeoCtubEditor
\ingroup Geometry_builder

Editor for a TGeoCtub.

\image html geom_ctub_pic.png

\image html geom_ctub_ed.png

*/

ClassImp(TGeoCtubEditor);

enum ETGeoCtubSegWid {
   kCTUB_THLO, kCTUB_PHLO, kCTUB_THHI, kCTUB_PHHI
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor for cut tube editor

TGeoCtubEditor::TGeoCtubEditor(const TGWindow *p, Int_t width,
                               Int_t height, UInt_t options, Pixel_t back)
  : TGeoTubeSegEditor(p, width, height, options, back)
{
   MakeTitle("Theta/phi low");
   TGTextEntry *nef;
   // Number entry for theta/phi of the lower normal
   TGCompositeFrame *compxyz = new TGCompositeFrame(this, 118, 30, kVerticalFrame | kRaisedFrame);
   TGCompositeFrame *f1 = new TGCompositeFrame(compxyz, 155, 30, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "TH_LO"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEThlo = new TGNumberEntry(f1, 0., 5, kCTUB_THLO);
   fEThlo->SetNumAttr(TGNumberFormat::kNEANonNegative);
   nef = (TGTextEntry*)fEThlo->GetNumberEntry();
   nef->SetToolTipText("Enter the theta angle of the lower plane normal");
   fEThlo->Associate(this);
   fEThlo->Connect("ValueSet(Long_t)", "TGeoCtubEditor", this, "DoThlo()");
   nef->Connect("TextChanged(const char *)","TGeoCtubEditor", this, "DoModified()");
   fEThlo->Resize(100,fEThlo->GetDefaultHeight());
   f1->AddFrame(fEThlo, new TGLayoutHints(kLHintsRight , 2, 2, 2, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));

   f1 = new TGCompositeFrame(compxyz, 155, 30, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "PH_LO"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEPhlo = new TGNumberEntry(f1, 0., 5, kCTUB_PHLO);
   fEPhlo->SetNumAttr(TGNumberFormat::kNEANonNegative);
   nef = (TGTextEntry*)fEPhlo->GetNumberEntry();
   nef->SetToolTipText("Enter the phi angle of the lower plane normal");
   fEPhlo->Associate(this);
   fEPhlo->Connect("ValueSet(Long_t)", "TGeoCtubEditor", this, "DoPhlo()");
   nef->Connect("TextChanged(const char *)","TGeoCtubEditor", this, "DoModified()");
   fEPhlo->Resize(100,fEPhlo->GetDefaultHeight());
   f1->AddFrame(fEPhlo, new TGLayoutHints(kLHintsRight , 2, 2, 2, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));
   AddFrame(compxyz, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));

   // Number entry for theta/phi of the lower normal
   MakeTitle("Theta/phi high");
   compxyz = new TGCompositeFrame(this, 118, 30, kVerticalFrame | kRaisedFrame);
   f1 = new TGCompositeFrame(compxyz, 155, 30, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "TH_HI"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEThhi = new TGNumberEntry(f1, 0., 5, kCTUB_THHI);
   fEThhi->SetNumAttr(TGNumberFormat::kNEANonNegative);
   nef = (TGTextEntry*)fEThhi->GetNumberEntry();
   nef->SetToolTipText("Enter the theta angle of the upper plane normal");
   fEThhi->Associate(this);
   fEThhi->Connect("ValueSet(Long_t)", "TGeoCtubEditor", this, "DoThhi()");
   nef->Connect("TextChanged(const char *)","TGeoCtubEditor", this, "DoModified()");
   fEThhi->Resize(100,fEThhi->GetDefaultHeight());
   f1->AddFrame(fEThhi, new TGLayoutHints(kLHintsRight , 2, 2, 2, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));

   f1 = new TGCompositeFrame(compxyz, 155, 30, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "PH_HI"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEPhhi = new TGNumberEntry(f1, 0., 5, kCTUB_PHHI);
   fEPhhi->SetNumAttr(TGNumberFormat::kNEANonNegative);
   nef = (TGTextEntry*)fEPhhi->GetNumberEntry();
   nef->SetToolTipText("Enter the phi angle of the upper plane normal");
   fEPhhi->Associate(this);
   fEPhhi->Connect("ValueSet(Long_t)", "TGeoCtubEditor", this, "DoPhhi()");
   nef->Connect("TextChanged(const char *)","TGeoCtubEditor", this, "DoModified()");
   fEPhhi->Resize(100,fEPhhi->GetDefaultHeight());
   f1->AddFrame(fEPhhi, new TGLayoutHints(kLHintsRight , 2, 2, 4, 4));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));
   AddFrame(compxyz, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
   TGeoTabManager::MoveFrame(fDFrame, this);
   TGeoTabManager::MoveFrame(fBFrame, this);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoCtubEditor::~TGeoCtubEditor()
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
/// Connect to the selected object.

void TGeoCtubEditor::SetModel(TObject* obj)
{
   if (obj == 0 || (obj->IsA()!=TGeoCtub::Class())) {
      SetActive(kFALSE);
      return;
   }
   fShape = (TGeoTube*)obj;
   fRmini = fShape->GetRmin();
   fRmaxi = fShape->GetRmax();
   fDzi = fShape->GetDz();
   fNamei = fShape->GetName();
   fPmini = ((TGeoTubeSeg*)fShape)->GetPhi1();
   fPmaxi = ((TGeoTubeSeg*)fShape)->GetPhi2();
   const Double_t *nlo = ((TGeoCtub*)fShape)->GetNlow();
   const Double_t *nhi = ((TGeoCtub*)fShape)->GetNhigh();
   fThlo = TMath::RadToDeg() * TMath::ACos(nlo[2]);
   fPhlo = TMath::RadToDeg() * TMath::ATan2(nlo[1], nlo[0]);
   fThhi = TMath::RadToDeg() * TMath::ACos(nhi[2]);
   fPhhi = TMath::RadToDeg() * TMath::ATan2(nhi[1], nhi[0]);

   fShapeName->SetText(fShape->GetName());
   fEPhi1->SetNumber(fPmini);
   fEPhi2->SetNumber(fPmaxi);
   fSPhi->SetPosition(fPmini,fPmaxi);
   fERmin->SetNumber(fRmini);
   fERmax->SetNumber(fRmaxi);
   fEDz->SetNumber(fDzi);
   fEThlo->SetNumber(fThlo);
   fEPhlo->SetNumber(fPhlo);
   fEThhi->SetNumber(fThhi);
   fEPhhi->SetNumber(fPhhi);
   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);

   if (fInit) ConnectSignals2Slots();
   SetActive();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for phi1.

void TGeoCtubEditor::DoThlo()
{
   Double_t thlo = fEThlo->GetNumber();
   if (thlo <= 90.) {thlo = 91.; fEThlo->SetNumber(thlo);}
   if (thlo > 180.) {thlo = 180.; fEThlo->SetNumber(thlo);}
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for phi1.

void TGeoCtubEditor::DoPhlo()
{
   Double_t phlo = fEPhlo->GetNumber();
   if (phlo >= 360.) {
      phlo = 0.;
      fEPhlo->SetNumber(phlo);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for phi1.

void TGeoCtubEditor::DoThhi()
{
   Double_t thhi = fEThhi->GetNumber();
   if (thhi >= 90.) {thhi = 89.; fEThhi->SetNumber(thhi);}
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for phi1.

void TGeoCtubEditor::DoPhhi()
{
   Double_t phhi = fEPhhi->GetNumber();
   if (phhi >= 360.) {
      phhi = 0.;
      fEPhhi->SetNumber(phhi);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for applying modifications.

void TGeoCtubEditor::DoApply()
{
   fApply->SetEnabled(kFALSE);
   const char *name = fShapeName->GetText();
   if (strcmp(name,fShape->GetName())) fShape->SetName(name);
   Double_t rmin = fERmin->GetNumber();
   Double_t rmax = fERmax->GetNumber();
   if (rmin<0 || rmax<rmin) return;
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
   Double_t thlo = TMath::DegToRad()*fEThlo->GetNumber();
   Double_t phlo = TMath::DegToRad()*fEPhlo->GetNumber();
   Double_t thhi = TMath::DegToRad()*fEThhi->GetNumber();
   Double_t phhi = TMath::DegToRad()*fEPhhi->GetNumber();
   Double_t lx = TMath::Sin(thlo)*TMath::Cos(phlo);
   Double_t ly = TMath::Sin(thlo)*TMath::Sin(phlo);
   Double_t lz = TMath::Cos(thlo);
   Double_t tx = TMath::Sin(thhi)*TMath::Cos(phhi);
   Double_t ty = TMath::Sin(thhi)*TMath::Sin(phhi);
   Double_t tz = TMath::Cos(thhi);
   ((TGeoCtub*)fShape)->SetCtubDimensions(rmin, rmax, dz, phi1, phi2, lx,ly,lz,tx,ty,tz);
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

void TGeoCtubEditor::DoUndo()
{
   fERmin->SetNumber(fRmini);
   fERmax->SetNumber(fRmaxi);
   fEDz->SetNumber(fDzi);
   fEPhi1->SetNumber(fPmini);
   fEPhi2->SetNumber(fPmaxi);
   fSPhi->SetPosition(fPmini,fPmaxi);
   fEThlo->SetNumber(fThlo);
   fEPhlo->SetNumber(fPhlo);
   fEThhi->SetNumber(fThhi);
   fEPhhi->SetNumber(fPhhi);

   DoApply();
   fUndo->SetEnabled(kFALSE);
   fApply->SetEnabled(kFALSE);
}
