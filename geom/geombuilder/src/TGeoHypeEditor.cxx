// @(#):$Id$
// Author: M.Gheata

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoHypeEditor
\ingroup Geometry_builder

Editor for a TGeoHype.

\image html geom_hype_pic.png

\image html geom_hype_ed.png

*/

#include "TGeoHypeEditor.h"
#include "TGeoTabManager.h"
#include "TGeoHype.h"
#include "TGeoManager.h"
#include "TVirtualGeoPainter.h"
#include "TVirtualPad.h"
#include "TView.h"
#include "TMath.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"

ClassImp(TGeoHypeEditor);

enum ETGeoHypeWid {
   kHYPE_NAME, kHYPE_RIN, kHYPE_ROUT,  kHYPE_DZ, kHYPE_STIN,
   kHYPE_STOUT, kHYPE_APPLY, kHYPE_UNDO
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor for Hype editor

TGeoHypeEditor::TGeoHypeEditor(const TGWindow *p, Int_t width,
                                   Int_t height, UInt_t options, Pixel_t back)
   : TGeoGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fShape   = 0;
   fRini = fRouti = fStIni = fStOuti = 0.0;
   fNamei = "";
   fIsModified = kFALSE;
   fIsShapeEditable = kTRUE;

   // TextEntry for shape name
   MakeTitle("Name");
   fShapeName = new TGTextEntry(this, new TGTextBuffer(50), kHYPE_NAME);
   fShapeName->Resize(135, fShapeName->GetDefaultHeight());
   fShapeName->SetToolTipText("Enter the hyperboloid name");
   fShapeName->Associate(this);
   AddFrame(fShapeName, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));

   TGTextEntry *nef;
   MakeTitle("Dimensions");
   // Number entry for Rin
   TGCompositeFrame *f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Rin"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fERin = new TGNumberEntry(f1, 0., 5, kHYPE_RIN);
   fERin->SetNumAttr(TGNumberFormat::kNEAPositive);
   fERin->Resize(100, fERin->GetDefaultHeight());
   nef = (TGTextEntry*)fERin->GetNumberEntry();
   nef->SetToolTipText("Enter the  inner radius ");
   fERin->Associate(this);
   f1->AddFrame(fERin, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // Number entry for Rout
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Rout"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fERout = new TGNumberEntry(f1, 0., 5, kHYPE_ROUT);
   fERout->SetNumAttr(TGNumberFormat::kNEAPositive);
   fERout->Resize(100, fERout->GetDefaultHeight());
   nef = (TGTextEntry*)fERout->GetNumberEntry();
   nef->SetToolTipText("Enter the outer radius");
   fERout->Associate(this);
   f1->AddFrame(fERout, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // Number entry for Dz
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "Dz"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEDz = new TGNumberEntry(f1, 0., 5, kHYPE_DZ);
   fEDz->SetNumAttr(TGNumberFormat::kNEAPositive);
   fEDz->Resize(100, fEDz->GetDefaultHeight());
   nef = (TGTextEntry*)fEDz->GetNumberEntry();
   nef->SetToolTipText("Enter the half-length in Dz");
   fEDz->Associate(this);
   f1->AddFrame(fEDz, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // Number entry for StIn.
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "StIn"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEStIn = new TGNumberEntry(f1, 0., 5, kHYPE_STIN);
   fEStIn->Resize(100, fEStIn->GetDefaultHeight());
   nef = (TGTextEntry*)fEStIn->GetNumberEntry();
   nef->SetToolTipText("Enter the stereo angle for inner surface");
   fEStIn->Associate(this);
   f1->AddFrame(fEStIn, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));

   // Number entry for StOut.
   f1 = new TGCompositeFrame(this, 155, 10, kHorizontalFrame | kFixedWidth);
   f1->AddFrame(new TGLabel(f1, "StOut"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEStOut = new TGNumberEntry(f1, 0., 5, kHYPE_STOUT);
   fEStOut->SetNumAttr(TGNumberFormat::kNEAPositive);
   fEStOut->Resize(100, fEStOut->GetDefaultHeight());
   nef = (TGTextEntry*)fEStOut->GetNumberEntry();
   nef->SetToolTipText("Enter the stereo angle for outer surface");
   fEStOut->Associate(this);
   f1->AddFrame(fEStOut, new TGLayoutHints(kLHintsRight, 2, 2, 4, 4));
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

TGeoHypeEditor::~TGeoHypeEditor()
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

void TGeoHypeEditor::ConnectSignals2Slots()
{
   fApply->Connect("Clicked()", "TGeoHypeEditor", this, "DoApply()");
   fUndo->Connect("Clicked()", "TGeoHypeEditor", this, "DoUndo()");
   fShapeName->Connect("TextChanged(const char *)", "TGeoHypeEditor", this, "DoModified()");
   fERin->Connect("ValueSet(Long_t)", "TGeoHypeEditor", this, "DoRin()");
   fERout->Connect("ValueSet(Long_t)", "TGeoHypeEditor", this, "DoRout()");
   fEDz->Connect("ValueSet(Long_t)", "TGeoHypeEditor", this, "DoDz()");
   fEStIn->Connect("ValueSet(Long_t)", "TGeoHypeEditor", this, "DoStIn()");
   fEStOut->Connect("ValueSet(Long_t)", "TGeoHypeEditor", this, "DoStOut()");
   fERin->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoHypeEditor", this, "DoModified()");
   fERout->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoHypeEditor", this, "DoModified()");
   fEDz->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoHypeEditor", this, "DoModified()");
   fEStIn->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoHypeEditor", this, "DoModified()");
   fEStOut->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoHypeEditor", this, "DoModified()");
   fInit = kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Connect to the selected object.

void TGeoHypeEditor::SetModel(TObject* obj)
{
   if (obj == 0 || (obj->IsA()!=TGeoHype::Class())) {
      SetActive(kFALSE);
      return;
   }
   fShape = (TGeoHype*)obj;
   fRini = fShape->GetRmin();
   fRouti = fShape->GetRmax();
   fDzi = fShape->GetDz();
   fStIni = fShape->GetStIn();
   fStOuti = fShape->GetStOut();
   const char *sname = fShape->GetName();
   if (!strcmp(sname, fShape->ClassName())) fShapeName->SetText("-no_name");
   else {
      fShapeName->SetText(sname);
      fNamei = sname;
   }
   fERin->SetNumber(fRini);
   fERout->SetNumber(fRouti);
   fEDz->SetNumber(fDzi);
   fEStIn->SetNumber(fStIni);
   fEStOut->SetNumber(fStOuti);
   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);

   if (fInit) ConnectSignals2Slots();
   SetActive();
}

////////////////////////////////////////////////////////////////////////////////
/// Check if shape drawing is delayed.

Bool_t TGeoHypeEditor::IsDelayed() const
{
   return (fDelayed->GetState() == kButtonDown);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for name.

void TGeoHypeEditor::DoName()
{
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for applying current settings.

void TGeoHypeEditor::DoApply()
{
   const char *name = fShapeName->GetText();
   if (strcmp(name,fShape->GetName())) fShape->SetName(name);
   Double_t rin = fERin->GetNumber();
   Double_t rout = fERout->GetNumber();
   Double_t dz = fEDz->GetNumber();
   Double_t stin = fEStIn->GetNumber();
   Double_t stout = fEStOut->GetNumber();
   Double_t tin = TMath::Tan(stin*TMath::DegToRad());
   Double_t tout = TMath::Tan(stout*TMath::DegToRad());
   if ((dz<=0) || (rin<0) || (rin>rout) ||
       (rin*rin+tin*tin*dz*dz > rout*rout+tout*tout*dz*dz)) {
      fUndo->SetEnabled();
      fApply->SetEnabled(kFALSE);
      return;
   }
   Double_t param[5];
   param[0] = dz;
   param[1] = rin;
   param[2] = stin;
   param[3] = rout;
   param[4] = stout;
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

void TGeoHypeEditor::DoModified()
{
   fApply->SetEnabled();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for undoing last operation.

void TGeoHypeEditor::DoUndo()
{
   fERin->SetNumber(fRini);
   fERout->SetNumber(fRouti);
   fEDz->SetNumber(fDzi);
   fEStIn->SetNumber(fStIni);
   fEStOut->SetNumber(fStOuti);
   DoApply();
   fUndo->SetEnabled(kFALSE);
   fApply->SetEnabled(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for Rin.

void TGeoHypeEditor::DoRin()
{
   Double_t rin = fERin->GetNumber();
   Double_t rout = fERout->GetNumber();
   Double_t dz = fEDz->GetNumber();
   Double_t stin = fEStIn->GetNumber();
   Double_t stout = fEStOut->GetNumber();
   Double_t tin = TMath::Tan(stin*TMath::DegToRad());
   Double_t tout = TMath::Tan(stout*TMath::DegToRad());
   if (rin<0) {
      rin = 0;
      fERin->SetNumber(rin);
   }
   Double_t rinmax = TMath::Sqrt((rout*rout+tout*tout*dz*dz)/(tin*tin*dz*dz));
   rinmax = TMath::Min(rinmax, rout);
   if (rin > rinmax) {
      rin = rinmax-1.e-6;
      fERin->SetNumber(rin);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for Rout.

void TGeoHypeEditor::DoRout()
{
   Double_t rin = fERin->GetNumber();
   Double_t rout = fERout->GetNumber();
   Double_t dz = fEDz->GetNumber();
   Double_t stin = fEStIn->GetNumber();
   Double_t stout = fEStOut->GetNumber();
   Double_t tin = TMath::Tan(stin*TMath::DegToRad());
   Double_t tout = TMath::Tan(stout*TMath::DegToRad());
   Double_t routmin = TMath::Sqrt((rin*rin+tin*tin*dz*dz)/(tout*tout*dz*dz));
   routmin = TMath::Max(routmin,rin);
   if (rout < routmin) {
      rout = routmin+1.e-6;
      fERout->SetNumber(rout);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for Z.

void TGeoHypeEditor::DoDz()
{
   Double_t rin = fERin->GetNumber();
   Double_t rout = fERout->GetNumber();
   Double_t dz = fEDz->GetNumber();
   Double_t stin = fEStIn->GetNumber();
   Double_t stout = fEStOut->GetNumber();
   if (TMath::Abs(stin-stout)<1.e-6) {
      stin = stout+1.;
      fEStIn->SetNumber(stin);
   }
   Double_t tin = TMath::Tan(stin*TMath::DegToRad());
   Double_t tout = TMath::Tan(stout*TMath::DegToRad());
   if (dz<=0) {
      dz = 0.1;
      fEDz->SetNumber(dz);
   }
   Double_t dzmax = TMath::Sqrt((rout*rout-rin*rin)/(tin*tin-tout*tout));
   if (dz>dzmax) {
      dz = dzmax;
      fEDz->SetNumber(dz);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for StIn.

void TGeoHypeEditor::DoStIn()
{
   Double_t rin = fERin->GetNumber();
   Double_t rout = fERout->GetNumber();
   Double_t dz = fEDz->GetNumber();
   Double_t stin = fEStIn->GetNumber();
   Double_t stout = fEStOut->GetNumber();
   if (stin >= 90) {
      stin = 89.;
      fEStIn->SetNumber(stin);
   }
   Double_t tin = TMath::Tan(stin*TMath::DegToRad());
   Double_t tout = TMath::Tan(stout*TMath::DegToRad());
   Double_t tinmax = TMath::Sqrt(tout*tout+(rout*rout-rin*rin)/(dz*dz));
   if (tin>tinmax) {
      tin = tinmax-1.e-6;
      stin = TMath::RadToDeg()*TMath::ATan(tin);
      fEStIn->SetNumber(stin);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for StOut.

void TGeoHypeEditor::DoStOut()
{
   Double_t rin = fERin->GetNumber();
   Double_t rout = fERout->GetNumber();
   Double_t dz = fEDz->GetNumber();
   Double_t stin = fEStIn->GetNumber();
   Double_t stout = fEStOut->GetNumber();
   if (stout > 90) {
      stout = 89;
      fEStOut->SetNumber(stout);
   }
   Double_t tin = TMath::Tan(stin*TMath::DegToRad());
   Double_t tout = TMath::Tan(stout*TMath::DegToRad());
   Double_t tinmin = TMath::Sqrt((rout*rout-rin*rin)/(dz*dz));
   if (tin < tinmin) {
      tin = tinmin;
      stin = TMath::RadToDeg()*TMath::ATan(tin);
      fEStIn->SetNumber(stin);
   }
   Double_t toutmin = TMath::Sqrt(tin*tin -tinmin*tinmin);
   if (tout < toutmin) {
      tout = toutmin+1.e-6;
      stout = TMath::RadToDeg()*TMath::ATan(tout);
      fEStOut->SetNumber(stout);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}
