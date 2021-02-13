// @(#):$Id$
// Author: M.Gheata

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoTranslationEditor
\ingroup Geometry_builder

Editor for a TGeoTranslation.

*/

#include "TGeoMatrixEditor.h"
#include "TGeoTabManager.h"
#include "TGeoMatrix.h"
#include "TVirtualPad.h"
#include "TGButton.h"
#include "TGButtonGroup.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"

ClassImp(TGeoTranslationEditor);

enum ETGeoMatrixWid {
   kMATRIX_NAME, kMATRIX_DX, kMATRIX_DY, kMATRIX_DZ,
   kMATRIX_PHI, kMATRIX_THETA, kMATRIX_PSI,
   kMATRIX_APPLY, kMATRIX_CANCEL, kMATRIX_UNDO
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor for translation editor

TGeoTranslationEditor::TGeoTranslationEditor(const TGWindow *p, Int_t width,
                                   Int_t height, UInt_t options, Pixel_t back)
   : TGeoGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fTranslation   = 0;
   fDxi = fDyi = fDzi = 0.0;
   fNamei = "";
   fIsModified = kFALSE;
   fIsEditable = kFALSE;

   // TextEntry for name
   MakeTitle("Name");
   fTransName = new TGTextEntry(this, new TGTextBuffer(50), kMATRIX_NAME);
   fTransName->Resize(135, fTransName->GetDefaultHeight());
   fTransName->SetToolTipText("Enter the translation name");
   fTransName->Associate(this);
   AddFrame(fTransName, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));

   TGTextEntry *nef;
   MakeTitle("Translations on axes");
   TGCompositeFrame *compxyz = new TGCompositeFrame(this, 118, 30, kVerticalFrame | kRaisedFrame | kDoubleBorder);
   // Number entry for dx
   TGCompositeFrame *f1 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f1->AddFrame(new TGLabel(f1, "DX"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fTransDx = new TGNumberEntry(f1, 0., 5, kMATRIX_DX);
   nef = (TGTextEntry*)fTransDx->GetNumberEntry();
   nef->SetToolTipText("Enter the translation on X");
   fTransDx->Associate(this);
   f1->AddFrame(fTransDx, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));

   // Number entry for dy
   TGCompositeFrame *f2 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f2->AddFrame(new TGLabel(f2, "DY"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fTransDy = new TGNumberEntry(f2, 0., 5, kMATRIX_DY);
   nef = (TGTextEntry*)fTransDy->GetNumberEntry();
   nef->SetToolTipText("Enter the translation on Y");
   fTransDy->Associate(this);
   f2->AddFrame(fTransDy, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));
   compxyz->AddFrame(f2, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));

   // Number entry for dx
   TGCompositeFrame *f3 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f3->AddFrame(new TGLabel(f3, "DZ"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fTransDz = new TGNumberEntry(f3, 0., 5, kMATRIX_DZ);
   nef = (TGTextEntry*)fTransDz->GetNumberEntry();
   nef->SetToolTipText("Enter the translation on Z");
   fTransDz->Associate(this);
   f3->AddFrame(fTransDz, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));
   compxyz->AddFrame(f3, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));

   compxyz->Resize(150,30);
   AddFrame(compxyz, new TGLayoutHints(kLHintsLeft, 6, 6, 2, 2));

   // Buttons
   TGCompositeFrame *f23 = new TGCompositeFrame(this, 118, 20, kHorizontalFrame | kSunkenFrame | kDoubleBorder);
   fApply = new TGTextButton(f23, "&Apply");
   f23->AddFrame(fApply, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));
   fApply->Associate(this);
   fCancel = new TGTextButton(f23, "&Cancel");
   f23->AddFrame(fCancel, new TGLayoutHints(kLHintsCenterX, 2, 2, 4, 4));
   fCancel->Associate(this);
   fUndo = new TGTextButton(f23, " &Undo ");
   f23->AddFrame(fUndo, new TGLayoutHints(kLHintsRight , 2, 2, 4, 4));
   fUndo->Associate(this);
   AddFrame(f23,  new TGLayoutHints(kLHintsLeft, 6, 6, 2, 2));
   fUndo->SetSize(fCancel->GetSize());
   fApply->SetSize(fCancel->GetSize());
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGeoTranslationEditor::~TGeoTranslationEditor()
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

void TGeoTranslationEditor::ConnectSignals2Slots()
{
   fApply->Connect("Clicked()", "TGeoTranslationEditor", this, "DoApply()");
   fCancel->Connect("Clicked()", "TGeoTranslationEditor", this, "DoCancel()");
   fUndo->Connect("Clicked()", "TGeoTranslationEditor", this, "DoUndo()");
   fTransName->Connect("TextChanged(const char *)", "TGeoTranslationEditor", this, "DoModified()");
   fTransDx->Connect("ValueSet(Long_t)", "TGeoTranslationEditor", this, "DoDx()");
   fTransDy->Connect("ValueSet(Long_t)", "TGeoTranslationEditor", this, "DoDy()");
   fTransDz->Connect("ValueSet(Long_t)", "TGeoTranslationEditor", this, "DoDz()");
   fTransDx->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTranslationEditor", this, "DoDx()");
   fTransDy->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTranslationEditor", this, "DoDy()");
   fTransDz->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTranslationEditor", this, "DoDz()");
   fInit = kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Connect to the new matrix.

void TGeoTranslationEditor::SetModel(TObject* obj)
{
   if (obj == 0 || (obj->IsA()!=TGeoTranslation::Class())) {
      SetActive(kFALSE);
      return;
   }
   fTranslation = (TGeoTranslation*)obj;
   fDxi = fTranslation->GetTranslation()[0];
   fDyi = fTranslation->GetTranslation()[1];
   fDzi = fTranslation->GetTranslation()[2];
   const char *sname = fTranslation->GetName();
   if (!strcmp(sname, fTranslation->ClassName())) fTransName->SetText("no_name");
   else {
      fTransName->SetText(sname);
      fNamei = sname;
   }
   fTransDx->SetNumber(fDxi);
   fTransDy->SetNumber(fDyi);
   fTransDz->SetNumber(fDzi);
   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);
   fCancel->SetEnabled(kFALSE);

   if (fInit) ConnectSignals2Slots();
   SetActive();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for name.

void TGeoTranslationEditor::DoName()
{
   const char *name = fTransName->GetText();
   if (!strcmp(name, "no_name") || !strcmp(name, fTranslation->GetName())) return;
   fTranslation->SetName(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for checking parameters.

Bool_t TGeoTranslationEditor::DoParameters()
{
   Double_t dx = fTransDx->GetNumber();
   Double_t dy = fTransDy->GetNumber();
   Double_t dz = fTransDz->GetNumber();
   Bool_t changed = kFALSE;
   if (dx != fTranslation->GetTranslation()[0] ||
       dy != fTranslation->GetTranslation()[1] ||
       dz != fTranslation->GetTranslation()[2]) changed = kTRUE;
   if (!changed) return kFALSE;
   fUndo->SetEnabled();
   fTranslation->SetTranslation(dx, dy, dz);
   if (fPad) {
      fPad->Modified();
      fPad->Update();
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for applying changes.

void TGeoTranslationEditor::DoApply()
{
   DoName();
   if (DoParameters()) {
      fUndo->SetEnabled();
      fCancel->SetEnabled(kFALSE);
      fApply->SetEnabled(kFALSE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for cancelling last modifications non-applied.

void TGeoTranslationEditor::DoCancel()
{
   if (!fNamei.Length()) fTransName->SetText("no_name");
   else fTransName->SetText(fNamei.Data());
   fTransDx->SetNumber(fDxi);
   fTransDy->SetNumber(fDyi);
   fTransDz->SetNumber(fDzi);
   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);
   fCancel->SetEnabled(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for notifying changes.

void TGeoTranslationEditor::DoModified()
{
   fApply->SetEnabled();
   if (fUndo->GetState()==kButtonDisabled) fCancel->SetEnabled();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for undoing last operation.

void TGeoTranslationEditor::DoUndo()
{
   DoCancel();
   DoParameters();
   fCancel->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);
   fApply->SetEnabled(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for dx.

void TGeoTranslationEditor::DoDx()
{
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for dx.

void TGeoTranslationEditor::DoDy()
{
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for dx.

void TGeoTranslationEditor::DoDz()
{
   DoModified();
}

/** \class TGeoRotationEditor
\ingroup Geometry_builder

Editor for a TGeoRotation.

*/

ClassImp(TGeoRotationEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor for rotation editor

TGeoRotationEditor::TGeoRotationEditor(const TGWindow *p, Int_t width,
                                   Int_t height, UInt_t options, Pixel_t back)
   : TGeoGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fRotation   = 0;
   fPhii = fThetai = fPsii = 0.0;
   fAngleX = fAngleY = fAngleZ = 0.0;
   fNamei = "";
   fIsModified = kFALSE;
   fIsEditable = kFALSE;

   // TextEntry for name
   MakeTitle("Name");
   fRotName = new TGTextEntry(this, new TGTextBuffer(50), kMATRIX_NAME);
   fRotName->Resize(135, fRotName->GetDefaultHeight());
   fRotName->SetToolTipText("Enter the rotation name");
   fRotName->Associate(this);
   AddFrame(fRotName, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));

   TGTextEntry *nef;
   MakeTitle("Euler angles");
   TGCompositeFrame *compxyz = new TGCompositeFrame(this, 140, 30, kVerticalFrame | kRaisedFrame | kDoubleBorder);
   // Number entry for phi angle
   TGCompositeFrame *f1 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f1->AddFrame(new TGLabel(f1, " PHI "), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fRotPhi = new TGNumberEntry(f1, 0., 5, kMATRIX_PHI);
   nef = (TGTextEntry*)fRotPhi->GetNumberEntry();
   nef->SetToolTipText("Modify the first rotation angle about Z");
   fRotPhi->Associate(this);
   fRotPhi->Resize(90, fRotPhi->GetDefaultHeight());
   f1->AddFrame(fRotPhi, new TGLayoutHints(kLHintsRight , 2, 2, 2, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));

   // Number entry for theta angle
   TGCompositeFrame *f2 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f2->AddFrame(new TGLabel(f2, "THETA"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fRotTheta = new TGNumberEntry(f2, 0., 5, kMATRIX_THETA);
   nef = (TGTextEntry*)fRotTheta->GetNumberEntry();
   nef->SetToolTipText("Modify the second rotation angle about the new X");
   fRotTheta->Associate(this);
   fRotTheta->Resize(90, fRotTheta->GetDefaultHeight());
   f2->AddFrame(fRotTheta, new TGLayoutHints(kLHintsRight , 2, 2, 2, 2));
   compxyz->AddFrame(f2, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));

   // Number entry for psi angle
   TGCompositeFrame *f3 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f3->AddFrame(new TGLabel(f3, " PSI "), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fRotPsi = new TGNumberEntry(f3, 0., 5, kMATRIX_PSI);
   nef = (TGTextEntry*)fRotPsi->GetNumberEntry();
   nef->SetToolTipText("Modify the third rotation angle about Z");
   fRotPsi->Associate(this);
   fRotPsi->Resize(90, fRotPsi->GetDefaultHeight());
   f3->AddFrame(fRotPsi, new TGLayoutHints(kLHintsRight , 2, 2, 2, 2));
   compxyz->AddFrame(f3, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));

   compxyz->Resize(150,compxyz->GetDefaultHeight());
   AddFrame(compxyz, new TGLayoutHints(kLHintsLeft, 6, 6, 2, 2));

   MakeTitle("Rotate about axis");
   compxyz = new TGCompositeFrame(this, 140, 30, kVerticalFrame | kRaisedFrame | kDoubleBorder);
   // Number entry for rotation angle about one axis
   f1 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f1->AddFrame(new TGLabel(f1, "ANGLE"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fRotAxis = new TGNumberEntry(f1, 0., 5, kMATRIX_DX);
   nef = (TGTextEntry*)fRotAxis->GetNumberEntry();
   nef->SetToolTipText("Enter the new rotation angle about the selected axis");
   fRotAxis->Associate(this);
   fRotAxis->Resize(90, fRotAxis->GetDefaultHeight());
   f1->AddFrame(fRotAxis, new TGLayoutHints(kLHintsRight , 2, 2, 2, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));

   // Radio buttons group for axis selection
   TGHButtonGroup *bg1 = new TGHButtonGroup(compxyz, " Axis ");
   fRotX = new TGRadioButton(bg1, " &X ", kMATRIX_DX);
   fRotY = new TGRadioButton(bg1, " &Y ", kMATRIX_DY);
   fRotZ = new TGRadioButton(bg1, " &Z ", kMATRIX_DZ);
   bg1->SetRadioButtonExclusive();
   bg1->Show();
   compxyz->AddFrame(bg1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));

   compxyz->Resize(150,compxyz->GetDefaultHeight());
   AddFrame(compxyz, new TGLayoutHints(kLHintsLeft, 6, 6, 2, 2));


   // Buttons
   TGCompositeFrame *f23 = new TGCompositeFrame(this, 118, 20, kHorizontalFrame | kSunkenFrame | kDoubleBorder);
   fApply = new TGTextButton(f23, "Apply");
   f23->AddFrame(fApply, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));
   fApply->Associate(this);
   fCancel = new TGTextButton(f23, "Cancel");
   f23->AddFrame(fCancel, new TGLayoutHints(kLHintsCenterX, 2, 2, 4, 4));
   fCancel->Associate(this);
   fUndo = new TGTextButton(f23, " Undo ");
   f23->AddFrame(fUndo, new TGLayoutHints(kLHintsRight , 2, 2, 4, 4));
   fUndo->Associate(this);
   AddFrame(f23,  new TGLayoutHints(kLHintsLeft, 6, 6, 2, 2));
   fUndo->SetSize(fCancel->GetSize());
   fApply->SetSize(fCancel->GetSize());
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoRotationEditor::~TGeoRotationEditor()
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

void TGeoRotationEditor::ConnectSignals2Slots()
{
   fApply->Connect("Clicked()", "TGeoRotationEditor", this, "DoApply()");
   fCancel->Connect("Clicked()", "TGeoRotationEditor", this, "DoCancel()");
   fUndo->Connect("Clicked()", "TGeoRotationEditor", this, "DoUndo()");
   fRotName->Connect("TextChanged(const char *)", "TGeoRotationEditor", this, "DoModified()");
   fRotPhi->Connect("ValueSet(Long_t)", "TGeoRotationEditor", this, "DoRotPhi()");
   fRotTheta->Connect("ValueSet(Long_t)", "TGeoRotationEditor", this, "DoRotTheta()");
   fRotPsi->Connect("ValueSet(Long_t)", "TGeoRotationEditor", this, "DoRotPsi()");
   fRotAxis->Connect("ValueSet(Long_t)", "TGeoRotationEditor", this, "DoRotAngle()");
   fInit = kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Connect to the selected rotation.

void TGeoRotationEditor::SetModel(TObject* obj)
{
   if (obj == 0 || (obj->IsA()!=TGeoRotation::Class())) {
      SetActive(kFALSE);
      return;
   }
   fRotation = (TGeoRotation*)obj;
   fRotation->GetAngles(fPhii, fThetai, fPsii);
   const char *sname = fRotation->GetName();
   if (!strcmp(sname, fRotation->ClassName())) fRotName->SetText("no_name");
   else {
      fRotName->SetText(sname);
      fNamei = sname;
   }
   fRotPhi->SetNumber(fPhii);
   fRotTheta->SetNumber(fThetai);
   fRotPsi->SetNumber(fPsii);
   fRotAxis->SetNumber(0.0);

   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);
   fCancel->SetEnabled(kFALSE);

   if (fInit) ConnectSignals2Slots();
   SetActive();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for name.

void TGeoRotationEditor::DoName()
{
   const char *name = fRotName->GetText();
   if (!strcmp(name, "no_name") || !strcmp(name, fRotation->GetName())) return;
   fRotation->SetName(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for phi (Euler X convention)

void TGeoRotationEditor::DoRotPhi()
{
   if (fRotPhi->GetNumber() < 0.) fRotPhi->SetNumber(fRotPhi->GetNumber()+360.);
   if (fRotPhi->GetNumber() >= 360.) fRotPhi->SetNumber(fRotPhi->GetNumber()-360.);
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for theta (Euler X convention)

void TGeoRotationEditor::DoRotTheta()
{
   if (fRotTheta->GetNumber() < 0.) fRotTheta->SetNumber(fRotTheta->GetNumber()+360.);
   if (fRotTheta->GetNumber() >= 360.) fRotTheta->SetNumber(fRotTheta->GetNumber()-360.);
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for psi (Euler X convention)

void TGeoRotationEditor::DoRotPsi()
{
   if (fRotPsi->GetNumber() < 0.) fRotPsi->SetNumber(fRotPsi->GetNumber()+360.);
   if (fRotPsi->GetNumber() >= 360.) fRotPsi->SetNumber(fRotPsi->GetNumber()-360.);
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for additional rotation about one axis.

void TGeoRotationEditor::DoRotAngle()
{
   if (fRotAxis->GetNumber() < 0.) fRotAxis->SetNumber(fRotAxis->GetNumber()+360.);
   if (fRotAxis->GetNumber() >= 360.) fRotAxis->SetNumber(fRotAxis->GetNumber()-360.);
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for checking parameters.

Bool_t TGeoRotationEditor::DoParameters()
{
   Double_t phi = fRotPhi->GetNumber();
   Double_t theta = fRotTheta->GetNumber();
   Double_t psi = fRotPsi->GetNumber();
   Double_t angle = fRotAxis->GetNumber();
   Double_t phi0 = 0., theta0 = 0., psi0 = 0.;
   fRotation->GetAngles(phi0,theta0,psi0);
   Bool_t changed = kFALSE;
   if (phi != psi0 || theta != theta0 || psi != psi0) changed = kTRUE;
   if (changed) fRotation->SetAngles(phi, theta, psi);
   // Check if we have to rotate about one axis
   if (angle != 0.0) {
      if (fRotX->IsOn()) {fRotation->RotateX(angle); changed = kTRUE;}
      if (fRotY->IsOn()) {fRotation->RotateY(angle); changed = kTRUE;}
      if (fRotZ->IsOn()) {fRotation->RotateZ(angle); changed = kTRUE;}
   }
   if (!changed) return kFALSE;
   fRotAxis->SetNumber(0.0);
   fUndo->SetEnabled();
   if (fPad) {
      fPad->Modified();
      fPad->Update();
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for applying modifications.

void TGeoRotationEditor::DoApply()
{
   DoName();
   if (DoParameters()) {
      fUndo->SetEnabled();
      fCancel->SetEnabled(kFALSE);
      fApply->SetEnabled(kFALSE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for cancelling last un-applied operations.

void TGeoRotationEditor::DoCancel()
{
   if (!fNamei.Length()) fRotName->SetText("no_name");
   else fRotName->SetText(fNamei.Data());
   fRotPhi->SetNumber(fPhii);
   fRotTheta->SetNumber(fThetai);
   fRotPsi->SetNumber(fPsii);
   fRotAxis->SetNumber(0.0);
   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);
   fCancel->SetEnabled(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for notifying changes.

void TGeoRotationEditor::DoModified()
{
   fApply->SetEnabled();
   if (fUndo->GetState()==kButtonDisabled) fCancel->SetEnabled();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for undoing last changes.

void TGeoRotationEditor::DoUndo()
{
   DoCancel();
   DoParameters();
   fCancel->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);
   fApply->SetEnabled(kFALSE);
}

/** \class TGeoCombiTransEditor
\ingroup Geometry_builder

Editor for a TGeoCombiTrans.

*/

ClassImp(TGeoCombiTransEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor for combi matrix editor

TGeoCombiTransEditor::TGeoCombiTransEditor(const TGWindow *p, Int_t width,
                                   Int_t height, UInt_t options, Pixel_t back)
   : TGeoGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fCombi   = 0;
   fPhii = fThetai = fPsii = 0.0;
   fDxi = fDyi = fDzi = 0.0;
   fAngleX = fAngleY = fAngleZ = 0.0;
   fNamei = "";
   fIsModified = kFALSE;
   fIsEditable = kFALSE;

   // TextEntry for name
   MakeTitle("Name");
   fRotName = new TGTextEntry(this, new TGTextBuffer(50), kMATRIX_NAME);
   fRotName->Resize(135, fRotName->GetDefaultHeight());
   fRotName->SetToolTipText("Enter the rotation name");
   fRotName->Associate(this);
   AddFrame(fRotName, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));

   TGTextEntry *nef;
   MakeTitle("Translations on axes");
   TGCompositeFrame *compxyz = new TGCompositeFrame(this, 118, 30, kVerticalFrame | kRaisedFrame | kDoubleBorder);
   // Number entry for dx
   TGCompositeFrame *f1 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f1->AddFrame(new TGLabel(f1, "DX"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fTransDx = new TGNumberEntry(f1, 0., 5, kMATRIX_DX);
   nef = (TGTextEntry*)fTransDx->GetNumberEntry();
   nef->SetToolTipText("Enter the translation on X");
   fTransDx->Associate(this);
   f1->AddFrame(fTransDx, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));

   // Number entry for dy
   TGCompositeFrame *f2 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f2->AddFrame(new TGLabel(f2, "DY"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fTransDy = new TGNumberEntry(f2, 0., 5, kMATRIX_DY);
   nef = (TGTextEntry*)fTransDy->GetNumberEntry();
   nef->SetToolTipText("Enter the translation on Y");
   fTransDy->Associate(this);
   f2->AddFrame(fTransDy, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));
   compxyz->AddFrame(f2, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));

   // Number entry for dx
   TGCompositeFrame *f3 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f3->AddFrame(new TGLabel(f3, "DZ"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fTransDz = new TGNumberEntry(f3, 0., 5, kMATRIX_DZ);
   nef = (TGTextEntry*)fTransDz->GetNumberEntry();
   nef->SetToolTipText("Enter the translation on Z");
   fTransDz->Associate(this);
   f3->AddFrame(fTransDz, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));
   compxyz->AddFrame(f3, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));

   compxyz->Resize(150,30);
   AddFrame(compxyz, new TGLayoutHints(kLHintsLeft, 6, 6, 2, 2));


   MakeTitle("Euler angles");
   compxyz = new TGCompositeFrame(this, 140, 30, kVerticalFrame | kRaisedFrame | kDoubleBorder);
   // Number entry for phi angle
   f1 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f1->AddFrame(new TGLabel(f1, " PHI "), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fRotPhi = new TGNumberEntry(f1, 0., 5, kMATRIX_PHI);
   nef = (TGTextEntry*)fRotPhi->GetNumberEntry();
   nef->SetToolTipText("Modify the first rotation angle about Z");
   fRotPhi->Associate(this);
   fRotPhi->Resize(90, fRotPhi->GetDefaultHeight());
   f1->AddFrame(fRotPhi, new TGLayoutHints(kLHintsRight , 2, 2, 2, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));

   // Number entry for theta angle
   f2 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f2->AddFrame(new TGLabel(f2, "THETA"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fRotTheta = new TGNumberEntry(f2, 0., 5, kMATRIX_THETA);
   nef = (TGTextEntry*)fRotTheta->GetNumberEntry();
   nef->SetToolTipText("Modify the second rotation angle about the new X");
   fRotTheta->Associate(this);
   fRotTheta->Resize(90, fRotTheta->GetDefaultHeight());
   f2->AddFrame(fRotTheta, new TGLayoutHints(kLHintsRight , 2, 2, 2, 2));
   compxyz->AddFrame(f2, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));

   // Number entry for psi angle
   f3 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f3->AddFrame(new TGLabel(f3, " PSI "), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fRotPsi = new TGNumberEntry(f3, 0., 5, kMATRIX_PSI);
   nef = (TGTextEntry*)fRotPsi->GetNumberEntry();
   nef->SetToolTipText("Modify the third rotation angle about Z");
   fRotPsi->Associate(this);
   fRotPsi->Resize(90, fRotPsi->GetDefaultHeight());
   f3->AddFrame(fRotPsi, new TGLayoutHints(kLHintsRight , 2, 2, 2, 2));
   compxyz->AddFrame(f3, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));

   compxyz->Resize(150,compxyz->GetDefaultHeight());
   AddFrame(compxyz, new TGLayoutHints(kLHintsLeft, 6, 6, 2, 2));

   MakeTitle("Rotate about axis");
   compxyz = new TGCompositeFrame(this, 140, 30, kVerticalFrame | kRaisedFrame | kDoubleBorder);
   // Number entry for rotation angle about one axis
   f1 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f1->AddFrame(new TGLabel(f1, "ANGLE"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fRotAxis = new TGNumberEntry(f1, 0., 5, kMATRIX_DX);
   nef = (TGTextEntry*)fRotAxis->GetNumberEntry();
   nef->SetToolTipText("Enter the new rotation angle about the selected axis");
   fRotAxis->Associate(this);
   fRotAxis->Resize(90, fRotAxis->GetDefaultHeight());
   f1->AddFrame(fRotAxis, new TGLayoutHints(kLHintsRight , 2, 2, 2, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));

   // Radio buttons group for axis selection
   TGHButtonGroup *bg1 = new TGHButtonGroup(compxyz, " Axis ");
   fRotX = new TGRadioButton(bg1, " &X ", kMATRIX_DX);
   fRotY = new TGRadioButton(bg1, " &Y ", kMATRIX_DY);
   fRotZ = new TGRadioButton(bg1, " &Z ", kMATRIX_DZ);
   bg1->SetRadioButtonExclusive();
   bg1->Show();
   compxyz->AddFrame(bg1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2));

   compxyz->Resize(150,compxyz->GetDefaultHeight());
   AddFrame(compxyz, new TGLayoutHints(kLHintsLeft, 6, 6, 2, 2));


   // Buttons
   TGCompositeFrame *f23 = new TGCompositeFrame(this, 118, 20, kHorizontalFrame | kSunkenFrame | kDoubleBorder);
   fApply = new TGTextButton(f23, "&Apply");
   f23->AddFrame(fApply, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));
   fApply->Associate(this);
   fCancel = new TGTextButton(f23, "&Cancel");
   f23->AddFrame(fCancel, new TGLayoutHints(kLHintsCenterX, 2, 2, 4, 4));
   fCancel->Associate(this);
   fUndo = new TGTextButton(f23, " &Undo ");
   f23->AddFrame(fUndo, new TGLayoutHints(kLHintsRight , 2, 2, 4, 4));
   fUndo->Associate(this);
   AddFrame(f23,  new TGLayoutHints(kLHintsLeft, 6, 6, 2, 2));
   fUndo->SetSize(fCancel->GetSize());
   fApply->SetSize(fCancel->GetSize());
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoCombiTransEditor::~TGeoCombiTransEditor()
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

void TGeoCombiTransEditor::ConnectSignals2Slots()
{
   fApply->Connect("Clicked()", "TGeoCombiTransEditor", this, "DoApply()");
   fCancel->Connect("Clicked()", "TGeoCombiTransEditor", this, "DoCancel()");
   fUndo->Connect("Clicked()", "TGeoCombiTransEditor", this, "DoUndo()");
   fRotName->Connect("TextChanged(const char *)", "TGeoCombiTransEditor", this, "DoModified()");
   fRotPhi->Connect("ValueSet(Long_t)", "TGeoCombiTransEditor", this, "DoRotPhi()");
   fRotTheta->Connect("ValueSet(Long_t)", "TGeoCombiTransEditor", this, "DoRotTheta()");
   fRotPsi->Connect("ValueSet(Long_t)", "TGeoCombiTransEditor", this, "DoRotPsi()");
   fRotAxis->Connect("ValueSet(Long_t)", "TGeoCombiTransEditor", this, "DoRotAngle()");
   fTransDx->Connect("ValueSet(Long_t)", "TGeoCombiTransEditor", this, "DoDx()");
   fTransDy->Connect("ValueSet(Long_t)", "TGeoCombiTransEditor", this, "DoDy()");
   fTransDz->Connect("ValueSet(Long_t)", "TGeoCombiTransEditor", this, "DoDz()");
   fTransDx->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoCombiTransEditor", this, "DoDx()");
   fTransDy->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoCombiTransEditor", this, "DoDy()");
   fTransDz->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoCombiTransEditor", this, "DoDz()");
   fInit = kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Connect to the selected combi matrix.

void TGeoCombiTransEditor::SetModel(TObject* obj)
{
   if (obj == 0 || (obj->IsA()!=TGeoCombiTrans::Class())) {
      SetActive(kFALSE);
      return;
   }
   fCombi = (TGeoCombiTrans*)obj;
   TGeoRotation *rot = fCombi->GetRotation();
   if (rot) rot->GetAngles(fPhii, fThetai, fPsii);
   const char *sname = fCombi->GetName();
   if (!strcmp(sname, fCombi->ClassName())) fRotName->SetText("no_name");
   else {
      fRotName->SetText(sname);
      fNamei = sname;
   }

   fDxi = fCombi->GetTranslation()[0];
   fDyi = fCombi->GetTranslation()[1];
   fDzi = fCombi->GetTranslation()[2];
   fTransDx->SetNumber(fDxi);
   fTransDy->SetNumber(fDyi);
   fTransDz->SetNumber(fDzi);

   fRotPhi->SetNumber(fPhii);
   fRotTheta->SetNumber(fThetai);
   fRotPsi->SetNumber(fPsii);
   fRotAxis->SetNumber(0.0);

   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);
   fCancel->SetEnabled(kFALSE);

   if (fInit) ConnectSignals2Slots();
   SetActive();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for name.

void TGeoCombiTransEditor::DoName()
{
   const char *name = fRotName->GetText();
   if (!strcmp(name, "no_name") || !strcmp(name, fCombi->GetName())) return;
   fCombi->SetName(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for phi (Euler X convention)

void TGeoCombiTransEditor::DoRotPhi()
{
   if (fRotPhi->GetNumber() < 0.) fRotPhi->SetNumber(fRotPhi->GetNumber()+360.);
   if (fRotPhi->GetNumber() >= 360.) fRotPhi->SetNumber(fRotPhi->GetNumber()-360.);
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for theta (Euler X convention)

void TGeoCombiTransEditor::DoRotTheta()
{
   if (fRotTheta->GetNumber() < 0.) fRotTheta->SetNumber(fRotTheta->GetNumber()+360.);
   if (fRotTheta->GetNumber() >= 360.) fRotTheta->SetNumber(fRotTheta->GetNumber()-360.);
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for psi (Euler X convention)

void TGeoCombiTransEditor::DoRotPsi()
{
   if (fRotPsi->GetNumber() < 0.) fRotPsi->SetNumber(fRotPsi->GetNumber()+360.);
   if (fRotPsi->GetNumber() >= 360.) fRotPsi->SetNumber(fRotPsi->GetNumber()-360.);
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for additional rotation about one axis.

void TGeoCombiTransEditor::DoRotAngle()
{
   if (fRotAxis->GetNumber() < 0.) fRotAxis->SetNumber(fRotAxis->GetNumber()+360.);
   if (fRotAxis->GetNumber() >= 360.) fRotAxis->SetNumber(fRotAxis->GetNumber()-360.);
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for checking parameters.

Bool_t TGeoCombiTransEditor::DoParameters()
{
   Double_t dx = fTransDx->GetNumber();
   Double_t dy = fTransDy->GetNumber();
   Double_t dz = fTransDz->GetNumber();
   Bool_t changedtr = kFALSE;
   if (dx != fCombi->GetTranslation()[0] ||
       dy != fCombi->GetTranslation()[1] ||
       dz != fCombi->GetTranslation()[2]) changedtr = kTRUE;
   if (changedtr) fCombi->SetTranslation(dx, dy, dz);
   Double_t phi = fRotPhi->GetNumber();
   Double_t theta = fRotTheta->GetNumber();
   Double_t psi = fRotPsi->GetNumber();
   Double_t angle = fRotAxis->GetNumber();
   Double_t phi0 = 0., theta0 = 0., psi0 = 0.;
   TGeoRotation *rot = fCombi->GetRotation();
   if (rot) rot->GetAngles(phi0,theta0,psi0);
   else {
      if (phi!=fPhii || theta!=fThetai || psi!=fPsii) {
         TGeoRotation r("rot",10.,0.,0.);
         fCombi->SetRotation(r);
         rot = fCombi->GetRotation();
         rot->SetAngles(0.,0.,0.);
      }
   }
   Bool_t changed = kFALSE;
   if (phi != psi0 || theta != theta0 || psi != psi0) changed = kTRUE;
   if (changed && rot) rot->SetAngles(phi, theta, psi);
   // Check if we have to rotate about one axis
   if (angle != 0.0) {
      if (fRotX->IsOn()) {fCombi->RotateX(angle); changed = kTRUE;}
      if (fRotY->IsOn()) {fCombi->RotateY(angle); changed = kTRUE;}
      if (fRotZ->IsOn()) {fCombi->RotateZ(angle); changed = kTRUE;}
   }
   if (changedtr) changed = kTRUE;
   if (!changed) return kFALSE;
   fRotAxis->SetNumber(0.0);
   fUndo->SetEnabled();
   if (fPad) {
      fPad->Modified();
      fPad->Update();
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for applying modifications.

void TGeoCombiTransEditor::DoApply()
{
   DoName();
   if (DoParameters()) {
      fUndo->SetEnabled();
      fCancel->SetEnabled(kFALSE);
      fApply->SetEnabled(kFALSE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for cancelling last un-applied operations.

void TGeoCombiTransEditor::DoCancel()
{
   if (!fNamei.Length()) fRotName->SetText("no_name");
   else fRotName->SetText(fNamei.Data());
   fTransDx->SetNumber(fDxi);
   fTransDy->SetNumber(fDyi);
   fTransDz->SetNumber(fDzi);
   fRotPhi->SetNumber(fPhii);
   fRotTheta->SetNumber(fThetai);
   fRotPsi->SetNumber(fPsii);
   fRotAxis->SetNumber(0.0);
   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);
   fCancel->SetEnabled(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for notifying changes.

void TGeoCombiTransEditor::DoModified()
{
   fApply->SetEnabled();
   if (fUndo->GetState()==kButtonDisabled) fCancel->SetEnabled();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for undoing last changes.

void TGeoCombiTransEditor::DoUndo()
{
   DoCancel();
   DoParameters();
   fCancel->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);
   fApply->SetEnabled(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for X.

void TGeoCombiTransEditor::DoDx()
{
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for Y.

void TGeoCombiTransEditor::DoDy()
{
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for Z.

void TGeoCombiTransEditor::DoDz()
{
   DoModified();
}

