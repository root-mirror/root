// @(#):$Id$
// Author: M.Gheata

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoTrd2Editor
\ingroup Geometry_builder

Editor for a TGeoTrd2.

\image html geom_trd2_pic.png

\image html geom_trd2_ed.png

*/

#include "TGeoTrd2Editor.h"
#include "TGeoTabManager.h"
#include "TGeoTrd2.h"
#include "TGeoManager.h"
#include "TVirtualGeoPainter.h"
#include "TVirtualPad.h"
#include "TView.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"

ClassImp(TGeoTrd2Editor);

enum ETGeoTrd2Wid {
   kTRD2_NAME, kTRD2_X1, kTRD2_X2,  kTRD2_Y1, kTRD2_Y2, kTRD2_Z,
   kTRD2_APPLY, kTRD2_UNDO
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor for trd2 editor

TGeoTrd2Editor::TGeoTrd2Editor(const TGWindow *p, Int_t width,
                                   Int_t height, UInt_t options, Pixel_t back)
   : TGeoGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fShape   = 0;
   fDxi1 = fDxi2 = fDyi1 = fDyi2 = fDzi = 0.0;
   fNamei = "";
   fIsModified = kFALSE;
   fIsShapeEditable = kFALSE;

   // TextEntry for shape name
   MakeTitle("Name");
   fShapeName = new TGTextEntry(this, new TGTextBuffer(50), kTRD2_NAME);
   fShapeName->Resize(135, fShapeName->GetDefaultHeight());
   fShapeName->SetToolTipText("Enter the box name");
   fShapeName->Associate(this);
   AddFrame(fShapeName, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));

   TGTextEntry *nef;
   MakeTitle("Trd2 dimensions");
   TGCompositeFrame *compxyz = new TGCompositeFrame(this, 118, 30, kVerticalFrame | kRaisedFrame | kDoubleBorder);

   // Number entry for dx1
   TGCompositeFrame *f1 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f1->AddFrame(new TGLabel(f1, "DX1"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEDx1 = new TGNumberEntry(f1, 0., 5, kTRD2_X1);
   fEDx1->SetNumAttr(TGNumberFormat::kNEAPositive);
   nef = (TGTextEntry*)fEDx1->GetNumberEntry();
   nef->SetToolTipText("Enter the half-length in X1");
   fEDx1->Associate(this);
   f1->AddFrame(fEDx1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));

   // Number entry for dx2
   f1 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f1->AddFrame(new TGLabel(f1, "DX2"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEDx2 = new TGNumberEntry(f1, 0., 5, kTRD2_X2);
   fEDx2->SetNumAttr(TGNumberFormat::kNEAPositive);
   nef = (TGTextEntry*)fEDx2->GetNumberEntry();
   nef->SetToolTipText("Enter the  half-length in X2");
   fEDx2->Associate(this);
   f1->AddFrame(fEDx2, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));

   // Number entry for dy1
   TGCompositeFrame *f2 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f2->AddFrame(new TGLabel(f2, "DY1"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEDy1 = new TGNumberEntry(f2, 0., 5, kTRD2_Y1);
   fEDy1->SetNumAttr(TGNumberFormat::kNEAPositive);
   nef = (TGTextEntry*)fEDy1->GetNumberEntry();
   nef->SetToolTipText("Enter the half-length in Y1");
   fEDy1->Associate(this);
   f2->AddFrame(fEDy1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));
   compxyz->AddFrame(f2, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));

   // Number entry for dy2
   f2 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f2->AddFrame(new TGLabel(f2, "DY2"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEDy2 = new TGNumberEntry(f2, 0., 5, kTRD2_Y2);
   fEDy2->SetNumAttr(TGNumberFormat::kNEAPositive);
   nef = (TGTextEntry*)fEDy2->GetNumberEntry();
   nef->SetToolTipText("Enter the half-length in Y2");
   fEDy2->Associate(this);
   f2->AddFrame(fEDy2, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));
   compxyz->AddFrame(f2, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));

   // Number entry for dz
   TGCompositeFrame *f3 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f3->AddFrame(new TGLabel(f3, "DZ"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEDz = new TGNumberEntry(f3, 0., 5, kTRD2_Z);
   fEDz->SetNumAttr(TGNumberFormat::kNEAPositive);
   nef = (TGTextEntry*)fEDz->GetNumberEntry();
   nef->SetToolTipText("Enter the  half-length in Z");
   fEDz->Associate(this);
   f3->AddFrame(fEDz, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));
   compxyz->AddFrame(f3, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));

   compxyz->Resize(150,30);
   AddFrame(compxyz, new TGLayoutHints(kLHintsLeft, 6, 6, 4, 4));

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

TGeoTrd2Editor::~TGeoTrd2Editor()
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

void TGeoTrd2Editor::ConnectSignals2Slots()
{
   fApply->Connect("Clicked()", "TGeoTrd2Editor", this, "DoApply()");
   fUndo->Connect("Clicked()", "TGeoTrd2Editor", this, "DoUndo()");
   fShapeName->Connect("TextChanged(const char *)", "TGeoTrd2Editor", this, "DoModified()");
   fEDx1->Connect("ValueSet(Long_t)", "TGeoTrd2Editor", this, "DoDx1()");
   fEDx2->Connect("ValueSet(Long_t)", "TGeoTrd2Editor", this, "DoDx2()");
   fEDy1->Connect("ValueSet(Long_t)", "TGeoTrd2Editor", this, "DoDy1()");
   fEDy2->Connect("ValueSet(Long_t)", "TGeoTrd2Editor", this, "DoDy2()");
   fEDz->Connect("ValueSet(Long_t)", "TGeoTrd2Editor", this, "DoDz()");
   fEDx1->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTrd2Editor", this, "DoModified()");
   fEDx2->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTrd2Editor", this, "DoModified()");
   fEDy1->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTrd2Editor", this, "DoModified()");
   fEDy2->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTrd2Editor", this, "DoModified()");
   fEDz->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTrd2Editor", this, "DoModified()");
   fInit = kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Connect to the selected object.

void TGeoTrd2Editor::SetModel(TObject* obj)
{
   if (obj == 0 || (obj->IsA()!=TGeoTrd2::Class())) {
      SetActive(kFALSE);
      return;
   }
   fShape = (TGeoTrd2*)obj;
   fDxi1 = fShape->GetDx1();
   fDxi2 = fShape->GetDx2();
   fDyi1 = fShape->GetDy1();
   fDyi2 = fShape->GetDy2();
   fDzi = fShape->GetDz();
   const char *sname = fShape->GetName();
   if (!strcmp(sname, fShape->ClassName())) fShapeName->SetText("-no_name");
   else {
      fShapeName->SetText(sname);
      fNamei = sname;
   }
   fEDx1->SetNumber(fDxi1);
   fEDx2->SetNumber(fDxi2);
   fEDy1->SetNumber(fDyi1);
   fEDy2->SetNumber(fDyi2);
   fEDz->SetNumber(fDzi);
   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);

   if (fInit) ConnectSignals2Slots();
   SetActive();
}

////////////////////////////////////////////////////////////////////////////////
/// Check if shape drawing is delayed.

Bool_t TGeoTrd2Editor::IsDelayed() const
{
   return (fDelayed->GetState() == kButtonDown);
}

////////////////////////////////////////////////////////////////////////////////
/// Perform name change.

void TGeoTrd2Editor::DoName()
{
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for applying modifications.

void TGeoTrd2Editor::DoApply()
{
   const char *name = fShapeName->GetText();
   if (strcmp(name,fShape->GetName())) fShape->SetName(name);
   Double_t dx1 = fEDx1->GetNumber();
   Double_t dx2 = fEDx2->GetNumber();
   Double_t dy1 = fEDy1->GetNumber();
   Double_t dy2 = fEDy2->GetNumber();
   Double_t dz = fEDz->GetNumber();
   Double_t param[5];
   param[0] = dx1;
   param[1] = dx2;
   param[2] = dy1;
   param[3] = dy2;
   param[4] = dz;
   fShape->SetDimensions(param);
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

void TGeoTrd2Editor::DoModified()
{
   fApply->SetEnabled();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for undoing last operation.

void TGeoTrd2Editor::DoUndo()
{
   fEDx1->SetNumber(fDxi1);
   fEDx2->SetNumber(fDxi2);
   fEDy1->SetNumber(fDyi1);
   fEDy2->SetNumber(fDyi2);
   fEDz->SetNumber(fDzi);
   DoApply();
   fUndo->SetEnabled(kFALSE);
   fApply->SetEnabled(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for dx1.

void TGeoTrd2Editor::DoDx1()
{
   Double_t dx1 = fEDx1->GetNumber();
   Double_t dx2 = fEDx2->GetNumber();
   if (dx1<0) {
      dx1 = 0;
      fEDx1->SetNumber(dx1);
   }
   if (dx1<1.e-6 && dx2<1.e-6) {
      dx1 = 0.1;
      fEDx1->SetNumber(dx1);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for dx2.

void TGeoTrd2Editor::DoDx2()
{
   Double_t dx1 = fEDx1->GetNumber();
   Double_t dx2 = fEDx2->GetNumber();
   if (dx2<0) {
      dx2 = 0;
      fEDx2->SetNumber(dx2);
   }
   if (dx1<1.e-6 && dx2<1.e-6) {
      dx2 = 0.1;
      fEDx2->SetNumber(dx2);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for dy1.

void TGeoTrd2Editor::DoDy1()
{
   Double_t dy1 = fEDy1->GetNumber();
   Double_t dy2 = fEDy2->GetNumber();
   if (dy1<0) {
      dy1 = 0;
      fEDy1->SetNumber(dy1);
   }
   if (dy1<1.e-6 && dy2<1.e-6) {
      dy1 = 0.1;
      fEDy1->SetNumber(dy1);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for dy2.

void TGeoTrd2Editor::DoDy2()
{
   Double_t dy1 = fEDy1->GetNumber();
   Double_t dy2 = fEDy2->GetNumber();
   if (dy2<0) {
      dy2 = 0;
      fEDy2->SetNumber(dy2);
   }
   if (dy1<1.e-6 && dy2<1.e-6) {
      dy2 = 0.1;
      fEDy2->SetNumber(dy2);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for dz.

void TGeoTrd2Editor::DoDz()
{
   Double_t dz = fEDz->GetNumber();
   if (dz<=0) {
      dz = 0.1;
      fEDz->SetNumber(dz);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}


