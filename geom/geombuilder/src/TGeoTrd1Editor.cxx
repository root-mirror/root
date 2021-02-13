// @(#):$Id$
// Author: M.Gheata

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoTrd1Editor
\ingroup Geometry_builder

Editor for a TGeoTrd1.

\image html geom_trd1_pic.png

\image html geom_trd1_ed.png

*/

#include "TGeoTrd1Editor.h"
#include "TGeoTabManager.h"
#include "TGeoTrd1.h"
#include "TGeoManager.h"
#include "TVirtualGeoPainter.h"
#include "TVirtualPad.h"
#include "TView.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"

ClassImp(TGeoTrd1Editor);

enum ETGeoTrd1Wid {
   kTRD1_NAME, kTRD1_X1, kTRD1_X2,  kTRD1_Y, kTRD1_Z,
   kTRD1_APPLY, kTRD1_UNDO
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor for trd1 editor

TGeoTrd1Editor::TGeoTrd1Editor(const TGWindow *p, Int_t width,
                                   Int_t height, UInt_t options, Pixel_t back)
   : TGeoGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fShape   = 0;
   fDxi1 = fDxi2 = fDyi = fDzi = 0.0;
   fNamei = "";
   fIsModified = kFALSE;
   fIsShapeEditable = kFALSE;

   // TextEntry for shape name
   MakeTitle("Name");
   fShapeName = new TGTextEntry(this, new TGTextBuffer(50), kTRD1_NAME);
   fShapeName->Resize(135, fShapeName->GetDefaultHeight());
   fShapeName->SetToolTipText("Enter the box name");
   fShapeName->Associate(this);
   AddFrame(fShapeName, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));

   TGTextEntry *nef;
   MakeTitle("Trd1 dimensions");
   TGCompositeFrame *compxyz = new TGCompositeFrame(this, 118, 30, kVerticalFrame | kRaisedFrame | kDoubleBorder);

   // Number entry for dx1
   TGCompositeFrame *f1 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f1->AddFrame(new TGLabel(f1, "DX1"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEDx1 = new TGNumberEntry(f1, 0., 5, kTRD1_X1);
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
   fEDx2 = new TGNumberEntry(f1, 0., 5, kTRD1_X2);
   fEDx2->SetNumAttr(TGNumberFormat::kNEAPositive);
   nef = (TGTextEntry*)fEDx2->GetNumberEntry();
   nef->SetToolTipText("Enter the  half-length in X2");
   fEDx2->Associate(this);
   f1->AddFrame(fEDx2, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));

   // Number entry for dy
   TGCompositeFrame *f2 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f2->AddFrame(new TGLabel(f2, "DY"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEDy = new TGNumberEntry(f2, 0., 5, kTRD1_Y);
   fEDy->SetNumAttr(TGNumberFormat::kNEAPositive);
   nef = (TGTextEntry*)fEDy->GetNumberEntry();
   nef->SetToolTipText("Enter the half-length in Y");
   fEDy->Associate(this);
   f2->AddFrame(fEDy, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));
   compxyz->AddFrame(f2, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));

   // Number entry for dz
   TGCompositeFrame *f3 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f3->AddFrame(new TGLabel(f3, "DZ"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fEDz = new TGNumberEntry(f3, 0., 5, kTRD1_Z);
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
/// Destructor

TGeoTrd1Editor::~TGeoTrd1Editor()
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

void TGeoTrd1Editor::ConnectSignals2Slots()
{
   fApply->Connect("Clicked()", "TGeoTrd1Editor", this, "DoApply()");
   fUndo->Connect("Clicked()", "TGeoTrd1Editor", this, "DoUndo()");
   fShapeName->Connect("TextChanged(const char *)", "TGeoTrd1Editor", this, "DoModified()");
   fEDx1->Connect("ValueSet(Long_t)", "TGeoTrd1Editor", this, "DoDx1()");
   fEDx2->Connect("ValueSet(Long_t)", "TGeoTrd1Editor", this, "DoDx2()");
   fEDy->Connect("ValueSet(Long_t)", "TGeoTrd1Editor", this, "DoDy()");
   fEDz->Connect("ValueSet(Long_t)", "TGeoTrd1Editor", this, "DoDz()");
   fEDx1->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTrd1Editor", this, "DoModified()");
   fEDx2->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTrd1Editor", this, "DoModified()");
   fEDy->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTrd1Editor", this, "DoModified()");
   fEDz->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoTrd1Editor", this, "DoModified()");
   fInit = kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Connect to the selected object.

void TGeoTrd1Editor::SetModel(TObject* obj)
{
   if (obj == 0 || (obj->IsA()!=TGeoTrd1::Class())) {
      SetActive(kFALSE);
      return;
   }
   fShape = (TGeoTrd1*)obj;
   fDxi1 = fShape->GetDx1();
   fDxi2 = fShape->GetDx2();
   fDyi = fShape->GetDy();
   fDzi = fShape->GetDz();
   const char *sname = fShape->GetName();
   if (!strcmp(sname, fShape->ClassName())) fShapeName->SetText("-no_name");
   else {
      fShapeName->SetText(sname);
      fNamei = sname;
   }
   fEDx1->SetNumber(fDxi1);
   fEDx2->SetNumber(fDxi2);
   fEDy->SetNumber(fDyi);
   fEDz->SetNumber(fDzi);
   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);


   if (fInit) ConnectSignals2Slots();
   SetActive();
}

////////////////////////////////////////////////////////////////////////////////
/// Check if shape drawing is delayed.

Bool_t TGeoTrd1Editor::IsDelayed() const
{
   return (fDelayed->GetState() == kButtonDown);
}

////////////////////////////////////////////////////////////////////////////////
/// Perform name change.

void TGeoTrd1Editor::DoName()
{
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for applying modifications.

void TGeoTrd1Editor::DoApply()
{
   const char *name = fShapeName->GetText();
   if (strcmp(name,fShape->GetName())) fShape->SetName(name);
   Double_t dx1 = fEDx1->GetNumber();
   Double_t dx2 = fEDx2->GetNumber();
   Double_t dy = fEDy->GetNumber();
   Double_t dz = fEDz->GetNumber();
   Double_t param[4];
   param[0] = dx1;
   param[1] = dx2;
   param[2] = dy;
   param[3] = dz;
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

void TGeoTrd1Editor::DoModified()
{
   fApply->SetEnabled();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for undoing last operation.

void TGeoTrd1Editor::DoUndo()
{
   fEDx1->SetNumber(fDxi1);
   fEDx2->SetNumber(fDxi2);
   fEDy->SetNumber(fDyi);
   fEDz->SetNumber(fDzi);
   DoApply();
   fUndo->SetEnabled(kFALSE);
   fApply->SetEnabled(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for dx1.

void TGeoTrd1Editor::DoDx1()
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

void TGeoTrd1Editor::DoDx2()
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
/// Slot for dy.

void TGeoTrd1Editor::DoDy()
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
/// Slot for dz.

void TGeoTrd1Editor::DoDz()
{
   Double_t dz = fEDz->GetNumber();
   if (dz<=0) {
      dz = 0.1;
      fEDz->SetNumber(dz);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}


