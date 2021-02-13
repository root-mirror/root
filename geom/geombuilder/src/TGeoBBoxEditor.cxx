// @(#):$Id$
// Author: M.Gheata

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoBBoxEditor
\ingroup Geometry_builder

Editor for a TGeoBBox.

\image html geom_box_pic.png

\image html geom_box_ed.png

*/

#include "TGeoBBoxEditor.h"
#include "TGeoTabManager.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TVirtualGeoPainter.h"
#include "TVirtualPad.h"
#include "TView.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"

ClassImp(TGeoBBoxEditor);

enum ETGeoBBoxWid {
   kBOX_NAME, kBOX_X, kBOX_Y, kBOX_Z,
   kBOX_OX, kBOX_OY, kBOX_OZ,
   kBOX_APPLY, kBOX_CANCEL, kBOX_UNDO
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor for volume editor.

TGeoBBoxEditor::TGeoBBoxEditor(const TGWindow *p, Int_t width,
                               Int_t height, UInt_t options, Pixel_t back)
   : TGeoGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fShape   = 0;
   fDxi = fDyi = fDzi = 0.0;
   memset(fOrigi, 0, 3*sizeof(Double_t));
   fNamei = "";
   fIsModified = kFALSE;
   fIsShapeEditable = kFALSE;

   // TextEntry for shape name
   MakeTitle("Name");
   fShapeName = new TGTextEntry(this, new TGTextBuffer(50), kBOX_NAME);
   fShapeName->Resize(135, fShapeName->GetDefaultHeight());
   fShapeName->SetToolTipText("Enter the box name");
   fShapeName->Associate(this);
   AddFrame(fShapeName, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));

   TGTextEntry *nef;
   MakeTitle("Box half-lengths");
   TGCompositeFrame *compxyz = new TGCompositeFrame(this, 118, 30, kVerticalFrame | kRaisedFrame | kDoubleBorder);
   // Number entry for dx
   TGCompositeFrame *f1 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f1->AddFrame(new TGLabel(f1, "DX"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fBoxDx = new TGNumberEntry(f1, 0., 5, kBOX_X);
   fBoxDx->SetNumAttr(TGNumberFormat::kNEAPositive);
   nef = (TGTextEntry*)fBoxDx->GetNumberEntry();
   nef->SetToolTipText("Enter the box half-length in X");
   fBoxDx->Associate(this);
   f1->AddFrame(fBoxDx, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));

   // Number entry for dy
   TGCompositeFrame *f2 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f2->AddFrame(new TGLabel(f2, "DY"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fBoxDy = new TGNumberEntry(f2, 0., 5, kBOX_Y);
   fBoxDy->SetNumAttr(TGNumberFormat::kNEAPositive);
   nef = (TGTextEntry*)fBoxDy->GetNumberEntry();
   nef->SetToolTipText("Enter the box half-length in Y");
   fBoxDy->Associate(this);
   f2->AddFrame(fBoxDy, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));
   compxyz->AddFrame(f2, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));

   // Number entry for dx
   TGCompositeFrame *f3 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f3->AddFrame(new TGLabel(f3, "DZ"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fBoxDz = new TGNumberEntry(f3, 0., 5, kBOX_Z);
   fBoxDz->SetNumAttr(TGNumberFormat::kNEAPositive);
   nef = (TGTextEntry*)fBoxDz->GetNumberEntry();
   nef->SetToolTipText("Enter the box half-length in Z");
   fBoxDz->Associate(this);
   f3->AddFrame(fBoxDz, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));
   compxyz->AddFrame(f3, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));

   compxyz->Resize(150,30);
   AddFrame(compxyz, new TGLayoutHints(kLHintsLeft, 6, 6, 4, 4));

   MakeTitle("Box origin");
   compxyz = new TGCompositeFrame(this, 118, 30, kVerticalFrame | kRaisedFrame | kDoubleBorder);
   // Number entry for dx
   f1 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                             kLHintsExpandX | kFixedWidth | kOwnBackground);
   f1->AddFrame(new TGLabel(f1, "OX"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fBoxOx = new TGNumberEntry(f1, 0., 5, kBOX_OX);
   nef = (TGTextEntry*)fBoxOx->GetNumberEntry();
   nef->SetToolTipText("Enter the box origin X coordinate");
   fBoxOx->Associate(this);
   f1->AddFrame(fBoxOx, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));

   // Number entry for dy
   f2 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                             kLHintsExpandX | kFixedWidth | kOwnBackground);
   f2->AddFrame(new TGLabel(f2, "OY"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fBoxOy = new TGNumberEntry(f2, 0., 5, kBOX_OY);
   nef = (TGTextEntry*)fBoxOy->GetNumberEntry();
   nef->SetToolTipText("Enter the box origin Y coordinate");
   fBoxOy->Associate(this);
   f2->AddFrame(fBoxOy, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));
   compxyz->AddFrame(f2, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));

   // Number entry for dx
   f3 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                             kLHintsExpandX | kFixedWidth | kOwnBackground);
   f3->AddFrame(new TGLabel(f3, "OZ"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fBoxOz = new TGNumberEntry(f3, 0., 5, kBOX_OZ);
   nef = (TGTextEntry*)fBoxOz->GetNumberEntry();
   nef->SetToolTipText("Enter the box origin Z coordinate");
   fBoxOz->Associate(this);
   f3->AddFrame(fBoxOz, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));
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

TGeoBBoxEditor::~TGeoBBoxEditor()
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

void TGeoBBoxEditor::ConnectSignals2Slots()
{
   fApply->Connect("Clicked()", "TGeoBBoxEditor", this, "DoApply()");
   fUndo->Connect("Clicked()", "TGeoBBoxEditor", this, "DoUndo()");
   fShapeName->Connect("TextChanged(const char *)", "TGeoBBoxEditor", this, "DoModified()");
   fBoxDx->Connect("ValueSet(Long_t)", "TGeoBBoxEditor", this, "DoDx()");
   fBoxDy->Connect("ValueSet(Long_t)", "TGeoBBoxEditor", this, "DoDy()");
   fBoxDz->Connect("ValueSet(Long_t)", "TGeoBBoxEditor", this, "DoDz()");
   fBoxDx->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoBBoxEditor", this, "DoModified()");
   fBoxDy->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoBBoxEditor", this, "DoModified()");
   fBoxDz->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoBBoxEditor", this, "DoModified()");
   fBoxOx->Connect("ValueSet(Long_t)", "TGeoBBoxEditor", this, "DoOx()");
   fBoxOy->Connect("ValueSet(Long_t)", "TGeoBBoxEditor", this, "DoOy()");
   fBoxOz->Connect("ValueSet(Long_t)", "TGeoBBoxEditor", this, "DoOz()");
   fBoxOx->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoBBoxEditor", this, "DoModified()");
   fBoxOy->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoBBoxEditor", this, "DoModified()");
   fBoxOz->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoBBoxEditor", this, "DoModified()");
   fInit = kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Update editor for a new selected box.

void TGeoBBoxEditor::SetModel(TObject* obj)
{
   if (obj == 0 || (obj->IsA()!=TGeoBBox::Class())) {
      SetActive(kFALSE);
      return;
   }
   fShape = (TGeoBBox*)obj;
   fDxi = fShape->GetDX();
   fDyi = fShape->GetDY();
   fDzi = fShape->GetDZ();
   memcpy(fOrigi, fShape->GetOrigin(), 3*sizeof(Double_t));
   const char *sname = fShape->GetName();
   if (!strcmp(sname, fShape->ClassName())) fShapeName->SetText("-no_name");
   else {
      fShapeName->SetText(sname);
      fNamei = sname;
   }
   fBoxDx->SetNumber(fDxi);
   fBoxDy->SetNumber(fDyi);
   fBoxDz->SetNumber(fDzi);
   fBoxOx->SetNumber(fOrigi[0]);
   fBoxOy->SetNumber(fOrigi[1]);
   fBoxOz->SetNumber(fOrigi[2]);
   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);

   if (fInit) ConnectSignals2Slots();
   SetActive();
}

////////////////////////////////////////////////////////////////////////////////
/// Check if shape drawing is delayed.

Bool_t TGeoBBoxEditor::IsDelayed() const
{
   return (fDelayed->GetState() == kButtonDown);
}

////////////////////////////////////////////////////////////////////////////////
///Slot for name.

void TGeoBBoxEditor::DoName()
{
   DoModified();
}

////////////////////////////////////////////////////////////////////////////////
///Slot for applying current parameters.

void TGeoBBoxEditor::DoApply()
{
   const char *name = fShapeName->GetText();
   if (strcmp(name,fShape->GetName())) fShape->SetName(name);
   Double_t dx = fBoxDx->GetNumber();
   Double_t dy = fBoxDy->GetNumber();
   Double_t dz = fBoxDz->GetNumber();
   Double_t orig[3];
   orig[0] = fBoxOx->GetNumber();
   orig[1] = fBoxOy->GetNumber();
   orig[2] = fBoxOz->GetNumber();
   fShape->SetBoxDimensions(dx, dy, dz, orig);
   fUndo->SetEnabled();
   fApply->SetEnabled(kFALSE);
   if (fPad) {
      if (gGeoManager && gGeoManager->GetPainter() && gGeoManager->GetPainter()->IsPaintingShape()) {
         TView *view = fPad->GetView();
         if (!view) {
            fShape->Draw();
            fPad->GetView()->ShowAxis();
         } else {
            const Double_t *origin = fShape->GetOrigin();
            view->SetRange(origin[0]-fShape->GetDX(), origin[1]-fShape->GetDY(), origin[2]-fShape->GetDZ(),
                           origin[0]+fShape->GetDX(), origin[1]+fShape->GetDY(), origin[2]+fShape->GetDZ());
            Update();
         }
      } else Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
///Slot for modifying current parameters.

void TGeoBBoxEditor::DoModified()
{
   fApply->SetEnabled();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for undoing last operation.

void TGeoBBoxEditor::DoUndo()
{
   fBoxDx->SetNumber(fDxi);
   fBoxDy->SetNumber(fDyi);
   fBoxDz->SetNumber(fDzi);
   fBoxOx->SetNumber(fOrigi[0]);
   fBoxOy->SetNumber(fOrigi[1]);
   fBoxOz->SetNumber(fOrigi[2]);
   DoApply();
   fUndo->SetEnabled(kFALSE);
   fApply->SetEnabled(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
///Slot for Dx modification.

void TGeoBBoxEditor::DoDx()
{
   Double_t dx = fBoxDx->GetNumber();
   if (dx<=0) {
      dx=0.1;
      fBoxDx->SetNumber(dx);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
///Slot for Dy modification.

void TGeoBBoxEditor::DoDy()
{
   Double_t dy = fBoxDy->GetNumber();
   if (dy<=0) {
      dy=0.1;
      fBoxDy->SetNumber(dy);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
///Slot for Dz modification.

void TGeoBBoxEditor::DoDz()
{
   Double_t dz = fBoxDz->GetNumber();
   if (dz<=0) {
      dz=0.1;
      fBoxDz->SetNumber(dz);
   }
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
///Slot for Ox modification.

void TGeoBBoxEditor::DoOx()
{
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
///Slot for Oy modification.

void TGeoBBoxEditor::DoOy()
{
   DoModified();
   if (!IsDelayed()) DoApply();
}

////////////////////////////////////////////////////////////////////////////////
///Slot for Oz modification.

void TGeoBBoxEditor::DoOz()
{
   DoModified();
   if (!IsDelayed()) DoApply();
}


