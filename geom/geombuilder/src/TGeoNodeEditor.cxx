// @(#):$Id: b94274447e46ae2cf700c562e57c3c38c6612b74 $
// Author: M.Gheata

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoNodeEditor
\ingroup Geometry_builder

Editor class for TGeoNode objects.

*/

#include "TGeoNodeEditor.h"
#include "TGedEditor.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TGTab.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"
#include "TGeoTabManager.h"

ClassImp(TGeoNodeEditor);

enum ETGeoNodeWid {
   kNODE_NAME, kNODE_ID, kNODE_VOLSEL, kNODE_MVOLSEL,
   kNODE_MATRIX, kNODE_EDIT_VOL, kNODE_EDIT_MATRIX
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor for node editor

TGeoNodeEditor::TGeoNodeEditor(const TGWindow *p, Int_t width,
                               Int_t height, UInt_t options, Pixel_t back)
   : TGeoGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fNode   = 0;
   fIsEditable = kTRUE;
   Pixel_t color;

   // TextEntry for medium name
   TGTextEntry *nef;
   MakeTitle("Name");
   TGCompositeFrame *f1 = new TGCompositeFrame(this, 140, 30, kHorizontalFrame | kRaisedFrame);
   fNodeName = new TGTextEntry(f1, new TGTextBuffer(50), kNODE_NAME);
   fNodeName->Resize(100, fNodeName->GetDefaultHeight());
   fNodeName->SetToolTipText("Enter the node name");
   fNodeName->Associate(this);
   f1->AddFrame(fNodeName, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));
   f1->AddFrame(new TGLabel(f1, "ID"), new TGLayoutHints(kLHintsLeft, 1, 1, 6, 0));
   fNodeNumber = new TGNumberEntry(f1, 0., 1, kNODE_ID);
   nef = (TGTextEntry*)fNodeNumber->GetNumberEntry();
   nef->SetToolTipText("Enter the node copy number");
   fNodeNumber->Associate(this);
   f1->AddFrame(fNodeNumber, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 3, 3, 2, 5));


// Mother volume selection
   MakeTitle("Mother volume");
   f1 = new TGCompositeFrame(this, 155, 30, kHorizontalFrame | kFixedWidth);
   fSelectedMother = 0;
   fLSelMother = new TGLabel(f1, "Select mother");
   gClient->GetColorByName("#0000ff", color);
   fLSelMother->SetTextColor(color);
   fLSelMother->ChangeOptions(kSunkenFrame | kDoubleBorder);
   f1->AddFrame(fLSelMother, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 1, 1, 2, 2));
   fBSelMother = new TGPictureButton(f1, fClient->GetPicture("rootdb_t.xpm"), kNODE_MVOLSEL);
   fBSelMother->SetToolTipText("Select one of the existing volumes");
   fBSelMother->Associate(this);
   f1->AddFrame(fBSelMother, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
   fEditMother = new TGTextButton(f1, "Edit");
   f1->AddFrame(fEditMother, new TGLayoutHints(kLHintsRight, 1, 1, 1, 1));
   fEditMother->Associate(this);
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 2));

// Volume selection
   MakeTitle("Volume");
   f1 = new TGCompositeFrame(this, 155, 30, kHorizontalFrame | kFixedWidth);
   fSelectedVolume = 0;
   fLSelVolume = new TGLabel(f1, "Select volume");
   gClient->GetColorByName("#0000ff", color);
   fLSelVolume->SetTextColor(color);
   fLSelVolume->ChangeOptions(kSunkenFrame | kDoubleBorder);
   f1->AddFrame(fLSelVolume, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 1, 1, 2, 2));
   fBSelVolume = new TGPictureButton(f1, fClient->GetPicture("rootdb_t.xpm"), kNODE_VOLSEL);
   fBSelVolume->SetToolTipText("Select one of the existing volumes");
   fBSelVolume->Associate(this);
   f1->AddFrame(fBSelVolume, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
   fEditVolume = new TGTextButton(f1, "Edit");
   f1->AddFrame(fEditVolume, new TGLayoutHints(kLHintsRight, 1, 1, 1, 1));
   fEditVolume->Associate(this);
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 2));

// Matrix selection
   MakeTitle("Matrix");
   f1 = new TGCompositeFrame(this, 155, 30, kHorizontalFrame | kFixedWidth);
   fSelectedMatrix = 0;
   fLSelMatrix = new TGLabel(f1, "Select matrix");
   gClient->GetColorByName("#0000ff", color);
   fLSelMatrix->SetTextColor(color);
   fLSelMatrix->ChangeOptions(kSunkenFrame | kDoubleBorder);
   f1->AddFrame(fLSelMatrix, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 1, 1, 2, 2));
   fBSelMatrix = new TGPictureButton(f1, fClient->GetPicture("rootdb_t.xpm"), kNODE_MATRIX);
   fBSelMatrix->SetToolTipText("Select one of the existing matrices");
   fBSelMatrix->Associate(this);
   f1->AddFrame(fBSelMatrix, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
   fEditMatrix = new TGTextButton(f1, "Edit");
   f1->AddFrame(fEditMatrix, new TGLayoutHints(kLHintsRight, 1, 1, 1, 1));
   fEditMatrix->Associate(this);
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 2));

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

TGeoNodeEditor::~TGeoNodeEditor()
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

void TGeoNodeEditor::ConnectSignals2Slots()
{
   fBSelMother->Connect("Clicked()", "TGeoNodeEditor", this, "DoSelectMother()");
   fBSelVolume->Connect("Clicked()", "TGeoNodeEditor", this, "DoSelectVolume()");
   fBSelMatrix->Connect("Clicked()", "TGeoNodeEditor", this, "DoSelectMatrix()");
   fApply->Connect("Clicked()", "TGeoNodeEditor", this, "DoApply()");
   fUndo->Connect("Clicked()", "TGeoNodeEditor", this, "DoUndo()");
   fEditMother->Connect("Clicked()", "TGeoNodeEditor", this, "DoEditMother()");
   fEditVolume->Connect("Clicked()", "TGeoNodeEditor", this, "DoEditVolume()");
   fEditMatrix->Connect("Clicked()", "TGeoNodeEditor", this, "DoEditMatrix()");
   fNodeName->Connect("TextChanged(const char *)", "TGeoNodeEditor", this, "DoNodeName()");
   fInit = kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Connect to a editable object.

void TGeoNodeEditor::SetModel(TObject* obj)
{
   if (obj == 0 || !obj->InheritsFrom(TGeoNode::Class())) {
      SetActive(kFALSE);
      return;
   }
   fNode = (TGeoNode*)obj;
   const char *sname = fNode->GetName();
   fNodeName->SetText(sname);

   fNodeNumber->SetNumber(fNode->GetNumber());

   fSelectedMother = fNode->GetMotherVolume();
   if (fSelectedMother) fLSelMother->SetText(fSelectedMother->GetName());
   fSelectedVolume = fNode->GetVolume();
   if (fSelectedVolume) fLSelVolume->SetText(fSelectedVolume->GetName());
   fSelectedMatrix = fNode->GetMatrix();
   if (fSelectedMatrix) fLSelMatrix->SetText(fSelectedMatrix->GetName());

   fApply->SetEnabled(kFALSE);
   fUndo->SetEnabled(kFALSE);

   if (fInit) ConnectSignals2Slots();
   SetActive();
}

////////////////////////////////////////////////////////////////////////////////
/// Select the mother volume.

void TGeoNodeEditor::DoSelectMother()
{
   TGeoVolume *vol = fSelectedMother;
   new TGeoVolumeDialog(fBSelMother, gClient->GetRoot(), 200,300);
   fSelectedMother = (TGeoVolume*)TGeoVolumeDialog::GetSelected();
   if (fSelectedMother) fLSelMother->SetText(fSelectedMother->GetName());
   else fSelectedMother = vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Select the volume.

void TGeoNodeEditor::DoSelectVolume()
{
   TGeoVolume *vol = fSelectedVolume;
   new TGeoVolumeDialog(fBSelVolume, gClient->GetRoot(), 200,300);
   fSelectedVolume = (TGeoVolume*)TGeoVolumeDialog::GetSelected();
   if (fSelectedVolume) fLSelVolume->SetText(fSelectedVolume->GetName());
   else fSelectedVolume = vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Select the matrix.

void TGeoNodeEditor::DoSelectMatrix()
{
   TGeoMatrix *matrix = fSelectedMatrix;
   new TGeoMatrixDialog(fBSelMatrix, gClient->GetRoot(), 200,300);
   fSelectedMatrix = (TGeoMatrix*)TGeoMatrixDialog::GetSelected();
   if (fSelectedMatrix) fLSelMatrix->SetText(fSelectedMatrix->GetName());
   else fSelectedMatrix = matrix;
}

////////////////////////////////////////////////////////////////////////////////
/// Edit the mother volume.

void TGeoNodeEditor::DoEditMother()
{
   if (!fSelectedMother) {
      fTabMgr->SetVolTabEnabled(kFALSE);
      return;
   }
   fTabMgr->SetVolTabEnabled();
   fTabMgr->GetVolumeEditor(fSelectedMother);
   fTabMgr->SetTab();
   fSelectedMother->Draw();
}

////////////////////////////////////////////////////////////////////////////////
/// Edit selected volume.

void TGeoNodeEditor::DoEditVolume()
{
   if (!fSelectedVolume) {
      fTabMgr->SetVolTabEnabled(kFALSE);
      return;
   }
   fTabMgr->SetVolTabEnabled();
   fTabMgr->GetVolumeEditor(fSelectedVolume);
   fTabMgr->SetTab();
   fSelectedVolume->Draw();
}

////////////////////////////////////////////////////////////////////////////////
/// Edit selected material.

void TGeoNodeEditor::DoEditMatrix()
{
   if (!fSelectedMatrix) return;
   fTabMgr->GetMatrixEditor(fSelectedMatrix);
}

////////////////////////////////////////////////////////////////////////////////
/// Change node name.

void TGeoNodeEditor::DoNodeName()
{
   const char *name = fNodeName->GetText();
   if (!name[0] || !strcmp(name, fNode->GetName())) return;
   fNode->SetName(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Change node copy number

void TGeoNodeEditor::DoNodeNumber()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for applying modifications.

void TGeoNodeEditor::DoApply()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for undoing last operation.

void TGeoNodeEditor::DoUndo()
{
}

