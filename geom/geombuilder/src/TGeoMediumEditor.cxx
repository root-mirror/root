// @(#):$Id: c816a00a89512fcc6cd4a75fb1343c76ebaa3c24 $
// Author: M.Gheata

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoMediumEditor
\ingroup Geometry_builder

Editor class for TGeo tracking media.

*/

#include "TGeoMediumEditor.h"
#include "TGeoTabManager.h"
#include "TGeoManager.h"
#include "TGeoMedium.h"
#include "TGeoMaterial.h"
#include "TGTab.h"
#include "TGComboBox.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"
#include "TG3DLine.h"

ClassImp(TGeoMediumEditor);

enum ETGeoMediumWid {
   kMED_NAME, kMED_ID, kMED_MATSEL,
   kMED_SENS, kMED_FLDOPT, kMED_EDIT_MAT,
   kMED_FIELDM, kMED_TMAX, kMED_STEMAX,
   kMED_DEEMAX, kMED_EPSIL, kMED_STMIN,
   kMED_APPLY, kMED_CANCEL, kMED_UNDO
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor for medium editor

TGeoMediumEditor::TGeoMediumEditor(const TGWindow *p, Int_t width,
                                   Int_t height, UInt_t options, Pixel_t back)
   : TGeoGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fMedium   = 0;
   fIsEditable = kFALSE;
   fIsModified = kFALSE;
   Pixel_t color;
   TGLabel *label;

   // TextEntry for medium name
   MakeTitle("Name");
   fMedName = new TGTextEntry(this, "", kMED_NAME);
   fMedName->SetDefaultSize(135, fMedName->GetDefaultHeight());
   fMedName->SetToolTipText("Enter the medium name");
   fMedName->Associate(this);
   AddFrame(fMedName, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 3, 1, 2, 2));

   TGTextEntry *nef;

// Composite frame for medium ID and sensitivity
   TGCompositeFrame *f1 = new TGCompositeFrame(this, 120, 30, kHorizontalFrame | kRaisedFrame);
   f1->AddFrame(new TGLabel(f1, "ID"), new TGLayoutHints(kLHintsLeft, 4, 1, 6, 0));
   fMedId = new TGNumberEntry(f1, 0., 1, kMED_ID);
   nef = (TGTextEntry*)fMedId->GetNumberEntry();
   nef->SetToolTipText("Enter the medium ID");
   fMedId->Associate(this);
   f1->AddFrame(fMedId, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 2, 2, 4, 4));
   fMedSensitive = new TGCheckButton(f1, "&Sens", kMED_SENS);
   fMedSensitive->Associate(this);
   f1->AddFrame(fMedSensitive, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));
   AddFrame(f1, new TGLayoutHints(kLHintsLeft, 3, 3, 2, 2));

   // Current material
   f1 = new TGCompositeFrame(this, 145, 10, kHorizontalFrame | kFixedWidth | kOwnBackground);
   f1->AddFrame(label = new TGLabel(f1, "Current material"), new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   f1->AddFrame(new TGHorizontal3DLine(f1), new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
   gClient->GetColorByName("#ff0000", color);
   label->SetTextColor(color);
   AddFrame(f1, new TGLayoutHints(kLHintsTop, 0, 0, 2, 0));
   f1 = new TGCompositeFrame(this, 155, 30, kHorizontalFrame);
   fSelectedMaterial = 0;
   fLSelMaterial = new TGLabel(f1, "Select material");
   gClient->GetColorByName("#0000ff", color);
   fLSelMaterial->SetTextColor(color);
   fLSelMaterial->ChangeOptions(kChildFrame | kSunkenFrame | kDoubleBorder);
   f1->AddFrame(fLSelMaterial, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 1, 1, 2, 2));
   fBSelMaterial = new TGPictureButton(f1, fClient->GetPicture("rootdb_t.xpm"), kMED_MATSEL);
   fBSelMaterial->SetToolTipText("Replace with one of the existing materials");
   fBSelMaterial->Associate(this);
   f1->AddFrame(fBSelMaterial, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
   fEditMaterial = new TGTextButton(f1, "Edit");
   f1->AddFrame(fEditMaterial, new TGLayoutHints(kLHintsLeft, 1, 1, 2, 2));
   fEditMaterial->SetToolTipText("Edit selected material");
   fEditMaterial->Associate(this);
   AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 2, 2, 0, 0));

// Combo box for magnetic field option
   f1 = new TGCompositeFrame(this, 145, 10, kHorizontalFrame | kLHintsExpandX | kFixedWidth | kOwnBackground);
   f1->AddFrame(label = new TGLabel(f1, "Mag. field option"), new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   f1->AddFrame(new TGHorizontal3DLine(f1), new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
   gClient->GetColorByName("#ff0000", color);
   label->SetTextColor(color);
   AddFrame(f1, new TGLayoutHints(kLHintsTop, 0, 0, 2, 0));
   fMagfldOption = new TGComboBox(this, kMED_FLDOPT);
   fMagfldOption->Resize(135, fMedName->GetDefaultHeight());
   AddFrame(fMagfldOption, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 2));


// Number entries for other settings
   f1 = new TGCompositeFrame(this, 145, 10, kHorizontalFrame | kFixedWidth | kOwnBackground);
   f1->AddFrame(label = new TGLabel(f1, "Medium cuts"), new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   f1->AddFrame(new TGHorizontal3DLine(f1), new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
   gClient->GetColorByName("#ff0000", color);
   label->SetTextColor(color);
   AddFrame(f1, new TGLayoutHints(kLHintsTop, 0, 0, 2, 0));
   TGCompositeFrame *compxyz = new TGCompositeFrame(this, 130, 30, kVerticalFrame | kRaisedFrame | kDoubleBorder);

   // Number entry for fieldm
   f1 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f1->AddFrame(new TGLabel(f1, "FIELDM"), new TGLayoutHints(kLHintsLeft, 1, 1, 4, 0));
   fMedFieldm = new TGNumberEntry(f1, 0., 5, kMED_FIELDM);
   nef = (TGTextEntry*)fMedFieldm->GetNumberEntry();
   nef->SetToolTipText("Maximum magnetic field [kilogauss]");
   fMedFieldm->Associate(this);
   fMedFieldm->Resize(90, fMedFieldm->GetDefaultHeight());
   f1->AddFrame(fMedFieldm, new TGLayoutHints(kLHintsRight | kFixedWidth , 2, 2, 2, 2));
   compxyz->AddFrame(f1, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 1, 1));

   // Number entry for tmaxfd
   TGCompositeFrame *f2 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f2->AddFrame(new TGLabel(f2, "TMAXFD"), new TGLayoutHints(kLHintsLeft, 1, 1, 4, 0));
   fMedTmaxfd = new TGNumberEntry(f2, 0., 5, kMED_TMAX);
   nef = (TGTextEntry*)fMedTmaxfd->GetNumberEntry();
   nef->SetToolTipText("Maximum angle per step due to field [deg]");
   fMedTmaxfd->Associate(this);
   fMedTmaxfd->Resize(90, fMedTmaxfd->GetDefaultHeight());
   f2->AddFrame(fMedTmaxfd, new TGLayoutHints(kLHintsRight | kFixedWidth , 2, 2, 2, 2));
   compxyz->AddFrame(f2, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 1, 1));

   // Number entry for stemax
   TGCompositeFrame *f3 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f3->AddFrame(new TGLabel(f3, "STEMAX"), new TGLayoutHints(kLHintsLeft, 1, 1, 4, 0));
   fMedStemax = new TGNumberEntry(f3, 0., 5, kMED_STEMAX);
   nef = (TGTextEntry*)fMedStemax->GetNumberEntry();
   nef->SetToolTipText("Maximum step allowed [cm]");
   fMedStemax->Associate(this);
   fMedStemax->Resize(90, fMedStemax->GetDefaultHeight());
   f3->AddFrame(fMedStemax, new TGLayoutHints(kLHintsRight | kFixedWidth , 2, 2, 2, 2));
   compxyz->AddFrame(f3, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 1, 1));

   // Number entry for deemax
   TGCompositeFrame *f4 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f4->AddFrame(new TGLabel(f4, "DEEMAX"), new TGLayoutHints(kLHintsLeft, 1, 1, 4, 0));
   fMedDeemax = new TGNumberEntry(f4, 0., 5, kMED_DEEMAX);
   nef = (TGTextEntry*)fMedDeemax->GetNumberEntry();
   nef->SetToolTipText("Maximum fraction of energy lost in a step");
   fMedDeemax->Associate(this);
   fMedDeemax->Resize(90, fMedDeemax->GetDefaultHeight());
   f4->AddFrame(fMedDeemax, new TGLayoutHints(kLHintsRight | kFixedWidth , 2, 2, 2, 2));
   compxyz->AddFrame(f4, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 1, 1));

   // Number entry for epsil
   TGCompositeFrame *f5 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f5->AddFrame(new TGLabel(f5, "EPSIL"), new TGLayoutHints(kLHintsLeft, 1, 1, 4, 0));
   fMedEpsil = new TGNumberEntry(f5, 0., 5, kMED_EPSIL);
   nef = (TGTextEntry*)fMedEpsil->GetNumberEntry();
   nef->SetToolTipText("Tracking precision [cm]");
   fMedEpsil->Associate(this);
   fMedEpsil->Resize(90, fMedEpsil->GetDefaultHeight());
   f5->AddFrame(fMedEpsil, new TGLayoutHints(kLHintsRight | kFixedWidth , 2, 2, 2, 2));
   compxyz->AddFrame(f5, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 1, 1));

   // Number entry for stmin
   TGCompositeFrame *f6 = new TGCompositeFrame(compxyz, 118, 10, kHorizontalFrame |
                                 kLHintsExpandX | kFixedWidth | kOwnBackground);
   f6->AddFrame(new TGLabel(f6, "STMIN"), new TGLayoutHints(kLHintsLeft, 1, 1, 4, 0));
   fMedStmin = new TGNumberEntry(f6, 0., 5, kMED_STMIN);
   nef = (TGTextEntry*)fMedStmin->GetNumberEntry();
   nef->SetToolTipText("Minimum step due to continuous processes [cm]");
   fMedStmin->Associate(this);
   fMedStmin->Resize(90, fMedStmin->GetDefaultHeight());
   f6->AddFrame(fMedStmin, new TGLayoutHints(kLHintsRight | kFixedWidth , 2, 2, 2, 2));
   compxyz->AddFrame(f6, new TGLayoutHints(kLHintsLeft | kLHintsExpandX , 2, 2, 1, 1));

   compxyz->Resize(160,50);
   AddFrame(compxyz, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));

   // Buttons
   TGCompositeFrame *f23 = new TGCompositeFrame(this, 118, 20, kHorizontalFrame | kSunkenFrame | kDoubleBorder);
   fApply = new TGTextButton(f23, "&Apply");
   f23->AddFrame(fApply, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));
   fApply->Associate(this);
   fUndo = new TGTextButton(f23, " &Undo ");
   f23->AddFrame(fUndo, new TGLayoutHints(kLHintsRight , 2, 2, 4, 4));
   fUndo->Associate(this);
   AddFrame(f23,  new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoMediumEditor::~TGeoMediumEditor()
{
   TGFrameElement *el;
   TIter next(GetList());
   while ((el = (TGFrameElement *)next())) {
      if (el->fFrame->IsA() == TGCompositeFrame::Class()  ||
          el->fFrame->IsA() == TGHorizontalFrame::Class() ||
          el->fFrame->IsA() == TGVerticalFrame::Class())
         TGeoTabManager::Cleanup((TGCompositeFrame*)el->fFrame);
   }
   Cleanup();
}

////////////////////////////////////////////////////////////////////////////////
/// Connect signals to slots.

void TGeoMediumEditor::ConnectSignals2Slots()
{
   fApply->Connect("Clicked()", "TGeoMediumEditor", this, "DoApply()");
   fUndo->Connect("Clicked()", "TGeoMediumEditor", this, "DoUndo()");
   fMedName->Connect("TextChanged(const char *)", "TGeoMediumEditor", this, "DoMedName()");
   fBSelMaterial->Connect("Clicked()", "TGeoMediumEditor", this, "DoSelectMaterial()");
   fEditMaterial->Connect("Clicked()", "TGeoMediumEditor", this, "DoEditMaterial()");
   fMedId->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoMediumEditor", this, "DoMedId()");
   fMedTmaxfd->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoMediumEditor", this, "DoTmaxfd()");
   fMedStemax->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoMediumEditor", this, "DoStemax()");
   fMedDeemax->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoMediumEditor", this, "DoDeemax()");
   fMedEpsil->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoMediumEditor", this, "DoEpsil()");
   fMedStmin->GetNumberEntry()->Connect("TextChanged(const char *)", "TGeoMediumEditor", this, "DoStmin()");
   fMedSensitive->Connect("Clicked()", "TGeoMediumEditor", this, "DoToggleSensitive()");
   fMagfldOption->Connect("Selected(Int_t)", "TGeoMediumEditor", this, "DoMagfldSelect(Int_t)");
   fInit = kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Connect to the selected object.

void TGeoMediumEditor::SetModel(TObject* obj)
{
   if (obj == 0 || !(obj->IsA()==TGeoMedium::Class())) {
      SetActive(kFALSE);
      return;
   }
   fMedium = (TGeoMedium*)obj;
   const char *sname = fMedium->GetName();
   if (!strcmp(sname, fMedium->ClassName())) fMedName->SetText("");
   else fMedName->SetText(sname);

   fMedId->SetNumber(fMedium->GetId());
   Int_t isvol = (Int_t)fMedium->GetParam(0);
   fMedSensitive->SetState((isvol==0)?kButtonUp:kButtonDown);

   fSelectedMaterial = fMedium->GetMaterial();
   if (fSelectedMaterial) fLSelMaterial->SetText(fSelectedMaterial->GetName());

   if (!fMagfldOption->GetNumberOfEntries()) {
      fMagfldOption->AddEntry("No field", 0);
      fMagfldOption->AddEntry("User decision", 1);
      fMagfldOption->AddEntry("Runge-Kutta", 2);
      fMagfldOption->AddEntry("Helix", 3);
      fMagfldOption->AddEntry("Helix3", 4);
      fMagfldOption->AddEntry("Unknown option", 5);
   }
   Int_t ifld = (Int_t)fMedium->GetParam(1);
   switch (ifld) {
      case 0:
         fMagfldOption->Select(0);
         break;
      case -1:
         fMagfldOption->Select(1);
         break;
      case 1:
         fMagfldOption->Select(2);
         break;
      case 2:
         fMagfldOption->Select(3);
         break;
      case 3:
         fMagfldOption->Select(4);
         break;
      default:
         fMagfldOption->Select(5);
         break;
   }

   fMedFieldm->SetNumber(fMedium->GetParam(2));
   fMedTmaxfd->SetNumber(fMedium->GetParam(3));
   fMedStemax->SetNumber(fMedium->GetParam(4));
   fMedDeemax->SetNumber(fMedium->GetParam(5));
   fMedEpsil->SetNumber(fMedium->GetParam(6));
   fMedStmin->SetNumber(fMedium->GetParam(7));

   fUndo->SetEnabled(kFALSE);
   fIsModified = kFALSE;


   if (fInit) ConnectSignals2Slots();
   SetActive();
}

////////////////////////////////////////////////////////////////////////////////
/// Edit selected material.

void TGeoMediumEditor::DoEditMaterial()
{
   fTabMgr->GetMaterialEditor(fMedium->GetMaterial());
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for medium name.

void TGeoMediumEditor::DoMedName()
{
   const char *name = fMedName->GetText();
   if (!name[0] || !strcmp(name, fMedium->GetName())) return;
   fMedium->SetName(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for medium id.

void TGeoMediumEditor::DoMedId()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Select the material component.

void TGeoMediumEditor::DoSelectMaterial()
{
   TGeoMaterial *material = fSelectedMaterial;
   new TGeoMaterialDialog(fBSelMaterial, gClient->GetRoot(), 200,300);
   fSelectedMaterial = (TGeoMaterial*)TGeoMaterialDialog::GetSelected();
   if (fSelectedMaterial) fLSelMaterial->SetText(fSelectedMaterial->GetName());
   else fSelectedMaterial = material;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for sensitivity.

void TGeoMediumEditor::DoToggleSensitive()
{
   fIsModified = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for mag. field.

void TGeoMediumEditor::DoMagfldSelect(Int_t)
{
   fIsModified = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for max field.

void TGeoMediumEditor::DoFieldm()
{
   fIsModified = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for tmaxfd.

void TGeoMediumEditor::DoTmaxfd()
{
   fIsModified = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for the max allowed step.

void TGeoMediumEditor::DoStemax()
{
   fIsModified = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for the maximum allowed dedx.

void TGeoMediumEditor::DoDeemax()
{
   fIsModified = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for tracking precision.

void TGeoMediumEditor::DoEpsil()
{
   fIsModified = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for min. step.

void TGeoMediumEditor::DoStmin()
{
   fIsModified = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for applying modifications.

void TGeoMediumEditor::DoApply()
{
   if (!fIsModified) return;
   Double_t isvol = (fMedSensitive->IsOn())?1:0;
   Double_t ifield = fMagfldOption->GetSelected();
   if (ifield>0) {
      ifield -= 1.;
      if (ifield < 1.) ifield -= 1.;
   }
   Double_t fieldm = fMedFieldm->GetNumber();
   Double_t tmaxfd = fMedTmaxfd->GetNumber();
   Double_t stemax = fMedStemax->GetNumber();
   Double_t deemax = fMedDeemax->GetNumber();
   Double_t epsil = fMedEpsil->GetNumber();
   Double_t stmin = fMedStmin->GetNumber();

   fMedium->SetParam(0,isvol);
   fMedium->SetParam(1,ifield);
   fMedium->SetParam(2,fieldm);
   fMedium->SetParam(3,tmaxfd);
   fMedium->SetParam(4,stemax);
   fMedium->SetParam(5,deemax);
   fMedium->SetParam(6,epsil);
   fMedium->SetParam(7,stmin);
   if (strcmp(fMedium->GetName(), fMedName->GetText())) fMedium->SetName(fMedName->GetText());
   if (fMedium->GetId() != fMedId->GetIntNumber()) fMedium->SetId(fMedId->GetIntNumber());
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for undoing last operation.

void TGeoMediumEditor::DoUndo()
{
}

