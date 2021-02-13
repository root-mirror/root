// @(#)root/eve:$Id$
// Authors: Alja & Matevz Tadel 2008

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveTextEditor.h"
#include "TEveText.h"
#include "TEveGValuators.h"

#include "TGLFontManager.h"

#include "TGLabel.h"
#include "TGComboBox.h"
#include "TGTextEntry.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TG3DLine.h"


/** \class TEveTextEditor
\ingroup TEve
GUI editor for TEveText.
*/

ClassImp(TEveTextEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveTextEditor::TEveTextEditor(const TGWindow *p, Int_t width, Int_t height,
                               UInt_t options, Pixel_t back) :
   TGedFrame(p, width, height, options | kVerticalFrame, back),
   fM(0),

   fSize(0),
   fFile(0),
   fMode(0),
   fExtrude(0),

   fLighting(0),
   fAutoLighting(0)
{
   MakeTitle("TEveText");

   // Text entry
   fText = new TGTextEntry(this);
   fText->Resize(135, fText->GetDefaultHeight());
   AddFrame(fText, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));
   fText->Connect("TextChanged(const char *)", "TEveTextEditor", this, "DoText(const char *)");

   // Face Size combo
   fSize = MakeLabeledCombo("Size:");
   Int_t* fsp = &TGLFontManager::GetFontSizeArray()->front();
   Int_t  nums = TGLFontManager::GetFontSizeArray()->size();
   for(Int_t i= 0; i< nums; i++)
   {
      fSize->AddEntry(Form("%-2d", fsp[i]), fsp[i]);
   }
   fSize->Connect("Selected(Int_t)", "TEveTextEditor", this, "DoFontSize()");

   // Font File combo
   fFile = MakeLabeledCombo("File:");
   TObjArray* farr = TGLFontManager::GetFontFileArray();
   TIter next(farr);
   TObjString* os;
   Int_t cnt = 0;
   while ((os = (TObjString*) next()) != 0)
   {
      fFile->AddEntry(Form("%s", os->GetString().Data()), cnt);
      cnt++;
   }
   fFile->Connect("Selected(Int_t)", "TEveTextEditor", this, "DoFontFile()");

   // Mode combo
   fMode = MakeLabeledCombo("Mode:");
   fMode->AddEntry("Bitmap",  TGLFont::kBitmap);
   fMode->AddEntry("Pixmap",  TGLFont::kPixmap);
   fMode->AddEntry("Texture", TGLFont::kTexture);
   fMode->AddEntry("Outline", TGLFont::kOutline);
   fMode->AddEntry("Polygon", TGLFont::kPolygon);
   fMode->AddEntry("Extrude", TGLFont::kExtrude);
   fMode->Connect("Selected(Int_t)", "TEveTextEditor", this, "DoFontMode()");

   fExtrude = new TEveGValuator(this, "Depth:", 90, 0);
   fExtrude->SetLabelWidth(45);
   fExtrude->SetNELength(5);
   // fExtrude->SetShowSlider(kFALSE);
   fExtrude->Build();
   fExtrude->SetLimits(0.01, 10, 100, TGNumberFormat::kNESRealTwo);
   fExtrude->SetToolTip("Extrusion depth.");
   fExtrude->Connect("ValueSet(Double_t)", "TEveTextEditor", this, "DoExtrude()");
   AddFrame(fExtrude, new TGLayoutHints(kLHintsTop, 4, 1, 1, 1));

   // GLConfig
   TGCompositeFrame *f1 = new TGCompositeFrame(this, 145, 10, kHorizontalFrame | kFitWidth | kFixedWidth );
   f1->AddFrame(new TGLabel(f1, "GLConfig"), new TGLayoutHints(kLHintsLeft, 1, 1, 0, 0));
   f1->AddFrame(new TGHorizontal3DLine(f1), new TGLayoutHints(kLHintsExpandX, 5, 5, 7, 7));
   AddFrame(f1, new TGLayoutHints(kLHintsTop, 0, 0, 8, 0));

   TGCompositeFrame *alf = new TGCompositeFrame(this, 145, 10, kHorizontalFrame );
   fAutoLighting  = new TGCheckButton(alf, "AutoLighting");
   alf->AddFrame(fAutoLighting, new TGLayoutHints(kLHintsLeft, 1,2,0,0));
   fAutoLighting->Connect("Toggled(Bool_t)", "TEveTextEditor", this, "DoAutoLighting()");
   fLighting  = new TGCheckButton(alf, "Lighting");
   alf->AddFrame(fLighting, new TGLayoutHints(kLHintsLeft, 1,2,0,0));
   fLighting->Connect("Toggled(Bool_t)", "TEveTextEditor", this, "DoLighting()");
   AddFrame(alf, new TGLayoutHints(kLHintsTop, 0, 0, 0, 0));
}

////////////////////////////////////////////////////////////////////////////////
/// Helper function. Creates TGComboBox with fixed size TGLabel.

TGComboBox* TEveTextEditor::MakeLabeledCombo(const char* name)
{
   UInt_t labelW = 45;
   UInt_t labelH = 20;
   TGHorizontalFrame* hf = new TGHorizontalFrame(this);
   // label
   TGCompositeFrame *labfr = new TGHorizontalFrame(hf, labelW, labelH, kFixedSize);
   TGLabel* label = new TGLabel(labfr, name);
   labfr->AddFrame(label, new TGLayoutHints(kLHintsLeft  | kLHintsBottom));
   hf->AddFrame(labfr, new TGLayoutHints(kLHintsLeft));
   // combo
   TGLayoutHints*  clh =  new TGLayoutHints(kLHintsLeft, 0,0,0,0);
   TGComboBox* combo = new TGComboBox(hf);
   combo->Resize(90, 20);
   hf->AddFrame(combo, clh);

   AddFrame(hf, new TGLayoutHints(kLHintsTop, 4, 1, 1, 1));
   return combo;
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

void TEveTextEditor::SetModel(TObject* obj)
{
   fM = dynamic_cast<TEveText*>(obj);
   if (strcmp(fM->GetText(), fText->GetText()))
      fText->SetText(fM->GetText());

   fSize->Select(fM->GetFontSize(), kFALSE);
   fFile->Select(fM->GetFontFile(), kFALSE);

   // mode
   fMode->Select(fM->GetFontMode(), kFALSE);

   // lightning
   fAutoLighting->SetState(fM->GetAutoLighting() ? kButtonDown : kButtonUp);
   if (fM->GetAutoLighting()) {
      fLighting->SetDisabledAndSelected(fM->GetLighting() ? kButtonDown : kButtonUp);
   } else {
      fLighting->SetEnabled();
      fLighting->SetState(fM->GetLighting() ? kButtonDown : kButtonUp);
   }

   // extrude
   if (fM->GetFontMode() == TGLFont::kExtrude)
   {
      ShowFrame(fExtrude);
      fExtrude->SetValue(fM->GetExtrude());
   }
   else
   {
      HideFrame(fExtrude);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting text.

void TEveTextEditor::DoText(const char* /*txt*/)
{
   fM->SetText(fText->GetText());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting FTGL attributes.

void TEveTextEditor::DoFontSize()
{
   fM->SetFontSize(fSize->GetSelected(), kFALSE);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting FTGL attributes.

void TEveTextEditor::DoFontFile()
{
   fM->SetFontFile(fFile->GetSelected());
   Update();
}
////////////////////////////////////////////////////////////////////////////////
/// Slot for setting FTGL attributes.

void TEveTextEditor::DoFontMode()
{
   fM->SetFontMode(fMode->GetSelected());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for setting an extrude depth.

void TEveTextEditor::DoExtrude()
{
   fM->SetExtrude(fExtrude->GetValue());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for enabling/disabling defaults.

void TEveTextEditor::DoAutoLighting()
{
   fM->SetAutoLighting(fAutoLighting->IsOn());
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for enabling/disabling GL lighting.

void TEveTextEditor::DoLighting()
{
   fM->SetLighting(fLighting->IsOn());
   Update();
}
