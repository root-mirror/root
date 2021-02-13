// @(#)root/ged:$Id$
// Author: Ilka Antcheva   10/05/04

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TAttLineEditor                                                      //
//                                                                      //
//  Implements GUI for editing line attributes.                         //
//           color, line width, line style                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
//Begin_Html
/*
<img src="gif/TAttLineEditor.gif">
*/
//End_Html


#include "TAttLineEditor.h"
#include "TGColorSelect.h"
#include "TGComboBox.h"
#include "TColor.h"
#include "TGraph.h"
#include "TGLabel.h"
#include "TGNumberEntry.h"
#include "TCanvas.h"
#include "TROOT.h"

ClassImp(TAttLineEditor);

////////////////////////////////////////////////////////////////////////////////
/// Constructor of line attributes GUI.

TAttLineEditor::TAttLineEditor(const TGWindow *p, Int_t width,
                               Int_t height, UInt_t options, Pixel_t back)
   : TGedFrame(p, width, height, options | kVerticalFrame, back)
{
   enum ELineWid {
      kCOLOR,
      kLINE_WIDTH,
      kLINE_STYLE,
      kALPHA,
      kALPHAFIELD
   };

   fPriority = 1;
   fAttLine = 0;

   MakeTitle("Line");

   TGCompositeFrame *f2 = new TGCompositeFrame(this, 80, 20, kHorizontalFrame);
   AddFrame(f2, new TGLayoutHints(kLHintsTop, 1, 1, 0, 0));

   fColorSelect = new TGColorSelect(f2, 0, kCOLOR);
   f2->AddFrame(fColorSelect, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
   fColorSelect->Associate(this);

   fStyleCombo = new TGLineStyleComboBox(this, kLINE_STYLE);
   fStyleCombo->Resize(137, 20);
   AddFrame(fStyleCombo, new TGLayoutHints(kLHintsLeft, 3, 1, 1, 1));
   fStyleCombo->Associate(this);

   fWidthCombo = new TGLineWidthComboBox(f2, kLINE_WIDTH);
   fWidthCombo->Resize(90, 20);
   f2->AddFrame(fWidthCombo, new TGLayoutHints(kLHintsLeft, 3, 1, 1, 1));
   fWidthCombo->Associate(this);

   TGLabel *AlphaLabel = new TGLabel(this,"Opacity");
   AddFrame(AlphaLabel,
            new TGLayoutHints(kLHintsLeft | kLHintsCenterY));
   TGHorizontalFrame *f2a = new TGHorizontalFrame(this);
   fAlpha = new TGHSlider(f2a,100,kSlider2|kScaleNo,kALPHA);
   fAlpha->SetRange(0,1000);
   f2a->AddFrame(fAlpha,new TGLayoutHints(kLHintsLeft | kLHintsCenterY));
   fAlphaField = new TGNumberEntryField(f2a, kALPHAFIELD, 0,
                                        TGNumberFormat::kNESReal,
                                        TGNumberFormat::kNEANonNegative);
   fAlphaField->Resize(40,20);
   if (!TCanvas::SupportAlpha()) {
      fAlpha->SetEnabled(kFALSE);
      AlphaLabel->Disable(kTRUE);
      fAlphaField->SetEnabled(kFALSE);
   }
   f2a->AddFrame(fAlphaField,new TGLayoutHints(kLHintsLeft | kLHintsCenterY));
   AddFrame(f2a, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor of line editor.

TAttLineEditor::~TAttLineEditor()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Connect signals to slots.

void TAttLineEditor::ConnectSignals2Slots()
{
   fColorSelect->Connect("ColorSelected(Pixel_t)", "TAttLineEditor", this, "DoLineColor(Pixel_t)");
   fColorSelect->Connect("AlphaColorSelected(ULong_t)", "TAttLineEditor", this, "DoLineAlphaColor(ULong_t)");
   fStyleCombo->Connect("Selected(Int_t)", "TAttLineEditor", this, "DoLineStyle(Int_t)");
   fWidthCombo->Connect("Selected(Int_t)", "TAttLineEditor", this, "DoLineWidth(Int_t)");
   fAlpha->Connect("Released()","TAttLineEditor", this, "DoAlpha()");
   fAlpha->Connect("PositionChanged(Int_t)","TAttLineEditor", this, "DoLiveAlpha(Int_t)");
   fAlphaField->Connect("ReturnPressed()","TAttLineEditor", this, "DoAlphaField()");
   fAlpha->Connect("Pressed()","TAttLineEditor", this, "GetCurAlpha()");

   fInit = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Pick up the used line attributes.

void TAttLineEditor::SetModel(TObject* obj)
{
   TAttLine *attline = dynamic_cast<TAttLine*>(obj);
   if (!attline) return;

   fAttLine = attline;
   fAvoidSignal = kTRUE;

   fStyleCombo->Select(fAttLine->GetLineStyle());

   if (obj->InheritsFrom(TGraph::Class())) {
      fWidthCombo->Select(TMath::Abs(fAttLine->GetLineWidth()%100));
   } else {
      fWidthCombo->Select(fAttLine->GetLineWidth());
   }

   Color_t c = fAttLine->GetLineColor();
   Pixel_t p = TColor::Number2Pixel(c);
   fColorSelect->SetColor(p);

   if (fInit) ConnectSignals2Slots();

   fAvoidSignal = kFALSE;

   if (TColor *color = gROOT->GetColor(fAttLine->GetLineColor())) {
      fAlpha->SetPosition((Int_t)(color->GetAlpha()*1000));
      fAlphaField->SetNumber(color->GetAlpha());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the line color.

void TAttLineEditor::DoLineColor(Pixel_t color)
{
   if (fAvoidSignal) return;
   fAttLine->SetLineColor(TColor::GetColor(color));

   if (TColor *tcolor = gROOT->GetColor(TColor::GetColor(color))) {
      fAlpha->SetPosition((Int_t)(tcolor->GetAlpha()*1000));
      fAlphaField->SetNumber(tcolor->GetAlpha());
   }

   Update();
}


////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the color with alpha.

void TAttLineEditor::DoLineAlphaColor(ULong_t p)
{
   TColor *color = (TColor *)p;

   if (fAvoidSignal) return;
   fAttLine->SetLineColor(color->GetNumber());
   fAlpha->SetPosition((Int_t)(color->GetAlpha()*1000));
   fAlphaField->SetNumber(color->GetAlpha());

   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the line style.

void TAttLineEditor::DoLineStyle(Int_t style)
{
   if (fAvoidSignal) return;
   fAttLine->SetLineStyle(style);
   Update();
}


////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the line width.

void TAttLineEditor::DoLineWidth(Int_t width)
{
   if (fAvoidSignal) return;
   if (dynamic_cast<TGraph*>(fAttLine)) {
      Int_t graphLineWidth = 100*Int_t(fAttLine->GetLineWidth()/100);
      if (graphLineWidth >= 0) {
         fAttLine->SetLineWidth(graphLineWidth+width);
      } else {
         fAttLine->SetLineWidth(-(TMath::Abs(graphLineWidth)+width));
      }
   } else {
      fAttLine->SetLineWidth(width);
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the alpha value from the entry field.

void TAttLineEditor::DoAlphaField()
{
   if (fAvoidSignal) return;

   if (TColor *color = gROOT->GetColor(fAttLine->GetLineColor())) {
      color->SetAlpha((Float_t)fAlphaField->GetNumber());
      fAlpha->SetPosition((Int_t)fAlphaField->GetNumber()*1000);
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the alpha value

void TAttLineEditor::DoAlpha()
{
   if (fAvoidSignal) return;

   if (TColor *color = gROOT->GetColor(fAttLine->GetLineColor())) {
      color->SetAlpha((Float_t)fAlpha->GetPosition()/1000);
      fAlphaField->SetNumber((Float_t)fAlpha->GetPosition()/1000);
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set alpha value online.

void TAttLineEditor::DoLiveAlpha(Int_t a)
{
   if (fAvoidSignal) return;
   fAlphaField->SetNumber((Float_t)a/1000);

   if (TColor *color = gROOT->GetColor(fAttLine->GetLineColor())) {
      // In case the color is not transparent a new color is created.
      if (color->GetAlpha() == 1.) {
         fAttLine->SetLineColor(TColor::GetColorTransparent(color->GetNumber(),0.99));
      } else {
         color->SetAlpha((Float_t)a/1000);
      }
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to update alpha value on click on Slider

void TAttLineEditor::GetCurAlpha()
{
   if (fAvoidSignal) return;

   if (TColor *color = gROOT->GetColor(fAttLine->GetLineColor())) {
      fAlpha->SetPosition((Int_t)(color->GetAlpha()*1000));
      fAlphaField->SetNumber(color->GetAlpha());
   }
   Update();
}
