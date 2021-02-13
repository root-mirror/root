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
//  TAttFillEditor                                                      //
//                                                                      //
//  Implements GUI for editing fill attributes.                         //                                             //
//             color and fill style                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
//Begin_Html
/*
<img src="gif/TAttFillEditor.gif">
*/
//End_Html

#include "TAttFillEditor.h"
#include "TGedPatternSelect.h"
#include "TGColorSelect.h"
#include "TColor.h"
#include "TGLabel.h"
#include "TGNumberEntry.h"
#include "TCanvas.h"
#include "TROOT.h"

ClassImp(TAttFillEditor);

enum EFillWid {
   kCOLOR,
   kPATTERN,
   kALPHA,
   kALPHAFIELD
};


////////////////////////////////////////////////////////////////////////////////
/// Constructor of fill attributes GUI.

TAttFillEditor::TAttFillEditor(const TGWindow *p, Int_t width,
                               Int_t height, UInt_t options, Pixel_t back)
   : TGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fPriority = 2;

   fAttFill = 0;

   MakeTitle("Fill");

   TGCompositeFrame *f2 = new TGCompositeFrame(this, 80, 20, kHorizontalFrame);
   fColorSelect = new TGColorSelect(f2, 0, kCOLOR);
   f2->AddFrame(fColorSelect, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
   fColorSelect->Associate(this);
   fPatternSelect = new TGedPatternSelect(f2, 1, kPATTERN);
   f2->AddFrame(fPatternSelect, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
   fPatternSelect->Associate(this);
   AddFrame(f2, new TGLayoutHints(kLHintsTop, 1, 1, 0, 0));

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
/// Destructor of fill editor.

TAttFillEditor::~TAttFillEditor()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Connect signals to slots.

void TAttFillEditor::ConnectSignals2Slots()
{
   fColorSelect->Connect("ColorSelected(Pixel_t)", "TAttFillEditor", this, "DoFillColor(Pixel_t)");
   fColorSelect->Connect("AlphaColorSelected(ULong_t)", "TAttFillEditor", this, "DoFillAlphaColor(ULong_t)");
   fPatternSelect->Connect("PatternSelected(Style_t)", "TAttFillEditor", this, "DoFillPattern(Style_t)");
   fAlpha->Connect("Released()","TAttFillEditor", this, "DoAlpha()");
   fAlpha->Connect("PositionChanged(Int_t)","TAttFillEditor", this, "DoLiveAlpha(Int_t)");
   fAlphaField->Connect("ReturnPressed()","TAttFillEditor", this, "DoAlphaField()");
   fAlpha->Connect("Pressed()","TAttFillEditor", this, "GetCurAlpha()");
   fInit = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Pick up the used fill attributes.

void TAttFillEditor::SetModel(TObject* obj)
{
   TAttFill *attfill = dynamic_cast<TAttFill *>(obj);
   if (!attfill) return;

   fAttFill = attfill;
   fAvoidSignal = kTRUE;

   Color_t c = fAttFill->GetFillColor();
   Pixel_t p = TColor::Number2Pixel(c);
   fColorSelect->SetColor(p, kFALSE);

   Style_t s = fAttFill->GetFillStyle();
   fPatternSelect->SetPattern(s, kFALSE);

   if (fInit) ConnectSignals2Slots();
   fAvoidSignal = kFALSE;

   if (TColor *color = gROOT->GetColor(fAttFill->GetFillColor())) {
      fAlpha->SetPosition((Int_t)(color->GetAlpha()*1000));
      fAlphaField->SetNumber(color->GetAlpha());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the fill area color with alpha.

void TAttFillEditor::DoFillAlphaColor(ULong_t p)
{
   TColor *color = (TColor *)p;

   if (fAvoidSignal) return;
   fAttFill->SetFillColor(color->GetNumber());
   fAlpha->SetPosition((Int_t)(color->GetAlpha()*1000));
   fAlphaField->SetNumber(color->GetAlpha());

   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the fill area color.

void TAttFillEditor::DoFillColor(Pixel_t color)
{
   if (fAvoidSignal) return;
   fAttFill->SetFillColor(TColor::GetColor(color));

   if (TColor *tcolor = gROOT->GetColor(TColor::GetColor(color))) {
      fAlpha->SetPosition((Int_t)(tcolor->GetAlpha()*1000));
      fAlphaField->SetNumber(tcolor->GetAlpha());
   }

   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the fill area pattern.

void TAttFillEditor::DoFillPattern(Style_t pattern)
{
   if (fAvoidSignal) return;
   fAttFill->SetFillStyle(pattern);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the alpha value from the entry field.

void TAttFillEditor::DoAlphaField()
{
   if (fAvoidSignal) return;

   if (TColor *color = gROOT->GetColor(fAttFill->GetFillColor())) {
      color->SetAlpha((Float_t)fAlphaField->GetNumber());
      fAlpha->SetPosition((Int_t)fAlphaField->GetNumber()*1000);
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the alpha value

void TAttFillEditor::DoAlpha()
{
   if (fAvoidSignal) return;

   if (TColor *color = gROOT->GetColor(fAttFill->GetFillColor())) {
      color->SetAlpha((Float_t)fAlpha->GetPosition()/1000);
      fAlphaField->SetNumber((Float_t)fAlpha->GetPosition()/1000);
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set alpha value online.

void TAttFillEditor::DoLiveAlpha(Int_t a)
{
   if (fAvoidSignal) return;
   fAlphaField->SetNumber((Float_t)a/1000);

   if (TColor *color = gROOT->GetColor(fAttFill->GetFillColor()))  {
      // In case the color is not transparent a new color is created.
      if (color->GetAlpha() == 1.) {
         fAttFill->SetFillColor(TColor::GetColorTransparent(color->GetNumber(),0.99));
      } else {
         color->SetAlpha((Float_t)a/1000);
      }
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to update alpha value on click on Slider

void TAttFillEditor::GetCurAlpha()
{
   if (fAvoidSignal) return;

   if (TColor *color = gROOT->GetColor(fAttFill->GetFillColor())) {
      fAlpha->SetPosition((Int_t)(color->GetAlpha()*1000));
      fAlphaField->SetNumber(color->GetAlpha());
   }
   Update();
}
