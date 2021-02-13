// @(#)root/ged:$Id$
// Author: Ilka Antcheva   11/05/04

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TAttMarkerEditor                                                    //
//                                                                      //
//  Implements GUI for editing marker attributes.                       //
//            color, style and size                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
//Begin_Html
/*
<img src="gif/TAttMarkerEditor.gif">
*/
//End_Html


#include "TAttMarkerEditor.h"
#include "TAttMarker.h"
#include "TGedMarkerSelect.h"
#include "TGColorSelect.h"
#include "TGNumberEntry.h"
#include "TColor.h"
#include "TGLabel.h"
#include "TCanvas.h"
#include "TROOT.h"


ClassImp(TAttMarkerEditor);

enum EMarkerWid {
   kCOLOR,
   kMARKER,
   kMARKER_SIZE,
   kALPHA,
   kALPHAFIELD
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor of marker attributes GUI.

TAttMarkerEditor::TAttMarkerEditor(const TGWindow *p, Int_t width,
                                   Int_t height,UInt_t options, Pixel_t back)
   : TGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fAttMarker = 0;
   fSizeForText = kFALSE;

   MakeTitle("Marker");

   TGCompositeFrame *f2 = new TGCompositeFrame(this, 80, 20, kHorizontalFrame);
   fColorSelect = new TGColorSelect(f2, 0, kCOLOR);
   f2->AddFrame(fColorSelect, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
   fColorSelect->Associate(this);

   fMarkerType = new TGedMarkerSelect(f2, 1, kMARKER);
   f2->AddFrame(fMarkerType, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
   fMarkerType->Associate(this);

   fMarkerSize = new TGNumberEntry(f2, 0., 4, kMARKER_SIZE,
                                   TGNumberFormat::kNESRealOne,
                                   TGNumberFormat::kNEANonNegative,
                                   TGNumberFormat::kNELLimitMinMax, 0.2, 5.0);
   fMarkerSize->GetNumberEntry()->SetToolTipText("Set marker size");
   f2->AddFrame(fMarkerSize, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
   fMarkerSize->Associate(this);
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
/// Destructor of marker editor.

TAttMarkerEditor::~TAttMarkerEditor()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Connect signals to slots.

void TAttMarkerEditor::ConnectSignals2Slots()
{
   fColorSelect->Connect("ColorSelected(Pixel_t)", "TAttMarkerEditor", this, "DoMarkerColor(Pixel_t)");
   fColorSelect->Connect("AlphaColorSelected(ULong_t)", "TAttMarkerEditor", this, "DoMarkerAlphaColor(ULong_t)");
   fMarkerType->Connect("MarkerSelected(Style_t)", "TAttMarkerEditor", this, "DoMarkerStyle(Style_t)");
   fMarkerSize->Connect("ValueSet(Long_t)", "TAttMarkerEditor", this, "DoMarkerSize()");
   (fMarkerSize->GetNumberEntry())->Connect("ReturnPressed()", "TAttMarkerEditor", this, "DoMarkerSize()");
   fAlpha->Connect("Released()","TAttMarkerEditor", this, "DoAlpha()");
   fAlpha->Connect("PositionChanged(Int_t)","TAttMarkerEditor", this, "DoLiveAlpha(Int_t)");
   fAlphaField->Connect("ReturnPressed()","TAttMarkerEditor", this, "DoAlphaField()");
   fAlpha->Connect("Pressed()","TAttMarkerEditor", this, "GetCurAlpha()");
   fInit = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Pick up the values of used marker attributes.

void TAttMarkerEditor::SetModel(TObject* obj)
{
   fAvoidSignal = kTRUE;

   fAttMarker = dynamic_cast<TAttMarker *>(obj);
   if (!fAttMarker) return;

   TString str = GetDrawOption();
   str.ToUpper();
   if (obj->InheritsFrom("TH2") && str.Contains("TEXT")) {
      fSizeForText = kTRUE;
   } else {
      fSizeForText = kFALSE;
   }
   Style_t marker = fAttMarker->GetMarkerStyle();
   if ((marker==1 || marker==6 || marker==7) && !fSizeForText) {
      fMarkerSize->SetNumber(1.);
      fMarkerSize->SetState(kFALSE);
   } else {
      Float_t s = fAttMarker->GetMarkerSize();
      fMarkerSize->SetState(kTRUE);
      fMarkerSize->SetNumber(s);
   }
   fMarkerType->SetMarkerStyle(marker);

   Color_t c = fAttMarker->GetMarkerColor();
   Pixel_t p = TColor::Number2Pixel(c);
   fColorSelect->SetColor(p);

   if (fInit) ConnectSignals2Slots();
   fAvoidSignal = kFALSE;

   if (TColor *color = gROOT->GetColor(fAttMarker->GetMarkerColor())) {
      fAlpha->SetPosition((Int_t)(color->GetAlpha()*1000));
      fAlphaField->SetNumber(color->GetAlpha());
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the marker color.

void TAttMarkerEditor::DoMarkerColor(Pixel_t color)
{
   if (fAvoidSignal) return;
   fAttMarker->SetMarkerColor(TColor::GetColor(color));

   if (TColor *tcolor = gROOT->GetColor(TColor::GetColor(color))) {
      fAlpha->SetPosition((Int_t)(tcolor->GetAlpha()*1000));
      fAlphaField->SetNumber(tcolor->GetAlpha());
   }

   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the color with alpha.

void TAttMarkerEditor::DoMarkerAlphaColor(ULong_t p)
{
   TColor *color = (TColor *)p;

   if (fAvoidSignal) return;
   fAttMarker->SetMarkerColor(color->GetNumber());
   fAlpha->SetPosition((Int_t)(color->GetAlpha()*1000));
   fAlphaField->SetNumber(color->GetAlpha());

   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the marker type.

void TAttMarkerEditor::DoMarkerStyle(Style_t marker)
{
   if (fAvoidSignal) return;
   if ((marker==1 || marker==6 || marker==7) && !fSizeForText) {
      fMarkerSize->SetNumber(1.);
      fMarkerSize->SetState(kFALSE);
   } else
      fMarkerSize->SetState(kTRUE);

   fAttMarker->SetMarkerStyle(marker);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the marker size.

void TAttMarkerEditor::DoMarkerSize()
{
   if (fAvoidSignal) return;
   Style_t marker = fAttMarker->GetMarkerStyle();
   if ((marker==1 || marker==6 || marker==7) && !fSizeForText) {
      fMarkerSize->SetNumber(1.);
      fMarkerSize->SetState(kFALSE);
   } else
      fMarkerSize->SetState(kTRUE);
   Float_t size = fMarkerSize->GetNumber();
   fAttMarker->SetMarkerSize(size);
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the alpha value from the entry field.

void TAttMarkerEditor::DoAlphaField()
{
   if (fAvoidSignal) return;

   if (TColor *color = gROOT->GetColor(fAttMarker->GetMarkerColor())) {
      color->SetAlpha((Float_t)fAlphaField->GetNumber());
      fAlpha->SetPosition((Int_t)fAlphaField->GetNumber()*1000);
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the alpha value

void TAttMarkerEditor::DoAlpha()
{
   if (fAvoidSignal) return;

   if (TColor *color = gROOT->GetColor(fAttMarker->GetMarkerColor())) {
      color->SetAlpha((Float_t)fAlpha->GetPosition()/1000);
      fAlphaField->SetNumber((Float_t)fAlpha->GetPosition()/1000);
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set alpha value online.

void TAttMarkerEditor::DoLiveAlpha(Int_t a)
{
   if (fAvoidSignal) return;
   fAlphaField->SetNumber((Float_t)a/1000);

   if (TColor *color = gROOT->GetColor(fAttMarker->GetMarkerColor())) {
      // In case the color is not transparent a new color is created.
      if (color->GetAlpha() == 1.) {
         fAttMarker->SetMarkerColor(TColor::GetColorTransparent(color->GetNumber(),0.99));
      } else {
         color->SetAlpha((Float_t)a/1000);
      }
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to update alpha value on click on Slider

void TAttMarkerEditor::GetCurAlpha()
{
   if (fAvoidSignal) return;

   if (TColor *color = gROOT->GetColor(fAttMarker->GetMarkerColor())) {
      fAlpha->SetPosition((Int_t)(color->GetAlpha()*1000));
      fAlphaField->SetNumber(color->GetAlpha());
   }
   Update();
}
