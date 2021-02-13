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
//  TAttTextEditor                                                      //
//                                                                      //
//  Implements GUI for editing text attributes.                         //
//          color, font, font size and alignment                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
//Begin_Html
/*
<img src="gif/TAttTextEditor.gif">
*/
//End_Html


#include "TAttTextEditor.h"
#include "TGedEditor.h"
#include "TGColorSelect.h"
#include "TGComboBox.h"
#include "TColor.h"
#include "TPaveLabel.h"
#include "TGLabel.h"
#include "TGNumberEntry.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "snprintf.h"

ClassImp(TAttTextEditor);

enum ETextWid {
   kCOLOR,
   kFONT_SIZE,
   kFONT_STYLE,
   kFONT_ALIGN,
   kALPHA,
   kALPHAFIELD
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor of text attributes GUI.

TAttTextEditor::TAttTextEditor(const TGWindow *p, Int_t width,
                               Int_t height, UInt_t options, Pixel_t back)
   : TGedFrame(p, width, height, options | kVerticalFrame, back)
{
   fPriority = 3;

   fAttText = 0;

   MakeTitle("Text");

   TGCompositeFrame *f2 = new TGCompositeFrame(this, 80, 20, kHorizontalFrame);
   fColorSelect = new TGColorSelect(f2, 0, kCOLOR);
   f2->AddFrame(fColorSelect, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
   fColorSelect->Associate(this);
   fSizeCombo = BuildFontSizeComboBox(f2, kFONT_SIZE);
   f2->AddFrame(fSizeCombo, new TGLayoutHints(kLHintsLeft, 3, 1, 1, 1));
   fSizeCombo->Resize(91, 20);
   fSizeCombo->Associate(this);
   AddFrame(f2, new TGLayoutHints(kLHintsTop, 1, 1, 0, 0));
   fTypeCombo = new TGFontTypeComboBox(this, kFONT_STYLE);
   fTypeCombo->Resize(137, 20);
   AddFrame(fTypeCombo, new TGLayoutHints(kLHintsLeft, 3, 1, 1, 1));
   fAlignCombo = BuildTextAlignComboBox(this, kFONT_ALIGN);
   fAlignCombo->Resize(137, 20);
   AddFrame(fAlignCombo, new TGLayoutHints(kLHintsLeft, 3, 1, 1, 1));

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
/// Destructor of text editor.

TAttTextEditor::~TAttTextEditor()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Connect signals to slots.

void TAttTextEditor::ConnectSignals2Slots()
{
   fAlpha->Connect("Released()","TAttTextEditor", this, "DoAlpha()");
   fAlpha->Connect("PositionChanged(Int_t)","TAttTextEditor", this, "DoLiveAlpha(Int_t)");
   fAlphaField->Connect("ReturnPressed()","TAttTextEditor", this, "DoAlphaField()");
   fAlpha->Connect("Pressed()","TAttTextEditor", this, "GetCurAlpha()");
   fColorSelect->Connect("ColorSelected(Pixel_t)", "TAttTextEditor", this, "DoTextColor(Pixel_t)");
   fColorSelect->Connect("AlphaColorSelected(ULong_t)", "TAttTextEditor", this, "DoTextAlphaColor(ULong_t)");
   fInit = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Pick up the values of used text attributes.

void TAttTextEditor::SetModel(TObject* obj)
{
   TAttText *atttext = dynamic_cast<TAttText *>(obj);
   if (!atttext) return;

   fAttText = atttext;
   fAvoidSignal = kTRUE;

   fTypeCombo->Select(fAttText->GetTextFont() / 10);

   Float_t s = fAttText->GetTextSize();
   Float_t dy;

   if (obj->InheritsFrom(TPaveLabel::Class())) {
      TBox *pl = (TBox*)obj;
      dy = s * (pl->GetY2() - pl->GetY1());
   }
   else
      dy = s * (fGedEditor->GetPad()->GetY2() - fGedEditor->GetPad()->GetY1());

   Int_t size = fGedEditor->GetPad()->YtoPixel(0.0) - fGedEditor->GetPad()->YtoPixel(dy);
   if (size > 50) size = 50;
   if (size < 0)  size = 0;
   fSizeCombo->Select(size, kFALSE);

   fAlignCombo->Select(fAttText->GetTextAlign(), kFALSE);

   Color_t c = fAttText->GetTextColor();
   Pixel_t p = TColor::Number2Pixel(c);
   fColorSelect->SetColor(p, kFALSE);

   if (fInit) ConnectSignals2Slots();
   fAvoidSignal = kFALSE;

   if (TColor *color = gROOT->GetColor(fAttText->GetTextColor())) {
      fAlpha->SetPosition((Int_t)color->GetAlpha()*1000);
      fAlphaField->SetNumber(color->GetAlpha());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the marker color.

void TAttTextEditor::DoTextColor(Pixel_t color)
{
   if (fAvoidSignal) return;
   fAttText->SetTextColor(TColor::GetColor(color));

   if (TColor *tcolor = gROOT->GetColor(TColor::GetColor(color))) {
      fAlpha->SetPosition((Int_t)(tcolor->GetAlpha()*1000));
      fAlphaField->SetNumber(tcolor->GetAlpha());
   }

   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot connected to the color with alpha.

void TAttTextEditor::DoTextAlphaColor(ULong_t p)
{
   TColor *color = (TColor *)p;

   if (fAvoidSignal) return;
   fAttText->SetTextColor(color->GetNumber());
   fAlpha->SetPosition((Int_t)(color->GetAlpha()*1000));
   fAlphaField->SetNumber(color->GetAlpha());

   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Process message.

Bool_t TAttTextEditor::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
   if(!fGedEditor || !fGedEditor->GetModel()) return kTRUE;

   Bool_t b = kFALSE;

   if (GET_MSG(msg) == kC_COLORSEL && GET_SUBMSG(msg) == kCOL_SELCHANGED) {
      if (parm1 != 0) fAttText->SetTextColor(TColor::GetColor(parm2));
      b = kTRUE;
      // SendMessage(fMsgWindow, msg, parm1, parm2);
   }

   if (GET_MSG(msg) == kC_COMMAND && GET_SUBMSG(msg) == kCM_COMBOBOX) {

      if (parm1 == kFONT_SIZE) {
         TVirtualPad* pad = fGedEditor->GetPad();
         Float_t dy = pad->AbsPixeltoY(0) - pad->AbsPixeltoY(parm2);
         Float_t textSize;

         if (fGedEditor->GetModel()->InheritsFrom(TPaveLabel::Class())) {
            TBox *pl = (TBox*)fGedEditor->GetModel();
            textSize = dy/(pl->GetY2() - pl->GetY1());
         }
         else
            textSize = dy/(pad->GetY2() - pad->GetY1());

         fAttText->SetTextSize(textSize);
         b = kTRUE;
      } else if (parm1 == kFONT_STYLE) {
         Int_t fontPrec = fAttText->GetTextFont()%10;
         fAttText->SetTextFont(parm2 * 10 + fontPrec);
         b = kTRUE;
      } else if (parm1 == kFONT_ALIGN) {
         fAttText->SetTextAlign(parm2);
         b = kTRUE;
      }
   }

   if (b && !fAvoidSignal) Update();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create text size combo box.

TGComboBox* TAttTextEditor::BuildFontSizeComboBox(TGFrame* parent, Int_t id)
{
   char a[100];
   TGComboBox *c = new TGComboBox(parent, id);

   c->AddEntry("Default", 0);
   for (int i = 1; i <= 50; i++) {
      snprintf(a, 99, "%d", i);
      c->AddEntry(a, i);
   }

   return c;
}

////////////////////////////////////////////////////////////////////////////////
/// Create text align combo box.

TGComboBox* TAttTextEditor::BuildTextAlignComboBox(TGFrame* parent, Int_t id)
{
   TGComboBox *c = new TGComboBox(parent, id);

   c->AddEntry("11 Bottom, Left",   11);
   c->AddEntry("21 Bottom, Middle", 21);
   c->AddEntry("31 Bottom, Right",  31);
   c->AddEntry("12 Middle, Left",   12);
   c->AddEntry("22 Middle, Middle", 22);
   c->AddEntry("32 Middle, Right",  32);
   c->AddEntry("13 Top, Left",   13);
   c->AddEntry("23 Top, Middle", 23);
   c->AddEntry("33 Top, Right",  33);

   return c;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the alpha value from the entry field.

void TAttTextEditor::DoAlphaField()
{
   if (fAvoidSignal) return;

   if (TColor *color = gROOT->GetColor(fAttText->GetTextColor())) {
      color->SetAlpha((Float_t)fAlphaField->GetNumber());
      fAlpha->SetPosition((Int_t)(fAlphaField->GetNumber()*1000));
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set the alpha value

void TAttTextEditor::DoAlpha()
{
   if (fAvoidSignal) return;

   if (TColor *color = gROOT->GetColor(fAttText->GetTextColor())) {
      color->SetAlpha((Float_t)fAlpha->GetPosition()/1000);
      fAlphaField->SetNumber((Float_t)fAlpha->GetPosition()/1000);
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to set alpha value online.

void TAttTextEditor::DoLiveAlpha(Int_t a)
{
   if (fAvoidSignal) return;
   fAlphaField->SetNumber((Float_t)a/1000);

   if (TColor *color = gROOT->GetColor(fAttText->GetTextColor())) {
      // In case the color is not transparent a new color is created.
      if (color->GetAlpha() == 1.) {
         fAttText->SetTextColor(TColor::GetColorTransparent(color->GetNumber(),0.99));
      } else {
         color->SetAlpha((Float_t)a/1000);
      }
   }
   Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot to update alpha value on click on Slider

void TAttTextEditor::GetCurAlpha()
{
   if (fAvoidSignal) return;

   if (TColor *color = gROOT->GetColor(fAttText->GetTextColor())) {
      fAlpha->SetPosition((Int_t)(color->GetAlpha()*1000));
      fAlphaField->SetNumber(color->GetAlpha());
   }
   Update();
}
