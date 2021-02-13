// @(#)root/gui:$Id: bcc7b53cab203466d39a6f17445d822d6fc19610 $
// Author: Bertrand Bellenot + Fons Rademakers + Valeriy Onuchin  23/04/03

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
/**************************************************************************

    This source is based on Xclass95, a Win95-looking GUI toolkit.
    Copyright (C) 1996, 1997 David Barth, Ricky Ralston, Hector Peraza.

    Xclass95 is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

**************************************************************************/

#include "TGFontDialog.h"
#include "TString.h"
#include "TError.h"
#include "TGColorSelect.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGListBox.h"
#include "TGComboBox.h"
#include "TVirtualX.h"
#include <cstdlib>


enum EFontDialog {
   kFDLG_OK          = 1,
   kFDLG_CANCEL,

   kFDLG_FONTNAMES   = 21,
   kFDLG_FONTSIZES,
   kFDLG_FONTSTYLES,
   kFDLG_FONTALIGNS,

   kFDLG_COLORSEL    = 31
};

static const char *gFDLG_DEFAULTSAMPLE = "AaBbCcYyZz 1234567890";


static const char *gFontSizes[] = {
    "8",  "9", "10", "11", "12", "13", "14", "16",
   "18", "20", "22", "24", "26", "28", "30", "32",
   "34", "36", "48", "72", 0
};

static const char *gFontStyles[] = {
   "Normal", "Bold", "Italic", "Bold Italic", 0
};

static TString gFontStylesReal[4];


static const char *gAlignTypes[] = {
    "top left", "top center", "top right",
    "middle left", "middle center", "middle right",
    "bottom left", "bottom center", "bottom right",
    0
};

static const Int_t gAlignValues[] = {
    kTextTop     | kTextLeft,
    kTextTop     | kTextCenterX,
    kTextTop     | kTextRight,
    kTextCenterY | kTextLeft,
    kTextCenterY | kTextCenterX,
    kTextCenterY | kTextRight,
    kTextBottom  | kTextLeft,
    kTextBottom  | kTextCenterX,
    kTextBottom  | kTextRight,
    0
};

/*
static const char *gFontList[] = {
   "Arial",
   "Comic Sans MS",
   "Courier New",
   "Georgia",
   "Impact",
   "Monotype",
   "Symbol",
   "Times New Roman",
   "Trebuchet MS",
   "Verdana",
   "Webdings",
   "Wingdings",
   0,
};

static const char *gFontList2[] = {
   "FreeSans",
   "FreeSerif",
   "FreeMono",
   "OpenSymbol",
   0
};
*/

//_____________________________________________________________________________
//
// TGFontDialog
//
// Font selection dialog, allowing to select one in the list of available fonts
// in the system.
//_____________________________________________________________________________

ClassImp(TGFontDialog);

////////////////////////////////////////////////////////////////////////////////
/// Create font dialog. When closed via OK button fontProp is set to
/// the newly selected font. If closed via Cancel button or WM close box
/// fontProp->fName == "".

TGFontDialog::TGFontDialog(const TGWindow *p, const TGWindow *t,
                           FontProp_t *fontProp, const TString &sample,
                           char **fontList, Bool_t wait) :
   TGTransientFrame(p, t, 100, 100), fFontNames(0), fFontSizes(0), fFontStyles(0),
   fTextAligns(0), fColorSelect(0), fFontProp(0), fItalic(0), fBold(0), fSize(0),
   fTextAlign(0), fTextColor(0), fNumberOfFonts(0)
{
   TGLabel *lbl;
   TGHorizontalFrame *hf, *hf2;
   TGVerticalFrame *vf;
   Int_t i, w;
   UInt_t width = 0, height = 0;

   fSampleTextGC = 0;
   fLabelFont    = 0;
   fSample       = 0;
   fHitOK        = kFALSE;
   fWaitFor      = wait;
   fInitFont     = 0;
   fInitColor    = 0;
   fInitAlign    = 0;

   if (!fontProp) {
      Error("TGFontDialog", "fontProp argument may not be 0");
      return;
   }
   SetCleanup(kDeepCleanup);

   TGLayoutHints *lh1 = new TGLayoutHints(kLHintsLeft | kLHintsExpandY);
   TGLayoutHints *lh2 = new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 0, 2, 2);

   fFontProp = fontProp;
   if (fontProp->fName == "") {
      fName      = fontList ? fontList[0] : "Arial";
      fItalic    = kFALSE;
      fBold      = kFALSE;
      fSize      = 12;
      fTextAlign = gAlignValues[6];   //kTextBottom  | kTextLeft
      fTextColor = TGFrame::GetBlackPixel();
   } else {
      fName      = fontProp->fName;
      fItalic    = fontProp->fItalic;
      fBold      = fontProp->fBold;
      fSize      = fontProp->fSize;
      fTextAlign = fontProp->fAlign;
      fTextColor = fontProp->fColor;
   }

   hf = new TGHorizontalFrame(this, 10, 10);
   AddFrame(hf, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5));

   //--------------------- font names listbox

   vf = new TGVerticalFrame(hf, 10, 10);

   lbl = new TGLabel(vf, new TGString("Font:"));
   vf->AddFrame(lbl, lh2);

   fFontNames = new TGListBox(vf, kFDLG_FONTNAMES);
   fFontNames->Resize(120, fFontNames->GetDefaultHeight());

   if (gVirtualX->InheritsFrom("TGX11") || gVirtualX->InheritsFrom("TGCocoa")) {
      fFontNames->Connect("Selected(char*)", "TGFontDialog", this, "UpdateStyleSize(char*)");
   }

   fFontNames->Associate(this);
   vf->AddFrame(fFontNames,  new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   hf->AddFrame(vf, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY,
                                      0, 10, 2, 2));

   //--------------------- font styles listbox

   vf = new TGVerticalFrame(hf, 10, 10);

   lbl = new TGLabel(vf, new TGString("Style:"));
   vf->AddFrame(lbl, lh2);

   fFontStyles = new TGListBox(vf, kFDLG_FONTSTYLES);
   fFontStyles->Resize(80, fFontStyles->GetDefaultHeight());
   fFontStyles->Associate(this);
   vf->AddFrame(fFontStyles, lh1);

   hf->AddFrame(vf, new TGLayoutHints(kLHintsLeft | kLHintsExpandY,
                                      0, 10, 2, 2));

   //--------------------- font sizes listbox

   vf = new TGVerticalFrame(hf, 10, 10);

   lbl = new TGLabel(vf, new TGString("Size:"));
   vf->AddFrame(lbl, lh2);

   fFontSizes = new TGListBox(vf, kFDLG_FONTSIZES);
   fFontSizes->Resize(50, fFontSizes->GetDefaultHeight());
   fFontSizes->Associate(this);
   vf->AddFrame(fFontSizes, lh1);

   hf->AddFrame(vf, new TGLayoutHints(kLHintsLeft | kLHintsExpandY,
                                      0, 0, 2, 2));

   //--------------------- Text align and color

   vf = new TGVerticalFrame(hf, 10, 10);

   lbl = new TGLabel(vf, new TGString("Text Align :"));
   vf->AddFrame(lbl, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));

   fTextAligns = new TGComboBox(vf, kFDLG_FONTALIGNS);
   fTextAligns->Associate(this);
   vf->AddFrame(fTextAligns, new TGLayoutHints(kLHintsCenterX | kLHintsTop |
                                               kLHintsExpandX, 5, 5, 0, 5));

   fTextAligns->Resize(110, 20);

   vf->Layout();
   vf->Resize(vf->GetDefaultSize());

   w = hf->GetDefaultWidth();

   hf2 = new TGHorizontalFrame(vf, 10, 10);
   vf->AddFrame(hf2, new TGLayoutHints(kLHintsNormal, 0, 0, 5, 5));

   lbl = new TGLabel(hf2, new TGString("Text Color"));
   hf2->AddFrame(lbl, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 5, 5, 20, 5));

   fColorSelect = new TGColorSelect(hf2, fTextColor, kFDLG_COLORSEL);
   fColorSelect->Associate(this);

   hf2->AddFrame(fColorSelect, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 5, 5, 20, 5));

   vf->Resize(vf->GetDefaultSize());

   hf->AddFrame(vf, new TGLayoutHints(kLHintsLeft | kLHintsExpandY | kLHintsTop,
                                      10, 0, lbl->GetDefaultHeight()+6, 0));

   //--------------------- initialize controls

   Resize(GetDefaultSize());

   Int_t cnt = 0;
   Bool_t own = kFALSE;

   if (!fontList) {
      fontList = gVirtualX->ListFonts("-*-*-*-*", 10000, cnt);
      own = kTRUE;
   }

   Build(fontList, cnt);

   for (i = 0; gAlignTypes[i] != 0; ++i) {
      fTextAligns->AddEntry(new TGString(gAlignTypes[i]), i);
   }

   for (i = 0; gAlignValues[i] != 0; ++i) {
      if (gAlignValues[i] == fTextAlign) {
         fTextAligns->Select(i);
      }
   }

   TGLBEntry *entry = fTextAligns->FindEntry("bottom center");
   if (entry)
      fTextAligns->Resize(entry->GetWidth() + 20, entry->GetHeight() > 20 ?
                          entry->GetHeight() : 20);

   vf->Layout();
   vf->Resize(vf->GetDefaultSize());

   //--------------------- sample box

   hf = new TGHorizontalFrame(this, 10, 10);
   AddFrame(hf, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 5, 5, 5, 5));

   TGGroupFrame *gf = new TGGroupFrame(hf, new TGString("Sample"));
   hf->AddFrame(gf, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 3, 3, 3, 3));

   TGCompositeFrame *cf = new TGCompositeFrame(gf, 200, 100, kSunkenFrame);
   gf->AddFrame(cf, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY |
                                      kLHintsExpandX | kLHintsExpandY,
                                      1, 1, 3, 1));

   if (sample == "")
      fSampleText = gFDLG_DEFAULTSAMPLE;
   else
      fSampleText = sample;

   for (i = 0; gFontSizes[i] != 0; ++i) {
      if (fSize == atoi(gFontSizes[i])) {
         fFontSizes->Select(i);
         break;
      }
   }

   if ((fBold == kFALSE) && (fItalic == kFALSE))
      fFontStyles->Select(0);
   if ((fBold == kTRUE) && (fItalic == kFALSE))
      fFontStyles->Select(1);
   if ((fBold == kFALSE) && (fItalic == kTRUE))
      fFontStyles->Select(2);
   if ((fBold == kTRUE) && (fItalic == kTRUE))
      fFontStyles->Select(3);

   GetFontName();   //sets fLName and other data members

   GCValues_t gcval;
   gcval.fMask = kGCForeground | kGCFont;
   gcval.fForeground = fTextColor;
   gcval.fFont = fLabelFont->GetFontHandle();
   fSampleTextGC = fClient->GetGC(&gcval, kTRUE);
   fSample = new TGLabel(cf, fSampleText); //, (*fSampleTextGC)(), (*fLabelFont)());
   fSample->SetTextJustify(gAlignValues[fTextAligns->GetSelected()]);
   cf->AddFrame(fSample, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY |
                                           kLHintsExpandX | kLHintsExpandY,
                                           1, 1, 1, 1));
   cf->Layout();
   if (own) gVirtualX->FreeFontNames(fontList);

   gf->Resize(w, 80);

   //--------------------- OK and Cancel buttons

   vf = new TGVerticalFrame(hf, 10, 10, kFixedWidth);

   TGTextButton *ok = new TGTextButton(vf, new TGHotString("&OK"), kFDLG_OK);
   ok->Associate(this);
   vf->AddFrame(ok, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 5, 5, 20, 5));

   TGTextButton *cancel = new TGTextButton(vf, new TGHotString("&Cancel"), kFDLG_CANCEL);
   cancel->Associate(this);
   vf->AddFrame(cancel, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 5, 5, 5, 5));

   vf->Resize(cancel->GetDefaultWidth()+70, vf->GetDefaultHeight());

   hf->AddFrame(vf, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));
   SetEditDisabled(kEditDisable);

   fInitAlign = fTextAlign;
   fInitColor = fTextColor;
   fInitFont = fLabelFont;

   fFontNames->GetScrollBar()->SetPosition(fFontNames->GetSelected()-3);
   fFontSizes->GetScrollBar()->SetPosition(fFontSizes->GetSelected()-3);
   fFontSizes->Layout();
   fFontNames->Layout();

   SetWindowName("Font Selection");
   SetIconName("Font Selection");
   SetClassHints("ROOT", "FontDialog");

   SetMWMHints(kMWMDecorAll | kMWMDecorResizeH  | kMWMDecorMaximize |
                              kMWMDecorMinimize | kMWMDecorMenu,
               kMWMFuncAll  | kMWMFuncResize    | kMWMFuncMaximize |
                              kMWMFuncMinimize,
               kMWMInputModeless);

   //---- set minimum size to the dialog box
   width  = GetDefaultWidth();
   height = GetDefaultHeight();
   MapSubwindows();
   Resize(width, height);

   //---- position the dialog relative to the parent's window
   CenterOnParent();

   // set minimum size
   SetWMSize(width, height);
   SetWMSizeHints(width, height, 10000, 10000, 0, 0);

   if (fWaitFor) {
      MapWindow();
      fClient->WaitForUnmap(this);
      DeleteWindow();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Delete all widgets.

TGFontDialog::~TGFontDialog()
{
   //fClient->FreeFont(fLabelFont);
   fClient->FreeGC(fSampleTextGC);
}

////////////////////////////////////////////////////////////////////////////////
/// Called when window is closed via window manager.

void TGFontDialog::CloseWindow()
{
   if (fWaitFor) {
      UnmapWindow();
      return;
   }

   if (!fHitOK) {
      fFontProp->fName = "";

      if (fInitFont) {
         SetFont(fInitFont);
         FontSelected((char*)fInitFont->GetName());
      }
      if (fInitColor) {
         SetColor(fInitColor);
         ColorSelected(fInitColor);
      }

      if (fInitAlign) {
         SetAlign(fInitAlign);
         AlignSelected(fInitAlign);
      }
   }
   fFontNames->Select(0);
   fFontStyles->Select(0);
   fFontSizes->Select(0);

   // don't call DeleteWindow() here since that will cause access
   // to the deleted dialog in the WaitFor() method (see ctor)
   UnmapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Handle dialog events.

Bool_t TGFontDialog::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
   switch (GET_MSG(msg)) {
      case kC_COMMAND:
         switch (GET_SUBMSG(msg)) {
            case kCM_BUTTON:
               switch (parm1) {

                  case kFDLG_OK:
                     fHitOK             = kTRUE;
                     fFontProp->fName   = fName;
                     fFontProp->fItalic = fItalic;
                     fFontProp->fBold   = fBold;
                     fFontProp->fSize   = fSize;
                     fFontProp->fColor  = fTextColor;
                     fFontProp->fAlign  = fTextAlign;
                     CloseWindow();
                     break;
                  case kFDLG_CANCEL:
                     fHitOK = kFALSE;
                     CloseWindow();
                     break;
               }
               break;

            case kCM_LISTBOX:
               GetFontName();
               break;

            case kCM_COMBOBOX:
               GetFontName();
               break;

            default:
               break;
         }
         break;

      case kC_COLORSEL:
         switch (GET_SUBMSG(msg)) {

            case kCOL_SELCHANGED:
               {
                  if (parm2 != (Int_t)fTextColor) {
                     fTextColor = parm2;
                     GetFontName();
                  }
               }
               break;

            default:
               break;
         }
         break;

      default:
         break;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Build font dialog.

Bool_t TGFontDialog::Build(char **fontList, Int_t cnt)
{
   TString family;
   TString font;

   Int_t  i;
   fNumberOfFonts = 1;
   Int_t n1, n2;

   for (i = 0; i < cnt; i++) {
      font = fontList[i];
      n1 = font.Index("-", 1);
      n2 = font.Index("-", n1+1);
      family = font(n1+1, n2-n1-1);

      if ((family.Length() == 0) || (family.BeginsWith("@"))) {
         continue;
      }

      if (!fFontNames->FindEntry(family.Data())) {
         fFontNames->AddEntry(family.Data(), fNumberOfFonts++);
      }
   }

   fFontNames->SortByName();

   TGTextLBEntry *le = (TGTextLBEntry*)fFontNames->FindEntry(fName.Data());

   if (le) {
      fFontNames->Select(le->EntryId());
   }

   UpdateStyleSize(fName.Data());

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Build font style and size list boxes

void TGFontDialog::UpdateStyleSize(const char *family)
{
   if (!family) {
      return;
   }

   TString font = family;
   Bool_t styles[4] = { kFALSE, kFALSE, kFALSE, kFALSE };
   Int_t cnt = 0;
   Int_t i = 0;

   TString fname;
   char **fontList = 0;

   fname = TString::Format("-*-%s-*-*", family);
   fontList = gVirtualX->ListFonts(fname.Data(), 1000, cnt);

   fFontSizes->RemoveEntries(0, 1000);
   fFontSizes->Layout();

   fFontStyles->RemoveEntries(0, 100);
   fFontStyles->Layout();

   if (!cnt || !fontList) {
      return;
   }

   TString style1;
   TString style2;
   TString sz;
   TString name;
   Int_t n1, n2;
   Bool_t x11 = gVirtualX->InheritsFrom("TGX11");
   Bool_t all_sizes = !x11;
   Bool_t all_styles = !x11;

   //
   if (gVirtualX->InheritsFrom("TGCocoa")) {
      all_sizes = kTRUE;
      all_styles = kFALSE;
   }

   int szn = 0;

   fFontSizes->AddEntry("12", szn++);

   for (i = 0; i < cnt; i++) {
      name = fontList[i];
      n1 = name.Index(family);

      if (n1 == kNPOS) {
         break;
      }
      n1 += font.Length() + 1;
      n2 = name.Index("-", n1);
      if (n2 == kNPOS) {
         break;
      }

      style1 = name(n1, n2 - n1);

      n1 = n2 + 1;
      n2 = name.Index("-", n1);
      if (n2 == kNPOS) {
         break;
      }
      style2 = name(n1, n2 - n1);

      if ((style1 == "normal") || (style1 == "medium")) {
         if (style2 == "r") {
            styles[0]  = kTRUE;
            gFontStylesReal[0] = style1 + "-" + style2;
         } else if (style2 == "i") {
            styles[2]  = kTRUE;
            gFontStylesReal[2] = style1 + "-" + style2;
         } else if (style2 == "o") {
            styles[2]  = kTRUE;
            gFontStylesReal[2] = style1 + "-" + style2;
         }
      } else if (style1 == "bold") {
         if (style2 == "r") {
            styles[1]  = kTRUE;
            gFontStylesReal[1] = style1 + "-" + style2;
         } else if (style2 == "i") {
            styles[3]  = kTRUE;
            gFontStylesReal[3] = style1 + "-" + style2;
         } else if (style2 == "o") {
            styles[3]  = kTRUE;
            gFontStylesReal[3] = style1 + "-" + style2;
         }
      } else if (style1 == "(null)") {
         styles[0]  = kTRUE;
         gFontStylesReal[0] = "normal-r";
         styles[1]  = kTRUE;
         gFontStylesReal[1] = "bold-r";
         styles[2]  = kTRUE;
         gFontStylesReal[2] = "normal-i";
         styles[3]  = kTRUE;
         gFontStylesReal[3] = "bold-i";
      }

      n1++;
      n2 = name.Index("-", n1);
      n1 = n2 + 1;
      n2 = name.Index("-", n1);
      n1 = n2 + 1;

      if (n2 != kNPOS) {
         n1 = n2 + 2;
         n2 = name.Index("-", n1);
         sz = name(n1, n2 - n1);
         if (!sz.IsDigit()) {
            continue;
         }

         all_sizes = (sz == "0") && !x11;
         if (!all_sizes) {
            sz.Strip();
            if (sz.Length() == 1) {
               sz = " " + sz;
            }
            if (!fFontSizes->FindEntry(sz.Data())) {
               fFontSizes->AddEntry(sz.Data(), szn++);
            }
         }
      }
   }
   gVirtualX->FreeFontNames(fontList);

   Bool_t nostyles = kTRUE;
   for (i = 0; gFontStyles[i] != 0; ++i) {
      if (all_styles || styles[i]) {
         nostyles = kFALSE;
         fFontStyles->AddEntry(new TGString(gFontStyles[i]), i);
      }
   }

   if (nostyles && x11) {
      fFontStyles->AddEntry(new TGString(gFontStyles[0]), 0);
      fBold = kFALSE;
      fItalic = kFALSE;
   }

   if (!fBold) {
      if (!fItalic) {
         fFontStyles->Select(0);
      } else {
         fFontStyles->Select(2);
      }
   } else {
      if (!fItalic) {
         fFontStyles->Select(1);
      } else {
         fFontStyles->Select(3);
      }
   }

//   if (nostyles) {
//      fFontNames->RemoveEntry(fFontNames->FindEntry(family)->EntryId());
//   }

   fFontStyles->MapSubwindows();
   fFontStyles->Layout();
//

   sz = TString::Format("%d", fSize);
   if (sz.Length() == 1) {
      sz = " " + sz;
   }

   for (i = 0; gFontSizes[i] != 0; ++i) {
      if (all_sizes && !fFontSizes->FindEntry(gFontSizes[i])) {
         fFontSizes->AddEntry(new TGString(gFontSizes[i]), i);
      }
      if (sz == gFontSizes[i]) {
         fFontSizes->Select(i);
      }
   }

   fFontSizes->SortByName();
   fFontSizes->MapSubwindows();
   fFontSizes->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Sets fLName and other data members.

void TGFontDialog::GetFontName()
{
   TGTextLBEntry *e;
   const char *size, *name;
   Int_t sel;
   Int_t sav = gErrorIgnoreLevel;
   gErrorIgnoreLevel = kFatal;

   TString oldName = fName;
   e = (TGTextLBEntry *) fFontNames->GetSelectedEntry();

   if (!e) {
      fFontNames->Select(1);
      e = (TGTextLBEntry *) fFontNames->GetSelectedEntry();
   }
   name = e ? e->GetText()->GetString() : "";
   fName = name;

   e = (TGTextLBEntry *) fFontSizes->GetSelectedEntry();
   size = e ? e->GetText()->GetString() : "0";
   fSize = atoi(size);

   sel = fFontStyles->GetSelected();
   if (sel < 0) sel = 0;

   switch(sel) {
      case 0:
         fItalic = kFALSE;
         fBold   = kFALSE;
         break;
      case 1:
         fItalic = kFALSE;
         fBold   = kTRUE;
         break;
      case 2:
         fItalic = kTRUE;
         fBold   = kFALSE;
         break;
      case 3:
         fItalic = kTRUE;
         fBold   = kTRUE;
         break;
      default:
         fItalic = kFALSE;
         fBold   = kFALSE;
         break;
   }

   const char *rgstry = "*";

   if ((fName == "Symbol") || (fName == "Webdings") || (fName == "Wingdings")) {
      rgstry = "microsoft";
   }

   TString oldFont = fLName;
   fLName = TString::Format("-*-%s-%s-*-*-%s-*-*-*-*-*-%s-*", name,
                            gFontStylesReal[sel].Data(), size, rgstry);

   if (oldFont != fLName) {
      if (fLabelFont) {
         //fClient->FreeFont(fLabelFont);
         fLabelFont = fClient->GetFont(fLName, kFALSE);

         if (!fLabelFont) {
            fLabelFont = fClient->GetFont("fixed");
         }
      } else {
         fLabelFont = fClient->GetFont("fixed");
      }

      if (!fLabelFont) {
         // should not happen
         fLName = oldFont;
      }
      else if (fSample) {
         fSample->SetTextFont(fLabelFont);
      }
   }

   Int_t oldAlign = fTextAlign;

   Int_t idx = fTextAligns->GetSelected();
   fTextAlign = gAlignValues[idx >= 0 ? idx : 6];

   if (fSample) {
      if (fTextAlign != oldAlign) {
         fSample->SetTextJustify(fTextAlign);
         AlignSelected(fTextAlign);
      }
      fSample->SetTextColor(fTextColor);
      fColorSelect->SetColor(fTextColor, kFALSE);
      ColorSelected(fTextColor);
   }
   FontSelected((char*)fLName.Data());
   fClient->NeedRedraw(this);
   gErrorIgnoreLevel = sav;
}

////////////////////////////////////////////////////////////////////////////////
/// Set font

void TGFontDialog::SetFont(TGFont *font)
{
   if (!font) {
      return;
   }
   TString name = font->GetName();

   if (name.Index("-", 1) == kNPOS) {
      return;
   }

   if (fSample) {
      fLabelFont = font;
      fSample->SetTextFont(fLabelFont);
   }
   fInitFont = font;

   TString style1;
   TString style2;
   TString sz;

   TString family;
   Int_t n1, n2;

   n1 = name.Index("-", 1);
   n2 = name.Index("-", n1 + 1);
   n1++;
   family = name(n1, n2 - n1);

   TGTextLBEntry *le = (TGTextLBEntry*)fFontNames->FindEntry(family.Data());

   if (le) {
      fFontNames->Select(le->EntryId());
   }

   n1 = n2 + 1;
   n2 = name.Index("-", n1);

   style1 = name(n1, n2 - n1);

   n1 = n2 + 1;
   n2 = name.Index("-", n1);
   if (n2 == kNPOS) {
      return;
   }
   style2 = name(n1, n2 - n1);

   if ((style1 == "normal") || (style1 == "medium")) {
      if (style2 == "r") {
         fFontStyles->Select(0);
      } else if (style2 == "i") {
         fFontStyles->Select(2);
      } else if (style2 == "o") {
         fFontStyles->Select(2);
      }
   } else if (style1 == "bold") {
      if (style2 == "r") {
         fFontStyles->Select(1);
      } else if (style2 == "i") {
         fFontStyles->Select(3);
      } else if (style2 == "o") {
         fFontStyles->Select(3);
      }
   }
   n1++;
   n2 = name.Index("-", n1);
   n1 = n2 + 1;
   n2 = name.Index("-", n1);
   n1 = n2 + 1;
   n2 = name.Index("-", n1);
   n1 = n2 + 1;
   if (n2 != kNPOS) {
      n1 = n2 + 1;
      n2 = name.Index("-", n1);
      sz = name(n1, n2 - n1);

      le = (TGTextLBEntry*)fFontSizes->FindEntry(sz.Data());
      if (le) {
         fFontSizes->Select(le->EntryId());
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set color

void TGFontDialog::SetColor(Pixel_t color)
{
   if (fSample) {
      fTextColor = color;
      fSample->SetTextColor(fTextColor);
   }
   fColorSelect->SetColor(color, kFALSE);
   fClient->NeedRedraw(fColorSelect);
   fInitColor = color;
}

////////////////////////////////////////////////////////////////////////////////
/// Set align

void TGFontDialog::SetAlign(Int_t align)
{
   if (fSample) {
      fTextAlign = align;
      fSample->SetTextJustify(fTextAlign);
   }

   for (int i = 0; gAlignValues[i] != 0; ++i) {
      if (gAlignValues[i] == align) {
         fTextAligns->Select(i);
      }
   }
   fInitAlign = align;
   fClient->NeedRedraw(fTextAligns);
}

////////////////////////////////////////////////////////////////////////////////
/// Enable/disable align combobox

void TGFontDialog::EnableAlign(Bool_t on)
{
   fTextAligns->SetEnabled(on);
}
