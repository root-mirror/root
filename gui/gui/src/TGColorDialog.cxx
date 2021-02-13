// @(#)root/gui:$Id: c1ecfb3a4b91442ae3b6fe3059ae838e463f2f56 $
// Author: Bertrand Bellenot + Fons Rademakers   22/08/02
// Author: Ilka Antcheva (color wheel support)   16/03/07

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
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

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGColorPalette, TGColorPick and TGColorDialog.                       //
//                                                                      //
// The TGColorPalette is a widget showing an matrix of color cells. The //
// colors can be set and selected.                                      //
//                                                                      //
// The TGColorPick is a widget which allows a color to be picked from   //
// HLS space. It consists of two elements: a color map window from      //
// where the user can select the hue and saturation level of a color,   //
// and a slider to select color's lightness.                            //
//                                                                      //
// Selecting a color in these two widgets will generate the event:      //
// kC_COLORSEL, kCOL_CLICK, widget id, 0.                               //
// and the signal:                                                      //
// ColorSelected(Pixel_t color)                                         //
//                                                                      //
// The TGColorDialog presents a full featured color selection dialog.   //
// It uses 2 TGColorPalette's and the TGColorPick widgets.              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <cstdlib>

#include "TGLabel.h"
#include "TGMsgBox.h"         // for ID_OK, ID_CANCEL
#include "TGLayout.h"
#include "TGGC.h"
#include "KeySymbols.h"
#include "TGColorDialog.h"
#include "TGTextEntry.h"
#include "TGButton.h"
#include "TGResourcePool.h"
#include "TColor.h"
#include "TColorWheel.h"
#include "TGColorSelect.h"
#include "TGTab.h"
#include "TRootEmbeddedCanvas.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"
#include "TVirtualX.h"
#include "snprintf.h"

ClassImp(TGColorPalette);
ClassImp(TGColorPick);
ClassImp(TGColorDialog);


// TODO:
// - implement "custom" colors.
// - optimize the code, specially the one handling the fColormap image
//   and dithering in pseudo-color modes; remove duplicated code.
// - improve the color allocation routine.
// - use a buffering pixmap for the fColormap image.

enum EColorDialog {
   kCDLG_OK       = 100,
   kCDLG_CANCEL,
   kCDLG_PREVIEW,
   kCDLG_ADD,

   kCDLG_SPALETTE = 200,
   kCDLG_CPALETTE,
   kCDLG_COLORPICK,

   kCDLG_HTE      = 300,
   kCDLG_LTE,
   kCDLG_STE,
   kCDLG_RTE,
   kCDLG_GTE,
   kCDLG_BTE,
   kCDLG_ALE
};

enum EColorPick {
   kCLICK_NONE,
   kCLICK_HS,
   kCLICK_L
};

enum EColorImage {
   kIMG_HS,
   kIMG_L
};

// "User" defined colors

static ULong_t gUcolor[24] = { 0xff000000 };


////////////////////////////////////////////////////////////////////////////////
/// TGColorPalette widget: this is just a grid of color cells of the
/// specified size. Colors can be selected by clicking on them or by
/// using the arrow keys.

TGColorPalette::TGColorPalette(const TGWindow *p, Int_t cols, Int_t rows, Int_t id) :
   TGFrame(p, 10, 10, kChildFrame)
{
   fWidgetId    = id;
   fWidgetFlags = kWidgetIsEnabled;
   fMsgWindow   = p;
   fDrawGC      = *fClient->GetResourcePool()->GetFrameGC();

   fCw = 20;
   fCh = 17;

   fRows = rows;
   fCols = cols;

   fCx = fCy = 0;

   fPixels = new ULong_t[fRows * fCols];

   for (Int_t i = 0; i < fRows * fCols; ++i) {
      fPixels[i] = TColor::RGB2Pixel(255, 255, 255);
   }

   gVirtualX->GrabButton(fId, kAnyButton, kAnyModifier,
                         kButtonPressMask | kButtonReleaseMask |
                         kPointerMotionMask, kNone, kNone);

   AddInput(kKeyPressMask | kEnterWindowMask | kLeaveWindowMask |
            kFocusChangeMask | kStructureNotifyMask);
   fEditDisabled = kEditDisable;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGColorPalette::~TGColorPalette()
{
   delete [] fPixels;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle button events in color palette

Bool_t TGColorPalette::HandleButton(Event_t *event)
{
   if (event->fCode != kButton1)
      return kFALSE;

   if ((event->fType == kButtonPress) && HasFocus())
      WantFocus();

   Int_t cx = event->fX / (fCw + 5);
   Int_t cy = event->fY / (fCh + 5);

   if (cx >= 0 && cx < fCols && cy >= 0 && cy < fRows) {

      DrawFocusHilite(kFALSE);

      fCx = cx;
      fCy = cy;

      DrawFocusHilite(kTRUE);

      SendMessage(fMsgWindow, MK_MSG(kC_COLORSEL, kCOL_CLICK), fWidgetId, 0);
      ColorSelected();
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse motion events in color palette.

Bool_t TGColorPalette::HandleMotion(Event_t *event)
{
   if (!IsEnabled())
      return kTRUE;

   Int_t cx = event->fX / (fCw + 5);
   Int_t cy = event->fY / (fCh + 5);

   if (cx >= 0 && cx < fCols && cy >= 0 && cy < fRows) {

      DrawFocusHilite(kFALSE);

      fCx = cx;
      fCy = cy;

      DrawFocusHilite(kTRUE);

      SendMessage(fMsgWindow, MK_MSG(kC_COLORSEL, kCOL_CLICK), fWidgetId, 0);
      ColorSelected();
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle keyboard events in color palette.

Bool_t TGColorPalette::HandleKey(Event_t *event)
{
   Char_t input[10];
   UInt_t keysym;

   if (event->fType == kGKeyPress) {

      gVirtualX->LookupString(event, input, sizeof(input), keysym);

      Int_t cx = fCx;
      Int_t cy = fCy;

      switch ((EKeySym)keysym) {
         case kKey_Left:
            if (cx > 0) --cx;
            break;

         case kKey_Right:
            if (cx < fCols - 1) ++cx;
            break;

         case kKey_Up:
            if (cy > 0) --cy;
            break;

         case kKey_Down:
            if (cy < fRows - 1) ++cy;
            break;

         case kKey_Home:
            cx = cy = 0;
            break;

         case kKey_End:
            cx = fCols - 1;
            cy = fRows - 1;
            break;

         default:
            break;
      }

      if (cx != fCx || cy != fCy) {

         DrawFocusHilite(kFALSE);

         fCx = cx;
         fCy = cy;

         DrawFocusHilite(kTRUE);

         SendMessage(fMsgWindow, MK_MSG(kC_COLORSEL, kCOL_CLICK), fWidgetId, 0);
         ColorSelected();
      }
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set color entries in color samples.

void TGColorPalette::SetColors(ULong_t colors[])
{
   for (Int_t i = 0; i < fRows * fCols; ++i)
      SetColor(i, colors[i]);
   gClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Set color at index ix of color entries.

void TGColorPalette::SetColor(Int_t ix, ULong_t color)
{
   fPixels[ix] = color;
   gClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Set current cell color.

void TGColorPalette::SetCurrentCellColor(ULong_t color)
{
   SetColor(fCy * fCols + fCx, color);
}

////////////////////////////////////////////////////////////////////////////////
/// Set color cell size.

void TGColorPalette::SetCellSize(Int_t w, Int_t h)
{
   fCw = w;
   fCh = h;
   gClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Return currently selected color value.

ULong_t TGColorPalette::GetCurrentColor() const
{
   if (fCx >= 0 && fCy >= 0)
      return GetColorByIndex(fCy * fCols + fCx);
   else
      return TColor::RGB2Pixel(0, 0, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw color palette.

void TGColorPalette::DoRedraw()
{
   Int_t i, j, k, x, y;

   k = 0;
   y = 2;
   for (i = 0; i < fRows; ++i) {
      x = 2;
      for (j = 0; j < fCols; ++j) {
         Draw3dRectangle(kSunkenFrame | kDoubleBorder, x, y, fCw, fCh);
         fDrawGC.SetForeground(fPixels[k++]);
         gVirtualX->FillRectangle(fId, fDrawGC(), x + 2, y + 2, fCw - 4, fCh - 4);
         x += fCw + 5;
      }
      y += fCh + 5;
   }

   DrawFocusHilite(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Add keyboard input.

void TGColorPalette::GotFocus()
{
   AddInput(kKeyPressMask | kKeyReleaseMask);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove keyboard input.

void TGColorPalette::LostFocus()
{
   RemoveInput(kKeyPressMask | kKeyReleaseMask);
   gClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a highlight rectangle around cell obtaining focus.

void TGColorPalette::DrawFocusHilite(Int_t onoff)
{
   if (fCx >= 0 && fCy >= 0) {
      GContext_t gc = onoff ? GetShadowGC()() : GetBckgndGC()();
      gVirtualX->DrawRectangle(fId, gc, fCx * (fCw + 5) + 0, fCy * (fCh + 5) + 0,
                               fCw + 3, fCh + 3);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// TGColorPick constructor.
/// TGColorPick is a widget which allows a color to be picked from HLS space.
/// It consists of two elements: a color map window from where the user can
/// select the hue and saturation level of a color, and a slider to select
/// color's lightness.

TGColorPick::TGColorPick(const TGWindow *p, Int_t w, Int_t h, Int_t id) :
   TGFrame(p, w, h, kChildFrame), fCursorGC(GetBlackGC())
{
   UInt_t iw, ih;

   fWidgetId    = id;
   fWidgetFlags = kWidgetIsEnabled;
   fMsgWindow   = p;

   fColormapRect.fX = 1;
   fColormapRect.fY = 1;
   fColormapRect.fWidth = w - 33 - 2;
   fColormapRect.fHeight = h - 2;
   fSliderRect.fX = w - 18 - 2;
   fSliderRect.fY = 1;
   fSliderRect.fWidth = 10;
   fSliderRect.fHeight = h - 2;

   fNColors = 0;

   if (!p) {
      MakeZombie();
      // coverity[uninit_member]
      return;
   }
   CreateImages();
   gVirtualX->GetImageSize(fLimage, iw, ih);

   fCx = 0;
   fCy = 0;
   fCz = (Int_t)ih / 2;

   fClick = kCLICK_NONE;

   UpdateCurrentColor();
   InitImages();

   gVirtualX->GrabButton(fId, kAnyButton, kAnyModifier,
                         kButtonPressMask | kButtonReleaseMask |
                         kPointerMotionMask, kNone, kNone);

   AddInput(kKeyPressMask | kEnterWindowMask | kLeaveWindowMask |
            kFocusChangeMask | kStructureNotifyMask);
   fEditDisabled = kEditDisable;
}

////////////////////////////////////////////////////////////////////////////////
/// TGColorPick destructor.

TGColorPick::~TGColorPick()
{
   if (IsZombie()) return;
   gVirtualX->DeleteImage(fHSimage);
   gVirtualX->DeleteImage(fLimage);
   FreeColors();
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button events in color pick widget.

Bool_t TGColorPick::HandleButton(Event_t *event)
{
   if (event->fCode != kButton1) return kFALSE;

   if (event->fType == kButtonPress) {
      if ((event->fX > fColormapRect.fX) && (event->fX < fColormapRect.fX + fColormapRect.fWidth) &&
          (event->fY > fColormapRect.fY) && (event->fY < fColormapRect.fY + fColormapRect.fHeight)) {

         fClick = kCLICK_HS;
         SetHScursor(event->fX - fColormapRect.fX, event->fY - fColormapRect.fY);

      } else if (event->fX > fSliderRect.fX) {

         fClick = kCLICK_L;
         SetLcursor(event->fY - fSliderRect.fY);

      }
   } else {  // ButtonRelease

      fClick = kCLICK_NONE;

   }

   UpdateCurrentColor();
   if (fClick == kCLICK_HS) SetSliderColor();

   SendMessage(fMsgWindow, MK_MSG(kC_COLORSEL, kCOL_CLICK), fWidgetId, kFALSE);
   ColorSelected();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse motion events in color pick widget.

Bool_t TGColorPick::HandleMotion(Event_t *event)
{
   if (!IsEnabled())
      return kTRUE;

   if (fClick == kCLICK_HS) {

      SetHScursor(event->fX - fColormapRect.fX, event->fY - fColormapRect.fY);

   } else if (fClick == kCLICK_L) {

      SetLcursor(event->fY - fSliderRect.fY);

   } else {

      return kTRUE;

   }

   UpdateCurrentColor();
   if (fClick == kCLICK_HS) SetSliderColor();

   SendMessage(fMsgWindow, MK_MSG(kC_COLORSEL, kCOL_CLICK), fWidgetId, kFALSE);
   ColorSelected();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create colormap and color slider images.

void TGColorPick::CreateImages()
{
   UInt_t width, height;

   width = fColormapRect.fWidth;
   height = fColormapRect.fHeight;
   fHSimage = gVirtualX->CreateImage(width, height);
   width = fSliderRect.fWidth;
   height = fSliderRect.fHeight;
   fLimage = gVirtualX->CreateImage(width, height);
}

////////////////////////////////////////////////////////////////////////////////
/// Try to allocate first a palette of 64 colors. Used by the dithered
/// version of the color maps.

void TGColorPick::AllocColors()
{
   ColorStruct_t color;
   Int_t i;

   for (i = 0; i < 64; ++i) {
      Int_t cc[4] = { 0, 21845, 43691, 65535 };
      color.fPixel = 0;
      color.fRed   = cc[i & 0x3];
      color.fGreen = cc[(i >> 2) & 0x3];
      color.fBlue  = cc[(i >> 4) & 0x3];
      if (gVirtualX->AllocColor(gVirtualX->GetColormap(), color) == 0)
         break;
      fColormap[i][0] = color.fRed / 256;
      fColormap[i][1] = color.fGreen / 256;
      fColormap[i][2] = color.fBlue / 256;
      fPixel[i] = color.fPixel;
   }

   fNColors = i;
   if (fNColors == 64) return;  // success

   // Failed, try a simpler 27-color.

   FreeColors();

   for (i = 0; i < 27; ++i) {
      Int_t cc[3] = { 0, 32768, 65535 };
      color.fPixel = 0;
      color.fRed   = cc[i % 3];
      color.fGreen = cc[(i / 3) % 3];
      color.fBlue  = cc[(i / 9) % 3];
      if (gVirtualX->AllocColor(gVirtualX->GetColormap(), color) == 0)
         break;
      fColormap[i][0] = color.fRed / 256;
      fColormap[i][1] = color.fGreen / 256;
      fColormap[i][2] = color.fBlue / 256;
      fPixel[i] = color.fPixel;
   }

   fNColors = i;
   if (fNColors == 27) return;  // success

   // Failed, try then a much simpler 8-color.

   FreeColors();

   for (i = 0; i < 8; ++i) {
      color.fPixel = 0;
      color.fRed   = (i & 1) * 65535;
      color.fGreen = ((i >> 1) & 1) * 65535;
      color.fBlue  = ((i >> 2) & 1) * 65535;
      if (gVirtualX->AllocColor(gVirtualX->GetColormap(), color) == 0)
         break;
      fColormap[i][0] = color.fRed / 256;
      fColormap[i][1] = color.fGreen / 256;
      fColormap[i][2] = color.fBlue / 256;
      fPixel[i] = color.fPixel;
   }

   fNColors = i;
   if (fNColors == 8) return;  // success

   // Failed, try to get at least 8 closest colors...
   // (TODO: search for closest colors in the colormap, right now we just
   // get as many as exact colors we can for the 8-color palette)

   FreeColors();

   for (i = 0; i < 8; ++i) {
      color.fPixel = 0;
      color.fRed   = (i & 1) * 65535;
      color.fGreen = ((i >> 1) & 1) * 65535;
      color.fBlue  = ((i >> 2) & 1) * 65535;
      if (gVirtualX->AllocColor(gVirtualX->GetColormap(), color) != 0) {
         fColormap[fNColors][0] = color.fRed / 256;
         fColormap[fNColors][1] = color.fGreen / 256;
         fColormap[fNColors][2] = color.fBlue / 256;
         fPixel[fNColors++] = color.fPixel;
      }
   }

   // continue with what we got...
}

////////////////////////////////////////////////////////////////////////////////
/// Free allocated colors.

void TGColorPick::FreeColors()
{
   for (Int_t i = 0; i < fNColors; i++)
      gVirtualX->FreeColor(gVirtualX->GetColormap(), fPixel[i]);
   fNColors = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a dithered version of the color map and lightness images for
/// display modes with reduced number of colors. The Floyd-Steinberg error
/// diffusion dithering algorithm is used.
/// This routine is called in PseudoColor modes only.

void TGColorPick::CreateDitheredImage(Pixmap_t image, Int_t which)
{
   const Int_t kWidth = 20;

   ColorStruct_t line[kWidth];
   struct { Int_t r, g, b; } ed[kWidth], ef;
   Int_t  x, y, c, v, e[4], nc = 0;
   Int_t  r, g, b;
   Int_t h, l, s;
   Long_t dist, sdist;
   Int_t iw, ih;

   gVirtualX->GetImageSize(image, (UInt_t&) iw, (UInt_t&) ih);

   for (x = 0; x < iw; ++x) {
      ed[x].r = ed[x].g = ed[x].b = 0;
   }

   if (fNColors == 0) AllocColors();

   for (y = 0; y < ih; ++y) {

      if (which == kIMG_HS) {

         for (x = 0; x < iw; ++x) {

            h = x * 255 / iw;
            l = 128;
            s = (ih - y) * 255 / ih;

            TColor::HLS2RGB(h, l, s, r, g, b);

            line[x].fRed   = r;
            line[x].fGreen = g;
            line[x].fBlue  = b;
         }

      } else if (which == kIMG_L) {

         TColor::Pixel2RGB(fCurrentColor, r, g, b);
         TColor::RGB2HLS(r, g, b, h, l, s);

         Int_t ll = (ih - y) * 255 / ih;

         TColor::HLS2RGB(h, ll, s, r, g, b);

         for (x = 0; x < iw; ++x) {
            line[x].fRed   = r;
            line[x].fGreen = g;
            line[x].fBlue  = b;
         }

      } else {

         return;

      }

      ef.r = ef.g = ef.b = 0;        // no forward error for first pixel

      for (x = 0; x < iw; ++x) {

         // add errors from previous line

         v = line[x].fRed + ed[x].r;
         if (v < 0) v = 0; else if (v > 255) v = 255;
         line[x].fRed = v;

         v = line[x].fGreen + ed[x].g;
         if (v < 0) v = 0; else if (v > 255) v = 255;
         line[x].fGreen = v;

         v = line[x].fBlue + ed[x].b;
         if (v < 0) v = 0; else if (v > 255) v = 255;
         line[x].fBlue = v;

      }

      for (x = 0; x < iw; ++x) {

         // add forward errors

         v = line[x].fRed + ef.r;
         if (v < 0) v = 0; else if (v > 255) v = 255;
         line[x].fRed = v;

         v = line[x].fGreen + ef.g;
         if (v < 0) v = 0; else if (v > 255) v = 255;
         line[x].fGreen = v;

         v = line[x].fBlue + ef.b;
         if (v < 0) v = 0; else if (v > 255) v = 255;
         line[x].fBlue = v;

         // find the nearest color in colormap[]

         sdist = 255L * 255L * 255L;
         for (c = 0; c < fNColors; ++c) {

            Int_t dr = line[x].fRed   - fColormap[c][0];
            Int_t dg = line[x].fGreen - fColormap[c][1];
            Int_t db = line[x].fBlue  - fColormap[c][2];

            dist = dr * dr + dg * dg + db * db;
            if (dist < sdist) {
               nc = c;
               sdist = dist;
            }
         }

         gVirtualX->PutPixel(image, x, y, fPixel[nc]);

#define FILTER(v) \
      e[0] = (7 * v) >> 4; \
      e[1] = v >> 4;       \
      e[2] = (5 * v) >> 4; \
      e[3] = (3 * v) >> 4;

         v = line[x].fRed - fColormap[nc][0];
         FILTER(v)

         ef.r = e[0];
         if (x < iw-1) ed[x+1].r = e[1];
         if (x == 0) ed[x].r = e[2]; else ed[x].r += e[2];
         if (x > 0) ed[x-1].r += e[3];

         v = line[x].fGreen - fColormap[nc][1];
         FILTER(v)

         ef.g = e[0];
         if (x < iw-1) ed[x+1].g = e[1];
         if (x == 0) ed[x].g = e[2]; else ed[x].g += e[2];
         if (x > 0) ed[x-1].g += e[3];

         v = line[x].fBlue - fColormap[nc][2];
         FILTER(v)

         ef.b = e[0];
         if (x < iw-1) ed[x+1].b = e[1];
         if (x == 0) ed[x].b = e[2]; else ed[x].b += e[2];
         if (x > 0) ed[x-1].b += e[3];

      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize color palette and slider images.

void TGColorPick::InitImages()
{
   Int_t width, height;
   Int_t h, l, s;
   Int_t r, g, b;

   gVirtualX->GetImageSize(fHSimage, (UInt_t&) width, (UInt_t&) height);

   // initialize fHSimage

   Int_t ncolors = gVirtualX->GetDepth();

   if (ncolors > 8) {
      for (Int_t y = 0; y < height; ++y) {
         for (Int_t x = 0; x < width; ++x) {

            r = g = b = 0;
            h = x * 255 / width;
            l = 128;
            s = (height - y) * 255 / height;

            TColor::HLS2RGB(h, l, s, r, g, b);

            ULong_t pixel = TColor::RGB2Pixel(r, g, b);
            gVirtualX->PutPixel(fHSimage, x, y, pixel);
         }
      }
   } else {
      CreateDitheredImage(fHSimage, kIMG_HS);
   }

   // initialize fLimage

   SetSliderColor();
}

////////////////////////////////////////////////////////////////////////////////
/// Set slider colors.

void TGColorPick::SetSliderColor()
{
   Int_t width, height;
   Int_t h, l, s;
   Int_t r, g, b;

   gVirtualX->GetImageSize(fLimage, (UInt_t&) width, (UInt_t&) height);

   Int_t ncolors = gVirtualX->GetDepth();

   if (ncolors > 8) {

      for (Int_t y = 0; y < height; ++y) {

         TColor::Pixel2RGB(fCurrentColor, r, g, b);
         TColor::RGB2HLS(r, g, b, h, l, s);

         l = (height - y) * 255 / height;

         TColor::HLS2RGB(h, l, s, r, g, b);

         ULong_t pixel = TColor::RGB2Pixel(r, g, b);

         for (Int_t x = 0; x < width; ++x) {
            gVirtualX->PutPixel(fLimage, x, y, pixel);
         }
      }
   } else {
      CreateDitheredImage(fLimage, kIMG_L);
   }

   gClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Position the slider cursor on right color position.

void TGColorPick::SetColor(ULong_t color)
{
   UInt_t width, height;
   Int_t h, l, s;
   Int_t r, g, b;

   gVirtualX->GetImageSize(fHSimage, width, height);

   fCurrentColor = color;

   TColor::Pixel2RGB(fCurrentColor, r, g, b);
   TColor::RGB2HLS(r, g, b, h, l, s);

   SetHScursor(h * (Int_t)width / 256, (255 - s) * (Int_t)height / 256);

   gVirtualX->GetImageSize(fLimage, width, height);

   SetLcursor((255 - l) * (Int_t)height / 256);

   SetSliderColor();
}

////////////////////////////////////////////////////////////////////////////////
/// Assign the current cursor position as currently selected color.

void TGColorPick::UpdateCurrentColor()
{
   UInt_t lwidth, lheight;
   UInt_t swidth, sheight;
   Int_t r, g, b;
   Int_t h, l, s;

   gVirtualX->GetImageSize(fLimage, lwidth, lheight);
   gVirtualX->GetImageSize(fHSimage, swidth, sheight);

   h = Int_t(fCx * 255 / swidth);
   l = Int_t((lheight - fCz) * 255 / lheight);
   s = Int_t((sheight - fCy) * 255 / sheight);

   TColor::HLS2RGB(h, l, s, r, g, b);
   fCurrentColor = TColor::RGB2Pixel(r, g, b);
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw the color pick widget.

void TGColorPick::DoRedraw()
{
   UInt_t lwidth, lheight;
   UInt_t swidth, sheight;

   gVirtualX->GetImageSize(fLimage, lwidth, lheight);
   gVirtualX->GetImageSize(fHSimage, swidth, sheight);

   DrawBorder();

   Draw3dRectangle(kSunkenFrame, fColormapRect.fX - 1, fColormapRect.fY - 1,
                   fColormapRect.fWidth + 2, fColormapRect.fHeight + 2);
   gVirtualX->PutImage(fId, GetBckgndGC()(), fHSimage,
                       fColormapRect.fX, fColormapRect.fY, 0, 0, swidth, sheight);

   Draw3dRectangle(kSunkenFrame, fSliderRect.fX - 1, fSliderRect.fY - 1,
                   fSliderRect.fWidth + 2, fSliderRect.fHeight + 2);
   gVirtualX->PutImage(fId, GetBckgndGC()(), fLimage,
                       fSliderRect.fX, fSliderRect.fY, 0, 0, lwidth, lheight);

   DrawHScursor(kTRUE);
   DrawLcursor(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Set hue / saturation cursor position.

void TGColorPick::SetHScursor(Int_t x, Int_t y)
{
   UInt_t width, height;

   gVirtualX->GetImageSize(fHSimage, width, height);

   DrawHScursor(kFALSE);

   fCx = x;
   fCy = y;

   if (fCx < 0)
      fCx = 0;
   else if (fCx >= (Int_t)width)
      fCx = (Int_t)width - 1;

   if (fCy < 0)
      fCy = 0;
   else if (fCy >= (Int_t)height)
      fCy = (Int_t)height - 1;

   DrawHScursor(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Set lightness slider cursor position.

void TGColorPick::SetLcursor(Int_t z)
{
   UInt_t width, height;

   gVirtualX->GetImageSize(fLimage, width, height);

   DrawLcursor(kFALSE);

   fCz = z - fSliderRect.fY;

   if (fCz < 0)
      fCz = 0;
   else if (fCz >= (Int_t)height)
      fCz = (Int_t)height - 1;

   DrawLcursor(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw hue / saturation cursor

void TGColorPick::DrawHScursor(Int_t onoff)
{
   UInt_t width, height;

   gVirtualX->GetImageSize(fHSimage, width, height);

   if (onoff) {
      Int_t x, y;
      Rectangle_t rect;

      x = fCx + fColormapRect.fX;
      y = fCy + fColormapRect.fY;

      rect.fX = fColormapRect.fX;
      rect.fY = fColormapRect.fX;
      rect.fWidth = fColormapRect.fWidth;
      rect.fHeight = fColormapRect.fHeight;
      gVirtualX->SetClipRectangles(fCursorGC(), 0, 0, &rect, 1);

      gVirtualX->FillRectangle(fId, fCursorGC(), x - 9, y - 1, 5, 3);
      gVirtualX->FillRectangle(fId, fCursorGC(), x - 1, y - 9, 3, 5);
      gVirtualX->FillRectangle(fId, fCursorGC(), x + 5, y - 1, 5, 3);
      gVirtualX->FillRectangle(fId, fCursorGC(), x - 1, y + 5, 3, 5);

   } else {
      Int_t  x, y;
      UInt_t w, h;

      x = fCx - 9; w = 19;
      y = fCy - 9; h = 19;

      if (x < 0) { w += x; x = 0; }
      if (y < 0) { h += y; y = 0; }

      if (x + w > width) w = width - x;
      if (y + h > height) h = height - y;

      gVirtualX->PutImage(fId, GetBckgndGC()(), fHSimage, x, y,
                          fColormapRect.fX + x, fColormapRect.fY + y, w, h);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw lightness slider cursor

void TGColorPick::DrawLcursor(Int_t onoff)
{
   Int_t l = fSliderRect.fX + fSliderRect.fWidth + 3;
   Int_t r = l + 5;
   Int_t t = fCz - 5 + fSliderRect.fY;
   Int_t b = t + 10;

   Point_t points[3];

   Int_t m = (t + b) >> 1;

   points[0].fX = r;
   points[0].fY = t;
   points[1].fX = r;
   points[1].fY = b;
   points[2].fX = l;
   points[2].fY = m;

   GContext_t gc = onoff ? GetShadowGC()() : GetBckgndGC()();

   gVirtualX->FillPolygon(fId, gc, points, 3);
}


////////////////////////////////////////////////////////////////////////////////
/// Color selection dialog constructor.
/// The TGColorDialog presents a full featured color selection dialog.
/// It uses 2 TGColorPalette's and the TGColorPick widgets.

TGColorDialog::TGColorDialog(const TGWindow *p, const TGWindow *m,
                             Int_t *retc, ULong_t *color, Bool_t wait) :
   TGTransientFrame(p, m, 200, 150)
{
   const Int_t kC_X = 175;  // Win95: 177
   const Int_t kC_Y = 180;  // Win95: 189

   Int_t i;

   fRetc = retc;
   fRetColor = 0;
   fRetTColor = 0;
   fInitColor = 0;
   if (color) {
      fRetColor = color;
      fRetTColor = gROOT->GetColor(TColor::GetColor(*color));
      fInitColor = *fRetColor;
   }
   fWaitFor = wait;

   if (fRetc) *fRetc = kMBCancel;

   TGHorizontalFrame *hftop = new TGHorizontalFrame(this, 10, 10);
   hftop->SetCleanup();
   AddFrame(hftop, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 10, 5));

   fTab = new TGTab(hftop, 300, 300);
   hftop->AddFrame(fTab);

   TGCompositeFrame *cf = new TGCompositeFrame(hftop, 10, 10);
   cf->SetCleanup();
   hftop->AddFrame(cf, new TGLayoutHints(kLHintsLeft | kLHintsTop, 5, 0, 30, 0));

   TGCompositeFrame *cf1 = new TGCompositeFrame(cf, 10, 10);
   cf1->SetCleanup();
   cf->AddFrame(cf1, new TGLayoutHints(kLHintsLeft | kLHintsTop, 5, 0, 30, 0));
   cf1->SetLayoutManager(new TGMatrixLayout(cf1, 0, 2, 4));

   cf1->AddFrame(new TGLabel(cf1, new TGHotString("Red:")));
   cf1->AddFrame(fRte = new TGTextEntry(cf1, fRtb = new TGTextBuffer(5), kCDLG_RTE),0);
   fRte->Resize(50, fRte->GetDefaultHeight());
   cf1->AddFrame(new TGLabel(cf1, new TGHotString("Green:")),0);
   cf1->AddFrame(fGte = new TGTextEntry(cf1, fGtb = new TGTextBuffer(5), kCDLG_GTE),0);
   fGte->Resize(50, fGte->GetDefaultHeight());
   cf1->AddFrame(new TGLabel(cf1, new TGHotString("Blue:")));
   cf1->AddFrame(fBte = new TGTextEntry(cf1, fBtb = new TGTextBuffer(5), kCDLG_BTE),0);
   fBte->Resize(50, fBte->GetDefaultHeight());
   cf1->AddFrame(new TGLabel(cf1, new TGHotString("Opacity:")),0);
   cf1->AddFrame(fAle = new TGTextEntry(cf1, fAlb = new TGTextBuffer(5), kCDLG_ALE),0);
   fAle->Resize(50, fAle->GetDefaultHeight());

   if (!TCanvas::SupportAlpha()) {
      fAle->SetEnabled(kFALSE);
   }

   TGCompositeFrame *cf2 = new TGCompositeFrame(cf, 10, 10);
   cf2->SetCleanup();
   cf->AddFrame(cf2, new TGLayoutHints(kLHintsLeft | kLHintsTop, 5, 0, 30, 0));
   cf2->SetLayoutManager(new TGMatrixLayout(cf2, 0, 2, 4));
   cf2->AddFrame(new TGLabel(cf2, new TGHotString("Hue:")),0);
   cf2->AddFrame(fHte = new TGTextEntry(cf2, fHtb = new TGTextBuffer(5), kCDLG_HTE),0);
   fHte->Resize(50, fHte->GetDefaultHeight());
   cf2->AddFrame(new TGLabel(cf2, new TGHotString("Sat:")),0);
   cf2->AddFrame(fSte = new TGTextEntry(cf2, fStb = new TGTextBuffer(5), kCDLG_STE),0);
   fSte->Resize(50, fSte->GetDefaultHeight());
   cf2->AddFrame(new TGLabel(cf2, new TGHotString("Lum:")),0);
   cf2->AddFrame(fLte = new TGTextEntry(cf2, fLtb = new TGTextBuffer(5), kCDLG_LTE),0);
   fLte->Resize(50, fLte->GetDefaultHeight());

   fHte->Associate(this);
   fLte->Associate(this);
   fSte->Associate(this);
   fRte->Associate(this);
   fGte->Associate(this);
   fBte->Associate(this);
   fAle->Associate(this);

   if (color) {
      UpdateRGBentries(color);
      UpdateHLSentries(color);
      UpdateAlpha(color);
      fCurrentColor = *color;
   } else {
      gClient->GetColorByName("red", fCurrentColor);
   }

   // color sample
   TGCompositeFrame *cf3 = new TGCompositeFrame(cf, 10, 10);
   cf3->SetCleanup();
   cf3->SetLayoutManager(new TGMatrixLayout(cf3, 0, 1, 0));
   cf3->AddFrame(fColorInfo = new TGLabel(cf3, new TGString("New: not set         ")),0);
   fColorInfo->SetTextJustify(kTextLeft);
   cf3->AddFrame(fSample = new TGFrame(cf3, 50, 25, kOwnBackground),0);
   cf3->AddFrame(fSampleOld = new TGFrame(cf3, 50, 25, kOwnBackground),0);
   cf3->AddFrame(new TGLabel(cf3, new TGString("Current")),0);
   cf->AddFrame(cf3, new TGLayoutHints(kLHintsLeft | kLHintsTop, 5, 5, 20, 0));
   fSample->SetBackgroundColor(fCurrentColor);
   fSampleOld->SetBackgroundColor(fCurrentColor);

   TGCompositeFrame *tf = fTab->AddTab("Color Wheel");
   TGCompositeFrame *tf1 = new TGCompositeFrame(tf, 60, 20, kHorizontalFrame);
   tf->AddFrame(tf1);
   fEcanvas = new TRootEmbeddedCanvas("wheel", tf1, 360, 360);
   tf1->AddFrame(fEcanvas);
   TCanvas *wcan = fEcanvas->GetCanvas();
   wcan->SetBit(kNoContextMenu);
   fColorWheel = new TColorWheel();
   fColorWheel->SetCanvas(wcan);
   fColorWheel->Draw();
   wcan->Update();
   wcan->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","TGColorDialog",this,
                 "SetColorInfo(Int_t,Int_t,Int_t,TObject*)");

   tf = fTab->AddTab("Basic Colors");
   TGCompositeFrame *tf2 = new TGCompositeFrame(tf, 60, 20, kHorizontalFrame);
   tf->AddFrame(tf2);

   TGVerticalFrame *vf1 = new TGVerticalFrame(tf2, 20, 20);
   vf1->SetCleanup();
   TGVerticalFrame *vf2 = new TGVerticalFrame(tf2, 20, 20);
   vf2->SetCleanup();

   tf2->AddFrame(vf1, new TGLayoutHints(kLHintsLeft | kLHintsExpandY));
   tf2->AddFrame(vf2, new TGLayoutHints(kLHintsLeft | kLHintsExpandY));

   //----------------------------------------- Left panel

   // basic colors

   fPalette = new TGColorPalette(vf1, 6, 8, kCDLG_SPALETTE);
   vf1->AddFrame(fPalette, new TGLayoutHints(kLHintsNormal, 5, 5, 15, 0));
   fPalette->Associate(this);

   for (i = 0; i < 48; ++i)
      fPalette->SetColor(i, TColor::Number2Pixel(i+10));  // root colors
      // the basic colors were set via bcolor
      //fPalette->SetColor(i, TColor::GetPixel(bcolor[i][0], bcolor[i][1], bcolor[i][2]));

   // add some default colors
   fPalette->SetColor(47, TGFrame::GetDefaultFrameBackground());

   Float_t r, g, b;

   r = 232./255;
   g = 232./255;
   b = 222./255;

   // Gui Builder background
   Pixel_t pixel = TColor::RGB2Pixel(r, g, b);
   fPalette->SetColor(46, pixel);

   r = 230./255;
   g = 230./255;
   b = 230./255;

   // a la MAC background
   pixel = TColor::RGB2Pixel(r, g, b);
   fPalette->SetColor(45, pixel);

   r = 172./255;
   g = 174./255;
   b = 205./255;

   // a la CDE background
   pixel = TColor::RGB2Pixel(r, g, b);
   fPalette->SetColor(44, pixel);

   r = 205./255;
   g = 195./255;
   b = 175./255;

   // a la FOX background
   pixel = TColor::RGB2Pixel(r, g, b);
   fPalette->SetColor(43, pixel);

   // custom colors

   vf1->AddFrame(new TGLabel(vf1, new TGHotString("&Custom Colors:")),
                 new TGLayoutHints(kLHintsNormal, 5, 0, 15, 2));

   fCpalette = new TGColorPalette(vf1, 6, 4, kCDLG_CPALETTE);
   vf1->AddFrame(fCpalette, new TGLayoutHints(kLHintsNormal, 5, 5, 5, 0));
   fCpalette->Associate(this);

   if (gUcolor[0] == 0xff000000) {
      for (i = 0; i < 24; i++)
         gUcolor[i] = TColor::RGB2Pixel(255, 255, 255);
   }
   fCpalette->SetColors(gUcolor);

   // button frame - OK, Cancel
   TGHorizontalFrame *hf = new TGHorizontalFrame(this, 10, 10, kFixedWidth);
   hf->SetCleanup();
   AddFrame(hf, new TGLayoutHints(kLHintsBottom | kLHintsRight, 5, 5, 10, 5));

   TGTextButton *ok = new TGTextButton(hf, new TGHotString("OK"), kCDLG_OK);
   TGTextButton *cancel = new TGTextButton(hf, new TGHotString("Cancel"), kCDLG_CANCEL);
   fPreview = new TGTextButton(hf, new TGHotString("&Preview"), kCDLG_PREVIEW);
   fPreview->Connect("Clicked()", "TGColorDialog", this, "DoPreview()");

   hf->AddFrame(ok, new TGLayoutHints(kLHintsBottom | kLHintsExpandX, 0, 3, 0, 0));
   hf->AddFrame(cancel, new TGLayoutHints(kLHintsBottom | kLHintsExpandX,3, 0, 0, 0));
   hf->AddFrame(fPreview, new TGLayoutHints(kLHintsBottom | kLHintsExpandX,3, 0, 0, 0));

   UInt_t w = ok->GetDefaultWidth();
   w = TMath::Max(w, cancel->GetDefaultWidth());
   hf->Resize(3 * (w + 30), hf->GetDefaultHeight());

   ok->Associate(this);
   cancel->Associate(this);

   //----------------------------------------- Right panel

   // fColormap frame

   fColors = new TGColorPick(vf2, kC_X + 23, kC_Y, kCDLG_COLORPICK);
   vf2->AddFrame(fColors, new TGLayoutHints(kLHintsLeft | kLHintsTop, 5, 0, 15, 5));
   fColors->Associate(this);

   if (color)
      fColors->SetColor(*color);

   TGTextButton *add = new TGTextButton(vf2, new TGHotString("&Add to Custom Colors"),
                                        kCDLG_ADD);
   vf2->AddFrame(add, new TGLayoutHints(kLHintsBottom | kLHintsExpandX,
                                        5, 10, 0, 5));
   add->Associate(this);

   MapSubwindows();
   Resize(GetDefaultSize());
   SetEditDisabled(kEditDisable);

   //---- make the message box non-resizable

   SetWMSize(fWidth, fHeight);
   SetWMSizeHints(fWidth, fHeight, fWidth, fHeight, 0, 0);

   SetWindowName("Color Selector");
   SetIconName("Color Selector");
   SetClassHints("ROOT", "ColorSelector");

   SetMWMHints(kMWMDecorAll | kMWMDecorResizeH | kMWMDecorMaximize |
                              kMWMDecorMinimize | kMWMDecorMenu,
               kMWMFuncAll  | kMWMFuncResize    | kMWMFuncMaximize |
                              kMWMFuncMinimize,
               kMWMInputModeless);


   //---- position relative to the parent's window

   if (fClient->IsEditable()) {
      const TGWindow *main = fMain;
      fMain = fClient->GetRoot();
      CenterOnParent(kTRUE, TGTransientFrame::kRight);
      fMain = main;
   } else {
      CenterOnParent();
   }

   if (fWaitFor) {
      MapWindow();
      fClient->WaitForUnmap(this);
      DeleteWindow();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// TGColorDialog destructor.

TGColorDialog::~TGColorDialog()
{
   fEcanvas->GetCanvas()->Disconnect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)");
   delete fEcanvas;
   Cleanup();
}

////////////////////////////////////////////////////////////////////////////////
/// Change current color.

void TGColorDialog::SetCurrentColor(Pixel_t col)
{
   if (fCurrentColor == col) {
      return;
   }
   fInitColor = *fRetColor = col;
   if((fRetTColor = gROOT->GetColor(TColor::GetColor(col)))) {};
   fCurrentColor = col;
   fColors->SetColor(col);
   fSample->SetBackgroundColor(col);
   ColorSelected(col);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit signal about selected color.

void TGColorDialog::ColorSelected(Pixel_t color)
{
   Emit("ColorSelected(Pixel_t)", color);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit signal about selected alpha and color.

void TGColorDialog::AlphaColorSelected(ULong_t color)
{
   Emit("AlphaColorSelected(ULong_t)", color);
}

////////////////////////////////////////////////////////////////////////////////
/// Called when window is closed via window manager.

void TGColorDialog::CloseWindow()
{
   // save user set colors
   for (Int_t i = 0; i < 24; ++i)
      gUcolor[i] = fCpalette->GetColorByIndex(i);

   if (*fRetc != kMBOk) {
      ColorSelected(fInitColor);
      ULong_t ptr;
      if((ptr = (ULong_t)gROOT->GetColor(TColor::GetColor(fInitColor)))) AlphaColorSelected(ptr);
   } else {
      ColorSelected(*fRetColor);
      AlphaColorSelected((ULong_t)fRetTColor);
   }
   // don't call DeleteWindow() here since that will cause access
   // to the deleted dialog in the WaitFor() method (see ctor)

   //OpenGL + XQuartz on Mac: gl context and related resources
   //must be deleted _before_ UnmapWindow.
   if (gVirtualX->InheritsFrom("TGX11") && fEcanvas->GetCanvas()->UseGL())
      fEcanvas->GetCanvas()->DeleteCanvasPainter();

   UnmapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Upadate Opacity text entry with alpha value of color c.

void TGColorDialog::UpdateAlpha(ULong_t *c)
{
   Char_t tmp[20];
   Double_t alpha;

   if (TColor *color = gROOT->GetColor(TColor::GetColor(*c))) {
      alpha = color->GetAlpha();
      snprintf(tmp, 20, "%.1f", alpha);
      fAlb->Clear();
      fAlb->AddText(0,tmp);
      gClient->NeedRedraw(fAle);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Update RGB text entries with RGB values of color c.

void TGColorDialog::UpdateRGBentries(ULong_t *c)
{
   Char_t tmp[20];

   Int_t r, g, b;
   TColor::Pixel2RGB(*c, r, g, b);

   snprintf(tmp, 20, "%d", r);
   fRtb->Clear();
   fRtb->AddText(0, tmp);
   gClient->NeedRedraw(fRte);

   snprintf(tmp, 20, "%d", g);
   fGtb->Clear();
   fGtb->AddText(0, tmp);
   gClient->NeedRedraw(fGte);

   snprintf(tmp, 20, "%d", b);
   fBtb->Clear();
   fBtb->AddText(0, tmp);
   gClient->NeedRedraw(fBte);
}

////////////////////////////////////////////////////////////////////////////////
/// Update HLS text entries with HLS values of color c.

void TGColorDialog::UpdateHLSentries(ULong_t *c)
{
   Char_t tmp[20];

   Int_t h, l, s;
   Int_t r, g, b;

   TColor::Pixel2RGB(*c, r, g, b);
   TColor::RGB2HLS(r, g, b, h, l, s);

   snprintf(tmp, 20, "%d", h);
   fHtb->Clear();
   fHtb->AddText(0, tmp);
   gClient->NeedRedraw(fHte);

   snprintf(tmp, 20, "%d", l);
   fLtb->Clear();
   fLtb->AddText(0, tmp);
   gClient->NeedRedraw(fLte);

   snprintf(tmp, 20, "%d", s);
   fStb->Clear();
   fStb->AddText(0, tmp);
   gClient->NeedRedraw(fSte);
}

////////////////////////////////////////////////////////////////////////////////
/// Process messages for the color selection dialog.

Bool_t TGColorDialog::ProcessMessage(Long_t msg, Long_t parm1, Long_t /*parm2*/)
{
   ULong_t color;
   Int_t h, l, s;
   Int_t r, g, b;

   switch (GET_MSG(msg)) {
      case kC_COMMAND:
         switch (GET_SUBMSG(msg)) {
            case kCM_BUTTON:
               switch(parm1) {
                  case kCDLG_ADD:
                     fCpalette->SetCurrentCellColor(fCurrentColor);
                     break;

                  case kCDLG_OK:
                     *fRetc = kMBOk;
                     *fRetColor = TColor::RGB2Pixel(atoi(fRtb->GetString()),
                                                    atoi(fGtb->GetString()),
                                                    atoi(fBtb->GetString()));
                     if ((fRetTColor = gROOT->GetColor(TColor::GetColor(*fRetColor)))) {
                        fRetTColor->SetAlpha(TMath::Max((Double_t)0, TMath::Min((Double_t)1,
                                             atof(fAlb->GetString()))));
                     }
                     CloseWindow();
                     break;
                  case kCDLG_CANCEL:
                     if (!fClient->IsEditable()) {
                        TGColorPopup *p = (TGColorPopup *)GetMain();
                        if (p && p->InheritsFrom("TGColorPopup"))
                           p->PreviewColor(fSampleOld->GetBackground());
                     }
                     CloseWindow();
                     break;
               }
               break;
         }
         break;
      case kC_COLORSEL:
         switch (GET_SUBMSG(msg)) {
            case kCOL_CLICK:
               switch (parm1) {
                  case kCDLG_SPALETTE:
                     color = fPalette->GetCurrentColor();
                     fSample->SetBackgroundColor(color);
                     ColorSelected(color);
                     gClient->NeedRedraw(fSample);
                     fCurrentColor = color;
                     fColors->SetColor(color);
                     UpdateRGBentries(&color);
                     UpdateHLSentries(&color);
                     UpdateAlpha(&color);
                     break;

                  case kCDLG_CPALETTE:
                     color = fCpalette->GetCurrentColor();
                     fSample->SetBackgroundColor(color);
                     ColorSelected(color);
                     gClient->NeedRedraw(fSample);
                     fCurrentColor = color;
                     fColors->SetColor(color);
                     UpdateRGBentries(&color);
                     UpdateHLSentries(&color);
                     UpdateAlpha(&color);
                     break;

                  case kCDLG_COLORPICK:
                     color = fColors->GetCurrentColor();
                     fSample->SetBackgroundColor(color);
                     ColorSelected(color);
                     gClient->NeedRedraw(fSample);
                     fCurrentColor = color;
                     UpdateRGBentries(&color);
                     UpdateHLSentries(&color);
                     UpdateAlpha(&color);
                     break;

               }
               break;
         }
         break;

      case kC_TEXTENTRY:
         switch (GET_SUBMSG(msg)) {
            case kTE_TEXTCHANGED:
               switch (parm1) {
                  case kCDLG_HTE:
                  case kCDLG_LTE:
                  case kCDLG_STE:

                     h = atoi(fHtb->GetString());
                     l = atoi(fLtb->GetString());
                     s = atoi(fStb->GetString());
                     TColor::HLS2RGB(h, l, s, r, g, b);

                     color = TColor::RGB2Pixel(r, g, b);
                     fSample->SetBackgroundColor(color);
                     ColorSelected(color);
                     gClient->NeedRedraw(fSample);
                     fCurrentColor = color;
                     fColors->SetColor(color);
                     UpdateRGBentries(&color);
                     break;

                  case kCDLG_RTE:
                  case kCDLG_GTE:
                  case kCDLG_BTE:
                     color = TColor::RGB2Pixel(atoi(fRtb->GetString()),
                                               atoi(fGtb->GetString()),
                                               atoi(fBtb->GetString()));
                     fSample->SetBackgroundColor(color);
                     ColorSelected(color);
                     gClient->NeedRedraw(fSample);
                     fCurrentColor = color;
                     fColors->SetColor(color);
                     UpdateHLSentries(&color);
                     break;

               }
               break;
         }
         break;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///  Set the color info in RGB and HLS parts

void TGColorDialog::SetColorInfo(Int_t event, Int_t px, Int_t py, TObject *object)
{
   if (object == fColorWheel) {
      Int_t n = fColorWheel->GetColor(px,py);
      if (n < 0) return;
      TColor *color = gROOT->GetColor(n);
      if (!color) return;
      ULong_t pcolor = color->GetPixel();
      if (event == kButton1Down) {
         UpdateRGBentries(&pcolor);
         UpdateHLSentries(&pcolor);
         UpdateAlpha(&pcolor);
         fSample->SetBackgroundColor(pcolor);
         fColorInfo->SetText(Form("New: %s",color->GetName()));
         gClient->NeedRedraw(fSample);
         gClient->NeedRedraw(fColorInfo);
         fCurrentColor = pcolor;
         fColors->SetColor(pcolor);
         ColorSelected(pcolor);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot method called when Preview button is clicked.

void TGColorDialog::DoPreview()
{
   TColor *tcolor;
   if ((tcolor = gROOT->GetColor(TColor::GetColor(fSample->GetBackground())))) {
      tcolor->SetAlpha(TMath::Max((Double_t)0, TMath::Min((Double_t)1, atof(fAlb->GetString()))));
   }

   if (fClient->IsEditable()) {
      ColorSelected(fSample->GetBackground());
      AlphaColorSelected((ULong_t)tcolor);
      return;
   }
   TGColorPopup *p = (TGColorPopup *)GetMain();
   if (p && p->InheritsFrom("TGColorPopup")) {
      if (tcolor) p->PreviewAlphaColor((ULong_t)tcolor);
      else p->PreviewColor(fSample->GetBackground());
   }
}
