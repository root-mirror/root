// $Id: TGHtml.cxx,v 1.4 2007/05/07 15:19:07 brun Exp $
// Author:  Valeriy Onuchin   03/05/2007

/*************************************************************************
 * Copyright (C) 1995-2001, Rene Brun, Fons Rademakers and Reiner Rohlfs *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**************************************************************************

    HTML widget for xclass. Based on tkhtml 1.28
    Copyright (C) 1997-2000 D. Richard Hipp <drh@acm.org>
    Copyright (C) 2002-2003 Hector Peraza.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Library General Public License for more details.

    You should have received a copy of the GNU Library General Public
    License along with this library; if not, write to the Free
    Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

**************************************************************************/

#include <cctype>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "TGHtml.h"
#include "THashTable.h"
#include "TObjString.h"
#include "TGIdleHandler.h"
#include "TImage.h"
#include "TGTextEntry.h"
#include "TGText.h"
#include "TGComboBox.h"
#include "TGListBox.h"
#include "TVirtualX.h"
#include "strlcpy.h"
#include "snprintf.h"

//_____________________________________________________________________________
//
// TGHtml
//
// The ROOT HTML widget. A derivate of TGView.
//_____________________________________________________________________________

ClassImp(TGHtml);

int HtmlTraceMask = 0; //HtmlTrace_Table1 | HtmlTrace_Table4;
int HtmlDepth = 0;

#define DEF_FRAME_BG_COLOR        "#c0c0c0"
#define DEF_FRAME_CURSOR          ""
#define DEF_BUTTON_FG             "black"
#define DEF_BUTTON_HIGHLIGHT_BG   "#d9d9d9"
#define DEF_BUTTON_HIGHLIGHT      "black"


////////////////////////////////////////////////////////////////////////////////
/// HTML Widget constructor.

TGHtml::TGHtml(const TGWindow *p, int w, int h, int id) : TGView(p, w, h, id)
{
   SetBit(kIsHtmlView);

   int i;

   fExiting = 0;
   fPFirst = 0;
   fPLast = 0;
   fNToken = 0;
   fLastSized = 0;
   fNextPlaced = 0;
   fFirstBlock = 0;
   fLastBlock = 0;
   fFirstInput = 0;
   fLastInput = 0;
   fNInput = 0;
   fNForm = 0;
   fVarId = 0;  // do we need this??
   fInputIdx = 0;
   fRadioIdx = 0;
   fSelBegin.fI = 0;
   fSelBegin.fP = 0;
   fSelEnd.fI = 0;
   fSelEnd.fP = 0;
   fPSelStartBlock = 0;
   fPSelEndBlock = 0;
   fInsOnTime = DEF_HTML_INSERT_ON_TIME;
   fInsOffTime = DEF_HTML_INSERT_OFF_TIME;
   fInsStatus = 0;
   fInsTimer = 0;
   fIns.fI = 0;
   fIns.fP = 0;
   fPInsBlock = 0;
   fInsIndex = 0;
   fZText = 0;
   fNText = 0;
   fNAlloc = 0;
   fNComplete = 0;
   fICol = 0;
   fIPlaintext = 0;
   fPScript = 0;
   fIdle = 0;
   fStyleStack = 0;
   fParaAlignment = ALIGN_None;
   fRowAlignment = ALIGN_None;
   fAnchorFlags = 0;
   fInDt = 0;
   fInTr = 0;
   fInTd = 0;
   fAnchorStart = 0;
   fFormStart = 0;
   fFormElemStart = 0;
   fFormElemLast = 0;
   fLoEndPtr = 0;
   fLoFormStart = 0;
   fInnerList = 0;
   ResetLayoutContext();
   fHighlightWidth = 0;
   fHighlightBgColorPtr = 0;
   fHighlightColorPtr = 0;
   for (i = 0; i < N_FONT; ++i) fAFont[i] = 0;
   memset(fFontValid, 0, sizeof(fFontValid));
   for (i = 0; i < N_COLOR; ++i) {
      fApColor[i] = 0;
      fIDark[i] = 0;
      fILight[i] = 0;
   }
   fFgColor = AllocColor("black");
   fBgColor = AllocColor("white"); //AllocColor("#c0c0c0");
   fNewLinkColor = AllocColor(DEF_HTML_UNVISITED);
   fOldLinkColor = AllocColor(DEF_HTML_VISITED);
   fSelectionColor = AllocColor(DEF_HTML_SELECTION_COLOR);

   fApColor[COLOR_Normal] = fFgColor;
   fApColor[COLOR_Visited] = fOldLinkColor;
   fApColor[COLOR_Unvisited] = fNewLinkColor;
   fApColor[COLOR_Selection] = fSelectionColor;
   fApColor[COLOR_Background] = fBgColor;

   fBgImage = 0;

   SetBackgroundColor(fApColor[COLOR_Background]->fPixel);
   SetBackgroundPixmap(0);  // force usage of solid color

   fColorUsed = 0;

   for (i = 0; i < N_CACHE_GC; ++i) {
      fAGcCache[i].fIndex = 0;
      fAGcCache[i].fColor = 0;
      fAGcCache[i].fFont  = 0;
      fAGcCache[i].fGc    = 0;
   }
   fLastGC = 0;
   fSelEndIndex =0;
   fSelStartIndex = 0;
   fGcNextToFree = 0;
   fImageList = 0;
   fZBaseHref = 0;
   fInnerList = 0;
   fFormPadding = 5;
   fOverrideFonts = 0;
   fOverrideColors = 0;
   fHasScript = 0;
   fHasFrames = 0;
   fAddEndTags = 0;
   fTableBorderMin = 0;
   fVarind = 0;
   fIdind = 0;
   fInParse = 0;
   fZGoto = 0;
   fExts = 0;
   fUnderlineLinks = kTRUE;
   fExportSelection = DEF_HTML_EXPORT_SEL;
   fTableRelief = HTML_RELIEF_RAISED;
   fRuleRelief = HTML_RELIEF_SUNKEN;
   fRulePadding = 5;
   fZBase = 0;
   fZBaseHref = 0;
   fCursor = kPointer;
   fMaxX = 0;
   fMaxY = 0;

   fXMargin = fYMargin = 0; //HTML_INDENT/4;

   fFlags = RESIZE_ELEMENTS | RELAYOUT;

   fDirtyLeft = LARGE_NUMBER;
   fDirtyRight = 0;
   fDirtyTop = LARGE_NUMBER;
   fDirtyBottom = 0;


   fVsb->SetAccelerated();
   fHsb->SetAccelerated();

   fLastUri = 0;

   AddInput(kExposureMask | kFocusChangeMask);
   AddInput(kButtonPressMask | kButtonReleaseMask | kPointerMotionMask);

   fUidTable = new THashTable(100);
}

////////////////////////////////////////////////////////////////////////////////
/// HTML widget destructor.

TGHtml::~TGHtml()
{
   int i;

   fExiting = 1;
   HClear();
   for (i = 0; i < N_FONT; i++) {
      if (fAFont[i] != 0) fClient->FreeFont(fAFont[i]);
   }
   if (fInsTimer) delete fInsTimer;
   if (fIdle) delete fIdle;

  // TODO: should also free colors!
}

////////////////////////////////////////////////////////////////////////////////
/// Start background update.

void TGHtml::UpdateBackgroundStart()
{
   //GCValues_t gcv;
   //unsigned int mask = GCTileStipXOrigin | GCTileStipYOrigin;
//
   //gcv.ts_x_origin = -fVvisible.fX;
   //gcv.ts_y_origin = -fVisibleStart.fY;
   //XChangeGC(GetDisplay(), _backGC, mask, &gcv);
}

////////////////////////////////////////////////////////////////////////////////
/// Free system color.

void TGHtml::FreeColor(ColorStruct_t *color)
{
   gVirtualX->FreeColor(gClient->GetDefaultColormap(), color->fPixel);
   delete color;
}

////////////////////////////////////////////////////////////////////////////////
/// Allocate system color by name.

ColorStruct_t *TGHtml::AllocColor(const char *name)
{
   ColorStruct_t *color = new ColorStruct_t;

   color->fPixel = 0;
   if (gVirtualX->ParseColor(fClient->GetDefaultColormap(), name, *color)) {
      if (!gVirtualX->AllocColor(fClient->GetDefaultColormap(), *color)) {
         // force allocation of pixel 0
         gVirtualX->QueryColor(fClient->GetDefaultColormap(), *color);
         gVirtualX->AllocColor(fClient->GetDefaultColormap(), *color);
      }
   }

   return color;
}

////////////////////////////////////////////////////////////////////////////////
/// Allocate system color by value.

ColorStruct_t *TGHtml::AllocColorByValue(ColorStruct_t *color)
{
   ColorStruct_t *c = new ColorStruct_t;
   *c = *color;

   if (!gVirtualX->AllocColor(gClient->GetDefaultColormap(), *c)) {
      // force allocation of pixel 0
      c->fPixel = 0;
      gVirtualX->QueryColor(gClient->GetDefaultColormap(), *c);
      gVirtualX->AllocColor(gClient->GetDefaultColormap(), *c);
   }

   return c;
}

////////////////////////////////////////////////////////////////////////////////
/// Erase all HTML from this widget and clear the screen. This is
/// typically done before loading a new document.

void TGHtml::Clear(Option_t *)
{
   HClear();
   TGView::Clear();
   fFlags |= REDRAW_TEXT | VSCROLL | HSCROLL;
   ScheduleRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Appends (or insert at the specified position) the given HTML text to the
/// end of any HTML text that may have been inserted by prior calls to this
/// command.  Then it runs the tokenizer, parser and layout engine as far as
/// possible with the text that is available. The display is updated
/// appropriately.

int TGHtml::ParseText(char *text, const char *index)
{
   SHtmlIndex_t iStart;
   TGHtmlElement *savePtr=0;

   iStart.fP = 0;
   iStart.fI = 0;

   fLoEndPtr = fPLast;

   if (index) {
      int rc = GetIndex(index, &iStart.fP, &iStart.fI);
      if (rc != 0) return kFALSE;  // malformed index
      if (iStart.fP) {
         savePtr = iStart.fP->fPNext;
         fPLast = iStart.fP;
         iStart.fP->fPNext = 0;
      }
   }

   TokenizerAppend(text);

   if (fLoEndPtr) {
      fFormStart = fLoFormStart;
      if (iStart.fP && savePtr) {
         AddStyle(fLoEndPtr);
         fPLast->fPNext = savePtr;
         savePtr->fPPrev = fPLast;
         fPLast = fLoEndPtr;
         fFlags |= REDRAW_TEXT | RELAYOUT;
         ScheduleRedraw();
      } else if (fLoEndPtr->fPNext) {
         AddStyle(fLoEndPtr->fPNext);
      }
   } else if (fPFirst) {
      fParaAlignment = ALIGN_None;
      fRowAlignment = ALIGN_None;
      fAnchorFlags = 0;
      fInDt = 0;
      fAnchorStart = 0;
      fFormStart = 0;
      fInnerList = 0;
      fNInput = 0;
      AddStyle(fPFirst);
   }
#if 1
   fLoEndPtr = fPLast;
   fLoFormStart = fFormStart;
#endif

   fFlags |= EXTEND_LAYOUT;
   ScheduleRedraw();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets relief mode of html table.

void TGHtml::SetTableRelief(int relief)
{
   if (fTableRelief != relief) {
      fTableRelief = relief;
      fFlags |= RELAYOUT;
      RedrawEverything();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Sets relief mode of html rule.

void TGHtml::SetRuleRelief(int relief)
{
   if (fRuleRelief != relief) {
      fRuleRelief = relief;
      fFlags |= RELAYOUT;
      RedrawEverything();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set/reset html links underline.

void TGHtml::UnderlineLinks(int onoff)
{
   if (fUnderlineLinks != onoff) {
      fUnderlineLinks = onoff;
//    fFlags |= RESIZE_ELEMENTS | RELAYOUT;
//    AddStyle(fPFirst);

      TGHtmlElement *p;
      SHtmlStyle_t style = GetCurrentStyle();
      for (p = fPFirst; p; p = p->fPNext) {
         if (p->fType == Html_A) {
            if (fAnchorStart) {
               style = PopStyleStack(Html_EndA);
               fAnchorStart = 0;
               fAnchorFlags = 0;
            }
            const char *z = p->MarkupArg("href", 0);
            if (z) {
               style.fColor = GetLinkColor(z);
               if (fUnderlineLinks) style.fFlags |= STY_Underline;
               fAnchorFlags |= STY_Anchor;
               PushStyleStack(Html_EndA, style);
               fAnchorStart = (TGHtmlAnchor *) p;
            }
         } else if (p->fType == Html_EndA) {
            if (fAnchorStart) {
               ((TGHtmlRef *)p)->fPOther = fAnchorStart;
               style = PopStyleStack(Html_EndA);
               fAnchorStart = 0;
               fAnchorFlags = 0;
            }
         }
         p->fStyle.fFlags &= ~STY_Underline;
         p->fStyle.fFlags |= (style.fFlags & STY_Underline);
      }

      RedrawEverything();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Sets base URI.

void TGHtml::SetBaseUri(const char *uri)
{
   if (fZBase) delete[] fZBase;
   fZBase = 0;
   if (uri) fZBase = StrDup(uri);
}

////////////////////////////////////////////////////////////////////////////////
/// Go to anchor position.

int TGHtml::GotoAnchor(const char *name)
{
   const char *z;
   TGHtmlElement *p;

   for (p = fPFirst; p; p = p->fPNext) {
      if (p->fType == Html_A) {
         z = p->MarkupArg("name", 0);
         if (z && strcmp(z, name) == 0) {
            ScrollToPosition(TGLongPosition(fVisible.fX, ((TGHtmlAnchor *)p)->fY));
            return kTRUE;
         }
      }
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Given a string, this procedure returns a unique identifier for the
/// string.
///
/// This procedure returns a pointer to a new char string corresponding to
/// the "string" argument. The new string has a value identical to string
/// (strcmp will return 0), but it's guaranteed that any other calls to this
/// procedure with a string equal to "string" will return exactly the same
/// result (i.e. can compare pointer *values* directly, without having to
/// call strcmp on what they point to).

const char *TGHtml::GetUid(const char *string)
{
   //int dummy;

   TObjString *obj = 0;
   obj = (TObjString*)fUidTable->FindObject(string);

   if (!obj) {
      obj = new TObjString(string);
      fUidTable->Add(obj);
   }

   return (const char *)obj->GetName();
}

////////////////////////////////////////////////////////////////////////////////
/// Computes virtual size of html area.

void TGHtml::ComputeVirtualSize()
{
   fVirtualSize = TGDimension(fMaxX, fMaxY);
}

////////////////////////////////////////////////////////////////////////////////
/// Clear the cache of GCs

void TGHtml::ClearGcCache()
{
   int i;

   for (i = 0; i < N_CACHE_GC; i++) {
      if (fAGcCache[i].fIndex) {
         gVirtualX->DeleteGC(fAGcCache[i].fGc);
         fAGcCache[i].fIndex = 0;
      }
   }
   fGcNextToFree = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the main layout context in the main widget.  This happens
/// before we redo the layout, or just before deleting the widget.

void TGHtml::ResetLayoutContext()
{
   fLayoutContext.Reset();
}

////////////////////////////////////////////////////////////////////////////////
/// This routine is invoked in order to redraw all or part of the HTML
/// widget. This might happen because the display has changed, or in
/// response to an expose event. In all cases, though, this routine
/// is called by an idle handler.

void TGHtml::Redraw()
{
   Pixmap_t pixmap;           // The buffer on which to render HTML
   int x, y, w, h;          // Virtual canvas coordinates of area to draw
   int hw;                  // highlight thickness
   int clipwinH, clipwinW;  // Width and height of the clipping window
   TGHtmlBlock *pBlock;      // For looping over blocks to be drawn
   int redoSelection = 0;   // kTRUE to recompute the selection

   // Don't do anything if we are in the middle of a parse.

   if (fInParse) {
      fFlags &= ~REDRAW_PENDING;
      return;
   }

   // Recompute the layout, if necessary or requested.
   //
   // Calling LayoutDoc() is tricky because LayoutDoc() may invoke one
   // or more user-overriden methods, and these methods could, in theory,
   // do nasty things. So we have to take precautions:
   //
   // *  Do not remove the REDRAW_PENDING flag until after LayoutDoc()
   //    has been called, to prevent a recursive call to Redraw().

   if ((fFlags & RESIZE_ELEMENTS) != 0 && (fFlags & STYLER_RUNNING) == 0) {
      TGHtmlImage *pImage;
      for (pImage = fImageList; pImage; pImage = pImage->fPNext) {
         pImage->fPList = 0;
      }
      fLastSized = 0;
      fFlags &= ~RESIZE_ELEMENTS;
      fFlags |= RELAYOUT;
   }

   // We used to make a distinction between RELAYOUT and EXTEND_LAYOUT.
   // RELAYOUT would be used when the widget was resized, but the
   // less compute-intensive EXTEND_LAYOUT would be used when new
   // text was appended.
   //
   // Unfortunately, EXTEND_LAYOUT has some problem that arise when
   // tables are used.  The quick fix is to make an EXTEND_LAYOUT do
   // a complete RELAYOUT.  Someday, we need to fix EXTEND_LAYOUT so
   // that it works right...

   if ((fFlags & (RELAYOUT | EXTEND_LAYOUT)) != 0
      && (fFlags & STYLER_RUNNING) == 0) {
      fNextPlaced = 0;
      //fNInput = 0;
      fVarId = 0;
      fMaxX = 0;
      fMaxY = 0;
      ResetLayoutContext();
      fFirstBlock = 0;
      fLastBlock = 0;
      redoSelection = 1;
      fFlags &= ~RELAYOUT;
      fFlags |= HSCROLL | VSCROLL | REDRAW_TEXT | EXTEND_LAYOUT;
   }

   if ((fFlags & EXTEND_LAYOUT) && fPFirst != 0) {
      LayoutDoc();
      fFlags &= ~EXTEND_LAYOUT;
      FormBlocks();
      MapControls();
      if (redoSelection && fSelBegin.fP && fSelEnd.fP) {
         UpdateSelection(1);
         UpdateInsert();
      }
   }
   fFlags &= ~REDRAW_PENDING;

   // No need to do any actual drawing if we aren't mapped

////  if (!IsMapped()) return;

   // Update the scrollbars.

   if ((fFlags & (HSCROLL | VSCROLL)) != 0) {
      ComputeVirtualSize();
      fFlags &= ~(HSCROLL | VSCROLL);

      if (fFlags & REDRAW_PENDING) return;
   }

   // Redraw the focus highlight, if requested
   hw = fHighlightWidth;
   if (fFlags & REDRAW_FOCUS) {
      if (hw > 0) {
#if 0
      unsigned long color;

      if (fFlags & GOT_FOCUS) {
         color = highlightColorPtr;
      } else {
         color = highlightBgColorPtr;
      }
      _DrawFocusHighlight(color);
#endif
      }
      fFlags &= ~REDRAW_FOCUS;
   }

   // If the styler is in a callback, abort further processing.
   // TODO: check this!

   if (fFlags & STYLER_RUNNING) {
      goto earlyOut;
   }

   MapControls();

   // Compute the virtual canvas coordinates corresponding to the
   // dirty region of the clipping window.

   clipwinW = fCanvas->GetWidth();
   clipwinH = fCanvas->GetHeight();
   if (fFlags & REDRAW_TEXT) {
      w = clipwinW;
      h = clipwinH;
      x = fVisible.fX;
      y = fVisible.fY;
      fDirtyLeft = 0;
      fDirtyTop = 0;
      fFlags &= ~REDRAW_TEXT;
   } else {
      if (fDirtyLeft < 0) fDirtyLeft = 0;
      if (fDirtyRight > clipwinW) fDirtyRight = clipwinW;
      if (fDirtyTop < 0) fDirtyTop = 0;
      if (fDirtyBottom > clipwinH) fDirtyBottom = clipwinH;
      w = fDirtyRight - fDirtyLeft;
      h = fDirtyBottom - fDirtyTop;
      x = fVisible.fX + fDirtyLeft;
      y = fVisible.fY + fDirtyTop;
   }

   // Skip the rest of the drawing process if the area to be refreshed is
   // less than zero
   if (w > 0 && h > 0) {
      GContext_t gcBg;
      TGRectangle xrec;
      // printf("Redraw %dx%d at %d,%d\n", w, h, x, y);

      // Allocate and clear a pixmap upon which to draw
      gcBg = GetGC(COLOR_Background, FONT_Any);
      pixmap = gVirtualX->CreatePixmap(fCanvas->GetId(), w, h);
      xrec.fX = 0;
      xrec.fY = 0;
      xrec.fW = w;
      xrec.fH = h;
#if 0

//old--    XFillRectangles(GetDisplay(), pixmap, gcBg, &xrec, 1);
//new--    if (fBgImage)
//           BGDraw(fVisible.fX, fVisible.fY, w, h, fBgImage);
#else

      fWhiteGC.SetTileStipXOrigin(-fVisible.fX - fDirtyLeft);
      fWhiteGC.SetTileStipYOrigin(-fVisible.fY - fDirtyTop);

      gVirtualX->FillRectangle(pixmap, fWhiteGC.GetGC(), 0, 0, w, h);
      UpdateBackgroundStart();  // back to original
#endif

      // Render all visible HTML onto the pixmap
      for (pBlock = fFirstBlock; pBlock; pBlock = pBlock->fBNext) {
         if (pBlock->fTop <= y+h && pBlock->fBottom >= y-10 &&
            pBlock->fLeft <= x+w && pBlock->fRight >= x-10) {
            BlockDraw(pBlock, pixmap, x, y, w, h, pixmap);
         }
      }

      // Finally, copy the pixmap onto the window and delete the pixmap
      gVirtualX->CopyArea(pixmap, fCanvas->GetId(),
                          gcBg, 0, 0, w, h, fDirtyLeft, fDirtyTop);
      gVirtualX->Update(kFALSE);

      gVirtualX->DeletePixmap(pixmap);
//    XFlush(GetDisplay());
   }

   // Redraw images, if requested
   if (fFlags & REDRAW_IMAGES) {
      TGHtmlImage *pImage;
      TGHtmlImageMarkup *pElem;
      int top, bottom, left, right;     // Coordinates of the clipping window
      int imageTop;                     // Top edge of image

      top = fVisible.fY;
      bottom = top + fCanvas->GetHeight();
      left = fVisible.fX;
      right = left + fCanvas->GetWidth();
      for (pImage = fImageList; pImage; pImage = pImage->fPNext) {
         for (pElem = pImage->fPList; pElem; pElem = pElem->fINext) {
            if (pElem->fRedrawNeeded == 0) continue;
            imageTop = pElem->fY - pElem->fAscent;
            if (imageTop > bottom || imageTop + pElem->fH < top
               || pElem->fX > right || pElem->fX + pElem->fW < left) continue;

            DrawImage(pElem, fCanvas->GetId(), left, top, right, bottom);
         }
      }
      fFlags &= ~(REDRAW_IMAGES | ANIMATE_IMAGES);
   }

   // Set the dirty region to the empty set.
earlyOut:
   fDirtyTop = LARGE_NUMBER;
   fDirtyLeft = LARGE_NUMBER;
   fDirtyBottom = 0;
   fDirtyRight = 0;

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure that a call to the Redraw() routine has been queued.

void TGHtml::ScheduleRedraw()
{
   if ((fFlags & REDRAW_PENDING) == 0 /*&& IsMapped()*/) {
      if (!fIdle) fIdle = new TGIdleHandler(this);
      fFlags |= REDRAW_PENDING;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handles idle event.

Bool_t TGHtml::HandleIdleEvent(TGIdleHandler *idle)
{
   if (idle != fIdle) return kFALSE;
   Redraw();
   delete fIdle;
   fIdle = NULL;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// If any part of the screen needs to be redrawn, then call this routine
/// with the values of a box (in window coordinates) that needs to be
/// redrawn. This routine will schedule an idle handler to do the redraw.
///
/// The box coordinates are relative to the clipping window (fCanvas).

void TGHtml::RedrawArea(int left, int top, int right, int bottom)
{
   if (bottom < 0) return;
   if (top > (int)fCanvas->GetHeight()) return;
   if (right < 0) return;
   if (left > (int)fCanvas->GetWidth()) return;
   if (fDirtyTop > top) fDirtyTop = top;
   if (fDirtyLeft > left) fDirtyLeft = left;
   if (fDirtyBottom < bottom) fDirtyBottom = bottom;
   if (fDirtyRight < right) fDirtyRight = right;
   ScheduleRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Draw region defined by [x,y] [w,h].

void TGHtml::DrawRegion(Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   TGView::DrawRegion(x, y, w, h);

#if 0
   RedrawArea(x, y, x + w + 1, y + h + 1);
#else
   int left = x;
   int top = y;
   int right = x + w + 1;
   int bottom = y + h + 1;
   if (bottom < 0) return;
   if (top > (int) fCanvas->GetHeight()) return;
   if (right < 0) return;
   if (left > (int)fCanvas->GetWidth()) return;
   if (fDirtyTop > top) fDirtyTop = top;
   if (fDirtyLeft > left) fDirtyLeft = left;
   if (fDirtyBottom < bottom) fDirtyBottom = bottom;
   if (fDirtyRight < right) fDirtyRight = right;

   fFlags |= REDRAW_PENDING;
   Redraw();
#endif
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Layout html widget.

Bool_t TGHtml::ItemLayout()
{
#if 0
   fFlags |= RELAYOUT | VSCROLL | HSCROLL;
   Redraw(); //RedrawEverything();
#else
   fNextPlaced = 0;
   //fNInput = 0;
   fVarId = 0;
   fMaxX = 0;
   fMaxY = 0;
   ResetLayoutContext();
   fFirstBlock = 0;
   fLastBlock = 0;
   if (fPFirst != 0) {
      LayoutDoc();
      FormBlocks();
      MapControls();
      if (fSelBegin.fP && fSelEnd.fP) {
         UpdateSelection(1);
         UpdateInsert();
      }
   }
   ComputeVirtualSize();
   ScheduleRedraw();
#endif
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw the TGHtmlBlock given.

void TGHtml::RedrawBlock(TGHtmlBlock *p)
{
   if (p) {
      RedrawArea(p->fLeft - fVisible.fX, p->fTop - fVisible.fY,
                 p->fRight - fVisible.fX + 1, p->fBottom - fVisible.fY);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Call this routine to force the entire widget to be redrawn.

void TGHtml::RedrawEverything()
{
   fFlags |= REDRAW_FOCUS | REDRAW_TEXT;
   ScheduleRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Call this routine to cause all of the rendered HTML at the
/// virtual canvas coordinate of Y and beyond to be redrawn.

void TGHtml::RedrawText(int y)
{
   int clipHeight;     // Height of the clipping window

   clipHeight = fCanvas->GetHeight();
   y -= fVisible.fY;
   if (y < clipHeight) {
      RedrawArea(0, y, LARGE_NUMBER, clipHeight);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Erase all data from the HTML widget. Bring it back to an empty screen.

void TGHtml::HClear()
{
   int i;
   TGHtmlElement *p, *fPNext;

   fXMargin = fYMargin = 0; //HTML_INDENT/4;

   DeleteControls();
   for (p = fPFirst; p; p = fPNext) {
      fPNext = p->fPNext;
      delete p;
   }
   fPFirst = 0;
   fPLast = 0;
   fNToken = 0;
   if (fZText) delete[] fZText;
   fZText = 0;
   fNText = 0;
   fNAlloc = 0;
   fNComplete = 0;
   fIPlaintext = 0;

   for (i = 0; i < N_COLOR; ++i) {
      if (fApColor[i] != 0) FreeColor(fApColor[i]);
      fApColor[i] = 0;
      fIDark[i] = 0;
      fILight[i] = 0;
   }

   if (!fExiting) {
      fFgColor = AllocColor("black");
      fBgColor = AllocColor("white"); //AllocColor("#c0c0c0");
      fNewLinkColor = AllocColor(DEF_HTML_UNVISITED);
      fOldLinkColor = AllocColor(DEF_HTML_VISITED);
      fSelectionColor = AllocColor(DEF_HTML_SELECTION_COLOR);

      fApColor[COLOR_Normal] = fFgColor;
      fApColor[COLOR_Visited] = fOldLinkColor;
      fApColor[COLOR_Unvisited] = fNewLinkColor;
      fApColor[COLOR_Selection] = fSelectionColor;
      fApColor[COLOR_Background] = fBgColor;

      SetBackgroundColor(fApColor[COLOR_Background]->fPixel);
      SetBackgroundPixmap(0);  // use solid color
   }

   fColorUsed = 0;
   while (fImageList) {
      TGHtmlImage *p2 = fImageList;
      fImageList = p2->fPNext;
      delete p2;
   }

   if (fBgImage) delete fBgImage;
   fBgImage = 0;

   while (fStyleStack) {
      SHtmlStyleStack_t *p2 = fStyleStack;
      fStyleStack = p2->fPNext;
      delete p2;
   }
   ClearGcCache();
   ResetLayoutContext();
//  if (fZBase) delete[] fZBase;
//  fZBase = 0;

   if (fZBaseHref) delete [] fZBaseHref;
   fZBaseHref = 0;
   fLastSized = 0;
   fNextPlaced = 0;
   fFirstBlock = 0;
   fLastBlock = 0;
   fNInput = 0;
   fNForm = 0;
   fVarId = 0;
   fParaAlignment = ALIGN_None;
   fRowAlignment = ALIGN_None;
   fAnchorFlags = 0;
   fInDt = 0;
   fAnchorStart = 0;
   fFormStart = 0;
   fInnerList = 0;
   fMaxX = 0;
   fMaxY = 0;
#if 0  // in OXView::Clear()
   fVisible = TGPosition(0, 0);
   _virtualSize = TGDimension(0, 0);
   ScrollTTGPosition(fVisible);
#endif
   fPInsBlock = 0;
   fIns.fP = 0;
   fSelBegin.fP = 0;
   fSelEnd.fP = 0;
   fPSelStartBlock = 0;
   fPSelEndBlock = 0;
   fHasScript = 0;
   fHasFrames = 0;
   fLastUri = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle timer event.

Bool_t TGHtml::HandleTimer(TTimer *t)
{
   if (t == fInsTimer) {
      if (fInsTimer) delete fInsTimer;
      fInsTimer = NULL;
      FlashCursor();
      return kTRUE;
   } else {
      TGHtmlImage *pImage;
      for (pImage = fImageList; pImage; pImage = pImage->fPNext) {
         if (pImage->fTimer == t) {
            AnimateImage(pImage);
            return kTRUE;
         }
      }
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Flash the insertion cursor.

void TGHtml::FlashCursor()
{
   if (fPInsBlock == 0 || fInsOnTime <= 0 || fInsOffTime <= 0) return;
   RedrawBlock(fPInsBlock);
   if ((fFlags & GOT_FOCUS) == 0) {
      fInsStatus = 0;
   } else if (fInsStatus) {
      fInsTimer = new TTimer(this, fInsOffTime);
      fInsStatus = 0;
   } else {
      fInsTimer = new TTimer(this, fInsOnTime);
      fInsStatus = 1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return a GC from the cache.  As many as N_CACHE_GCs are kept valid
/// at any one time.  They are replaced using an LRU algorithm.
///
/// A value of FONT_Any (-1) for the font means "don't care".

GContext_t TGHtml::GetGC(int color, int font)
{
   int i, j;
   GcCache_t *p = fAGcCache;
   GCValues_t gcValues;
   TGFont *xfont;

   // Check for an existing GC.

   if (color < 0 || color >= N_COLOR) color = 0;
   if (font < FONT_Any || font >= N_FONT) font = FONT_Default;

   for (i = 0; i < N_CACHE_GC; i++, p++) {
      if (p->fIndex == 0) continue;
      if ((font < 0 || p->fFont == font) && p->fColor == color) {
         if (p->fIndex > 1) {
            for (j = 0; j < N_CACHE_GC; j++) {
               if (fAGcCache[j].fIndex && fAGcCache[j].fIndex < p->fIndex ) {
                  fAGcCache[j].fIndex++;
               }
            }
            p->fIndex = 1;
         }
         return fAGcCache[i].fGc;
      }
   }

   // No GC matches. Find a place to allocate a new GC.

   p = fAGcCache;
   for (i = 0; i < N_CACHE_GC; i++, p++) {
      if (p->fIndex == 0 || p->fIndex == N_CACHE_GC) break;
   }
   if (i >= N_CACHE_GC) {  // No slot, so free one (round-robin)
      p = fAGcCache;
      for (i = 0; i < N_CACHE_GC && i < fGcNextToFree; ++i, ++p) {}
         fGcNextToFree = (fGcNextToFree + 1) % N_CACHE_GC;
         gVirtualX->DeleteGC(p->fGc);
   }
   gcValues.fForeground = fApColor[color]->fPixel;
   gcValues.fGraphicsExposures = kTRUE;
   gcValues.fMask = kGCForeground | kGCGraphicsExposures;

   if (font < 0) font = FONT_Default;
   xfont = GetFont(font);

   if (xfont) {
      gcValues.fFont = xfont->GetFontHandle();
      gcValues.fMask |= kGCFont;
   }

   p->fGc = gVirtualX->CreateGC(fId, &gcValues);

   if (p->fIndex == 0) p->fIndex = N_CACHE_GC + 1;
   for (j = 0; j < N_CACHE_GC; j++) {
      if (fAGcCache[j].fIndex && fAGcCache[j].fIndex < p->fIndex) {
         fAGcCache[j].fIndex++;
      }
   }
   p->fIndex = 1;
   p->fFont = font;
   p->fColor = color;

   return p->fGc;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve any valid GC. The font and color don't matter since the
/// GC will only be used for copying.

GContext_t TGHtml::GetAnyGC()
{
   int i;
   GcCache_t *p = fAGcCache;

   for (i = 0; i < N_CACHE_GC; i++, p++) {
      if (p->fIndex) return p->fGc;
   }

   return GetGC(COLOR_Normal, FONT_Default);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle focus change event.

Bool_t TGHtml::HandleFocusChange(Event_t *event)
{
   if (event->fType == kFocusIn) {
      fFlags |= GOT_FOCUS | REDRAW_FOCUS;
      ScheduleRedraw();
      UpdateInsert();
   } else {  // FocusOut
      fFlags &= ~GOT_FOCUS;
      fFlags |= REDRAW_FOCUS;
      ScheduleRedraw();
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// This routine searchs for a hyperlink beneath the coordinates x,y
/// and returns a pointer to the HREF for that hyperlink. The text
/// is held in one of the markup argv[] fields of the <a> markup.

TGHtmlInput *TGHtml::GetInputElement(int x, int y)
{
   TGHtmlInput *p;     // For looping over all controls
   int vx, vy, vw, vh;    // Part of the virtual canvas that is visible

   vx = fVisible.fX;
   vy = fVisible.fY;
   vw = fCanvas->GetWidth();
   vh = fCanvas->GetHeight();
   for (p = fFirstInput; p; p = p->fINext) {
      if (p->fFrame == 0) continue;
      if (p->fY < vy + vh && p->fY + p->fH > vy &&
          p->fX < vx + vw && p->fX + p->fW > vx) {
         if ((x > p->fX) && (y > p->fY) && (x < (p->fX + p->fW)) &&
             (y < (p->fY + p->fH)) ) {
            return p;
         }
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle html input (button, checkbox, ...) event.

Bool_t TGHtml::HandleHtmlInput(TGHtmlInput *pr, Event_t *event)
{
   Window_t childdum;
   Event_t eventSt;
   eventSt.fType      = event->fType;
   eventSt.fWindow    = event->fWindow;
   eventSt.fTime      = event->fTime;
   eventSt.fX         = 2;
   eventSt.fY         = 2;
   eventSt.fXRoot     = event->fXRoot;
   eventSt.fYRoot     = event->fYRoot;
   eventSt.fCode      = event->fCode;
   eventSt.fState     = event->fState;
   eventSt.fWidth     = event->fWidth;
   eventSt.fHeight    = event->fHeight;
   eventSt.fCount     = event->fCount;
   eventSt.fSendEvent = event->fSendEvent;
   eventSt.fHandle    = event->fHandle;
   eventSt.fFormat    = event->fFormat;
   eventSt.fUser[0]   = event->fUser[0];
   eventSt.fUser[1]   = event->fUser[1];
   eventSt.fUser[2]   = event->fUser[2];
   eventSt.fUser[3]   = event->fUser[3];
   eventSt.fUser[4]   = event->fUser[4];
   gVirtualX->TranslateCoordinates(GetId(), pr->fFrame->GetId(),
                                   event->fX, event->fY, eventSt.fX,
                                   eventSt.fY, childdum);

   const char *name = pr->MarkupArg("name", 0);
   const char *val = pr->MarkupArg("value", 0);
   switch (pr->fItype) {
      case INPUT_TYPE_Submit:
      case INPUT_TYPE_Button: {
         TGButton *b = (TGButton *) pr->fFrame;
         Bool_t was = !b->IsDown();
         b->HandleButton(&eventSt);
         Bool_t now = !b->IsDown();
         if (!was && now) {
            if (pr->fItype == INPUT_TYPE_Submit)
               SubmitClicked(val);   // emit SubmitClicked
            else
               ButtonClicked(name, val);   // emit ButtonClicked
         }
         break;
      }
      case INPUT_TYPE_Radio: {
         TGRadioButton *rb = (TGRadioButton *) pr->fFrame;
         Bool_t was = !rb->IsDown();
         rb->HandleButton(&eventSt);
         Bool_t now = !rb->IsDown();
         if ((!was && now) || (was && !now)) {
            HandleRadioButton(pr);
            RadioChanged(name, val);      // emit RadioChanged
         }
         break;
      }
      case INPUT_TYPE_Checkbox: {
         TGCheckButton *cb = (TGCheckButton *) pr->fFrame;
         Bool_t was = !cb->IsDown();
         cb->HandleButton(&eventSt);
         Bool_t now = !cb->IsDown();
         if ((!was && now) || (was && !now))
            CheckToggled(name, !now, val);   // emit CheckToggled
         break;
      }
      case INPUT_TYPE_Text:
      case INPUT_TYPE_Password: {
         TGTextEntry *te = (TGTextEntry *) pr->fFrame;
         te->SetFocus();
         break;
      }
      case INPUT_TYPE_Select: {
         RemoveInput(kButtonPressMask | kButtonReleaseMask | kPointerMotionMask);
         eventSt.fUser[0] = childdum;
         if (pr->fFrame->InheritsFrom("TGComboBox"))
            ((TGComboBox *)pr->fFrame)->HandleButton(&eventSt);
         else if (pr->fFrame->InheritsFrom("TGListBox"))
            ((TGListBox *)pr->fFrame)->HandleButton(&eventSt);
         InputSelected(name, val); // emit InputSelected
         AddInput(kButtonPressMask | kButtonReleaseMask | kPointerMotionMask);
         break;
      }
      default:
         break;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle radio button event.

Bool_t TGHtml::HandleRadioButton(TGHtmlInput *p)
{
   TGHtmlInput *pr;
   for (pr = fFirstInput; pr; pr = pr->fINext) {
      if ((pr->fPForm == p->fPForm) && (pr->fItype == INPUT_TYPE_Radio)) {
         if (pr != p) {
            if (strcmp(pr->MarkupArg("name", ""), p->MarkupArg("name", "")) == 0) {
               ((TGRadioButton *)pr->fFrame)->SetState(kButtonUp);
            }
         }
      }
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Emit ButtonClicked() signal.

void TGHtml::ButtonClicked(const char *name, const char *val)
{
   Long_t args[2];

   args[0] = (Long_t)name;
   args[1] = (Long_t)val;

   Emit("ButtonClicked(char*,char*)", args);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit CheckToggled() signal.

void TGHtml::CheckToggled(const char *name, Bool_t on, const char *val)
{
   Long_t args[3];

   args[0] = (Long_t)name;
   args[1] = on;
   args[2] = (Long_t)val;

   Emit("CheckToggled(char*,Bool_t,char*)", args);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit RadioChanged() signal.

void TGHtml::RadioChanged(const char *name, const char *val)
{
   Long_t args[2];

   args[0] = (Long_t)name;
   args[1] = (Long_t)val;

   Emit("RadioChanged(char*,char*)", args);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit Selected() signal.

void TGHtml::InputSelected(const char *name, const char *val)
{
   Long_t args[2];

   args[0] = (Long_t)name;
   args[1] = (Long_t)val;

   Emit("InputSelected(char*,char*)", args);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit SubmitClicked() signal.

void TGHtml::SubmitClicked(const char *val)
{
   Emit("SubmitClicked(char*)", val);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button event.

Bool_t TGHtml::HandleButton(Event_t *event)
{
   int amount, ch;

   ch = fCanvas->GetHeight();
   amount = fScrollVal.fY * TMath::Max(ch/6, 1);

   int ix = event->fX + fVisible.fX;
   int iy = event->fY + fVisible.fY;
   TGHtmlInput *pr = GetInputElement(ix, iy);
   if (pr) {
      HandleHtmlInput(pr, event);
   }
   if ((event->fType == kButtonPress) && (event->fCode == kButton1)) {
      int x = event->fX + fVisible.fX;
      int y = event->fY + fVisible.fY;
      const char *uri = GetHref(x, y);

#if 0  // insertion cursor test
      char ix[20];
      sprintf(ix, "begin");
      SetInsert(ix);
#endif

      if (uri) {
         uri = ResolveUri(uri);
         if (uri) {
            MouseDown(uri);
            //!!delete[] uri;
         }
      }
   } else if (event->fCode == kButton4) {
      ScrollToPosition(TGLongPosition(fVisible.fX, fVisible.fY / fScrollVal.fY - amount));
   } else if (event->fCode == kButton5) {
      ScrollToPosition(TGLongPosition(fVisible.fX, fVisible.fY / fScrollVal.fY + amount));
   } else {
      return TGView::HandleButton(event);
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// handle mouse motion events

Bool_t TGHtml::HandleMotion(Event_t *event)
{
   int x = event->fX + fVisible.fX;
   int y = event->fY + fVisible.fY;
   const char *uri = GetHref(x, y);

   if (uri) {
      gVirtualX->SetCursor(fId, gVirtualX->CreateCursor(kHand));
   } else {
      gVirtualX->SetCursor(fId, gVirtualX->CreateCursor(kPointer));
   }

   if (uri != fLastUri) {
      fLastUri = uri;
      if (uri) uri = ResolveUri(uri);
      MouseOver(uri);
      //!!if (uri) delete [] uri;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// The rendering and layout routines should call this routine in order to
/// get a font structure. The iFont parameter specifies which of the N_FONT
/// fonts should be obtained. The font is allocated if necessary.

TGFont *TGHtml::GetFont(int iFont)
{
   TGFont *toFree = 0;

   if (iFont < 0) iFont = 0;
   if (iFont >= N_FONT) { iFont = N_FONT - 1; CANT_HAPPEN; }

   // If the font has previously been allocated, but the "fontValid" bitmap
   // shows it is no longer valid, then mark it for freeing later. We use
   // a policy of allocate-before-free because xclass' font cache operates
   // much more efficiently that way.

   if (!FontIsValid(iFont) && fAFont[iFont] != 0) {
      toFree = fAFont[iFont];
      fAFont[iFont] = 0;
   }

   // If we need to allocate a font, first construct the font name then
   // allocate it.

   if (fAFont[iFont] == 0) {
      char name[200];         // Name of the font
      const char *familyStr = "";
      int iFamily;
      int iSize;
      int size;

      iFamily = FontFamily(iFont) >> 3;
      iSize = FontSize(iFont) + 1;

      switch (iFamily) {
//#define TIMES
#ifdef TIMES
         case 0:  familyStr = "times -%d";                 break;
         case 1:  familyStr = "times -%d bold";            break;
         case 2:  familyStr = "times -%d italic";          break;
         case 3:  familyStr = "times -%d bold italic";     break;
         case 4:  familyStr = "courier -%d";               break;
         case 5:  familyStr = "courier -%d bold";          break;
         case 6:  familyStr = "courier -%d italic";        break;
         case 7:  familyStr = "courier -%d bold italic";   break;
         default: familyStr = "times -16";                 CANT_HAPPEN;
#else
         case 0:  familyStr = "helvetica -%d";             break;
         case 1:  familyStr = "helvetica -%d bold";        break;
         case 2:  familyStr = "helvetica -%d italic";      break;
         case 3:  familyStr = "helvetica -%d bold italic"; break;
         case 4:  familyStr = "courier -%d";               break;
         case 5:  familyStr = "courier -%d bold";          break;
         case 6:  familyStr = "courier -%d italic";        break;
         case 7:  familyStr = "courier -%d bold italic";   break;
         default: familyStr = "helvetica -14";             CANT_HAPPEN;
#endif
      }
#if 0
      switch (iSize) {
         case 1:  size = 6+finc/*8*/;   break;
         case 2:  size = 10+finc/*10*/;  break;
         case 3:  size = 12+finc/*12*/;  break;
         case 4:  size = 14+finc/*14*/;  break;
         case 5:  size = 20+finc/*16*/;  break;
         case 6:  size = 24+finc/*18*/;  break;
         case 7:  size = 30+finc/*24*/;  break;
         default: size = 14+finc/*14*/;  CANT_HAPPEN;
      }
#else
      switch (iSize) {
         case 1:  size = 8;   break;
         case 2:  size = 10;  break;
         case 3:  size = 12;  break;
         case 4:  size = 14;  break;
         case 5:  size = 16;  break;
         case 6:  size = 18;  break;
         case 7:  size = 24;  break;
         default: size = 14;  CANT_HAPPEN;
      }
#endif
#ifdef TIMES
      if (iFamily < 4) size += 2;
#endif

      snprintf(name, 199, familyStr, size);

      // Get the named font
      fAFont[iFont] = fClient->GetFont(name);

      if (fAFont[iFont] == 0) {
         fprintf(stderr, "TGHtml: could not get font \"%s\", trying fixed\n",
                 name);
         fAFont[iFont] = fClient->GetFont("fixed");
      }
      if (fAFont[iFont] == 0 ){
         fprintf(stderr, "TGHtml: could not get font \"fixed\", trying "
                 "\"helvetica -12\"\n");
         fAFont[iFont] = fClient->GetFont("helvetica -12");
      }
      FontSetValid(iFont);
   }

   // Free the expired font, if any.

   if (toFree) fClient->FreeFont(toFree);

   return fAFont[iFont];
}

////////////////////////////////////////////////////////////////////////////////
/// Only support rect and circles for now

int TGHtml::InArea(TGHtmlMapArea *p, int left, int top, int x, int y)
{
   int *ip = p->fCoords;
   if (!ip) return 0;

   if (p->fMType == HTML_MAP_RECT) {
      return ((left + ip[0]) <= x && (left + ip[2]) >= x &&
              (top  + ip[1]) <= y && (top  + ip[3]) >= y);
   } else if (p->fMType == HTML_MAP_CIRCLE) {
      int dx = left + ip[0] - x;
      int dy = top + ip[1] - y;
      return (dx * dx + dy * dy <= ip[2] * ip[2]);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns html map element.

TGHtmlElement *TGHtml::GetMap(const char *name)
{
   TGHtmlElement *p = fPFirst;
   const char *z, *zb;

   while (p) {
      if (p->fType == Html_MAP) {
         z = p->MarkupArg("name", 0);
         zb = p->MarkupArg("shape", 0);
         if (zb && *zb != 'r') return 0;
         if (z && !strcmp(z, name)) return p;
      }
      p = p->fPNext;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the squared distance between two colors

float TGHtml::ColorDistance(ColorStruct_t *pA, ColorStruct_t *pB)
{
   float x, y, z;

   x = 0.30 * (pA->fRed - pB->fRed);
   y = 0.61 * (pA->fGreen - pB->fGreen);
   z = 0.11 * (pA->fBlue - pB->fBlue);

   return x*x + y*y + z*z;
}

////////////////////////////////////////////////////////////////////////////////
/// This routine returns an index between 0 and N_COLOR-1 which indicates
/// which ColorStruct_t structure in the fApColor[] array should be used to describe
/// the color specified by the given name.

int TGHtml::GetColorByName(const char *zColor)
{
   ColorStruct_t *pNew;
   int iColor;
   const char *name;  // unique!
   int i, n;
   char zAltColor[16];

   // Netscape accepts color names that are just HEX values, without
   // the # up front.  This isn't valid HTML, but we support it for
   // compatibility.

   n = strlen(zColor);
   if (n == 6 || n == 3 || n == 9 || n == 12) {
      for (i = 0; i < n; i++) {
         if (!isxdigit(zColor[i])) break;
      }
      if (i == n) {
         snprintf(zAltColor, 15, "#%s", zColor);
      } else {
         strlcpy(zAltColor, zColor, sizeof(zAltColor));
      }
      name = GetUid(zAltColor);
   } else {
      name = GetUid(zColor);
   }

   pNew = AllocColor(name);
   if (pNew == 0) {
      return 0;      // Color 0 is always the default
   }

   iColor = GetColorByValue(pNew);
   FreeColor(pNew);

   return iColor;
}


// Macros used in the computation of appropriate shadow colors.

#define MAX_COLOR    65535
#define MAX(A,B)     ((A)<(B)?(B):(A))
#define MIN(A,B)     ((A)<(B)?(A):(B))

////////////////////////////////////////////////////////////////////////////////
/// Check to see if the given color is too dark to be easily distinguished
/// from black.

int TGHtml::IsDarkColor(ColorStruct_t *p)
{
   float x, y, z;

   x = 0.50 * p->fRed;
   y = 1.00 * p->fGreen;
   z = 0.28 * p->fBlue;
   return (x*x + y*y + z*z) < (0.05 * MAX_COLOR * MAX_COLOR);
}

////////////////////////////////////////////////////////////////////////////////
/// Given that the background color is iBgColor, figure out an
/// appropriate color for the dark part of a 3D shadow.

int TGHtml::GetDarkShadowColor(int iBgColor)
{
   if (fIDark[iBgColor] == 0) {
      ColorStruct_t *pRef, val;
      val.fMask = kDoRed | kDoGreen | kDoBlue;
      val.fPixel = 0;
      pRef = fApColor[iBgColor];
      if (IsDarkColor(pRef)) {
         int t1, t2;
         t1 = (int) MIN(MAX_COLOR, pRef->fRed * 1.2);
         t2 = (pRef->fRed * 3 + MAX_COLOR) / 4;
         val.fRed = MAX(t1, t2);
         t1 = (int) MIN(MAX_COLOR, pRef->fGreen * 1.2);
         t2 = (pRef->fGreen * 3 + MAX_COLOR) / 4;
         val.fGreen = MAX(t1, t2);
         t1 = (int) MIN(MAX_COLOR, pRef->fBlue * 1.2);
         t2 = (pRef->fBlue * 3 + MAX_COLOR) / 4;
         val.fBlue = MAX(t1, t2);
      } else {
         val.fRed = (unsigned short) (pRef->fRed * 0.6);
         val.fGreen = (unsigned short) (pRef->fGreen * 0.6);
         val.fBlue = (unsigned short) (pRef->fBlue * 0.6);
      }
      fIDark[iBgColor] = GetColorByValue(&val) + 1;
   }

   return fIDark[iBgColor] - 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Check to see if the given color is too light to be easily distinguished
/// from white.

int TGHtml::IsLightColor(ColorStruct_t *p)
{
   return p->fGreen >= 0.85 * MAX_COLOR;
}

////////////////////////////////////////////////////////////////////////////////
/// Given that the background color is iBgColor, figure out an
/// appropriate color for the bright part of the 3D shadow.

int TGHtml::GetLightShadowColor(int iBgColor)
{
   if (fILight[iBgColor] == 0) {
      ColorStruct_t *pRef, val;
      val.fMask = kDoRed | kDoGreen | kDoBlue;
      val.fPixel = 0;
      pRef = fApColor[iBgColor];
      if (IsLightColor(pRef)) {
         val.fRed = (unsigned short) (pRef->fRed * 0.9);
         val.fGreen = (unsigned short) (pRef->fGreen * 0.9);
         val.fBlue = (unsigned short) (pRef->fBlue * 0.9);
      } else {
         int t1, t2;
         t1 = (int) MIN(MAX_COLOR, pRef->fGreen * 1.4);
         t2 = (pRef->fGreen + MAX_COLOR) / 2;
         val.fGreen = MAX(t1, t2);
         t1 = (int) MIN(MAX_COLOR, pRef->fRed * 1.4);
         t2 = (pRef->fRed + MAX_COLOR) / 2;
         val.fRed = MAX(t1, t2);
         t1 = (int) MIN(MAX_COLOR, pRef->fBlue * 1.4);
         t2 = (pRef->fBlue + MAX_COLOR) / 2;
         val.fBlue = MAX(t1, t2);
      }
      fILight[iBgColor] = GetColorByValue(&val) + 1;
   }

   return fILight[iBgColor] - 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Find a color integer for the color whose color components
/// are given by pRef.

int TGHtml::GetColorByValue(ColorStruct_t *pRef)
{
   int i;
   float dist;
   float closestDist;
   int closest;
   int r, g, b;
# define COLOR_MASK  0xf800

   // Search for an exact match
   r = pRef->fRed & COLOR_MASK;
   g = pRef->fGreen & COLOR_MASK;
   b = pRef->fBlue & COLOR_MASK;
   for (i = 0; i < N_COLOR; i++) {
      ColorStruct_t *p = fApColor[i];
      if (p &&
         ((p->fRed & COLOR_MASK) == r) &&
         ((p->fGreen & COLOR_MASK) == g) &&
         ((p->fBlue & COLOR_MASK) == b)) {
         fColorUsed |= (1<<i);
         return i;
      }
   }

   // No exact matches. Look for a completely unused slot
   for (i = N_PREDEFINED_COLOR; i < N_COLOR; i++) {
      if (fApColor[i] == 0) {
         fApColor[i] = AllocColorByValue(pRef);
         fColorUsed |= (1<<i);
         return i;
      }
   }

   // No empty slots. Look for a slot that contains a color that
   // isn't currently in use.
   for (i = N_PREDEFINED_COLOR; i < N_COLOR; i++) {
      if (((fColorUsed >> i) & 1) == 0) {
         FreeColor(fApColor[i]);
         fApColor[i] = AllocColorByValue(pRef);
         fColorUsed |= (1<<i);
         return i;
      }
   }

   // Ok, find the existing color that is closest to the color requested
   // and use it.
   closest = 0;
   closestDist = ColorDistance(pRef, fApColor[0]);
   for (i = 1; i < N_COLOR; i++) {
      dist = ColorDistance(pRef, fApColor[i]);
      if (dist < closestDist) {
         closestDist = dist;
         closest = i;
      }
   }

   return closest;
}

////////////////////////////////////////////////////////////////////////////////
/// This routine searchs for a hyperlink beneath the coordinates x,y
/// and returns a pointer to the HREF for that hyperlink. The text
/// is held in one of the markup argv[] fields of the <a> markup.

const char *TGHtml::GetHref(int x, int y, const char **target)
{
   TGHtmlBlock *pBlock;
   TGHtmlElement *pElem;

   for (pBlock = fFirstBlock; pBlock; pBlock = pBlock->fBNext) {
      if (pBlock->fTop > y || pBlock->fBottom < y ||
          pBlock->fLeft > x || pBlock->fRight < x) continue;
      pElem = pBlock->fPNext;
      if (pElem->fType == Html_IMG) {
         TGHtmlImageMarkup *image = (TGHtmlImageMarkup *) pElem;
         if (image->fPMap) {
            pElem = image->fPMap->fPNext;
            while (pElem && pElem->fType != Html_EndMAP) {
               if (pElem->fType == Html_AREA) {
                  if (InArea((TGHtmlMapArea *) pElem, pBlock->fLeft, pBlock->fTop, x, y)) {
                     if (target) *target = pElem->MarkupArg("target", 0);
                     return pElem->MarkupArg("href", 0);
                  }
               }
               pElem = pElem->fPNext;
            }
            continue;
         }
      }
      if ((pElem->fStyle.fFlags & STY_Anchor) == 0) continue;
      switch (pElem->fType) {
         case Html_Text:
         case Html_Space:
         case Html_IMG:
            while (pElem && pElem->fType != Html_A) pElem = pElem->fPPrev;
            if (pElem == 0 || pElem->fType != Html_A) break;
            if (target) *target = pElem->MarkupArg("target", 0);
            return pElem->MarkupArg("href", 0);

            default:
               break;
      }
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return coordinates of item

int TGHtml::ElementCoords(TGHtmlElement *p, int /*i*/, int pct, int *coords)
{
   TGHtmlBlock *pBlock;

   while (p && p->fType != Html_Block) p = p->fPPrev;
   if (!p) return 1;

   pBlock = (TGHtmlBlock *) p;
   if (pct) {
      TGHtmlElement *pEnd = fPLast;
      TGHtmlBlock *pb2;
      while (pEnd && pEnd->fType != Html_Block) pEnd = pEnd->fPPrev;
      pb2 = (TGHtmlBlock *) pEnd;
#define HGCo(dir) (pb2 && pb2->dir) ? pBlock->dir * 100 / pb2->dir : 0
      coords[0] = HGCo(fLeft);
      coords[1] = HGCo(fTop);
      coords[3] = HGCo(fRight);
      coords[4] = HGCo(fBottom);
   } else {
      coords[0] = pBlock->fLeft;
      coords[1] = pBlock->fTop;
      coords[2] = pBlock->fRight;
      coords[3] = pBlock->fBottom;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns html element matching attribute name and value.

TGHtmlElement *TGHtml::AttrElem(const char *name, char *value)
{
   TGHtmlElement *p;
   const char *z;

   for (p = fPFirst; p; p = p->fPNext) {
      if (p->fType != Html_A) continue;
      z = p->MarkupArg(name, 0);
      if (z && (strcmp(z, value) == 0)) return p;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Given the selection end-points in fSelBegin and fSelEnd, recompute
/// pSelBeginBlock and fPSelEndBlock, then call UpdateSelectionDisplay()
/// to update the display.
///
/// This routine should be called whenever the selection changes or
/// whenever the set of TGHtmlBlock structures change.

void TGHtml::UpdateSelection(int forceUpdate)
{
   TGHtmlBlock *pBlock;
   int index;
   int needUpdate = forceUpdate;
   int temp;

   if (fSelEnd.fP == 0) fSelBegin.fP = 0;

   IndexToBlockIndex(fSelBegin, &pBlock, &index);
   if (needUpdate || pBlock != fPSelStartBlock) {
      needUpdate = 1;
      RedrawBlock(fPSelStartBlock);
      fPSelStartBlock = pBlock;
      fSelStartIndex = index;
   } else if (index != fSelStartIndex) {
      RedrawBlock(pBlock);
      fSelStartIndex = index;
   }

   if (fSelBegin.fP == 0) fSelEnd.fP = 0;

   IndexToBlockIndex(fSelEnd, &pBlock, &index);
   if (needUpdate || pBlock != fPSelEndBlock) {
      needUpdate = 1;
      RedrawBlock(fPSelEndBlock);
      fPSelEndBlock = pBlock;
      fSelEndIndex = index;
   } else if (index != fSelEndIndex) {
      RedrawBlock(pBlock);
      fSelEndIndex = index;
   }

   if (fPSelStartBlock && fPSelStartBlock == fPSelEndBlock &&
       fSelStartIndex > fSelEndIndex) {
      temp = fSelStartIndex;
      fSelStartIndex = fSelEndIndex;
      fSelEndIndex = temp;
   }

   if (needUpdate) {
      fFlags |= ANIMATE_IMAGES;
      UpdateSelectionDisplay();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// The fPSelStartBlock and fPSelEndBlock values have been changed.
/// This routine's job is to loop over all TGHtmlBlocks and either
/// set or clear the HTML_Selected bits in the .fFlags field
/// as appropriate.  For every TGHtmlBlock where the bit changes,
/// mark that block for redrawing.

void TGHtml::UpdateSelectionDisplay()
{
   int selected = 0;
   SHtmlIndex_t tempIndex;
   TGHtmlBlock *pTempBlock;
   int temp;
   TGHtmlBlock *p;

   for (p = fFirstBlock; p; p = p->fBNext) {
      if (p == fPSelStartBlock) {
         selected = 1;
         RedrawBlock(p);
      } else if (!selected && p == fPSelEndBlock) {
         selected = 1;
         tempIndex = fSelBegin;
         fSelBegin = fSelEnd;
         fSelEnd = tempIndex;
         pTempBlock = fPSelStartBlock;
         fPSelStartBlock = fPSelEndBlock;
         fPSelEndBlock = pTempBlock;
         temp = fSelStartIndex;
         fSelStartIndex = fSelEndIndex;
         fSelEndIndex = temp;
         RedrawBlock(p);
      }
      if (p->fFlags & HTML_Selected) {
         if (!selected) {
            p->fFlags &= ~HTML_Selected;
            RedrawBlock(p);
         }
      } else {
         if (selected) {
            p->fFlags |= HTML_Selected;
            RedrawBlock(p);
         }
      }
      if (p == fPSelEndBlock) {
         selected = 0;
         RedrawBlock(p);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Clear selection.

void TGHtml::LostSelection()
{
   if (fExportSelection) {
      // clear selection
      fPSelStartBlock = 0;
      fPSelEndBlock = 0;
      fSelBegin.fP = 0;
      fSelEnd.fP = 0;
      UpdateSelectionDisplay();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set selection.

int TGHtml::SelectionSet(const char *startIx, const char *endIx)
{
   SHtmlIndex_t sBegin, sEnd;
   int bi, ei;

   if (GetIndex(startIx, &sBegin.fP, &sBegin.fI)) {
      // malformed start index
      return kFALSE;
   }

   if (GetIndex(endIx, &sEnd.fP, &sEnd.fI)) {
      // malformed end index
      return kFALSE;
   }

   bi = TokenNumber(sBegin.fP);
   ei = TokenNumber(sEnd.fP);

   if (!(sBegin.fP && sEnd.fP)) return kTRUE;

   if (bi < ei || (bi == ei && sBegin.fI <= sEnd.fI)) {
      fSelBegin = sBegin;
      fSelEnd = sEnd;
   } else {
      fSelBegin = sEnd;
      fSelEnd = sBegin;
   }

   UpdateSelection(0);
   if (fExportSelection) {
      // TODO:
      // get selection ownership ... fId, XA_PRIMARY
      // selection lost handler must directly call LostSelection()
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Recompute the position of the insertion cursor based on the
/// position in fIns.

void TGHtml::UpdateInsert()
{
   IndexToBlockIndex(fIns, &fPInsBlock, &fInsIndex);
   RedrawBlock(fPInsBlock);
   if (fInsTimer == 0) {
      fInsStatus = 0;
      FlashCursor();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the position of the insertion cursor.

int TGHtml::SetInsert(const char *insIx)
{
   SHtmlIndex_t i;

   if (!insIx) {
      RedrawBlock(fPInsBlock);
      fInsStatus = 0;
      fPInsBlock = 0;
      fIns.fP = 0;
   } else {
      if (GetIndex(insIx, &i.fP, &i.fI)) {
         // malformed index
         return kFALSE;
      }
      RedrawBlock(fPInsBlock);
      fIns = i;
      UpdateInsert();
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a html widget as a C++ statement(s) on output stream out.

void TGHtml::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   out << "   TGHtml *";
   out << GetName() << " = new TGHtml(" << fParent->GetName()
       << "," << GetWidth() << "," << GetHeight()
       << ");"<< std::endl;
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   if (fCanvas->GetBackground() != TGFrame::GetWhitePixel()) {
      out << "   " << GetName() << "->ChangeBackground(" << fCanvas->GetBackground() << ");" << std::endl;
   }

   TString fn;
   TGText txt(GetText());
   fn.Form("Html%s.htm", GetName()+5);
   txt.Save(fn.Data());
   out << "   " << "FILE *f = fopen(\"" << fn.Data() << "\", \"r\");" << std::endl;
   out << "   " << "if (f) {" << std::endl;
   out << "      " << GetName() << "->Clear();" << std::endl;
   out << "      " << GetName() << "->Layout();" << std::endl;
   out << "      " << GetName() << "->SetBaseUri(\"\");" << std::endl;
   out << "      " << "char *buf = (char *)calloc(4096, sizeof(char));" << std::endl;
   out << "      " << "while (fgets(buf, 4096, f)) {" << std::endl;
   out << "         " << GetName() << "->ParseText(buf);" << std::endl;
   out << "      " << "}" << std::endl;
   out << "      " << "free(buf);" << std::endl;
   out << "      " << "fclose(f);" << std::endl;
   out << "   " << "}" << std::endl;
   out << "   " << GetName() << "->Layout();" << std::endl;
}
