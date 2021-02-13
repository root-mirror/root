// @(#)root/gui:$Id$
// Author: Fons Rademakers   30/6/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
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
// TGView                                                               //
//                                                                      //
// A TGView provides the infrastructure for text viewer and editor      //
// widgets. It provides a canvas (TGViewFrame) and (optionally) a       //
// vertical and horizontal scrollbar and methods for marking and        //
// scrolling.                                                           //
//                                                                      //
// The TGView (and derivatives) will generate the following             //
// event messages:                                                      //
// kC_TEXTVIEW, kTXT_ISMARKED, widget id, [true|false]                  //
// kC_TEXTVIEW, kTXT_DATACHANGE, widget id, 0                           //
// kC_TEXTVIEW, kTXT_CLICK2, widget id, position (y << 16) | x)         //
// kC_TEXTVIEW, kTXT_CLICK3, widget id, position (y << 16) | x)         //
// kC_TEXTVIEW, kTXT_F3, widget id, true                                //
// kC_TEXTVIEW, kTXT_OPEN, widget id, 0                                 //
// kC_TEXTVIEW, kTXT_CLOSE, widget id, 0                                //
// kC_TEXTVIEW, kTXT_SAVE, widget id, 0                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGView.h"
#include "TGScrollBar.h"
#include "TGResourcePool.h"
#include "TMath.h"
#include "KeySymbols.h"
#include "RConfigure.h"
#include "TVirtualX.h"

ClassImp(TGViewFrame);

////////////////////////////////////////////////////////////////////////////////
/// Create a editor frame.

TGViewFrame::TGViewFrame(TGView *v, UInt_t w, UInt_t h, UInt_t options,
                         ULong_t back) :
   TGCompositeFrame(v, w, h, options | kOwnBackground, back)
{
   fView = v;

   SetBackgroundColor(back);

   gVirtualX->GrabButton(fId, kAnyButton, kAnyModifier,
                         kButtonPressMask | kButtonReleaseMask |
                         kButtonMotionMask, kNone, kNone);

   AddInput(kKeyPressMask | kEnterWindowMask | kLeaveWindowMask |
            kFocusChangeMask);

   SetWindowAttributes_t wattr;
   wattr.fMask = kWAWinGravity | kWABitGravity;
   wattr.fBitGravity = 1; // NorthWestGravity
   wattr.fWinGravity = 1;
   gVirtualX->ChangeWindowAttributes(fId, &wattr);

   // guibuiding settings
   fEditDisabled = kEditDisableGrab | kEditDisableKeyEnable | kEditDisableBtnEnable;
}


ClassImp(TGView);

////////////////////////////////////////////////////////////////////////////////
/// Create an editor view, containing an TGEditorFrame and (optionally)
/// a horizontal and vertical scrollbar.

TGView::TGView(const TGWindow *p, UInt_t w, UInt_t h, Int_t id,
               UInt_t xMargin, UInt_t yMargin, UInt_t options,
               UInt_t sboptions, ULong_t back)
       : TGCompositeFrame(p, w, h, options, GetDefaultFrameBackground())
{
   fWidgetId    = id;
   fMsgWindow   = p;
   fWidgetFlags = kWidgetWantFocus;

   fXMargin = xMargin;
   fYMargin = yMargin;
   fScrollVal.fX = 1;
   fScrollVal.fY = 1;
   fExposedRegion.Empty();

   fClipboard = fClient->GetResourcePool()->GetClipboard();

   fCanvas = new TGViewFrame(this, 10, 10, kChildFrame | kOwnBackground, back);
   AddFrame(fCanvas);

   if (!(sboptions & kNoHSB)) {
      fHsb = new TGHScrollBar(this, 10, 10, kChildFrame);
      AddFrame(fHsb);
      fHsb->Associate(this);
   } else {
      fHsb = 0;
   }

   if (!(sboptions & kNoVSB)) {
      fVsb = new TGVScrollBar(this, 10, 10, kChildFrame);
      AddFrame(fVsb);
      fVsb->Associate(this);
   } else {
      fVsb = 0;
   }

   fWhiteGC.SetGraphicsExposures(kTRUE);
   fWhiteGC.SetBackground(back);

   // sets for guibuilding
   if (fVsb) {
      fVsb->SetEditDisabled(kEditDisableGrab | kEditDisableBtnEnable);
   }
   if (fHsb) {
      fHsb->SetEditDisabled(kEditDisableGrab  | kEditDisableBtnEnable);
   }

   fEditDisabled = kEditDisableLayout;

   // layout manager is not used
   delete fLayoutManager;
   fLayoutManager = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete view.

TGView::~TGView()
{
   if (!MustCleanup()) {
      delete fCanvas;
      delete fHsb;
      delete fVsb;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Clear view.

void TGView::Clear(Option_t *)
{
   fScrolling = -1;

   fMousePos.fX = fMousePos.fY = -1;
   fVisible.fX  = fVisible.fY = 0;
   UpdateBackgroundStart();
   fVirtualSize = TGDimension(0, 0);

   gVirtualX->ClearArea(fCanvas->GetId(), 0, 0,
                        fCanvas->GetWidth(), fCanvas->GetHeight());
   Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Scroll view in specified direction to make newTop the visible location.

void TGView::SetVisibleStart(Int_t newTop, Int_t direction)
{
   if (direction == kHorizontal) {
      if (newTop / fScrollVal.fX == fVisible.fX / fScrollVal.fX) {
         return;
      }
      ScrollCanvas(newTop, kHorizontal);
   } else {
      if (newTop / fScrollVal.fY == fVisible.fY / fScrollVal.fY) {
         return;
      }
      ScrollCanvas(newTop, kVertical);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw region.

void TGView::DrawRegion(Int_t, Int_t, UInt_t, UInt_t)
{
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// update a part of view

void TGView::UpdateRegion(Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   x = x < 0 ? 0 : x;
   y = y < 0 ? 0 : y;

   w = x + w > fCanvas->GetWidth() ? fCanvas->GetWidth() - x : w;
   h = y + h > fCanvas->GetHeight() ? fCanvas->GetHeight() - y : h;

   if (fExposedRegion.IsEmpty()) {
      fExposedRegion.fX = x;
      fExposedRegion.fY = y;
      fExposedRegion.fW = w;
      fExposedRegion.fH = h;
   } else {
      TGRectangle r(x, y, w, h);
      fExposedRegion.Merge(r);
   }

   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// set some gc values

void TGView::UpdateBackgroundStart()
{
   fWhiteGC.SetTileStipXOrigin(-fVisible.fX);
   fWhiteGC.SetTileStipYOrigin(-fVisible.fY);
}

////////////////////////////////////////////////////////////////////////////////
/// handle button

Bool_t TGView::HandleButton(Event_t *event)
{
   if (event->fType == kButtonPress) {
      int amount, ch;

      ch = fCanvas->GetHeight();

      if (fScrollVal.fY == 1) {
         amount = fScrollVal.fY * TMath::Max(ch/6, 1);
      } else {
         amount = fScrollVal.fY * 5;
      }

      if (event->fState & kKeyShiftMask) {
         amount = fScrollVal.fY;
      } else if (event->fState & kKeyControlMask) {
         amount = ch - TMath::Max(ch / 20, 1);
      }

      if (event->fCode == kButton4) {
         ScrollDown(amount);
         return kTRUE;
      } else if (event->fCode == kButton5) {
         ScrollUp(amount);
         return kTRUE;
      }
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// redraw

void TGView::DoRedraw()
{
   DrawBorder();

   if (!fExposedRegion.IsEmpty()) {
      DrawRegion(fExposedRegion.fX, fExposedRegion.fY,
                 fExposedRegion.fW, fExposedRegion.fH);
      fExposedRegion.Empty();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle expose events.

Bool_t TGView::HandleExpose(Event_t *event)
{
   if (event->fWindow == fCanvas->GetId()) {

      TGPosition pos(event->fX, event->fY);
      TGDimension dim(event->fWidth, event->fHeight);
      TGRectangle rect(pos, dim);

      if (fExposedRegion.IsEmpty()) {
         fExposedRegion = rect;
      } else {
         if (((!rect.fX && !fExposedRegion.fY) ||
              (!rect.fY && !fExposedRegion.fX)) &&
             ((rect.fX >= (int)fExposedRegion.fW) ||
              (rect.fY >= (int)fExposedRegion.fH))) {
            DrawRegion(rect.fX, rect.fY, rect.fW, rect.fY);
         } else {
            fExposedRegion.Merge(rect);
         }
      }

      fClient->NeedRedraw(this);
   } else {
      return TGCompositeFrame::HandleExpose(event);
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Process scrollbar messages.

Bool_t TGView::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
   switch(GET_MSG(msg)) {
      case kC_HSCROLL:
         switch(GET_SUBMSG(msg)) {
            case kSB_SLIDERTRACK:
            case kSB_SLIDERPOS:
               SetVisibleStart(Int_t(parm1 * fScrollVal.fX), kHorizontal);
               break;
         }
         break;

      case kC_VSCROLL:
         switch(GET_SUBMSG(msg)) {
            case kSB_SLIDERTRACK:
            case kSB_SLIDERPOS:
               SetVisibleStart(Int_t(parm1 * fScrollVal.fY), kVertical);
               break;
         }
         break;

      default:
         break;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// layout view

void TGView::Layout()
{
   Bool_t need_vsb, need_hsb;
   Int_t cw, ch;

   need_vsb = need_hsb = kFALSE;

   // test whether we need scrollbars
   cw = fWidth - (fBorderWidth << 1) - fXMargin - 1;
   ch = fHeight - (fBorderWidth << 1) - fYMargin - 1;

   fCanvas->SetWidth(cw);
   fCanvas->SetHeight(ch);
   ItemLayout();

   if ((Int_t)fVirtualSize.fWidth > cw) {
      if (fHsb) {
         need_hsb = kTRUE;
         if (fVsb) ch -= fVsb->GetDefaultWidth();
         if (ch < 0) ch = 0;
         fCanvas->SetHeight(ch);
         ItemLayout();
      }
   }

   if ((Int_t)fVirtualSize.fHeight > ch) {
      if (fVsb) {
         need_vsb = kTRUE;
         if (fHsb) cw -= fHsb->GetDefaultHeight();
         if (cw < 0) cw = 0;
         fCanvas->SetWidth(cw);
         ItemLayout();
      }
   }

   // re-check again (putting the scrollbar could have changed things)

   if ((Int_t)fVirtualSize.fWidth > cw) {
      if (!need_hsb) {
         need_hsb = kTRUE;
         if (fVsb) ch -= fVsb->GetDefaultWidth();
         if (ch < 0) ch = 0;
         fCanvas->SetHeight(ch);
         ItemLayout();
      }
   }

   if (fHsb) {
      if (need_hsb) {
         fHsb->MoveResize(fBorderWidth + fXMargin, ch + fBorderWidth + fYMargin,
                          cw, fHsb->GetDefaultHeight());
         fHsb->MapRaised();
      } else {
         fHsb->UnmapWindow();
         fHsb->SetPosition(0);
      }
   }

   if (fVsb) {
      if (need_vsb) {
         fVsb->MoveResize(cw + fBorderWidth + fXMargin,  fBorderWidth + fYMargin,
                          fVsb->GetDefaultWidth(), ch);
         fVsb->MapWindow();
      } else {
         fVsb->UnmapWindow();
         fVsb->SetPosition(0);
      }
   }
   fCanvas->MoveResize(fBorderWidth + fXMargin, fBorderWidth + fYMargin, cw, ch);

   if (fHsb) {
      fHsb->SetRange(fVirtualSize.fWidth / fScrollVal.fX, fCanvas->GetWidth() / fScrollVal.fX);
   }

   if (fVsb) {
      fVsb->SetRange(fVirtualSize.fHeight / fScrollVal.fY, fCanvas->GetHeight() / fScrollVal.fY);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the border of the text edit widget.

void TGView::DrawBorder()
{
   switch (fOptions & (kSunkenFrame | kRaisedFrame | kDoubleBorder)) {
      case kSunkenFrame | kDoubleBorder:
         if (gClient->GetStyle() < 2) {
            gVirtualX->DrawLine(fId, GetShadowGC()(), 0, 0, fWidth-2, 0);
            gVirtualX->DrawLine(fId, GetShadowGC()(), 0, 0, 0, fHeight-2);
            gVirtualX->DrawLine(fId, GetBlackGC()(), 1, 1, fWidth-3, 1);
            gVirtualX->DrawLine(fId, GetBlackGC()(), 1, 1, 1, fHeight-3);

            gVirtualX->DrawLine(fId, GetHilightGC()(), 0, fHeight-1, fWidth-1, fHeight-1);
            gVirtualX->DrawLine(fId, GetHilightGC()(), fWidth-1, fHeight-1, fWidth-1, 0);
            gVirtualX->DrawLine(fId, GetBckgndGC()(),  1, fHeight-2, fWidth-2, fHeight-2);
            gVirtualX->DrawLine(fId, GetBckgndGC()(),  fWidth-2, 1, fWidth-2, fHeight-2);
            break;
         }
      default:
         TGFrame::DrawBorder();
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Scroll the canvas to pos.

void TGView::ScrollToPosition(TGLongPosition pos)
{
   if (pos.fX < 0) pos.fX = 0;
   if (pos.fY < 0) pos.fY = 0;
   if (pos.fX != fHsb->GetPosition()) fHsb->SetPosition(pos.fX / fScrollVal.fX);
   if (pos.fY != fVsb->GetPosition()) fVsb->SetPosition(pos.fY / fScrollVal.fY);
}

////////////////////////////////////////////////////////////////////////////////
/// Scroll the canvas to new_top in the kVertical or kHorizontal direction.

void TGView::ScrollCanvas(Int_t new_top, Int_t direction)
{
   Point_t points[4];
   Int_t xsrc, ysrc, xdest, ydest, cpyheight, cpywidth;

   if (new_top < 0) {
      return;
   }

   if (direction == kVertical) {
      if (new_top == fVisible.fY) {
         return;
      }

      points[0].fX = points[3].fX = 0;
      points[1].fX = points[2].fX = fCanvas->GetWidth();
      xsrc = xdest = 0;
      cpywidth = 0;
      if (new_top < fVisible.fY) {
         ysrc = 0;
         ydest = Int_t(fVisible.fY - new_top);
         cpyheight = ydest;
         if (ydest > (Int_t)fCanvas->GetHeight()) {
            ydest = fCanvas->GetHeight();
         }

         points[1].fY = points[0].fY = 0;
         points[3].fY = points[2].fY = ydest; // -1;
      } else {
         ydest = 0;
         ysrc = Int_t(new_top - fVisible.fY);
         cpyheight= ysrc;
         if (ysrc > (Int_t)fCanvas->GetHeight()) {
            ysrc = fCanvas->GetHeight();
         }
         points[1].fY = points[0].fY = fCanvas->GetHeight()-ysrc; // +1;
         points[3].fY = points[2].fY = fCanvas->GetHeight();
      }
      fVisible.fY = new_top;

      if (fVisible.fY < 0) {
         fVisible.fY = 0;
      }
   } else {
      if (new_top == fVisible.fX) {
         return;
      }

      points[0].fY = points[1].fY = 0;
      points[2].fY = points[3].fY = fCanvas->GetHeight();
      ysrc = ydest = 0;
      cpyheight = 0;

      if (new_top < fVisible.fX) {
         xsrc = 0;
         xdest = Int_t(fVisible.fX - new_top);
         cpywidth = xdest;
         if (xdest < 0) {
            xdest = fCanvas->GetWidth();
         }
         points[0].fX = points[3].fX = 0;
         points[1].fX = points[2].fX = xdest ; // -1;
      } else {
         xdest = 0;
         xsrc =  Int_t(new_top - fVisible.fX);
         cpywidth = xsrc;
         if (xsrc > (Int_t)fCanvas->GetWidth()) {
            xsrc = fCanvas->GetWidth();
         }
         points[0].fX = points[3].fX = fCanvas->GetWidth()-xsrc; // +1;
         points[1].fX = points[2].fX = fCanvas->GetWidth();
      }
      fVisible.fX = new_top;
      if (fVisible.fX < 0) {
         fVisible.fX = 0;
      }
   }

   UpdateBackgroundStart();

#ifdef R__HAS_COCOA
   //With QuartzView it's quite tough to copy window's pixels to window.
   //TODO: non-optimal solution.
   DrawRegion(0, 0, GetWidth(), GetHeight());
#else
   // Copy the scrolled region to its new position
   gVirtualX->CopyArea(fCanvas->GetId(), fCanvas->GetId(), fWhiteGC(),
                       xsrc, ysrc, fCanvas->GetWidth()-cpywidth,
                       fCanvas->GetHeight()-cpyheight, xdest, ydest);

   UInt_t xdiff = points[2].fX - points[0].fX;
   UInt_t ydiff = points[2].fY - points[0].fY;

   // under windows we need to redraw larger area (why?)
#ifdef WIN32
   xdiff = xdiff << 1;
   ydiff = ydiff << 1;
#endif

   DrawRegion(points[0].fX, points[0].fY, xdiff, ydiff);
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Change background color of the canvas frame.

void TGView::ChangeBackground(Pixel_t col)
{
   fCanvas->SetBackgroundColor(col);
   fWhiteGC.SetBackground(col);
   fWhiteGC.SetForeground(col);
   DrawRegion(0, 0, fCanvas->GetWidth(), fCanvas->GetHeight());
}

////////////////////////////////////////////////////////////////////////////////
/// Set background color of the canvas frame.

void TGView::SetBackgroundColor(Pixel_t col)
{
   fCanvas->SetBackgroundColor(col);
   fWhiteGC.SetBackground(col);
   fWhiteGC.SetForeground(col);
}

////////////////////////////////////////////////////////////////////////////////
/// Set backgound  pixmap

void TGView::SetBackgroundPixmap(Pixmap_t p)
{
   fCanvas->SetBackgroundPixmap(p);
}
