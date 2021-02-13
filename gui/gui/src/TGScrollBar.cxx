// @(#)root/gui:$Id$
// Author: Fons Rademakers   10/01/98

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
// TGScrollBar and TGScrollBarElement                                   //
//                                                                      //
// The classes in this file implement scrollbars. Scrollbars can be     //
// either placed horizontal or vertical. A scrollbar contains three     //
// TGScrollBarElements: The "head", "tail" and "slider". The head and   //
// tail are fixed at either end and have the typical arrows in them.    //
//                                                                      //
// The TGHScrollBar will generate the following event messages:         //
// kC_HSCROLL, kSB_SLIDERPOS, position, 0                               //
// kC_HSCROLL, kSB_SLIDERTRACK, position, 0                             //
//                                                                      //
// The TGVScrollBar will generate the following event messages:         //
// kC_VSCROLL, kSB_SLIDERPOS, position, 0                               //
// kC_VSCROLL, kSB_SLIDERTRACK, position, 0                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGScrollBar.h"
#include "TGResourcePool.h"
#include "TGPicture.h"
#include "TImage.h"
#include "TSystem.h"
#include "TTimer.h"
#include "TEnv.h"
#include "TVirtualX.h"

#include <iostream>


Pixmap_t TGScrollBar::fgBckgndPixmap = 0;
Int_t    TGScrollBar::fgScrollBarWidth = kDefaultScrollBarWidth;

ClassImp(TGScrollBarElement);
ClassImp(TGScrollBar);
ClassImp(TGHScrollBar);
ClassImp(TGVScrollBar);


////////////////////////////////////////////////////////////////////////////////

class TSBRepeatTimer : public TTimer {
private:
   TGScrollBar   *fScrollBar;   // scroll bar
   Int_t          fSmallInc;    // step
public:
   TSBRepeatTimer(TGScrollBar *s, Long_t ms, Int_t inc) : TTimer(ms, kTRUE)
      { fScrollBar = s;  fSmallInc = inc; }

   Bool_t Notify();
   Int_t  GetSmallInc() const { return fSmallInc; }
};

////////////////////////////////////////////////////////////////////////////////
/// Notify when timer times out and reset the timer.

Bool_t TSBRepeatTimer::Notify()
{
   fScrollBar->HandleTimer(this);
   Reset();
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGScrollBarElement::TGScrollBarElement(const TGWindow *p, const TGPicture *pic,
                                       UInt_t w, UInt_t h, UInt_t options, Pixel_t back) :
                                       TGFrame(p, w, h, options | kOwnBackground, back)
{
   fPic = fPicN = pic;
   fState = kButtonUp;
   fPicD = 0;
   fStyle = 0;
   if ((gClient->GetStyle() > 1) || (p && p->InheritsFrom("TGScrollBar")))
      fStyle = gClient->GetStyle();

   fBgndColor = fBackground;
   fHighColor = gClient->GetResourcePool()->GetHighLightColor();
   AddInput(kEnterWindowMask | kLeaveWindowMask);
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGScrollBarElement::~TGScrollBarElement()
{
   if (fPicD) {
      fClient->FreePicture(fPicD);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Change state of scrollbar element (either up or down).

void TGScrollBarElement::SetState(Int_t state)
{
   if (state != fState) {
      switch (state) {
         case kButtonDown:
            fOptions &= ~kRaisedFrame;
            fOptions |= kSunkenFrame;
            break;
         case kButtonUp:
         case kButtonDisabled:
            fOptions &= ~kSunkenFrame;
            fOptions |= kRaisedFrame;
            break;
      }
      fState = state;
      fClient->NeedRedraw(this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Enable/Disable scroll bar button chaging the state.

void TGScrollBarElement::SetEnabled(Bool_t on)
{
   if (on) {
      if (fState == kButtonUp) {
         return;
      }
      SetState(kButtonUp);
      fPic = fPicN;
   } else {
      if (fState == kButtonDisabled) {
         return;
      }
      SetState(kButtonDisabled);

      if (!fPicD) {
         TImage *img = TImage::Create();
         if (!img) return;
         TImage *img2 = TImage::Create();
         if (!img2) {
            if (img) delete img;
            return;
         }
         TString back = gEnv->GetValue("Gui.BackgroundColor", "#c0c0c0");
         img2->FillRectangle(back.Data(), 0, 0, fPic->GetWidth(), fPic->GetHeight());
         img->SetImage(fPicN->GetPicture(), fPicN->GetMask());
         Pixmap_t mask = img->GetMask();
         img2->Merge(img, "overlay");

         TString name = "disbl_";
         name += fPic->GetName();
         fPicD = fClient->GetPicturePool()->GetPicture(name.Data(), img2->GetPixmap(),
                                                       mask);
         delete img;
         delete img2;
      }
      fPic = fPicD;
   }
   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw border around scollbar element.

void TGScrollBarElement::DrawBorder()
{
   switch (fOptions & (kSunkenFrame | kRaisedFrame)) {
      case kSunkenFrame: // pressed
         gVirtualX->DrawLine(fId, GetBlackGC()(), 0, 0, fWidth-2, 0);
         gVirtualX->DrawLine(fId, GetBlackGC()(), 0, 0, 0, fHeight-2);
         gVirtualX->DrawLine(fId, GetShadowGC()(),  1, 1, fWidth-3, 1);
         gVirtualX->DrawLine(fId, GetShadowGC()(),  1, 1, 1, fHeight-3);

         gVirtualX->DrawLine(fId, GetWhiteGC()(),  0, fHeight-1, fWidth-1, fHeight-1);
         gVirtualX->DrawLine(fId, GetWhiteGC()(),  fWidth-1, fHeight-1, fWidth-1, 1);
         gVirtualX->DrawLine(fId, GetBckgndGC()(),  1, fHeight-2, fWidth-2, fHeight-2);
         gVirtualX->DrawLine(fId, GetBckgndGC()(),  fWidth-2, fHeight-2, fWidth-2, 2);

         if (fPic) {
            int x = (fWidth - fPic->GetWidth()) >> 1;
            int y = (fHeight - fPic->GetHeight()) >> 1;
            fPic->Draw(fId, GetBckgndGC()(), x+1, y+1); // 3, 3
         }
         break;

      case kRaisedFrame: // normal
      case kButtonDisabled:
         if (fStyle > 0) {
            // new modern (flat) version
            if (fBackground == fHighColor) {
               gVirtualX->DrawRectangle(fId, GetShadowGC()(), 0, 0, fWidth-1, fHeight-1);
            }
            else {
               if (fPic == 0)
                  gVirtualX->DrawRectangle(fId, GetShadowGC()(), 0, 0, fWidth-1, fHeight-1);
               else
                  gVirtualX->DrawRectangle(fId, GetBckgndGC()(), 0, 0, fWidth-1, fHeight-1);
            }
            if (fParent && fParent->InheritsFrom("TGHScrollBar")) {
               if (fWidth > 20) {
                  gVirtualX->DrawLine(fId, GetShadowGC()(), (fWidth/2)-3, 4, (fWidth/2)-3, fHeight-5);
                  gVirtualX->DrawLine(fId, GetShadowGC()(), (fWidth/2),   4, (fWidth/2),   fHeight-5);
                  gVirtualX->DrawLine(fId, GetShadowGC()(), (fWidth/2)+3, 4, (fWidth/2)+3, fHeight-5);
               }
            }
            else if (fParent && fParent->InheritsFrom("TGVScrollBar")) {
               if (fHeight > 20) {
                  gVirtualX->DrawLine(fId, GetShadowGC()(), 4, (fHeight/2)-3, fWidth-5, (fHeight/2)-3);
                  gVirtualX->DrawLine(fId, GetShadowGC()(), 4, (fHeight/2)  , fWidth-5, (fHeight/2));
                  gVirtualX->DrawLine(fId, GetShadowGC()(), 4, (fHeight/2)+3, fWidth-5, (fHeight/2)+3);
               }
            }
            else { // not used in a scroll bar (e.g. in a combo box)
               gVirtualX->DrawRectangle(fId, GetShadowGC()(), 0, 0, fWidth-1, fHeight-1);
            }
         }
         else {
            gVirtualX->DrawLine(fId, GetBckgndGC()(), 0, 0, fWidth-2, 0);
            gVirtualX->DrawLine(fId, GetBckgndGC()(), 0, 0, 0, fHeight-2);
            gVirtualX->DrawLine(fId, GetHilightGC()(), 1, 1, fWidth-3, 1);
            gVirtualX->DrawLine(fId, GetHilightGC()(), 1, 1, 1, fHeight-3);

            gVirtualX->DrawLine(fId, GetShadowGC()(),  1, fHeight-2, fWidth-2, fHeight-2);
            gVirtualX->DrawLine(fId, GetShadowGC()(),  fWidth-2, fHeight-2, fWidth-2, 1);
            gVirtualX->DrawLine(fId, GetBlackGC()(), 0, fHeight-1, fWidth-1, fHeight-1);
            gVirtualX->DrawLine(fId, GetBlackGC()(), fWidth-1, fHeight-1, fWidth-1, 0);
         }
         if (fPic) {
            int x = (fWidth - fPic->GetWidth()) >> 1;
            int y = (fHeight - fPic->GetHeight()) >> 1;
            fPic->Draw(fId, GetBckgndGC()(), x, y); // 2, 2
         }
         break;

      default:
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse crossing event.

Bool_t TGScrollBarElement::HandleCrossing(Event_t *event)
{
   if (fStyle > 0) {
      TGScrollBarElement *el = 0;
      TGScrollBar *bar = 0;
      if ((event->fType == kEnterNotify) && (fState != kButtonDisabled)) {
         fBgndColor = fHighColor;
      } else {
         fBgndColor = fBackground;
      }
      if (event->fType == kLeaveNotify) {
         fBgndColor = fBackground;
      }
      gVirtualX->SetWindowBackground(fId, fBgndColor);
      TGFrame::DoRedraw();
      DrawBorder();
      if (fParent && fParent->InheritsFrom("TGScrollBar")) {
         bar = (TGScrollBar *)fParent;
         if ((el = bar->GetHead()) != this) {
            el->ChangeBackground(fBgndColor);
            el->DrawBorder();
         }
         if ((el = bar->GetTail()) != this) {
            el->ChangeBackground(fBgndColor);
            el->DrawBorder();
         }
         if ((el = bar->GetSlider()) != this) {
            el->ChangeBackground(fBgndColor);
            el->DrawBorder();
         }
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGScrollBar::TGScrollBar(const TGWindow *p, UInt_t w, UInt_t h,
                         UInt_t options, Pixel_t back) :
   TGFrame(p, w, h, options | kOwnBackground, back),
   fX0(0), fY0(0), fXp(0), fYp(0), fDragging(kFALSE), fGrabPointer(kTRUE),
   fRange(0), fPsize(0), fPos(0), fSliderSize(0), fSliderRange(0),
   fSmallInc(1), fHead(0), fTail(0), fSlider(0), fHeadPic(0),
   fTailPic(0), fRepeat(0), fSubw()
{
   fAccelerated = kFALSE;

   fBgndColor = fBackground;
   fHighColor = gClient->GetResourcePool()->GetHighLightColor();

   fMsgWindow = p;
   if (gClient->GetStyle() == 0)
      SetBackgroundPixmap(GetBckgndPixmap());
   SetWindowName();
   AddInput(kEnterWindowMask | kLeaveWindowMask);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete a scrollbar (either horizontal or vertical).

TGScrollBar::~TGScrollBar()
{
   delete fHead;
   delete fTail;
   delete fSlider;
   if (fHeadPic) fClient->FreePicture(fHeadPic);
   if (fTailPic) fClient->FreePicture(fTailPic);
   if (fRepeat) { delete fRepeat; fRepeat = 0; }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse crossing event.

Bool_t TGScrollBar::HandleCrossing(Event_t *event)
{
   if (gClient->GetStyle() > 0) {
      if (event->fType == kEnterNotify) {
         fBgndColor = fHighColor;
      } else {
         fBgndColor = fBackground;
      }
      if (event->fType == kLeaveNotify) {
         fBgndColor = fBackground;
      }
      fHead->ChangeBackground(fBgndColor);
      fTail->ChangeBackground(fBgndColor);
      fSlider->ChangeBackground(fBgndColor);
      fHead->DrawBorder();
      fTail->DrawBorder();
      fSlider->DrawBorder();
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle repeat timer for horizontal or vertical scrollbar. Every time
/// timer times out we move slider.

Bool_t TGScrollBar::HandleTimer(TTimer *t)
{
   // shorten time out time to 50 milli seconds
   t->SetTime(50);

   Window_t  dum1, dum2;
   Event_t   ev;

   ev.fCode    = kButton1;
   ev.fType    = kButtonPress;
   ev.fUser[0] = fSubw;

   if (IsAccelerated()) {
      ++fSmallInc;
      if (fSmallInc > 100) fSmallInc = 100;
   }

   gVirtualX->QueryPointer(fId, dum1, dum2, ev.fXRoot, ev.fYRoot, ev.fX, ev.fY,
                      ev.fState);

   HandleButton(&ev);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning scrollbar background pixmap.

Pixmap_t TGScrollBar::GetBckgndPixmap()
{
   static Bool_t init = kFALSE;
   if (!init) {
      fgBckgndPixmap = gClient->GetResourcePool()->GetCheckeredPixmap();
      init = kTRUE;
   }
   return fgBckgndPixmap;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning the scrollbar width.

Int_t TGScrollBar::GetScrollBarWidth()
{
   return fgScrollBarWidth;
}

////////////////////////////////////////////////////////////////////////////////
/// Change background color

void TGScrollBar::ChangeBackground(Pixel_t back)
{
   TGFrame::ChangeBackground(back);
   fHead->ChangeBackground(back);
   fTail->ChangeBackground(back);
   fSlider->ChangeBackground(back);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a horizontal scrollbar widget.

TGHScrollBar::TGHScrollBar(const TGWindow *p, UInt_t w, UInt_t h,
                           UInt_t options, ULong_t back) :
    TGScrollBar(p, w, h, options, back)
{
   fHeadPic = fClient->GetPicture("arrow_left.xpm");
   fTailPic = fClient->GetPicture("arrow_right.xpm");

   if (!fHeadPic || !fTailPic) {
      Error("TGHScrollBar", "arrow_*.xpm not found");
      return;
   }
   fHead   = new TGScrollBarElement(this, fHeadPic, fgScrollBarWidth, fgScrollBarWidth,
                                    kRaisedFrame);
   fTail   = new TGScrollBarElement(this, fTailPic, fgScrollBarWidth, fgScrollBarWidth,
                                    kRaisedFrame);
   fSlider = new TGScrollBarElement(this, 0, fgScrollBarWidth, 50,
                                    kRaisedFrame);

   gVirtualX->GrabButton(fId, kAnyButton, kAnyModifier, kButtonPressMask |
                         kButtonReleaseMask | kPointerMotionMask, kNone, kNone);

   fDragging = kFALSE;
   fX0 = fY0 = (fgScrollBarWidth = TMath::Max(fgScrollBarWidth, 5));
   fPos = 0;

   fRange = TMath::Max((Int_t) w - (fgScrollBarWidth << 1), 1);
   fPsize = fRange >> 1;

   fSliderSize  = 50;
   fSliderRange = 1;

   fHead->SetEditDisabled(kEditDisable | kEditDisableGrab);
   fTail->SetEditDisabled(kEditDisable | kEditDisableGrab);
   fSlider->SetEditDisabled(kEditDisable | kEditDisableGrab);
   fEditDisabled = kEditDisableLayout | kEditDisableHeight | kEditDisableBtnEnable;
}

////////////////////////////////////////////////////////////////////////////////
/// Layout and move horizontal scrollbar components.

void TGHScrollBar::Layout()
{
   // Should also recalculate the slider size and range, etc.
   fHead->Move(0, 0);
   fHead->Resize(fgScrollBarWidth, fgScrollBarWidth);
   fTail->Move(fWidth-fgScrollBarWidth, 0);
   fTail->Resize(fgScrollBarWidth, fgScrollBarWidth);

   if (fSlider->GetX() != fX0) {
      fSlider->Move(fX0, 0);
      fSlider->Resize(50, fgScrollBarWidth);
      fClient->NeedRedraw(fSlider);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle a mouse button event in a horizontal scrolbar.

Bool_t TGHScrollBar::HandleButton(Event_t *event)
{
   Int_t newpos;

   if (event->fCode == kButton4) {
      if (!fHead->IsEnabled()) {
         return kFALSE;
      }
      //scroll left
      newpos = fPos - fPsize;
      if (newpos<0) newpos = 0;
      SetPosition(newpos);
      return kTRUE;
   }
   if (event->fCode == kButton5) {
      if (!fTail->IsEnabled()) {
         return kFALSE;
      }
      // scroll right
      newpos = fPos + fPsize;
      SetPosition(newpos);
      return kTRUE;
   }

   if (event->fType == kButtonPress) {
      if (event->fCode == kButton3) {
         fX0 = event->fX - fSliderSize/2;
         fX0 = TMath::Max(fX0, fgScrollBarWidth);
         fX0 = TMath::Min(fX0, fgScrollBarWidth + fSliderRange);
         ULong_t pos = (ULong_t)(fX0 - fgScrollBarWidth) * (ULong_t)(fRange-fPsize) / (ULong_t)fSliderRange;
         fPos = (Int_t)pos;

         fPos = TMath::Max(fPos, 0);
         fPos = TMath::Min(fPos, fRange-fPsize);
         fSlider->Move(fX0, 0);

         SendMessage(fMsgWindow, MK_MSG(kC_HSCROLL, kSB_SLIDERTRACK), fPos, 0);
         PositionChanged(fPos);
         return kTRUE;
      }

      // fUser[0] contains the subwindow (child) in which the event occured
      // (see GX11Gui.cxx)
      Window_t subw = (Window_t)event->fUser[0];

      if (subw == fSlider->GetId()) {
         fXp = event->fX - fX0;
         fYp = event->fY - fY0;
         fDragging = kTRUE;

      } else {

         if (!fRepeat)
            fRepeat = new TSBRepeatTimer(this, 400, fSmallInc); // 500
         fRepeat->Reset();
         gSystem->AddTimer(fRepeat);
         fSubw = subw;

         if (subw == fHead->GetId()) {
            //if (!fHead->IsEnabled()) {
             //  return kFALSE;
            //}
            fHead->SetState(kButtonDown);
            fPos -= fSmallInc;
         } else if (subw == fTail->GetId()) {
            //if (!fTail->IsEnabled()) {
           //    return kFALSE;
           // }
            fTail->SetState(kButtonDown);
            fPos += fSmallInc;
         } else if (event->fX > fgScrollBarWidth && event->fX < fX0)
            fPos -= fPsize;
         else if (event->fX > fX0+fSliderSize && event->fX < (Int_t)fWidth-fgScrollBarWidth)
            fPos += fPsize;

         fPos = TMath::Max(fPos, 0);
         fPos = TMath::Min(fPos, fRange-fPsize);

         fX0 = fgScrollBarWidth + fPos * fSliderRange / TMath::Max(fRange-fPsize, 1);

         fX0 = TMath::Max(fX0, fgScrollBarWidth);
         fX0 = TMath::Min(fX0, fgScrollBarWidth + fSliderRange);

         fSlider->Move(fX0, 0);

         SendMessage(fMsgWindow, MK_MSG(kC_HSCROLL, kSB_SLIDERTRACK), fPos, 0);
         PositionChanged(fPos);
      }

      // last argument kFALSE forces all specified events to this window
      if (fGrabPointer && !fClient->IsEditable())
         gVirtualX->GrabPointer(fId, kButtonPressMask | kButtonReleaseMask |
                                kPointerMotionMask, kNone, kNone,
                                kTRUE, kFALSE);
   } else {
      fHead->SetState(kButtonUp);
      fTail->SetState(kButtonUp);

      if (fRepeat) {
         fRepeat->Remove();
         fRepeat->SetTime(400);  // might have been shortened in HandleTimer()
         fSmallInc = ((TSBRepeatTimer*)fRepeat)->GetSmallInc();
      }

      fDragging = kFALSE;

      fPos = TMath::Max(fPos, 0);
      fPos = TMath::Min(fPos, fRange-fPsize);

      SendMessage(fMsgWindow, MK_MSG(kC_HSCROLL, kSB_SLIDERPOS), fPos, 0);
      PositionChanged(fPos);

      if (fGrabPointer)
         gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);  // ungrab pointer
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse motion event in a horizontal scrollbar.

Bool_t TGHScrollBar::HandleMotion(Event_t *event)
{
   if (fDragging) {
      fX0 = event->fX - fXp;
      fY0 = event->fY - fYp;

      fX0 = TMath::Max(fX0, fgScrollBarWidth);
      fX0 = TMath::Min(fX0, fgScrollBarWidth + fSliderRange);
      fSlider->Move(fX0, 0);
      ULong_t pos = (ULong_t)(fX0 - fgScrollBarWidth) * (ULong_t)(fRange-fPsize) / (ULong_t)fSliderRange;
      fPos = (Int_t)pos;

      fPos = TMath::Max(fPos, 0);
      fPos = TMath::Min(fPos, fRange-fPsize);

      SendMessage(fMsgWindow, MK_MSG(kC_HSCROLL, kSB_SLIDERTRACK), fPos, 0);
      PositionChanged(fPos);
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set range of horizontal scrollbar.

void TGHScrollBar::SetRange(Int_t range, Int_t page_size)
{
   fRange = TMath::Max(range, 1);
   fPsize = TMath::Max(page_size, 0);
   fPos = TMath::Max(fPos, 0);
   fPos = TMath::Min(fPos, fRange-fPsize);

   fSliderSize = TMath::Max(fPsize * (fWidth - (fgScrollBarWidth << 1)) /
                            fRange, (UInt_t) 6);
   fSliderSize = TMath::Min((UInt_t)fSliderSize, fWidth - (fgScrollBarWidth << 1));

   fSliderRange = TMath::Max(fWidth - (fgScrollBarWidth << 1) - fSliderSize,
                             (UInt_t) 1);

   fX0 = fgScrollBarWidth + fPos * fSliderRange / TMath::Max(fRange-fPsize, 1);
   fX0 = TMath::Max(fX0, fgScrollBarWidth);
   fX0 = TMath::Min(fX0, fgScrollBarWidth + fSliderRange);

   fSlider->Move(fX0, 0);
   fSlider->Resize(fSliderSize, fgScrollBarWidth);
   fClient->NeedRedraw(fSlider);

   //  fPos = (fX0 - fgScrollBarWidth) * (fRange-fPsize) / fSliderRange;

   SendMessage(fMsgWindow, MK_MSG(kC_HSCROLL, kSB_SLIDERPOS), fPos, 0);
   PositionChanged(fPos);
   RangeChanged(fRange);
   PageSizeChanged(fPsize);
}

////////////////////////////////////////////////////////////////////////////////
/// Set logical slider position of horizontal scrollbar.

void TGHScrollBar::SetPosition(Int_t pos)
{
   fPos = TMath::Max(pos, 0);
   fPos = TMath::Min(pos, fRange-fPsize);

   fX0 = fgScrollBarWidth + fPos * fSliderRange / TMath::Max(fRange-fPsize, 1);
   fX0 = TMath::Max(fX0, fgScrollBarWidth);
   fX0 = TMath::Min(fX0, fgScrollBarWidth + fSliderRange);

   fSlider->Move(fX0, 0);
   fSlider->Resize(fSliderSize, fgScrollBarWidth);
   fClient->NeedRedraw(fSlider);

   SendMessage(fMsgWindow, MK_MSG(kC_HSCROLL, kSB_SLIDERPOS), fPos, 0);
   PositionChanged(fPos);
}


////////////////////////////////////////////////////////////////////////////////
/// Create a vertical scrollbar.

TGVScrollBar::TGVScrollBar(const TGWindow *p, UInt_t w, UInt_t h,
                           UInt_t options, ULong_t back) :
    TGScrollBar(p, w, h, options, back)
{
   fHeadPic = fClient->GetPicture("arrow_up.xpm");
   fTailPic = fClient->GetPicture("arrow_down.xpm");

   if (!fHeadPic || !fTailPic) {
      Error("TGVScrollBar", "arrow_*.xpm not found");
      return;
   }
   fHead   = new TGScrollBarElement(this, fHeadPic, fgScrollBarWidth, fgScrollBarWidth,
                                    kRaisedFrame);
   fTail   = new TGScrollBarElement(this, fTailPic, fgScrollBarWidth, fgScrollBarWidth,
                                    kRaisedFrame);
   fSlider = new TGScrollBarElement(this, 0, fgScrollBarWidth, 50,
                                    kRaisedFrame);

   gVirtualX->GrabButton(fId, kAnyButton, kAnyModifier, kButtonPressMask |
                         kButtonReleaseMask | kPointerMotionMask, kNone, kNone);

   fDragging = kFALSE;
   fX0 = fY0 = (fgScrollBarWidth = TMath::Max(fgScrollBarWidth, 5));
   fPos = 0;

   fRange = TMath::Max((Int_t) h - (fgScrollBarWidth << 1), 1);
   fPsize = fRange >> 1;

   fSliderSize  = 50;
   fSliderRange = 1;

   fHead->SetEditDisabled(kEditDisable | kEditDisableGrab);
   fTail->SetEditDisabled(kEditDisable | kEditDisableGrab);
   fSlider->SetEditDisabled(kEditDisable | kEditDisableGrab);
   fEditDisabled = kEditDisableLayout | kEditDisableWidth | kEditDisableBtnEnable;
}

////////////////////////////////////////////////////////////////////////////////
/// Layout and move vertical scrollbar components.

void TGVScrollBar::Layout()
{
   // Should recalculate fSliderSize, fSliderRange, fX0, fY0, etc. too...
   fHead->Move(0, 0);
   fHead->Resize(fgScrollBarWidth, fgScrollBarWidth);
   fTail->Move(0, fHeight-fgScrollBarWidth);
   fTail->Resize(fgScrollBarWidth, fgScrollBarWidth);

   if (fSlider->GetY() != fY0) {
      fSlider->Move(0, fY0);
      fSlider->Resize(fgScrollBarWidth, 50);
      fClient->NeedRedraw(fSlider);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button event in vertical scrollbar.

Bool_t TGVScrollBar::HandleButton(Event_t *event)
{
   Int_t newpos;

   if (event->fCode == kButton4) {
      if (!fHead->IsEnabled()) {
         return kFALSE;
      }
      //scroll up
      newpos = fPos - fPsize;
      if (newpos<0) newpos = 0;
      SetPosition(newpos);
      return kTRUE;
   }
   if (event->fCode == kButton5) {
      if (!fTail->IsEnabled()) {
         return kFALSE;
      }

      // scroll down
      newpos = fPos + fPsize;
      SetPosition(newpos);
      return kTRUE;
   }

   if (event->fType == kButtonPress) {
      if (event->fCode == kButton3) {
         fY0 = event->fY - fSliderSize/2;
         fY0 = TMath::Max(fY0, fgScrollBarWidth);
         fY0 = TMath::Min(fY0, fgScrollBarWidth + fSliderRange);
         ULong_t pos = (ULong_t)(fY0 - fgScrollBarWidth) * (ULong_t)(fRange-fPsize) / (ULong_t)fSliderRange;
         fPos = (Int_t)pos;

         fPos = TMath::Max(fPos, 0);
         fPos = TMath::Min(fPos, fRange-fPsize);
         fSlider->Move(0, fY0);

         SendMessage(fMsgWindow, MK_MSG(kC_VSCROLL, kSB_SLIDERTRACK), fPos, 0);
         PositionChanged(fPos);
         return kTRUE;
      }

      // fUser[0] contains the subwindow (child) in which the event occured
      // (see GX11Gui.cxx)
      Window_t subw = (Window_t)event->fUser[0];

      if (subw == fSlider->GetId()) {
         fXp = event->fX - fX0;
         fYp = event->fY - fY0;
         fDragging = kTRUE;

      } else {

         if (!fRepeat)
            fRepeat = new TSBRepeatTimer(this, 400, fSmallInc); // 500
         fRepeat->Reset();
         gSystem->AddTimer(fRepeat);
         fSubw = subw;

         if (subw == fHead->GetId()) {
            //if (!fHead->IsEnabled()) {
            //   return kFALSE;
           // }
            fHead->SetState(kButtonDown);
            fPos -= fSmallInc;
         } else if (subw == fTail->GetId()) {
            //if (!fTail->IsEnabled()) {
            //   return kFALSE;
            //}
            fTail->SetState(kButtonDown);
            fPos += fSmallInc;
         } else if (event->fY > fgScrollBarWidth && event->fY < fY0)
            fPos -= fPsize;
         else if (event->fY > fY0+fSliderSize && event->fY < (Int_t)fHeight-fgScrollBarWidth)
            fPos += fPsize;

         fPos = TMath::Max(fPos, 0);
         fPos = TMath::Min(fPos, fRange-fPsize);

         ULong_t y0 = (ULong_t)fgScrollBarWidth + (ULong_t)fPos * (ULong_t)fSliderRange / (ULong_t)TMath::Max(fRange-fPsize, 1);
         fY0 = (Int_t)y0;

         fY0 = TMath::Max(fY0, fgScrollBarWidth);
         fY0 = TMath::Min(fY0, fgScrollBarWidth + fSliderRange);

         fSlider->Move(0, fY0);

         SendMessage(fMsgWindow, MK_MSG(kC_VSCROLL, kSB_SLIDERTRACK), fPos, 0);
         PositionChanged(fPos);
      }

      // last argument kFALSE forces all specified events to this window
      if (fGrabPointer && !fClient->IsEditable())
         gVirtualX->GrabPointer(fId, kButtonPressMask | kButtonReleaseMask |
                                kPointerMotionMask, kNone, kNone,
                                kTRUE, kFALSE);
   } else {
      fHead->SetState(kButtonUp);
      fTail->SetState(kButtonUp);

      if (fRepeat) {
         fRepeat->Remove();
         fRepeat->SetTime(400);  // might have been shortened in HandleTimer()
         fSmallInc = ((TSBRepeatTimer*)fRepeat)->GetSmallInc();
      }

      fDragging = kFALSE;

      fPos = TMath::Max(fPos, 0);
      fPos = TMath::Min(fPos, fRange-fPsize);

      SendMessage(fMsgWindow, MK_MSG(kC_VSCROLL, kSB_SLIDERPOS), fPos, 0);
      PositionChanged(fPos);

      if (fGrabPointer) {
         gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);  // ungrab pointer
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse motion in a vertical scrollbar.

Bool_t TGVScrollBar::HandleMotion(Event_t *event)
{
   if (fDragging) {
      fX0 = event->fX - fXp;
      fY0 = event->fY - fYp;

      fY0 = TMath::Max(fY0, fgScrollBarWidth);
      fY0 = TMath::Min(fY0, fgScrollBarWidth + fSliderRange);
      fSlider->Move(0, fY0);
      ULong_t pos = (ULong_t)(fY0 - fgScrollBarWidth) * (ULong_t)(fRange-fPsize) / fSliderRange;
      fPos = (Int_t)pos;

      fPos = TMath::Max(fPos, 0);
      fPos = TMath::Min(fPos, fRange-fPsize);

      SendMessage(fMsgWindow, MK_MSG(kC_VSCROLL, kSB_SLIDERTRACK), fPos, 0);
      PositionChanged(fPos);
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set range of vertical scrollbar.

void TGVScrollBar::SetRange(Int_t range, Int_t page_size)
{
   fRange = TMath::Max(range, 1);
   fPsize = TMath::Max(page_size, 0);
   fPos = TMath::Max(fPos, 0);
   fPos = TMath::Min(fPos, fRange-fPsize);

   fSliderSize = TMath::Max(fPsize * (fHeight - (fgScrollBarWidth << 1)) /
                            fRange, (UInt_t) 6);
   fSliderSize = TMath::Min((UInt_t)fSliderSize, fHeight - (fgScrollBarWidth << 1));

   fSliderRange = TMath::Max(fHeight - (fgScrollBarWidth << 1) - fSliderSize,
                             (UInt_t)1);

   ULong_t y0 = (ULong_t)fgScrollBarWidth + (ULong_t)fPos * (ULong_t)fSliderRange / (ULong_t)TMath::Max(fRange-fPsize, 1);
   fY0 = (Int_t)y0;
   fY0 = TMath::Max(fY0, fgScrollBarWidth);
   fY0 = TMath::Min(fY0, fgScrollBarWidth + fSliderRange);

   fSlider->Move(0, fY0);
   fSlider->Resize(fgScrollBarWidth, fSliderSize);
   fClient->NeedRedraw(fSlider);

   //  fPos = (fY0 - fgScrollBarWidth) * (fRange-fPsize) / fSliderRange;


   SendMessage(fMsgWindow, MK_MSG(kC_VSCROLL, kSB_SLIDERPOS), fPos, 0);
   PositionChanged(fPos);
   RangeChanged(fRange);
   PageSizeChanged(fPsize);
}

////////////////////////////////////////////////////////////////////////////////
/// Set logical slider position of vertical scrollbar.

void TGVScrollBar::SetPosition(Int_t pos)
{
   fPos = TMath::Max(pos, 0);
   fPos = TMath::Min(pos, fRange-fPsize);

   ULong_t y0 = (ULong_t)fgScrollBarWidth + (ULong_t)fPos * (ULong_t)fSliderRange / (ULong_t)TMath::Max(fRange-fPsize, 1);
   fY0 = (Int_t)y0;
   fY0 = TMath::Max(fY0, fgScrollBarWidth);
   fY0 = TMath::Min(fY0, fgScrollBarWidth + fSliderRange);

   fSlider->Move(0, fY0);
   fSlider->Resize(fgScrollBarWidth, fSliderSize);
   fClient->NeedRedraw(fSlider);

   SendMessage(fMsgWindow, MK_MSG(kC_VSCROLL, kSB_SLIDERPOS), fPos, 0);
   PositionChanged(fPos);
}

////////////////////////////////////////////////////////////////////////////////
/// Save an horizontal scrollbar as a C++ statement(s) on output stream out.

void TGHScrollBar::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   if (fBackground != GetDefaultFrameBackground()) SaveUserColor(out, option);

   out <<"   TGHScrollBar *";
   out << GetName() << " = new TGHScrollBar(" << fParent->GetName()
       << "," << GetWidth() << "," << GetHeight();

   if (fBackground == GetDefaultFrameBackground()) {
      if (!GetOptions()) {
         out <<");" << std::endl;
      } else {
         out << "," << GetOptionString() <<");" << std::endl;
      }
   } else {
      out << "," << GetOptionString() << ",ucolor);" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   out << "   " << GetName() <<"->SetRange(" << GetRange() << "," << GetPageSize() << ");" << std::endl;
   out << "   " << GetName() <<"->SetPosition(" << GetPosition() << ");" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Save an vertical scrollbar as a C++ statement(s) on output stream out.

void TGVScrollBar::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   if (fBackground != GetDefaultFrameBackground()) SaveUserColor(out, option);

   out<<"   TGVScrollBar *";
   out << GetName() <<" = new TGVScrollBar("<< fParent->GetName()
       << "," << GetWidth() << "," << GetHeight();

   if (fBackground == GetDefaultFrameBackground()) {

      if (!GetOptions()) {
         out <<");" << std::endl;
      } else {
         out << "," << GetOptionString() <<");" << std::endl;
      }
   } else {
      out << "," << GetOptionString() << ",ucolor);" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   out << "   " << GetName() <<"->SetRange(" << GetRange() << "," << GetPageSize() << ");" << std::endl;
   out << "   " << GetName() <<"->SetPosition(" << GetPosition() << ");" << std::endl;
}
