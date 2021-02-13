// @(#)root/gui:$Id$
// Author: Fons Rademakers   14/01/98

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
// TGSlider, TGVSlider and TGHSlider                                    //
//                                                                      //
// Slider widgets allow easy selection of a range.                      //
// Sliders can be either horizontal or vertical oriented and there is   //
// a choice of two different slider types and three different types     //
// of tick marks.                                                       //
//                                                                      //
// TGSlider is an abstract base class. Use the concrete TGVSlider and   //
// TGHSlider.                                                           //
//                                                                      //
// Dragging the slider will generate the event:                         //
// kC_VSLIDER, kSL_POS, slider id, position  (for vertical slider)      //
// kC_HSLIDER, kSL_POS, slider id, position  (for horizontal slider)    //
//                                                                      //
// Pressing the mouse will generate the event:                          //
// kC_VSLIDER, kSL_PRESS, slider id, 0  (for vertical slider)           //
// kC_HSLIDER, kSL_PRESS, slider id, 0  (for horizontal slider)         //
//                                                                      //
// Releasing the mouse will generate the event:                         //
// kC_VSLIDER, kSL_RELEASE, slider id, 0  (for vertical slider)         //
// kC_HSLIDER, kSL_RELEASE, slider id, 0  (for horizontal slider)       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGSlider.h"
#include "TGPicture.h"
#include "TImage.h"
#include "TEnv.h"
#include "TVirtualX.h"

#include <iostream>

ClassImp(TGSlider);
ClassImp(TGVSlider);
ClassImp(TGHSlider);


////////////////////////////////////////////////////////////////////////////////
/// Slider constructor.

TGSlider::TGSlider(const TGWindow *p, UInt_t w, UInt_t h, UInt_t type, Int_t id,
                   UInt_t options, ULong_t back)
   : TGFrame(p, w, h, options, back)
{
   fDisabledPic = 0;
   fWidgetId    = id;
   fWidgetFlags = kWidgetWantFocus | kWidgetIsEnabled;
   fMsgWindow   = p;

   fType     = type;
   fScale    = 10;
   fDragging = kFALSE;
   fPos = fRelPos = 0;
   fVmax = fVmin = 0;
   fSliderPic = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates disabled picture.

void TGSlider::CreateDisabledPicture()
{
   if (!fSliderPic) return;

   TImage *img = TImage::Create();
   if (!img) return;
   TImage *img2 = TImage::Create();
   if (!img2) {
      if (img) delete img;
      return;
   }
   TString back = gEnv->GetValue("Gui.BackgroundColor", "#c0c0c0");
   img2->FillRectangle(back.Data(), 0, 0, fSliderPic->GetWidth(),
                       fSliderPic->GetHeight());
   img->SetImage(fSliderPic->GetPicture(), fSliderPic->GetMask());
   Pixmap_t mask = img->GetMask();
   img2->Merge(img, "overlay");

   TString name = "disbl_";
   name += fSliderPic->GetName();
   fDisabledPic = fClient->GetPicturePool()->GetPicture(name.Data(),
                                             img2->GetPixmap(), mask);
   delete img;
   delete img2;
}

////////////////////////////////////////////////////////////////////////////////
/// Set state of widget. If kTRUE=enabled, kFALSE=disabled.

void TGSlider::SetState(Bool_t state)
{
   if (state) {
      SetFlags(kWidgetIsEnabled);
   } else {
      ClearFlags(kWidgetIsEnabled);
   }
   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a vertical slider widget.

TGVSlider::TGVSlider(const TGWindow *p, UInt_t h, UInt_t type, Int_t id,
                     UInt_t options, ULong_t back) :
   TGSlider(p, kSliderWidth, h, type, id, options, back)
{
   if ((fType & kSlider1))
      fSliderPic = fClient->GetPicture("slider1h.xpm");
   else
      fSliderPic = fClient->GetPicture("slider2h.xpm");

   if (!fSliderPic)
      Error("TGVSlider", "slider?h.xpm not found");

   CreateDisabledPicture();

   gVirtualX->GrabButton(fId, kAnyButton, kAnyModifier,
                         kButtonPressMask | kButtonReleaseMask |
                         kPointerMotionMask, kNone, kNone);

   AddInput(kStructureNotifyMask);
   // set initial values
   fPos = h/2; fVmin = 0; fVmax = h; fYp = 0;
   fEditDisabled = kEditDisableWidth;

   if (!p && fClient->IsEditable()) {
      Resize(GetDefaultWidth(), 100);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Delete vertical slider widget.

TGVSlider::~TGVSlider()
{
   if (fSliderPic) fClient->FreePicture(fSliderPic);
   if (fDisabledPic) fClient->FreePicture(fDisabledPic);
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw vertical slider widget.

void TGVSlider::DoRedraw()
{
   // cleanup the drawable
   gVirtualX->ClearWindow(fId);

   GContext_t drawGC = IsEnabled() ? GetBlackGC()() : GetShadowGC()();

   gVirtualX->DrawLine(fId, GetShadowGC()(), fWidth/2, 8, fWidth/2-1, 8);
   gVirtualX->DrawLine(fId, GetShadowGC()(), fWidth/2-1, 8, fWidth/2-1, fHeight-9);
   gVirtualX->DrawLine(fId, GetHilightGC()(), fWidth/2+1, 8, fWidth/2+1, fHeight-8);
   gVirtualX->DrawLine(fId, GetHilightGC()(), fWidth/2+1, fHeight-8, fWidth/2, fHeight-8);
   gVirtualX->DrawLine(fId, drawGC, fWidth/2, 9, fWidth/2, fHeight-9);

   // check scale
   if (fScale == 1) fScale++;
   if (fScale * 2 > (int)fHeight) fScale = 0;
   if (fScale > 0 && !(fType & kScaleNo)) {
      int lines = ((int)fHeight-16) / fScale;
      int remain = ((int)fHeight-16) % fScale;
      if (lines < 1) lines = 1;
      for (int i = 0; i <= lines; i++) {
         int y = i * fScale + (i * remain) / lines;
         gVirtualX->DrawLine(fId, drawGC, fWidth/2+8, y+7, fWidth/2+10, y+7);
         if ((fType & kSlider2) && (fType & kScaleBoth))
            gVirtualX->DrawLine(fId, drawGC, fWidth/2-9, y+7, fWidth/2-11, y+7);
      }
   }
   if (fPos < fVmin) fPos = fVmin;
   if (fPos > fVmax) fPos = fVmax;

   // calc slider-picture position
   fRelPos = (((int)fHeight-16) * (fPos - fVmin)) / (fVmax - fVmin) + 8;
   const TGPicture *pic = fSliderPic;
   if (!IsEnabled()) {
      if (!fDisabledPic) CreateDisabledPicture();
      pic = fDisabledPic ? fDisabledPic : fSliderPic;
   }
   if (pic) pic->Draw(fId, GetBckgndGC()(), fWidth/2-7, fRelPos-6);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button event in vertical slider.

Bool_t TGVSlider::HandleButton(Event_t *event)
{
   if (!IsEnabled()) return kTRUE;
   if (event->fCode == kButton4 || event->fCode == kButton5) {
      Int_t oldPos = fPos;
      int m = (fVmax - fVmin) / (fWidth-16);
      if (event->fCode == kButton4)
         fPos -= ((m) ? m : 1);
      else if (event->fCode == kButton5)
         fPos += ((m) ? m : 1);
      if (fPos > fVmax) fPos = fVmax;
      if (fPos < fVmin) fPos = fVmin;
      SendMessage(fMsgWindow, MK_MSG(kC_HSLIDER, kSL_POS),
                  fWidgetId, fPos);
      fClient->ProcessLine(fCommand, MK_MSG(kC_HSLIDER, kSL_POS),
                           fWidgetId, fPos);
      if (fPos != oldPos) {
         PositionChanged(fPos);
         fClient->NeedRedraw(this);
      }
      return kTRUE;
   }
   if (event->fType == kButtonPress) {
      // constrain to the slider width
      if (event->fX < (Int_t)fWidth/2-7 || event->fX > (Int_t)fWidth/2+7) {
         return kTRUE;
      }
      // last argument kFALSE forces all specified events to this window
      gVirtualX->GrabPointer(fId, kButtonPressMask | kButtonReleaseMask |
                             kPointerMotionMask, kNone, kNone,
                             kTRUE, kFALSE);

      if (event->fY >= fRelPos - 7 && event->fY <= fRelPos + 7) {
         // slider selected
         fDragging = kTRUE;
         fYp = event->fY - (fRelPos-7);
         SendMessage(fMsgWindow, MK_MSG(kC_VSLIDER, kSL_PRESS), fWidgetId, 0);
         fClient->ProcessLine(fCommand, MK_MSG(kC_VSLIDER, kSL_PRESS), fWidgetId, 0);
         Pressed();
      } else {
         if (event->fCode == kButton1) {
            // scroll up or down
            int m = (fVmax - fVmin) / (fHeight-16);
            if (event->fY < fRelPos) {
               fPos -= ((m) ? m : 1);
            }
            if (event->fY > fRelPos) {
               fPos += ((m) ? m : 1);
            }
         } else if (event->fCode == kButton2) {
            // set absolute position
            fPos = ((fVmax - fVmin) * event->fY) / (fHeight-16) + fVmin;
         }
         if (fPos > fVmax) fPos = fVmax;
         if (fPos < fVmin) fPos = fVmin;
         SendMessage(fMsgWindow, MK_MSG(kC_VSLIDER, kSL_POS),
                     fWidgetId, fPos);
         fClient->ProcessLine(fCommand, MK_MSG(kC_VSLIDER, kSL_POS),
                              fWidgetId, fPos);
         PositionChanged(fPos);
      }
      fClient->NeedRedraw(this);

   } else {
      // ButtonRelease
      fDragging = kFALSE;
      gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);  // ungrab pointer

      SendMessage(fMsgWindow, MK_MSG(kC_VSLIDER, kSL_RELEASE), fWidgetId, 0);
      fClient->ProcessLine(fCommand, MK_MSG(kC_VSLIDER, kSL_RELEASE), fWidgetId, 0);
      Released();
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse motion event in vertical slider.

Bool_t TGVSlider::HandleMotion(Event_t *event)
{
   if (fDragging) {
      int old = fPos;
      fPos = ((fVmax - fVmin) * (event->fY - fYp)) / ((int)fHeight-16) + fVmin;
      if (fPos > fVmax) fPos = fVmax;
      if (fPos < fVmin) fPos = fVmin;

      // check if position changed
      if (old != fPos) {
         fClient->NeedRedraw(this);
         SendMessage(fMsgWindow, MK_MSG(kC_VSLIDER, kSL_POS),
                     fWidgetId, fPos);
         fClient->ProcessLine(fCommand, MK_MSG(kC_VSLIDER, kSL_POS),
                              fWidgetId, fPos);
         PositionChanged(fPos);
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handles resize events for this widget.

Bool_t TGVSlider::HandleConfigureNotify(Event_t* event)
{
   TGFrame::HandleConfigureNotify(event);
   fClient->NeedRedraw(this);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create horizontal slider widget.

TGHSlider::TGHSlider(const TGWindow *p, UInt_t w, UInt_t type, Int_t id,
                     UInt_t options, ULong_t back) :
   TGSlider(p, w, kSliderHeight, type, id, options, back)
{
   if ((fType & kSlider1))
      fSliderPic = fClient->GetPicture("slider1v.xpm");
   else
      fSliderPic = fClient->GetPicture("slider2v.xpm");

   if (!fSliderPic)
      Error("TGHSlider", "slider?v.xpm not found");

   CreateDisabledPicture();

   gVirtualX->GrabButton(fId, kAnyButton, kAnyModifier,
                         kButtonPressMask | kButtonReleaseMask |
                         kPointerMotionMask, kNone, kNone);

   AddInput(kStructureNotifyMask);
   // set initial values
   fPos = w/2; fVmin = 0; fVmax = w; fXp = 0;
   fEditDisabled = kEditDisableHeight;

   if (!p && fClient->IsEditable()) {
      Resize(100, GetDefaultHeight());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Delete a horizontal slider widget.

TGHSlider::~TGHSlider()
{
   if (fSliderPic) fClient->FreePicture(fSliderPic);
   if (fDisabledPic) fClient->FreePicture(fDisabledPic);
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw horizontal slider widget.

void TGHSlider::DoRedraw()
{
   // cleanup drawable
   gVirtualX->ClearWindow(fId);

   GContext_t drawGC = IsEnabled() ? GetBlackGC()() : GetShadowGC()();

   gVirtualX->DrawLine(fId, GetShadowGC()(), 8, fHeight/2, 8, fHeight/2-1);
   gVirtualX->DrawLine(fId, GetShadowGC()(), 8, fHeight/2-1, fWidth-9, fHeight/2-1);
   gVirtualX->DrawLine(fId, GetHilightGC()(), 8, fHeight/2+1, fWidth-8, fHeight/2+1);
   gVirtualX->DrawLine(fId, GetHilightGC()(), fWidth-8, fHeight/2+1, fWidth-8, fHeight/2);
   gVirtualX->DrawLine(fId, drawGC, 9, fHeight/2, fWidth-9, fHeight/2);

   if (fScale == 1) fScale++;
   if (fScale * 2 > (int)fWidth) fScale = 0;
   if (fScale > 0 && !(fType & kScaleNo)) {
      int lines = ((int)fWidth-16) / fScale;
      int remain = ((int)fWidth-16) % fScale;
      if (lines < 1) lines = 1;
      for (int i = 0; i <= lines; i++) {
         int x = i * fScale + (i * remain) / lines;
         gVirtualX->DrawLine(fId, drawGC, x+7, fHeight/2+8, x+7, fHeight/2+10);
         if ((fType & kSlider2) && (fType & kScaleBoth))
            gVirtualX->DrawLine(fId, drawGC, x+7, fHeight/2-9, x+7, fHeight/2-11);
      }
   }
   if (fPos < fVmin) fPos = fVmin;
   if (fPos > fVmax) fPos = fVmax;

   // calc slider-picture position
   fRelPos = (((int)fWidth-16) * (fPos - fVmin)) / (fVmax - fVmin) + 8;
   const TGPicture *pic = fSliderPic;
   if (!IsEnabled()) {
      if (!fDisabledPic) CreateDisabledPicture();
      pic = fDisabledPic ? fDisabledPic : fSliderPic;
   }
   if (pic) pic->Draw(fId, GetBckgndGC()(), fRelPos-6, fHeight/2-7);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button event in horizontal slider widget.

Bool_t TGHSlider::HandleButton(Event_t *event)
{
   if (!IsEnabled()) return kTRUE;
   if (event->fCode == kButton4 || event->fCode == kButton5) {
      Int_t oldPos = fPos;
      int m = (fVmax - fVmin) / (fWidth-16);
      if (event->fCode == kButton4)
         fPos += ((m) ? m : 1);
      else if (event->fCode == kButton5)
         fPos -= ((m) ? m : 1);
      if (fPos > fVmax) fPos = fVmax;
      if (fPos < fVmin) fPos = fVmin;
      SendMessage(fMsgWindow, MK_MSG(kC_HSLIDER, kSL_POS),
                  fWidgetId, fPos);
      fClient->ProcessLine(fCommand, MK_MSG(kC_HSLIDER, kSL_POS),
                           fWidgetId, fPos);
      if (fPos != oldPos) {
         PositionChanged(fPos);
         fClient->NeedRedraw(this);
      }
      return kTRUE;
   }
   if (event->fType == kButtonPress) {
      // constrain to the slider height
      if (event->fY < (Int_t)fHeight/2-7 || event->fY > (Int_t)fHeight/2+7) {
         return kTRUE;
      }
      if (event->fX >= fRelPos - 7 && event->fX <= fRelPos + 7) {
         // slider selected
         fDragging = kTRUE;
         fXp = event->fX - (fRelPos-7);
         SendMessage(fMsgWindow, MK_MSG(kC_HSLIDER, kSL_PRESS), fWidgetId, 0);
         fClient->ProcessLine(fCommand, MK_MSG(kC_HSLIDER, kSL_PRESS), fWidgetId, 0);
         Pressed();
      } else {
         if (event->fCode == kButton1) {
            int m = (fVmax - fVmin) / (fWidth-16);
            if (event->fX < fRelPos) {
               fPos -= ((m) ? m : 1);
            }
            if (event->fX > fRelPos) {
               fPos += ((m) ? m : 1);
            }
         } else if (event->fCode == kButton2) {
            fPos = ((fVmax - fVmin) * event->fX) / (fWidth-16) + fVmin;
         }
         if (fPos > fVmax) fPos = fVmax;
         if (fPos < fVmin) fPos = fVmin;
         SendMessage(fMsgWindow, MK_MSG(kC_HSLIDER, kSL_POS),
                     fWidgetId, fPos);
         fClient->ProcessLine(fCommand, MK_MSG(kC_HSLIDER, kSL_POS),
                              fWidgetId, fPos);
         PositionChanged(fPos);
      }
      fClient->NeedRedraw(this);

      // last argument kFALSE forces all specified events to this window
      gVirtualX->GrabPointer(fId, kButtonPressMask | kButtonReleaseMask | kPointerMotionMask,
                             kNone, kNone, kTRUE, kFALSE);
   } else {
      // ButtonRelease
      fDragging = kFALSE;
      gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);  // ungrab pointer

      SendMessage(fMsgWindow, MK_MSG(kC_HSLIDER, kSL_RELEASE), fWidgetId, 0);
      fClient->ProcessLine(fCommand, MK_MSG(kC_HSLIDER, kSL_RELEASE), fWidgetId, 0);
      Released();
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse motion event in horizontal slide widget.

Bool_t TGHSlider::HandleMotion(Event_t *event)
{
   if (fDragging) {
      int old = fPos;
      fPos = ((fVmax - fVmin) * (event->fX - fXp)) / ((int)fWidth-16) + fVmin;
      if (fPos > fVmax) fPos = fVmax;
      if (fPos < fVmin) fPos = fVmin;

      // check if position changed
      if (old != fPos) {
         fClient->NeedRedraw(this);
         SendMessage(fMsgWindow, MK_MSG(kC_HSLIDER, kSL_POS),
                     fWidgetId, fPos);
         fClient->ProcessLine(fCommand, MK_MSG(kC_HSLIDER, kSL_POS),
                              fWidgetId, fPos);
         PositionChanged(fPos);
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handles resize events for this widget.

Bool_t TGHSlider::HandleConfigureNotify(Event_t* event)
{
   TGFrame::HandleConfigureNotify(event);
   fClient->NeedRedraw(this);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the slider type as a string - used in SavePrimitive().

TString TGSlider::GetTypeString() const
{
   TString stype;

   if (fType) {
      if (fType & kSlider1)  {
         if (stype.Length() == 0) stype  = "kSlider1";
         else                     stype += " | kSlider1";
      }
      if (fType & kSlider2)  {
         if (stype.Length() == 0) stype  = "kSlider2";
         else                     stype += " | kSlider2";
      }
      if (fType & kScaleNo)  {
         if (stype.Length() == 0) stype  = "kScaleNo";
         else                     stype += " | kScaleNo";
      }
      if (fType & kScaleDownRight) {
         if (stype.Length() == 0) stype  = "kScaleDownRight";
         else                     stype += " | kScaleDownRight";
      }
      if (fType & kScaleBoth) {
         if (stype.Length() == 0) stype  = "kScaleBoth";
         else                     stype += " | kScaleBoth";
      }
   }
   return stype;
}

////////////////////////////////////////////////////////////////////////////////
/// Save an horizontal slider as a C++ statement(s) on output stream out.

void TGHSlider::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   if (fBackground != GetDefaultFrameBackground()) SaveUserColor(out, option);

   out <<"   TGHSlider *";
   out << GetName() << " = new TGHSlider(" << fParent->GetName()
       << "," << GetWidth() << ",";
   out << GetTypeString() << "," << WidgetId();

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

   if (fVmin != 0 || fVmax != (Int_t)fWidth)
      out << "   " << GetName() <<"->SetRange(" << fVmin << "," << fVmax << ");" << std::endl;

   if (fPos != (Int_t)fWidth/2)
      out << "   " << GetName() <<"->SetPosition(" << GetPosition() << ");" << std::endl;

   if (fScale != 10)
      out << "   " << GetName() <<"->SetScale(" << fScale << ");" << std::endl;

   if (!IsEnabled())
      out << "   " << GetName() <<"->SetState(kFALSE);" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Save an horizontal slider as a C++ statement(s) on output stream out.

void TGVSlider::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   if (fBackground != GetDefaultFrameBackground()) SaveUserColor(out, option);

   out<<"   TGVSlider *";
   out << GetName() <<" = new TGVSlider("<< fParent->GetName()
       << "," << GetHeight() << ",";
   out << GetTypeString() << "," << WidgetId();

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

   if (fVmin != 0 || fVmax != (Int_t)fHeight)
      out << "   " << GetName() <<"->SetRange(" << fVmin << "," << fVmax << ");" << std::endl;

   if (fPos != (Int_t)fHeight/2)
      out << "   " << GetName() <<"->SetPosition(" << GetPosition() << ");" << std::endl;

   if (fScale != 10)
      out << "   " << GetName() <<"->SetScale(" << fScale << ");" << std::endl;

   if (!IsEnabled())
      out << "   " << GetName() <<"->SetState(kFALSE);" << std::endl;
}
