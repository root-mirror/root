// @(#)root/gui:$Id$
// Author: Bertrand Bellenot   20/01/06

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
// TGTripleVSlider and TGTripleHSlider                                  //
//                                                                      //
// TripleSlider inherit from DoubleSlider widgets and allow easy        //
// selection of a min, max and pointer value out of a range.            //
// The pointer position can be constrained to edges of slider and / or  //
// can be relative to the slider position.                              //
//                                                                      //
// To change the min value press the mouse near to the left / bottom    //
// edge of the slider.                                                  //
// To change the max value press the mouse near to the right / top      //
// edge of the slider.                                                  //
// To change both values simultaneously press the mouse near to the     //
// center of the slider.                                                //
// To change pointer value press the mouse on the pointer and drag it   //
// to the desired position                                              //
//                                                                      //
// Dragging the slider will generate the event:                         //
// kC_VSLIDER, kSL_POS, slider id, 0  (for vertical slider)             //
// kC_HSLIDER, kSL_POS, slider id, 0  (for horizontal slider)           //
//                                                                      //
// Pressing the mouse will generate the event:                          //
// kC_VSLIDER, kSL_PRESS, slider id, 0  (for vertical slider)           //
// kC_HSLIDER, kSL_PRESS, slider id, 0  (for horizontal slider)         //
//                                                                      //
// Releasing the mouse will generate the event:                         //
// kC_VSLIDER, kSL_RELEASE, slider id, 0  (for vertical slider)         //
// kC_HSLIDER, kSL_RELEASE, slider id, 0  (for horizontal slider)       //
//                                                                      //
// Moving the pointer will generate the event:                          //
// kC_VSLIDER, kSL_POINTER, slider id, 0  (for vertical slider)         //
// kC_HSLIDER, kSL_POINTER, slider id, 0  (for horizontal slider)       //
//                                                                      //
// Use the functions GetMinPosition(), GetMaxPosition() and             //
// GetPosition() to retrieve the position of the slider.                //
// Use the function GetPointerPosition() to retrieve the position of    //
// the pointer                                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGTripleSlider.h"
#include "TGPicture.h"
#include "TSystem.h"
#include "TVirtualX.h"

#include <iostream>
#include <cstdlib>

ClassImp(TGTripleVSlider);
ClassImp(TGTripleHSlider);

////////////////////////////////////////////////////////////////////////////////
/// Create a vertical slider widget.

TGTripleVSlider::TGTripleVSlider(const TGWindow *p, UInt_t h, UInt_t type, Int_t id,
                                 UInt_t options, ULong_t back,
                                 Bool_t reversed, Bool_t mark_ends,
                                 Bool_t constrained, Bool_t relative)
    : TGDoubleVSlider(p, h, type, id, options, back, reversed, mark_ends)
{
   fPointerPic = fClient->GetPicture("slider1h.xpm");
   if (!fPointerPic)
      Error("TGTripleVSlider", "slider1h.xpm not found");
   fConstrained = constrained;
   fRelative = relative;
   fCz = 0;
   fSCz = 0;
   AddInput(kStructureNotifyMask);
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete vertical slider widget.

TGTripleVSlider::~TGTripleVSlider()
{
   if (fPointerPic) fClient->FreePicture(fPointerPic);
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw vertical slider widget.

void TGTripleVSlider::DoRedraw()
{
   TGDoubleVSlider::DoRedraw();
   // Draw Pointer
   DrawPointer();
}

////////////////////////////////////////////////////////////////////////////////
/// Draw slider pointer

void TGTripleVSlider::DrawPointer()
{
   if (fPointerPic) fPointerPic->Draw(fId, GetBckgndGC()(), fWidth/2-7, fCz-5);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button event in vertical slider.

Bool_t TGTripleVSlider::HandleButton(Event_t *event)
{
   if (event->fType == kButtonPress && event->fCode == kButton1) {
      // constrain to the slider width
      if (event->fX < (Int_t)fWidth/2-7 || event->fX > (Int_t)fWidth/2+7) {
         return kTRUE;
      }
      fPressPoint = event->fY;
      fPressSmin  = fSmin;
      fPressSmax  = fSmax;

      int relMin = (int)((fHeight-16) * (fSmin - fVmin) / (fVmax - fVmin)) + 1;
      int relMax = (int)((fHeight-16) * (fSmax - fVmin) / (fVmax - fVmin) + 15);
      if (fPressPoint > (fCz - 5) && fPressPoint < (fCz + 5) &&
          event->fX > ((Int_t)fWidth / 2) - 7 && event->fX < ((Int_t)fWidth / 2) + 5)
         // move pointer
         fMove = 4;
      else if (fPressPoint < (relMax - relMin) / 4 + relMin)
         // move only min value
         fMove = 1;
      else if (fPressPoint > (relMax - relMin) / 4 * 3 + relMin)
         // move only max value
         fMove = 2;
      else
         // move min and max value
         fMove = 3;

      SendMessage(fMsgWindow, MK_MSG(kC_VSLIDER, kSL_PRESS), fWidgetId, 0);
      fClient->ProcessLine(fCommand, MK_MSG(kC_VSLIDER, kSL_PRESS), fWidgetId, 0);
      Pressed();

      // last argument kFALSE forces all specified events to this window
      gVirtualX->GrabPointer(fId, kButtonPressMask | kButtonReleaseMask |
                             kPointerMotionMask, kNone, kNone,
                             kTRUE, kFALSE);
   } else if (event->fType == kButtonRelease && event->fCode == kButton1) {
      SendMessage(fMsgWindow, MK_MSG(kC_VSLIDER, kSL_RELEASE), fWidgetId, 0);
      fClient->ProcessLine(fCommand, MK_MSG(kC_VSLIDER, kSL_RELEASE), fWidgetId, 0);
      Released();
      fMove = 0;
      gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);  // ungrab pointer
   } else
      fMove = 0;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handles resize events for this widget.

Bool_t TGTripleVSlider::HandleConfigureNotify(Event_t* event)
{
   TGFrame::HandleConfigureNotify(event);
   SetPointerPosition(fSCz);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse motion event in vertical slider.

Bool_t TGTripleVSlider::HandleMotion(Event_t *event)
{
   if (fMove < 3) {
      // if the mouse pointer is on the cursor,
      // and we are not moving anything,
      // set the cursor shape as Pointer
      if (event->fY > (fCz - 5) && event->fY < (fCz + 5) &&
          event->fX > ((Int_t)fWidth / 2) - 7 &&
          event->fX < ((Int_t)fWidth / 2) + 5 &&
          fMove == 0)
         gVirtualX->SetCursor(fId, kNone);
      else
         ChangeCursor(event);
   }
   static int oldDiff = 0;
   static Long64_t was = gSystem->Now();
   Long64_t now = gSystem->Now();

   if (fMove == 0)  return kTRUE;
   if ((now-was) < 50) return kTRUE;
   was = now;

   int     diff;
   Double_t oldMin, oldMax;

   diff    = event->fY - fPressPoint;
   oldMin  = fSmin;
   oldMax  = fSmax;

   if (fMove == 1) {
      // change of min value
      oldDiff = 0;
      fSmin = fPressSmin + diff * (fVmax - fVmin) / (fHeight-16);
      if (fSmin < fVmin) fSmin = fVmin;
      if (fSmin > fSmax) fSmin = fSmax;
   } else if (fMove == 2) {
      // change of max value
      oldDiff = 0;
      fSmax = fPressSmax + diff * (fVmax - fVmin) / (fHeight-16);
      if (fSmax > fVmax) fSmax = fVmax;
      if (fSmax < fSmin) fSmax = fSmin;
   } else if (fMove == 3) {
      // change of min and of max value
      Double_t logicalDiff;
      logicalDiff = diff * (fVmax - fVmin) / (fHeight-16);
      if (fPressSmax + logicalDiff > fVmax)
         logicalDiff = fVmax - fPressSmax;
      if (fPressSmin + logicalDiff < fVmin)
         logicalDiff = fVmin - fPressSmin;
      fSmax = fPressSmax + logicalDiff;
      fSmin = fPressSmin + logicalDiff;
      if (fRelative) {
         if (abs(diff) < 3) oldDiff = diff;
         SetPointerPos(diff - oldDiff, 3);
         oldDiff = diff;
      }
   }
   else if (fMove == 4) {
      // change pointer position
      oldDiff = 0;
      SetPointerPos(event->fY, 1);
   }
   if (fMove != 4){
      SetPointerPos(0, 2);
   }
   // check if position has changed
   if (fMove != 0 && (fSmax != oldMax || fSmin != oldMin)) {
      fClient->NeedRedraw(this);
      SendMessage(fMsgWindow, MK_MSG(kC_VSLIDER, kSL_POS), fWidgetId, 0);
      fClient->ProcessLine(fCommand, MK_MSG(kC_VSLIDER, kSL_POS), fWidgetId, 0);
      PositionChanged();
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set pointer position constrained in the slider range.

void TGTripleVSlider::SetConstrained(Bool_t on)
{
   fConstrained = on;

   if (fConstrained) {
      if (GetPointerPosition() <= GetMinPosition())
         SetPointerPos((Int_t)GetMinPosition(), 3);
      else if (GetPointerPosition() >= GetMaxPosition())
         SetPointerPos((Int_t)GetMaxPosition(), 3);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set slider pointer position in pixel value.

void TGTripleVSlider::SetPointerPos(Int_t z, Int_t opt)
{
   static Long64_t was = gSystem->Now();
   Bool_t lcheck = (opt == 1);
   Int_t oldPos = fCz;

   if (opt < 2) {
      fCz = z;

      if (fCz < 7)
         fCz = 7;
      else if (fCz >= (Int_t)fHeight - 7)
         fCz = (Int_t)fHeight - 7;
   }
   if (opt == 3) {
      lcheck = kTRUE;
      fCz += z;
      if (fCz < 7)
         fCz = 7;
      else if (fCz >= (Int_t)fHeight-7)
         fCz = (Int_t)fHeight - 7;
   }
   if (fConstrained) {
      int relMin = (int)((fHeight-16) * (fSmin - fVmin) / (fVmax - fVmin)) + 1;
      int relMax = (int)((fHeight-16) * (fSmax - fVmin) / (fVmax - fVmin) + 15);
      if(fCz < relMin+7) {
         fCz = relMin+7;
         lcheck = kTRUE;
      }
      if(fCz > relMax-7) {
         fCz = relMax-7;
         lcheck = kTRUE;
      }
   }
   if (lcheck)
      fSCz = fVmin + ((Double_t)(fCz-8) * (fVmax - fVmin) / (Double_t)(fHeight-16));
   if(fSCz < fVmin) fSCz = fVmin;
   if(fSCz > fVmax) fSCz = fVmax;
   if (fConstrained) {
      if(fSCz < fSmin) fSCz = fSmin;
      if(fSCz > fSmax) fSCz = fSmax;
   }

   DrawPointer();
   fClient->NeedRedraw(this);
   if (fCz != oldPos) {
      Long64_t now = gSystem->Now();
      if ((fMove != 4) && ((now-was) < 150)) return;
      was = now;
      SendMessage(fMsgWindow, MK_MSG(kC_VSLIDER, kSL_POINTER), fWidgetId, 0);
      fClient->ProcessLine(fCommand, MK_MSG(kC_VSLIDER, kSL_POINTER), fWidgetId, 0);
      PointerPositionChanged();
      fClient->NeedRedraw(this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set pointer position in scaled (real) value

void TGTripleVSlider::SetPointerPosition(Double_t pos)
{
   if (fReversedScale) {
      fSCz = fVmin + fVmax - pos;
   }
   else {
      fSCz = pos;
   }
   Double_t absPos = (fSCz - fVmin) * (fHeight-16) / (fVmax - fVmin);
   SetPointerPos((int)(absPos+5.0), 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Create horizontal slider widget.

TGTripleHSlider::TGTripleHSlider(const TGWindow *p, UInt_t w, UInt_t type, Int_t id,
                                 UInt_t options, ULong_t back,
                                 Bool_t reversed, Bool_t mark_ends,
                                 Bool_t constrained, Bool_t relative)
    : TGDoubleHSlider(p, w, type, id, options, back, reversed, mark_ends)
{
   fPointerPic = fClient->GetPicture("slider1v.xpm");
   if (!fPointerPic)
      Error("TGTripleVSlider", "slider1v.xpm not found");
   fConstrained = constrained;
   fRelative = relative;
   fCz = 0;
   fSCz = 0;
   AddInput(kStructureNotifyMask);
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete a horizontal slider widget.

TGTripleHSlider::~TGTripleHSlider()
{
   if (fPointerPic) fClient->FreePicture(fPointerPic);
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw horizontal slider widget.

void TGTripleHSlider::DoRedraw()
{
   TGDoubleHSlider::DoRedraw();
   // Draw Pointer
   DrawPointer();
}

////////////////////////////////////////////////////////////////////////////////
/// Draw slider pointer

void TGTripleHSlider::DrawPointer()
{
   if (fPointerPic) fPointerPic->Draw(fId, GetBckgndGC()(), fCz-5, fHeight/2-7);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button event in horizontal slider widget.

Bool_t TGTripleHSlider::HandleButton(Event_t *event)
{
   if (event->fType == kButtonPress && event->fCode == kButton1) {
      // constrain to the slider height
      if (event->fY < (Int_t)fHeight/2-7 || event->fY > (Int_t)fHeight/2+7) {
         return kTRUE;
      }
      fPressPoint = event->fX;
      fPressSmin  = fSmin;
      fPressSmax  = fSmax;

      int relMin = (int)((fWidth-16) * (fSmin - fVmin) / (fVmax - fVmin)) + 1;
      int relMax = (int)((fWidth-16) * (fSmax - fVmin) / (fVmax - fVmin) + 15);
      if (fPressPoint > (fCz - 5) && fPressPoint < (fCz + 5) &&
          event->fY > ((Int_t)fHeight / 2) - 7 && event->fY < ((Int_t)fHeight / 2) + 5)
         // move pointer
         fMove = 4;
      else if (fPressPoint < (relMax - relMin) / 4 + relMin)
         // move only min value
         fMove = 1;
      else if (fPressPoint > (relMax - relMin) / 4 * 3 + relMin)
         // move only max value
         fMove = 2;
      else
         // move min and max value
         fMove = 3;

      SendMessage(fMsgWindow, MK_MSG(kC_HSLIDER, kSL_PRESS), fWidgetId, 0);
      fClient->ProcessLine(fCommand, MK_MSG(kC_HSLIDER, kSL_PRESS), fWidgetId, 0);
      Pressed();

      // last argument kFALSE forces all specified events to this window
      gVirtualX->GrabPointer(fId, kButtonPressMask | kButtonReleaseMask |
                             kPointerMotionMask, kNone, kNone,
                             kTRUE, kFALSE);
   } else if (event->fType == kButtonRelease && event->fCode == kButton1) {
      SendMessage(fMsgWindow, MK_MSG(kC_HSLIDER, kSL_RELEASE), fWidgetId, 0);
      fClient->ProcessLine(fCommand, MK_MSG(kC_HSLIDER, kSL_RELEASE), fWidgetId, 0);
      Released();
      fMove = 0;
      gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);  // ungrab pointer
   } else
      fMove = 0;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handles resize events for this widget.

Bool_t TGTripleHSlider::HandleConfigureNotify(Event_t* event)
{
   TGFrame::HandleConfigureNotify(event);
   SetPointerPosition(fSCz);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse motion event in horizontal slide widget.

Bool_t TGTripleHSlider::HandleMotion(Event_t *event)
{
   if (fMove < 3) {
      // if the mouse pointer is on the cursor,
      // and we are not moving anything,
      // set the cursor shape as Pointer
      if (event->fX > (fCz - 5) && event->fX < (fCz + 5) &&
          event->fY > ((Int_t)fHeight / 2) - 7 &&
          event->fY < ((Int_t)fHeight / 2) + 5 &&
          fMove == 0)
         gVirtualX->SetCursor(fId, kNone);
      else
         ChangeCursor(event);
   }
   static int oldDiff = 0;
   static Long64_t was = gSystem->Now();
   Long64_t now = gSystem->Now();

   if (fMove == 0)  return kTRUE;
   if ((now-was) < 50) return kTRUE;
   was = now;

   int     diff;
   Double_t oldMin, oldMax;

   diff    = event->fX - fPressPoint;
   oldMin  = fSmin;
   oldMax  = fSmax;

   if (fMove == 1) {
      // change of min value
      oldDiff = 0;
      fSmin = fPressSmin + diff * (fVmax - fVmin) / (fWidth-16);
      if (fSmin < fVmin) fSmin = fVmin;
      if (fSmin > fSmax) fSmin = fSmax;
   } else if (fMove == 2) {
      // change of max value
      oldDiff = 0;
      fSmax = fPressSmax + diff * (fVmax - fVmin) / (fWidth-16);
      if (fSmax > fVmax) fSmax = fVmax;
      if (fSmax < fSmin) fSmax = fSmin;
   } else if (fMove == 3) {
      // change of min and of max value
      Double_t logicalDiff;
      logicalDiff = diff * (fVmax - fVmin) / (fWidth-16);
      if (fPressSmax + logicalDiff > fVmax)
         logicalDiff = fVmax - fPressSmax;
      if (fPressSmin + logicalDiff < fVmin)
         logicalDiff = fVmin - fPressSmin;
      fSmax = fPressSmax + logicalDiff;
      fSmin = fPressSmin + logicalDiff;
      if (fRelative) {
         if (abs(diff) < 3) oldDiff = diff;
         SetPointerPos(diff - oldDiff, 3);
         oldDiff = diff;
      }
   }
   else if (fMove == 4) {
      // change pointer position
      oldDiff = 0;
      SetPointerPos(event->fX, 1);
   }
   if (fMove != 4) {
      SetPointerPos(0, 2);
   }
   // check if position has changed
   if (fMove != 0 && (fSmax != oldMax || fSmin != oldMin)) {
      fClient->NeedRedraw(this);
      SendMessage(fMsgWindow, MK_MSG(kC_HSLIDER, kSL_POS), fWidgetId, 0);
      fClient->ProcessLine(fCommand, MK_MSG(kC_HSLIDER, kSL_POS), fWidgetId, 0);
      PositionChanged();
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set pointer position constrained in the slider range.

void TGTripleHSlider::SetConstrained(Bool_t on)
{
   fConstrained = on;

   if (fConstrained) {
      if (GetPointerPosition() <= GetMinPosition())
         SetPointerPos((Int_t)GetMinPosition(), 3);
      else if (GetPointerPosition() >= GetMaxPosition())
         SetPointerPos((Int_t)GetMaxPosition(), 3);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set slider pointer position in pixel value.

void TGTripleHSlider::SetPointerPos(Int_t z, Int_t opt)
{
   static Long64_t was = gSystem->Now();
   Bool_t lcheck = (opt == 1);
   Int_t oldPos = fCz;

   if (opt < 2) {
      fCz = z;

      if (fCz < 7)
         fCz = 7;
      else if (fCz >= (Int_t)fWidth-7)
         fCz = (Int_t)fWidth-7;
   }
   if (opt == 3) {
      lcheck = kTRUE;
      fCz += z;
      if (fCz < 7)
         fCz = 7;
      else if (fCz >= (Int_t)fWidth-7)
         fCz = (Int_t)fWidth-7;
   }
   if (fConstrained) {
      int relMin = (int)((fWidth-16) * (fSmin - fVmin) / (fVmax - fVmin)) + 1;
      int relMax = (int)((fWidth-16) * (fSmax - fVmin) / (fVmax - fVmin) + 15);
      if(fCz < relMin+7) {
         fCz = relMin+7;
         lcheck = kTRUE;
      }
      if(fCz > relMax-7) {
         fCz = relMax-7;
         lcheck = kTRUE;
      }
   }
   if (lcheck)
      fSCz = fVmin + ((Double_t)(fCz-8) * (fVmax - fVmin) / (Double_t)(fWidth-16));
   if(fSCz < fVmin) fSCz = fVmin;
   if(fSCz > fVmax) fSCz = fVmax;
   if (fConstrained) {
      if(fSCz < fSmin) fSCz = fSmin;
      if(fSCz > fSmax) fSCz = fSmax;
   }

   DrawPointer();
   fClient->NeedRedraw(this);
   if (fCz != oldPos) {
      Long64_t now = gSystem->Now();
      if ((fMove != 4) && ((now-was) < 150)) return;
      was = now;
      SendMessage(fMsgWindow, MK_MSG(kC_HSLIDER, kSL_POINTER), fWidgetId, 0);
      fClient->ProcessLine(fCommand, MK_MSG(kC_HSLIDER, kSL_POINTER), fWidgetId, 0);
      PointerPositionChanged();
      fClient->NeedRedraw(this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set pointer position in scaled (real) value

void TGTripleHSlider::SetPointerPosition(Double_t pos)
{
   if (fReversedScale) {
      fSCz = fVmin + fVmax - pos;
   }
   else {
      fSCz = pos;
   }
   Double_t absPos = (fSCz - fVmin) * (fWidth-16) / (fVmax - fVmin);
   SetPointerPos((int)(absPos+5.0), 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Save an horizontal slider as a C++ statement(s) on output stream out.

void TGTripleHSlider::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   SaveUserColor(out, option);

   out <<"   TGTripleHSlider *";
   out << GetName() << " = new TGTripleHSlider(" << fParent->GetName()
       << "," << GetWidth() << ",";
   out << GetSString() << "," << WidgetId() << ",";
   out << GetOptionString() << ",ucolor";
   if (fMarkEnds) {
      if (fReversedScale)
         out << ",kTRUE,kTRUE";
      else
         out << ",kFALSE,kTRUE";
   } else if (fReversedScale) {
      out << ",kTRUE,kFALSE";
   } else {
      out << ",kFALSE,kFALSE";
   }
   if (!fConstrained) {
      if (fRelative)
         out << ",kFALSE,kTRUE);" << std::endl;
      else
         out << ",kFALSE,kFALSE);" << std::endl;
   }
   else if (fRelative) {
      out << ",kTRUE);" << std::endl;
   }
   else {
      out << ");" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   if (fVmin != 0 || fVmax != (Int_t)fWidth)
      out << "   " << GetName() << "->SetRange(" << fVmin << "," << fVmax
          << ");" << std::endl;

   if (fSmin != fWidth/8*3 || fSmax != fWidth/8*5)
      out << "   " << GetName() << "->SetPosition(" << GetMinPosition()
          << "," << GetMaxPosition() << ");" << std::endl;

   if (fScale != 10)
      out << "   " << GetName() << "->SetScale(" << fScale << ");" << std::endl;

   out << "   " << GetName() << "->SetPointerPosition(" << fSCz << ");" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Save an horizontal slider as a C++ statement(s) on output stream out.

void TGTripleVSlider::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   SaveUserColor(out, option);

   out<<"   TGTripleVSlider *";
   out << GetName() << " = new TGTripleVSlider("<< fParent->GetName()
       << "," << GetHeight() << ",";
   out << GetSString() << "," << WidgetId() << ",";
   out << GetOptionString() << ",ucolor";
   if (fMarkEnds) {
      if (fReversedScale)
         out << ",kTRUE,kTRUE";
      else
         out << ",kFALSE,kTRUE";
   } else if (fReversedScale) {
      out << ",kTRUE,kFALSE";
   } else {
      out << ",kFALSE,kFALSE";
   }
   if (!fConstrained) {
      if (fRelative)
         out << ",kFALSE,kTRUE);" << std::endl;
      else
         out << ",kFALSE,kFALSE);" << std::endl;
   }
   else if (fRelative) {
      out << ",kTRUE);" << std::endl;
   }
   else {
      out << ");" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   if (fVmin != 0 || fVmax != (Int_t)fHeight)
      out << "   " << GetName() <<"->SetRange(" << fVmin << "," << fVmax
          << ");" << std::endl;

   if (fSmin != fHeight/8*3 || fSmax != fHeight/8*5)
      out << "   " << GetName() << "->SetPosition(" << GetMinPosition()
          << "," << GetMaxPosition() << ");" << std::endl;

   if (fScale != 10)
      out << "   " << GetName() << "->SetScale(" << fScale << ");" << std::endl;

   out << "   " << GetName() << "->SetPointerPosition(" << fSCz << ");" << std::endl;
}
