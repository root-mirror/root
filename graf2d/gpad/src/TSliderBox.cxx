// @(#)root/gpad:$Id$
// Author: Rene Brun   23/11/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"
#include "TSlider.h"
#include "TSliderBox.h"
#include "TVirtualX.h"

#include <cstring>

ClassImp(TSliderBox);


/** \class TSliderBox
\ingroup gpad

The moving box in a TSlider
*/

////////////////////////////////////////////////////////////////////////////////
/// SliderBox default constructor.

TSliderBox::TSliderBox(): TWbox()
{
   fSlider = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// SliderBox normal constructor.

TSliderBox::TSliderBox(Double_t x1, Double_t y1,Double_t x2, Double_t  y2, Color_t color, Short_t bordersize, Short_t bordermode)
           :TWbox(x1,y1,x2,y2,color,bordersize,bordermode)
{
   fSlider = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// SliderBox default destructor.

TSliderBox::~TSliderBox()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Interaction with a slider.

void TSliderBox::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (!gPad) return;

   const Int_t kMaxDiff = 5;
   const Int_t kMinSize = 20;

   static Int_t px1, px2, py1, py2, pxl, pyl, pxt, pyt, pxold, pyold;
   static Int_t px1p, px2p, py1p, py2p;
   static Bool_t  pL, pR, pTop, pBot, pINSIDE;
   Int_t  wx, wy;
   Bool_t doing_again = kFALSE;
   Bool_t opaque  = gPad->OpaqueMoving();
   Bool_t ropaque = gPad->OpaqueResizing();

   TVirtualPad  *parent = gPad;

   Int_t border   = parent->GetBorderSize();
   Int_t pxpadmin = parent->XtoAbsPixel(parent->GetX1()) + border;
   Int_t pxpadmax = parent->XtoAbsPixel(parent->GetX2()) - border;
   Int_t pypadmin = parent->YtoAbsPixel(parent->GetY1()) - border;
   Int_t pypadmax = parent->YtoAbsPixel(parent->GetY2()) + border;

   Bool_t vertical = kTRUE;
   if (pxpadmax-pxpadmin > pypadmin-pypadmax) vertical = kFALSE;

again:

   switch (event) {

   case kButton1Down:

      gVirtualX->SetLineColor(-1);
      TAttLine::Modify();  //Change line attributes only if necessary
      if (GetFillColor())
         gVirtualX->SetLineColor(GetFillColor());
      else
         gVirtualX->SetLineColor(1);
      gVirtualX->SetLineWidth(2);

      // No break !!!

   case kMouseMotion:

      px1 = gPad->XtoAbsPixel(GetX1());
      py1 = gPad->YtoAbsPixel(GetY1());
      px2 = gPad->XtoAbsPixel(GetX2());
      py2 = gPad->YtoAbsPixel(GetY2());

      if (px1 < px2) {
         pxl = px1;
         pxt = px2;
      } else {
         pxl = px2;
         pxt = px1;
      }
      if (py1 < py2) {
         pyl = py1;
         pyt = py2;
      } else {
         pyl = py2;
         pyt = py1;
      }

      px1p = pxpadmin;
      py1p = pypadmin;
      px2p = pxpadmax;
      py2p = pypadmax;

      pL = pR  = pTop = pBot = pINSIDE = kFALSE;

      if (vertical && (px > pxl+kMaxDiff && px < pxt-kMaxDiff) &&
          TMath::Abs(py - pyl) < kMaxDiff) {             // top edge
         pxold = pxl; pyold = pyl; pTop = kTRUE;
         gPad->SetCursor(kTopSide);
      }

      if (vertical && (px > pxl+kMaxDiff && px < pxt-kMaxDiff) &&
          TMath::Abs(py - pyt) < kMaxDiff) {             // bottom edge
         pxold = pxt; pyold = pyt; pBot = kTRUE;
         gPad->SetCursor(kBottomSide);
      }

      if (!vertical && (py > pyl+kMaxDiff && py < pyt-kMaxDiff) &&
          TMath::Abs(px - pxl) < kMaxDiff) {             // left edge
         pxold = pxl; pyold = pyl; pL = kTRUE;
         gPad->SetCursor(kLeftSide);
      }

      if (!vertical && (py > pyl+kMaxDiff && py < pyt-kMaxDiff) &&
          TMath::Abs(px - pxt) < kMaxDiff) {             // right edge
         pxold = pxt; pyold = pyt; pR = kTRUE;
         gPad->SetCursor(kRightSide);
      }

      if ((px > pxl+kMaxDiff && px < pxt-kMaxDiff) &&
          (py > pyl+kMaxDiff && py < pyt-kMaxDiff)) {    // inside box
         pxold = px; pyold = py; pINSIDE = kTRUE;
         if (event == kButton1Down)
            gPad->SetCursor(kMove);
         else
            gPad->SetCursor(kCross);
      }

      fResizing = kFALSE;
      if ( pL || pR || pTop || pBot)
         fResizing = kTRUE;

      if ( !pL && !pR && !pTop && !pBot && !pINSIDE)
         gPad->SetCursor(kCross);

      break;

   case kButton1Motion:

      wx = wy = 0;

      if (pTop) {
         if (!ropaque) gVirtualX->DrawBox(px1, py1, px2, py2, TVirtualX::kHollow);
         py2 += py - pyold;
         if (py2 > py1-kMinSize) { py2 = py1-kMinSize; wy = py2; }
         if (py2 < py2p) { py2 = py2p; wy = py2; }
         if (!ropaque) gVirtualX->DrawBox(px1, py1, px2, py2, TVirtualX::kHollow);
      }
      if (pBot) {
         if (!ropaque) gVirtualX->DrawBox(px1, py1, px2, py2, TVirtualX::kHollow);
         py1 += py - pyold;
         if (py1 < py2+kMinSize) { py1 = py2+kMinSize; wy = py1; }
         if (py1 > py1p) { py1 = py1p; wy = py1; }
         if (!ropaque) gVirtualX->DrawBox(px1, py1, px2, py2, TVirtualX::kHollow);
      }
      if (pL) {
         if (!ropaque) gVirtualX->DrawBox(px1, py1, px2, py2, TVirtualX::kHollow);
         px1 += px - pxold;
         if (px1 > px2-kMinSize) { px1 = px2-kMinSize; wx = px1; }
         if (px1 < px1p) { px1 = px1p; wx = px1; }
         if (!ropaque) gVirtualX->DrawBox(px1, py1, px2, py2, TVirtualX::kHollow);
      }
      if (pR) {
         if (!ropaque) gVirtualX->DrawBox(px1, py1, px2, py2, TVirtualX::kHollow);
         px2 += px - pxold;
         if (px2 < px1+kMinSize) { px2 = px1+kMinSize; wx = px2; }
         if (px2 > px2p) { px2 = px2p; wx = px2; }
         if (!ropaque) gVirtualX->DrawBox(px1, py1, px2, py2, TVirtualX::kHollow);
      }
      if (pINSIDE) {
         if (!opaque) gVirtualX->DrawBox(px1, py1, px2, py2, TVirtualX::kHollow);  // draw the old box
         Int_t dx = px - pxold;
         Int_t dy = py - pyold;
         px1 += dx; py1 += dy; px2 += dx; py2 += dy;
         if (px1 < px1p) { dx = px1p - px1; px1 += dx; px2 += dx; wx = px+dx; }
         if (px2 > px2p) { dx = px2 - px2p; px1 -= dx; px2 -= dx; wx = px-dx; }
         if (py1 > py1p) { dy = py1 - py1p; py1 -= dy; py2 -= dy; wy = py-dy; }
         if (py2 < py2p) { dy = py2p - py2; py1 += dy; py2 += dy; wy = py+dy; }
         if (!opaque) gVirtualX->DrawBox(px1, py1, px2, py2, TVirtualX::kHollow);  // draw the new box
      }

      if (wx || wy) {
         if (wx) px = wx;
         if (wy) py = wy;
         gVirtualX->Warp(px, py);
      }

      pxold = px;
      pyold = py;

      if ((pINSIDE && opaque) || (fResizing && ropaque)) {
         event = kButton1Up;
         doing_again = kTRUE;
         goto again;
      }

      break;

   case kButton1Up:

      if (pTop || pBot || pL || pR || pINSIDE) {
         fX1 = gPad->AbsPixeltoX(px1);
         fY1 = gPad->AbsPixeltoY(py1);
         fX2 = gPad->AbsPixeltoX(px2);
         fY2 = gPad->AbsPixeltoY(py2);
      }

      if (pINSIDE) {
         // if it was not a pad that was moved then it must have been
         // a box or something like that so we have to redraw the pad
         if (parent == gPad) gPad->Modified(kTRUE);
         if (!doing_again) gPad->SetCursor(kCross);
      }

      if (pTop || pBot ||  pL || pR )
         gPad->Modified(kTRUE);

      // Restore original event type
      if (doing_again)
         event = kButton1Motion;
      else {
         gVirtualX->SetLineColor(-1);
         gVirtualX->SetLineWidth(-1);
      }

      break;
   }


   // Give control to object using the slider

   Float_t xpmin,xpmax;
   if (vertical) {
      xpmin = Float_t(pypadmin-py1)/Float_t(pypadmin-pypadmax);
      xpmax = Float_t(pypadmin-py2)/Float_t(pypadmin-pypadmax);
   } else {  //vertical slider
      xpmin = Float_t(px1-pxpadmin)/Float_t(pxpadmax-pxpadmin);
      xpmax = Float_t(px2-pxpadmin)/Float_t(pxpadmax-pxpadmin);
   }
   fSlider->SetMinimum(xpmin);
   fSlider->SetMaximum(xpmax);

   //A user method to execute?
   Int_t lenMethod = strlen(fSlider->GetMethod());
   if (event == kButton1Up && lenMethod > 0 ) {
      gPad->SetCursor(kWatch);
      gROOT->ProcessLine(fSlider->GetMethod());
      return;
   }

   //An object connected to this slider?
   TObject *obj = fSlider->GetObject();
   if (obj) {
      obj->ExecuteEvent(event,0,0);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out

void TSliderBox::SavePrimitive(std::ostream &, Option_t * /*= ""*/)
{
}
