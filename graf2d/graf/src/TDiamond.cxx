// @(#)root/graf:$Id$
// Author: Rene Brun   22/06/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <stdlib.h>

#include <iostream>
#include "TBufferFile.h"
#include "TROOT.h"
#include "TDiamond.h"
#include "TVirtualPad.h"
#include "TVirtualX.h"
#include "TMath.h"

ClassImp(TDiamond);

/** \class TDiamond
\ingroup BasicGraphics

Draw a Diamond.

A diamond is defined by:

- Its central left coordinates x1,y1
- Its top central coordinates x2,y2

A diamond has line attributes (see TAttLine) and fill area attributes (see TAttFill).

Like for the class TPaveText, a TDiamond may have one or more line(s) of text inside.

Begin_Macro(source)
../../../tutorials/graphics/diamond.C
End_Macro
*/

////////////////////////////////////////////////////////////////////////////////
/// Diamond default constructor.

TDiamond::TDiamond(): TPaveText()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Diamond standard constructor.

TDiamond::TDiamond(Double_t x1, Double_t y1,Double_t x2, Double_t  y2)
     :TPaveText(x1,y1,x2,y2)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Diamond destructor.

TDiamond::~TDiamond()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TDiamond::TDiamond(const TDiamond &diamond) : TPaveText()
{
   TBufferFile b(TBuffer::kWrite);
   TDiamond *p = (TDiamond*)(&diamond);
   p->Streamer(b);
   b.SetReadMode();
   b.SetBufferOffset(0);
   Streamer(b);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from point px,py to a diamond.
///
///  Compute the closest distance of approach from point px,py to the
///  edges of this diamond.
///  The distance is computed in pixels units.

Int_t TDiamond::DistancetoPrimitive(Int_t px, Int_t py)
{
   return TPaveText::DistancetoPrimitive(px,py);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this diamond with its current attributes.

void TDiamond::Draw(Option_t *option)
{
   AppendPad(option);

}

////////////////////////////////////////////////////////////////////////////////
/// Execute action corresponding to one event.
///
/// This member function is called when a Diamond object is clicked.
///
/// If the mouse is clicked inside the diamond, the diamond is moved.
///
/// If the mouse is clicked on the 4 tops (pL,pR,pTop,pBot), the diamond is
/// rescaled.

void TDiamond::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (!gPad) return;

   const Int_t kMaxDiff = 5;
   const Int_t kMinSize = 20;

   static Int_t px1, px2, py1, py2, pxl, pyl, pxt, pyt, pxold, pyold;
   static Int_t px1p, px2p, py1p, py2p;
   static Int_t pTx,pTy,pLx,pLy,pRx,pRy,pBx,pBy;
   static Double_t x1c,x2c,x3c,x4c;
   static Bool_t pTop, pL, pR, pBot, pINSIDE;
   static Int_t i,x[5], y[5];
   Int_t  wx, wy;
   TVirtualPad  *parent;
   Bool_t opaque  = gPad->OpaqueMoving();
   Bool_t ropaque = gPad->OpaqueResizing();

   if (!gPad->IsEditable()) return;

   parent = gPad;

   switch (event) {

   case kArrowKeyPress:
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

      px1p = parent->XtoAbsPixel(parent->GetX1()) + parent->GetBorderSize();
      py1p = parent->YtoAbsPixel(parent->GetY1()) - parent->GetBorderSize();
      px2p = parent->XtoAbsPixel(parent->GetX2()) - parent->GetBorderSize();
      py2p = parent->YtoAbsPixel(parent->GetY2()) + parent->GetBorderSize();

      pTx = pBx = (pxl+pxt)/2;
      pLy = pRy = (pyl+pyt)/2;
      pTy = pyl;
      pBy = pyt;
      pLx = pxl;
      pRx = pxt;

      pTop = pL = pR = pBot = pINSIDE = kFALSE;

      if ((TMath::Abs(px-(pxl+pxt)/2) < kMaxDiff) &&
          (TMath::Abs(py - pyl) < kMaxDiff)) {             // top edge
         pxold = pxl; pyold = pyl; pTop = kTRUE;
         gPad->SetCursor(kTopSide);
      }

      if ((TMath::Abs(px-(pxl+pxt)/2) < kMaxDiff) &&
          (TMath::Abs(py - pyt) < kMaxDiff)) {             // bottom edge
         pxold = pxt; pyold = pyt; pBot = kTRUE;
         gPad->SetCursor(kBottomSide);
      }

      if ((TMath::Abs(py-(pyl+pyt)/2) < kMaxDiff) &&
          (TMath::Abs(px - pxl) < kMaxDiff)) {             // left edge
         pxold = pxl; pyold = pyl; pL = kTRUE;
         gPad->SetCursor(kLeftSide);
      }

      if ((TMath::Abs(py-(pyl+pyt)/2) < kMaxDiff) &&
          (TMath::Abs(px - pxt) < kMaxDiff)) {             // right edge
         pxold = pxt; pyold = pyt; pR = kTRUE;
         gPad->SetCursor(kRightSide);
      }

      x1c = (py-pTy)*(pTx-pLx)/(pTy-pLy)+pTx;
      x2c = (py-pTy)*(pRx-pTx)/(pRy-pTy)+pTx;
      x3c = (py-pRy)*(pRx-pBx)/(pRy-pBy)+pRx;
      x4c = (py-pBy)*(pBx-pLx)/(pBy-pLy)+pBx;

      if (px > x1c+kMaxDiff && px < x2c-kMaxDiff &&
          px > x4c+kMaxDiff && px < x3c-kMaxDiff) {    // inside box
         pxold = px; pyold = py; pINSIDE = kTRUE;
         if (event == kButton1Down)
            gPad->SetCursor(kMove);
         else
            gPad->SetCursor(kCross);
      }

      fResizing = kFALSE;
      if (pTop || pL || pR || pBot)
         fResizing = kTRUE;

      if (!pTop && !pL && !pR && !pBot && !pINSIDE)
         gPad->SetCursor(kCross);

      break;

   case kArrowKeyRelease:
   case kButton1Motion:

      wx = wy = 0;
      x[0] = x[2] = x[4] = (px1+px2)/2;
      x[1] = px2;
      x[3] = px1;
      y[0] = y[4] = py1;
      y[2] = py2;
      y[1] = y[3] = (py1+py2)/2;
      if (pTop) {
         for (i=0;i<4;i++) gVirtualX->DrawLine(x[i], y[i], x[i+1], y[i+1]);
         py2 += py - pyold;
         if (py2 > py1-kMinSize) { py2 = py1-kMinSize; wy = py2; }
         if (py2 < py2p) { py2 = py2p; wy = py2; }
         y[2] = py2;
         y[1] = y[3] = (py1+py2)/2;
         for (i=0;i<4;i++) gVirtualX->DrawLine(x[i], y[i], x[i+1], y[i+1]);
      }
      if (pBot) {
         for (i=0;i<4;i++) gVirtualX->DrawLine(x[i], y[i], x[i+1], y[i+1]);
         py1 += py - pyold;
         if (py1 < py2+kMinSize) { py1 = py2+kMinSize; wy = py1; }
         if (py1 > py1p) { py1 = py1p; wy = py1; }
         y[0] = y[4] = py1;
         y[1] = y[3] = (py1+py2)/2;
         for (i=0;i<4;i++) gVirtualX->DrawLine(x[i], y[i], x[i+1], y[i+1]);
      }
      if (pL) {
         for (i=0;i<4;i++) gVirtualX->DrawLine(x[i], y[i], x[i+1], y[i+1]);
         px1 += px - pxold;
         if (px1 > px2-kMinSize) { px1 = px2-kMinSize; wx = px1; }
         if (px1 < px1p) { px1 = px1p; wx = px1; }
         x[3] = px1;
         x[0] = x[2] = x[4] = (px1+px2)/2;
         for (i=0;i<4;i++) gVirtualX->DrawLine(x[i], y[i], x[i+1], y[i+1]);
      }
      if (pR) {
         for (i=0;i<4;i++) gVirtualX->DrawLine(x[i], y[i], x[i+1], y[i+1]);
         px2 += px - pxold;
         if (px2 < px1+kMinSize) { px2 = px1+kMinSize; wx = px2; }
         if (px2 > px2p) { px2 = px2p; wx = px2; }
         x[1] = px2;
         x[0] = x[2] = x[4] = (px1+px2)/2;
         for (i=0;i<4;i++) gVirtualX->DrawLine(x[i], y[i], x[i+1], y[i+1]);
      }
      if (pINSIDE) {
         for (i=0;i<4;i++) gVirtualX->DrawLine(x[i], y[i], x[i+1], y[i+1]);
         Int_t dx = px - pxold;
         Int_t dy = py - pyold;
         px1 += dx; py1 += dy; px2 += dx; py2 += dy;
         if (px1 < px1p) { dx = px1p - px1; px1 += dx; px2 += dx; wx = px+dx; }
         if (px2 > px2p) { dx = px2 - px2p; px1 -= dx; px2 -= dx; wx = px-dx; }
         if (py1 > py1p) { dy = py1 - py1p; py1 -= dy; py2 -= dy; wy = py-dy; }
         if (py2 < py2p) { dy = py2p - py2; py1 += dy; py2 += dy; wy = py+dy; }
         x[0] = x[2] = x[4] = (px1+px2)/2;
         x[1] = px2;
         x[3] = px1;
         y[0] = y[4] = py1;
         y[2] = py2;
         y[1] = y[3] = (py1+py2)/2;
         for (i=0;i<4;i++) gVirtualX->DrawLine(x[i], y[i], x[i+1], y[i+1]);
      }

      if (wx || wy) {
         if (wx) px = wx;
         if (wy) py = wy;
         gVirtualX->Warp(px, py);
      }

      pxold = px;
      pyold = py;

      if ((pINSIDE && opaque) || (fResizing && ropaque)) {
         if (pTop || pBot || pL || pR) {
            fX1 = gPad->AbsPixeltoX(px1);
            fY1 = gPad->AbsPixeltoY(py1);
            fX2 = gPad->AbsPixeltoX(px2);
            fY2 = gPad->AbsPixeltoY(py2);
         }
         if (pINSIDE) {
            fX1 = gPad->AbsPixeltoX(px1);
            fY1 = gPad->AbsPixeltoY(py1);
            fX2 = gPad->AbsPixeltoX(px2);
            fY2 = gPad->AbsPixeltoY(py2);
            // if it was not a pad that was moved then it must have been
            // a box or something like that so we have to redraw the pad
            if (parent == gPad) gPad->Modified(kTRUE);
         }

         if (pINSIDE) gPad->ShowGuidelines(this, event, 'i', true);
         if (pTop) gPad->ShowGuidelines(this, event, 't', true);
         if (pBot) gPad->ShowGuidelines(this, event, 'b', true);
         if (pL) gPad->ShowGuidelines(this, event, 'l', true);
         if (pR) gPad->ShowGuidelines(this, event, 'r', true);

         if (pTop || pL || pR || pBot)
            gPad->Modified(kTRUE);
      }

      break;

   case kButton1Up:

      if (opaque) {
         gPad->ShowGuidelines(this, event);
      } else {
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
         }
      }

      if (pTop || pL || pR || pBot) gPad->Modified(kTRUE);

      if (!opaque) {
         gVirtualX->SetLineColor(-1);
         gVirtualX->SetLineWidth(-1);
      }

      break;

   case kButton1Locate:

      ExecuteEvent(kButton1Down, px, py);

      while (1) {
         px = py = 0;
         event = gVirtualX->RequestLocator(1, 1, px, py);

         ExecuteEvent(kButton1Motion, px, py);

         if (event != -1) {                     // button is released
            ExecuteEvent(kButton1Up, px, py);
            return;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint this diamond with its current attributes.

void TDiamond::Paint(Option_t *)
{
   Double_t x[7],y[7],depx,depy;
   Double_t x1 = fX1;
   Double_t y1 = fY1;
   Double_t x2 = fX2;
   Double_t y2 = fY2;
   Int_t fillstyle = GetFillStyle();
   Int_t fillcolor = GetFillColor();
   Int_t linecolor = GetLineColor();
   if (fBorderSize) {
      Double_t wy = gPad->PixeltoY(0) - gPad->PixeltoY(fBorderSize);
      Double_t wx = gPad->PixeltoX(fBorderSize) - gPad->PixeltoX(0);
      // Draw the frame top right
      if (y2-y1>x2-x1) {
         depx = wx;
         depy = 0;
         }
      else if (y2-y1<x2-x1) {
         depx = 0;
         depy = -wy;
         }
      else {
         depx = wx;
         depy = -wy;
      }
      x[0] = x[2] = (x1+x2)/2+depx;
      x[1] = x2+depx;
      x[3] = x1+depx;
      y[0] = y2+depy;
      y[2] = y1+depy;
      y[1] = y[3] =(y1+y2)/2+depy;
      x[4] = x[0]; y[4] = y[0];
      SetFillStyle(fillstyle);
      SetFillColor(linecolor);
      TAttFill::Modify();  //Change fill area attributes only if necessary
      gPad->PaintFillArea(4,x,y);
   }
   x[0] = x[2] = (x1+x2)/2;
   x[1] = x2;
   x[3] = x1;
   y[0] = y2;
   y[2] = y1;
   y[1] = y[3] = (y1+y2)/2;
   x[4] = x[0]; y[4] =y[0];
   SetLineColor(linecolor);
   SetFillColor(fillcolor);
   TAttLine::Modify();  //Change line attributes only if necessary
   TAttFill::Modify();  //Change fill area attributes only if necessary
   gPad->PaintFillArea(4,x,y);
   gPad->PaintPolyLine(5,x,y);

   // Paint list of primitives (test,etc)
   PaintPrimitives(kDiamond);
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out.

void TDiamond::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   Bool_t saved = gROOT->ClassSaved(TDiamond::Class());
   if (saved) {
      out<<"   ";
   } else {
      out<<"   TDiamond *";
   }
   out<<"diamond = new TDiamond("<<fX1<<","<<fY1<<","<<fX2<<","<<fY2<<");"<<std::endl;

   SaveFillAttributes(out,"diamond",0,1001);
   SaveLineAttributes(out,"diamond",1,1,1);
   SaveTextAttributes(out,"diamond",11,0,1,62,0.05);

   SaveLines(out,"diamond",saved);
   out<<"   diamond->Draw();"<<std::endl;
}
