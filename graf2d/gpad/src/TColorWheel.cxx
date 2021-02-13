// @(#)root/gpad:$Id$
// Author: Rene Brun   10/03/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"
#include "TColorWheel.h"
#include "TCanvas.h"
#include "TText.h"
#include "TGraph.h"
#include "TArc.h"
#include "TLine.h"
#include "TColor.h"
#include "TMath.h"
#include "snprintf.h"

ClassImp(TColorWheel);

/** \class TColorWheel
\ingroup gpad

Draw the ROOT Color Wheel.

The wheel contains the recommended 216 colors to be used in web applications.
The colors in the Color Wheel are created by TColor::CreateColorWheel.

Using this color set for your text, background or graphics will give your
application a consistent appearance across different platforms and browsers.

Colors are grouped by hue, the aspect most important in human perception.
Touching color chips have the same hue, but with different brightness and vividness.

Colors of slightly different hues __clash__. If you intend to display
colors of the same hue together, you should pick them from the same group.

Each color chip is identified by a mnemonic (e.g. kYellow) and a number.
The keywords, kRed, kBlue, kYellow, kPink, etc are defined in the header file __Rtypes.h__
that is included in all ROOT other header files. We strongly recommend to use these keywords
in your code instead of hardcoded color numbers, e.g.:
~~~ {.cpp}
   myObject.SetFillColor(kRed);
   myObject.SetFillColor(kYellow-10);
   myLine.SetLineColor(kMagenta+2);
~~~

Begin_Macro
{
   TColorWheel *w = new TColorWheel();
   cw = new TCanvas("cw","cw",0,0,400,400);
   w->SetCanvas(cw);
   w->Draw();
}
End_Macro
*/

////////////////////////////////////////////////////////////////////////////////
/// constructor

TColorWheel::TColorWheel() :TNamed("wheel","ROOT Color Wheel")
{
   fCanvas = 0;
   fArc    = 0;
   fLine   = 0;
   fText   = 0;
   fGraph  = 0;
   fRmin   = 2.1;
   fRmax   = 9.5;
   fR0     = 4;
   fDr     = 1;
   fRgray  = 1.8;
   fX[ 0]  = fR0-0.2*fDr; fY[ 0] =  0;
   fX[ 1]  = fR0+fDr;     fY[ 1] =  0.75*fDr;
   fX[ 2]  = fR0+fDr;     fY[ 2] = -0.75*fDr;
   fX[ 3]  = fR0+2.2*fDr; fY[ 3] =  1.5*fDr;
   fX[ 4]  = fR0+2.2*fDr; fY[ 4] =  0;
   fX[ 5]  = fR0+2.2*fDr; fY[ 5] = -1.5*fDr;
   fX[ 6]  = fR0+3.4*fDr; fY[ 6] =  2.2*fDr;
   fX[ 7]  = fR0+3.4*fDr; fY[ 7] =  0.7*fDr;
   fX[ 8]  = fR0+3.4*fDr; fY[ 8] = -0.7*fDr;
   fX[ 9]  = fR0+3.4*fDr; fY[ 9] = -2.2*fDr;
   fX[10]  = fR0+4.6*fDr; fY[10] =  2.8*fDr;
   fX[11]  = fR0+4.6*fDr; fY[11] =  1.4*fDr;
   fX[12]  = fR0+4.6*fDr; fY[12] =  0;
   fX[13]  = fR0+4.6*fDr; fY[13] = -1.4*fDr;
   fX[14]  = fR0+4.6*fDr; fY[14] = -2.8*fDr;

   SetBit(kCanDelete);
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TColorWheel::~TColorWheel()
{
   //delete fCanvas;  please don't do that
   delete fArc;
   delete fLine;
   delete fText;
   delete fGraph;
}

////////////////////////////////////////////////////////////////////////////////
/// always return the color wheel

Int_t TColorWheel::DistancetoPrimitive(Int_t px, Int_t py)
{
   if (px+py < 0) return 1;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint the color wheel

void TColorWheel::Draw(Option_t *option)
{
   if (!fCanvas) {
      fCanvas = new TCanvas("wheel","ROOT Color Wheel",10,10,400,400);
      fCanvas->ToggleEventStatus();
   }
   fCanvas->Range(-10.5,-10.5,10.5,10.5);
   fCanvas->SetBorderMode(0);
   fCanvas->SetFillColor(TColor::GetColor(243,241,174));
   AppendPad(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the color number pointed by the mouse

Int_t TColorWheel::GetColor(Int_t px, Int_t py) const
{
   Double_t x = fCanvas->AbsPixeltoX(px);
   Double_t y = fCanvas->AbsPixeltoY(py);
   Int_t n = 0;

   n = InGray(x,y);                   if (n >= 0) return n;
   n = InCircles   (x,y,kMagenta, 0); if (n >= 0) return n;
   n = InRectangles(x,y,kPink,   30); if (n >= 0) return n;
   n = InCircles   (x,y,kRed,    60); if (n >= 0) return n;
   n = InRectangles(x,y,kOrange, 90); if (n >= 0) return n;
   n = InCircles   (x,y,kYellow,120); if (n >= 0) return n;
   n = InRectangles(x,y,kSpring,150); if (n >= 0) return n;
   n = InCircles   (x,y,kGreen, 180); if (n >= 0) return n;
   n = InRectangles(x,y,kTeal,  210); if (n >= 0) return n;
   n = InCircles   (x,y,kCyan,  240); if (n >= 0) return n;
   n = InRectangles(x,y,kAzure, 270); if (n >= 0) return n;
   n = InCircles   (x,y,kBlue,  300); if (n >= 0) return n;
   n = InRectangles(x,y,kViolet,330); if (n >= 0) return n;
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the color number pointed by the mouse

char  *TColorWheel::GetObjectInfo(Int_t px, Int_t py) const
{
   static char info[50];
   info[0] = 0;

   Int_t n = GetColor(px,py);
   if (n < 0) return info;
   TColor *color = gROOT->GetColor(n);
   if (!color) return info;
   Int_t r = (Int_t)(255.01*color->GetRed());
   Int_t g = (Int_t)(255.01*color->GetGreen());
   Int_t b = (Int_t)(255.01*color->GetBlue());
   int res = snprintf(info,sizeof(info),"col %d, %s, r=%3d, g=%3d, b=%3d",n,color->GetName(),r,g,b);
   // check improbable error condition, suppress gcc9 warnings
   if ((res < 0) || (res >= (int) sizeof(info)))
      info[0] = 0;
   return info;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the color number when the mouse point to a circle

Int_t TColorWheel::InCircles(Double_t x, Double_t y, Int_t coffset, Double_t angle) const
{
   Double_t ang = angle*TMath::DegToRad();
   Double_t u,v;
   Rotate(x,y,u,v,ang);
   Double_t r2 = 0.7*0.7*fDr*fDr;
   for (Int_t i=0;i<15;i++) {
      Double_t dx = u-fX[i];
      Double_t dy = v-fY[i];
      if (dx*dx+dy*dy < r2) return coffset+i-10;
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the color number when the mouse point to the gray circle

Int_t TColorWheel::InGray(Double_t x, Double_t y) const
{
   if (x*x+y*y > fRgray*fRgray) return -1;
   Double_t ang = TMath::ATan2(y,x)*TMath::RadToDeg();
   if (ang < 0) ang += 360;
   if (ang <  60) return kWhite;
   if (ang < 120) return kGray;
   if (ang < 180) return kGray+1;
   if (ang < 240) return kGray+2;
   if (ang < 300) return kGray+3;
   return kBlack;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the color number when the mouse point to a rectangle

Int_t TColorWheel::InRectangles(Double_t x, Double_t y, Int_t coffset, Double_t angle) const
{
   Double_t ang = angle*TMath::DegToRad();
   Double_t u,v;
   Rotate(x,y,u,v,ang);
   if (TMath::Abs(v) > 1) return -1;
   if (u < fRmin || u > fRmax) return -1;
   Int_t div = (Int_t)(10*(u-fRmin)/(fRmax-fRmin));
   if (v > 0) return coffset + div+1;
   return coffset+div-9;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint the color wheel

void TColorWheel::Paint(Option_t * /*option*/)
{
   if (!fArc) {
      fArc   = new TArc;
      fLine  = new TLine;
      fText  = new TText();
      fGraph = new TGraph();
   }
   PaintGray();
   PaintCircles   (kMagenta, 0);
   PaintRectangles(kPink,   30);
   PaintCircles   (kRed,    60);
   PaintRectangles(kOrange, 90);
   PaintCircles   (kYellow,120);
   PaintRectangles(kSpring,150);
   PaintCircles   (kGreen, 180);
   PaintRectangles(kTeal,  210);
   PaintCircles   (kCyan,  240);
   PaintRectangles(kAzure, 270);
   PaintCircles   (kBlue,  300);
   PaintRectangles(kViolet,330);

   fText->SetTextFont(72);
   fText->SetTextColor(kBlue);
   fText->SetTextAlign(11);
   fText->SetTextSize(0.03);
   fText->SetTextAngle(0);
   fText->PaintText(-10.2,-10.2,"ROOT Color Wheel");
}

////////////////////////////////////////////////////////////////////////////////
/// Draw one color of type circle

void TColorWheel::PaintCircle(Int_t coffset,Int_t n, Double_t x, Double_t y, Double_t ang) const
{
   Double_t u,v;
   Rotate(x,y,u,v,ang);
   Int_t colorn = coffset+n;
   TColor *color = gROOT->GetColor(colorn);
   if (!color) return;
   fArc->SetFillColor(colorn);
   fArc->SetLineColor(14);
   Double_t r = 0.7*fDr;
   fArc->PaintEllipse(u,v,r,r,0,360,0);
   fText->SetTextSize(0.03);
   fText->SetTextAlign(22);
   if (255*color->GetLight() <150 && n != 0) fText->SetTextColor(0);
   if (n>0) fText->PaintText(u,v,Form("+%d",n));
   else     fText->PaintText(u,v,Form("%d", n));
}

////////////////////////////////////////////////////////////////////////////////
/// Draw all colors of type circle

void TColorWheel::PaintCircles(Int_t coffset, Double_t angle) const
{
   Double_t ang = TMath::DegToRad()*angle;
   Double_t u,v,u0,v0;
   Rotate(fR0+4.6*fDr,2.8*fDr,u0,v0,ang);
   Rotate(fR0+5.8*fDr,2.1*fDr,u,v,ang);
   fLine->PaintLine(u,v,u0,v0);
   fText->SetTextAlign(22);
   fText->SetTextFont(72);
   fText->SetTextColor(1);
   fText->SetTextSize(0.03);
   Double_t tangle = angle-90;
   if (angle == 240) tangle = -30;
   if (angle == 300) tangle =  30;
   TColor *col = gROOT->GetColor(coffset);
   if (!col) return;
   fText->SetTextAngle(tangle);
   fText->PaintText(u,v,col->GetName());

   for (Int_t i=0;i<15;i++) {
      PaintCircle(coffset,i-10, fX[i], fY[i],    ang);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw all colors of type rectangle

void TColorWheel::PaintRectangles(Int_t coffset, Double_t angle) const
{
   Double_t ang = TMath::DegToRad()*angle;
   Double_t rmin = fRmin, rmax=fRmax;
   Double_t dr = (rmax-rmin)/10;
   Double_t dy = -1.0;

   Double_t u,v,u0,v0;
   Rotate(rmax+0.62*dr,0,u,v,ang);
   Rotate(rmax-dr,0.9*dy,u0,v0,ang);
   fLine->PaintLine(u,v,u0,v0);
   fText->SetTextAlign(22);
   fText->SetTextFont(72);
   fText->SetTextColor(1);
   fText->SetTextSize(0.03);
   Double_t tangle = angle+90;
   if (angle == 30)  tangle = -60;
   if (angle == 90)  tangle = 0;
   if (angle == 150) tangle = 60;
   if (angle == 210) tangle = -60;
   if (angle == 270) tangle = 0;
   fText->SetTextAngle(tangle);
   TColor *color = gROOT->GetColor(coffset);
   if (!color) return;
   fText->PaintText(u,v,color->GetName());

   Double_t x[5],y[5];
   Int_t n=-10;
   for (Int_t j=0;j<2;j++) {
      for (Int_t i=0;i<10;i++) {
         n++;
         Int_t colorn = coffset +n;
         color = gROOT->GetColor(colorn);
         Rotate(rmin+i*dr,    0,x[0],y[0],ang);
         Rotate(rmin+i*dr,   dy,x[1],y[1],ang);
         Rotate(rmin+i*dr+dr,dy,x[2],y[2],ang);
         Rotate(rmin+i*dr+dr, 0,x[3],y[3],ang);
         fGraph->SetFillColor(colorn);
         fGraph->PaintGraph(4,x,y,"f");
         Rotate(rmin+i*dr+0.5*dr,0.5*dy,x[0],y[0],ang);
         fText->SetTextSize(0.03);
         fText->SetTextAlign(22);
         if (color) {
            if (255*color->GetLight() <110) fText->SetTextColor(0);
         }
         Double_t tang = angle-90;
         if (angle > 180) tang -=180;
         fText->SetTextAngle(tang);
         if (n > 0) fText->PaintText(x[0],y[0],Form("+%d",n));
         else       fText->PaintText(x[0],y[0],Form("%d",n));
      }
      dy=1;
   }

   Rotate(rmin,-dy,x[0],y[0],ang);
   Rotate(rmax,-dy,x[1],y[1],ang);
   Rotate(rmax, dy,x[2],y[2],ang);
   Rotate(rmin, dy,x[3],y[3],ang);
   Rotate(rmin,-dy,x[4],y[4],ang);
   fGraph->SetLineColor(1);
   fGraph->SetLineWidth(1);
   fGraph->PaintGraph(5,x,y,"l");
   fLine->SetLineWidth(1);
   Rotate(rmin+3*dr,-dy,x[0],y[0],ang);
   Rotate(rmin+3*dr, dy,x[1],y[1],ang);
   fLine->PaintLine(x[0],y[0],x[1],y[1]);
   Rotate(rmin+6*dr,-dy,x[0],y[0],ang);
   Rotate(rmin+6*dr, dy,x[1],y[1],ang);
   fLine->PaintLine(x[0],y[0],x[1],y[1]);
   Rotate(rmin+9*dr,-dy,x[0],y[0],ang);
   Rotate(rmin+9*dr, dy,x[1],y[1],ang);
   fLine->PaintLine(x[0],y[0],x[1],y[1]);
   Rotate(rmin+7*dr,-dy,x[0],y[0],ang);
   Rotate(rmin+7*dr, dy,x[1],y[1],ang);
   fLine->PaintLine(x[0],y[0],x[1],y[1]);
   Rotate(rmin+6*dr,0,x[0],y[0],ang);
   Rotate(rmax, 0,x[1],y[1],ang);
   fLine->PaintLine(x[0],y[0],x[1],y[1]);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the gray colors + white + black

void TColorWheel::PaintGray() const
{
   Double_t r = fRgray;
   fArc->SetFillColor(kWhite);
   fArc->PaintEllipse(0,0,r,r,0,60,0);
   fArc->SetFillColor(kGray);
   fArc->PaintEllipse(0,0,r,r,60,120,0);
   fArc->SetFillColor(kGray+1);
   fArc->PaintEllipse(0,0,r,r,120,180,0);
   fArc->SetFillColor(kGray+2);
   fArc->PaintEllipse(0,0,r,r,180,240,0);
   fArc->SetFillColor(kGray+3);
   fArc->PaintEllipse(0,0,r,r,240,300,0);
   fArc->SetFillColor(kBlack);
   fArc->PaintEllipse(0,0,r,r,300,360,0);

   fText->SetTextAlign(22);
   fText->SetTextFont(62);
   fText->SetTextColor(1);
   fText->SetTextSize(0.02);
   fText->SetTextAngle(40);
   fText->PaintText(0.5*r,0.3*r,"kWhite");
   fText->SetTextAngle(0);
   fText->PaintText(0,0.8*r,"kGray");
   fText->SetTextColor(10);
   fText->SetTextFont(72);
   fText->SetTextSize(0.03);
   fText->PaintText(-0.6*r, 0.3*r,"+1");
   fText->PaintText(-0.6*r,-0.3*r,"+2");
   fText->PaintText(0,-0.6*r,"+3");
   fText->SetTextAngle(-40);
   fText->SetTextSize(0.02);
   fText->SetTextFont(62);
   fText->PaintText(0.5*r,-0.35*r,"kBlack");
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate point x,y with an angle=ang

void TColorWheel::Rotate(Double_t x, Double_t y, Double_t &u, Double_t &v, Double_t ang) const
{
   u = x*TMath::Cos(ang)  + y*TMath::Sin(ang);
   v = x*TMath::Sin(ang)  - y*TMath::Cos(ang);
}
