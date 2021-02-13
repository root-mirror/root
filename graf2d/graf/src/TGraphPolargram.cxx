// @(#)root/graf:$Id$
// Author: Sebastian Boser, Mathieu Demaret 02/02/06

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGraphPolargram
\ingroup BasicGraphics

To draw polar axis

TGraphPolargram draw the polar axis of the TGraphPolar.

Example:

Begin_Macro(source)
{
   TCanvas * CPol = new TCanvas("CPol","TGraphPolar Examples",500,500);

   Double_t rmin=0;
   Double_t rmax=TMath::Pi()*2;
   Double_t r[1000];
   Double_t theta[1000];

   TF1 * fp1 = new TF1("fplot","cos(x)",rmin,rmax);
   for (Int_t ipt = 0; ipt < 1000; ipt++) {
      r[ipt] = ipt*(rmax-rmin)/1000+rmin;
      theta[ipt] = fp1->Eval(r[ipt]);
   }
   TGraphPolar * grP1 = new TGraphPolar(1000,r,theta);
   grP1->SetTitle("");
   grP1->SetLineColor(2);
   grP1->Draw("AOL");
}
End_Macro
*/

#include "TGraphPolar.h"
#include "TGraphPolargram.h"
#include "TGaxis.h"
#include "THLimitsFinder.h"
#include "TVirtualPad.h"
#include "TLatex.h"
#include "TEllipse.h"
#include "TMath.h"

ClassImp(TGraphPolargram);

////////////////////////////////////////////////////////////////////////////////
/// TGraphPolargram Constructor.

TGraphPolargram::TGraphPolargram(const char* name, Double_t rmin, Double_t rmax,
                                 Double_t tmin, Double_t tmax):
                                 TNamed(name,"Polargram")
{
   Init();
   fNdivRad          = 508;
   fNdivPol          = 508;
   fPolarLabels      = NULL;
   fRwrmax           = rmax;
   fRwrmin           = rmin;
   fRwtmin           = tmin;
   fRwtmax           = tmax;
}

////////////////////////////////////////////////////////////////////////////////
/// Short constructor used in the case of a spider plot.

TGraphPolargram::TGraphPolargram(const char* name):
                                 TNamed(name,"Polargram")
{
   Init();
   fNdivRad     = 0;
   fNdivPol     = 0;
   fPolarLabels = NULL;
   fRwrmax      = 1;
   fRwrmin      = 0;
   fRwtmax      = 0;
   fRwtmin      = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphPolargram destructor.

TGraphPolargram::~TGraphPolargram()
{
   if (fPolarLabels != NULL) delete [] fPolarLabels;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the Polar range.
/// \param[in] tmin   the start number.
/// \param[in] tmax   the end number.

void TGraphPolargram::ChangeRangePolar(Double_t tmin, Double_t tmax)
{
   if (tmin < tmax) {
      fRwtmin = tmin;
      fRwtmax = tmax;
   }
   if (gPad) gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// Everything within the circle belongs to the TGraphPolargram.

Int_t TGraphPolargram::DistancetoPrimitive(Int_t px, Int_t py)
{
   Int_t i;
   Double_t x = gPad->AbsPixeltoX(px);
   Double_t y = gPad->AbsPixeltoY(py);

   // Check if close to a (major) radial line.
   Double_t rad = TMath::Sqrt(x*x+y*y);
   Int_t div    = (Int_t)rad*(fNdivRad%100);
   Double_t dr  = TMath::Min(TMath::Abs(rad-div*1./(fNdivRad%100)),
                             TMath::Abs(rad-(div+1)*1./(fNdivRad%100)));
   Int_t drad   = gPad->XtoPixel(dr)-gPad->XtoPixel(0);

   // Check if close to a (major) Polar line.
   // This is not a proper calculation, but rather fast.
   Int_t dt = kMaxPixel;
   for (i=0; i<(fNdivPol%100); i++) {
      Double_t theta = i*2*TMath::Pi()/(fNdivPol%100);

      // Attention: px,py in pixel units, line given in user coordinates.
      Int_t dthis = DistancetoLine(px,py,0.,0.,TMath::Cos(theta),
                                               TMath::Sin(theta));

      // Fails if we are outside box described by the line.
      // (i.e for all hor/vert lines)
      if (dthis==9999) {

         // Outside -> Get distance to endpoint of line.
         if (rad>1) {
            dthis = (Int_t)TMath::Sqrt(
                    TMath::Power(px-gPad->XtoPixel(TMath::Cos(theta)),2)+
                    TMath::Power(py-gPad->YtoPixel(TMath::Sin(theta)),2));
         } else {

            // Check for horizontal line.
            if (((TMath::Abs(theta-TMath::Pi())<0.1) &&
                ((px-gPad->XtoPixel(0))<0))          ||
                ((TMath::Abs(theta)<0.1)             &&
                ((px-gPad->XtoPixel(0))>0))) {
               dthis = TMath::Abs(py-gPad->YtoPixel(0.));
            }

            //Check for vertical line.
            if (((TMath::Abs(theta-TMath::PiOver2())<0.1)   &&
               ((py-gPad->YtoPixel(0))>0))                 ||
                ((TMath::Abs(theta-3*TMath::PiOver2())<0.1) &&
                                  (py-gPad->YtoPixel(0))<0)) {
               dthis = TMath::Abs(px-gPad->XtoPixel(0.));
            }
            if (dthis==9999) {

               // Inside, but out of box for nonorthogonal line ->
               // get distance to start point.
               dthis = (Int_t)TMath::Sqrt(
                       TMath::Power(px-gPad->XtoPixel(0.),2)+
                       TMath::Power(py-gPad->YtoPixel(0.),2));
            }
         }
      }

      // Take distance to closes line.
      dt = TMath::Min(dthis,dt);
   }
   return TMath::Min(drad, dt);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw Polargram.

void TGraphPolargram::Draw(Option_t* options)
{
   Paint(options);
   AppendPad(options);
}

////////////////////////////////////////////////////////////////////////////////
/// Indicate that there is something to click here.

void TGraphPolargram::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (!gPad) return;

   Int_t kMaxDiff = 20;
   static Int_t d1, d2, d3, px1, py1, px3, py3;
   static Bool_t p1, p2, p3, p4, p5, p6, p7, p8;
   Double_t px2, py2;
   p2 = p3 = p4 = p5 = p6 = p7 = p8 = kFALSE;
   if (!gPad->IsEditable()) return;
   switch (event) {
      case kMouseMotion:
         px1 = gPad->XtoAbsPixel(TMath::Cos(GetAngle()));
         py1 = gPad->YtoAbsPixel(TMath::Sin(GetAngle()));
         d1  = TMath::Abs(px1 - px) + TMath::Abs(py1-py); //simply take sum of pixels differences
         p1  = kFALSE;
         px2 = gPad->XtoAbsPixel(-1);
         py2 = gPad->YtoAbsPixel(1);
         d2  = (Int_t)(TMath::Abs(px2 - px) + TMath::Abs(py2 - py)) ;
         px3 = gPad->XtoAbsPixel(-1);
         py3 = gPad->YtoAbsPixel(-1);
         d3  = TMath::Abs(px3 - px) + TMath::Abs(py3 - py) ; //simply take sum of pixels differences
         // check if point is close to the radial axis
         if (d1 < kMaxDiff) {
            gPad->SetCursor(kMove);
            p1 = kTRUE;
         }
         // check if point is close to the left high axis
         if ( d2 < kMaxDiff) {
            gPad->SetCursor(kHand);
            p7 = kTRUE;
         }
         // check if point is close to the left down axis
         if ( d3 < kMaxDiff) {
            gPad->SetCursor(kHand);
            p8 = kTRUE;
         }
         // check if point is close to a main circle
         if (!p1 && !p7 ) {
            p6 = kTRUE;
            gPad->SetCursor(kHand);
         }
         break;

      case kButton1Down:
         // Record initial coordinates
         //px4 = px;
         //py4 = py;

      case kButton1Motion:
         if (p1) {
            px2 = gPad->AbsPixeltoX(px);
            py2 = gPad->AbsPixeltoY(py);
            if ( px2 < 0 && py2 < 0)  {p2 = kTRUE;};
            if ( px2 < 0 && py2 > 0 ) {p3 = kTRUE;};
            if ( px2 > 0 && py2 > 0 ) {p4 = kTRUE;};
            if ( px2 > 0 && py2 < 0 ) {p5 = kTRUE;};
            px2 = TMath::ACos(TMath::Abs(px2));
            py2 = TMath::ASin(TMath::Abs(py2));
            if (p2) {
               fAxisAngle = TMath::Pi()+(px2+py2)/2;
               p2 = kFALSE;
            };
            if (p3) {
               fAxisAngle = TMath::Pi()-(px2+py2)/2;
               p3 = kFALSE;
            };
            if (p4) {
               fAxisAngle = (px2+py2)/2;
               p4 = kFALSE;
            };
            if (p5) {
               fAxisAngle = -(px2+py2)/2;
               p5 = kFALSE;
            };
         }
         break;

      case kButton1Up:
         Paint();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Find the alignement rule to apply for TText::SetTextAlign(Short_t).

Int_t TGraphPolargram::FindAlign(Double_t angle)
{
   Double_t pi = TMath::Pi();

   while(angle < 0 || angle > 2*pi){
      if(angle < 0) angle+=2*pi;
      if(angle > 2*pi) angle-=2*pi;
   }
   if(!TestBit(TGraphPolargram::kLabelOrtho)){
      if(angle > 0 && angle < pi/2) return 11;
      else if(angle > pi/2 && angle < pi) return 31;
      else if(angle > pi && angle < 3*pi/2) return 33;
      else if(angle > 3*pi/2 && angle < 2*pi) return 13;
      else if(angle == 0 || angle == 2*pi) return 12;
      else if(angle == pi/2) return 21;
      else if(angle == pi) return 32;
      else if(angle == 3*pi/2) return 23;
      else return 0;
   }
   else{
      if(angle >= 0 && angle <= pi/2) return 12;
      else if((angle > pi/2 && angle <= pi) || (angle > pi && angle <= 3*pi/2)) return 32;
      else if(angle > 3*pi/2 && angle <= 2*pi) return 12;
      else return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Determine the orientation of the polar labels according to their angle.

Double_t TGraphPolargram::FindTextAngle(Double_t angle)
{
   Double_t pi = TMath::Pi();
   Double_t convraddeg = 180.0/pi;

   while(angle < 0 || angle > 2*pi){
      if(angle < 0) angle+=2*pi;
      if(angle > 2*pi) angle-=2*pi;
   }

   if(angle >= 0 && angle <= pi/2) return angle*convraddeg;
   else if(angle > pi/2 && angle <= pi) return (angle + pi)*convraddeg;
   else if(angle > pi && angle <= 3*pi/2) return (angle - pi)*convraddeg;
   else if(angle > 3*pi/2 && angle <= 2*pi) return angle*convraddeg;
   else return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize some of the fields of TGraphPolargram.

void TGraphPolargram::Init()
{
   fAxisAngle        = 0;
   fCutRadial        = 0;
   fDegree           = kFALSE;
   fGrad             = kFALSE;
   fLineStyle        = 3;
   fPolarLabelColor  = 1;
   fPolarLabelFont   = 62;
   fPolarOffset      = 0.04;
   fPolarTextSize    = 0.04;
   fRadialOffset     = 0.025;
   fRadian           = kTRUE;
   fRadialLabelColor = 1;
   fRadialLabelFont  = 62;
   fRadialTextSize   = 0.035;
   fTickpolarSize    = 0.02;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint TGraphPolargram.

void TGraphPolargram::Paint(Option_t * chopt)
{
   Int_t optionpoldiv, optionraddiv;
   Bool_t optionLabels = kTRUE;

   TString opt = chopt;
   opt.ToUpper();

   if(opt.Contains('P')) optionpoldiv=1; else optionpoldiv=0;
   if(opt.Contains('R')) optionraddiv=1; else optionraddiv=0;
   if(opt.Contains('O')) SetBit(TGraphPolargram::kLabelOrtho);
   else ResetBit(TGraphPolargram::kLabelOrtho);
   if(!opt.Contains('P') && !opt.Contains('R')) optionpoldiv=optionraddiv=1;
   if(opt.Contains('N')) optionLabels = kFALSE;

   if(optionraddiv) PaintRadialDivisions(kTRUE);
   else PaintRadialDivisions(kFALSE);
   if(optionpoldiv) PaintPolarDivisions(optionLabels);
}

////////////////////////////////////////////////////////////////////////////////
/// This is simplified from TEllipse::PaintEllipse.
/// Draw this ellipse with new coordinates.

void TGraphPolargram::PaintCircle(Double_t x1, Double_t y1, Double_t r,
                            Double_t phimin, Double_t phimax, Double_t theta)
{
   Int_t i;
   const Int_t np = 200; // Number of point to draw circle
   static Double_t x[np+3], y[np+3];

   // Set number of points approximatively proportional to the ellipse
   // circumference.

   Double_t circ = TMath::Pi()*2*r*(phimax-phimin)/36;
   Int_t n = (Int_t)(np*circ/((gPad->GetX2()-gPad->GetX1())+
                              (gPad->GetY2()-gPad->GetY1())));
   if (n < 8) n  = 8;
   if (n > np) n = np;
   Double_t angle,dx,dy;
   Double_t dphi = (phimax-phimin)*TMath::Pi()/(180*n);
   Double_t ct   = TMath::Cos(TMath::Pi()*theta/180);
   Double_t st   = TMath::Sin(TMath::Pi()*theta/180);
   for (i=0; i<=n; i++) {
      angle = phimin*TMath::Pi()/180 + Double_t(i)*dphi;
      dx    = r*TMath::Cos(angle);
      dy    = r*TMath::Sin(angle);
      x[i]  =  x1 + dx*ct - dy*st;
      y[i]  =  y1 + dx*st + dy*ct;
   }
   gPad->PaintPolyLine(n+1,x,y);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw Polar divisions.
/// Check for editable pad or create default.

void TGraphPolargram::PaintPolarDivisions(Bool_t optionLabels)
{
   Int_t i, j, rnum, rden, first, last;
   if (!gPad) return ;

   gPad->RangeAxis(-1,-1,1,1);
   gPad->Range(-1.25,-1.25,1.25,1.25);
   Int_t ndivMajor = fNdivPol%100;
   Int_t ndivMinor = fNdivPol/100;

   if (!gPad->GetLogy()) {
      for (i=0; i<ndivMajor; i++) {
         Double_t txtval    = fRwtmin + i*(fRwtmax-fRwtmin)/ndivMajor;
         Double_t theta     = i*2*TMath::Pi()/ndivMajor;
         Double_t costheta  = TMath::Cos(theta);
         Double_t sintheta  = TMath::Sin(theta);
         Double_t tantheta  = TMath::Tan(theta);
         Double_t costhetas = (1+fPolarOffset)*costheta;
         Double_t sinthetas = (1+fPolarOffset)*sintheta;
         Double_t corr = 0.01;

         TLatex *textangular = new TLatex();
         textangular->SetTextColor(GetPolarColorLabel());
         textangular->SetTextFont(GetPolarLabelFont());

         const char* form = (char *)" ";
         TGaxis axis;
         if (TestBit(TGraphPolargram::kLabelOrtho)) {
            // Polar numbers are aligned with their axis.
            if(fPolarLabels == NULL && optionLabels){;
               if (fRadian) {
                  // Radian case.
                  ReduceFraction(2*i, ndivMajor, rnum, rden); // Reduces the fraction.
                  if (rnum == 0)                       form = Form("%d",rnum);
                  if (rnum == 1 && rden == 1)          form = Form("#pi");
                  if (rnum == 1 && rden != 1)          form = Form("#frac{#pi}{%d}",rden);
                  if (rnum != 1 && rden == 1 && i !=0) form= Form("%d#pi",rnum);
                  if (rnum != 1 && rden != 1)          form = Form("#frac{%d#pi}{%d}",rnum,rden);
                  textangular->SetTextAlign(FindAlign(theta));
                  textangular->PaintLatex(costhetas,
                                          sinthetas, FindTextAngle(theta),
                                          GetPolarLabelSize(), form);
               } else {
                  // Any other cases: numbers are aligned with their axis.
                  form = Form("%5.3g",txtval);
                  axis.LabelsLimits(form,first,last);
                  TString s = Form("%s",form);
                  if (first != 0) s.Remove(0, first);
                  textangular->SetTextAlign(FindAlign(theta));
                  textangular->PaintLatex(costhetas,
                                          sinthetas, FindTextAngle(theta),
                                          GetPolarLabelSize(), s);
               }
            } else if (fPolarLabels){
               // print the specified polar labels
               textangular->SetTextAlign(FindAlign(theta));
               textangular->PaintLatex(costhetas,sinthetas,FindTextAngle(theta),
                                       GetPolarLabelSize(), fPolarLabels[i]);
            }
         } else {
            // Polar numbers are shown horizontally.
            if(fPolarLabels == NULL && optionLabels){
               if (fRadian) {
               // Radian case
                  ReduceFraction(2*i, ndivMajor, rnum, rden);
                  if (rnum == 0) form = Form("%d",rnum);
                  if (rnum == 1 && rden == 1)          form = Form("#pi");
                  if (rnum == 1 && rden != 1)          form = Form("#frac{#pi}{%d}",rden);
                  if (rnum != 1 && rden == 1 && i !=0) form = Form("%d#pi",rnum);
                  if (rnum != 1 && rden != 1)          form = Form("#frac{%d#pi}{%d}",rnum,rden);
                  if(theta >= 3*TMath::Pi()/12.0 && theta < 2*TMath::Pi()/3.0) corr=0.04;
                  textangular->SetTextAlign(FindAlign(theta));
                  textangular->PaintLatex(costhetas,corr+sinthetas,0,
                                          GetPolarLabelSize(),form);
               } else {
               // Any other cases where numbers are shown horizontally.
                  form = Form("%5.3g",txtval);
                  axis.LabelsLimits(form,first,last);
                  TString s = Form("%s",form);
                  if (first != 0) s.Remove(0, first);
                  if(theta >= 3*TMath::Pi()/12.0 && theta < 2*TMath::Pi()/3.0) corr=0.04;
                  textangular->SetTextAlign(FindAlign(theta));
                  textangular->PaintLatex(costhetas, //j'ai efface des offset la
                                          corr+sinthetas,0,GetPolarLabelSize(),s);
               }
            } else if (fPolarLabels) {
               // print the specified polar labels
               textangular->SetTextAlign(FindAlign(theta));
               textangular->PaintText(costhetas,sinthetas,fPolarLabels[i]);
            }
         }
         TAttLine::Modify();
       //Check if SetTickPolar is activated, and draw tick marks
         Bool_t issettickpolar = gPad->GetTicky();

         if (issettickpolar) {
            if (theta != 0 && theta !=TMath::Pi()) {
               gPad->PaintLine((sintheta-GetTickpolarSize())/tantheta,sintheta-GetTickpolarSize(),
               (sintheta+GetTickpolarSize())/tantheta,sintheta+GetTickpolarSize());
            }
            if (theta == 0 || theta ==TMath::Pi()) {
               gPad->PaintLine(1-GetTickpolarSize(),0,1+GetTickpolarSize(),0);
               gPad->PaintLine(-1+GetTickpolarSize(),0,-1-GetTickpolarSize(),0);
            }
         }
         TAttLine::SetLineStyle(1);
         TAttLine::Modify();
         gPad->PaintLine(0.,0.,costheta,sintheta);
         delete textangular;
       // Add minor lines w/o text.
         Int_t oldLineStyle = GetLineStyle();
         TAttLine::SetLineStyle(2);  //Minor lines always in this style.
         TAttLine::Modify();  //Changes line attributes apart from style.
         for (j=1; j<ndivMinor; j++) {
            Double_t thetamin = theta+j*2*TMath::Pi()/(ndivMajor*ndivMinor);
            gPad->PaintLine(0.,0.,TMath::Cos(thetamin),TMath::Sin(thetamin));
         }
         TAttLine::SetLineStyle(oldLineStyle);
         TAttLine::Modify();
      }
   } else {
         Int_t big = (Int_t)fRwtmax;
         Int_t test= 1;
         while (big >= 10) {
            big = big/10;
            test++;
         }
         for (i=1; i<=test; i++) {
         Double_t txtval    = pow((double)10,(double)(i-1));
         Double_t theta     = (i-1)*2*TMath::Pi()/(double)(test);
         Double_t costheta  = TMath::Cos(theta);
         Double_t sintheta  = TMath::Sin(theta);
         Double_t tantheta  = TMath::Tan(theta);
         Double_t costhetas = (1+fPolarOffset)*costheta;
         Double_t sinthetas = (1+fPolarOffset)*sintheta;
         Double_t corr      = 0.01;

         TLatex *textangular = new TLatex();
         textangular->SetTextColor(GetPolarColorLabel());
         textangular->SetTextFont(GetPolarLabelFont());

         const char* form = (char *)" ";
         TGaxis axis;

         if (TestBit(TGraphPolargram::kLabelOrtho)) {
            if(fPolarLabels==NULL && optionLabels){
            // Polar numbers are aligned with their axis.
               form = Form("%5.3g",txtval);
               axis.LabelsLimits(form,first,last);
               TString s = Form("%s",form);
               if (first != 0) s.Remove(0, first);
               textangular->SetTextAlign(FindAlign(theta));
               textangular->PaintLatex(costhetas,
                                       sinthetas, FindTextAngle(theta), GetPolarLabelSize(), s);
            }
            else if (fPolarLabels){
               // print the specified polar labels
               textangular->SetTextAlign(FindAlign(theta));
               textangular->PaintText(costhetas,sinthetas,fPolarLabels[i]);
            }

         } else {
            if(fPolarLabels==NULL && optionLabels){
            // Polar numbers are shown horizontally.
               form = Form("%5.3g",txtval);
               axis.LabelsLimits(form,first,last);
               TString s = Form("%s",form);
               if (first != 0) s.Remove(0, first);
               if(theta >= 3*TMath::Pi()/12.0 && theta < 2*TMath::Pi()/3.0) corr=0.04;
               textangular->SetTextAlign(FindAlign(theta));
               textangular->PaintLatex(costhetas,
                                       corr+sinthetas,0,GetPolarLabelSize(),s);
            } else if (fPolarLabels){
               // print the specified polar labels
               textangular->SetTextAlign(FindAlign(theta));
               textangular->PaintText(costhetas,sinthetas,fPolarLabels[i]);
            }
         }

         TAttLine::Modify();
         //Check if SetTickPolar is activated, and draw tick marks
         Bool_t issettickpolar = gPad->GetTicky();
         if (issettickpolar) {
            if (theta != 0 && theta !=TMath::Pi()) {
               gPad->PaintLine((sintheta-GetTickpolarSize())/tantheta,sintheta-GetTickpolarSize(),
               (sintheta+GetTickpolarSize())/tantheta,sintheta+GetTickpolarSize());
            }
            if (theta == 0 || theta ==TMath::Pi()) {
               gPad->PaintLine(1-GetTickpolarSize(),0,1+GetTickpolarSize(),0);
               gPad->PaintLine(-1+GetTickpolarSize(),0,-1-GetTickpolarSize(),0);
            }
         }
         TAttLine::SetLineStyle(1);
         TAttLine::Modify();
         gPad->PaintLine(0.,0.,costheta,sintheta);
         delete textangular;
         // Add minor lines w/o text.
         Int_t oldLineStyle = GetLineStyle();
         TAttLine::SetLineStyle(2);  //Minor lines always in this style.
         TAttLine::Modify();  //Changes line attributes apart from style.
         Double_t a=0;
         Double_t b,c,d;
         b = TMath::Log(10)*test;
         d= 2*TMath::Pi()/(double)test;
         for (j=1; j<9; j++) {
            a=TMath::Log(j+1)-TMath::Log(j)+a;
            c=a/b*6.28+d*(i-1);
            gPad->PaintLine(0.,0.,TMath::Cos(c),TMath::Sin(c));
         }
         TAttLine::SetLineStyle(oldLineStyle);
         TAttLine::Modify();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint radial divisions.
/// Check for editable pad or create default.

void TGraphPolargram::PaintRadialDivisions(Bool_t drawaxis)
{
   static char chopt[8] = "";
   Int_t i,j;
   Int_t ndiv      = TMath::Abs(fNdivRad);
   Int_t ndivMajor = ndiv%100;
   Int_t ndivMinor = ndiv/100;
   Int_t ndivmajor = 0;
   Double_t frwrmin = 0., frwrmax = 0., binWidth = 0;

   THLimitsFinder::Optimize(fRwrmin,fRwrmax,ndivMajor,frwrmin,
                               frwrmax, ndivmajor,binWidth,"");

   if (!gPad) return ;
   if (!gPad->GetLogx()) {
      gPad->RangeAxis(-1,-1,1,1);
      gPad->Range(-1.25,-1.25,1.25,1.25);
      Double_t umin  = fRwrmin;
      Double_t umax  = fRwrmax;
      Double_t rmajmin  = (frwrmin-fRwrmin)/(fRwrmax-fRwrmin);
      Double_t rmajmax  = (frwrmax-fRwrmin)/(fRwrmax-fRwrmin);
      Double_t dist  = (rmajmax-rmajmin)/ndivmajor;
      Int_t ndivminor = 0;

      chopt[0] = 0;
      strncat(chopt, "SDH", 4);
      if (fNdivRad < 0) strncat(chopt, "N",2);
      if(drawaxis){
      // Paint axis.
         TGaxis axis;
         axis.SetLabelSize(GetRadialLabelSize());
         axis.SetLabelColor(GetRadialColorLabel());
         axis.SetLabelFont(GetRadialLabelFont());
         axis.SetLabelOffset(GetRadialOffset());
         axis.PaintAxis(0, 0, TMath::Cos(GetAngle()), TMath::Sin(GetAngle()),
                                        umin, umax,  ndiv, chopt, 0., kFALSE);
      }

      // Paint Circles.
      // First paint main circle.
      PaintCircle(0.,0.,1,0.,360,0);
      // Optimised case.
      if (fNdivRad>0 ) {
         Double_t frwrmini = 0., frwrmaxi = 0., binWidth2 =0;
         THLimitsFinder::Optimize(frwrmin,frwrmin+binWidth,ndivMinor,frwrmini,
                                    frwrmaxi, ndivminor,binWidth2,"");
         Double_t dist2 = dist/(ndivminor);
         // Paint major circles.
         for (i=1; i<=ndivmajor+2; i++) {
            TAttLine::SetLineStyle(1);
            TAttLine::Modify();
            PaintCircle(0.,0.,rmajmin,0.,360,0);

            //Paint minor circles.
            TAttLine::SetLineStyle(2);
            TAttLine::Modify();
            for (j=1; j<ndivminor+1; j++) {
               if (rmajmin+j*dist2<=1) PaintCircle(0.,0.,rmajmin+j*dist2,0.,360,0);
            }
            rmajmin = (frwrmin-fRwrmin)/(fRwrmax-fRwrmin)+(i-1)*dist;
         }
      // Non-optimized case.
      } else {

         // Paint major circles.
         for (i=1; i<=ndivMajor; i++) {
            TAttLine::SetLineStyle(1);
            TAttLine::Modify();
            Double_t rmaj = i*1./ndivMajor;
            PaintCircle(0.,0.,rmaj,0.,360,0);

            // Paint minor circles.
            for (j=1; j<ndivMinor; j++) {
               TAttLine::SetLineStyle(2);
               TAttLine::Modify();
               PaintCircle(0.,0.,rmaj- j*1./(ndivMajor*ndivMinor),0.,360,0);
            }
         }
      }
   } else {
   // Draw Log scale on radial axis if option activated.
      Int_t big = (Int_t)fRwrmax;
      Int_t test= 0;
      while (big >= 10) {
         big = big/10;
         test++;
      }
      for (i=1; i<=test; i++) {
         TAttLine::SetLineStyle(1);
         TAttLine::Modify();
         Double_t ecart;
         ecart = ((double) i)/ ((double) test);
         PaintCircle(0.,0.,ecart,0,360,0);
         TAttLine::SetLineStyle(GetLineStyle());
         TAttLine::Modify();
         Double_t a=0;
         Double_t b,c,d;
         b = TMath::Log(10)*test;
         d = 1/(double)test;
         for (j=1; j<9; j++) {
            a = TMath::Log(j+1)-TMath::Log(j)+a;
            c = a/b+d*(i-1);
            PaintCircle(0,0.,c,0.,360,0);
         }
      }
   }
   TAttLine::SetLineStyle(1);
   TAttLine::Modify();
}

////////////////////////////////////////////////////////////////////////////////
/// Reduce fractions.

void TGraphPolargram::ReduceFraction(Int_t num, Int_t den, Int_t &rnum, Int_t &rden)
{
   Int_t a = 0;
   Int_t b = 0;
   Int_t i = 0;
   Int_t j = 0;
   a = den;
   b = num;
   if (b > a) {
      j = b;
   } else {
      j = a;
   }
   for (i=j; i > 1; i--) {
      if ((a % i == 0) && (b % i == 0)) {
         a = a/i;
         b = b/i;
      }
   }
   rden = a;
   rnum = b;
}

////////////////////////////////////////////////////////////////////////////////
/// Set axis angle.

void TGraphPolargram::SetAxisAngle(Double_t angle)
{
   fAxisAngle = angle/180*TMath::Pi();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the number of Polar divisions: enter a number ij with 0<i<99 and 0<j<99
/// - i sets the major Polar divisions.
/// - j sets the minor Polar divisions.

void TGraphPolargram::SetNdivPolar(Int_t ndiv)
{
   if (ndiv > 0)
      fNdivPol = ndiv;
   if (gPad) gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the number of radial divisions: enter a number ij with 0<i<99 and 0<j<99
/// - i sets the major radial divisions.
/// - j sets the minor radial divisions.

void TGraphPolargram::SetNdivRadial(Int_t ndiv)
{
   fNdivRad = ndiv;
   if (gPad) gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// Set some specified polar labels, used in the case of a spider plot.

void TGraphPolargram::SetPolarLabel(Int_t div, const TString & label)
{
   if(fPolarLabels == NULL)
      fPolarLabels = new TString[fNdivPol];
   fPolarLabels[div]=label;
   if (gPad) gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// Set Polar labels color.

void TGraphPolargram::SetPolarLabelColor(Color_t tcolorangular )
{
   fPolarLabelColor = tcolorangular;
}

////////////////////////////////////////////////////////////////////////////////

void TGraphPolargram::SetPolarLabelFont(Font_t tfontangular)
{
   // Set Polar label font.

   fPolarLabelFont = tfontangular;
}

////////////////////////////////////////////////////////////////////////////////
/// Set angular labels size.

void TGraphPolargram::SetPolarLabelSize(Double_t angularsize )
{
   fPolarTextSize = angularsize;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the labels offset.

void TGraphPolargram::SetPolarOffset(Double_t angularOffset)
{
   fPolarOffset = angularOffset;
   if (gPad) gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// Set radial labels color.

void TGraphPolargram::SetRadialLabelColor(Color_t tcolorradial )
{
   fRadialLabelColor = tcolorradial;
}

////////////////////////////////////////////////////////////////////////////////
/// Set radial label font.

void TGraphPolargram::SetRadialLabelFont(Font_t tfontradial)
{
   fRadialLabelFont = tfontradial;
}

////////////////////////////////////////////////////////////////////////////////
/// Set radial labels size.

void TGraphPolargram::SetRadialLabelSize(Double_t radialsize )
{
   fRadialTextSize = radialsize;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the labels offset.

void TGraphPolargram::SetRadialOffset(Double_t radialOffset)
{
   fRadialOffset = radialOffset;
   if (gPad) gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// Allows to change range Polar.
/// \param[in] tmin   the start number.
/// \param[in] tmax   the end number.

void TGraphPolargram::SetRangePolar(Double_t tmin, Double_t tmax)
{
   fDegree = kFALSE;
   fGrad   = kFALSE;
   fRadian = kFALSE;

   if (tmin < tmax) {
      fRwtmin = tmin;
      fRwtmax = tmax;
   }
   if (gPad) gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the radial range.
/// \param[in] rmin   radius at center of the circle.
/// \param[in] rmax   radius at the intersection of the right X axis part and the circle.

void TGraphPolargram::SetRangeRadial(Double_t rmin, Double_t rmax)
{
   if (rmin < rmax) {
      fRwrmin = rmin;
      fRwrmax = rmax;
   }
   if (gPad) gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// Set polar ticks size.

void TGraphPolargram::SetTickpolarSize(Double_t tickpolarsize)
{
   fTickpolarSize = tickpolarsize;
}

////////////////////////////////////////////////////////////////////////////////
/// The Polar circle is labelled using degree.

void TGraphPolargram::SetToDegree()
{
   fDegree = kTRUE;
   fGrad   = kFALSE;
   fRadian = kFALSE;
   ChangeRangePolar(0,360);
}

////////////////////////////////////////////////////////////////////////////////
/// The Polar circle is labelled using gradian.

void TGraphPolargram::SetToGrad()
{
   fGrad   = kTRUE;
   fRadian = kFALSE;
   fDegree = kFALSE;
   ChangeRangePolar(0,200);
}

////////////////////////////////////////////////////////////////////////////////
/// The Polar circle is labelled using radian.

void TGraphPolargram::SetToRadian()
{
   fRadian = kTRUE;
   fGrad   = kFALSE;
   fDegree = kFALSE;
   ChangeRangePolar(0,2*TMath::Pi());
}

////////////////////////////////////////////////////////////////////////////////
/// Set range from 0 to 2*pi

void TGraphPolargram::SetTwoPi()
{
   SetRangePolar(0,2*TMath::Pi());
}
