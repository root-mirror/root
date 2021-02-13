// @(#)root/spectrumpainter:$Id: TSpectrum2Painter.cxx,v 1.00
// Author: Miroslav Morhac 29/09/06

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


/** \class TSpectrum2Painter
    \ingroup Spectrumpainter

Two-dimensional graphics function

TSpectrum2Painter is a set of graphical functions developed by Miroslav
Morhac to paint 2D-histograms in three dimensions. This package is accessed
via THistPainter in a transparent way. For the ROOT user it is enough to use
the "SPEC" option to draw a 2D-Histogram. This option offers many
functionalities detailed in the header of the PaintSpectrum function.

Reference:

Morhac M., Kliman J., Matousek V., Turzo I.: Sophisticated visualization
algorithms for analysis of multidimensional experimental nuclear data. Acta
Pysica Slovaca Vol. 54/ 4 (2004), pp. 385-400.
*/

#include <climits>

#include "TROOT.h"
#include "TColor.h"
#include "TMath.h"
#include "TLine.h"
#include "TEllipse.h"
#include "TVirtualPad.h"
#include "TBox.h"
#include "TF1.h"
#include "TH2.h"
#include "TGaxis.h"
#include "THLimitsFinder.h"
#include "TSpectrum2Painter.h"
#include "strlcpy.h"

ClassImp (TSpectrum2Painter)


////////////////////////////////////////////////////////////////////////////////
/// TSpectrum2Painter normal constructor

TSpectrum2Painter::TSpectrum2Painter(TH2* h2, Int_t bs)
   : TNamed ("Spectrum Painter2","Miroslav Morhac Painter")
{
   int i, j;
   double val;
   gPad->Range(0, 0, 1 ,1);
   fXmin = 0;
   fXmax = h2->GetNbinsX() - 1;
   fYmin = 0;
   fYmax = h2->GetNbinsY() - 1;
   fZmin = 0, fZmax = 0;
   fMaximumXScreenResolution = bs;

   for (i = 0;i <= fXmax; i++) {
      for (j = 0;j <= fYmax; j++) {
         val = h2->GetBinContent(i + 1,j + 1);
         if (val > fZmax) fZmax = val;
      }
   }

   fBx1 = gPad->XtoPixel(0.1); //axis positions
   fBx2 = gPad->XtoPixel(0.99);
   fBy1 = gPad->YtoPixel(0.99);
   fBy2 = gPad->YtoPixel(0.05);

   fModeGroup = kModeGroupLightHeight;

   fDisplayMode = kDisplayModeSurface;

   fZscale = kZScaleLinear; // Scale linear, log.

   fNodesx = fXmax-fXmin+1; // Number of nodes in x dimension of grid
   fNodesy = fYmax-fYmin+1; // Number of nodes in y dimension of grid

   fContWidth = 50; // Width between contours,
                    // applies only for contours display mode.
   fAlpha     = 20; // Angles of display,alfa+beta must be less or equal to 90,
                    // alpha- angle between base line of Canvas and left lower
                    // edge of picture picture base plane
   fBeta      = 60; // Angle between base line of Canvas and right lower edge
                    // of picture base plane
   fViewAngle = 0;  // Rotation angle of the view,
                    // it can be 0, 90, 180, 270 degrees.

   fLevels       = 256; // Number of color levels for rainbowed display modes,
                        // It does not apply for simple display modes
                        // algorithm group
   fRainbow1Step = 1;   // Determines the first component step for
                        // neighbouring color levels, applies only for
                        // rainbowed display modes, it does not apply for
                        // simple display modes algorithm group.
   fRainbow2Step = 1;   // Determines the second component step for
                        // neighbouring color levels, applies only for
                        // rainbowed display modes, it does not apply for
                        // simple display modes algorithm group.
   fRainbow3Step = 1;   // Determines the third component step for
                        // neighbouring color levels, applies only for
                        // rainbowed display modes, it does not apply for
                        // simple display modes algorithm group.

   fColorAlg = kColorAlgRgbSmooth; // Applies only for rainbowed display modes
                                   // (rgb smooth algorithm, rgb modulo color
                                   // component, cmy smooth algorithm, cmy
                                   // modulo color component, cie smooth
                                   // algorithm, cie modulo color component,
                                   // yiq smooth algorithm, yiq modulo color
                                   // component, hsv smooth algorithm, hsv
                                   // modulo color component, it does not
                                   // apply for simple display modes
                                   // algorithm group.

   fLHweight = 0.5; // Weight between shading according to fictive light
                    // source and according to channels counts, applies only
                    // for kModeGroupLightHeight modes group.

   fXlight   = 1000; // X position of fictive light source, applies only for
                     // rainbowed display modes with shading according to light.
   fYlight   = 1000; // Y position of fictive light source, applies only for
                     // rainbowed display modes with shading according to light.
   fZlight   = 1000; // Z position of fictive light source, applies only for
                     // rainbowed display modes with shading according to light.

   fShadow   = kShadowsNotPainted; // Determines whether shadow will be drawn
                                   // (no shadow, shadow), for rainbowed
                                   // display modes with shading according to
                                   // light.

   fShading  = kShaded; // Determines whether the picture will shaded,
                        // smoothed (no shading, shading), for rainbowed
                        // display modes only.

   fBezier   = kNoBezierInterpol; // Determines Bezier interpolation (applies
                                  // only for simple display modes group for
                                  // grid, x_lines, y_lines display modes).

   fPenColor = kBlack;         // Color of spectrum.
   fPenWidth = 1;              // Width of line.
   fPenDash  = kPenStyleSolid; // Style of pen.

   fChanmarkEnDis  = kChannelMarksNotDrawn; // Decides whether the channel
                                            // marks are shown.
   fChanmarkColor  = kBlue;                 // Color of channel marks.
   fChanmarkWidth  = 8;                     // Width of channel marks.
   fChanmarkHeight = 8;                     // Height of channel marks.
   fChanmarkStyle  = kChannelMarksStyleDot; // Style of channel marks.

   fChanlineEnDis   = kChannelGridNotDrawn; // Decides whether the channel lines
                                            // (grid) are shown.
   fChanlineColor   = kRed;                 // Color of channel marks.
   fNewColor        = 0;
   fEnvelope        = new Short_t [fMaximumXScreenResolution];
   fEnvelopeContour = new Short_t [fMaximumXScreenResolution];
   for (i=0;i<fMaximumXScreenResolution;i++) {
      fEnvelope[i]        = fBy2;
      fEnvelopeContour[i] = fBy2;
   }
   fH2 = h2;
}


////////////////////////////////////////////////////////////////////////////////
/// TSpectrum2Painter destructor

TSpectrum2Painter::~TSpectrum2Painter()
{
   TColor* col;
   for (int i=0; i<256; i++) {
      col = gROOT->GetColor(250+i);
      if (col) delete col;
   }
   if (fEnvelope) delete [] fEnvelope;
   if (fEnvelopeContour) delete [] fEnvelopeContour;
}


////////////////////////////////////////////////////////////////////////////////
/// Reads out the value from histogram and calculates screen coordinates
///
/// Parameters:
///    - it - node in x- direction
///    - jt - node in y- direction
///    - zmt - control variable

void TSpectrum2Painter::Transform(Int_t it,Int_t jt,Int_t zmt)
{
   Int_t lxt,lyt,ix,iy;
   Double_t zf = 0;
   Double_t p1,p2;
   p1        = fXmin+fKx*(Double_t)it;
   p2        = fYmin+fKy*(Double_t)jt;
   ix        = (Int_t)p1;
   iy        = (Int_t)p2;
   fDxspline = p1;
   fDyspline = p2;
   if ((zmt==0)||(zmt==-3)||(zmt==-4)) {
      zf = fH2->GetBinContent(ix+1,iy+1);
   } else if (zmt==-2) zf = fZPresetValue;
   if (zf<fZmin) zf = fZmin;
   fZeq = zf;
   switch (fZscale) {
      case kZScaleLog:
         if (zf>=1.0) zf = log(zf);
         else         zf = 0;
         break;
      case kZScaleSqrt:
         if (zf>0) zf = sqrt(zf);
         else      zf = 0;
         break;
   }
   lxt = (Int_t)(fTxx*(Double_t)it+fTxy*(Double_t)jt+fVx);
   lyt = (Int_t)(fTyx*(Double_t)it+fTyy*(Double_t)jt+fTyz*zf+fVy);
   if (lxt<fBx1) lxt = fBx1;
   if (lxt>fBx2) lxt = fBx2;
   if (lyt<fBy1) lyt = fBy1;
   if (lyt>fBy2) lyt = fBy2;
   fXt = lxt;
   fYt = lyt;
   fZ  = zf;
   return;
}


////////////////////////////////////////////////////////////////////////////////
/// Calculates and returns color value for the surface triangle
/// given by function parameters:
///    -dx1,dy1,z1 coordinates of the first point in 3d space
///    -dx2,dy2,z2 coordinates of the second point in 3d space
///    -dx3,dy3,z3 coordinates of the third point in 3d space

Double_t TSpectrum2Painter::ColorCalculation(
         Double_t dx1, Double_t dy1, Double_t z1,
         Double_t dx2, Double_t dy2, Double_t z2,
         Double_t dx3, Double_t dy3, Double_t z3)
{
   Double_t da,db,dc=0,dd,dl,dm,dn,xtaz,ytaz,ztaz,v=0,v1;
   Double_t pi=3.1415927;
   Int_t i;
   switch (fZscale) {
      case kZScaleLog:
         if (z1>900) z1 = 900;
         z1 = exp(z1);
         if (z2>900) z2 = 900;
         z2 = exp(z2);
         if (z3>900) z3 = 900;
         z3 = exp(z3);
         break;
      case kZScaleSqrt:
         z1 = z1*z1;
         z2 = z2*z2;
         z3 = z3*z3;
         break;
   }
   i = fViewAngle;
   i = i/90;
   if ((i==1)||(i==3)) {
      da  = dx1;
      dx1 = dx2;
      dx2 = da;
      da  = dy1;
      dy1 = dy2;
      dy2 = da;
      da  = z1;
      z1  = z2;
      z2  = da;
   }
   xtaz = (dx1+dx2+dx3)/3;
   ytaz = (dy1+dy2+dy3)/3;
   ztaz = (z1+z2+z3)/3;
   if (fModeGroup==kModeGroupLight) {
      dn = (Double_t)fZlight-ztaz;
      dm = (Double_t)fYlight-ytaz;
      dl = (Double_t)fXlight-xtaz;
      da = (dy2-dy1)*(z3-z1)-(dy3-dy1)*(z2-z1);
      db = (z2-z1)*(dx3-dx1)-(z3-z1)*(dx2-dx1);
      dc = (dx2-dx1)*(dy3-dy1)-(dx3-dx1)*(dy2-dy1);
      dd = (da*da+db*db+dc*dc)*(dl*dl+dm*dm+dn*dn);
      dd = sqrt(dd);
      if (dd!=0) v = (da*dl+db*dm+dc*dn)/dd;
      else       v = 0;
      if (v<-1) v=-1;
      if (v>1) v=1;
      v = asin(v);
      v = v+pi/2;
      v = v/pi;
   } else if (fModeGroup==kModeGroupHeight) {
      da = fZmax-fZmin;
      if (ztaz<fZmin) ztaz=fZmin;
      if (ztaz>=fZmax) ztaz=fZmax-1;
      db = ztaz-fZmin;
      if (da!=0) {
         switch (fZscale) {
            case kZScaleLinear:
               dc = db/da;
               break;
            case kZScaleLog:
               if (da>=1) da=log(da);
               if (db>=1) db=log(db);
               if (da!=0) dc=db/da;
               else dc=0;
               break;
            case kZScaleSqrt:
               da = sqrt(da);
               db = sqrt(db);
               dc = db/da;
               break;
         }
      } else {
         dc=0;
      }
      i = (Int_t)dc;
      v = dc-i;
   } else if (fModeGroup==kModeGroupLightHeight) {
      dn = (Double_t)fZlight-ztaz;
      dm = (Double_t)fYlight-ytaz;
      dl = (Double_t)fXlight-xtaz;
      da = (dy2-dy1)*(z3-z1)-(dy3-dy1)*(z2-z1);
      db = (z2-z1)*(dx3-dx1)-(z3-z1)*(dx2-dx1);
      dc = (dx2-dx1)*(dy3-dy1)-(dx3-dx1)*(dy2-dy1);
      dd = (da*da+db*db+dc*dc)*(dl*dl+dm*dm+dn*dn);
      dd = sqrt(dd);
      if (dd!=0) v = (da*dl+db*dm+dc*dn)/dd;
      else v = 0;
      if (v<-1) v=-1;
      if (v>1) v=1;
      v  = asin(v);
      v  = v+pi/2;
      v  = v/pi;
      da = fZmax-fZmin;
      if (ztaz<fZmin)  ztaz = fZmin;
      if (ztaz>=fZmax) ztaz = fZmax-1;
      db = ztaz-fZmin;
      if (da!=0) {
         switch (fZscale) {
            case kZScaleLinear:
               dc = db/da;
               break;
            case kZScaleLog:
               if (da>=1) da = log(da);
               if (db>=1) db = log(db);
               if (da!=0) dc = db/da;
               else       dc = 0;
               break;
            case kZScaleSqrt:
               da = sqrt(da);
               db = sqrt(db);
               dc = db/da;
               break;
         }
      } else {
         dc = 0;
      }
      i  = (Int_t)dc;
      v1 = dc-i;
      v  = fLHweight*v+(1-fLHweight)*v1;
   }
   if (fShadow==kShadowsNotPainted) {
      da = 1.0/(Double_t)fLevels;
      if (v<da) v = da;
   } else {
      da = 2.0/(Double_t)fLevels;
      if (v<da) v = da;
   }
   return(v);
}


////////////////////////////////////////////////////////////////////////////////
/// Determines whether the center of the triangle in 3-d space
/// given by function parameters:
///   - xtaz,ytaz,ztaz
///     is in shadow or not. If yes it return 1 otherwise it returns 0.

Double_t TSpectrum2Painter::ShadowColorCalculation(Double_t xtaz, Double_t ytaz,
                                                   Double_t ztaz,
                                                   Double_t shad_noise)
{
   Int_t sx2,sy2,sz1,sz2,skrokx,skroky,xmax,ymax;
   Double_t sx1,sy1;
   Double_t pom1,pom2,sdx1=0,sdx2=0,sdy1,sdy2,spriz;
   switch (fZscale) {
      case kZScaleLog:
         if (ztaz>900) ztaz = 900;
         ztaz = exp(ztaz);
         if (ztaz>32767) ztaz = 32767;
         break;
      case kZScaleSqrt:
         ztaz = ztaz*ztaz;
         break;
   }
   spriz = 0;
   sx1   = xtaz;
   sy1   = ytaz;
   sz1   = (Int_t)ztaz;
   sx2   = fXlight;
   sy2   = fYlight;
   sz2   = fZlight;
   xmax  = fXmax;
   ymax  = fYmax;
   if (sx1!=sx2) {
      if (sx1<sx2) skrokx =  1;
      else         skrokx = -1;
      if (sy1<sy2) skroky =  1;
      else         skroky = -1;
      pom1 = sx2-sx1;
      pom2 = sy2-sy1;
      if (TMath::Abs(pom1)>0.0000001) sdx1 = pom2/pom1;
      pom1 = sx1;
      pom2 = sy1;
      sdy1 = pom2-sdx1*pom1;
      pom1 = sx2-sx1;
      pom2 = sz2-sz1;
      if (TMath::Abs(pom1)>0.0000001) sdx2 = pom2/pom1;
      pom1  = sx1;
      pom2  = sz1;
      sdy2  = pom2-sdx2*pom1;
      spriz = 0;
      pom1  = sx1;
      pom2  = pom1*sdx1+sdy1;
      sy1   = pom2;
      for (;(sx1>(fXmin-skrokx)) && (sx1<(xmax-skrokx)) &&
            (sy1>(fYmin-skroky)) && (sy1<(ymax-skroky)) &&
            (spriz==0);sx1+=skrokx) {
         pom1 = sx1;
         pom2 = pom1*sdx1+sdy1;
         sy1  = pom2+skroky;
         if ((sy1>=fYmin)&&(sy1<=fYmax)) {
            sz1  = (Int_t)(fH2->GetBinContent((Int_t)sx1+1,(Int_t)sy1+1));
            pom2 = pom1*sdx2+sdy2;
            sz2  = (Int_t)(pom2+shad_noise);
            if (sz1>sz2) spriz = 1;
         }
      }
   } else if (sy1!=sy2) {
      if (sy1<sy2) skroky =  1;
      else         skroky = -1;
      pom1 = sy2-sy1;
      pom2 = sz2-sz1;
      if (TMath::Abs(pom1)>0.0000001) sdx2 = pom2/pom1;
      pom1  = sy1;
      pom2  = sz1;
      sdy2  = pom2-sdx2*pom1;
      spriz = 0;
      for (;(sy1>(fYmin-skroky)) && (sy1<(ymax-skroky)) &&
            (spriz==0);sy1+=skroky) {
         sz1  = (Int_t)(fH2->GetBinContent((Int_t)sx1+1,(Int_t)sy1+1));
         pom1 = sy1;
         pom2 = pom1*sdx2+sdy2;
         sz2  = (Int_t)(pom2+shad_noise);
         if (sz1>sz2) spriz=1;
      }
   }
   return(spriz);
}


////////////////////////////////////////////////////////////////////////////////
/// This function calculates color for one palette entry
/// given by function parameter ui. Other parameters
///    -ui1,ui2,ui3
///    represent r, g, b color components of the basic pen color.

void TSpectrum2Painter::ColorModel(unsigned ui, unsigned ui1, unsigned ui2,
                                   unsigned ui3)
{
   unsigned long uinc1=0,uinc2=0,uinc3=0,upom,i;
   Double_t a,b,c,d,h,v,s,f;
   Int_t j,iv=ui;
   Double_t red=0,green=0,blue=0;
   if (iv<0)        iv = 0;
   else if (iv>255) iv = 255;
   if (gROOT->GetColor(250+iv)) {
      fNewColorIndex = 250+iv;
      return;
   }
   if (fColorAlg%2==0) {
      a     = fRainbow1Step;
      a     = ui*a;
      a     = ui1+a;
      if (a >= UINT_MAX) uinc1 = UINT_MAX;
      else               uinc1 = (unsigned)a;
      upom  = uinc1%256;
      i     = (uinc1-upom)/256;
      if ((i%2)==0) uinc1 = upom;
      else          uinc1 = 255-upom;
      b     = fRainbow2Step;
      b     = ui*b;
      b     = ui2+b;
      uinc2 = (Int_t)b;
      upom  = uinc2%256;
      i     = (uinc2-upom)/256;
      if ((i%2)==0) uinc2 = upom;
      else          uinc2 = 255-upom;
      c     = fRainbow3Step;
      c     = ui*c;
      c     = ui3+c;
      uinc3 = (Int_t)c;
      upom  = uinc3%256;
      i     = (uinc3-upom)/256;
      if ((i%2)==0) uinc3 = upom;
      else          uinc3 = 255-upom;
      if (fColorAlg==kColorAlgCmySmooth) {
         uinc1 = 255-uinc1;
         uinc2 = 255-uinc2;
         uinc3 = 255-uinc3;
      } else if (fColorAlg==kColorAlgCieSmooth) {
         a = uinc1;
         b = uinc2;
         c = uinc3;
         d = a+b+c;
         if (d!=0) {
            a = a/d;
            b = b/d;
            c = c/d;
         }
         red   = a*255;
         green = b*255;
         blue  = c*255;
         uinc1 = (Int_t)red;
         uinc2 = (Int_t)green;
         uinc3 = (Int_t)blue;
      } else if (fColorAlg==kColorAlgYiqSmooth) {
         a     = uinc1;
         b     = uinc2;
         c     = uinc3;
         a     = a/256;
         b     = b/256;
         c     = c/256;
         red   = a+0.956*b+0.62*c;
         green = a-0.272*b-0.647*c;
         blue  = a-1.108*b+1.705*c;
         if (red>=2)       red   = red-2;
         else if (red>=1)  red   = 2-red;
         if (green<0)      green = -green;
         if (blue>=2)      blue  = blue-2;
         else if (blue>=1) blue  = 2-blue;
         else if (blue<-1) blue  = 2+blue;
         else if (blue<0)  blue  = -blue;
         red   = red*255;
         green = green*255;
         blue  = blue*255;
         uinc1 = (Int_t)red;
         uinc2 = (Int_t)green;
         uinc3 = (Int_t)blue;
      } else if (fColorAlg==kColorAlgHvsSmooth) {
         h = uinc1;
         v = uinc2;
         s = uinc3;
         h = h/256;
         v = v/256;
         s = s/256;
         if (s==0) {
            red   = v;
            green = v;
            blue  = v;
         } else {
            if (h==1.0) h=0;
            h = 6.0*h;
            j = (Int_t)h;
            f = h-j;
            a = v*(1-s);
            b = v*(1-s*f);
            c = v*(1-s*(1-f));
            switch (j) {
               case 0:
                  red   = v;
                  green = c;
                  blue  = a;
                  break;
               case 1:
                  red   = b;
                  green = v;
                  blue  = a;
                  break;
               case 2:
                  red   = a;
                  green = v;
                  blue  = c;
                  break;
               case 3:
                  red   = a;
                  green = b;
                  blue  = v;
                  break;
               case 4:
                  red   = c;
                  green = a;
                  blue  = v;
                  break;
               case 5:
                  red   = v;
                  green = a;
                  blue  = b;
                  break;
            }
         }
         red   = red*255;
         green = green*255;
         blue  = blue*255;
         uinc1 = (Int_t)red;
         uinc2 = (Int_t)green;
         uinc3 = (Int_t)blue;
      }
      ui = uinc1+uinc2*256+uinc3*65536;
   } else if (fColorAlg%2==1) {
      a     = fRainbow1Step;
      a     = ui*a;
      a     = ui1/2+a;
      uinc1 = (Int_t)a;
      uinc1 = uinc1%256;
      b     = fRainbow2Step;
      b     = ui*b;
      b     = ui2/2+b;
      uinc2 = (Int_t)b;
      uinc2 = uinc2%256;
      c     = fRainbow3Step;
      c     = ui*c;
      c     = ui3/2+c;
      uinc3 = (Int_t)c;
      uinc3 = uinc3%256;
      if (fColorAlg==kColorAlgCmyModulo) {
         uinc1 = 255-uinc1;
         uinc2 = 255-uinc2;
         uinc3 = 255-uinc3;
      } else if (fColorAlg==kColorAlgCieModulo) {
         a = uinc1;
         b = uinc2;
         c = uinc3;
         d = a+b+c;
         if (d!=0) {
            a = a/d;
            b = b/d;
            c = c/d;
         }
         red   = a*255;
         green = b*255;
         blue  = c*255;
         uinc1 = (Int_t)red;
         uinc2 = (Int_t)green;
         uinc3 = (Int_t)blue;
      } else if (fColorAlg==kColorAlgYiqModulo) {
         a     = uinc1;
         b     = uinc2;
         c     = uinc3;
         a     = a/256;
         b     = b/256;
         c     = c/256;
         red   = a+0.956*b+0.62*c;
         green = a-0.272*b-0.647*c;
         blue  = a-1.108*b+1.705*c;
         if (red>=2)       red   = red-2;
         else if (red>=1)  red   = red-1;
         if (green<0)      green = 1+green;
         if (blue>=2)      blue  = blue-2;
         else if (blue>=1) blue  = blue-1;
         else if (blue<-1) blue  = 2+blue;
         else if (blue<0)  blue  = 1+blue;
         red   = red*255;
         green = green*255;
         blue  = blue*255;
         uinc1 = (Int_t)red;
         uinc2 = (Int_t)green;
         uinc3 = (Int_t)blue;
      } else if (fColorAlg==kColorAlgHvsModulo) {
         h = uinc1;
         v = uinc2;
         s = uinc3;
         h = h/256;
         v = v/256;
         s = s/256;
         if (s==0) {
            red   = v;
            green = v;
            blue  = v;
         } else {
            if (h==1.0) h = 0;
            h = 6.0*h;
            j = (Int_t)h;
            f = h-j;
            a = v*(1-s);
            b = v*(1-s*f);
            c = v*(1-s*(1-f));
            switch (j) {
               case 0:
                  red   = v;
                  green = c;
                  blue  = a;
                  break;
               case 1:
                  red   = b;
                  green = v;
                  blue  = a;
                  break;
               case 2:
                  red   = a;
                  green = v;
                  blue  = c;
                  break;
               case 3:
                  red   = a;
                  green = b;
                  blue  = v;
                  break;
               case 4:
                  red   = c;
                  green = a;
                  blue  = v;
                  break;
               case 5:
                  red   = v;
                  green = a;
                  blue  = b;
                  break;
            }
         }
         red   = red*255;
         green = green*255;
         blue  = blue*255;
         uinc1 = (Int_t)red;
         uinc2 = (Int_t)green;
         uinc3 = (Int_t)blue;
      }
      ui = uinc1+uinc2*256+uinc3*65536;
   }
   red            = uinc1;
   green          = uinc2;
   blue           = uinc3;
   red            = red/255.0;
   green          = green/255.0;
   blue           = blue/255.0;
   fNewColor      = new TColor(250+iv,red,green,blue);
   fNewColorIndex = 250+iv;
   return;
}


////////////////////////////////////////////////////////////////////////////////
/// This function is called from BezierBlend function.

Int_t TSpectrum2Painter::BezC(Int_t i)
{
   Int_t j,a;
   a = 1;
   for (j=i+1;j<=3;j++) a = a*j;
   for (j=1;j<=3-i;j++) a = a/j;
   return a;
}


////////////////////////////////////////////////////////////////////////////////
/// This function calculates Bezier approximation.

Double_t TSpectrum2Painter::BezierBlend(Int_t i,Double_t bezf)
{
   Int_t j;
   Double_t v;
   v = BezC(i);
   for (j=1;j<=i;j++)   v = v*bezf;
   for (j=1;j<=3-i;j++) v = v*(1-bezf);
   return v;
}


////////////////////////////////////////////////////////////////////////////////
/// Calculates screen coordinates of the smoothed point.
/// Parameter bezf changes within the interval 0 to 1 in 0.1 steps.

void TSpectrum2Painter::BezierSmoothing(Double_t bezf)
{
   Int_t i;
   Double_t b;
   fGbezx = 0;
   fGbezy = 0;
   for (i=0;i<4;i++) {
      b       = BezierBlend(i,bezf);
      fGbezx += fBzX[i]*b;
      fGbezy += fBzY[i]*b;
   }
   return;
}


////////////////////////////////////////////////////////////////////////////////
/// Ensures hidden surface removal.

void TSpectrum2Painter::Envelope(Int_t x1,Int_t y1,Int_t x2,Int_t y2)
{
   Int_t x,y,krok,xold=0,yold=0,prvy,yprv=0;
   Double_t fx,fy,fx1,fy1;
   if (y1<fBy1) y1 = fBy1;
   if (y2<fBy1) y2 = fBy1;
   if (x1==x2) {
      if ((y1>=fEnvelope[x1]) && (y2>=fEnvelope[x1])) {
         if (x1>0) {
            if (y1<=fEnvelope[x1-1]||y2<=fEnvelope[x1-1]) {
               if (y1>fEnvelope[x1-1]) y1 = fEnvelope[x1-1];
               if (y2>fEnvelope[x1-1]) y2 = fEnvelope[x1-1];
               fLine = 2;
               fXs   = x1;
               fYs   = y1;
               fXe   = x2;
               fYe   = y2;
               return;
            }
         }
         if (x1<fBx2) {
            if (y1<=fEnvelope[x1+1]||y2<=fEnvelope[x1+1]) {
               if (y1>fEnvelope[x1+1]) y1 = fEnvelope[x1+1];
               if (y2>fEnvelope[x1+1]) y2 = fEnvelope[x1+1];
               fLine = 2;
               fXs   = x1;
               fYs   = y1;
               fXe   = x2;
               fYe   = y2;
               return;
            }
         }
         fLine=0;
         return;
      }
      if ((y1<fEnvelope[x1]) && (y2<fEnvelope[x1])) {
         fLine = 2;
         fXs   = x1;
         fYs   = y1;
         fXe   = x2;
         fYe   = y2;
         if (y1<y2) fEnvelope[x1] = y1;
         else       fEnvelope[x1] = y2;
         return;
      }
      if (y1<y2) {
         fLine = 2;
         fXs   = x1;
         fYs   = y1;
         fXe   = x1;
         fYe   = fEnvelope[x1];
         fEnvelope[x1] = y1;
         return;
      } else {
         fLine = 2;
         fXs   = x1;
         fYs   = y2;
         fXe   = x1;
         fYe   = fEnvelope[x1];
         fEnvelope[x1] = y2;
         return;
      }
   }
   krok  = (x1<x2)? 1:-1;
   fLine = 0;
   prvy  = 0;
   x     = x1;
   y     = y1;
l1:
   if (y<=fEnvelope[x]) {
      xold = x;
      yold = y;
      if (fLine==0) {
         fLine = 1;
         if (prvy==1) {
            if (yprv<=fEnvelope[x]) fYs = yprv;
            else                    fYs = fEnvelope[x];
            fXs = x;
         } else {
            fXs = x;
            fYs = y;
         }
      }
      if (x!=x2) fEnvelope[x] = y;
   } else {
      prvy = 1;
      yprv = y;
      if (fLine==1) {
         fLine = 2;
         fXe   = xold;
         fYe   = yold;
      }
   }
   if (x1==x2) {
      if (y1!=y2) y += (y1<y2)? +1:-1;
      if (y!=y2) goto l1;
   } else {
      x  += krok;
      fy1 = y2-y1;
      fx1 = x2-x1;
      fx  = x-x1;
      fy  = fy1*fx/fx1;
      y   = (Int_t)(y1+fy);
      if (((x<=x2)&&(x1<x2)) || ((x>=x2)&&(x1>x2))) goto l1;
   }
   return;
}


////////////////////////////////////////////////////////////////////////////////
/// Ensures hidden surface removal for Bars, BarsX and BarsY
/// display modes.

void TSpectrum2Painter::EnvelopeBars(Int_t x1,Int_t y1,Int_t x2,Int_t y2)
{
   Int_t x,y,krok,xold=0,yold=0,prvy,xprv,yprv=0;
   Double_t fx,fy,fx1,fy1;
   if (x1==x2) {
      if ((y1>=fEnvelope[x1]) && (y2>=fEnvelope[x1])) {
         fLine = 0;
         return;
      }
      if ((y1<fEnvelope[x1]) && (y2<fEnvelope[x1])) {
         fLine = 2;
         fXs   = x1;
         fYs   = y1;
         fXe   = x2;
         fYe   = y2;
         if (y1<y2) fEnvelope[x1] = y1;
         else       fEnvelope[x1] = y2;
         return;
      }
      if (y1<y2) {
         fLine = 2;
         fXs   = x1;
         fYs   = y1;
         fXe   = x1;
         fYe   = fEnvelope[x1];
         fEnvelope[x1] = y1;
         return;
      } else {
         fLine = 2;
         fXs   = x1;
         fYs   = y2;
         fXe   = x1;
         fYe   = fEnvelope[x1];
         fEnvelope[x1] = y2;
         return;
      }
   }
   krok  = (x1<x2)? 1:-1;
   fLine = 0;
   prvy  = 0;
   x     = x1;
   y     = y1;
l1:
   if (y<=fEnvelope[x]) {
      xold = x;
      yold = y;
      if (fLine==0) {
         fLine = 1;
         if (prvy==1) {
            xprv = x;
            fXs  = xprv;
            fYs  = yprv;
         } else {
            fXs = x;
            fYs = y;
         }
      }
      if (x!=x2) fEnvelope[x] = y;
   } else {
      prvy = 1;
      xprv = x;
      yprv = y;
      if (fLine==1) {
         fLine = 2;
         fXe   = xold;
         fYe   = yold;
      }
   }
   if (x1==x2) {
      if (y1!=y2) y+=(y1<y2)? +1:-1;
      if (y!=y2) goto l1;
   } else {
      x  += krok;
      fy1 = y2-y1;
      fx1 = x2-x1;
      fx  = x-x1;
      fy  = fy1*fx/fx1;
      y   = (Int_t)(y1+fy);
      if (((x<=x2)&&(x1<x2)) || ((x>=x2)&&(x1>x2))) goto l1;
   }
   return;
}


////////////////////////////////////////////////////////////////////////////////
/// Draws channel mark at the screen coordinates x, y. Width of
/// the mark is w, height is h and the type of the mark is determined by the
/// parameter type.

void TSpectrum2Painter::DrawMarker(Int_t x,Int_t y,Int_t w,Int_t h,Int_t type)
{
   TLine *line=new TLine();
   TEllipse *ellipse=new TEllipse();
   line->SetLineColor(fChanmarkColor);
   line->SetLineWidth(1);
   line->SetLineStyle(kPenStyleSolid);
   ellipse->SetLineColor(fChanmarkColor);
   ellipse->SetLineWidth(1);
   ellipse->SetLineStyle(kPenStyleSolid);
   switch (type) {
      case kChannelMarksStyleDot:
         ellipse->SetX1(gPad->PixeltoX(x));
         ellipse->SetY1(gPad->PixeltoY(y)+1);
         ellipse->SetR1(gPad->PixeltoX(w/2));
         ellipse->SetR2(gPad->PixeltoY(h/2));
         ellipse->SetPhimin(0);
         ellipse->SetPhimax(360);
         ellipse->SetTheta(0);
         ellipse->Paint("");
         break;
      case kChannelMarksStyleCross:
         line->PaintLine(gPad->PixeltoX(x-w/2),gPad->PixeltoY(y)+1,
                         gPad->PixeltoX(x+w/2),gPad->PixeltoY(y)+1);
         line->PaintLine(gPad->PixeltoX(x)    ,gPad->PixeltoY(y-h/2)+1,
                         gPad->PixeltoX(x)    ,gPad->PixeltoY(y+h/2+1)+1);
         break;
      case kChannelMarksStyleStar:
         line->PaintLine(gPad->PixeltoX(x-w/2)  ,gPad->PixeltoY(y)+1,
                         gPad->PixeltoX(x+w/2+1),gPad->PixeltoY(y)+1);
         line->PaintLine(gPad->PixeltoX(x)      ,gPad->PixeltoY(y-h/2)+1,
                         gPad->PixeltoX(x)      ,gPad->PixeltoY(y+h/2+1)+1);
         line->PaintLine(gPad->PixeltoX(x-w/2)  ,gPad->PixeltoY(y-h/2)+1,
                         gPad->PixeltoX(x+w/2+1),gPad->PixeltoY(y+h/2+1)+1);
         line->PaintLine(gPad->PixeltoX(x-w/2)  ,gPad->PixeltoY(y+h/2)+1,
                         gPad->PixeltoX(x+w/2+1),gPad->PixeltoY(y-h/2-1)+1);
         break;
      case kChannelMarksStyleRectangle:
         line->PaintLine(gPad->PixeltoX(x-w/2),gPad->PixeltoY(y-h/2)+1,
                         gPad->PixeltoX(x-w/2),gPad->PixeltoY(y+h/2)+1);
         line->PaintLine(gPad->PixeltoX(x-w/2),gPad->PixeltoY(y+h/2)+1,
                         gPad->PixeltoX(x+w/2),gPad->PixeltoY(y+h/2)+1);
         line->PaintLine(gPad->PixeltoX(x+w/2),gPad->PixeltoY(y+h/2)+1,
                         gPad->PixeltoX(x+w/2),gPad->PixeltoY(y-h/2)+1);
         line->PaintLine(gPad->PixeltoX(x+w/2),gPad->PixeltoY(y-h/2)+1,
                         gPad->PixeltoX(x-w/2),gPad->PixeltoY(y-h/2)+1);
         break;
      case kChannelMarksStyleX:
         line->PaintLine(gPad->PixeltoX(x-w/2)  ,gPad->PixeltoY(y-h/2)+1,
                         gPad->PixeltoX(x+w/2+1),gPad->PixeltoY(y+h/2+1)+1);
         line->PaintLine(gPad->PixeltoX(x-w/2)  ,gPad->PixeltoY(y+h/2)+1,
                         gPad->PixeltoX(x+w/2+1),gPad->PixeltoY(y-h/2-1)+1);
         break;
      case kChannelMarksStyleDiamond:
         line->PaintLine(gPad->PixeltoX(x)    ,gPad->PixeltoY(y-h/2)+1,
                         gPad->PixeltoX(x-w/2),gPad->PixeltoY(y)+1);
         line->PaintLine(gPad->PixeltoX(x-w/2),gPad->PixeltoY(y)+1,
                         gPad->PixeltoX(x)    ,gPad->PixeltoY(y+h/2)+1);
         line->PaintLine(gPad->PixeltoX(x)    ,gPad->PixeltoY(y+h/2)+1,
                         gPad->PixeltoX(x+w/2),gPad->PixeltoY(y)+1);
         line->PaintLine(gPad->PixeltoX(x+w/2),gPad->PixeltoY(y)+1,
                         gPad->PixeltoX(x)    ,gPad->PixeltoY(y-h/2)+1);
         break;
      case kChannelMarksStyleTriangle:
         line->PaintLine(gPad->PixeltoX(x)    ,gPad->PixeltoY(y-h/2)+1,
                         gPad->PixeltoX(x-w/2),gPad->PixeltoY(y+h/2)+1);
         line->PaintLine(gPad->PixeltoX(x-w/2),gPad->PixeltoY(y+h/2)+1,
                         gPad->PixeltoX(x+w/2),gPad->PixeltoY(y+h/2)+1);
         line->PaintLine(gPad->PixeltoX(x+w/2),gPad->PixeltoY(y+h/2)+1,
                         gPad->PixeltoX(x)    ,gPad->PixeltoY(y-h/2)+1);
         break;
   }
   delete line;
   delete ellipse;
   return;
}


////////////////////////////////////////////////////////////////////////////////
/// Calculates screen coordinates of the line given by two
/// nodes for contours display mode. The line is given by two points
/// xr, yr, xs, ys. Finally it draws the line.

void TSpectrum2Painter::Slice(Double_t xr, Double_t yr, Double_t xs,
                              Double_t ys, TLine *line)
{
   Int_t krok,xi,yi,xj,yj,a,b,as,bs,pr,ae,be;
   Double_t fx,fy,fx1,fy1;
   xi = (Int_t)(fTxx*(xr-fXmin)/fKx+fTxy*(yr-fYmin)/fKy+fVx);
   xj = (Int_t)(fTxx*(xs-fXmin)/fKx+fTxy*(ys-fYmin)/fKy+fVx);
   yi = (Int_t)(fTyx*(xr-fXmin)/fKx+fTyy*(yr-fYmin)/fKy+fTyz*fZ+fVy);
   yj = (Int_t)(fTyx*(xs-fXmin)/fKx+fTyy*(ys-fYmin)/fKy+fTyz*fZ+fVy);
   as = xi;
   bs = yi;
   ae = xj;
   be = yj;
   a  = xi;
   b  = yi;
   pr = 0;
   krok = (xi<xj)? 1:-1;
l1:
   if (b<=fEnvelope[a]) {
      fEnvelopeContour[a] = b;
      if (pr==0) {
         pr = 1;
         as = a;
         bs = b;
      }
   } else {
      if (pr==1) {
         pr = 2;
         ae = a;
         be = b;
      }
   }
   if (xi==xj) {
      if (yi!=yj) b += (yi<yj)? +1:-1;
      if (b!=yj)  goto l1;
   } else {
      a  += krok;
      fy1 = yj-yi;
      fx1 = xj-xi;
      fx  = a-xi;
      fy  = fy1*fx/fx1;
      b   = (Int_t)(yi+fy);
      if (a!=xj) goto l1;
   }
   if (pr!=0) {
      if (pr==1) {
         ae = xj;
         be = yj;
      }
      line->PaintLine(gPad->PixeltoX(as),gPad->PixeltoY(bs)+1,
                      gPad->PixeltoX(ae),gPad->PixeltoY(be)+1);
   }
   return;
}


////////////////////////////////////////////////////////////////////////////////
/// Copies envelope vector, which ensures hidden surface removal for the
/// contours display mode.

void TSpectrum2Painter::CopyEnvelope(Double_t xr, Double_t xs, Double_t yr,
                                     Double_t ys)
{
   Int_t xi,xj,a;
   xi = (Int_t)(fTxx*(xr-fXmin)/fKx+fTxy*(yr-fYmin)/fKy+fVx);
   xj = (Int_t)(fTxx*(xs-fXmin)/fKx+fTxy*(ys-fYmin)/fKy+fVx);
   if (xi<xj) {
      for (a=xi;a<=xj;a++) {
         if (fEnvelopeContour[a]<fEnvelope[a])
            fEnvelope[a] = fEnvelopeContour[a];
         fEnvelopeContour[a] = fBy2;
      }
   } else if (xj<xi) {
      for (a=xj;a<=xi;a++) {
         if (fEnvelopeContour[a]<fEnvelope[a])
            fEnvelope[a] = fEnvelopeContour[a];
         fEnvelopeContour[a] = fBy2;
      }
   }
   return;
}


////////////////////////////////////////////////////////////////////////////////
/// Paints histogram according to preset parameters.
/// ### Visualization
/// #### Goal: to present 2-dimensional spectra in suitable visual form
/// This package has several display mode groups and display modes, which can be
/// employed for the presentation of 2-dimensional histograms
/// #### Display modes groups:
///
///   - `kModeGroupSimple`      - it covers simple display modes using one
///                                       color only
///   - `kModeGroupLight`       - in this group the shading is carried out
///                                       according to the position of the fictive
///                                       light source
///   - `kModeGroupHeight`      - in this group the shading is carried out
///                                       according to the channel contents
///   - `kModeGroupLightHeight` - combination of two previous shading
///                                       algorithms. One can control the weight
///                                       between both algorithms.
///
/// #### Display modes:
///
///   - `kDisplayModePoints, `
///   - `kDisplayModeGrid, `
///   - `kDisplayModeContours,`
///   - `kDisplayModeBars,`
///   - `kDisplayModeLinesX,`
///   - `kDisplayModeLinesY,`
///   - `kDisplayModeBarsX,`
///   - `kDisplayModeBarsY,`
///   - `kDisplayModeNeedles,`
///   - `kDisplayModeSurface,`
///   - `kDisplayModeTriangles.`
///
/// one can combine the above given modes groups and display modes. The meaningful
/// combinations (denoted by x) are given in the next table.
///
/// |           | Simple | Light | Height | Light-Height |
/// |-----------|--------|-------|--------|--------------|
/// | Points    |    X   |   X   |    X   |       X      |
/// | Grid      |    X   |   X   |    X   |       X      |
/// | Contours  |    X   |   -   |    X   |       -      |
/// | Bars      |    X   |   -   |    X   |       -      |
/// | LinesX    |    X   |   X   |    X   |       X      |
/// | LinesY    |    X   |   X   |    X   |       X      |
/// | BarsX     |    X   |   -   |    X   |       -      |
/// | BarsY     |    X   |   -   |    X   |       -      |
/// | Needles   |    X   |   -   |    -   |       -      |
/// | Surface   |    -   |   X   |    X   |       X      |
/// | Triangles |    X   |   X   |    X   |       X      |
///
/// #### Function: void TSpectrum2Painter::SetDisplayMode (Int_t modeGroup, Int_t displayMode)
///
/// This function controls the display mode group and display mode of the
/// histogram drawing. To illustrate the possible effects of the various display
/// modes we introduce a set of examples. Default values:
///
///   - `modeGroup = kModeGroupLightHeight `
///   - `displayMode = kDisplayModeSurface `
///
/// \image html spectrumpainter001.jpg
///
/// Simple modes group, display mode = points, 256 x 256 channels.
/// \image html spectrumpainter002.jpg
///
/// Simple modes group, display mode = grid, 64 x 64 channels.
/// \image html spectrumpainter003.jpg
///
/// Simple modes group, display mode = contours, 64 x 64 channels.
/// \image html spectrumpainter004.jpg
///
/// Simple modes group, display mode = bars, 64 x 64 channels.
/// \image html spectrumpainter005.jpg
///
/// Simple modes group, display mode = linesX, 64 x 64 channels.
/// \image html spectrumpainter006.jpg
///
/// Simple modes group, display mode = linesY, 64 x 64 channels.
/// \image html spectrumpainter007.jpg
///
/// Simple modes group, display mode = barsX, 64 x 64 channels.
/// \image html spectrumpainter008.jpg
///
/// Simple modes group, display mode = barsY, 64 x 64 channels.
/// \image html spectrumpainter009.jpg
///
/// Simple modes group, display mode = needles, 64 x 64 channels.
/// \image html spectrumpainter010.jpg
///
/// Simple modes group, display mode = triangles, 64 x 64 channels.
/// \image html spectrumpainter011.jpg
///
/// Light modes group, display mode = points, 256 x 256 channels.
/// \image html spectrumpainter012.jpg
///
/// Light modes group, display mode = grid, 256 x 256 channels.
/// \image html spectrumpainter013.jpg
///
/// Light modes group, display mode = surface, 64 x 64 channels.
/// \image html spectrumpainter014.jpg
///
/// Light modes group, display mode = triangles, 64 x 64 channels.
/// \image html spectrumpainter015.jpg
///
/// Height modes group, display mode = points, 256 x 256 channels.
/// \image html spectrumpainter016.jpg
///
/// Height modes group, display mode = grid, 256 x 256 channels.
/// \image html spectrumpainter017.jpg
///
/// Height modes group, display mode = contours, 64 x 64 channels.
/// \image html spectrumpainter018.jpg
///
/// Height modes group, display mode = bars, 64 x 64 channels.
/// \image html spectrumpainter019.jpg
///
/// Height modes group, display mode = surface, 64 x 64 channels.
/// \image html spectrumpainter020.jpg
///
/// Height modes group, display mode = triangles, 64 x 64 channels.
/// \image html spectrumpainter021.jpg
///
/// Light - height modes group, display mode = surface, 64 x 64 channels. The weight
/// between both shading algorithms is set to 0.5. One can observe the influence of
/// both shadings.
///
/// #### Function: TSpectrum2Painter::SetPenAttr(Int_t color,Int_t style,Int_t width)
///
/// Using this function one can change pen color, pen style and pen width.
/// Possible pen styles are:
///
///   - ` kPenStyleSolid,`
///   - ` kPenStyleDash,`
///   - ` kPenStyleDot,`
///   - ` kPenStyleDashDot.`
///
/// Default values:
///
///   - ` color = kBlack`
///   - ` style = kPenStyleSolid`
///   - ` width = 1`
///
/// \image html spectrumpainter022.jpg
///
/// Simple modes group, display mode = linesX, 64 x 64 channels. Pen width = 3.
///
/// #### Function: TSpectrum2Painter::SetNodes(Int_t nodesx,Int_t nodesy)
///
/// Sometimes the displayed region is rather large. When displaying all channels
/// pictures become very dense and complicated. It is very difficult to understand
/// overall shape of the data. Therefore in the package we have implemented the
/// possibility to change the density of displayed channels. Only channels
/// coinciding with given nodes are displayed. In the next figure we introduce the
/// example of the above presented spectrum with number of nodes set to 64x64.
///
/// Default values:
///
///   - ` nodesx = Xmax-Xmin+1`
///   - ` nodesy = Ymax-Ymin+1`
///
/// \image html spectrumpainter023.jpg
///
/// Simple modes group, display mode = grid, 256 x 256 channels.
/// Number of nodes is 64x64.
///
/// #### Function: void TSpectrum2Painter::SetAngles (Int_t alpha,Int_t beta, Int_t view)
///
/// One can change the angles of the position of 3-d space and to rotate the
/// space. Alpha parameter defines the angle between bottom horizontal screen line
/// and the displayed space on the right side of the picture and beta on the left
/// side, respectively. One can rotate the 3-d space around vertical axis going
/// through the center of it employing the view parameter.
/// Allowed values are 0, 90, 180 and 270 degrees respectively.
///
/// Default values:
///
///   - ` alpha = 20`
///   - ` beta = 60`
///   - ` view = 0`
///
/// \image html spectrumpainter024.jpg
///
/// Light modes group, display mode = surface, 256 x 256 channels. Angles are
/// set as follows: alpha=40, beta=30, view=0.
/// \image html spectrumpainter025.jpg
///
/// Light modes group, display mode = surface, 256 x 256 channels. Angles are
/// set as follows: alpha=30, beta=30, view=90.
///
/// #### Function: TSpectrum2Painter::SetZScale(Int_t scale)
///
/// One can change the scale of z-axis. Possible values are:
///
///   - ` kZScaleLinear`
///   - ` kZScaleLog`
///   - ` kZScaleSqrt`
///
/// Default value is:
///
///   - ` scale = kZScaleLinear`
///
/// \image html spectrumpainter026.jpg
///
/// Height modes group, display mode = surface, 64 x 64 channels, log scale.
///
/// #### Function: TSpectrum2Painter::SetColorIncrements(Double_t r,Double_t g,Double_t b);
///
/// For sophisticated shading (in kModeGroupLight, kModeGroupHeight
/// and kModeGroupLightHeight display modes groups) the color palette starts
/// from the basic pen color (see SetPenAttr function). There is a predefined number
/// of color levels (256). Color in every level is calculated by adding the
/// increments of the r, g, b components to the previous level. Using this function
/// one can change the color increments between two neighbouring color levels. The
/// function does not apply for kModeGroupSimple display modes group.
/// Default values: r=1, g=1, b=1;
/// \image html spectrumpainter027.jpg
///
/// Light modes group, display mode = surface, 64 x 64 channels, color increments
/// r=1, g=2, b=3.
/// \image html spectrumpainter028.jpg
///
/// Light modes group, display mode = surface, 64 x 64 channels, color
/// increments r=4, g=2, b=1.
///
/// #### Function: TSpectrum2Painter::SetColorAlgorithm(Int_t colorAlgorithm)
///
/// To define the colors one can employ one of the following color algorithms
/// (rgb, cmy, cie, yiq, hvs models [1], [2]). When the level of a component
/// achieves the limit value one can choose either smooth transition (by decreasing
/// the limit value) or sharp - modulo transition (continuing with 0 value). This
/// makes possible to realize various visual effects. One can choose from the
/// following set of the algorithms:
///
///   - ` kColorAlgRgbSmooth `
///   - ` kColorAlgRgbModulo `
///   - ` kColorAlgCmySmooth `
///   - ` kColorAlgCmyModulo `
///   - ` kColorAlgCieSmooth `
///   - ` kColorAlgCieModulo `
///   - ` kColorAlgYiqSmooth `
///   - ` kColorAlgYiqModulo `
///   - ` kColorAlgHvsSmooth `
///   - ` kColorAlgHvsModulo `
///
/// The function does not apply for kModeGroupSimple display modes group.
/// Default value is:
///
///   - ` colorAlgorithm = kColorAlgRgbSmooth`
///
/// \image html spectrumpainter029.jpg
///
/// Light modes group, display mode = surface, 64 x 64 channels, color algorithm
/// is cmy smooth.
/// \image html spectrumpainter030.jpg
///
/// Light modes group, display mode = surface, 64 x 64 channels, color algorithm
/// is hvs smooth.
/// \image html spectrumpainter031.jpg
///
/// Light modes group, display mode = surface, 64 x 64 channels, color algorithm
/// is yiq smooth.
/// \image html spectrumpainter032.jpg
///
/// Light modes group, display mode = surface, 64 x 64 channels, color algorithm
/// is rgb modulo.
/// \image html spectrumpainter033.jpg
///
/// Height modes group, display mode = surface, 256 x 256 channels, color
/// algorithm is rgb modulo, increments r=5, g=5, b=5, angles alpha=0, beta=90,
/// view=0.
///
/// #### Function: TSpectrum2Painter::SetLightPosition(Int_t x, Int_t y, Int_t z)
///
/// In kModeGroupLight and kModeGroupLightHeight display modes
/// groups the color palette is calculated according to the fictive light source
/// position in 3-d space. Using this function one can change the position of the
/// source and thus to achieve various graphical effects. The function does not
/// apply for kModeGroupSimple and kModeGroupHeight display modes
/// groups. Default values are: x=1000, y=1000, z=1000.
/// \image html spectrumpainter034.jpg
///
/// Light modes group, display mode = surface, 64 x 64 channels. Position of the
/// light source was set to x=0, y=1000, z=1000.
///
/// #### Function: TSpectrum2Painter::SetShading(Int_t shading,Int_t shadow)
///
/// Surface of the picture is composed of triangles. If desired the edges of the
/// neighbouring triangles can be smoothed (shaded). If desired the display of the
/// shadow can be painted as well. The function does not apply for
/// kModeGroupSimple display modes group.
///
/// Possible values for shading are:
///
///   - ` kNotShaded`
///   - ` kShaded.`
///
/// Possible values for shadow are:
///
///   - ` kShadowsNotPainted`
///   - ` kShadowsPainted`
///
/// Default values:
///
///   - ` shading = kShaded`
///   - ` shadow = kShadowsNotPainted`
///
/// \image html spectrumpainter035.jpg
///
/// Light modes group, display mode = surface, 64 x 64 channels, not shaded.
/// \image html spectrumpainter036.jpg
///
/// Light modes group, display mode = surface, 64 x 64 channels, shaded, with
/// shadow.
///
/// #### Function: TSpectrum2Painter::SetBezier(Int_t bezier)
///
/// For kModeGroupSimple display modes group and for kDisplayModeGrid,
/// kDisplayModeLinesX >and kDisplayModeLinesY display modes one
/// can smooth data using Bezier smoothing algorithm. The function does not apply
/// for other display modes groups and display modes. Possible values are:
///
///   - ` kNoBezierInterpol`
///   - ` kBezierInterpol`
///
/// Default value is:
///
///   - ` bezier = kNoBezierInterpol.`
///
/// \image html spectrumpainter005.jpg
///
/// Simple modes group, display mode = linesX, 64 x 64 channels with Bezier
/// smoothing.
///
/// #### Function: TSpectrum2Painter::SetContourWidth(Int_t width)
///
/// This function applies only for kDisplayModeContours display mode.
/// One can change the width between horizontal slices and thus their density.
/// Default value: width=50.
/// \image html spectrumpainter037.jpg
///
/// Simple modes group, display mode = contours, 64 x 64 channels. Width between
/// slices was set to 30.
///
/// #### Function: TSpectrum2Painter::SetLightHeightWeight(Double_t weight)
///
/// For kModeGroupLightHeight display modes group one can change the
/// weight between both shading algorithm. The function does not apply for other
/// display modes groups. Default value is: weight=0.5.
/// \image html spectrumpainter038.jpg
///
/// Light - height modes group, display mode = surface, 64 x 64 channels.
/// The weight between both shading algorithms is set to 0.7.
///
/// #### Function: TSpectrum2Painter::SetChanMarks(Int_t enable,Int_t color,Int_t width,Int_t height,Int_t style)
/// In addition to the surface drawn using any above given algorithm one can display
/// channel marks. One can control the color as well as the width, height
/// (in pixels) and the style of the marks. The parameter enable can be set to:
///
///   - `kChannelMarksNotDrawn`
///   - `kChannelMarksDrawn.`
///
/// The possible styles can be chosen from the set:
///
///   - ` kChannelMarksStyleDot`
///   - ` kChannelMarksStyleCross`
///   - ` kChannelMarksStyleStar`
///   - ` kChannelMarksStyleRectangle`
///   - ` kChannelMarksStyleX`
///   - ` kChannelMarksStyleDiamond`
///   - ` kChannelMarksStyleTriangle.`
///
/// \image html spectrumpainter039.jpg
///
/// Light modes group, display mode = surface, 64 x 64 channels,
/// with marks (red circles).</p>
///
/// #### Function: TSpectrum2Painter::SetChanGrid(Int_t enable,Int_t color)
///
/// In addition to the surface drawn using any above given algorithm one can
/// display grid using the color parameter. The parameter enable can be set to:
///
///   - ` kChannelGridNotDrawn`
///   - ` kChannelGridDrawn.`
///
/// \image html spectrumpainter040.jpg
///
/// Height modes group, display mode = surface, 64 x 64 channels, with blue grid.
/// \image html spectrumpainter041.jpg
///
/// Height modes group, display mode = surface, 64 x 64 channels, with marks
/// (red circles) and blue grid.
/// #### References:
///
///  [1] Morhac M., Kliman J., Matouoek V., Turzo I.,
/// Sophisticated visualization algorithms for analysis of multidimensional
/// experimental nuclear data, Acta Physica Slovaca 54 (2004) 385.
///
///  [2] D. Hearn, M. P. Baker: Computer Graphics, Prentice Hall International,
/// Inc.  1994.
/// #### Script:
///
/// Example to draw source spectrum (class TSpectrum2Painter).
/// To execute this example, do:
/// ~~~
/// root > .x VisA.C
/// ~~~
/// ~~~ {.cpp}
/// #include "TSpectrum2Painter.h"
///
/// void VisA() {
///    TFile *f = new TFile("TSpectrum2.root");
///    TH2F *graph=(TH2F*) f->Get("graph2;1");
///    TCanvas *Graph2 = new TCanvas("Graph2","Illustration of 2D graphics",10,10,1000,700);
///    graph->Draw("SPEC");
/// }
/// ~~~

void TSpectrum2Painter::Paint(Option_t * /*option*/)
{


   Int_t turni,turnj,w1,w2,x,y;
   Int_t q1=0,q2=0,qv=0,smer=0,flag=0,i=0,j=0,x1=0,y1=0,x2=0,y2=0,x3=0,y3=0,x4=0,y4=0,uhl=0,xp1=0,yp1=0,xp2=0,yp2=0;
   Int_t ix5,iy5,x6,y6,x7,y7,y8,x1d,y1d,x2d=0,y2d=0;
   Int_t i1=0,i2=0,i3=0,i4=0,j1=0,j2=0,j3=0,j4=0;
   Int_t s1=0,s2=0,s3=0,s4=0,t1=0,t2=0,t3=0,t4=0;
   Double_t dx1,dx2,dx3,dx4,dy1,dy2,dy3,dy4,z1,z2,z3,z4,zl,zh;
   Double_t xa,xb=0,ya,yb=0,x5=0,y5=0;
   Double_t da=0,db=0,dc=0,dd=0,xtaz,ytaz,ztaz,v,shad_noise;
   Int_t iv=0,ekv,stvor,sx1,sx2,sx3,sx4,sx5,sy1,sy2,sy3,sy4,sy5;
   Double_t pom1,pom2,sdx1,sdy1,sdx2=0,sdy2,sdx3,sdy3,sdy4,spriz;
   Int_t sr1=0,sr2=0,sr3=0,sr4=0,sr5=0,sr6=0,sr7=0,sr8=0;
   Int_t tr1=0,tr2=0,tr3=0,tr4=0,tr5=0,tr6=0,tr7=0,tr8=0;
   Int_t il,iv1=0,iv2=0,iv3=0,iv4=0;
   Double_t v1=0,v2=0,v3=0,v4=0,dxr1,dxr2,dyr1,dyr2,zr1,zr2,bezf;
   Double_t dcount_reg,z1l,z2l,z3l,z4l,sdx2p,sdy2p,dap,dbp,dcp,ddp;
   Int_t sx1p,sy1p,sx3p,uip=0;
   Double_t bezx1,bezy1,bezx2,bezy2;
   Double_t p000x,p000y,p100x,p100y,p010x,p010y,p110x,p110y;
   Double_t p001x,p001y,p101x,p101y,p011x,p011y,p111x,p111y;
   Int_t ibezx1=0,ibezy1=0,ibezx2,ibezy2;
   unsigned ui1,ui2,ui3;
   Double_t fi,alfa,beta,x3max,y3max,mul,movx,movy;
   Double_t xmin,xmax,ymin,ymax,zmin,zmax,mx,my,mz;
   Double_t mxx,mxy,myx,myy,myz,px,py,kx,ky;
   Double_t bxl,bxh,byl,byh,xd,yd,a,b,rotx,roty;
   TLine  *line = new TLine();
   TBox   *box  = new TBox();
   TColor *pen_col;
   pen_col = (TColor*)(gROOT->GetListOfColors()->At(fPenColor));
   ui1 = (Int_t)(256*pen_col->GetRed());
   ui2 = (Int_t)(256*pen_col->GetGreen());
   ui3 = (Int_t)(256*pen_col->GetBlue());

   if (fBx2>=fMaximumXScreenResolution) {
      printf("The canvas size exceed the maximum X screen resolution.\n");
      printf("Use the option bf() to increase the buffer size (it should be greater than %d).\n",fBx2);
      return;
   }

   for (i=fBx1;i<fBx2;i++) {
      fEnvelope[i]        = fBy2;
      fEnvelopeContour[i] = fBy2;
   }

//   gPad->Range(0, 0, 1 ,1);

   // Set the histogram's parameters.
   fBx1      = gPad->XtoPixel(0.1);
   fBx2      = gPad->XtoPixel(0.99);
   fBy1      = gPad->YtoPixel(0.99);
   fBy2      = gPad->YtoPixel(0.05);
   fXmin     = fH2->GetXaxis()->GetFirst();
   fXmax     = fH2->GetXaxis()->GetLast();
   fYmin     = fH2->GetYaxis()->GetFirst();
   fYmax     = fH2->GetYaxis()->GetLast();
   fZmax     = fH2->GetMaximum();
   fZmin     = fH2->GetMinimum();

   // Calculation of display parameters.
   xmin  = fXmin;
   xmax  = fXmax;
   ymin  = fYmin;
   ymax  = fYmax;
   zmin  = fZmin;
   zmax  = fZmax;
   xd    = (xmax-xmin)/2;
   yd    = (ymax-ymin)/2;
   a     = (xmax+xmin)/2;
   b     = (ymax+ymin)/2;
   fi    = (fViewAngle*3.1415927)/180;
   alfa  = (fAlpha*3.1415927)/180;
   beta  = (fBeta*3.1415927)/180;
   rotx  = (-1)*a*cos(fi)+b*sin(fi)+xd*TMath::Abs(cos(fi))+yd*TMath::Abs(sin(fi));
   roty  = (-1)*a*sin(fi)-b*cos(fi)+xd*TMath::Abs(sin(fi))+yd*TMath::Abs(cos(fi));
   x3max = (xmax-xmin)*TMath::Abs(cos(fi))+(ymax-ymin)*TMath::Abs(sin(fi));
   y3max = (xmax-xmin)*TMath::Abs(sin(fi))+(ymax-ymin)*TMath::Abs(cos(fi));
   bxl   = fBx1;
   bxh   = fBx2;
   byl   = fBy1;
   byh   = fBy2;
   mx    = (bxh-bxl)/(x3max*(cos(alfa)+cos(beta)));
   my    = (bxh-bxl)/(y3max*(cos(alfa)+cos(beta)));
   mul   = (byh-byl)/(bxh-bxl);
   movx  = bxl+my*cos(alfa)*y3max;
   mxx   = mx*cos(beta)*cos(fi)-my*cos(alfa)*sin(fi);
   mxy   = (-1)*mx*cos(beta)*sin(fi)-my*cos(alfa)*cos(fi);
   myx   = mul*(mx*sin(beta)*cos(fi)+my*sin(alfa)*sin(fi));
   myy   = mul*((-1)*mx*sin(beta)*sin(fi)+my*sin(alfa)*cos(fi));
   px    = rotx*mx*cos(beta)-roty*my*cos(alfa)+movx;
   kx    = (xmax-xmin)/(fNodesx-1);
   ky    = (ymax-ymin)/(fNodesy-1);
   fKx   = kx;
   fKy   = ky;
   fMxx  = mxx;
   fMxy  = mxy;
   fMyx  = myx;
   fMyy  = myy;
   fTxx  = mxx*kx;
   fTxy  = mxy*ky;
   fTyx  = myx*kx;
   fTyy  = myy*ky;
   fVx   = mxx*xmin+mxy*ymin+px;
   if (fZscale==kZScaleLinear) {
      mz     = (bxh-bxl)*(cos(alfa)+cos(beta)-sin(alfa)-sin(beta));
      mz     = mz/((zmax-zmin)*(cos(alfa)+cos(beta)));
      movy   = byl+mul*mz*zmax;
      myz    = (-1)*mz*mul;
      py     = mul*(rotx*mx*sin(beta)+roty*my*sin(alfa))+movy;
      fTyz   = myz;
      fVy    = myx*xmin+myy*ymin+py;
      fNuSli = (zmax-zmin)/(Double_t)fContWidth;
   } else if (fZscale==kZScaleLog) {
      if (zmin>=1) zmin = log(zmin);
      else         zmin = 0;
      if (zmax>=1) zmax = log(zmax);
      else         zmax = 0;
      if ((zmax-zmin)<0.000001) zmax = zmin+0.000001;
      mz     = (bxh-bxl)*(cos(alfa)+cos(beta)-sin(alfa)-sin(beta));
      mz     = mz/((zmax-zmin)*(cos(alfa)+cos(beta)));
      movy   = byl+mul*mz*zmax;
      myz    = (-1)*mz*mul;
      py     = mul*(rotx*mx*sin(beta)+roty*my*sin(alfa))+movy;
      fTyz   = myz;
      fVy    = myx*xmin+myy*ymin+py;
      fNuSli = (zmax-zmin)/(Double_t)fContWidth;
   } else if (fZscale==kZScaleSqrt) {
      if (zmin>=1) zmin = sqrt(zmin);
      else         zmin = 0;
      if (zmax>=1) zmax = sqrt(zmax);
      else         zmax = 0;
      if ((zmax-zmin)<0.000001) zmax = zmin+0.000001;
      mz     = (bxh-bxl)*(cos(alfa)+cos(beta)-sin(alfa)-sin(beta));
      mz     = mz/((zmax-zmin)*(cos(alfa)+cos(beta)));
      movy   = byl+mul*mz*zmax;
      myz    = (-1)*mz*mul;
      py     = mul*(rotx*mx*sin(beta)+roty*my*sin(alfa))+movy;
      fTyz   = myz;
      fVy    = myx*xmin+myy*ymin+py;
      fNuSli = (zmax-zmin)/(Double_t)fContWidth;
   }

   // End of calculations of display parameters.
   dcount_reg=fContWidth;
   switch (fZscale) {
      case kZScaleLog:
         dcount_reg=log(dcount_reg);
         break;
      case kZScaleSqrt:
         dcount_reg=sqrt(dcount_reg);
         break;
   }
   shad_noise  = fZmax;
   shad_noise /= 100.;
   w1          = fNodesx-1;
   w2          = fNodesy-1;

   // Drawing axis in backplanes.
   Transform(0,0,-1);
   p000x = gPad->PixeltoX(fXt);
   p000y = gPad->PixeltoY(fYt)+1;
   Transform(w1,0,-1);
   p100x = gPad->PixeltoX(fXt);
   p100y = gPad->PixeltoY(fYt)+1;
   Transform(0,w2,-1);
   p010x = gPad->PixeltoX(fXt);
   p010y = gPad->PixeltoY(fYt)+1;
   Transform(w1,w2,-1);
   p110x = gPad->PixeltoX(fXt);
   p110y = gPad->PixeltoY(fYt)+1;
   fZPresetValue = fZmax;
   Transform(0,0,-2);
   p001x = gPad->PixeltoX(fXt);
   p001y = gPad->PixeltoY(fYt)+1;
   Transform(w1,0,-2);
   p101x = gPad->PixeltoX(fXt);
   p101y = gPad->PixeltoY(fYt)+1;
   Transform(0,w2,-2);
   p011x = gPad->PixeltoX(fXt);
   p011y = gPad->PixeltoY(fYt)+1;
   Transform(w1,w2,-2);
   p111x = gPad->PixeltoX(fXt);
   p111y = gPad->PixeltoY(fYt)+1;
   Double_t bmin, bmax, binLow, binHigh, binWidth;
   Double_t axisLevel, gridDist, gridY1, gridY2;
   Int_t ndivx = 0, ndivy, ndivz, nbins;
   TGaxis *axis  = new TGaxis();
   TGaxis *xaxis = new TGaxis();
   TGaxis *yaxis = new TGaxis();
   TGaxis *zaxis = new TGaxis();
   line->SetLineStyle(kPenStyleDot);
   if (fViewAngle==0) {
      axis->PaintAxis(p000x, p000y, p100x, p100y, bmin, bmax, ndivx, "");
      axis->PaintAxis(p000x, p000y, p010x, p010y, bmin, bmax, ndivx, "");
      if(fAlpha+fBeta<90)
         axis->PaintAxis(p000x, p000y, p001x, p001y, bmin, bmax, ndivx, "");
      if(fAlpha+fBeta<90)
         axis->PaintAxis(p100x, p100y, p101x, p101y, bmin, bmax, ndivx, "");
      axis->PaintAxis(p101x, p101y, p001x, p001y, bmin, bmax, ndivx, "");
      axis->PaintAxis(p001x, p001y, p011x, p011y, bmin, bmax, ndivx, "");
      if (fZscale==kZScaleLinear) {
         bmin  = fZmin;
         bmax  = fZmax;
         ndivz = 10;
         THLimitsFinder::Optimize(bmin, bmax, ndivz, binLow, binHigh,
                                  nbins, binWidth, " ");
         for (i = 0; i < nbins + 1; i++) {
            axisLevel = binLow+i*binWidth;
            gridDist  = (axisLevel-bmin)*(p001y-p000y)/(bmax-bmin);
            gridY1    = p000y + gridDist, gridY2 = p100y + gridDist;
            line->PaintLine(p000x,gridY1,p100x,gridY2);
            gridY2    = p010y + gridDist;
            line->PaintLine(p000x,gridY1,p010x,gridY2);
         }
      }
   } else if (fViewAngle==90) {
      axis->PaintAxis(p010x, p010y, p000x, p000y, bmin, bmax, ndivx, "");
      axis->PaintAxis(p010x, p010y, p110x, p110y, bmin, bmax, ndivx, "");
      if(fAlpha+fBeta<90)
         axis->PaintAxis(p010x, p010y, p011x, p011y, bmin, bmax, ndivx, "");
      if(fAlpha+fBeta<90)
         axis->PaintAxis(p000x, p000y, p001x, p001y, bmin, bmax, ndivx, "");
      axis->PaintAxis(p001x, p001y, p011x, p011y, bmin, bmax, ndivx, "");
      axis->PaintAxis(p011x, p011y, p111x, p111y, bmin, bmax, ndivx, "");
      if (fZscale==kZScaleLinear) {
         bmin  = fZmin;
         bmax  = fZmax;
         ndivz = 10;
         THLimitsFinder::Optimize(bmin, bmax, ndivz, binLow, binHigh,
                                  nbins, binWidth, " ");
         for (i = 0; i < nbins + 1; i++) {
            axisLevel = binLow+i*binWidth;
            gridDist  = (axisLevel-bmin)*(p011y-p010y)/(bmax-bmin);
            gridY1    = p010y + gridDist, gridY2 = p000y + gridDist;
            line->PaintLine(p010x,gridY1,p000x,gridY2);
            gridY2    = p110y + gridDist;
            line->PaintLine(p010x,gridY1,p110x,gridY2);
         }
      }
   } else if (fViewAngle==180) {
      axis->PaintAxis(p110x, p110y, p010x, p010y, bmin, bmax, ndivx, "");
      axis->PaintAxis(p110x, p110y, p100x, p100y, bmin, bmax, ndivx, "");
      if(fAlpha+fBeta<90)
         axis->PaintAxis(p110x, p110y, p111x, p111y, bmin, bmax, ndivx, "");
      if(fAlpha+fBeta<90)
         axis->PaintAxis(p010x, p010y, p011x, p011y, bmin, bmax, ndivx, "");
      axis->PaintAxis(p011x, p011y, p111x, p111y, bmin, bmax, ndivx, "");
      axis->PaintAxis(p111x, p111y, p101x, p101y, bmin, bmax, ndivx, "");
      if (fZscale==kZScaleLinear) {
         bmin  = fZmin;
         bmax  = fZmax;
         ndivz = 10;
         THLimitsFinder::Optimize(bmin, bmax, ndivz, binLow, binHigh,
                                  nbins, binWidth, " ");
         for (i = 0; i < nbins + 1; i++) {
            axisLevel = binLow+i*binWidth;
            gridDist  = (axisLevel-bmin)*(p111y-p110y)/(bmax-bmin);
            gridY1    = p110y + gridDist, gridY2 = p010y + gridDist;
            line->PaintLine(p110x,gridY1,p010x,gridY2);
            gridY2    = p100y + gridDist;
            line->PaintLine(p110x,gridY1,p100x,gridY2);
         }
      }
   } else if (fViewAngle==270) {
      axis->PaintAxis(p100x, p100y, p110x, p110y, bmin, bmax, ndivx, "");
      axis->PaintAxis(p100x, p100y, p000x, p000y, bmin, bmax, ndivx, "");
      if(fAlpha+fBeta<90)
         axis->PaintAxis(p100x, p100y, p101x, p101y, bmin, bmax, ndivx, "");
      if(fAlpha+fBeta<90)
         axis->PaintAxis(p110x, p110y, p111x, p111y, bmin, bmax, ndivx, "");
      axis->PaintAxis(p111x, p111y, p101x, p101y, bmin, bmax, ndivx, "");
      axis->PaintAxis(p101x, p101y, p001x, p001y, bmin, bmax, ndivx, "");
      if (fZscale==kZScaleLinear) {
         bmin  = fZmin;
         bmax  = fZmax;
         ndivz = 10;
         THLimitsFinder::Optimize(bmin, bmax, ndivz, binLow, binHigh,
                                  nbins, binWidth, " ");
         for (i = 0; i < nbins + 1; i++) {
            axisLevel = binLow+i*binWidth;
            gridDist  = (axisLevel-bmin)*(p101y-p100y)/(bmax-bmin);
            gridY1    = p100y + gridDist, gridY2 = p110y + gridDist;
            line->PaintLine(p100x,gridY1,p110x,gridY2);
            gridY2    = p000y + gridDist;
            line->PaintLine(p100x,gridY1,p000x,gridY2);
         }
      }
   }

   // End.
   line->ResetAttLine("");
   line->SetLineColor(fPenColor);
   line->SetLineWidth(fPenWidth);
   line->SetLineStyle(fPenDash);
   turni = 0;
   turnj = 0;
   Transform(w1,0,0);
   x1    = fXt;
   Transform(0,0,0);
   x2    = fXt;
   Transform(0,w2,0);
   x3    = fXt;
   if (x2>=x1) turnj = 1;
   if (x3>=x2) turni = 1;
   q1 = 1;
   q2 = 0;
   qv = 1;
   do {
      uhl  = 0;
      smer = 0;
      flag = 0;
l2:
      if (turni==1) {
         i = q1;
      } else {
         i = w1-q1;
      }
      if (turnj==1) {
         j = q2;
      } else {
         j = w2-q2;
      }
      Transform(i,j,0);
      x1 = fXt;
      y1 = fYt;
      Transform(i,j,-1);
      x1d = fXt;
      y1d = fYt;
      do {
         if (flag==0) {
            flag = 1;
            if (smer==0) q1 -= 1;
            else         q2 -= 1;
         } else {
            flag = 0;
            if (smer==0) q2 += 1;
            else         q1 += 1;
         }
         if (turni==1) {
            i = q1;
         } else {
            i = w1-q1;
         }
         if (turnj==1) {
            j = q2;
         } else {
            j = w2-q2;
         }
         Transform(i,j,0);
         x2 = fXt;
         y2 = fYt;
         if (flag==1) {
            x  = x1;
            y  = y1;
            x1 = x2;
            y1 = y2;
            x2 = x;
            y2 = y;
         }
         switch (fDisplayMode) {
            case kDisplayModePoints:
               if (fModeGroup==kModeGroupSimple) {
                  Envelope(x1,y1,x2,y2);
                  if (y1<=fEnvelope[x1]) {
                     line->PaintLine(gPad->PixeltoX(x1)  ,gPad->PixeltoY(y1)+1,
                                     gPad->PixeltoX(x1+1),gPad->PixeltoY(y1)+1);
                  }
                  if (y2<=fEnvelope[x2]) {
                     line->PaintLine(gPad->PixeltoX(x2)  ,gPad->PixeltoY(y2)+1,
                                     gPad->PixeltoX(x2+1),gPad->PixeltoY(y2)+1);
                  }
               } else {
                  if ((q1!=q2||smer!=0) && flag==1) {
                     s1 = q1+1;
                     t1 = q2;
                     s2 = q1;
                     t2 = q2;
                     s3 = q1;
                     t3 = q2+1;
                     s4 = q1+1;
                     t4 = q2+1;
                     if (fShading==kShaded) {
                        sr1 = s1;
                        tr1 = (Int_t)TMath::Max(t1-1,0);
                        sr2 = s2;
                        tr2 = (Int_t)TMath::Max(t2-1,0);
                        sr3 = (Int_t)TMath::Max(s2-1,0);
                        tr3 = t2;
                        sr4 = (Int_t)TMath::Max(s3-1,0);
                        tr4 = t3;
                        sr5 = s3;
                        tr5 = t3+1;
                        sr6 = s4;
                        tr6 = t4+1;
                        sr7 = s4+1;
                        tr7 = t4;
                        sr8 = s1+1;
                        tr8 = t1;
                     }
                     if (turni==1) {
                        i1 = s1;
                        i2 = s2;
                        i3 = s3;
                        i4 = s4;
                     } else {
                        i1 = (Int_t)TMath::Max(w1-s1,0);
                        i2 = (Int_t)TMath::Max(w1-s2,0);
                        i3 = (Int_t)TMath::Max(w1-s3,0);
                        i4 = (Int_t)TMath::Max(w1-s4,0);
                        if (fShading==kShaded) {
                           sr1 = (Int_t)TMath::Max(w1-sr1,0);
                           sr2 = (Int_t)TMath::Max(w1-sr2,0);
                           sr3 = (Int_t)TMath::Max(w1-sr3,0);
                           sr4 = (Int_t)TMath::Max(w1-sr4,0);
                           sr5 = (Int_t)TMath::Max(w1-sr5,0);
                           sr6 = (Int_t)TMath::Max(w1-sr6,0);
                           sr7 = (Int_t)TMath::Max(w1-sr7,0);
                           sr8 = (Int_t)TMath::Max(w1-sr8,0);
                        }
                     }
                     if (turnj==1) {
                        j1 = t1;
                        j2 = t2;
                        j3 = t3;
                        j4 = t4;
                     } else {
                        j1 = (Int_t)TMath::Max(w2-t1,0);
                        j2 = (Int_t)TMath::Max(w2-t2,0);
                        j3 = (Int_t)TMath::Max(w2-t3,0);
                        j4 = (Int_t)TMath::Max(w2-t4,0);
                        if (fShading==kShaded) {
                           tr1 = (Int_t)TMath::Max(w2-tr1,0);
                           tr2 = (Int_t)TMath::Max(w2-tr2,0);
                           tr3 = (Int_t)TMath::Max(w2-tr3,0);
                           tr4 = (Int_t)TMath::Max(w2-tr4,0);
                           tr5 = (Int_t)TMath::Max(w2-tr5,0);
                           tr6 = (Int_t)TMath::Max(w2-tr6,0);
                           tr7 = (Int_t)TMath::Max(w2-tr7,0);
                           tr8 = (Int_t)TMath::Max(w2-tr8,0);
                        }
                     }
                     Transform(i1,j1,0);
                     x1   = fXt;
                     y1   = fYt;
                     dx1  = fDxspline;
                     dy1  = fDyspline;
                     z1   = fZ;
                     Transform(i2,j2,0);
                     x2   = fXt;
                     y2   = fYt;
                     dx2  = fDxspline;
                     dy2  = fDyspline;
                     z2   = fZ;
                     Transform(i3,j3,0);
                     x3   = fXt;
                     y3   = fYt;
                     dx3  = fDxspline;
                     dy3  = fDyspline;
                     z3   = fZ;
                     Transform(i4,j4,0);
                     x4   = fXt;
                     y4   = fYt;
                     dx4  = fDxspline;
                     dy4  = fDyspline;
                     z4   = fZ;
                     Envelope(x1,y1,x2,y2);
                     Envelope(x2,y2,x3,y3);
                     xtaz = (dx1+dx2+dx4)/3;
                     ytaz = (dy1+dy2+dy4)/3;
                     ztaz = (z1+z2+z4)/3;
                     v    = ColorCalculation(dx1,dy1,z1,dx2,dy2,z2,dx4,dy4,z4);
                     if (fShading==kShaded) {
                        if (fShadow==kShadowsNotPainted) {
                           if (sr1<0||sr1>w1||tr1<0||tr1>w2) Transform(sr1,tr1,-1);
                           else                              Transform(sr1,tr1,0);
                           dxr1  = fDxspline;
                           dyr1  = fDyspline;
                           zr1   = fZ;
                           if (sr8<0||sr8>w1||tr8<0||tr8>w2) Transform(sr8,tr8,-1);
                           else                              Transform(sr8,tr8,0);
                           dxr2  = fDxspline;
                           dyr2  = fDyspline;
                           zr2   = fZ;
                           v     = v+ColorCalculation(dxr1,dyr1,zr1,dx2,dy2,z2,dx1,dy1,z1);
                           v     = v+ColorCalculation(dxr1,dyr1,zr1,dx1,dy1,z1,dxr2,dyr2,zr2);
                           v     = v+ColorCalculation(dxr2,dyr2,zr2,dx1,dy1,z1,dx4,dy4,z4);
                           v1    = v/4;
                           if (sr3<0||sr3>w1||tr3<0||tr3>w2) Transform(sr3,tr3,-1);
                           else                              Transform(sr3,tr3,0);
                           dxr1  = fDxspline;
                           dyr1  = fDyspline;
                           zr1   = fZ;
                           if (sr2<0||sr2>w1||tr2<0||tr2>w2) Transform(sr2,tr2,-1);
                           else                              Transform(sr2,tr2,0);
                           dxr2  = fDxspline;
                           dyr2  = fDyspline;
                           zr2   = fZ;
                           v     = ColorCalculation(dx1,dy1,z1,dx2,dy2,z2,dx3,dy3,z3);
                           v     = v+ColorCalculation(dx2,dy2,z2,dxr1,dyr1,zr1,dx3,dy3,z3);
                           v     = v+ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dxr1,dyr1,zr1);
                           v     = v+ColorCalculation(dxr2,dyr2,zr2,dx2,dy2,z2,dx1,dy1,z1);
                           v2    = v/4;
                           if (sr5<0||sr5>w1||tr5<0||tr5>w2) Transform(sr5,tr5,-1);
                           else                              Transform(sr5,tr5,0);
                           dxr1  = fDxspline;
                           dyr1  = fDyspline;
                           zr1   = fZ;
                           if (sr4<0||sr4>w1||tr4<0||tr4>w2) Transform(sr4,tr4,-1);
                           else                              Transform(sr4,tr4,0);
                           dxr2  = fDxspline;
                           dyr2  = fDyspline;
                           zr2   = fZ;
                           v     = ColorCalculation(dx2,dy2,z2,dx3,dy3,z3,dx4,dy4,z4);
                           v     = v+ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr1,dyr1,zr1);
                           v     = v+ColorCalculation(dx3,dy3,z3,dxr2,dyr2,zr2,dxr1,dyr1,zr1);
                           v     = v+ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dx3,dy3,z3);
                           v3    = v/4;
                           if (sr7<0||sr7>w1||tr7<0||tr7>w2) Transform(sr7,tr7,-1);
                           else                              Transform(sr7,tr7,0);
                           dxr1  = fDxspline;
                           dyr1  = fDyspline;
                           zr1   = fZ;
                           if (sr6<0||sr6>w1||tr6<0||tr6>w2) Transform(sr6,tr6,-1);
                           else                              Transform(sr6,tr6,0);
                           dxr2  = fDxspline;
                           dyr2  = fDyspline;
                           zr2   = fZ;
                           v     = ColorCalculation(dx1,dy1,z1,dx3,dy3,z3,dx4,dy4,z4);
                           v     = v+ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr2,dyr2,zr2);
                           v     = v+ColorCalculation(dx4,dy4,z4,dxr2,dyr2,zr2,dxr1,dyr1,zr1);
                           v     = v+ColorCalculation(dx1,dy1,z1,dx4,dy4,z4,dxr1,dyr1,zr1);
                           v4    = v/4;
                        } else {
                           spriz = ShadowColorCalculation(xtaz,ytaz,ztaz,shad_noise);
                           v     = v+spriz;
                           v     = v/2;
                           if (sr1<0||sr1>w1||tr1<0||tr1>w2) Transform(sr1,tr1,-1);
                           else                              Transform(sr1,tr1,0);
                           dxr1  = fDxspline;
                           dyr1  = fDyspline;
                           zr1   = fZ;
                           if (sr8<0||sr8>w1||tr8<0||tr8>w2) Transform(sr8,tr8,-1);
                           else                              Transform(sr8,tr8,0);
                           dxr2  = fDxspline;
                           dyr2  = fDyspline;
                           zr2   = fZ;
                           da    = (dxr1+dx2+dx1)/3;
                           db    = (dyr1+dy2+dy1)/3;
                           dc    = (zr1+z2+z1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dxr1,dyr1,zr1,dx2,dy2,z2,dx1,dy1,z1)+spriz)/2;
                           da    = (dxr1+dxr2+dx1)/3;
                           db    = (dyr1+dyr2+dy1)/3;
                           dc    = (zr1+zr2+z1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dxr1,dyr1,zr1,dx1,dy1,z1,dxr2,dyr2,zr2)+spriz)/2;
                           da    = (dxr2+dx1+dx4)/3;
                           db    = (dyr2+dy1+dy4)/3;
                           dc    = (zr2+z1+z4)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dxr2,dyr2,zr2,dx1,dy1,z1,dx4,dy4,z4)+spriz)/2;
                           v1    = v/4;
                           if (sr3<0||sr3>w1||tr3<0||tr3>w2) Transform(sr3,tr3,-1);
                           else                              Transform(sr3,tr3,0);
                           dxr1  = fDxspline;
                           dyr1  = fDyspline;
                           zr1   = fZ;
                           if (sr2<0||sr2>w1||tr2<0||tr2>w2) Transform(sr2,tr2,-1);
                           else                              Transform(sr2,tr2,0);
                           dxr2  = fDxspline;
                           dyr2  = fDyspline;
                           zr2   = fZ;
                           da    = (dx1+dx2+dx3)/3;
                           db    = (dy1+dy2+dy3)/3;
                           dc    = (z1+z2+z3)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = (ColorCalculation(dx1,dy1,z1,dx2,dy2,z2,dx3,dy3,z3)+spriz)/2;
                           da    = (dx2+dxr1+dx3)/3;
                           db    = (dy2+dyr1+dy3)/3;
                           dc    = (z2+zr1+z3)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx2,dy2,z2,dxr1,dyr1,zr1,dx3,dy3,z3)+spriz)/2;
                           da    = (dx2+dxr2+dxr1)/3;
                           db    = (dy2+dyr2+dyr1)/3;
                           dc    = (z2+zr2+zr1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dxr1,dyr1,zr1)+spriz)/2;
                           da    = (dxr2+dx2+dx1)/3;
                           db    = (dyr2+dy2+dy1)/3;
                           dc    = (zr2+z2+z1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dxr2,dyr2,zr2,dx2,dy2,z2,dx1,dy1,z1)+spriz)/2;
                           v2    = v/4;
                           if (sr5<0||sr5>w1||tr5<0||tr5>w2) Transform(sr5,tr5,-1);
                           else                              Transform(sr5,tr5,0);
                           dxr1  = fDxspline;
                           dyr1  = fDyspline;
                           zr1   = fZ;
                           if (sr4<0||sr4>w1||tr4<0||tr4>w2) Transform(sr4,tr4,-1);
                           else                              Transform(sr4,tr4,0);
                           dxr2  = fDxspline;
                           dyr2  = fDyspline;
                           zr2   = fZ;
                           da    = (dx2+dx3+dx4)/3;
                           db    = (dy2+dy3+dy4)/3;
                           dc    = (z2+z3+z4)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = (ColorCalculation(dx2,dy2,z2,dx3,dy3,z3,dx4,dy4,z4)+spriz)/2;
                           da    = (dx4+dx3+dxr1)/3;
                           db    = (dy4+dy3+dyr1)/3;
                           dc    = (z4+z3+zr1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr1,dyr1,zr1)+spriz)/2;
                           da    = (dx3+dxr2+dxr1)/3;
                           db    = (dy3+dyr2+dyr1)/3;
                           dc    = (z3+zr2+zr1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx3,dy3,z3,dxr2,dyr2,zr2,dxr1,dyr1,zr1)+spriz)/2;
                           da    = (dx2+dxr2+dx3)/3;
                           db    = (dy2+dyr2+dy3)/3;
                           dc    = (z2+zr2+z3)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dx3,dy3,z3)+spriz)/2;
                           v3    = v/4;
                           if (sr7<0||sr7>w1||tr7<0||tr7>w2) Transform(sr7,tr7,-1);
                           else                              Transform(sr7,tr7,0);
                           dxr1  = fDxspline;
                           dyr1  = fDyspline;
                           zr1   = fZ;
                           if (sr6<0||sr6>w1||tr6<0||tr6>w2) Transform(sr6,tr6,-1);
                           else                              Transform(sr6,tr6,0);
                           dxr2  = fDxspline;
                           dyr2  = fDyspline;
                           zr2   = fZ;
                           da    = (dx1+dx3+dx4)/3;
                           db    = (dy1+dy3+dy4)/3;
                           dc    = (z1+z3+z4)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = (ColorCalculation(dx1,dy1,z1,dx3,dy3,z3,dx4,dy4,z4)+spriz)/2;
                           da    = (dx4+dx3+dxr2)/3;
                           db    = (dy4+dy3+dyr2)/3;
                           dc    = (z4+z3+zr2)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr2,dyr2,zr2)+spriz)/2;
                           da    = (dx4+dxr2+dxr1)/3;
                           db    = (dy4+dyr2+dyr1)/3;
                           dc    = (z4+zr2+zr1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx4,dy4,z4,dxr2,dyr2,zr2,dxr1,dyr1,zr1)+spriz)/2;
                           da    = (dx1+dx4+dxr1)/3;
                           db    = (dy1+dy4+dyr1)/3;
                           dc    = (z1+z4+zr1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx1,dy1,z1,dx4,dy4,z4,dxr1,dyr1,zr1)+spriz)/2;
                           v4    = v/4;
                        }
                     }
                     spriz = 0;
                     if (fShadow==kShadowsNotPainted) {
                        if (fShading==kNotShaded) {
                           v  = v*fLevels+0.5;
                           iv = fLevels-(Int_t)v;
                        } else {
                           v1  = v1*fLevels;
                           iv1 = fLevels-(Int_t)v1;
                           v2  = v2*fLevels;
                           iv2 = fLevels-(Int_t)v2;
                           v4 = v4*fLevels;
                           iv4 = fLevels-(Int_t)v4;
                        }
                     } else {
                        spriz = ShadowColorCalculation(xtaz,ytaz,ztaz,shad_noise);
                        if (fShading==kNotShaded) {
                           v  = v*fLevels/2.0;
                           iv = fLevels-(Int_t)(v+0.5);
                        } else {
                           v1  = v1*fLevels;
                           iv1 = fLevels-(Int_t)v1;
                           v2  = v2*fLevels;
                           iv2 = fLevels-(Int_t)v2;
                           v4  = v4*fLevels;
                           iv4 = fLevels-(Int_t)v4;
                        }
                     }
                     if (fShading==kNotShaded) {
                        ColorModel(iv,ui1,ui2,ui3);
                        line->SetLineColor(fNewColorIndex);
                        if (fEnvelope[x1]>=y1) {
                           line->PaintLine(gPad->PixeltoX(x1),gPad->PixeltoY(y1)+1,gPad->PixeltoX(x1+1),gPad->PixeltoY(y1)+1);
                           fEnvelope[x1] = y1;
                        }
                        if (fEnvelope[x2]>=y2) {
                           line->PaintLine(gPad->PixeltoX(x2),gPad->PixeltoY(y2)+1,gPad->PixeltoX(x2+1),gPad->PixeltoY(y2)+1);
                           fEnvelope[x2] = y2;
                        }
                        if (fEnvelope[x4]>=y4) {
                           line->PaintLine(gPad->PixeltoX(x4),gPad->PixeltoY(y4)+1,gPad->PixeltoX(x4+1),gPad->PixeltoY(y4)+1);
                           fEnvelope[x4] = y4;
                        }
                     } else {
                        if (fEnvelope[x1]>=y1) {
                           iv = iv1;
                           ColorModel(iv,ui1,ui2,ui3);
                           line->SetLineColor(fNewColorIndex);
                           line->PaintLine(gPad->PixeltoX(x1),gPad->PixeltoY(y1)+1,gPad->PixeltoX(x1+1),gPad->PixeltoY(y1)+1);
                           fEnvelope[x1] = y1;
                        }
                        if (fEnvelope[x2]>=y2) {
                           iv = iv2;
                           ColorModel(iv,ui1,ui2,ui3);
                           line->SetLineColor(fNewColorIndex);
                           line->PaintLine(gPad->PixeltoX(x2),gPad->PixeltoY(y2)+1,gPad->PixeltoX(x2+1),gPad->PixeltoY(y2)+1);
                           fEnvelope[x2]=y2;
                        }
                        if (fEnvelope[x4]>=y4) {
                           iv = iv4;
                           ColorModel(iv,ui1,ui2,ui3);
                           line->SetLineColor(fNewColorIndex);
                           line->PaintLine(gPad->PixeltoX(x4),gPad->PixeltoY(y4)+1,gPad->PixeltoX(x4+1),gPad->PixeltoY(y4)+1);
                           fEnvelope[x4] = y4;
                        }
                     }
                     xtaz = (dx3+dx2+dx4)/3;
                     ytaz = (dy3+dy2+dy4)/3;
                     ztaz = (z3+z2+z4)/3;
                     if (fShading==kNotShaded) v=ColorCalculation(dx2,dy2,z2,dx3,dy3,z3,dx4,dy4,z4);
                     spriz = 0;
                     if (fShadow==kShadowsNotPainted) {
                        if (fShading==kNotShaded) {
                           v  = v*fLevels;
                           iv = fLevels-(Int_t)v;
                        } else {
                           v3  = v3*fLevels;
                           iv3 = fLevels-(Int_t)v3;
                        }
                     } else {
                        spriz = ShadowColorCalculation(xtaz,ytaz,ztaz,shad_noise);
                        if (fShading==kNotShaded) {
                           v  = v*fLevels/2;
                           iv = fLevels-(Int_t)v;
                           iv = (Int_t)(iv-fLevels*spriz/2);
                        } else {
                           v3  = v3*fLevels;
                           iv3 = fLevels-(Int_t)v3;
                        }
                     }
                     if (fShading==kNotShaded) {
                        ColorModel(iv,ui1,ui2,ui3);
                        line->ResetAttLine("");
                        line->SetLineColor(fNewColorIndex);
                        if (fEnvelope[x3]>=y3) {
                           line->PaintLine(gPad->PixeltoX(x3),gPad->PixeltoY(y3)+1,gPad->PixeltoX(x3+1),gPad->PixeltoY(y3)+1);
                           fEnvelope[x3] = y3;
                        }
                     } else {
                        if (fEnvelope[x3]>=y3) {
                           iv = iv3;
                           ColorModel(iv,ui1,ui2,ui3);
                           line->ResetAttLine("");
                           line->SetLineColor(fNewColorIndex);
                           line->PaintLine(gPad->PixeltoX(x3),gPad->PixeltoY(y3)+1,gPad->PixeltoX(x3+1),gPad->PixeltoY(y3)+1);
                           fEnvelope[x3]=y3;
                        }
                     }
                  }
               }
               break;
            case kDisplayModeGrid:
               if (fBezier==kNoBezierInterpol) {
                  if (fModeGroup==kModeGroupSimple) {
                        Envelope(x1,y1,x2,y2);
                        if (fLine!=0) {
                           if (fLine==1) {
                              fXe = x2;
                              fYe = y2;
                           }
                           line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                        }
                  } else {
                     if ((q1!=q2||smer!=0)&&flag==1) {
                        s1 = q1+1;
                        t1 = q2;
                        s2 = q1;
                        t2 = q2;
                        s3 = q1;
                        t3 = q2+1;
                        s4 = q1+1;
                        t4 = q2+1;
                        if (fShading==kShaded) {
                           sr1 = s1;
                           tr1 = (Int_t)TMath::Max(t1-1,0);
                           sr2 = s2;
                           tr2 = (Int_t)TMath::Max(t2-1,0);
                           sr3 = (Int_t)TMath::Max(s2-1,0);
                           tr3 = t2;
                           sr4 = (Int_t)TMath::Max(s3-1,0);
                           tr4 = t3;
                           sr5 = s3;
                           tr5 = t3+1;
                           sr6 = s4;
                           tr6 = t4+1;
                           sr7 = s4+1;
                           tr7 = t4;
                           sr8 = s1+1;
                           tr8 = t1;
                        }
                        if (turni==1) {
                           i1 = s1;
                           i2 = s2;
                           i3 = s3;
                           i4 = s4;
                        } else {
                           i1 = (Int_t)TMath::Max(w1-s1,0);
                           i2 = (Int_t)TMath::Max(w1-s2,0);
                           i3 = (Int_t)TMath::Max(w1-s3,0);
                           i4 = (Int_t)TMath::Max(w1-s4,0);
                           if (fShading==kShaded) {
                              sr1 = (Int_t)TMath::Max(w1-sr1,0);
                              sr2 = (Int_t)TMath::Max(w1-sr2,0);
                              sr3 = (Int_t)TMath::Max(w1-sr3,0);
                              sr4 = (Int_t)TMath::Max(w1-sr4,0);
                              sr5 = (Int_t)TMath::Max(w1-sr5,0);
                              sr6 = (Int_t)TMath::Max(w1-sr6,0);
                              sr7 = (Int_t)TMath::Max(w1-sr7,0);
                              sr8 = (Int_t)TMath::Max(w1-sr8,0);
                           }
                        }
                        if (turnj==1) {
                           j1 = t1;
                           j2 = t2;
                           j3 = t3;
                           j4 = t4;
                        } else {
                           j1 = (Int_t)TMath::Max(w2-t1,0);
                           j2 = (Int_t)TMath::Max(w2-t2,0);
                           j3 = (Int_t)TMath::Max(w2-t3,0);
                           j4 = (Int_t)TMath::Max(w2-t4,0);
                           if (fShading==kShaded) {
                              tr1 = (Int_t)TMath::Max(w2-tr1,0);
                              tr2 = (Int_t)TMath::Max(w2-tr2,0);
                              tr3 = (Int_t)TMath::Max(w2-tr3,0);
                              tr4 = (Int_t)TMath::Max(w2-tr4,0);
                              tr5 = (Int_t)TMath::Max(w2-tr5,0);
                              tr6 = (Int_t)TMath::Max(w2-tr6,0);
                              tr7 = (Int_t)TMath::Max(w2-tr7,0);
                              tr8 = (Int_t)TMath::Max(w2-tr8,0);
                           }
                        }
                        Transform(i1,j1,0);
                        x1  = fXt;
                        y1  = fYt;
                        dx1 = fDxspline;
                        dy1 = fDyspline;
                        z1  = fZ;
                        Transform(i2,j2,0);
                        x2  = fXt;
                        y2  = fYt;
                        dx2 = fDxspline;
                        dy2 = fDyspline;
                        z2  = fZ;
                        Transform(i3,j3,0);
                        x3  = fXt;
                        y3  = fYt;
                        dx3 = fDxspline;
                        dy3 = fDyspline;
                        z3  = fZ;
                        Transform(i4,j4,0);
                        x4  = fXt;
                        y4  = fYt;
                        dx4 = fDxspline;
                        dy4 = fDyspline;
                        z4  = fZ;
                        Envelope(x1,y1,x2,y2);
                        Envelope(x2,y2,x3,y3);
                        xtaz = (dx1+dx2+dx4)/3;
                        ytaz = (dy1+dy2+dy4)/3;
                        ztaz = (z1+z2+z4)/3;
                        v    = ColorCalculation(dx1,dy1,z1,dx2,dy2,z2,dx4,dy4,z4);
                        if (fShading==kShaded) {
                           if (fShadow==kShadowsNotPainted) {
                              if (sr1<0||sr1>w1||tr1<0||tr1>w2) Transform(sr1,tr1,-1);
                              else                              Transform(sr1,tr1,0);
                              dxr1  = fDxspline;
                              dyr1  = fDyspline;
                              zr1   = fZ;
                              if (sr8<0||sr8>w1||tr8<0||tr8>w2) Transform(sr8,tr8,-1);
                              else                              Transform(sr8,tr8,0);
                              dxr2  = fDxspline;
                              dyr2  = fDyspline;
                              zr2   = fZ;
                              v     = v+ColorCalculation(dxr1,dyr1,zr1,dx2,dy2,z2,dx1,dy1,z1);
                              v     = v+ColorCalculation(dxr1,dyr1,zr1,dx1,dy1,z1,dxr2,dyr2,zr2);
                              v     = v+ColorCalculation(dxr2,dyr2,zr2,dx1,dy1,z1,dx4,dy4,z4);
                              v1    = v/4;
                              if (sr3<0||sr3>w1||tr3<0||tr3>w2) Transform(sr3,tr3,-1);
                              else                              Transform(sr3,tr3,0);
                              dxr1  = fDxspline;
                              dyr1  = fDyspline;
                              zr1   = fZ;
                              if (sr2<0||sr2>w1||tr2<0||tr2>w2) Transform(sr2,tr2,-1);
                              else                              Transform(sr2,tr2,0);
                              dxr2  = fDxspline;
                              dyr2  = fDyspline;
                              zr2   = fZ;
                              v     = ColorCalculation(dx1,dy1,z1,dx2,dy2,z2,dx3,dy3,z3);
                              v     = v+ColorCalculation(dx2,dy2,z2,dxr1,dyr1,zr1,dx3,dy3,z3);
                              v     = v+ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dxr1,dyr1,zr1);
                              v     = v+ColorCalculation(dxr2,dyr2,zr2,dx2,dy2,z2,dx1,dy1,z1);
                              v2    = v/4;
                              if (sr5<0||sr5>w1||tr5<0||tr5>w2) Transform(sr5,tr5,-1);
                              else                              Transform(sr5,tr5,0);
                              dxr1  = fDxspline;
                              dyr1  = fDyspline;
                              zr1   = fZ;
                              if (sr4<0||sr4>w1||tr4<0||tr4>w2) Transform(sr4,tr4,-1);
                              else                              Transform(sr4,tr4,0);
                              dxr2  = fDxspline;
                              dyr2  = fDyspline;
                              zr2   = fZ;
                              v     = ColorCalculation(dx2,dy2,z2,dx3,dy3,z3,dx4,dy4,z4);
                              v     = v+ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr1,dyr1,zr1);
                              v     = v+ColorCalculation(dx3,dy3,z3,dxr2,dyr2,zr2,dxr1,dyr1,zr1);
                              v     = v+ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dx3,dy3,z3);
                              v3    = v/4;
                              if (sr7<0||sr7>w1||tr7<0||tr7>w2) Transform(sr7,tr7,-1);
                              else                              Transform(sr7,tr7,0);
                              dxr1  = fDxspline;
                              dyr1  = fDyspline;
                              zr1   = fZ;
                              if (sr6<0||sr6>w1||tr6<0||tr6>w2) Transform(sr6,tr6,-1);
                              else                              Transform(sr6,tr6,0);
                              dxr2  = fDxspline;
                              dyr2  = fDyspline;
                              zr2   = fZ;
                              v     = ColorCalculation(dx1,dy1,z1,dx3,dy3,z3,dx4,dy4,z4);
                              v     = v+ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr2,dyr2,zr2);
                              v     = v+ColorCalculation(dx4,dy4,z4,dxr2,dyr2,zr2,dxr1,dyr1,zr1);
                              v     = v+ColorCalculation(dx1,dy1,z1,dx4,dy4,z4,dxr1,dyr1,zr1);
                              v4    = v/4;
                           } else {
                              spriz = ShadowColorCalculation(xtaz,ytaz,ztaz,shad_noise);
                              v     = v+spriz;
                              v     = v/2;
                              if (sr1<0||sr1>w1||tr1<0||tr1>w2) Transform(sr1,tr1,-1);
                              else                              Transform(sr1,tr1,0);
                              dxr1  = fDxspline;
                              dyr1  = fDyspline;
                              zr1   = fZ;
                              if (sr8<0||sr8>w1||tr8<0||tr8>w2) Transform(sr8,tr8,-1);
                              else                              Transform(sr8,tr8,0);
                              dxr2  = fDxspline;
                              dyr2  = fDyspline;
                              zr2   = fZ;
                              da    = (dxr1+dx2+dx1)/3;
                              db    = (dyr1+dy2+dy1)/3;
                              dc    = (zr1+z2+z1)/3;
                              spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                              v     = v+(ColorCalculation(dxr1,dyr1,zr1,dx2,dy2,z2,dx1,dy1,z1)+spriz)/2;
                              da    = (dxr1+dxr2+dx1)/3;
                              db    = (dyr1+dyr2+dy1)/3;
                              dc    = (zr1+zr2+z1)/3;
                              spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                              v     = v+(ColorCalculation(dxr1,dyr1,zr1,dx1,dy1,z1,dxr2,dyr2,zr2)+spriz)/2;
                              da    = (dxr2+dx1+dx4)/3;
                              db    = (dyr2+dy1+dy4)/3;
                              dc    = (zr2+z1+z4)/3;
                              spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                              v     = v+(ColorCalculation(dxr2,dyr2,zr2,dx1,dy1,z1,dx4,dy4,z4)+spriz)/2;
                              v1    = v/4;
                              if (sr3<0||sr3>w1||tr3<0||tr3>w2) Transform(sr3,tr3,-1);
                              else                              Transform(sr3,tr3,0);
                              dxr1  = fDxspline;
                              dyr1  = fDyspline;
                              zr1   = fZ;
                              if (sr2<0||sr2>w1||tr2<0||tr2>w2) Transform(sr2,tr2,-1);
                              else                              Transform(sr2,tr2,0);
                              dxr2  = fDxspline;
                              dyr2  = fDyspline;
                              zr2   = fZ;
                              da    = (dx1+dx2+dx3)/3;
                              db    = (dy1+dy2+dy3)/3;
                              dc    = (z1+z2+z3)/3;
                              spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                              v     = (ColorCalculation(dx1,dy1,z1,dx2,dy2,z2,dx3,dy3,z3)+spriz)/2;
                              da    = (dx2+dxr1+dx3)/3;
                              db    = (dy2+dyr1+dy3)/3;
                              dc    = (z2+zr1+z3)/3;
                              spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                              v     = v+(ColorCalculation(dx2,dy2,z2,dxr1,dyr1,zr1,dx3,dy3,z3)+spriz)/2;
                              da    = (dx2+dxr2+dxr1)/3;
                              db    = (dy2+dyr2+dyr1)/3;
                              dc    = (z2+zr2+zr1)/3;
                              spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                              v     = v+(ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dxr1,dyr1,zr1)+spriz)/2;
                              da    = (dxr2+dx2+dx1)/3;
                              db    = (dyr2+dy2+dy1)/3;
                              dc    = (zr2+z2+z1)/3;
                              spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                              v     = v+(ColorCalculation(dxr2,dyr2,zr2,dx2,dy2,z2,dx1,dy1,z1)+spriz)/2;
                              v2    = v/4;
                              if (sr5<0||sr5>w1||tr5<0||tr5>w2) Transform(sr5,tr5,-1);
                              else                              Transform(sr5,tr5,0);
                              dxr1  = fDxspline;
                              dyr1  = fDyspline;
                              zr1   = fZ;
                              if (sr4<0||sr4>w1||tr4<0||tr4>w2) Transform(sr4,tr4,-1);
                              else                              Transform(sr4,tr4,0);
                              dxr2  = fDxspline;
                              dyr2  = fDyspline;
                              zr2   = fZ;
                              da    = (dx2+dx3+dx4)/3;
                              db    = (dy2+dy3+dy4)/3;
                              dc    = (z2+z3+z4)/3;
                              spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                              v     = (ColorCalculation(dx2,dy2,z2,dx3,dy3,z3,dx4,dy4,z4)+spriz)/2;
                              da    = (dx4+dx3+dxr1)/3;
                              db    = (dy4+dy3+dyr1)/3;
                              dc    = (z4+z3+zr1)/3;
                              spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                              v     = v+(ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr1,dyr1,zr1)+spriz)/2;
                              da    = (dx3+dxr2+dxr1)/3;
                              db    = (dy3+dyr2+dyr1)/3;
                              dc    = (z3+zr2+zr1)/3;
                              spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                              v     = v+(ColorCalculation(dx3,dy3,z3,dxr2,dyr2,zr2,dxr1,dyr1,zr1)+spriz)/2;
                              da    = (dx2+dxr2+dx3)/3;
                              db    = (dy2+dyr2+dy3)/3;
                              dc    = (z2+zr2+z3)/3;
                              spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                              v     = v+(ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dx3,dy3,z3)+spriz)/2;
                              v3    = v/4;
                              if (sr7<0||sr7>w1||tr7<0||tr7>w2) Transform(sr7,tr7,-1);
                              else                              Transform(sr7,tr7,0);
                              dxr1  = fDxspline;
                              dyr1  = fDyspline;
                              zr1   = fZ;
                              if (sr6<0||sr6>w1||tr6<0||tr6>w2) Transform(sr6,tr6,-1);
                              else                              Transform(sr6,tr6,0);
                              dxr2  = fDxspline;
                              dyr2  = fDyspline;
                              zr2   = fZ;
                              da    = (dx1+dx3+dx4)/3;
                              db    = (dy1+dy3+dy4)/3;
                              dc    = (z1+z3+z4)/3;
                              spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                              v     = (ColorCalculation(dx1,dy1,z1,dx3,dy3,z3,dx4,dy4,z4)+spriz)/2;
                              da    = (dx4+dx3+dxr2)/3;
                              db    = (dy4+dy3+dyr2)/3;
                              dc    = (z4+z3+zr2)/3;
                              spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                              v     = v+(ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr2,dyr2,zr2)+spriz)/2;
                              da    = (dx4+dxr2+dxr1)/3;
                              db    = (dy4+dyr2+dyr1)/3;
                              dc    = (z4+zr2+zr1)/3;
                              spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                              v     = v+(ColorCalculation(dx4,dy4,z4,dxr2,dyr2,zr2,dxr1,dyr1,zr1)+spriz)/2;
                              da    = (dx1+dx4+dxr1)/3;
                              db    = (dy1+dy4+dyr1)/3;
                              dc    = (z1+z4+zr1)/3;
                              spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                              v     = v+(ColorCalculation(dx1,dy1,z1,dx4,dy4,z4,dxr1,dyr1,zr1)+spriz)/2;
                              v4    = v/4;
                           }
                        }
                        spriz = 0;
                        if (fShadow==kShadowsNotPainted) {
                           if (fShading==kNotShaded) {
                              v  = v*fLevels+0.5;
                              iv = fLevels-(Int_t)v;
                           } else {
                              v1  = v1*fLevels;
                              iv1 = fLevels-(Int_t)v1;
                              v2  = v2*fLevels;
                              iv2 = fLevels-(Int_t)v2;
                              v4  = v4*fLevels;
                              iv4 = fLevels-(Int_t)v4;
                           }
                        } else {
                           spriz = ShadowColorCalculation(xtaz,ytaz,ztaz,shad_noise);
                           if (fShading==kNotShaded) {
                              v  = v*fLevels/2.0;
                              iv = fLevels-(Int_t)(v+0.5);
                           } else {
                              v1  = v1*fLevels;
                              iv1 = fLevels-(Int_t)v1;
                              v2  = v2*fLevels;
                              iv2 = fLevels-(Int_t)v2;
                              v4  = v4*fLevels;
                              iv4 = fLevels-(Int_t)v4;
                           }
                        }
                        if (fShading==kNotShaded) {
                           ColorModel(iv,ui1,ui2,ui3);
                           line->SetLineColor(fNewColorIndex);
                        } else {
                           dx1 = x1;
                           dy1 = y1;
                           dx2 = x2;
                           dy2 = y2;
                           dx3 = x4;
                           dy3 = y4;
                           z1  = iv1;
                           z2  = iv2;
                           z3  = iv4;
                           da  = (dy2-dy1)*(z3-z1)-(dy3-dy1)*(z2-z1);
                           db  = (z2-z1)*(dx3-dx1)-(z3-z1)*(dx2-dx1);
                           dc  = (dx2-dx1)*(dy3-dy1)-(dx3-dx1)*(dy2-dy1);
                           dd  = -da*dx1-db*dy1-dc*z1;
                        }
                        sx1 = x1;
                        sy1 = y1;
                        sx2 = x2;
                        sy2 = y2;
                        if (sx2<sx1) {
                           sx4 = sx1;
                           sy4 = sy1;
                           sx1 = sx2;
                           sy1 = sy2;
                           sx2 = sx4;
                           sy2 = sy4;
                        }
                        sdx1 = 0;
                        pom1 = sy2-sy1;
                        pom2 = sx2-sx1;
                        if (pom2!=0) sdx1 = pom1/pom2;
                        pom1 = sy1;
                        pom2 = sx1;
                        sdy1 = pom1-sdx1*pom2;
                        for (sx4=sx1,sx5=sx1,sy5=sy1;sx4<=sx2;sx4++) {
                           pom1 = sx4;
                           sdy4 = sdx1*pom1+sdy1;
                           sy4  = (Int_t)(sdy4);
                           if (sy4<=fEnvelope[sx4]) {
                              fEnvelope[sx4] = sy4;
                              if (fShading==kNotShaded) {
                                 line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                              } else {
                                 dx1 = sx4;
                                 dy1 = sy4;
                                 if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                                 else       v = (iv1+iv2+iv4)/3;
                                 iv = (Int_t)v;
                                 ColorModel(iv,ui1,ui2,ui3);
                                 line->SetLineColor(fNewColorIndex);
                                 line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                              }
                              sy5 = sy4;
                           } else {
                              sy4 = fEnvelope[sx4];
                              if (fShading==kNotShaded&&sy5<=fEnvelope[sx5]) {
                                 line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                              } else if (sy5<=fEnvelope[sx5]) {
                                 dx1 = sx4;
                                 dy1 = sy4;
                              if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                              else       v = (iv1+iv2+iv4)/3;
                              iv = (Int_t)v;
                              ColorModel(iv,ui1,ui2,ui3);
                              line->SetLineColor(fNewColorIndex);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                              }
                              sy5 = fEnvelope[sx4];
                           }
                           sx5 = sx4;
                        }
                        sx1 = x1;
                        sy1 = y1;
                        sx3 = x4;
                        sy3 = y4;
                        if (sx3<sx1) {
                           sx4 = sx1;
                           sy4 = sy1;
                           sx1 = sx3;
                           sy1 = sy3;
                           sx3 = sx4;
                           sy3 = sy4;
                        }
                        pom1  = sy3-sy1;
                        pom2  = sx3-sx1;
                        if (pom2!=0) sdx2 = pom1/pom2;
                        pom1  = sy1;
                        pom2  = sx1;
                        sdy2  = pom1-sdx2*pom2;
                        sx1p  = sx1;
                        sy1p  = sy1;
                        sx3p  = sx3;
                        sdx2p = sdx2;
                        sdy2p = sdy2;
                        dap   = da;
                        dbp   = db;
                        dcp   = dc;
                        ddp   = dd;
                        uip   = fNewColorIndex;
                        xtaz  = (dx3+dx2+dx4)/3;
                        ytaz  = (dy3+dy2+dy4)/3;
                        ztaz  = (z3+z2+z4)/3;
                        if (fShading==kNotShaded) v = ColorCalculation(dx2,dy2,z2,dx3,dy3,z3,dx4,dy4,z4);
                        spriz = 0;
                        if (fShadow==kShadowsNotPainted) {
                           if (fShading==kNotShaded) {
                              v  = v*fLevels;
                              iv = fLevels-(Int_t)v;
                           } else {
                              v3  = v3*fLevels;
                              iv3 = fLevels-(Int_t)v3;
                           }
                        } else {
                           spriz = ShadowColorCalculation(xtaz,ytaz,ztaz,shad_noise);
                           if (fShading==kNotShaded) {
                              v  = v*fLevels/2;
                              iv = fLevels-(Int_t)v;
                              iv = (Int_t)(iv-fLevels*spriz/2);
                           } else {
                              v3  = v3*fLevels;
                              iv3 = fLevels-(Int_t)v3;
                           }
                        }
                        if (fShading==kNotShaded) {
                           ColorModel(iv,ui1,ui2,ui3);
                           line->SetLineColor(fNewColorIndex);
                        } else {
                           dx1 = x2;
                           dy1 = y2;
                           dx2 = x3;
                           dy2 = y3;
                           dx3 = x4;
                           dy3 = y4;
                           z1  = iv2;
                           z2  = iv3;
                           z3  = iv4;
                           da  = (dy2-dy1)*(z3-z1)-(dy3-dy1)*(z2-z1);
                           db  = (z2-z1)*(dx3-dx1)-(z3-z1)*(dx2-dx1);
                           dc  = (dx2-dx1)*(dy3-dy1)-(dx3-dx1)*(dy2-dy1);
                           dd  = -da*dx1-db*dy1-dc*z1;
                        }
                        sx1 = x2;
                        sy1 = y2;
                        sx2 = x3;
                        sy2 = y3;
                        if (sx2<sx1) {
                           sx4 = sx1;
                           sy4 = sy1;
                           sx1 = sx2;
                           sy1 = sy2;
                           sx2 = sx4;
                           sy2 = sy4;
                        }
                        pom1 = sy2-sy1;
                        pom2 = sx2-sx1;
                        sdx1 = 0;
                        if (pom2!=0) sdx1 = pom1/pom2;
                        pom1 = sy1;
                        pom2 = sx1;
                        sdy1 = pom1-sdx1*pom2;
                        for (sx4=sx1,sx5=sx1,sy5=sy1;sx4<=sx2;sx4++) {
                           pom1 = sx4;
                           sdy4 = sdx1*pom1+sdy1;
                           sy4  = (Int_t)sdy4;
                           if (sy4<=fEnvelope[sx4]) {
                              fEnvelope[sx4] = sy4;
                              if (fShading==kNotShaded) {
                                 line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                              } else {
                                 dx1 = sx4;
                                 dy1 = sy4;
                                 if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                                 else       v = (iv1+iv2+iv4)/3;
                                 iv = (Int_t)v;
                                 ColorModel(iv,ui1,ui2,ui3);
                                 line->SetLineColor(fNewColorIndex);
                                 line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                              }
                              sy5 = sy4;
                           } else {
                              sy4 = fEnvelope[sx4];
                              if (fShading==kNotShaded&&sy5<=fEnvelope[sx5]) {
                                 line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                              } else if (sy5<=fEnvelope[sx5]) {
                                 dx1 = sx4;
                                 dy1 = sy4;
                                 if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                                 else       v = (iv1+iv2+iv4)/3;
                                 iv = (Int_t)v;
                                 ColorModel(iv,ui1,ui2,ui3);
                                 line->SetLineColor(fNewColorIndex);
                                 line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                              }
                              sy5 = fEnvelope[sx4];
                           }
                           sx5 = sx4;
                        }
                        for (sx4=sx1p,sx5=sx1p,sy5=sy1p;sx4<=sx3p;sx4++) {
                           pom1 = sx4;
                           sdy4 = sdx2p*pom1+sdy2p;
                           sy4  = (Int_t)sdy4;
                           if (sy4<=fEnvelope[sx4]) {
                              fEnvelope[sx4]=sy4;
                              if (fShading==kNotShaded) {
                                 line->SetLineColor(uip);
                                 line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                              } else {
                                 dx1 = sx4;
                                 dy1 = sy4;
                                 if (dcp!=0) v = (-ddp-dap*dx1-dbp*dy1)/dcp;
                                 else        v = (iv1+iv2+iv4)/3;
                                 iv = (Int_t)v;
                                 ColorModel(iv,ui1,ui2,ui3);
                                 line->SetLineColor(fNewColorIndex);
                                 line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                              }
                              sy5 = sy4;
                           } else {
                              sy4 = fEnvelope[sx4];
                              if (fShading==kNotShaded&&sy5<=fEnvelope[sx5]) {
                                 line->SetLineColor(uip);
                                 line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                              } else if (sy5<=fEnvelope[sx5]) {
                                 dx1 = sx4;
                                 dy1 = sy4;
                                 if (dcp!=0) v = (-ddp-dap*dx1-dbp*dy1)/dcp;
                                 else        v = (iv1+iv2+iv4)/3;
                                 iv = (Int_t)v;
                                 ColorModel(iv,ui1,ui2,ui3);
                                 line->SetLineColor(fNewColorIndex);
                                 line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                              }
                              sy5 = fEnvelope[sx4];
                           }
                           sx5 = sx4;
                        }
                        sx2 = x3;
                        sy2 = y3;
                        sx3 = x4;
                        sy3 = y4;
                        if (sx3<sx2) {
                           sx4 = sx2;
                           sy4 = sy2;
                           sx2 = sx3;
                           sy2 = sy3;
                           sx3 = sx4;
                           sy3 = sy4;
                        }
                        sdx2 = 0;
                        pom1 = sy3-sy2;
                        pom2 = sx3-sx2;
                        if (pom2!=0) sdx2 = pom1/pom2;
                        pom1 = sy2;
                        pom2 = sx2;
                        sdy2 = pom1-sdx2*pom2;
                        for (sx4=sx2,sx5=sx2,sy5=sy2;sx4<=sx3;sx4++) {
                           pom1 = sx4;
                           sdy4 = sdx2*pom1+sdy2;
                           sy4  = (Int_t)sdy4;
                           if (sy4<=fEnvelope[sx4]) {
                              fEnvelope[sx4] = sy4;
                              if (fShading==kNotShaded) {
                                 line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                              } else {
                                 dx1 = sx4;
                                 dy1 = sy4;
                                 if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                                 else       v =(iv1+iv2+iv4)/3;
                                 iv = (Int_t)v;
                                 ColorModel(iv,ui1,ui2,ui3);
                                 line->SetLineColor(fNewColorIndex);
                                 line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                              }
                              sy5 = sy4;
                           } else {
                              sy4 = fEnvelope[sx4];
                              if (fShading==kNotShaded&&sy5<=fEnvelope[sx5]) {
                                 line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                              } else if (sy5<=fEnvelope[sx5]) {
                                 dx1 = sx4;
                                 dy1 = sy4;
                                 if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                                 else       v =(iv1+iv2+iv4)/3;
                                 iv = (Int_t)v;
                                 ColorModel(iv,ui1,ui2,ui3);
                                 line->SetLineColor(fNewColorIndex);
                                 line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                              }
                              sy5 = fEnvelope[sx4];
                           }
                           sx5 = sx4;
                        }
                     }
                  }
               } else {
                  if (((flag==0)&&(smer==0))||((flag!=0)&&(smer!=0))) {
                     s1 = q1;
                     t1 = (Int_t)TMath::Max(q2-1,0);
                     s2 = q1;
                     t2 = (Int_t)TMath::Min(q2+2,w2);
                  } else if (((flag!=0)&&(smer==0))||((flag==0)&&(smer!=0))) {
                     s1 = (Int_t)TMath::Max(q1-1,0);
                     t1 = q2;
                     s2 = (Int_t)TMath::Min(q1+2,w1);
                     t2 = q2;
                  }
                  if (turni==1) {
                     i1 = s1;
                     i2 = s2;
                  } else {
                     i1 = w1-s1;
                     i2 = w1-s2;
                  }
                  if (turnj==1) {
                     j1 = t1;
                     j2 = t2;
                  } else {
                     j1 = w2-t1;
                     j2 = w2-t2;
                  }
                  Transform(i1,j1,0);
                  x3 = fXt;
                  y3 = fYt;
                  Transform(i2,j2,0);
                  x4    = fXt;
                  y4    = fYt;
                  bezx1 = x1+(x2-x1)/3;
                  bezx2 = x1+2*(x2-x1)/3;
                  bezy1 = y1+(y2-y3)/6;
                  bezy2 = y2-(y4-y1)/6;
                  if (x1<=x2) {
                     if (bezx1<=x1) {
                        bezx1 = x1;
                        bezy1 = y1;
                     }
                     if (bezx1>=x2) {
                        bezx1 = x2;
                        bezy1 = y2;
                     }
                     if (bezx2<=x1) {
                        bezx2 = x1;
                        bezy2 = y1;
                     }
                     if (bezx2>=x2) {
                        bezx2 = x2;
                        bezy2 = y2;
                     }
                     fBzX[0] = x1;
                     fBzY[0] = y1;
                     fBzX[1] = (Int_t)bezx1;
                     fBzY[1] = (Int_t)bezy1;
                     fBzX[2] = (Int_t)bezx2;
                     fBzY[2] = (Int_t)bezy2;
                     fBzX[3] = x2;
                     fBzY[3] = y2;
                     for (bezf=0;bezf<1.01;bezf+=0.1) {
                        BezierSmoothing(bezf);
                        if (bezf==0) {
                           ibezx1 = (Int_t)(fGbezx+0.5);
                           ibezy1 = (Int_t)(fGbezy+0.5);
                        } else {
                           ibezx2 = ibezx1;
                           ibezy2 = ibezy1;
                           ibezx1 = (Int_t)(fGbezx+0.5);
                           ibezy1 = (Int_t)(fGbezy+0.5);
                           Envelope(ibezx2,ibezy2,ibezx1,ibezy1);
                           if (fLine!=0) {
                              if (fLine==1) {
                                 fXe = ibezx1;
                                 fYe = ibezy1;
                              }
                              line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                           }
                        }
                     }
                  } else if (x1>x2) {
                     if (bezx1>=x1) {
                        bezx1 = x1;
                        bezy1 = y1;
                     }
                     if (bezx1<=x2) {
                        bezx1 = x2;
                        bezy1 = y2;
                     }
                     if (bezx2>=x1) {
                        bezx2 = x1;
                        bezy2 = y1;
                     }
                     if (bezx2<=x2) {
                        bezx2 = x2;
                        bezy2 = y2;
                     }
                     fBzX[0] = x1;
                     fBzY[0] = y1;
                     fBzX[1] = (Int_t)bezx1;
                     fBzY[1] = (Int_t)bezy1;
                     fBzX[2] = (Int_t)bezx2;
                     fBzY[2] = (Int_t)bezy2;
                     fBzX[3] = x2;
                     fBzY[3] = y2;
                     for (bezf=0;bezf<1.01;bezf+=0.1) {
                        BezierSmoothing(bezf);
                        if (bezf==0) {
                           ibezx1 = (Int_t)(fGbezx+0.5);
                           ibezy1 = (Int_t)(fGbezy+0.5);
                        } else {
                           ibezx2 = ibezx1;
                           ibezy2 = ibezy1;
                           ibezx1 = (Int_t)(fGbezx+0.5);
                           ibezy1 = (Int_t)(fGbezy+0.5);
                           Envelope(ibezx1,ibezy1,ibezx2,ibezy2);
                           if (fLine!=0) {
                              if (fLine==1) {
                                 fXe = ibezx2;
                                 fYe = ibezy2;
                              }
                              line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                           }
                        }
                     }
                  }
               }
               break;
            case kDisplayModeContours:
               if ((q1!=q2||smer!=0)&&flag==1) {
                  s1 = q1+1;
                  t1 = q2;
                  s2 = q1;
                  t2 = q2;
                  s3 = q1;
                  t3 = q2+1;
                  s4 = q1+1;
                  t4 = q2+1;
                  if (turni==1) {
                     i1 = (Int_t)TMath::Min(w1,s1);
                     i2 = (Int_t)TMath::Min(w1,s2);
                     i3 = (Int_t)TMath::Min(w1,s3);
                     i4 = (Int_t)TMath::Min(w1,s4);
                  } else {
                     i1 = (Int_t)TMath::Max(w1-s1,0);
                     i2 = (Int_t)TMath::Max(w1-s2,0);
                     i3 = (Int_t)TMath::Max(w1-s3,0);
                     i4 = (Int_t)TMath::Max(w1-s4,0);
                  }
                  if (turnj==1) {
                     j1 = (Int_t)TMath::Min(w2,t1);
                     j2 = (Int_t)TMath::Min(w2,t2);
                     j3 = (Int_t)TMath::Min(w2,t3);
                     j4 = (Int_t)TMath::Min(w2,t4);
                  } else {
                     j1 = (Int_t)TMath::Max(w2-t1,0);
                     j2 = (Int_t)TMath::Max(w2-t2,0);
                     j3 = (Int_t)TMath::Max(w2-t3,0);
                     j4 = (Int_t)TMath::Max(w2-t4,0);
                  }
                  Transform(i1,j1,0);
                  dx1 = fDxspline;
                  dy1 = fDyspline;
                  z1  = fZ;
                  z1l = fZeq;
                  Transform(i2,j2,0);
                  dx2 = fDxspline;
                  dy2 = fDyspline;
                  z2  = fZ;
                  z2l = fZeq;
                  Transform(i3,j3,0);
                  dx3 = fDxspline;
                  dy3 = fDyspline;
                  z3  = fZ;
                  z3l = fZeq;
                  Transform(i4,j4,0);
                  dx4 = fDxspline;
                  dy4 = fDyspline;
                  z4  = fZ;
                  z4l = fZeq;
                  zh  = (Double_t)TMath::Max(z1,z2);
                  zh  = (Double_t)TMath::Max(zh,z3);
                  zh  = (Double_t)TMath::Max(zh,z4);
                  zl  = (Double_t)TMath::Min(z1l,z2l);
                  zl  = (Double_t)TMath::Min(zl,z3l);
                  zl  = (Double_t)TMath::Min(zl,z4l);
                  i1 = (Int_t)(zl/dcount_reg+1);
                  if (z1!=z2||z2!=z3||z3!=z4) {
                     do {
                        fZ = i1*dcount_reg;
                        switch (fZscale) {
                           case kZScaleLog:
                              if (fZ>=1.0) fZ = log(fZ);
                              else         fZ = 0;
                              break;
                           case kZScaleSqrt:
                              if (fZ>0) fZ = sqrt(fZ);
                              else      fZ = 0;
                              break;
                        }
                        if (fModeGroup!=kModeGroupSimple) {
                           v  = ColorCalculation(dx1,dy1,fZ,dx2,dy2,fZ,dx4,dy4,fZ);
                           v  = v*fLevels+0.5;
                           iv = fLevels-(Int_t)v;
                           ColorModel(iv,ui1,ui2,ui3);
                           line->SetLineColor(fNewColorIndex);
                        }
                        if (fZ>zh) goto eqend;
                        i1   += 1;
                        ekv   = 0;
                        stvor = 0;
                        if ((z2<=fZ&&fZ<z1)||(z2<fZ&&fZ<=z1)) {
                           xb = (fZ-z2)*(dx1-dx2)/(z1-z2)+dx2;
                           goto ekvi1;
                        }
                        if ((z1<=fZ&&fZ<z2)||(z1<fZ&&fZ<=z2)) {
                           xb = (fZ-z1)*(dx2-dx1)/(z2-z1)+dx1;
                           goto ekvi1;
                        }
                        if (z2==fZ&&fZ==z1) {
                           xb     = dx2;
ekvi1:
                           yb     = dy2;
                           ekv    = 1;
                           x5     = xb;
                           y5     = yb;
                           stvor += 1;
                        }
                        if ((z1<=fZ&&fZ<z4)||(z1<fZ&&fZ<=z4)) {
                           ya = (fZ-z1)*(dy4-dy1)/(z4-z1)+dy1;
                           goto ekvi2;
                        }
                        if ((z4<=fZ&&fZ<z1)||(z4<fZ&&fZ<=z1)) {
                           ya = (fZ-z4)*(dy1-dy4)/(z1-z4)+dy4;
                           goto ekvi2;
                        }
                        if (z4==fZ&&fZ==z1) {
                           ya = dy1;
ekvi2:
                           xa = dx1;
                           if (ekv==1) {
                              Slice(xa,ya,xb,yb,line);
                              stvor += 1;
                           }
                           xb  = xa;
                           yb  = ya;
                           ekv = 1;
                        }
                        if ((z3<=fZ&&fZ<z4)||(z3<fZ&&fZ<=z4)) {
                           xa = (fZ-z3)*(dx4-dx3)/(z4-z3)+dx3;
                           goto ekvi3;
                        }
                        if ((z4<=fZ&&fZ<z3)||(z4<fZ&&fZ<=z3)) {
                           xa = (fZ-z4)*(dx3-dx4)/(z3-z4)+dx4;
                           goto ekvi3;
                        }
                        if (z4==fZ&&fZ==z3) {
                           xa = dx4;
ekvi3:
                           ya = dy4;
                           if (ekv==1) {
                              Slice(xa,ya,xb,yb,line);
                              stvor += 1;
                           }
                           xb  = xa;
                           yb  = ya;
                           ekv = 1;
                        }
                        if ((z2<=fZ&&fZ<z3)||(z2<fZ&&fZ<=z3)) {
                           ya = (fZ-z2)*(dy3-dy2)/(z3-z2)+dy2;
                           goto ekvi4;
                        }
                        if ((z3<=fZ&&fZ<z2)||(z3<fZ&&fZ<=z2)) {
                           ya = (fZ-z3)*(dy2-dy3)/(z2-z3)+dy3;
                           goto ekvi4;
                        }
                        if (z3==fZ&&fZ==z2) {
                           ya = dy3;
ekvi4:
                           xa = dx3;
                           if (ekv==1) {
                              Slice(xa,ya,xb,yb,line);
                              stvor += 1;
                           }
                           if (stvor==4) Slice(xa,ya,x5,y5,line);
                        }
                     } while (fZ<=zh);
eqend:
                     CopyEnvelope(dx1,dx3,dy1,dy3);
                  }
               }
               break;
            case kDisplayModeBars:
            case kDisplayModeBarsX:
            case kDisplayModeBarsY:
               if ((q1!=q2||smer!=0)&&flag==1) {
                  s1 = q1+1;
                  t1 = q2;
                  s2 = q1;
                  t2 = q2;
                  s3 = q1;
                  t3 = q2+1;
                  s4 = q1+1;
                  t4 = q2+1;
               }
               if (turni==1) {
                  if (fDisplayMode==kDisplayModeBarsX) {
                     if (s1<=w1&&s2<=w1&&s3<=w1&&s4<=w1) {
                        i1 = s1;
                        i2 = s2;
                        i3 = s3;
                        i4 = s4;
                     }
                  } else {
                     i1 = (Int_t)TMath::Min(w1,s1);
                     i2 = (Int_t)TMath::Min(w1,s2);
                     i3 = (Int_t)TMath::Min(w1,s3);
                     i4 = (Int_t)TMath::Min(w1,s4);
                  }
               } else {
                  if (fDisplayMode==kDisplayModeBarsX) {
                     if (s1<=w1&&s2<=w1&&s3<=w1&&s4<=w1) {
                        i1 = w1-s1;
                        i2 = w1-s2;
                        i3 = w1-s3;
                        i4 = w1-s4;
                     }
                  } else {
                     i1 = (Int_t)TMath::Max(w1-s1,0);
                     i2 = (Int_t)TMath::Max(w1-s2,0);
                     i3 = (Int_t)TMath::Max(w1-s3,0);
                     i4 = (Int_t)TMath::Max(w1-s4,0);
                  }
               }
               if (turnj==1) {
                  if (fDisplayMode==kDisplayModeBarsY) {
                     if (t1<=w2&&t2<=w2&&t3<=w2&&t4<=w2) {
                        j1 = t1;
                        j2 = t2;
                        j3 = t3;
                        j4 = t4;
                     }
                  } else {
                     j1 = (Int_t)TMath::Min(w2,t1);
                     j2 = (Int_t)TMath::Min(w2,t2);
                     j3 = (Int_t)TMath::Min(w2,t3);
                     j4 = (Int_t)TMath::Min(w2,t4);
                  }
               } else {
                  if (fDisplayMode==kDisplayModeBarsY) {
                     if (t1<=w2&&t2<=w2&&t3<=w2&&t4<=w2) {
                        j1 = w2-t1;
                        j2 = w2-t2;
                        j3 = w2-t3;
                        j4 = w2-t4;
                     }
                  } else {
                     j1 = (Int_t)TMath::Max(w2-t1,0);
                     j2 = (Int_t)TMath::Max(w2-t2,0);
                     j3 = (Int_t)TMath::Max(w2-t3,0);
                     j4 = (Int_t)TMath::Max(w2-t4,0);
                  }
               }
               Transform(i1,j1,0);
               x1  = fXt;
               dx1 = fDxspline;
               dy1 = fDyspline;
               z1  = fZ;
               Transform(i2,j2,0);
               x2  = fXt;
               dx2 = fDxspline;
               dy2 = fDyspline;
               z2  = fZ;
               Transform(i3,j3,0);
               x3  = fXt;
               dx3 = fDxspline;
               dy3 = fDyspline;
               z3  = fZ;
               Transform(i4,j4,0);
               x4  = fXt;
               y4  = fYt;
               dx4 = fDxspline;
               dy4 = fDyspline;
               z4  = fZ;
               Transform(i1,j1,-1);
               ix5 = fXt;
               iy5 = fYt;
               Transform(i2,j2,-1);
               x6  = fXt;
               y6  = fYt;
               Transform(i3,j3,-1);
               x7  = fXt;
               y7  = fYt;
               Transform(i4,j4,-1);
               y8  = fYt;
               y1  = iy5+(y4-y8);
               y2  = y6+(y4-y8);
               y3  = y7+(y4-y8);
               if ((fDisplayMode==kDisplayModeBars)&&(q1!=q2||smer!=0)&&(flag==1)) {
                  EnvelopeBars(ix5,iy5,x6,y6);
                  if (fLine!=0) {
                     if (fLine==1) {
                        fXe = x6;
                        fYe = y6;
                     }
                     line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                  }
                  EnvelopeBars(x6,y6,x7,y7);
                  if (fLine!=0) {
                     if (fLine==1) {
                        fXe = x7;
                        fYe = y7;
                     }
                     line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                  }
                  EnvelopeBars(ix5,iy5,x1,y1);
                  if (fLine!=0) {
                     if (fLine==1) {
                        fXe = x1;
                        fYe = y1;
                     }
                     line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                  }
                  EnvelopeBars(x6,y6,x2,y2);
                  if (fLine!=0) {
                     if (fLine==1) {
                        fXe = x2;
                        fYe = y2;
                     }
                     line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                  }
                  EnvelopeBars(x7,y7,x3,y3);
                  if (fLine!=0) {
                     if (fLine==1) {
                        fXe = x3;
                        fYe = y3;
                     }
                     line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                  }
                  if (fModeGroup!=kModeGroupSimple) {
                     v   = ColorCalculation(dx1,dy1,z4,dx2,dy2,z4,dx4,dy4,z4);
                     v   = v*fLevels+0.5;
                     iv  = fLevels-(Int_t)v;
                     uip = fNewColorIndex;
                     ColorModel(iv,ui1,ui2,ui3);
                     line->SetLineColor(fNewColorIndex);
                     sx1 = x1;
                     sy1 = y1;
                     sx2 = x2;
                     sy2 = y2;
                     sx3 = x4;
                     sy3 = y4;
                     if (sx2<sx1) {
                        sx4 = sx1;
                        sy4 = sy1;
                        sx1 = sx2;
                        sy1 = sy2;
                        sx2 = sx4;
                        sy2 = sy4;
                     }
                     if (sx3<sx1) {
                        sx4 = sx1;
                        sy4 = sy1;
                        sx1 = sx3;
                        sy1 = sy3;
                        sx3 = sx4;
                        sy3 = sy4;
                     }
                     if (sy2<sy3) {
                        sx4 = sx2;
                        sy4 = sy2;
                        sx2 = sx3;
                        sy2 = sy3;
                        sx3 = sx4;
                        sy3 = sy4;
                     }
                     sdx1 = 0;
                     sdx2 = 0;
                     sdx3 = 0;
                     pom1 = sy2-sy1;
                     pom2 = sx2-sx1;
                     if (pom2!=0) sdx1 = pom1/pom2;
                     pom1 = sy1;
                     pom2 = sx1;
                     sdy1 = pom1-sdx1*pom2;
                     pom1 = sy3-sy1;
                     pom2 = sx3-sx1;
                     if (pom2!=0) sdx2 = pom1/pom2;
                     pom1 = sy1;
                     pom2 = sx1;
                     sdy2 = pom1-sdx2*pom2;
                     pom1 = sy3-sy2;
                     pom2 = sx3-sx2;
                     if (pom2!=0) sdx3 = pom1/pom2;
                     pom1 = sy2;
                     pom2 = sx2;
                     sdy3 = pom1-sdx3*pom2;
                     if (sx2<sx3) {
                        if (sx1!=sx2) {
                           for (sx4=sx1;sx4<=sx2;sx4++) {
                              pom1 = sx4;
                              sdy4 = sdx1*pom1+sdy1;
                              sy4  = (Int_t)sdy4;
                              if (sx3!=sx1) {
                                 sdy4 = sdx2*pom1+sdy2;
                                 sy5  = (Int_t)sdy4;
                                 y5   = fEnvelope[sx4];
                                 if (sy4<sy5) {
                                    pom1 = sy4;
                                    sy4 = sy5;
                                    sy5 = (Int_t)pom1;
                                 }
                                 if ((sy4<=y5)||(sy5<y5)) {
                                    sy4 = (Int_t)TMath::Min(sy4,(Int_t)y5);
                                    sy5 = (Int_t)TMath::Min(sy5,(Int_t)y5);
                                    line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx4),gPad->PixeltoY(sy5)+1);
                                 }
                              }
                           }
                        }
                        if (sx2!=sx3) {
                           for (sx4=sx2;sx4<=sx3;sx4++) {
                              pom1 = sx4;
                              sdy4 = sdx3*pom1+sdy3;
                              sy4  = (Int_t)sdy4;
                              if (sx3!=sx1) {
                                 sdy4 = sdx2*pom1+sdy2;
                                 sy5  = (Int_t)sdy4;
                                 y5   = fEnvelope[sx4];
                                 if (sy4<sy5) {
                                    pom1 = sy4;
                                    sy4  = sy5;
                                    sy5  = (Int_t)pom1;
                                 }
                                 if ((sy4<=y5)||(sy5<y5)) {
                                    sy4 = (Int_t)TMath::Min(sy4,(Int_t)y5);
                                    sy5 = (Int_t)TMath::Min(sy5,(Int_t)y5);
                                    line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx4),gPad->PixeltoY(sy5)+1);
                                 }
                              }
                           }
                        }
                     } else {
                        if (sx3!=sx1) {
                           for (sx4=sx1;sx4<=sx3;sx4++) {
                              pom1 = sx4;
                              sdy4 = sdx2*pom1+sdy2;
                              sy4  = (Int_t)sdy4;
                              if (sx2!=sx1) {
                                 sdy4 = sdx1*pom1+sdy1;
                                 sy5  = (Int_t)sdy4;
                                 y5   = fEnvelope[sx4];
                                 if (sy4<sy5) {
                                    pom1 = sy4;
                                    sy4  = sy5;
                                    sy5  = (Int_t)pom1;
                                 }
                                 if ((sy4<=y5)||(sy5<y5)) {
                                    sy4 = (Int_t)TMath::Min(sy4,(Int_t)y5);
                                    sy5 = (Int_t)TMath::Min(sy5,(Int_t)y5);
                                    line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx4),gPad->PixeltoY(sy5)+1);
                                 }
                              }
                           }
                        }
                        if (sx2!=sx3) {
                           for (sx4=sx3;sx4<=sx2;sx4++) {
                              pom1 = sx4;
                              sdy4 = sdx3*pom1+sdy3;
                              sy4  = (Int_t)sdy4;
                              if (sx2!=sx1) {
                                 sdy4 = sdx1*pom1+sdy1;
                                 sy5  = (Int_t)sdy4;
                                 y5   = fEnvelope[sx4];
                                 if (sy4<sy5) {
                                    pom1 = sy4;
                                    sy4  = sy5;
                                    sy5  = (Int_t)pom1;
                                 }
                                 if ((sy4<=y5)||(sy5<y5)) {
                                    sy4 = (Int_t)TMath::Min(sy4,(Int_t)y5);
                                    sy5 = (Int_t)TMath::Min(sy5,(Int_t)y5);
                                    line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx4),gPad->PixeltoY(sy5)+1);
                                 }
                              }
                           }
                        }
                     }
                     sx1 = x2;
                     sy1 = y2;
                     sx2 = x3;
                     sy2 = y3;
                     sx3 = x4;
                     sy3 = y4;
                     if (sx2<sx1) {
                        sx4 = sx1;
                        sy4 = sy1;
                        sx1 = sx2;
                        sy1 = sy2;
                        sx2 = sx4;
                        sy2 = sy4;
                     }
                     if (sx3<sx1) {
                        sx4 = sx1;
                        sy4 = sy1;
                        sx1 = sx3;
                        sy1 = sy3;
                        sx3 = sx4;
                        sy3 = sy4;
                     }
                     if (sy2<sy3) {
                        sx4 = sx2;
                        sy4 = sy2;
                        sx2 = sx3;
                        sy2 = sy3;
                        sx3 = sx4;
                        sy3 = sy4;
                     }
                     sdx1 = 0;
                     sdx2 = 0;
                     sdx3 = 0;
                     pom1 = sy2-sy1;
                     pom2 = sx2-sx1;
                     if (pom2!=0) sdx1 = pom1/pom2;
                     pom1 = sy1;
                     pom2 = sx1;
                     sdy1 = pom1-sdx1*pom2;
                     pom1 = sy3-sy1;
                     pom2 = sx3-sx1;
                     if (pom2!=0) sdx2 = pom1/pom2;
                     pom1 = sy1;
                     pom2 = sx1;
                     sdy2 = pom1-sdx2*pom2;
                     pom1 = sy3-sy2;
                     pom2 = sx3-sx2;
                     if (pom2!=0) sdx3 = pom1/pom2;
                     pom1 = sy2;
                     pom2 = sx2;
                     sdy3 = pom1-sdx3*pom2;
                     if (sx2<sx3) {
                        if (sx1!=sx2) {
                           for (sx4=sx1;sx4<=sx2;sx4++) {
                              pom1 = sx4;
                              sdy4 = sdx1*pom1+sdy1;
                              sy4  = (Int_t)sdy4;
                              if (sx3!=sx1) {
                                 sdy4 = sdx2*pom1+sdy2;
                                 sy5  = (Int_t)sdy4;
                                 y5   = fEnvelope[sx4];
                                 if (sy4<sy5) {
                                    pom1 = sy4;
                                    sy4  = sy5;
                                    sy5  = (Int_t)pom1;
                                 }
                                 if ((sy4<=y5)||(sy5<y5)) {
                                    sy4 = (Int_t)TMath::Min(sy4,(Int_t)y5);
                                    sy5 = (Int_t)TMath::Min(sy5,(Int_t)y5);
                                    line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx4),gPad->PixeltoY(sy5)+1);
                                 }
                              }
                           }
                        }
                        if (sx2!=sx3) {
                           for (sx4=sx2;sx4<=sx3;sx4++) {
                              pom1 = sx4;
                              sdy4 = sdx3*pom1+sdy3;
                              sy4  = (Int_t)sdy4;
                              if (sx3!=sx1) {
                                 sdy4 = sdx2*pom1+sdy2;
                                 sy5  = (Int_t)sdy4;
                                 y5   = fEnvelope[sx4];
                                 if (sy4<sy5) {
                                    pom1 = sy4;
                                    sy4  = sy5;
                                    sy5  = (Int_t)pom1;
                                 }
                                 if ((sy4<=y5)||(sy5<y5)) {
                                    sy4 = (Int_t)TMath::Min(sy4,(Int_t)y5);
                                    sy5 = (Int_t)TMath::Min(sy5,(Int_t)y5);
                                    line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx4),gPad->PixeltoY(sy5)+1);
                                 }
                              }
                           }
                        }
                     } else {
                        if (sx3!=sx1) {
                           for (sx4=sx1;sx4<=sx3;sx4++) {
                              pom1 = sx4;
                              sdy4 = sdx2*pom1+sdy2;
                              sy4  = (Int_t)sdy4;
                              if (sx2!=sx1) {
                                 sdy4 = sdx1*pom1+sdy1;
                                 sy5  = (Int_t)sdy4;
                                 y5 = fEnvelope[sx4];
                                 if (sy4<sy5) {
                                    pom1 = sy4;
                                    sy4  = sy5;
                                    sy5  = (Int_t)pom1;
                                 }
                                 if ((sy4<=y5)||(sy5<y5)) {
                                    sy4 = (Int_t)TMath::Min(sy4,(Int_t)y5);
                                    sy5 = (Int_t)TMath::Min(sy5,(Int_t)y5);
                                    line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx4),gPad->PixeltoY(sy5)+1);
                                 }
                              }
                           }
                        }
                        if (sx2!=sx3) {
                           for (sx4=sx3;sx4<=sx2;sx4++) {
                              pom1 = sx4;
                              sdy4 = sdx3*pom1+sdy3;
                              sy4  = (Int_t)sdy4;
                              if (sx2!=sx1) {
                                 sdy4 = sdx1*pom1+sdy1;
                                 sy5  = (Int_t)sdy4;
                                 y5   = fEnvelope[sx4];
                                 if (sy4<sy5) {
                                    pom1 = sy4;
                                    sy4  = sy5;
                                    sy5  = (Int_t)pom1;
                                 }
                                 if ((sy4<=y5)||(sy5<y5)) {
                                    sy4  = (Int_t)TMath::Min(sy4,(Int_t)y5);
                                    sy5  = (Int_t)TMath::Min(sy5,(Int_t)y5);
                                    line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx4),gPad->PixeltoY(sy5)+1);
                                 }
                              }
                           }
                        }
                     }
                     line->SetLineColor(uip);
                  }
                  EnvelopeBars(x1,y1,x2,y2);
                  if (fLine!=0) {
                     if (fLine==1) {
                        fXe = x2;
                        fYe = y2;
                     }
                     line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                  }
                  EnvelopeBars(x2,y2,x3,y3);
                  if (fLine!=0) {
                     if (fLine==1) {
                        fXe = x3;
                        fYe = y3;
                     }
                     line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                  }
                  EnvelopeBars(x1,y1,x4,y4);
                  if (fLine!=0) {
                     if (fLine==1) {
                        fXe = x4;
                        fYe = y4;
                     }
                     line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                  }
                  EnvelopeBars(x4,y4,x3,y3);
                  if (fLine!=0) {
                     if (fLine==1) {
                        fXe = x3;
                        fYe = y3;
                     }
                     line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                  }
               } else if ((fDisplayMode==kDisplayModeBarsY)&&(((flag!=0)&&(smer==0))||((flag==0)&&(smer!=0)))) {
                  EnvelopeBars(ix5,iy5,x6,y6);
                  if (fLine!=0) {
                     if (fLine==1) {
                        fXe = x6;
                        fYe = y6;
                     }
                     line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                  }
                  EnvelopeBars(x1,y1,ix5,iy5);
                  if (fLine!=0) {
                     if (fLine==1) {
                        fXe = ix5;
                        fYe = iy5;
                     }
                     line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                  }
                  EnvelopeBars(x2,y2,x6,y6);
                  if (fLine!=0) {
                     if (fLine==1) {
                        fXe = x6;
                        fYe = y6;
                     }
                     line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                  }
                  if (fModeGroup!=kModeGroupSimple) {
                     v   = ColorCalculation(dx1,dy1,z4,dx2,dy2,z4,dx4,dy4,z4);
                     v   = v*fLevels+0.5;
                     iv  = fLevels-(Int_t)v;
                     uip = fNewColorIndex;
                     ColorModel(iv,ui1,ui2,ui3);
                     line->SetLineColor(fNewColorIndex);
                  }
                  EnvelopeBars(x1,y1,x2,y2);
                  if (fLine!=0) {
                     if (fLine==1) {
                        fXe = x2;
                        fYe = y2;
                     }
                     line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                  }
                  if (fModeGroup!=kModeGroupSimple) {
                     line->SetLineColor(uip);
                  }
               } else if ((fDisplayMode==kDisplayModeBarsX)&&(((flag==0)&&(smer==0))||((flag!=0)&&(smer!=0)))) {
                  EnvelopeBars(x7,y7,x6,y6);
                  if (fLine!=0) {
                     if (fLine==1) {
                        fXe = x6;
                        fYe = y6;
                     }
                     line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                  }
                  EnvelopeBars(x2,y2,x6,y6);
                  if (fLine!=0) {
                     if (fLine==1) {
                        fXe = x6;
                        fYe = y6;
                     }
                     line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                  }
                  EnvelopeBars(x3,y3,x7,y7);
                  if (fLine!=0) {
                     if (fLine==1) {
                        fXe = x7;
                        fYe = y7;
                     }
                     line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                  }
                  if (fModeGroup!=kModeGroupSimple) {
                     v   = ColorCalculation(dx1,dy1,z4,dx2,dy2,z4,dx4,dy4,z4);
                     v   = v*fLevels+0.5;
                     iv  = fLevels-(Int_t)v;
                     uip = fNewColorIndex;
                     ColorModel(iv,ui1,ui2,ui3);
                     line->SetLineColor(fNewColorIndex);
                  }
                  EnvelopeBars(x3,y3,x2,y2);
                  if (fLine!=0) {
                     if (fLine==1) {
                        fXe = x2;
                        fYe = y2;
                     }
                     line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                  }
                  if (fModeGroup!=kModeGroupSimple) {
                     line->SetLineColor(uip);
                  }
               }
               break;
            case kDisplayModeLinesX:
               if (fModeGroup==kModeGroupSimple) {
                  if (((flag==0)&&(smer==0))||((flag!=0)&&(smer!=0))) {
                     if (fBezier==kNoBezierInterpol) {
                        Envelope(x1,y1,x2,y2);
                        if (fLine!=0) {
                           if (fLine==1) {
                              fXe = x2;
                              fYe = y2;
                           }
                           line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                        }
                     } else {
                        s1 = q1;
                        t1 = (Int_t)TMath::Max(q2-1,0);
                        s2 = q1;
                        t2 = (Int_t)TMath::Min(q2+2,w2);
                        if (turni==1) {
                           i1 = s1;
                           i2 = s2;
                        } else {
                           i1 = (Int_t)TMath::Max(w1-s1,0);
                           i2 = (Int_t)TMath::Max(w1-s2,0);
                        }
                        if (turnj==1) {
                           j1 = t1;
                           j2 = t2;
                        } else {
                           j1 = (Int_t)TMath::Max(w2-t1,0);
                           j2 = (Int_t)TMath::Max(w2-t2,0);
                        }
                        Transform(i1,j1,0);
                        x3 = fXt;
                        y3 = fYt;
                        Transform(i2,j2,0);
                        x4    = fXt;
                        y4    = fYt;
                        bezx1 = x1+(x2-x1)/3;
                        bezx2 = x1+2*(x2-x1)/3;
                        bezy1 = y1+(y2-y3)/6;
                        bezy2 = y2-(y4-y1)/6;
                        if (x1<=x2) {
                           if (bezx1<=x1) {
                              bezx1 = x1;
                              bezy1 = y1;
                           }
                           if (bezx1>=x2) {
                              bezx1 = x2;
                              bezy1 = y2;
                           }
                           if (bezx2<=x1) {
                              bezx2 = x1;
                              bezy2 = y1;
                           }
                           if (bezx2>=x2) {
                              bezx2 = x2;
                              bezy2 = y2;
                           }
                           fBzX[0] = x1;
                           fBzY[0] = y1;
                           fBzX[1] = (Int_t)bezx1;
                           fBzY[1] = (Int_t)bezy1;
                           fBzX[2] = (Int_t)bezx2;
                           fBzY[2] = (Int_t)bezy2;
                           fBzX[3] = x2;
                           fBzY[3] = y2;
                           for (bezf=0;bezf<1.01;bezf+=0.1) {
                              BezierSmoothing(bezf);
                              if (bezf==0) {
                                 ibezx1 = (Int_t)(fGbezx+0.5);
                                 ibezy1 = (Int_t)(fGbezy+0.5);
                              } else {
                                 ibezx2 = ibezx1;
                                 ibezy2 = ibezy1;
                                 ibezx1 = (Int_t)(fGbezx+0.5);
                                 ibezy1 = (Int_t)(fGbezy+0.5);
                                 Envelope(ibezx2,ibezy2,ibezx1,ibezy1);
                                 if (fLine!=0) {
                                    if (fLine==1) {
                                       fXe = ibezx1;
                                       fYe = ibezy1;
                                    }
                                    line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                                 }
                              }
                           }
                        } else if (x1>x2) {
                           if (bezx1>=x1) {
                              bezx1 = x1;
                              bezy1 = y1;
                           }
                           if (bezx1<=x2) {
                              bezx1 = x2;
                              bezy1 = y2;
                           }
                           if (bezx2>=x1) {
                              bezx2 = x1;
                              bezy2 = y1;
                           }
                           if (bezx2<=x2) {
                              bezx2 = x2;
                              bezy2 = y2;
                           }
                           fBzX[0] = x1;
                           fBzY[0] = y1;
                           fBzX[1] = (Int_t)bezx1;
                           fBzY[1] = (Int_t)bezy1;
                           fBzX[2] = (Int_t)bezx2;
                           fBzY[2] = (Int_t)bezy2;
                           fBzX[3] = x2;
                           fBzY[3] = y2;
                           for (bezf=0;bezf<1.01;bezf+=0.1) {
                              BezierSmoothing(bezf);
                              if (bezf==0) {
                                 ibezx1 = (Int_t)(fGbezx+0.5);
                                 ibezy1 = (Int_t)(fGbezy+0.5);
                              } else {
                                 ibezx2 = ibezx1;
                                 ibezy2 = ibezy1;
                                 ibezx1 = (Int_t)(fGbezx+0.5);
                                 ibezy1 = (Int_t)(fGbezy+0.5);
                                 Envelope(ibezx1,ibezy1,ibezx2,ibezy2);
                                 if (fLine!=0) {
                                    if (fLine==1) {
                                       fXe = ibezx2;
                                       fYe = ibezy2;
                                    }
                                    line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                                 }
                              }
                           }
                        }
                     }
                  }
               } else {
                  if ((q1!=q2||smer!=0)&&flag==1) {
                     s1 = q1+1;
                     t1 = q2;
                     s2 = q1;
                     t2 = q2;
                     s3 = q1;
                     t3 = q2+1;
                     s4 = q1+1;
                     t4 = q2+1;
                     if (fShading==kShaded) {
                        sr1 = s1;
                        tr1 = (Int_t)TMath::Max(t1-1,0);
                        sr2 = s2;
                        tr2 = (Int_t)TMath::Max(t2-1,0);
                        sr3 = (Int_t)TMath::Max(s2-1,0);
                        tr3 = t2;
                        sr4 = (Int_t)TMath::Max(s3-1,0);
                        tr4 = t3;
                        sr5 = s3;
                        tr5 = t3+1;
                        sr6 = s4;
                        tr6 = t4+1;
                        sr7 = s4+1;
                        tr7 = t4;
                        sr8 = s1+1;
                        tr8 = t1;
                     }
                     if (turni==1) {
                        i1 = s1;
                        i2 = s2;
                        i3 = s3;
                        i4 = s4;
                     } else {
                        i1 = (Int_t)TMath::Max(w1-s1,0);
                        i2 = (Int_t)TMath::Max(w1-s2,0);
                        i3 = (Int_t)TMath::Max(w1-s3,0);
                        i4 = (Int_t)TMath::Max(w1-s4,0);
                        if (fShading==kShaded) {
                           sr1 = (Int_t)TMath::Max(w1-sr1,0);
                           sr2 = (Int_t)TMath::Max(w1-sr2,0);
                           sr3 = (Int_t)TMath::Max(w1-sr3,0);
                           sr4 = (Int_t)TMath::Max(w1-sr4,0);
                           sr5 = (Int_t)TMath::Max(w1-sr5,0);
                           sr6 = (Int_t)TMath::Max(w1-sr6,0);
                           sr7 = (Int_t)TMath::Max(w1-sr7,0);
                           sr8 = (Int_t)TMath::Max(w1-sr8,0);
                        }
                     }
                     if (turnj==1) {
                        j1 = t1;
                        j2 = t2;
                        j3 = t3;
                        j4 = t4;
                     } else {
                        j1 = (Int_t)TMath::Max(w2-t1,0);
                        j2 = (Int_t)TMath::Max(w2-t2,0);
                        j3 = (Int_t)TMath::Max(w2-t3,0);
                        j4 = (Int_t)TMath::Max(w2-t4,0);
                        if (fShading==kShaded) {
                           tr1 = (Int_t)TMath::Max(w2-tr1,0);
                           tr2 = (Int_t)TMath::Max(w2-tr2,0);
                           tr3 = (Int_t)TMath::Max(w2-tr3,0);
                           tr4 = (Int_t)TMath::Max(w2-tr4,0);
                           tr5 = (Int_t)TMath::Max(w2-tr5,0);
                           tr6 = (Int_t)TMath::Max(w2-tr6,0);
                           tr7 = (Int_t)TMath::Max(w2-tr7,0);
                           tr8 = (Int_t)TMath::Max(w2-tr8,0);
                        }
                     }
                     Transform(i1,j1,0);
                     x1  = fXt;
                     y1  = fYt;
                     dx1 = fDxspline;
                     dy1 = fDyspline;
                     z1  = fZ;
                     Transform(i2,j2,0);
                     x2  = fXt;
                     y2  = fYt;
                     dx2 = fDxspline;
                     dy2 = fDyspline;
                     z2  = fZ;
                     Transform(i3,j3,0);
                     x3  = fXt;
                     y3  = fYt;
                     dx3 = fDxspline;
                     dy3 = fDyspline;
                     z3  = fZ;
                     Transform(i4,j4,0);
                     x4  = fXt;
                     y4  = fYt;
                     dx4 = fDxspline;
                     dy4 = fDyspline;
                     z4  = fZ;
                     Envelope(x1,y1,x2,y2);
                     Envelope(x2,y2,x3,y3);
                     xtaz = (dx1+dx2+dx4)/3;
                     ytaz = (dy1+dy2+dy4)/3;
                     ztaz = (z1+z2+z4)/3;
                     v    = ColorCalculation(dx1,dy1,z1,dx2,dy2,z2,dx4,dy4,z4);
                     if (fShading==kShaded) {
                        if (fShadow==kShadowsNotPainted) {
                           if (sr1<0||sr1>w1||tr1<0||tr1>w2) Transform(sr1,tr1,-1);
                           else                              Transform(sr1,tr1,0);
                           dxr1 = fDxspline;
                           dyr1 = fDyspline;
                           zr1  = fZ;
                           if (sr8<0||sr8>w1||tr8<0||tr8>w2) Transform(sr8,tr8,-1);
                           else                              Transform(sr8,tr8,0);
                           dxr2 = fDxspline;
                           dyr2 = fDyspline;
                           zr2  = fZ;
                           v    = v+ColorCalculation(dxr1,dyr1,zr1,dx2,dy2,z2,dx1,dy1,z1);
                           v    = v+ColorCalculation(dxr1,dyr1,zr1,dx1,dy1,z1,dxr2,dyr2,zr2);
                           v    = v+ColorCalculation(dxr2,dyr2,zr2,dx1,dy1,z1,dx4,dy4,z4);
                           v1   = v/4;
                           if (sr3<0||sr3>w1||tr3<0||tr3>w2) Transform(sr3,tr3,-1);
                           else                              Transform(sr3,tr3,0);
                           dxr1 = fDxspline;
                           dyr1 = fDyspline;
                           zr1  = fZ;
                           if (sr2<0||sr2>w1||tr2<0||tr2>w2) Transform(sr2,tr2,-1);
                           else                              Transform(sr2,tr2,0);
                           dxr2 = fDxspline;
                           dyr2 = fDyspline;
                           zr2  = fZ;
                           v    = ColorCalculation(dx1,dy1,z1,dx2,dy2,z2,dx3,dy3,z3);
                           v    = v+ColorCalculation(dx2,dy2,z2,dxr1,dyr1,zr1,dx3,dy3,z3);
                           v    = v+ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dxr1,dyr1,zr1);
                           v    = v+ColorCalculation(dxr2,dyr2,zr2,dx2,dy2,z2,dx1,dy1,z1);
                           v2   = v/4;
                           if (sr5<0||sr5>w1||tr5<0||tr5>w2) Transform(sr5,tr5,-1);
                           else                              Transform(sr5,tr5,0);
                           dxr1 = fDxspline;
                           dyr1 = fDyspline;
                           zr1  = fZ;
                           if (sr4<0||sr4>w1||tr4<0||tr4>w2) Transform(sr4,tr4,-1);
                           else                              Transform(sr4,tr4,0);
                           dxr2 = fDxspline;
                           dyr2 = fDyspline;
                           zr2  = fZ;
                           v    = ColorCalculation(dx2,dy2,z2,dx3,dy3,z3,dx4,dy4,z4);
                           v    = v+ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr1,dyr1,zr1);
                           v    = v+ColorCalculation(dx3,dy3,z3,dxr2,dyr2,zr2,dxr1,dyr1,zr1);
                           v    = v+ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dx3,dy3,z3);
                           v3   = v/4;
                           if (sr7<0||sr7>w1||tr7<0||tr7>w2) Transform(sr7,tr7,-1);
                           else                              Transform(sr7,tr7,0);
                           dxr1 = fDxspline;
                           dyr1 = fDyspline;
                           zr1  = fZ;
                           if (sr6<0||sr6>w1||tr6<0||tr6>w2) Transform(sr6,tr6,-1);
                           else                              Transform(sr6,tr6,0);
                           dxr2 = fDxspline;
                           dyr2 = fDyspline;
                           zr2  = fZ;
                           v    = ColorCalculation(dx1,dy1,z1,dx3,dy3,z3,dx4,dy4,z4);
                           v    = v+ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr2,dyr2,zr2);
                           v    = v+ColorCalculation(dx4,dy4,z4,dxr2,dyr2,zr2,dxr1,dyr1,zr1);
                           v    = v+ColorCalculation(dx1,dy1,z1,dx4,dy4,z4,dxr1,dyr1,zr1);
                           v4   = v/4;
                        } else {
                           spriz = ShadowColorCalculation(xtaz,ytaz,ztaz,shad_noise);
                           v     = v+spriz;
                           v     = v/2;
                           if (sr1<0||sr1>w1||tr1<0||tr1>w2) Transform(sr1,tr1,-1);
                           else                              Transform(sr1,tr1,0);
                           dxr1  = fDxspline;
                           dyr1  = fDyspline;
                           zr1   = fZ;
                           if (sr8<0||sr8>w1||tr8<0||tr8>w2) Transform(sr8,tr8,-1);
                           else                              Transform(sr8,tr8,0);
                           dxr2  = fDxspline;
                           dyr2  = fDyspline;
                           zr2   = fZ;
                           da    = (dxr1+dx2+dx1)/3;
                           db    = (dyr1+dy2+dy1)/3;
                           dc    = (zr1+z2+z1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dxr1,dyr1,zr1,dx2,dy2,z2,dx1,dy1,z1)+spriz)/2;
                           da    = (dxr1+dxr2+dx1)/3;
                           db    = (dyr1+dyr2+dy1)/3;
                           dc    = (zr1+zr2+z1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dxr1,dyr1,zr1,dx1,dy1,z1,dxr2,dyr2,zr2)+spriz)/2;
                           da    = (dxr2+dx1+dx4)/3;
                           db    = (dyr2+dy1+dy4)/3;
                           dc    = (zr2+z1+z4)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dxr2,dyr2,zr2,dx1,dy1,z1,dx4,dy4,z4)+spriz)/2;
                           v1    = v/4;
                           if (sr3<0||sr3>w1||tr3<0||tr3>w2) Transform(sr3,tr3,-1);
                           else                              Transform(sr3,tr3,0);
                           dxr1  = fDxspline;
                           dyr1  = fDyspline;
                           zr1   = fZ;
                           if (sr2<0||sr2>w1||tr2<0||tr2>w2) Transform(sr2,tr2,-1);
                           else                              Transform(sr2,tr2,0);
                           dxr2  = fDxspline;
                           dyr2  = fDyspline;
                           zr2   = fZ;
                           da    = (dx1+dx2+dx3)/3;
                           db    = (dy1+dy2+dy3)/3;
                           dc    = (z1+z2+z3)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = (ColorCalculation(dx1,dy1,z1,dx2,dy2,z2,dx3,dy3,z3)+spriz)/2;
                           da    = (dx2+dxr1+dx3)/3;
                           db    = (dy2+dyr1+dy3)/3;
                           dc    = (z2+zr1+z3)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx2,dy2,z2,dxr1,dyr1,zr1,dx3,dy3,z3)+spriz)/2;
                           da    = (dx2+dxr2+dxr1)/3;
                           db    = (dy2+dyr2+dyr1)/3;
                           dc    = (z2+zr2+zr1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dxr1,dyr1,zr1)+spriz)/2;
                           da    = (dxr2+dx2+dx1)/3;
                           db    = (dyr2+dy2+dy1)/3;
                           dc    = (zr2+z2+z1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dxr2,dyr2,zr2,dx2,dy2,z2,dx1,dy1,z1)+spriz)/2;
                           v2    = v/4;
                           if (sr5<0||sr5>w1||tr5<0||tr5>w2) Transform(sr5,tr5,-1);
                           else                              Transform(sr5,tr5,0);
                           dxr1  = fDxspline;
                           dyr1  = fDyspline;
                           zr1   = fZ;
                           if (sr4<0||sr4>w1||tr4<0||tr4>w2) Transform(sr4,tr4,-1);
                           else                              Transform(sr4,tr4,0);
                           dxr2  = fDxspline;
                           dyr2  = fDyspline;
                           zr2   = fZ;
                           da    = (dx2+dx3+dx4)/3;
                           db    = (dy2+dy3+dy4)/3;
                           dc    = (z2+z3+z4)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = (ColorCalculation(dx2,dy2,z2,dx3,dy3,z3,dx4,dy4,z4)+spriz)/2;
                           da    = (dx4+dx3+dxr1)/3;
                           db    = (dy4+dy3+dyr1)/3;
                           dc    = (z4+z3+zr1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr1,dyr1,zr1)+spriz)/2;
                           da    = (dx3+dxr2+dxr1)/3;
                           db    = (dy3+dyr2+dyr1)/3;
                           dc    = (z3+zr2+zr1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx3,dy3,z3,dxr2,dyr2,zr2,dxr1,dyr1,zr1)+spriz)/2;
                           da    = (dx2+dxr2+dx3)/3;
                           db    = (dy2+dyr2+dy3)/3;
                           dc    = (z2+zr2+z3)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dx3,dy3,z3)+spriz)/2;
                           v3    = v/4;
                           if (sr7<0||sr7>w1||tr7<0||tr7>w2) Transform(sr7,tr7,-1);
                           else                              Transform(sr7,tr7,0);
                           dxr1  = fDxspline;
                           dyr1  = fDyspline;
                           zr1   = fZ;
                           if (sr6<0||sr6>w1||tr6<0||tr6>w2) Transform(sr6,tr6,-1);
                           else                              Transform(sr6,tr6,0);
                           dxr2  = fDxspline;
                           dyr2  = fDyspline;
                           zr2   = fZ;
                           da    = (dx1+dx3+dx4)/3;
                           db    = (dy1+dy3+dy4)/3;
                           dc    = (z1+z3+z4)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = (ColorCalculation(dx1,dy1,z1,dx3,dy3,z3,dx4,dy4,z4)+spriz)/2;
                           da    = (dx4+dx3+dxr2)/3;
                           db    = (dy4+dy3+dyr2)/3;
                           dc    = (z4+z3+zr2)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr2,dyr2,zr2)+spriz)/2;
                           da    = (dx4+dxr2+dxr1)/3;
                           db    = (dy4+dyr2+dyr1)/3;
                           dc    = (z4+zr2+zr1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx4,dy4,z4,dxr2,dyr2,zr2,dxr1,dyr1,zr1)+spriz)/2;
                           da    = (dx1+dx4+dxr1)/3;
                           db    = (dy1+dy4+dyr1)/3;
                           dc    = (z1+z4+zr1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx1,dy1,z1,dx4,dy4,z4,dxr1,dyr1,zr1)+spriz)/2;
                           v4    = v/4;
                        }
                     }
                     spriz = 0;
                     if (fShadow==kShadowsNotPainted) {
                        if (fShading==kNotShaded) {
                           v  = v*fLevels+0.5;
                           iv = fLevels-(Int_t)v;
                        } else {
                           v1  = v1*fLevels;
                           iv1 = fLevels-(Int_t)v1;
                           v2  = v2*fLevels;
                           iv2 = fLevels-(Int_t)v2;
                           v4  = v4*fLevels;
                           iv4 = fLevels-(Int_t)v4;
                        }
                     } else {
                     spriz = ShadowColorCalculation(xtaz,ytaz,ztaz,shad_noise);
                        if (fShading==kNotShaded) {
                           v  = v*fLevels/2.0;
                           iv = fLevels-(Int_t)(v+0.5);
                        } else {
                           v1  = v1*fLevels;
                           iv1 = fLevels-(Int_t)v1;
                           v2  = v2*fLevels;
                           iv2 = fLevels-(Int_t)v2;
                           v4  = v4*fLevels;
                           iv4 = fLevels-(Int_t)v4;
                        }
                     }
                     if (fShading==kNotShaded) {
                        ColorModel(iv,ui1,ui2,ui3);
                        line->SetLineColor(fNewColorIndex);
                     } else {
                        dx1 = x1;
                        dy1 = y1;
                        dx2 = x2;
                        dy2 = y2;
                        dx3 = x4;
                        dy3 = y4;
                        z1  = iv1;
                        z2  = iv2;
                        z3  = iv4;
                        da  = (dy2-dy1)*(z3-z1)-(dy3-dy1)*(z2-z1);
                        db  = (z2-z1)*(dx3-dx1)-(z3-z1)*(dx2-dx1);
                        dc  = (dx2-dx1)*(dy3-dy1)-(dx3-dx1)*(dy2-dy1);
                        dd  = -da*dx1-db*dy1-dc*z1;
                     }
                     sx1 = x1;
                     sy1 = y1;
                     sx3 = x4;
                     sy3 = y4;
                     if (sx3<sx1) {
                        sx4 = sx1;
                        sy4 = sy1;
                        sx1 = sx3;
                        sy1 = sy3;
                        sx3 = sx4;
                        sy3 = sy4;
                     }
                     pom1 = sy3-sy1;
                     pom2 = sx3-sx1;
                     if (pom2!=0) sdx2 = pom1/pom2;
                     pom1  = sy1;
                     pom2  = sx1;
                     sdy2  = pom1-sdx2*pom2;
                     sx1p  = sx1;
                     sy1p  = sy1;
                     sx3p  = sx3;
                     sdx2p = sdx2;
                     sdy2p = sdy2;
                     dap   = da;
                     dbp   = db;
                     dcp   = dc;
                     ddp   = dd;
                     uip   = fNewColorIndex;
                     xtaz  = (dx3+dx2+dx4)/3;
                     ytaz  = (dy3+dy2+dy4)/3;
                     ztaz  = (z3+z2+z4)/3;
                     if (fShading==kNotShaded) v=ColorCalculation(dx2,dy2,z2,dx3,dy3,z3,dx4,dy4,z4);
                     spriz = 0;
                     if (fShadow==kShadowsNotPainted) {
                        if (fShading==kNotShaded) {
                           v  = v*fLevels;
                           iv = fLevels-(Int_t)v;
                        } else {
                           v3  = v3*fLevels;
                           iv3 = fLevels-(Int_t)v3;
                        }
                     } else {
                        spriz = ShadowColorCalculation(xtaz,ytaz,ztaz,shad_noise);
                        if (fShading==kNotShaded) {
                           v  = v*fLevels/2;
                           iv = fLevels-(Int_t)v;
                           iv = (Int_t)(iv-fLevels*spriz/2);
                        } else {
                           v3  = v3*fLevels;
                           iv3 = fLevels-(Int_t)v3;
                        }
                     }
                     if (fShading==kNotShaded) {
                        ColorModel(iv,ui1,ui2,ui3);
                        line->SetLineColor(fNewColorIndex);
                     } else {
                        dx1 = x2;
                        dy1 = y2;
                        dx2 = x3;
                        dy2 = y3;
                        dx3 = x4;
                        dy3 = y4;
                        z1  = iv2;
                        z2  = iv3;
                        z3  = iv4;
                        da  = (dy2-dy1)*(z3-z1)-(dy3-dy1)*(z2-z1);
                        db  = (z2-z1)*(dx3-dx1)-(z3-z1)*(dx2-dx1);
                        dc  = (dx2-dx1)*(dy3-dy1)-(dx3-dx1)*(dy2-dy1);
                        dd  = -da*dx1-db*dy1-dc*z1;
                     }
                     sx1 = x2;
                     sy1 = y2;
                     sx2 = x3;
                     sy2 = y3;
                     if (sx2<sx1) {
                        sx4 = sx1;
                        sy4 = sy1;
                        sx1 = sx2;
                        sy1 = sy2;
                        sx2 = sx4;
                        sy2 = sy4;
                     }
                     pom1 = sy2-sy1;
                     pom2 = sx2-sx1;
                     sdx1 = 0;
                     if (pom2!=0) sdx1 = pom1/pom2;
                     pom1 = sy1;
                     pom2 = sx1;
                     sdy1 = pom1-sdx1*pom2;
                     for (sx4=sx1,sx5=sx1,sy5=sy1;sx4<=sx2;sx4++) {
                        pom1 = sx4;
                        sdy4 = sdx1*pom1+sdy1;
                        sy4  = (Int_t)sdy4;
                        if (sy4<=fEnvelope[sx4]) {
                           fEnvelope[sx4] = sy4;
                           if (fShading==kNotShaded) {
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           } else {
                              dx1 = sx4;
                              dy1 = sy4;
                              if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                              else       v =(iv1+iv2+iv4)/3;
                              iv = (Int_t)v;
                              ColorModel(iv,ui1,ui2,ui3);
                              line->SetLineColor(fNewColorIndex);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           }
                           sy5 = sy4;
                        } else {
                           sy4 = fEnvelope[sx4];
                           if (fShading==kNotShaded&&sy5<=fEnvelope[sx5]) {
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           } else if (sy5<=fEnvelope[sx5]) {
                              dx1 = sx4;
                              dy1 = sy4;
                              if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                              else       v =(iv1+iv2+iv4)/3;
                              iv = (Int_t)v;
                              ColorModel(iv,ui1,ui2,ui3);
                              line->SetLineColor(fNewColorIndex);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           }
                           sy5 = fEnvelope[sx4];
                        }
                        sx5 = sx4;
                     }
                     for (sx4=sx1p,sx5=sx1p,sy5=sy1p;sx4<=sx3p;sx4++) {
                        pom1 = sx4;
                        sdy4 = sdx2p*pom1+sdy2p;
                        sy4  = (Int_t)sdy4;
                        if (sy4<=fEnvelope[sx4]) {
                           fEnvelope[sx4]=sy4;
                           if (fShading==kNotShaded) {
                              line->SetLineColor(uip);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           } else {
                              dx1 = sx4;
                              dy1 = sy4;
                              if (dcp!=0) v = (-ddp-dap*dx1-dbp*dy1)/dcp;
                              else        v = (iv1+iv2+iv4)/3;
                              iv = (Int_t)v;
                              ColorModel(iv,ui1,ui2,ui3);
                              line->SetLineColor(fNewColorIndex);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           }
                           sy5 = sy4;
                        } else {
                           sy4 = fEnvelope[sx4];
                           if (fShading==kNotShaded&&sy5<=fEnvelope[sx5]) {
                              line->SetLineColor(uip);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           } else if (sy5<=fEnvelope[sx5]) {
                              dx1 = sx4;
                              dy1 = sy4;
                              if (dcp!=0) v = (-ddp-dap*dx1-dbp*dy1)/dcp;
                              else        v = (iv1+iv2+iv4)/3;
                              iv = (Int_t)v;
                              ColorModel(iv,ui1,ui2,ui3);
                              line->SetLineColor(fNewColorIndex);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           }
                           sy5 = fEnvelope[sx4];
                        }
                        sx5 = sx4;
                     }
                  }
               }
               break;
            case kDisplayModeLinesY:
               if (fModeGroup==kModeGroupSimple) {
                  if (((flag!=0)&&(smer==0))||((flag==0)&&(smer!=0))) {
                     if (fBezier==kNoBezierInterpol) {
                        Envelope(x1,y1,x2,y2);
                        if (fLine!=0) {
                           if (fLine==1) {
                              fXe = x2;
                              fYe = y2;
                           }
                           line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                        }
                     } else {
                        s1 = (Int_t)TMath::Max(q1-1,0);
                        t1 = q2;
                        s2 = (Int_t)TMath::Min(q1+2,w1);
                        t2 = q2;
                        if (turni==1) {
                           i1 = s1;
                           i2 = s2;
                        } else {
                           i1 = w1-s1;
                           i2 = w1-s2;
                        }
                        if (turnj==1) {
                           j1 = t1;
                           j2 = t2;
                        } else {
                           j1 = w2-t1;
                           j2 = w2-t2;
                        }
                        Transform(i1,j1,0);
                        x3 = fXt;
                        y3 = fYt;
                        Transform(i2,j2,0);
                        x4 = fXt;
                        y4 = fYt;
                        bezx1 = x1+(x2-x1)/3;
                        bezx2 = x1+2*(x2-x1)/3;
                        bezy1 = y1+(y2-y3)/6;
                        bezy2 = y2-(y4-y1)/6;
                        if (x1<=x2) {
                           if (bezx1<=x1) {
                              bezx1 = x1;
                              bezy1 = y1;
                           }
                           if (bezx1>=x2) {
                              bezx1 = x2;
                              bezy1 = y2;
                           }
                           if (bezx2<=x1) {
                              bezx2 = x1;
                              bezy2 = y1;
                           }
                           if (bezx2>=x2) {
                              bezx2 = x2;
                              bezy2 = y2;
                           }
                           fBzX[0] = x1;
                           fBzY[0] = y1;
                           fBzX[1] = (Int_t)bezx1;
                           fBzY[1] = (Int_t)bezy1;
                           fBzX[2] = (Int_t)bezx2;
                           fBzY[2] = (Int_t)bezy2;
                           fBzX[3] = x2;
                           fBzY[3] = y2;
                           for (bezf=0;bezf<1.01;bezf+=0.1) {
                           BezierSmoothing(bezf);
                              if (bezf==0) {
                                 ibezx1 = (Int_t)(fGbezx+0.5);
                                 ibezy1 = (Int_t)(fGbezy+0.5);
                              } else {
                                 ibezx2 = ibezx1;
                                 ibezy2 = ibezy1;
                                 ibezx1 = (Int_t)(fGbezx+0.5);
                                 ibezy1 = (Int_t)(fGbezy+0.5);
                                 Envelope(ibezx2,ibezy2,ibezx1,ibezy1);
                                 if (fLine!=0) {
                                    if (fLine==1) {
                                       fXe = ibezx1;
                                       fYe = ibezy1;
                                    }
                                    line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                                 }
                              }
                           }
                        } else if (x1>x2) {
                           if (bezx1>=x1) {
                              bezx1 = x1;
                              bezy1 = y1;
                           }
                           if (bezx1<=x2) {
                              bezx1 = x2;
                              bezy1 = y2;
                           }
                           if (bezx2>=x1) {
                              bezx2 = x1;
                              bezy2 = y1;
                           }
                           if (bezx2<=x2) {
                              bezx2 = x2;
                              bezy2 = y2;
                           }
                           fBzX[0] = x1;
                           fBzY[0] = y1;
                           fBzX[1] = (Int_t)bezx1;
                           fBzY[1] = (Int_t)bezy1;
                           fBzX[2] = (Int_t)bezx2;
                           fBzY[2] = (Int_t)bezy2;
                           fBzX[3] = x2;
                           fBzY[3] = y2;
                           for (bezf=0;bezf<1.01;bezf+=0.1) {
                              BezierSmoothing(bezf);
                              if (bezf==0) {
                                 ibezx1 = (Int_t)(fGbezx+0.5);
                                 ibezy1 = (Int_t)(fGbezy+0.5);
                              } else {
                              ibezx2 = ibezx1;
                              ibezy2 = ibezy1;
                              ibezx1 = (Int_t)(fGbezx+0.5);
                              ibezy1 = (Int_t)(fGbezy+0.5);
                              Envelope(ibezx1,ibezy1,ibezx2,ibezy2);
                              if (fLine!=0) {
                                 if (fLine==1) {
                                    fXe = ibezx2;
                                    fYe = ibezy2;
                                 }
                                 line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                                 }
                              }
                           }
                        }
                     }
                  }
               } else {
                  if ((q1!=q2||smer!=0)&&flag==1) {
                     s1 = q1+1;
                     t1 = q2;
                     s2 = q1;
                     t2 = q2;
                     s3 = q1;
                     t3 = q2+1;
                     s4 = q1+1;
                     t4 = q2+1;
                     if (fShading==kShaded) {
                        sr1 = s1;
                        tr1 = (Int_t)TMath::Max(t1-1,0);
                        sr2 = s2;
                        tr2 = (Int_t)TMath::Max(t2-1,0);
                        sr3 = (Int_t)TMath::Max(s2-1,0);
                        tr3 = t2;
                        sr4 = (Int_t)TMath::Max(s3-1,0);
                        tr4 = t3;
                        sr5 = s3;
                        tr5 = t3+1;
                        sr6 = s4;
                        tr6 = t4+1;
                        sr7 = s4+1;
                        tr7 = t4;
                        sr8 = s1+1;
                        tr8 = t1;
                     }
                        if (turni==1) {
                        i1 = s1;
                        i2 = s2;
                        i3 = s3;
                        i4 = s4;
                     } else {
                        i1 = (Int_t)TMath::Max(w1-s1,0);
                        i2 = (Int_t)TMath::Max(w1-s2,0);
                        i3 = (Int_t)TMath::Max(w1-s3,0);
                        i4 = (Int_t)TMath::Max(w1-s4,0);
                        if (fShading==kShaded) {
                           sr1 = (Int_t)TMath::Max(w1-sr1,0);
                           sr2 = (Int_t)TMath::Max(w1-sr2,0);
                           sr3 = (Int_t)TMath::Max(w1-sr3,0);
                           sr4 = (Int_t)TMath::Max(w1-sr4,0);
                           sr5 = (Int_t)TMath::Max(w1-sr5,0);
                           sr6 = (Int_t)TMath::Max(w1-sr6,0);
                           sr7 = (Int_t)TMath::Max(w1-sr7,0);
                           sr8 = (Int_t)TMath::Max(w1-sr8,0);
                        }
                     }
                     if (turnj==1) {
                        j1 = t1;
                        j2 = t2;
                        j3 = t3;
                        j4 = t4;
                     } else {
                        j1 = (Int_t)TMath::Max(w2-t1,0);
                        j2 = (Int_t)TMath::Max(w2-t2,0);
                        j3 = (Int_t)TMath::Max(w2-t3,0);
                        j4 = (Int_t)TMath::Max(w2-t4,0);
                        if (fShading==kShaded) {
                           tr1 = (Int_t)TMath::Max(w2-tr1,0);
                           tr2 = (Int_t)TMath::Max(w2-tr2,0);
                           tr3 = (Int_t)TMath::Max(w2-tr3,0);
                           tr4 = (Int_t)TMath::Max(w2-tr4,0);
                           tr5 = (Int_t)TMath::Max(w2-tr5,0);
                           tr6 = (Int_t)TMath::Max(w2-tr6,0);
                           tr7 = (Int_t)TMath::Max(w2-tr7,0);
                           tr8 = (Int_t)TMath::Max(w2-tr8,0);
                        }
                     }
                     Transform(i1,j1,0);
                     x1  = fXt;
                     y1  = fYt;
                     dx1 = fDxspline;
                     dy1 = fDyspline;
                     z1  = fZ;
                     Transform(i2,j2,0);
                     x2  = fXt;
                     y2  = fYt;
                     dx2 = fDxspline;
                     dy2 = fDyspline;
                     z2  = fZ;
                     Transform(i3,j3,0);
                     x3  = fXt;
                     y3  = fYt;
                     dx3 = fDxspline;
                     dy3 = fDyspline;
                     z3  = fZ;
                     Transform(i4,j4,0);
                     x4  = fXt;
                     y4  = fYt;
                     dx4 = fDxspline;
                     dy4 = fDyspline;
                     z4  = fZ;
                     Envelope(x1,y1,x2,y2);
                     Envelope(x2,y2,x3,y3);
                     xtaz = (dx1+dx2+dx4)/3;
                     ytaz = (dy1+dy2+dy4)/3;
                     ztaz = (z1+z2+z4)/3;
                     v    = ColorCalculation(dx1,dy1,z1,dx2,dy2,z2,dx4,dy4,z4);
                     if (fShading==kShaded) {
                        if (fShadow==kShadowsNotPainted) {
                           if (sr1<0||sr1>w1||tr1<0||tr1>w2) Transform(sr1,tr1,-1);
                           else                              Transform(sr1,tr1,0);
                           dxr1 = fDxspline;
                           dyr1 = fDyspline;
                           zr1  = fZ;
                           if (sr8<0||sr8>w1||tr8<0||tr8>w2) Transform(sr8,tr8,-1);
                           else                              Transform(sr8,tr8,0);
                           dxr2 = fDxspline;
                           dyr2 = fDyspline;
                           zr2  = fZ;
                           v    = v+ColorCalculation(dxr1,dyr1,zr1,dx2,dy2,z2,dx1,dy1,z1);
                           v    = v+ColorCalculation(dxr1,dyr1,zr1,dx1,dy1,z1,dxr2,dyr2,zr2);
                           v    = v+ColorCalculation(dxr2,dyr2,zr2,dx1,dy1,z1,dx4,dy4,z4);
                           v1   = v/4;
                           if (sr3<0||sr3>w1||tr3<0||tr3>w2) Transform(sr3,tr3,-1);
                           else                              Transform(sr3,tr3,0);
                           dxr1 = fDxspline;
                           dyr1 = fDyspline;
                           zr1  = fZ;
                           if (sr2<0||sr2>w1||tr2<0||tr2>w2) Transform(sr2,tr2,-1);
                           else                              Transform(sr2,tr2,0);
                           dxr2 = fDxspline;
                           dyr2 = fDyspline;
                           zr2  = fZ;
                           v    = ColorCalculation(dx1,dy1,z1,dx2,dy2,z2,dx3,dy3,z3);
                           v    = v+ColorCalculation(dx2,dy2,z2,dxr1,dyr1,zr1,dx3,dy3,z3);
                           v    = v+ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dxr1,dyr1,zr1);
                           v    = v+ColorCalculation(dxr2,dyr2,zr2,dx2,dy2,z2,dx1,dy1,z1);
                           v2   = v/4;
                           if (sr5<0||sr5>w1||tr5<0||tr5>w2) Transform(sr5,tr5,-1);
                           else                              Transform(sr5,tr5,0);
                           dxr1 = fDxspline;
                           dyr1 = fDyspline;
                           zr1  = fZ;
                           if (sr4<0||sr4>w1||tr4<0||tr4>w2) Transform(sr4,tr4,-1);
                           else                              Transform(sr4,tr4,0);
                           dxr2 = fDxspline;
                           dyr2 = fDyspline;
                           zr2  = fZ;
                           v    = ColorCalculation(dx2,dy2,z2,dx3,dy3,z3,dx4,dy4,z4);
                           v    = v+ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr1,dyr1,zr1);
                           v    = v+ColorCalculation(dx3,dy3,z3,dxr2,dyr2,zr2,dxr1,dyr1,zr1);
                           v    = v+ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dx3,dy3,z3);
                           v3   = v/4;
                           if (sr7<0||sr7>w1||tr7<0||tr7>w2) Transform(sr7,tr7,-1);
                           else                              Transform(sr7,tr7,0);
                           dxr1 = fDxspline;
                           dyr1 = fDyspline;
                           zr1  = fZ;
                           if (sr6<0||sr6>w1||tr6<0||tr6>w2) Transform(sr6,tr6,-1);
                           else                              Transform(sr6,tr6,0);
                           dxr2 = fDxspline;
                           dyr2 = fDyspline;
                           zr2  = fZ;
                           v    = ColorCalculation(dx1,dy1,z1,dx3,dy3,z3,dx4,dy4,z4);
                           v    = v+ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr2,dyr2,zr2);
                           v    = v+ColorCalculation(dx4,dy4,z4,dxr2,dyr2,zr2,dxr1,dyr1,zr1);
                           v    = v+ColorCalculation(dx1,dy1,z1,dx4,dy4,z4,dxr1,dyr1,zr1);
                           v4   = v/4;
                        } else {
                           spriz = ShadowColorCalculation(xtaz,ytaz,ztaz,shad_noise);
                           v     = v+spriz;
                           v     = v/2;
                           if (sr1<0||sr1>w1||tr1<0||tr1>w2) Transform(sr1,tr1,-1);
                           else                              Transform(sr1,tr1,0);
                           dxr1 = fDxspline;
                           dyr1 = fDyspline;
                           zr1  = fZ;
                           if (sr8<0||sr8>w1||tr8<0||tr8>w2) Transform(sr8,tr8,-1);
                           else                              Transform(sr8,tr8,0);
                           dxr2  = fDxspline;
                           dyr2  = fDyspline;
                           zr2   = fZ;
                           da    = (dxr1+dx2+dx1)/3;
                           db    = (dyr1+dy2+dy1)/3;
                           dc    = (zr1+z2+z1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dxr1,dyr1,zr1,dx2,dy2,z2,dx1,dy1,z1)+spriz)/2;
                           da    = (dxr1+dxr2+dx1)/3;
                           db    = (dyr1+dyr2+dy1)/3;
                           dc    = (zr1+zr2+z1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dxr1,dyr1,zr1,dx1,dy1,z1,dxr2,dyr2,zr2)+spriz)/2;
                           da    = (dxr2+dx1+dx4)/3;
                           db    = (dyr2+dy1+dy4)/3;
                           dc    = (zr2+z1+z4)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dxr2,dyr2,zr2,dx1,dy1,z1,dx4,dy4,z4)+spriz)/2;
                           v1    = v/4;
                           if (sr3<0||sr3>w1||tr3<0||tr3>w2) Transform(sr3,tr3,-1);
                           else                              Transform(sr3,tr3,0);
                           dxr1  = fDxspline;
                           dyr1  = fDyspline;
                           zr1   = fZ;
                           if (sr2<0||sr2>w1||tr2<0||tr2>w2) Transform(sr2,tr2,-1);
                           else                              Transform(sr2,tr2,0);
                           dxr2  = fDxspline;
                           dyr2  = fDyspline;
                           zr2   = fZ;
                           da    = (dx1+dx2+dx3)/3;
                           db    = (dy1+dy2+dy3)/3;
                           dc    = (z1+z2+z3)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = (ColorCalculation(dx1,dy1,z1,dx2,dy2,z2,dx3,dy3,z3)+spriz)/2;
                           da    = (dx2+dxr1+dx3)/3;
                           db    = (dy2+dyr1+dy3)/3;
                           dc    = (z2+zr1+z3)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx2,dy2,z2,dxr1,dyr1,zr1,dx3,dy3,z3)+spriz)/2;
                           da    = (dx2+dxr2+dxr1)/3;
                           db    = (dy2+dyr2+dyr1)/3;
                           dc    = (z2+zr2+zr1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dxr1,dyr1,zr1)+spriz)/2;
                           da    = (dxr2+dx2+dx1)/3;
                           db    = (dyr2+dy2+dy1)/3;
                           dc    = (zr2+z2+z1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dxr2,dyr2,zr2,dx2,dy2,z2,dx1,dy1,z1)+spriz)/2;
                           v2    = v/4;
                           if (sr5<0||sr5>w1||tr5<0||tr5>w2) Transform(sr5,tr5,-1);
                           else                              Transform(sr5,tr5,0);
                           dxr1  = fDxspline;
                           dyr1  = fDyspline;
                           zr1   = fZ;
                           if (sr4<0||sr4>w1||tr4<0||tr4>w2) Transform(sr4,tr4,-1);
                           else                              Transform(sr4,tr4,0);
                           dxr2  = fDxspline;
                           dyr2  = fDyspline;
                           zr2   = fZ;
                           da    = (dx2+dx3+dx4)/3;
                           db    = (dy2+dy3+dy4)/3;
                           dc    = (z2+z3+z4)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = (ColorCalculation(dx2,dy2,z2,dx3,dy3,z3,dx4,dy4,z4)+spriz)/2;
                           da    = (dx4+dx3+dxr1)/3;
                           db    = (dy4+dy3+dyr1)/3;
                           dc    = (z4+z3+zr1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr1,dyr1,zr1)+spriz)/2;
                           da    = (dx3+dxr2+dxr1)/3;
                           db    = (dy3+dyr2+dyr1)/3;
                           dc    = (z3+zr2+zr1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx3,dy3,z3,dxr2,dyr2,zr2,dxr1,dyr1,zr1)+spriz)/2;
                           da    = (dx2+dxr2+dx3)/3;
                           db    = (dy2+dyr2+dy3)/3;
                           dc    = (z2+zr2+z3)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dx3,dy3,z3)+spriz)/2;
                           v3    = v/4;
                           if (sr7<0||sr7>w1||tr7<0||tr7>w2) Transform(sr7,tr7,-1);
                           else                              Transform(sr7,tr7,0);
                           dxr1  = fDxspline;
                           dyr1  = fDyspline;
                           zr1   = fZ;
                           if (sr6<0||sr6>w1||tr6<0||tr6>w2) Transform(sr6,tr6,-1);
                           else                              Transform(sr6,tr6,0);
                           dxr2  = fDxspline;
                           dyr2  = fDyspline;
                           zr2   = fZ;
                           da    = (dx1+dx3+dx4)/3;
                           db    = (dy1+dy3+dy4)/3;
                           dc    = (z1+z3+z4)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = (ColorCalculation(dx1,dy1,z1,dx3,dy3,z3,dx4,dy4,z4)+spriz)/2;
                           da    = (dx4+dx3+dxr2)/3;
                           db    = (dy4+dy3+dyr2)/3;
                           dc    = (z4+z3+zr2)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr2,dyr2,zr2)+spriz)/2;
                           da    = (dx4+dxr2+dxr1)/3;
                           db    = (dy4+dyr2+dyr1)/3;
                           dc    = (z4+zr2+zr1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx4,dy4,z4,dxr2,dyr2,zr2,dxr1,dyr1,zr1)+spriz)/2;
                           da    = (dx1+dx4+dxr1)/3;
                           db    = (dy1+dy4+dyr1)/3;
                           dc    = (z1+z4+zr1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx1,dy1,z1,dx4,dy4,z4,dxr1,dyr1,zr1)+spriz)/2;
                           v4    = v/4;
                        }
                     }
                     spriz = 0;
                     if (fShadow==kShadowsNotPainted) {
                        if (fShading==kNotShaded) {
                           v  = v*fLevels+0.5;
                           iv = fLevels-(Int_t)v;
                        } else {
                           v1  = v1*fLevels;
                           iv1 = fLevels-(Int_t)v1;
                           v2  = v2*fLevels;
                           iv2 = fLevels-(Int_t)v2;
                           v4  = v4*fLevels;
                           iv4 = fLevels-(Int_t)v4;
                        }
                     } else {
                        spriz = ShadowColorCalculation(xtaz,ytaz,ztaz,shad_noise);
                        if (fShading==kNotShaded) {
                           v  = v*fLevels/2.0;
                           iv = fLevels-(Int_t)(v+0.5);
                        } else {
                           v1  = v1*fLevels;
                           iv1 = fLevels-(Int_t)v1;
                           v2  = v2*fLevels;
                           iv2 = fLevels-(Int_t)v2;
                           v4  = v4*fLevels;
                           iv4 = fLevels-(Int_t)v4;
                        }
                     }
                     if (fShading==kNotShaded) {
                        ColorModel(iv,ui1,ui2,ui3);
                        line->SetLineColor(fNewColorIndex);
                     } else {
                        dx1 = x1;
                        dy1 = y1;
                        dx2 = x2;
                        dy2 = y2;
                        dx3 = x4;
                        dy3 = y4;
                        z1  = iv1;
                        z2  = iv2;
                        z3  = iv4;
                        da  = (dy2-dy1)*(z3-z1)-(dy3-dy1)*(z2-z1);
                        db  = (z2-z1)*(dx3-dx1)-(z3-z1)*(dx2-dx1);
                        dc  = (dx2-dx1)*(dy3-dy1)-(dx3-dx1)*(dy2-dy1);
                        dd  = -da*dx1-db*dy1-dc*z1;
                     }
                     sx1 = x1;
                     sy1 = y1;
                     sx2 = x2;
                     sy2 = y2;
                     if (sx2<sx1) {
                        sx4 = sx1;
                        sy4 = sy1;
                        sx1 = sx2;
                        sy1 = sy2;
                        sx2 = sx4;
                        sy2 = sy4;
                     }
                     sdx1 = 0;
                     pom1 = sy2-sy1;
                     pom2 = sx2-sx1;
                     if (pom2!=0) sdx1 = pom1/pom2;
                     pom1 = sy1;
                     pom2 = sx1;
                     sdy1 = pom1-sdx1*pom2;
                     for (sx4=sx1,sx5=sx1,sy5=sy1;sx4<=sx2;sx4++) {
                        pom1 = sx4;
                        sdy4 = sdx1*pom1+sdy1;
                        sy4  = (Int_t)sdy4;
                        if (sy4<=fEnvelope[sx4]) {
                           fEnvelope[sx4] = sy4;
                           if (fShading==kNotShaded) {
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           } else {
                              dx1 = sx4;
                              dy1 = sy4;
                              if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                              else       v =(iv1+iv2+iv4)/3;
                              iv = (Int_t)v;
                              ColorModel(iv,ui1,ui2,ui3);
                              line->SetLineColor(fNewColorIndex);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           }
                           sy5 = sy4;
                        } else {
                           sy4 = fEnvelope[sx4];
                           if (fShading==kNotShaded&&sy5<=fEnvelope[sx5]) {
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           } else if (sy5<=fEnvelope[sx5]) {
                              dx1 = sx4;
                              dy1 = sy4;
                              if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                              else       v =(iv1+iv2+iv4)/3;
                              iv = (Int_t)v;
                              ColorModel(iv,ui1,ui2,ui3);
                              line->SetLineColor(fNewColorIndex);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           }
                           sy5 = fEnvelope[sx4];
                        }
                        sx5 = sx4;
                     }
                     xtaz = (dx3+dx2+dx4)/3;
                     ytaz = (dy3+dy2+dy4)/3;
                     ztaz = (z3+z2+z4)/3;
                     if (fShading==kNotShaded) v=ColorCalculation(dx2,dy2,z2,dx3,dy3,z3,dx4,dy4,z4);
                     spriz = 0;
                     if (fShadow==kShadowsNotPainted) {
                        if (fShading==kNotShaded) {
                           v  = v*fLevels;
                           iv = fLevels-(Int_t)v;
                        } else {
                           v3  = v3*fLevels;
                           iv3 = fLevels-(Int_t)v3;
                        }
                     } else {
                        spriz = ShadowColorCalculation(xtaz,ytaz,ztaz,shad_noise);
                        if (fShading==kNotShaded) {
                           v  = v*fLevels/2;
                           iv = fLevels-(Int_t)v;
                           iv = (Int_t)(iv-fLevels*spriz/2);
                        } else {
                           v3  = v3*fLevels;
                           iv3 = fLevels-(Int_t)v3;
                        }
                     }
                     if (fShading==kNotShaded) {
                        ColorModel(iv,ui1,ui2,ui3);
                        line->SetLineColor(fNewColorIndex);
                     } else {
                        dx1 = x2;
                        dy1 = y2;
                        dx2 = x3;
                        dy2 = y3;
                        dx3 = x4;
                        dy3 = y4;
                        z1  = iv2;
                        z2  = iv3;
                        z3  = iv4;
                        da  = (dy2-dy1)*(z3-z1)-(dy3-dy1)*(z2-z1);
                        db  = (z2-z1)*(dx3-dx1)-(z3-z1)*(dx2-dx1);
                        dc  = (dx2-dx1)*(dy3-dy1)-(dx3-dx1)*(dy2-dy1);
                        dd  = -da*dx1-db*dy1-dc*z1;
                     }
                     sx2 = x3;
                     sy2 = y3;
                     sx3 = x4;
                     sy3 = y4;
                     if (sx3<sx2) {
                        sx4 = sx2;
                        sy4 = sy2;
                        sx2 = sx3;
                        sy2 = sy3;
                        sx3 = sx4;
                        sy3 = sy4;
                     }
                     sdx2 = 0;
                     pom1 = sy3-sy2;
                     pom2 = sx3-sx2;
                     if (pom2!=0) sdx2 = pom1/pom2;
                     pom1 = sy2;
                     pom2 = sx2;
                     sdy2 = pom1-sdx2*pom2;
                     for (sx4=sx2,sx5=sx2,sy5=sy2;sx4<=sx3;sx4++) {
                        pom1 = sx4;
                        sdy4 = sdx2*pom1+sdy2;
                        sy4  = (Int_t)sdy4;
                        if (sy4<=fEnvelope[sx4]) {
                           fEnvelope[sx4] = sy4;
                           if (fShading==kNotShaded) {
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           } else {
                              dx1 = sx4;
                              dy1 = sy4;
                              if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                              else       v =(iv1+iv2+iv4)/3;
                              iv = (Int_t)v;
                              ColorModel(iv,ui1,ui2,ui3);
                              line->SetLineColor(fNewColorIndex);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           }
                           sy5 = sy4;
                        } else {
                           sy4 = fEnvelope[sx4];
                           if (fShading==kNotShaded&&sy5<=fEnvelope[sx5]) {
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           } else if (sy5<=fEnvelope[sx5]) {
                              dx1 = sx4;
                              dy1 = sy4;
                              if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                              else       v =(iv1+iv2+iv4)/3;
                              iv = (Int_t)v;
                              ColorModel(iv,ui1,ui2,ui3);
                              line->SetLineColor(fNewColorIndex);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           }
                           sy5 = fEnvelope[sx4];
                        }
                        sx5 = sx4;
                     }
                  }
               }
               break;
            case kDisplayModeNeedles:
               Transform(i,j,-1);
               x2d = fXt;
               y2d = fYt;
               if (flag==1) {
                  x   = x1d;
                  y   = y1d;
                  x1d = x2d;
                  y1d = y2d;
                  x2d = x;
                  y2d = y;
               }
               line->PaintLine(gPad->PixeltoX(x1),gPad->PixeltoY(y1)+1,gPad->PixeltoX(x1d),gPad->PixeltoY(y1d)+1);
               line->PaintLine(gPad->PixeltoX(x2),gPad->PixeltoY(y2)+1,gPad->PixeltoX(x2d),gPad->PixeltoY(y2d)+1);
               break;
            case kDisplayModeSurface:
               box->SetFillStyle(1001);
               if ((q1!=q2||smer!=0)&&flag==1) {
                  s1 = q1+1;
                  t1 = q2;
                  s2 = q1;
                  t2 = q2;
                  s3 = q1;
                  t3 = q2+1;
                  s4 = q1+1;
                  t4 = q2+1;
                  if (fShading==kShaded) {
                     sr1 = s1;
                     tr1 = (Int_t)TMath::Max(t1-1,0);
                     sr2 = s2;
                     tr2 = (Int_t)TMath::Max(t2-1,0);
                     sr3 = (Int_t)TMath::Max(s2-1,0);
                     tr3 = t2;
                     sr4 = (Int_t)TMath::Max(s3-1,0);
                     tr4 = t3;
                     sr5 = s3;
                     tr5 = t3+1;
                     sr6 = s4;
                     tr6 = t4+1;
                     sr7 = s4+1;
                     tr7 = t4;
                     sr8 = s1+1;
                     tr8 = t1;
                  }
                  if (turni==1) {
                     i1 = s1;
                     i2 = s2;
                     i3 = s3;
                     i4 = s4;
                  } else {
                     i1 = (Int_t)TMath::Max(w1-s1,0);
                     i2 = (Int_t)TMath::Max(w1-s2,0);
                     i3 = (Int_t)TMath::Max(w1-s3,0);
                     i4 = (Int_t)TMath::Max(w1-s4,0);
                     if (fShading==kShaded) {
                     sr1 = (Int_t)TMath::Max(w1-sr1,0);
                     sr2 = (Int_t)TMath::Max(w1-sr2,0);
                     sr3 = (Int_t)TMath::Max(w1-sr3,0);
                     sr4 = (Int_t)TMath::Max(w1-sr4,0);
                     sr5 = (Int_t)TMath::Max(w1-sr5,0);
                     sr6 = (Int_t)TMath::Max(w1-sr6,0);
                     sr7 = (Int_t)TMath::Max(w1-sr7,0);
                     sr8 = (Int_t)TMath::Max(w1-sr8,0);
                     }
                  }
                  if (turnj==1) {
                     j1 = t1;
                     j2 = t2;
                     j3 = t3;
                     j4 = t4;
                  } else {
                     j1 = (Int_t)TMath::Max(w2-t1,0);
                     j2 = (Int_t)TMath::Max(w2-t2,0);
                     j3 = (Int_t)TMath::Max(w2-t3,0);
                     j4 = (Int_t)TMath::Max(w2-t4,0);
                     if (fShading==kShaded) {
                        tr1 = (Int_t)TMath::Max(w2-tr1,0);
                        tr2 = (Int_t)TMath::Max(w2-tr2,0);
                        tr3 = (Int_t)TMath::Max(w2-tr3,0);
                        tr4 = (Int_t)TMath::Max(w2-tr4,0);
                        tr5 = (Int_t)TMath::Max(w2-tr5,0);
                        tr6 = (Int_t)TMath::Max(w2-tr6,0);
                        tr7 = (Int_t)TMath::Max(w2-tr7,0);
                        tr8 = (Int_t)TMath::Max(w2-tr8,0);
                     }
                  }
                  Transform(i1,j1,0);
                  x1  = fXt;
                  y1  = fYt;
                  dx1 = fDxspline;
                  dy1 = fDyspline;
                  z1  = fZ;
                  Transform(i2,j2,0);
                  x2  = fXt;
                  y2  = fYt;
                  dx2 = fDxspline;
                  dy2 = fDyspline;
                  z2  = fZ;
                  Transform(i3,j3,0);
                  x3  = fXt;
                  y3  = fYt;
                  dx3 = fDxspline;
                  dy3 = fDyspline;
                  z3  = fZ;
                  Transform(i4,j4,0);
                  x4  = fXt;
                  y4  = fYt;
                  dx4 = fDxspline;
                  dy4 = fDyspline;
                  z4  = fZ;
                  Envelope(x1,y1,x2,y2);
                  Envelope(x2,y2,x3,y3);
                  xtaz = (dx1+dx2+dx4)/3;
                  ytaz = (dy1+dy2+dy4)/3;
                  ztaz = (z1+z2+z4)/3;
                  v    = ColorCalculation(dx1,dy1,z1,dx2,dy2,z2,dx4,dy4,z4);
                  if (fShading==kShaded) {
                     if (fShadow==kShadowsNotPainted) {
                        if (sr1<0||sr1>w1||tr1<0||tr1>w2) Transform(sr1,tr1,-1);
                        else                              Transform(sr1,tr1,0);
                        dxr1 = fDxspline;
                        dyr1 = fDyspline;
                        zr1  = fZ;
                        if (sr8<0||sr8>w1||tr8<0||tr8>w2) Transform(sr8,tr8,-1);
                        else                              Transform(sr8,tr8,0);
                        dxr2 = fDxspline;
                        dyr2 = fDyspline;
                        zr2  = fZ;
                        v    = v+ColorCalculation(dxr1,dyr1,zr1,dx2,dy2,z2,dx1,dy1,z1);
                        v    = v+ColorCalculation(dxr1,dyr1,zr1,dx1,dy1,z1,dxr2,dyr2,zr2);
                        v    = v+ColorCalculation(dxr2,dyr2,zr2,dx1,dy1,z1,dx4,dy4,z4);
                        v1   = v/4;
                        if (sr3<0||sr3>w1||tr3<0||tr3>w2) Transform(sr3,tr3,-1);
                        else                              Transform(sr3,tr3,0);
                        dxr1 = fDxspline;
                        dyr1 = fDyspline;
                        zr1  = fZ;
                        if (sr2<0||sr2>w1||tr2<0||tr2>w2) Transform(sr2,tr2,-1);
                        else                              Transform(sr2,tr2,0);
                        dxr2 = fDxspline;
                        dyr2 = fDyspline;
                        zr2  = fZ;
                        v    = ColorCalculation(dx1,dy1,z1,dx2,dy2,z2,dx3,dy3,z3);
                        v    = v+ColorCalculation(dx2,dy2,z2,dxr1,dyr1,zr1,dx3,dy3,z3);
                        v    = v+ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dxr1,dyr1,zr1);
                        v    = v+ColorCalculation(dxr2,dyr2,zr2,dx2,dy2,z2,dx1,dy1,z1);
                        v2   = v/4;
                        if (sr5<0||sr5>w1||tr5<0||tr5>w2) Transform(sr5,tr5,-1);
                        else                              Transform(sr5,tr5,0);
                        dxr1 = fDxspline;
                        dyr1 = fDyspline;
                        zr1  = fZ;
                        if (sr4<0||sr4>w1||tr4<0||tr4>w2) Transform(sr4,tr4,-1);
                        else                              Transform(sr4,tr4,0);
                        dxr2 = fDxspline;
                        dyr2 = fDyspline;
                        zr2  = fZ;
                        v    = ColorCalculation(dx2,dy2,z2,dx3,dy3,z3,dx4,dy4,z4);
                        v    = v+ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr1,dyr1,zr1);
                        v    = v+ColorCalculation(dx3,dy3,z3,dxr2,dyr2,zr2,dxr1,dyr1,zr1);
                        v    = v+ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dx3,dy3,z3);
                        v3   = v/4;
                        if (sr7<0||sr7>w1||tr7<0||tr7>w2) Transform(sr7,tr7,-1);
                        else                              Transform(sr7,tr7,0);
                        dxr1 = fDxspline;
                        dyr1 = fDyspline;
                        zr1  = fZ;
                        if (sr6<0||sr6>w1||tr6<0||tr6>w2) Transform(sr6,tr6,-1);
                        else                              Transform(sr6,tr6,0);
                        dxr2 = fDxspline;
                        dyr2 = fDyspline;
                        zr2  = fZ;
                        v    = ColorCalculation(dx1,dy1,z1,dx3,dy3,z3,dx4,dy4,z4);
                        v    = v+ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr2,dyr2,zr2);
                        v    = v+ColorCalculation(dx4,dy4,z4,dxr2,dyr2,zr2,dxr1,dyr1,zr1);
                        v    = v+ColorCalculation(dx1,dy1,z1,dx4,dy4,z4,dxr1,dyr1,zr1);
                        v4   = v/4;
                     } else {
                        spriz = ShadowColorCalculation(xtaz,ytaz,ztaz,shad_noise);
                        v     = v+spriz;
                        v     = v/2;
                        if (sr1<0||sr1>w1||tr1<0||tr1>w2) Transform(sr1,tr1,-1);
                        else                              Transform(sr1,tr1,0);
                        dxr1 = fDxspline;
                        dyr1 = fDyspline;
                        zr1  = fZ;
                        if (sr8<0||sr8>w1||tr8<0||tr8>w2) Transform(sr8,tr8,-1);
                        else                              Transform(sr8,tr8,0);
                        dxr2  = fDxspline;
                        dyr2  = fDyspline;
                        zr2   = fZ;
                        da    = (dxr1+dx2+dx1)/3;
                        db    = (dyr1+dy2+dy1)/3;
                        dc    = (zr1+z2+z1)/3;
                        spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                        v     = v+(ColorCalculation(dxr1,dyr1,zr1,dx2,dy2,z2,dx1,dy1,z1)+spriz)/2;
                        da    = (dxr1+dxr2+dx1)/3;
                        db    = (dyr1+dyr2+dy1)/3;
                        dc    = (zr1+zr2+z1)/3;
                        spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                        v     = v+(ColorCalculation(dxr1,dyr1,zr1,dx1,dy1,z1,dxr2,dyr2,zr2)+spriz)/2;
                        da    = (dxr2+dx1+dx4)/3;
                        db    = (dyr2+dy1+dy4)/3;
                        dc    = (zr2+z1+z4)/3;
                        spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                        v     = v+(ColorCalculation(dxr2,dyr2,zr2,dx1,dy1,z1,dx4,dy4,z4)+spriz)/2;
                        v1    = v/4;
                        if (sr3<0||sr3>w1||tr3<0||tr3>w2) Transform(sr3,tr3,-1);
                        else                              Transform(sr3,tr3,0);
                        dxr1  = fDxspline;
                        dyr1  = fDyspline;
                        zr1   = fZ;
                        if (sr2<0||sr2>w1||tr2<0||tr2>w2) Transform(sr2,tr2,-1);
                        else                              Transform(sr2,tr2,0);
                        dxr2  = fDxspline;
                        dyr2  = fDyspline;
                        zr2   = fZ;
                        da    = (dx1+dx2+dx3)/3;
                        db    = (dy1+dy2+dy3)/3;
                        dc    = (z1+z2+z3)/3;
                        spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                        v     = (ColorCalculation(dx1,dy1,z1,dx2,dy2,z2,dx3,dy3,z3)+spriz)/2;
                        da    = (dx2+dxr1+dx3)/3;
                        db    = (dy2+dyr1+dy3)/3;
                        dc    = (z2+zr1+z3)/3;
                        spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                        v     = v+(ColorCalculation(dx2,dy2,z2,dxr1,dyr1,zr1,dx3,dy3,z3)+spriz)/2;
                        da    = (dx2+dxr2+dxr1)/3;
                        db    = (dy2+dyr2+dyr1)/3;
                        dc    = (z2+zr2+zr1)/3;
                        spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                        v     = v+(ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dxr1,dyr1,zr1)+spriz)/2;
                        da    = (dxr2+dx2+dx1)/3;
                        db    = (dyr2+dy2+dy1)/3;
                        dc    = (zr2+z2+z1)/3;
                        spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                        v     = v+(ColorCalculation(dxr2,dyr2,zr2,dx2,dy2,z2,dx1,dy1,z1)+spriz)/2;
                        v2    = v/4;
                        if (sr5<0||sr5>w1||tr5<0||tr5>w2) Transform(sr5,tr5,-1);
                        else                              Transform(sr5,tr5,0);
                        dxr1  = fDxspline;
                        dyr1  = fDyspline;
                        zr1   = fZ;
                        if (sr4<0||sr4>w1||tr4<0||tr4>w2) Transform(sr4,tr4,-1);
                        else                              Transform(sr4,tr4,0);
                        dxr2  = fDxspline;
                        dyr2  = fDyspline;
                        zr2   = fZ;
                        da    = (dx2+dx3+dx4)/3;
                        db    = (dy2+dy3+dy4)/3;
                        dc    = (z2+z3+z4)/3;
                        spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                        v     = (ColorCalculation(dx2,dy2,z2,dx3,dy3,z3,dx4,dy4,z4)+spriz)/2;
                        da    = (dx4+dx3+dxr1)/3;
                        db    = (dy4+dy3+dyr1)/3;
                        dc    = (z4+z3+zr1)/3;
                        spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                        v     = v+(ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr1,dyr1,zr1)+spriz)/2;
                        da    = (dx3+dxr2+dxr1)/3;
                        db    = (dy3+dyr2+dyr1)/3;
                        dc    = (z3+zr2+zr1)/3;
                        spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                        v     = v+(ColorCalculation(dx3,dy3,z3,dxr2,dyr2,zr2,dxr1,dyr1,zr1)+spriz)/2;
                        da    = (dx2+dxr2+dx3)/3;
                        db    = (dy2+dyr2+dy3)/3;
                        dc    = (z2+zr2+z3)/3;
                        spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                        v     = v+(ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dx3,dy3,z3)+spriz)/2;
                        v3    = v/4;
                        if (sr7<0||sr7>w1||tr7<0||tr7>w2) Transform(sr7,tr7,-1);
                        else                              Transform(sr7,tr7,0);
                        dxr1  = fDxspline;
                        dyr1  = fDyspline;
                        zr1   = fZ;
                        if (sr6<0||sr6>w1||tr6<0||tr6>w2) Transform(sr6,tr6,-1);
                        else                              Transform(sr6,tr6,0);
                        dxr2  = fDxspline;
                        dyr2  = fDyspline;
                        zr2   = fZ;
                        da    = (dx1+dx3+dx4)/3;
                        db    = (dy1+dy3+dy4)/3;
                        dc    = (z1+z3+z4)/3;
                        spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                        v     = (ColorCalculation(dx1,dy1,z1,dx3,dy3,z3,dx4,dy4,z4)+spriz)/2;
                        da    = (dx4+dx3+dxr2)/3;
                        db    = (dy4+dy3+dyr2)/3;
                        dc    = (z4+z3+zr2)/3;
                        spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                        v     = v+(ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr2,dyr2,zr2)+spriz)/2;
                        da    = (dx4+dxr2+dxr1)/3;
                        db    = (dy4+dyr2+dyr1)/3;
                        dc    = (z4+zr2+zr1)/3;
                        spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                        v     = v+(ColorCalculation(dx4,dy4,z4,dxr2,dyr2,zr2,dxr1,dyr1,zr1)+spriz)/2;
                        da    = (dx1+dx4+dxr1)/3;
                        db    = (dy1+dy4+dyr1)/3;
                        dc    = (z1+z4+zr1)/3;
                        spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                        v     = v+(ColorCalculation(dx1,dy1,z1,dx4,dy4,z4,dxr1,dyr1,zr1)+spriz)/2;
                        v4    = v/4;
                     }
                  }
                  spriz = 0;
                  if (fShadow==kShadowsNotPainted) {
                     if (fShading==kNotShaded) {
                        v  = v*fLevels+0.5;
                        iv = fLevels-(Int_t)v;
                     } else {
                        v1  = v1*fLevels+0.5;
                        iv1 = fLevels-(Int_t)v1;
                        v2  = v2*fLevels+0.5;
                        iv2 = fLevels-(Int_t)v2;
                        v4  = v4*fLevels+0.5;
                        iv4 = fLevels-(Int_t)v4;
                     }
                  } else {
                     spriz = ShadowColorCalculation(xtaz,ytaz,ztaz,shad_noise);
                     if (fShading==kNotShaded) {
                        v  = v*fLevels/2.0;
                        iv = fLevels-(Int_t)(v+0.5);
                     } else {
                        v1  = v1*fLevels;
                        iv1 = fLevels-(Int_t)v1;
                        v2  = v2*fLevels;
                        iv2 = fLevels-(Int_t)v2;
                        v4  = v4*fLevels;
                        iv4 = fLevels-(Int_t)v4;
                     }
                  }
                  if (fShading==kNotShaded) {
                     ColorModel(iv,ui1,ui2,ui3);
                     box->SetFillColor(fNewColorIndex);
                  } else {
                     dx1 = x1;
                     dy1 = y1;
                     dx2 = x2;
                     dy2 = y2;
                     dx3 = x4;
                     dy3 = y4;
                     z1  = iv1;
                     z2  = iv2;
                     z3  = iv4;
                     da  = (dy2-dy1)*(z3-z1)-(dy3-dy1)*(z2-z1);
                     db  = (z2-z1)*(dx3-dx1)-(z3-z1)*(dx2-dx1);
                     dc  = (dx2-dx1)*(dy3-dy1)-(dx3-dx1)*(dy2-dy1);
                     dd  = -da*dx1-db*dy1-dc*z1;
                  }
                  sx1 = x1;
                  sy1 = y1;
                  sx2 = x2;
                  sy2 = y2;
                  sx3 = x4;
                  sy3 = y4;
                  if (sx2<sx1) {
                     sx4 = sx1;
                     sy4 = sy1;
                     sx1 = sx2;
                     sy1 = sy2;
                     sx2 = sx4;
                     sy2 = sy4;
                  }
                  if (sx3<sx1) {
                     sx4 = sx1;
                     sy4 = sy1;
                     sx1 = sx3;
                     sy1 = sy3;
                     sx3 = sx4;
                     sy3 = sy4;
                  }
                  if (sy2<sy3) {
                     sx4 = sx2;
                     sy4 = sy2;
                     sx2 = sx3;
                     sy2 = sy3;
                     sx3 = sx4;
                     sy3 = sy4;
                  }
                  sdx1 = 0;
                  sdx2 = 0;
                  sdx3 = 0;
                  pom1 = sy2-sy1;
                  pom2 = sx2-sx1;
                  if (pom2!=0) sdx1 = pom1/pom2;
                  pom1 = sy1;
                  pom2 = sx1;
                  sdy1 = pom1-sdx1*pom2;
                  pom1 = sy3-sy1;
                  pom2 = sx3-sx1;
                  if (pom2!=0) sdx2 = pom1/pom2;
                  pom1 = sy1;
                  pom2 = sx1;
                  sdy2 = pom1-sdx2*pom2;
                  pom1 = sy3-sy2;
                  pom2 = sx3-sx2;
                  if (pom2!=0) sdx3 = pom1/pom2;
                  pom1 = sy2;
                  pom2 = sx2;
                  sdy3 = pom1-sdx3*pom2;
                  if (sx2<sx3) {
                     if (sx1!=sx2) {
                        for (sx4=sx1;sx4<=sx2;sx4++) {
                           pom1 = sx4;
                           sdy4 = sdx1*pom1+sdy1;
                           sy4  = (Int_t)sdy4;
                           if (sx3!=sx1) {
                              sdy4 = sdx2*pom1+sdy2;
                              sy5  = (Int_t)sdy4;
                              y5 = fEnvelope[sx4];
                              if (sy4<sy5) {
                                 pom1 = sy4;
                                 sy4  = sy5;
                                 sy5  = (Int_t)pom1;
                              }
                              if ((sy4<=y5)||(sy5<y5)) {
                                 sy4 = (Int_t)TMath::Min(sy4,(Int_t)y5);
                                 sy5 = (Int_t)TMath::Min(sy5,(Int_t)y5);
                                 if (fShading==kNotShaded) {
                                    box->PaintBox(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx4+1),gPad->PixeltoY(sy5-1)+1);
                                 } else {
                                    for (il=sy5;il<=sy4+1;il++) {
                                       dx1 = sx4;
                                       if(il<=sy4) dy1 = il;
                                       else dy1 = sy4;
                                       if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                                       else       v =(iv1+iv2+iv4)/3;
                                       iv = (Int_t)v;
                                       ColorModel(iv,ui1,ui2,ui3);
                                       box->SetFillColor(fNewColorIndex);
                                       box->PaintBox(gPad->PixeltoX(sx4),gPad->PixeltoY(il)+1,gPad->PixeltoX(sx4+1),gPad->PixeltoY(il-1)+1);
                                    }
                                 }
                              }
                           }
                        }
                     }
                     if (sx2!=sx3) {
                        for (sx4=sx2;sx4<=sx3;sx4++) {
                           pom1 = sx4;
                           sdy4 = sdx3*pom1+sdy3;
                           sy4  = (Int_t)sdy4;
                           if (sx3!=sx1) {
                              sdy4 = sdx2*pom1+sdy2;
                              sy5  = (Int_t)sdy4;
                              y5  = fEnvelope[sx4];
                              if (sy4<sy5) {
                                 pom1 = sy4;
                                 sy4  = sy5;
                                 sy5  = (Int_t)pom1;
                              }
                              if ((sy4<=y5)||(sy5<y5)) {
                                 sy4 = (Int_t)TMath::Min(sy4,(Int_t)y5);
                                 sy5 = (Int_t)TMath::Min(sy5,(Int_t)y5);
                                 if (fShading==kNotShaded) {
                                    box->PaintBox(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx4+1),gPad->PixeltoY(sy5-1)+1);
                                 } else {
                                    for (il=sy5;il<=sy4+1;il++) {
                                       dx1 = sx4;
                                       if(il<=sy4) dy1 = il;
                                       else dy1 = sy4;
                                       if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                                       else       v =(iv1+iv2+iv4)/3;
                                       iv = (Int_t)v;
                                       ColorModel(iv,ui1,ui2,ui3);
                                       box->SetFillColor(fNewColorIndex);
                                       box->PaintBox(gPad->PixeltoX(sx4),gPad->PixeltoY(il)+1,gPad->PixeltoX(sx4+1),gPad->PixeltoY(il-1)+1);
                                    }
                                 }
                              }
                           }
                        }
                     }
                  } else {
                     if (sx3!=sx1) {
                        for (sx4=sx1;sx4<=sx3;sx4++) {
                           pom1 = sx4;
                           sdy4 = sdx2*pom1+sdy2;
                           sy4  = (Int_t)sdy4;
                           if (sx2!=sx1) {
                              sdy4 = sdx1*pom1+sdy1;
                              sy5  = (Int_t)sdy4;
                              y5   = fEnvelope[sx4];
                              if (sy4<sy5) {
                                 pom1 = sy4;
                                 sy4  = sy5;
                                 sy5  = (Int_t)pom1;
                              }
                              if ((sy4<=y5)||(sy5<y5)) {
                                 sy4 = (Int_t)TMath::Min(sy4,(Int_t)y5);
                                 sy5 = (Int_t)TMath::Min(sy5,(Int_t)y5);
                                 if (fShading==kNotShaded) {
                                    box->PaintBox(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx4+1),gPad->PixeltoY(sy5-1)+1);
                                 } else {
                                    for (il=sy5;il<=sy4+1;il++) {
                                       dx1 = sx4;
                                       if(il<=sy4) dy1 = il;
                                       else dy1 = sy4;
                                       if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                                       else       v =(iv1+iv2+iv4)/3;
                                       iv = (Int_t)v;
                                       ColorModel(iv,ui1,ui2,ui3);
                                       box->SetFillColor(fNewColorIndex);
                                       box->PaintBox(gPad->PixeltoX(sx4),gPad->PixeltoY(il)+1,gPad->PixeltoX(sx4+1),gPad->PixeltoY(il-1)+1);
                                    }
                                 }
                              }
                           }
                        }
                     }
                     if (sx2!=sx3) {
                        for (sx4=sx3;sx4<=sx2;sx4++) {
                           pom1 = sx4;
                           sdy4 = sdx3*pom1+sdy3;
                           sy4  = (Int_t)sdy4;
                           if (sx2!=sx1) {
                              sdy4 = sdx1*pom1+sdy1;
                              sy5  = (Int_t)sdy4;
                              y5 = fEnvelope[sx4];
                              if (sy4<sy5) {
                                 pom1 = sy4;
                                 sy4  = sy5;
                                 sy5  = (Int_t)pom1;
                              }
                              if ((sy4<=y5)||(sy5<y5)) {
                                 sy4 = (Int_t)TMath::Min(sy4,(Int_t)y5);
                                 sy5 = (Int_t)TMath::Min(sy5,(Int_t)y5);
                                 if (fShading==kNotShaded) {
                                    box->PaintBox(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx4+1),gPad->PixeltoY(sy5-1)+1);
                                 } else {
                                    for (il=sy5;il<=sy4+1;il++) {
                                       dx1 = sx4;
                                       if(il<=sy4) dy1 = il;
                                       else dy1 = sy4;
                                       if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                                       else       v =(iv1+iv2+iv4)/3;
                                       iv = (Int_t)v;
                                       ColorModel(iv,ui1,ui2,ui3);
                                       box->SetFillColor(fNewColorIndex);
                                       box->PaintBox(gPad->PixeltoX(sx4),gPad->PixeltoY(il)+1,gPad->PixeltoX(sx4+1),gPad->PixeltoY(il-1)+1);
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
                  xtaz = (dx3+dx2+dx4)/3;
                  ytaz = (dy3+dy2+dy4)/3;
                  ztaz = (z3+z2+z4)/3;
                  if (fShading==kNotShaded) v=ColorCalculation(dx2,dy2,z2,dx3,dy3,z3,dx4,dy4,z4);
                  spriz = 0;
                  if (fShadow==kShadowsNotPainted) {
                     if (fShading==kNotShaded) {
                        v  = v*fLevels;
                        iv = fLevels-(Int_t)v;
                     } else {
                        v3  = v3*fLevels;
                        iv3 = fLevels-(Int_t)v3;
                     }
                  } else {
                     spriz = ShadowColorCalculation(xtaz,ytaz,ztaz,shad_noise);
                     if (fShading==kNotShaded) {
                        v  = v*fLevels/2;
                        iv = fLevels-(Int_t)v;
                        iv = (Int_t)(iv-fLevels*spriz/2);
                     } else {
                        v3  = v3*fLevels;
                        iv3 = fLevels-(Int_t)v3;
                     }
                  }
                  if (fShading==kNotShaded) {
                     ColorModel(iv,ui1,ui2,ui3);
                     box->SetFillColor(fNewColorIndex);
                  } else {
                     dx1 = x2;
                     dy1 = y2;
                     dx2 = x3;
                     dy2 = y3;
                     dx3 = x4;
                     dy3 = y4;
                     z1  = iv2;
                     z2  = iv3;
                     z3  = iv4;
                     da  = (dy2-dy1)*(z3-z1)-(dy3-dy1)*(z2-z1);
                     db  = (z2-z1)*(dx3-dx1)-(z3-z1)*(dx2-dx1);
                     dc  = (dx2-dx1)*(dy3-dy1)-(dx3-dx1)*(dy2-dy1);
                     dd  = -da*dx1-db*dy1-dc*z1;
                  }
                  sx1 = x2;
                  sy1 = y2;
                  sx2 = x3;
                  sy2 = y3;
                  sx3 = x4;
                  sy3 = y4;
                  if (sx2<sx1) {
                     sx4 = sx1;
                     sy4 = sy1;
                     sx1 = sx2;
                     sy1 = sy2;
                     sx2 = sx4;
                     sy2 = sy4;
                  }
                  if (sx3<sx1) {
                     sx4 = sx1;
                     sy4 = sy1;
                     sx1 = sx3;
                     sy1 = sy3;
                     sx3 = sx4;
                     sy3 = sy4;
                  }
                  if (sy2<sy3) {
                     sx4 = sx2;
                     sy4 = sy2;
                     sx2 = sx3;
                     sy2 = sy3;
                     sx3 = sx4;
                     sy3 = sy4;
                  }
                  pom1 = sy2-sy1;
                  pom2 = sx2-sx1;
                  sdx1 = 0;
                  sdx2 = 0;
                  sdx3 = 0;
                  if (pom2!=0) sdx1 = pom1/pom2;
                  pom1 = sy1;
                  pom2 = sx1;
                  sdy1 = pom1-sdx1*pom2;
                  pom1 = sy3-sy1;
                  pom2 = sx3-sx1;
                  if (pom2!=0) sdx2 = pom1/pom2;
                  pom1 = sy1;
                  pom2 = sx1;
                  sdy2 = pom1-sdx2*pom2;
                  pom1 = sy3-sy2;
                  pom2 = sx3-sx2;
                  if (pom2!=0) sdx3 = pom1/pom2;
                  pom1 = sy2;
                  pom2 = sx2;
                  sdy3 = pom1-sdx3*pom2;
                  if (sx2<sx3) {
                     if (sx1!=sx2) {
                        for (sx4=sx1;sx4<=sx2;sx4++) {
                           pom1 = sx4;
                           sdy4 = sdx1*pom1+sdy1;
                           sy4  = (Int_t)sdy4;
                           if (sx3!=sx1) {
                              sdy4 = sdx2*pom1+sdy2;
                              sy5  = (Int_t)sdy4;
                              y5 = fEnvelope[sx4];
                              if (sy4<sy5) {
                                 pom1 = sy4;
                                 sy4  = sy5;
                                 sy5  = (Int_t)pom1;
                              }
                              if ((sy4<=y5)||(sy5<y5)) {
                                 sy4 = (Int_t)TMath::Min(sy4,(Int_t)y5);
                                 sy5 = (Int_t)TMath::Min(sy5,(Int_t)y5);
                                 if (fShading==kNotShaded) {
                                    box->PaintBox(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx4+1),gPad->PixeltoY(sy5-1)+1);
                                 } else {
                                    for (il=sy5;il<=sy4+1;il++) {
                                       dx1 = sx4;
                                       if(il<=sy4) dy1 = il;
                                       else dy1 = sy4;
                                       if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                                       else       v =(iv2+iv3+iv4)/3;
                                       iv = (Int_t)v;
                                       ColorModel(iv,ui1,ui2,ui3);
                                       box->SetFillColor(fNewColorIndex);
                                       box->PaintBox(gPad->PixeltoX(sx4),gPad->PixeltoY(il)+1,gPad->PixeltoX(sx4+1),gPad->PixeltoY(il-1)+1);
                                    }
                                 }
                              }
                           }
                        }
                     }
                     if (sx2!=sx3) {
                        for (sx4=sx2;sx4<=sx3;sx4++) {
                           pom1 = sx4;
                           sdy4 = sdx3*pom1+sdy3;
                           sy4  = (Int_t)sdy4;
                           if (sx3!=sx1) {
                              sdy4 = sdx2*pom1+sdy2;
                              sy5  = (Int_t)sdy4;
                              y5 = fEnvelope[sx4];
                              if (sy4<sy5) {
                                 pom1 = sy4;
                                 sy4  = sy5;
                                 sy5  = (Int_t)pom1;
                              }
                              if ((sy4<=y5)||(sy5<y5)) {
                                 sy4 = (Int_t)TMath::Min(sy4,(Int_t)y5);
                                 sy5 = (Int_t)TMath::Min(sy5,(Int_t)y5);
                                 if (fShading==kNotShaded) {
                                    box->PaintBox(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx4+1),gPad->PixeltoY(sy5-1)+1);
                                 } else {
                                    for (il=sy5;il<=sy4+1;il++) {
                                       dx1 = sx4;
                                       if(il<=sy4) dy1 = il;
                                       else dy1 = sy4;
                                       if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                                       else       v =(iv2+iv3+iv4)/3;
                                       iv = (Int_t)v;
                                       ColorModel(iv,ui1,ui2,ui3);
                                       box->SetFillColor(fNewColorIndex);
                                       box->PaintBox(gPad->PixeltoX(sx4),gPad->PixeltoY(il)+1,gPad->PixeltoX(sx4+1),gPad->PixeltoY(il-1)+1);
                                    }
                                 }
                              }
                           }
                        }
                     }
                  } else {
                     if (sx3!=sx1) {
                        for (sx4=sx1;sx4<=sx3;sx4++) {
                           pom1 = sx4;
                           sdy4 = sdx2*pom1+sdy2;
                           sy4  = (Int_t)sdy4;
                           if (sx2!=sx1) {
                              sdy4 = sdx1*pom1+sdy1;
                              sy5  = (Int_t)sdy4;
                              y5   = fEnvelope[sx4];
                              if (sy4<sy5) {
                                 pom1 = sy4;
                                 sy4  = sy5;
                                 sy5  = (Int_t)pom1;
                              }
                              if ((sy4<=y5)||(sy5<y5)) {
                                 sy4 = (Int_t)TMath::Min(sy4,(Int_t)y5);
                                 sy5 = (Int_t)TMath::Min(sy5,(Int_t)y5);
                                 if (fShading==kNotShaded) {
                                    box->PaintBox(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx4+1),gPad->PixeltoY(sy5-1)+1);
                                 } else {
                                    for (il=sy5;il<=sy4+1;il++) {
                                       dx1 = sx4;
                                       if(il<=sy4) dy1 = il;
                                       else dy1 = sy4;
                                       if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                                       else       v =(iv2+iv3+iv4)/3;
                                       iv = (Int_t)v;
                                       ColorModel(iv,ui1,ui2,ui3);
                                       box->SetFillColor(fNewColorIndex);
                                       box->PaintBox(gPad->PixeltoX(sx4),gPad->PixeltoY(il)+1,gPad->PixeltoX(sx4+1),gPad->PixeltoY(il-1)+1);
                                    }
                                 }
                              }
                           }
                        }
                     }
                     if (sx2!=sx3) {
                        for (sx4=sx3;sx4<=sx2;sx4++) {
                           pom1 = sx4;
                           sdy4 = sdx3*pom1+sdy3;
                           sy4  = (Int_t)sdy4;
                           if (sx2!=sx1) {
                              sdy4 = sdx1*pom1+sdy1;
                              sy5  = (Int_t)sdy4;
                              y5 = fEnvelope[sx4];
                              if (sy4<sy5) {
                                 pom1 = sy4;
                                 sy4  = sy5;
                                 sy5  = (Int_t)pom1;
                              }
                              if ((sy4<=y5)||(sy5<y5)) {
                                 sy4 = (Int_t)TMath::Min(sy4,(Int_t)y5);
                                 sy5 = (Int_t)TMath::Min(sy5,(Int_t)y5);
                                 if (fShading==kNotShaded) {
                                    box->PaintBox(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx4+1),gPad->PixeltoY(sy5-1)+1);
                                 } else {
                                    for (il=sy5;il<=sy4+1;il++) {
                                       dx1 = sx4;
                                       if(il<=sy4) dy1 = il;
                                       else dy1 = sy4;
                                       if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                                       else       v =(iv2+iv3+iv4)/3;
                                       iv = (Int_t)v;
                                       ColorModel(iv,ui1,ui2,ui3);
                                       box->SetFillColor(fNewColorIndex);
                                       box->PaintBox(gPad->PixeltoX(sx4),gPad->PixeltoY(il)+1,gPad->PixeltoX(sx4+1),gPad->PixeltoY(il-1)+1);
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
               break;
            case kDisplayModeTriangles:
               if (fModeGroup==kModeGroupSimple) {
                  if ((q1!=q2||smer!=0)&&flag==1) {
                     s1 = q1+1;
                     t1 = q2;
                     s2 = q1;
                     t2 = q2;
                     s3 = q1;
                     t3 = q2+1;
                     s4 = q1+1;
                     t4 = q2+1;
                  }
                  if (turni==1) {
                     i1 = (Int_t)TMath::Min(w1,s1);
                     i2 = (Int_t)TMath::Min(w1,s2);
                     i3 = (Int_t)TMath::Min(w1,s3);
                     i4 = (Int_t)TMath::Min(w1,s4);
                  } else {
                     i1 = (Int_t)TMath::Max(w1-s1,0);
                     i2 = (Int_t)TMath::Max(w1-s2,0);
                     i3 = (Int_t)TMath::Max(w1-s3,0);
                     i4 = (Int_t)TMath::Max(w1-s4,0);
                  }
                  if (turnj==1) {
                     j1 = (Int_t)TMath::Min(w2,t1);
                     j2 = (Int_t)TMath::Min(w2,t2);
                     j3 = (Int_t)TMath::Min(w2,t3);
                     j4 = (Int_t)TMath::Min(w2,t4);
                  } else {
                     j1 = (Int_t)TMath::Max(w2-t1,0);
                     j2 = (Int_t)TMath::Max(w2-t2,0);
                     j3 = (Int_t)TMath::Max(w2-t3,0);
                     j4 = (Int_t)TMath::Max(w2-t4,0);
                  }
                  Transform(i1,j1,0);
                  x1 = fXt;
                  y1 = fYt;
                  Transform(i2,j2,0);
                  x2 = fXt;
                  y2 = fYt;
                  Transform(i3,j3,0);
                  x3 = fXt;
                  y3 = fYt;
                  Transform(i4,j4,0);
                  x4 = fXt;
                  y4 = fYt;
                  if ((q1!=q2||smer!=0)&&flag==1) {
                     Envelope(x1,y1,x2,y2);
                     if (fLine!=0) {
                        if (fLine==1) {
                           fXe = x2;
                           fYe = y2;
                        }
                        line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                     }
                     Envelope(x2,y2,x3,y3);
                     if (fLine!=0) {
                        if (fLine==1) {
                           fXe = x3;
                           fYe = y3;
                        }
                        line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                     }
                     Envelope(x2,y2,x4,y4);
                     if (fLine!=0) {
                        if (fLine==1) {
                           fXe = x4;
                           fYe = y4;
                        }
                        line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                     }
                     Envelope(x1,y1,x4,y4);
                     if (fLine!=0) {
                        if (fLine==1) {
                           fXe = x4;
                           fYe = y4;
                        }
                        line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                     }
                     Envelope(x3,y3,x4,y4);
                     if (fLine!=0) {
                        if (fLine==1) {
                           fXe = x4;
                           fYe = y4;
                        }
                        line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
                     }
                  }
               } else {
                  if ((q1!=q2||smer!=0)&&flag==1) {
                     s1 = q1+1;
                     t1 = q2;
                     s2 = q1;
                     t2 = q2;
                     s3 = q1;
                     t3 = q2+1;
                     s4 = q1+1;
                     t4 = q2+1;
                     if (fShading==kShaded) {
                        sr1 = s1;
                        tr1 = (Int_t)TMath::Max(t1-1,0);
                        sr2 = s2;
                        tr2 = (Int_t)TMath::Max(t2-1,0);
                        sr3 = (Int_t)TMath::Max(s2-1,0);
                        tr3 = t2;
                        sr4 = (Int_t)TMath::Max(s3-1,0);
                        tr4 = t3;
                        sr5 = s3;
                        tr5 = t3+1;
                        sr6 = s4;
                        tr6 = t4+1;
                        sr7 = s4+1;
                        tr7 = t4;
                        sr8 = s1+1;
                        tr8 = t1;
                     }
                     if (turni==1) {
                        i1 = s1;
                        i2 = s2;
                        i3 = s3;
                        i4 = s4;
                     } else {
                        i1 = (Int_t)TMath::Max(w1-s1,0);
                        i2 = (Int_t)TMath::Max(w1-s2,0);
                        i3 = (Int_t)TMath::Max(w1-s3,0);
                        i4 = (Int_t)TMath::Max(w1-s4,0);
                        if (fShading==kShaded) {
                           sr1 = (Int_t)TMath::Max(w1-sr1,0);
                           sr2 = (Int_t)TMath::Max(w1-sr2,0);
                           sr3 = (Int_t)TMath::Max(w1-sr3,0);
                           sr4 = (Int_t)TMath::Max(w1-sr4,0);
                           sr5 = (Int_t)TMath::Max(w1-sr5,0);
                           sr6 = (Int_t)TMath::Max(w1-sr6,0);
                           sr7 = (Int_t)TMath::Max(w1-sr7,0);
                           sr8 = (Int_t)TMath::Max(w1-sr8,0);
                        }
                     }
                     if (turnj==1) {
                        j1 = t1;
                        j2 = t2;
                        j3 = t3;
                        j4 = t4;
                     } else {
                        j1 = (Int_t)TMath::Max(w2-t1,0);
                        j2 = (Int_t)TMath::Max(w2-t2,0);
                        j3 = (Int_t)TMath::Max(w2-t3,0);
                        j4 = (Int_t)TMath::Max(w2-t4,0);
                        if (fShading==kShaded) {
                           tr1 = (Int_t)TMath::Max(w2-tr1,0);
                           tr2 = (Int_t)TMath::Max(w2-tr2,0);
                           tr3 = (Int_t)TMath::Max(w2-tr3,0);
                           tr4 = (Int_t)TMath::Max(w2-tr4,0);
                           tr5 = (Int_t)TMath::Max(w2-tr5,0);
                           tr6 = (Int_t)TMath::Max(w2-tr6,0);
                           tr7 = (Int_t)TMath::Max(w2-tr7,0);
                           tr8 = (Int_t)TMath::Max(w2-tr8,0);
                        }
                     }
                     Transform(i1,j1,0);
                     x1  = fXt;
                     y1  = fYt;
                     dx1 = fDxspline;
                     dy1 = fDyspline;
                     z1  = fZ;
                     Transform(i2,j2,0);
                     x2  = fXt;
                     y2  = fYt;
                     dx2 = fDxspline;
                     dy2 = fDyspline;
                     z2  = fZ;
                     Transform(i3,j3,0);
                     x3  = fXt;
                     y3  = fYt;
                     dx3 = fDxspline;
                     dy3 = fDyspline;
                     z3  = fZ;
                     Transform(i4,j4,0);
                     x4  = fXt;
                     y4  = fYt;
                     dx4 = fDxspline;
                     dy4 = fDyspline;
                     z4  = fZ;
                     Envelope(x1,y1,x2,y2);
                     Envelope(x2,y2,x3,y3);
                     xtaz = (dx1+dx2+dx4)/3;
                     ytaz = (dy1+dy2+dy4)/3;
                     ztaz = (z1+z2+z4)/3;
                     v    = ColorCalculation(dx1,dy1,z1,dx2,dy2,z2,dx4,dy4,z4);
                     if (fShading==kShaded) {
                        if (fShadow==kShadowsNotPainted) {
                           if (sr1<0||sr1>w1||tr1<0||tr1>w2) Transform(sr1,tr1,-1);
                           else                              Transform(sr1,tr1,0);
                           dxr1 = fDxspline;
                           dyr1 = fDyspline;
                           zr1  = fZ;
                           if (sr8<0||sr8>w1||tr8<0||tr8>w2) Transform(sr8,tr8,-1);
                           else                              Transform(sr8,tr8,0);
                           dxr2 = fDxspline;
                           dyr2 = fDyspline;
                           zr2  = fZ;
                           v    = v+ColorCalculation(dxr1,dyr1,zr1,dx2,dy2,z2,dx1,dy1,z1);
                           v    = v+ColorCalculation(dxr1,dyr1,zr1,dx1,dy1,z1,dxr2,dyr2,zr2);
                           v    = v+ColorCalculation(dxr2,dyr2,zr2,dx1,dy1,z1,dx4,dy4,z4);
                           v1   = v/4;
                           if (sr3<0||sr3>w1||tr3<0||tr3>w2) Transform(sr3,tr3,-1);
                           else                              Transform(sr3,tr3,0);
                           dxr1 = fDxspline;
                           dyr1 = fDyspline;
                           zr1  = fZ;
                           if (sr2<0||sr2>w1||tr2<0||tr2>w2) Transform(sr2,tr2,-1);
                           else                              Transform(sr2,tr2,0);
                           dxr2 = fDxspline;
                           dyr2 = fDyspline;
                           zr2  = fZ;
                           v    = ColorCalculation(dx1,dy1,z1,dx2,dy2,z2,dx3,dy3,z3);
                           v    = v+ColorCalculation(dx2,dy2,z2,dxr1,dyr1,zr1,dx3,dy3,z3);
                           v    = v+ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dxr1,dyr1,zr1);
                           v    = v+ColorCalculation(dxr2,dyr2,zr2,dx2,dy2,z2,dx1,dy1,z1);
                           v2   = v/4;
                           if (sr5<0||sr5>w1||tr5<0||tr5>w2) Transform(sr5,tr5,-1);
                           else                              Transform(sr5,tr5,0);
                           dxr1 = fDxspline;
                           dyr1 = fDyspline;
                           zr1  = fZ;
                           if (sr4<0||sr4>w1||tr4<0||tr4>w2) Transform(sr4,tr4,-1);
                           else                              Transform(sr4,tr4,0);
                           dxr2 = fDxspline;
                           dyr2 = fDyspline;
                           zr2  = fZ;
                           v    = ColorCalculation(dx2,dy2,z2,dx3,dy3,z3,dx4,dy4,z4);
                           v    = v+ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr1,dyr1,zr1);
                           v    = v+ColorCalculation(dx3,dy3,z3,dxr2,dyr2,zr2,dxr1,dyr1,zr1);
                           v    = v+ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dx3,dy3,z3);
                           v3   = v/4;
                           if (sr7<0||sr7>w1||tr7<0||tr7>w2) Transform(sr7,tr7,-1);
                           else                              Transform(sr7,tr7,0);
                           dxr1 = fDxspline;
                           dyr1 = fDyspline;
                           zr1  = fZ;
                           if (sr6<0||sr6>w1||tr6<0||tr6>w2) Transform(sr6,tr6,-1);
                           else                              Transform(sr6,tr6,0);
                           dxr2 = fDxspline;
                           dyr2 = fDyspline;
                           zr2  = fZ;
                           v    = ColorCalculation(dx1,dy1,z1,dx3,dy3,z3,dx4,dy4,z4);
                           v    = v+ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr2,dyr2,zr2);
                           v    = v+ColorCalculation(dx4,dy4,z4,dxr2,dyr2,zr2,dxr1,dyr1,zr1);
                           v    = v+ColorCalculation(dx1,dy1,z1,dx4,dy4,z4,dxr1,dyr1,zr1);
                           v4   = v/4;
                        } else {
                           spriz = ShadowColorCalculation(xtaz,ytaz,ztaz,shad_noise);
                           v     = v+spriz;
                           v     = v/2;
                           if (sr1<0||sr1>w1||tr1<0||tr1>w2) Transform(sr1,tr1,-1);
                           else                              Transform(sr1,tr1,0);
                           dxr1  = fDxspline;
                           dyr1  = fDyspline;
                           zr1   = fZ;
                           if (sr8<0||sr8>w1||tr8<0||tr8>w2) Transform(sr8,tr8,-1);
                           else                              Transform(sr8,tr8,0);
                           dxr2  = fDxspline;
                           dyr2  = fDyspline;
                           zr2   = fZ;
                           da    = (dxr1+dx2+dx1)/3;
                           db    = (dyr1+dy2+dy1)/3;
                           dc    = (zr1+z2+z1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dxr1,dyr1,zr1,dx2,dy2,z2,dx1,dy1,z1)+spriz)/2;
                           da    = (dxr1+dxr2+dx1)/3;
                           db    = (dyr1+dyr2+dy1)/3;
                           dc    = (zr1+zr2+z1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dxr1,dyr1,zr1,dx1,dy1,z1,dxr2,dyr2,zr2)+spriz)/2;
                           da    = (dxr2+dx1+dx4)/3;
                           db    = (dyr2+dy1+dy4)/3;
                           dc    = (zr2+z1+z4)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dxr2,dyr2,zr2,dx1,dy1,z1,dx4,dy4,z4)+spriz)/2;
                           v1    = v/4;
                           if (sr3<0||sr3>w1||tr3<0||tr3>w2) Transform(sr3,tr3,-1);
                           else                              Transform(sr3,tr3,0);
                           dxr1  = fDxspline;
                           dyr1  = fDyspline;
                           zr1   = fZ;
                           if (sr2<0||sr2>w1||tr2<0||tr2>w2) Transform(sr2,tr2,-1);
                           else                              Transform(sr2,tr2,0);
                           dxr2  = fDxspline;
                           dyr2  = fDyspline;
                           zr2   = fZ;
                           da    = (dx1+dx2+dx3)/3;
                           db    = (dy1+dy2+dy3)/3;
                           dc    = (z1+z2+z3)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = (ColorCalculation(dx1,dy1,z1,dx2,dy2,z2,dx3,dy3,z3)+spriz)/2;
                           da    = (dx2+dxr1+dx3)/3;
                           db    = (dy2+dyr1+dy3)/3;
                           dc    = (z2+zr1+z3)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx2,dy2,z2,dxr1,dyr1,zr1,dx3,dy3,z3)+spriz)/2;
                           da    = (dx2+dxr2+dxr1)/3;
                           db    = (dy2+dyr2+dyr1)/3;
                           dc    = (z2+zr2+zr1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dxr1,dyr1,zr1)+spriz)/2;
                           da    = (dxr2+dx2+dx1)/3;
                           db    = (dyr2+dy2+dy1)/3;
                           dc    = (zr2+z2+z1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dxr2,dyr2,zr2,dx2,dy2,z2,dx1,dy1,z1)+spriz)/2;
                           v2    = v/4;
                           if (sr5<0||sr5>w1||tr5<0||tr5>w2) Transform(sr5,tr5,-1);
                           else                              Transform(sr5,tr5,0);
                           dxr1  = fDxspline;
                           dyr1  = fDyspline;
                           zr1   = fZ;
                           if (sr4<0||sr4>w1||tr4<0||tr4>w2) Transform(sr4,tr4,-1);
                           else                              Transform(sr4,tr4,0);
                           dxr2  = fDxspline;
                           dyr2  = fDyspline;
                           zr2   = fZ;
                           da    = (dx2+dx3+dx4)/3;
                           db    = (dy2+dy3+dy4)/3;
                           dc    = (z2+z3+z4)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = (ColorCalculation(dx2,dy2,z2,dx3,dy3,z3,dx4,dy4,z4)+spriz)/2;
                           da    = (dx4+dx3+dxr1)/3;
                           db    = (dy4+dy3+dyr1)/3;
                           dc    = (z4+z3+zr1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr1,dyr1,zr1)+spriz)/2;
                           da    = (dx3+dxr2+dxr1)/3;
                           db    = (dy3+dyr2+dyr1)/3;
                           dc    = (z3+zr2+zr1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx3,dy3,z3,dxr2,dyr2,zr2,dxr1,dyr1,zr1)+spriz)/2;
                           da    = (dx2+dxr2+dx3)/3;
                           db    = (dy2+dyr2+dy3)/3;
                           dc    = (z2+zr2+z3)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx2,dy2,z2,dxr2,dyr2,zr2,dx3,dy3,z3)+spriz)/2;
                           v3    = v/4;
                           if (sr7<0||sr7>w1||tr7<0||tr7>w2) Transform(sr7,tr7,-1);
                           else                              Transform(sr7,tr7,0);
                           dxr1  = fDxspline;
                           dyr1  = fDyspline;
                           zr1   = fZ;
                           if (sr6<0||sr6>w1||tr6<0||tr6>w2) Transform(sr6,tr6,-1);
                           else                              Transform(sr6,tr6,0);
                           dxr2  = fDxspline;
                           dyr2  = fDyspline;
                           zr2   = fZ;
                           da    = (dx1+dx3+dx4)/3;
                           db    = (dy1+dy3+dy4)/3;
                           dc    = (z1+z3+z4)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = (ColorCalculation(dx1,dy1,z1,dx3,dy3,z3,dx4,dy4,z4)+spriz)/2;
                           da    = (dx4+dx3+dxr2)/3;
                           db    = (dy4+dy3+dyr2)/3;
                           dc    = (z4+z3+zr2)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx4,dy4,z4,dx3,dy3,z3,dxr2,dyr2,zr2)+spriz)/2;
                           da    = (dx4+dxr2+dxr1)/3;
                           db    = (dy4+dyr2+dyr1)/3;
                           dc    = (z4+zr2+zr1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx4,dy4,z4,dxr2,dyr2,zr2,dxr1,dyr1,zr1)+spriz)/2;
                           da    = (dx1+dx4+dxr1)/3;
                           db    = (dy1+dy4+dyr1)/3;
                           dc    = (z1+z4+zr1)/3;
                           spriz = ShadowColorCalculation(da,db,dc,shad_noise);
                           v     = v+(ColorCalculation(dx1,dy1,z1,dx4,dy4,z4,dxr1,dyr1,zr1)+spriz)/2;
                           v4    = v/4;
                        }
                     }
                     spriz = 0;
                     if (fShadow==kShadowsNotPainted) {
                        if (fShading==kNotShaded) {
                           v  = v*fLevels+0.5;
                           iv = fLevels-(Int_t)v;
                        } else {
                           v1  = v1*fLevels;
                           iv1 = fLevels-(Int_t)v1;
                           v2  = v2*fLevels;
                           iv2 = fLevels-(Int_t)v2;
                           v4  = v4*fLevels;
                           iv4 = fLevels-(Int_t)v4;
                        }
                     } else {
                        spriz = ShadowColorCalculation(xtaz,ytaz,ztaz,shad_noise);
                        if (fShading==kNotShaded) {
                           v  = v*fLevels/2.0;
                           iv = fLevels-(Int_t)(v+0.5);
                        } else {
                           v1  = v1*fLevels;
                           iv1 = fLevels-(Int_t)v1;
                           v2  = v2*fLevels;
                           iv2 = fLevels-(Int_t)v2;
                           v4  = v4*fLevels;
                           iv4 = fLevels-(Int_t)v4;
                        }
                     }
                     if (fShading==kNotShaded) {
                        ColorModel(iv,ui1,ui2,ui3);
                        line->SetLineColor(fNewColorIndex);
                     } else {
                        dx1 = x1;
                        dy1 = y1;
                        dx2 = x2;
                        dy2 = y2;
                        dx3 = x4;
                        dy3 = y4;
                        z1  = iv1;
                        z2  = iv2;
                        z3  = iv4;
                        da  = (dy2-dy1)*(z3-z1)-(dy3-dy1)*(z2-z1);
                        db  = (z2-z1)*(dx3-dx1)-(z3-z1)*(dx2-dx1);
                        dc  = (dx2-dx1)*(dy3-dy1)-(dx3-dx1)*(dy2-dy1);
                        dd  = -da*dx1-db*dy1-dc*z1;
                     }
                     sx1 = x1;
                     sy1 = y1;
                     sx2 = x2;
                     sy2 = y2;
                     if (sx2<sx1) {
                        sx4 = sx1;
                        sy4 = sy1;
                        sx1 = sx2;
                        sy1 = sy2;
                        sx2 = sx4;
                        sy2 = sy4;
                     }
                     sdx1 = 0;
                     pom1 = sy2-sy1;
                     pom2 = sx2-sx1;
                     if (pom2!=0) sdx1 = pom1/pom2;
                     pom1 = sy1;
                     pom2 = sx1;
                     sdy1 = pom1-sdx1*pom2;
                     for (sx4=sx1,sx5=sx1,sy5=sy1;sx4<=sx2;sx4++) {
                        pom1 = sx4;
                        sdy4 = sdx1*pom1+sdy1;
                        sy4  = (Int_t)sdy4;
                        if (sy4<=fEnvelope[sx4]) {
                           fEnvelope[sx4] = sy4;
                           if (fShading==kNotShaded) {
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           } else {
                              dx1 = sx4;
                              dy1 = sy4;
                              if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                              else v = (iv1+iv2+iv4)/3;
                              iv = (Int_t)v;
                              ColorModel(iv,ui1,ui2,ui3);
                              line->SetLineColor(fNewColorIndex);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           }
                           sy5 = sy4;
                        } else {
                        sy4 = fEnvelope[sx4];
                           if (fShading==kNotShaded&&sy5<=fEnvelope[sx5]) {
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           } else if (sy5<=fEnvelope[sx5]) {
                              dx1 = sx4;
                              dy1 = sy4;
                              if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                              else v = (iv1+iv2+iv4)/3;
                              iv = (Int_t)v;
                              ColorModel(iv,ui1,ui2,ui3);
                              line->SetLineColor(fNewColorIndex);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           }
                           sy5 = fEnvelope[sx4];
                        }
                        sx5 = sx4;
                     }
                     sx1 = x1;
                     sy1 = y1;
                     sx3 = x4;
                     sy3 = y4;
                     if (sx3<sx1) {
                        sx4 = sx1;
                        sy4 = sy1;
                        sx1 = sx3;
                        sy1 = sy3;
                        sx3 = sx4;
                        sy3 = sy4;
                     }
                     pom1 = sy3-sy1;
                     pom2 = sx3-sx1;
                     if (pom2!=0) sdx2 = pom1/pom2;
                     pom1  = sy1;
                     pom2  = sx1;
                     sdy2  = pom1-sdx2*pom2;
                     sx1p  = sx1;
                     sy1p  = sy1;
                     sx3p  = sx3;
                     sdx2p = sdx2;
                     sdy2p = sdy2;
                     dap   = da;
                     dbp   = db;
                     dcp   = dc;
                     ddp   = dd;
                     uip   = fNewColorIndex;
                     xtaz  = (dx3+dx2+dx4)/3;
                     ytaz  = (dy3+dy2+dy4)/3;
                     ztaz  = (z3+z2+z4)/3;
                     if (fShading==kNotShaded) v = ColorCalculation(dx2,dy2,z2,dx3,dy3,z3,dx4,dy4,z4);
                     spriz = 0;
                     if (fShadow==kShadowsNotPainted) {
                        if (fShading==kNotShaded) {
                           v  =  v*fLevels;
                           iv = fLevels-(Int_t)v;
                        } else {
                           v3 = v3*fLevels;
                           iv3 = fLevels-(Int_t)v3;
                        }
                     } else {
                        spriz = ShadowColorCalculation(xtaz,ytaz,ztaz,shad_noise);
                        if (fShading==kNotShaded) {
                           v  =  v*fLevels/2;
                           iv = fLevels-(Int_t)v;
                           iv = (Int_t)(iv-fLevels*spriz/2);
                        } else {
                           v3  = v3*fLevels;
                           iv3 = fLevels-(Int_t)v3;
                        }
                     }
                     if (fShading==kNotShaded) {
                        ColorModel(iv,ui1,ui2,ui3);
                        line->SetLineColor(fNewColorIndex);
                     } else {
                        dx1 = x2;
                        dy1 = y2;
                        dx2 = x3;
                        dy2 = y3;
                        dx3 = x4;
                        dy3 = y4;
                        z1 = iv2;
                        z2 = iv3;
                        z3 = iv4;
                        da = (dy2-dy1)*(z3-z1)-(dy3-dy1)*(z2-z1);
                        db = (z2-z1)*(dx3-dx1)-(z3-z1)*(dx2-dx1);
                        dc = (dx2-dx1)*(dy3-dy1)-(dx3-dx1)*(dy2-dy1);
                        dd = -da*dx1-db*dy1-dc*z1;
                     }
                     sx1 = x2;
                     sy1 = y2;
                     sx2 = x3;
                     sy2 = y3;
                     if (sx2<sx1) {
                        sx4 = sx1;
                        sy4 = sy1;
                        sx1 = sx2;
                        sy1 = sy2;
                        sx2 = sx4;
                        sy2 = sy4;
                     }
                     pom1 = sy2-sy1;
                     pom2 = sx2-sx1;
                     sdx1 = 0;
                     if (pom2!=0) sdx1 = pom1/pom2;
                     pom1 = sy1;
                     pom2 = sx1;
                     sdy1 = pom1-sdx1*pom2;
                     for (sx4=sx1,sx5=sx1,sy5=sy1;sx4<=sx2;sx4++) {
                        pom1 = sx4;
                        sdy4 = sdx1*pom1+sdy1;
                        sy4  = (Int_t)sdy4;
                        if (sy4<=fEnvelope[sx4]) {
                           fEnvelope[sx4] = sy4;
                           if (fShading==kNotShaded) {
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           } else {
                              dx1 = sx4;
                              dy1 = sy4;
                              if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                              else v = (iv1+iv2+iv4)/3;
                              iv = (Int_t)v;
                              ColorModel(iv,ui1,ui2,ui3);
                              line->SetLineColor(fNewColorIndex);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           }
                           sy5 = sy4;
                        } else {
                           sy4 = fEnvelope[sx4];
                           if (fShading==kNotShaded&&sy5<=fEnvelope[sx5]) {
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           } else if (sy5<=fEnvelope[sx5]) {
                              dx1 = sx4;
                              dy1 = sy4;
                              if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                              else v = (iv1+iv2+iv4)/3;
                              iv = (Int_t)v;
                              ColorModel(iv,ui1,ui2,ui3);
                              line->SetLineColor(fNewColorIndex);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           }
                           sy5 = fEnvelope[sx4];
                        }
                        sx5 = sx4;
                     }
                     sx1 = x2;
                     sy1 = y2;
                     sx2 = x4;
                     sy2 = y4;
                     if (sx2<sx1) {
                        sx4 = sx1;
                        sy4 = sy1;
                        sx1 = sx2;
                        sy1 = sy2;
                        sx2 = sx4;
                        sy2 = sy4;
                     }
                     pom1 = sy2-sy1;
                     pom2 = sx2-sx1;
                     sdx1 = 0;
                     if (pom2!=0) sdx1 = pom1/pom2;
                     pom1 = sy1;
                     pom2 = sx1;
                     sdy1 = pom1-sdx1*pom2;
                     for (sx4=sx1,sx5=sx1,sy5=sy1;sx4<=sx2;sx4++) {
                        pom1 = sx4;
                        sdy4 = sdx1*pom1+sdy1;
                        sy4  = (Int_t)sdy4;
                        if (sy4<=fEnvelope[sx4]) {
                           fEnvelope[sx4] = sy4;
                           if (fShading==kNotShaded) {
                           line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                        } else {
                           dx1 = sx4;
                           dy1 = sy4;
                           if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                           else v = (iv1+iv2+iv4)/3;
                           iv = (Int_t)v;
                           ColorModel(iv,ui1,ui2,ui3);
                           line->SetLineColor(fNewColorIndex);
                           line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                        }
                        sy5 = sy4;
                        } else {
                           sy4 = fEnvelope[sx4];
                           if (fShading==kNotShaded&&sy5<=fEnvelope[sx5]) {
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           } else if (sy5<=fEnvelope[sx5]) {
                              dx1 = sx4;
                              dy1 = sy4;
                              if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                              else v = (iv1+iv2+iv4)/3;
                              iv = (Int_t)v;
                              ColorModel(iv,ui1,ui2,ui3);
                              line->SetLineColor(fNewColorIndex);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           }
                           sy5 = fEnvelope[sx4];
                        }
                        sx5 = sx4;
                     }
                     for (sx4=sx1p,sx5=sx1p,sy5=sy1p;sx4<=sx3p;sx4++) {
                        pom1 = sx4;
                        sdy4 = sdx2p*pom1+sdy2p;
                        sy4  = (Int_t)sdy4;
                        if (sy4<=fEnvelope[sx4]) {
                           fEnvelope[sx4] = sy4;
                           if (fShading==kNotShaded) {
                              line->SetLineColor(uip);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           } else {
                              dx1 = sx4;
                              dy1 = sy4;
                              if (dcp!=0) v = (-ddp-dap*dx1-dbp*dy1)/dcp;
                              else v = (iv1+iv2+iv4)/3;
                              iv = (Int_t)v;
                              ColorModel(iv,ui1,ui2,ui3);
                              line->SetLineColor(fNewColorIndex);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           }
                           sy5 = sy4;
                        } else {
                           sy4 = fEnvelope[sx4];
                           if (fShading==kNotShaded&&sy5<=fEnvelope[sx5]) {
                              line->SetLineColor(uip);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           } else if (sy5<=fEnvelope[sx5]) {
                              dx1 = sx4;
                              dy1 = sy4;
                              if (dcp!=0) v = (-ddp-dap*dx1-dbp*dy1)/dcp;
                              else v = (iv1+iv2+iv4)/3;
                              iv = (Int_t)v;
                              ColorModel(iv,ui1,ui2,ui3);
                              line->SetLineColor(fNewColorIndex);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           }
                           sy5 = fEnvelope[sx4];
                        }
                        sx5 = sx4;
                     }
                     sx2 = x3;
                     sy2 = y3;
                     sx3 = x4;
                     sy3 = y4;
                     if (sx3<sx2) {
                        sx4 = sx2;
                        sy4 = sy2;
                        sx2 = sx3;
                        sy2 = sy3;
                        sx3 = sx4;
                        sy3 = sy4;
                     }
                     sdx2 = 0;
                     pom1 = sy3-sy2;
                     pom2 = sx3-sx2;
                     if (pom2!=0) sdx2 = pom1/pom2;
                     pom1 = sy2;
                     pom2 = sx2;
                     sdy2 = pom1-sdx2*pom2;
                     for (sx4=sx2,sx5=sx2,sy5=sy2;sx4<=sx3;sx4++) {
                        pom1 = sx4;
                        sdy4 = sdx2*pom1+sdy2;
                        sy4  = (Int_t)sdy4;
                        if (sy4<=fEnvelope[sx4]) {
                           fEnvelope[sx4] = sy4;
                           if (fShading==kNotShaded) {
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           } else {
                              dx1 = sx4;
                              dy1 = sy4;
                              if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                              else v = (iv1+iv2+iv4)/3;
                              iv = (Int_t)v;
                              ColorModel(iv,ui1,ui2,ui3);
                              line->SetLineColor(fNewColorIndex);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           }
                           sy5 = sy4;
                        } else {
                           sy4 = fEnvelope[sx4];
                           if (fShading==kNotShaded&&sy5<=fEnvelope[sx5]) {
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           } else if (sy5<=fEnvelope[sx5]) {
                              dx1 = sx4;
                              dy1 = sy4;
                              if (dc!=0) v = (-dd-da*dx1-db*dy1)/dc;
                              else v = (iv1+iv2+iv4)/3;
                              iv = (Int_t)v;
                              ColorModel(iv,ui1,ui2,ui3);
                              line->SetLineColor(fNewColorIndex);
                              line->PaintLine(gPad->PixeltoX(sx4),gPad->PixeltoY(sy4)+1,gPad->PixeltoX(sx5),gPad->PixeltoY(sy5)+1);
                           }
                           sy5 = fEnvelope[sx4];
                        }
                        sx5 = sx4;
                     }
                  }
               }
               break;
         }
         if (flag==1) {
            x   = x1;
            y   = y1;
            x1  = x2;
            y1  = y2;
            x2  = x;
            y2  = y;
            x   = x1d;
            y   = y1d;
            x1d = x2d;
            y1d = y2d;
            x2d = x;
            y2d = y;
            if (smer==0) q1 += 1;
            else         q2 += 1;
         } else {
            x1  = x2;
            y1  = y2;
            x1d = x2d;
            y1d = y2d;
         }
      } while ((q1!=qv||(q2!=(qv-1) && q2!=w2)||smer!=0||flag!=1) &&
               ((q1!=(qv-1) && q1!=w1)||q2!=qv||smer!=1||flag!=1) &&
                uhl!=1);
      if (qv<=w2&&qv<=w1) {
         if (uhl==0) {
            if (smer==0) {
               smer = 1;
               q1   = 0;
               flag = 0;
               q2   = qv;
               xp1  = x1;
               yp1  = y1;
               goto l2;
            } else {
               smer = 0;
               uhl  = 1;
               q1   = qv;
               q2   = qv-1;
               xp2  = x1;
               yp2  = y1;
               x1   = xp1;
               y1   = yp1;
               flag = 1;
               goto l2;
            }
         } else {
            if (smer==0) {
               smer = 1;
               x1   = xp2;
               y1   = yp2;
               q1   = qv-1;
               q2   = qv;
               flag = 1;
               goto l2;
            }
         }
      }
      q2 = 0;
      qv += 1;
      q1 = qv;
      if (q1>w1) {
         q2   = qv;
         q1   = 0;
         smer = 1;
         flag = 0;
         uhl  = 0;
         if (q2<=w2) goto l2;
      }
   } while (q1<=w1&&q2<=w2);
   if (fChanmarkEnDis==kChannelMarksDrawn ||
       fChanlineEnDis==kChannelGridDrawn) {
      line->SetLineWidth(1);
      for (i=fBx1;i<=fBx2;i++) fEnvelope[i] = fBy2;
      turni = 0;
      turnj = 0;
      Transform(w1,0,0);
      x1 = fXt;
      Transform(0,0,0);
      x2 = fXt;
      Transform(0,w2,0);
      x3 = fXt;
      if (x2>=x1) turnj = 1;
      if (x3>=x2) turni = 1;
      q1 = 1;
      q2 = 0;
      qv = 1;
      do {
         uhl  = 0;
         smer = 0;
         flag = 0;
lc2:
         if (turni==1) {
            i = q1;
         } else {
            i = w1-q1;
         }
         if (turnj==1) {
            j = q2;
         } else {
            j = w2-q2;
         }
         Transform(i,j,0);
         x1 = fXt;
         y1 = fYt;
         Transform(i,j,-1);
         x1d = fXt;
         y1d = fYt;
         do {
            if (flag==0) {
               flag = 1;
               if (smer==0) q1 -= 1;
               else         q2 -= 1;
            } else {
               flag = 0;
               if (smer==0) q2 += 1;
               else         q1 += 1;
            }
            if (turni==1) {
               i = q1;
            } else {
               i = w1-q1;
            }
            if (turnj==1) {
               j = q2;
            } else {
               j = w2-q2;
            }
            Transform(i,j,0);
            x2 = fXt;
            y2 = fYt;
            if (flag==1) {
               x  = x1;
               y  = y1;
               x1 = x2;
               y1 = y2;
               x2 = x;
               y2 = y;
            }
            Envelope(x1,y1,x2,y2);
            if (fChanlineEnDis==kChannelGridDrawn) {
               if (fLine!=0) {
                  if (fLine==1) {
                     fXe = x2;
                     fYe = y2;
                  }
                  line->SetLineColor(fChanlineColor);
                  line->PaintLine(gPad->PixeltoX(fXs),gPad->PixeltoY(fYs)+1,gPad->PixeltoX(fXe),gPad->PixeltoY(fYe)+1);
               }
            }
            if (fChanmarkEnDis==kChannelMarksDrawn) {
               if (y1<=fEnvelope[x1]) {
                  DrawMarker(x1,y1,fChanmarkWidth,fChanmarkHeight,fChanmarkStyle);
               }
            }
            if (flag==1) {
               x   = x1;
               y   = y1;
               x1  = x2;
               y1  = y2;
               x2  = x;
               y2  = y;
               x   = x1d;
               y   = y1d;
               x1d = x2d;
               y1d = y2d;
               x2d = x;
               y2d = y;
               if (smer==0) q1 += 1;
               else         q2 += 1;
            } else {
               x1  = x2;
               y1  = y2;
               x1d = x2d;
               y1d = y2d;
            }
         } while ((q1!=qv||(q2!=(qv-1)&&q2!=w2)||smer!=0||flag!=1) &&
                  ((q1!=(qv-1)&&q1!=w1)||q2!=qv||smer!=1||flag!=1) &&
                   uhl!=1);
         if (qv<=w2&&qv<=w1) {
            if (uhl==0) {
               if (smer==0) {
                  smer = 1;
                  q1   = 0;
                  flag = 0;
                  q2   = qv;
                  xp1  = x1;
                  yp1  = y1;
                  goto lc2;
               } else {
                  smer = 0;
                  uhl  = 1;
                  q1   = qv;
                  q2   = qv-1;
                  xp2  = x1;
                  yp2  = y1;
                  x1   = xp1;
                  y1   = yp1;
                  flag = 1;
                  goto lc2;
               }
            } else {
               if (smer==0) {
                  smer = 1;
                  x1   = xp2;
                  y1   = yp2;
                  q1   = qv-1;
                  q2   = qv;
                  flag = 1;
                  goto lc2;
               }
            }
         }
         q2  = 0;
         qv += 1;
         q1  = qv;
         if (q1>w1) {
            q2   = qv;
            q1   = 0;
            smer = 1;
            flag = 0;
            uhl  = 0;
            if (q2<=w2) goto lc2;
         }
      } while (q1<=w1&&q2<=w2);
   }

   // Paint axis.
   static char chopt[10] = "";
   if (fViewAngle==0) {
      axis->PaintAxis(p101x, p101y, p111x, p111y, bmin, bmax, ndivx, "");
      axis->PaintAxis(p011x, p011y, p111x, p111y, bmin, bmax, ndivx, "");
      if(fAlpha+fBeta<90)
         axis->PaintAxis(p110x, p110y, p111x, p111y, bmin, bmax, ndivx, "");
      ndivx = fH2->GetXaxis()->GetNdivisions();
      bmin  = fH2->GetXaxis()->GetXmin();
      bmax  = fH2->GetXaxis()->GetXmax();
      xaxis->SetLabelOffset(xaxis->GetLabelOffset()-xaxis->GetTickSize());
      chopt[0] = 0; strlcat(chopt, "SDH-",10);
      if (ndivx < 0) {
         strlcat(chopt, "N",10);
         ndivx = -ndivx;
      }
      xaxis->PaintAxis(p010x, p010y, p110x, p110y, bmin, bmax, ndivx, chopt);
      ndivy = fH2->GetYaxis()->GetNdivisions();
      bmin  = fH2->GetYaxis()->GetXmin();
      bmax  = fH2->GetYaxis()->GetXmax();
      yaxis->SetLabelOffset(yaxis->GetLabelOffset()+yaxis->GetTickSize());
      chopt[0] = 0; strlcat(chopt, "SDH+",10);
      if (ndivy < 0) {
         strlcat(chopt, "N",10);
         ndivy = -ndivy;
      }
      yaxis->PaintAxis(p100x, p100y, p110x, p110y, bmin, bmax, ndivy, chopt);
      if(fAlpha+fBeta<90){
         ndivz = 510;
         bmin  = fZmin;
         bmax  = fZmax;
         zaxis->SetLabelOffset(zaxis->GetLabelOffset()-zaxis->GetTickSize());
         if (fZscale==kZScaleLog) {
            ndivz = 510;
            if (bmin <= 0) bmin=0.001*bmax;
            zaxis->PaintAxis(p010x, p010y, p011x, p011y, bmin, bmax, ndivz, "G+");
         } else if (fZscale==kZScaleSqrt) {
            TF1 *f1=new TF1("f1","sqrt(x)",bmin,bmax);
            TGaxis *a1 = new TGaxis(p010x, p010y, p011x, p011y, "f1", ndivz, "SDH+");
            a1->SetLabelOffset(a1->GetLabelOffset()-a1->GetTickSize());
            a1->Paint();
            delete f1;
            delete a1;
         } else {
            chopt[0] = 0; strlcat(chopt, "SDH+",10);
            if (ndivz < 0) {
               strlcat(chopt, "N",10);
               ndivz = -ndivz;
            }
            zaxis->PaintAxis(p010x, p010y, p011x, p011y, bmin, bmax, ndivz, chopt);
         }
      }
   } else if (fViewAngle==90) {
      axis->PaintAxis(p001x, p001y, p101x, p101y, bmin, bmax, ndivx, "");
      axis->PaintAxis(p111x, p111y, p101x, p101y, bmin, bmax, ndivx, "");
      if(fAlpha+fBeta<90)
         axis->PaintAxis(p100x, p100y, p101x, p101y, bmin, bmax, ndivx, "");
      ndivx = fH2->GetXaxis()->GetNdivisions();
      bmin  = fH2->GetXaxis()->GetXmin();
      bmax  = fH2->GetXaxis()->GetXmax();
      xaxis->SetLabelOffset(xaxis->GetLabelOffset()+xaxis->GetTickSize());
      chopt[0] = 0; strlcat(chopt, "SDH+",10);
      if (ndivx < 0) {
         strlcat(chopt, "N",10);
         ndivx = -ndivx;
      }
      xaxis->PaintAxis(p000x, p000y, p100x, p100y, bmin, bmax, ndivx, chopt);
      ndivy = fH2->GetYaxis()->GetNdivisions();
      bmin  = fH2->GetYaxis()->GetXmin();
      bmax  = fH2->GetYaxis()->GetXmax();
      yaxis->SetLabelOffset(yaxis->GetLabelOffset()+yaxis->GetTickSize());
      chopt[0] = 0; strlcat(chopt, "SDH+",10);
      if (ndivy < 0) {
         strlcat(chopt, "N",10);
         ndivy = -ndivy;
      }
      yaxis->PaintAxis(p100x, p100y, p110x, p110y, bmin, bmax, ndivy, chopt);
      if(fAlpha+fBeta<90){
         ndivz = 510;
         bmin  = fZmin;
         bmax  = fZmax;
         zaxis->SetLabelOffset(zaxis->GetLabelOffset()-zaxis->GetTickSize());
         if (fZscale==kZScaleLog) {
            ndivz = 510;
            if (bmin <= 0) bmin=0.001*bmax;
            zaxis->PaintAxis(p110x, p110y, p111x, p111y, bmin, bmax, ndivz, "G+");
         } else if (fZscale==kZScaleSqrt) {
            TF1 *f1=new TF1("f1","sqrt(x)",bmin,bmax);
            TGaxis *a1 = new TGaxis(p110x, p110y, p111x, p111y, "f1", ndivz, "SDH+");
            a1->SetLabelOffset(a1->GetLabelOffset()-a1->GetTickSize());
            a1->Paint();
            delete f1;
            delete a1;
         } else {
            chopt[0] = 0; strlcat(chopt, "SDH+",10);
            if (ndivz < 0) {
               strlcat(chopt, "N",10);
               ndivz = -ndivz;
            }
            zaxis->PaintAxis(p110x, p110y, p111x, p111y, bmin, bmax, ndivz, chopt);
         }
      }
   } else if (fViewAngle==180) {
      axis->PaintAxis(p011x, p011y, p001x, p001y, bmin, bmax, ndivx, "");
      axis->PaintAxis(p101x, p101y, p001x, p001y, bmin, bmax, ndivx, "");
      if(fAlpha+fBeta<90)
         axis->PaintAxis(p000x, p000y, p001x, p001y, bmin, bmax, ndivx, "");
      ndivx = fH2->GetXaxis()->GetNdivisions();
      bmin  = fH2->GetXaxis()->GetXmin();
      bmax  = fH2->GetXaxis()->GetXmax();
      xaxis->SetLabelOffset(xaxis->GetLabelOffset()+xaxis->GetTickSize());
      chopt[0] = 0; strlcat(chopt, "SDH+",10);
      if (ndivx < 0) {
         strlcat(chopt, "N",10);
         ndivx = -ndivx;
      }
      xaxis->PaintAxis(p000x, p000y, p100x, p100y, bmin, bmax, ndivx, chopt);
      ndivy = fH2->GetYaxis()->GetNdivisions();
      bmin  = fH2->GetYaxis()->GetXmin();
      bmax  = fH2->GetYaxis()->GetXmax();
      yaxis->SetLabelOffset(yaxis->GetLabelOffset()-yaxis->GetTickSize());
      chopt[0] = 0; strlcat(chopt, "SDH-",10);
      if (ndivy < 0) {
         strlcat(chopt, "N",10);
         ndivy = -ndivy;
      }
      yaxis->PaintAxis(p000x, p000y, p010x, p010y, bmin, bmax, ndivy, chopt);
      if(fAlpha+fBeta<90){
         ndivz = 510;
         bmin  = fZmin;
         bmax  = fZmax;
         zaxis->SetLabelOffset(zaxis->GetLabelOffset()-zaxis->GetTickSize());
         if (fZscale==kZScaleLog) {
            ndivz=510;
            if (bmin <= 0) bmin=0.001*bmax;
            zaxis->PaintAxis(p100x, p100y, p101x, p101y, bmin, bmax, ndivz, "G+");
         } else if (fZscale==kZScaleSqrt) {
            TF1 *f1=new TF1("f1","sqrt(x)",bmin,bmax);
            TGaxis *a1 = new TGaxis(p100x, p100y, p101x, p101y, "f1", ndivz, "SDH+");
            a1->SetLabelOffset(a1->GetLabelOffset()-a1->GetTickSize());
            a1->Paint();
            delete f1;
            delete a1;
         } else {
            chopt[0] = 0; strlcat(chopt, "SDH+",10);
            if (ndivz < 0) {
               strlcat(chopt, "N",10);
               ndivz = -ndivz;
            }
            zaxis->PaintAxis(p100x, p100y, p101x, p101y, bmin, bmax, ndivz, chopt);
         }
      }
   } else if (fViewAngle==270) {
      axis->PaintAxis(p111x, p111y, p011x, p011y, bmin, bmax, ndivx, "");
      axis->PaintAxis(p001x, p001y, p011x, p011y, bmin, bmax, ndivx, "");
      if(fAlpha+fBeta<90)
         axis->PaintAxis(p010x, p010y, p011x, p011y, bmin, bmax, ndivx, "");
      ndivx = fH2->GetXaxis()->GetNdivisions();
      bmin  = fH2->GetXaxis()->GetXmin();
      bmax  = fH2->GetXaxis()->GetXmax();
      xaxis->SetLabelOffset(xaxis->GetLabelOffset()-xaxis->GetTickSize());
      chopt[0] = 0; strlcat(chopt, "SDH-",10);
      if (ndivx < 0) {
         strlcat(chopt, "N",10);
         ndivx = -ndivx;
      }
      xaxis->PaintAxis(p010x, p010y, p110x, p110y, bmin, bmax, ndivx, chopt);
      ndivy = fH2->GetYaxis()->GetNdivisions();
      bmin  = fH2->GetYaxis()->GetXmin();
      bmax  = fH2->GetYaxis()->GetXmax();
      yaxis->SetLabelOffset(yaxis->GetLabelOffset()-yaxis->GetTickSize());
      chopt[0] = 0; strlcat(chopt, "SDH-",10);
      if (ndivy < 0) {
         strlcat(chopt, "N",10);
         ndivy = -ndivy;
      }
      yaxis->PaintAxis(p000x, p000y, p010x, p010y, bmin, bmax, ndivy, chopt);
      if(fAlpha+fBeta<90){
         ndivz = 510;
         bmin  = fZmin;
         bmax  = fZmax;
         zaxis->SetLabelOffset(zaxis->GetLabelOffset()-zaxis->GetTickSize());
         if (fZscale==kZScaleLog) {
            ndivz = 510;
            if (bmin <= 0) bmin=0.001*bmax;
            zaxis->PaintAxis(p000x, p000y, p001x, p001y, bmin, bmax, ndivz, "G+");
         } else if (fZscale==kZScaleSqrt) {
            TF1 *f1=new TF1("f1","sqrt(x)",bmin,bmax);
            TGaxis *a1 = new TGaxis(p000x, p000y, p001x, p001y, "f1", ndivz, "SDH+");
            a1->SetLabelOffset(a1->GetLabelOffset()-a1->GetTickSize());
            a1->Paint();
            delete f1;
            delete a1;
         } else {
            chopt[0] = 0; strlcat(chopt, "SDH+",10);
            if (ndivz < 0) {
               strlcat(chopt, "N",10);
               ndivz = -ndivz;
            }
            zaxis->PaintAxis(p000x, p000y, p001x, p001y, bmin, bmax, ndivz, "SDH+");
         }
      }
   }

   // End.
   delete axis;
   delete xaxis;
   delete yaxis;
   delete zaxis;
}


////////////////////////////////////////////////////////////////////////////////
/// Sets display group mode and display mode:
///    - modeGroup - the following group modes can be set: simple modes-kPicture2ModeGroupSimple, modes with shading according to light-kPicture2ModeGroupLight, modes with shading according to channels counts-kPicture2ModeGroupHeight, modes of combination of shading according to light and to channels counts-kPicture2ModeGroupLightHeight
///    - displayMode - posible display modes are: points, grid, contours, bars, x_lines, y_lines, bars_x, bars_y, needles, surface, triangles

void TSpectrum2Painter::SetDisplayMode(Int_t modeGroup,Int_t displayMode)
{
   if (modeGroup>=kModeGroupSimple&&modeGroup<=kModeGroupLightHeight) {
      if (displayMode>=kDisplayModePoints&&displayMode<=kDisplayModeTriangles) {
         fModeGroup   = modeGroup;
         fDisplayMode = displayMode;
      }
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Sets pen attributes:
///    - color - color of spectrum
///    - style - style of pen (solid, dash, dot, dash-dot)
///    - width - width of pen in pixels

void TSpectrum2Painter::SetPenAttr(Int_t color,Int_t style,Int_t width)
{
   if (color>=0 && style >=kPenStyleSolid && style <= kPenStyleDashDot && width > 0) {
      fPenColor = color;
      fPenDash  = style;
      fPenWidth = width;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Sets nodes in both directions:
///    - nodesx, nodesy, only the bins at the nodes points are displayed

void TSpectrum2Painter::SetNodes(Int_t nodesx,Int_t nodesy)
{
   if (nodesx>1&&nodesy>1) {
      fNodesx = nodesx;
      fNodesy = nodesy;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Sets angles of the view:
///    - alpha - angles of display,alpha+beta must be less or equal to 90, alpha- angle between base line of Canvas and right lower edge of picture base plane
///    - beta - angle between base line of Canvas and left lower edge of picture base plane
///    - view - rotation angle of the view, it can be 0, 90, 180, 270 degrees

void TSpectrum2Painter::SetAngles(Int_t alpha,Int_t beta,Int_t view)
{
   if (alpha>=0&&alpha<=90&&beta>=0&&beta<=90&&alpha+beta<=90) {
      fAlpha = alpha;
      fBeta  = beta;
   }
   view = view/90;
   if (view>=0&&view<=3) fViewAngle = view*90;
}


////////////////////////////////////////////////////////////////////////////////
/// Sets z-axis scale:
///    - scale - linear, sqrt or log

void TSpectrum2Painter::SetZScale(Int_t scale)
{
   if (scale>=kZScaleLinear&&scale<=kZScaleSqrt) {
      fZscale = scale;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Sets shading color algorithm:
///    - olorAlgorithm - applies only for rainbowed display modes
///    (rgb smooth algorithm, rgb modulo color component, cmy smooth algorithm,
///    - cmy modulo color component, cie smooth algorithm,
///    - cie modulo color component, yiq smooth algorithm,
///    - yiq modulo color component, hsv smooth algorithm,
///    - hsv modulo color component, it does not apply for simple display modes
///      algorithm group

void TSpectrum2Painter::SetColorAlgorithm(Int_t colorAlgorithm)
{
   if (fModeGroup!=kModeGroupSimple) {
      if (colorAlgorithm>=kColorAlgRgbSmooth&&colorAlgorithm<=kColorAlgHvsModulo) fColorAlg = colorAlgorithm;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Sets color increments between two color levels for r, g, b components:
///    - r, g, b - color increments between two color levels

void TSpectrum2Painter::SetColorIncrements(Double_t r,Double_t g,Double_t b)
{
   if (r>=0&&r<=255&&g>=0&&g<=255&&b>=0&&b<=255) {
      fRainbow1Step = r;
      fRainbow2Step = g;
      fRainbow3Step = b;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Sets position of fictive light source in 3D space:
///    - x, y, z

void TSpectrum2Painter::SetLightPosition(Int_t x,Int_t y,Int_t z)
{
   if (x>=0&&y>=0&&z>=0) {
      fXlight = x;
      fYlight = y;
      fZlight = z;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Sets on/off shading and shadow switches:
///   - shading - determines whether the picture will shaded, smoothed (no shading, shading), for rainbowed display modes only
///   - shadow - determines whether shadow will be drawn, for rainbowed display modes with shading according to light

void TSpectrum2Painter::SetShading(Int_t shading,Int_t shadow)
{
   if (fModeGroup!=kModeGroupSimple) {
      if (shading==kNotShaded||shading==kShaded) fShading = shading;
      if (shadow==kShadowsNotPainted||shadow==kShadowsPainted) fShadow = shadow;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Sets on/off Bezier smoothing:
///   - bezier - determines Bezier interpolation (applies only for simple
///     display modes group for grid, x_lines, y_lines display modes)

void TSpectrum2Painter::SetBezier(Int_t bezier)
{
   if (fDisplayMode==kDisplayModeGrid || fDisplayMode==kDisplayModeLinesX ||
       fDisplayMode==kDisplayModeLinesY) {
      if (bezier==kBezierInterpol||bezier==kNoBezierInterpol) fBezier = bezier;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Sets width between horizontal slices:
///   - width - width between contours, applies only for contours display mode

void TSpectrum2Painter::SetContourWidth(Int_t width)
{
   if (width>=1) fContWidth = width;
}


////////////////////////////////////////////////////////////////////////////////
/// Sets weight between shading according to fictive light source and according to channels counts:
///   - weight - weight between shading according to fictive light source and according to channels counts, applies only for kPicture2ModeGroupLightHeight modes group

void TSpectrum2Painter::SetLightHeightWeight(Double_t weight)
{
   if (fModeGroup==kModeGroupLightHeight) {
      if (weight>=0&&weight<=1) fLHweight = weight;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Sets enables/disables drawing of channel marks and sets their attributes:
///   - enable - decides whether the channel marks are shown or not
///   - color - color of channel marks
///   - width - width of channel marks in pixels
///   - height - height of channel marks in pixels
///   - style - style of channel marks (dot, cross, star, rectangle, X, diamond, triangle)

void TSpectrum2Painter::SetChanMarks(Int_t enable,Int_t color,Int_t width,Int_t height,Int_t style)
{
   if (enable==kChannelMarksNotDrawn||enable==kChannelMarksDrawn) {
      if (enable==kChannelMarksDrawn) {
         if (style>=kChannelMarksStyleDot&&style<=kChannelMarksStyleTriangle) {
            fChanmarkStyle  = style;
            fChanmarkColor  = color;
            if (width>=4) {
               fChanmarkWidth  = width;
            }
            else fChanmarkWidth  = 4;
            if (height>=4) {
               fChanmarkHeight = height;
            }
            else fChanmarkHeight  = 4;
         }
      }
      fChanmarkEnDis = enable;
   }
}


////////////////////////////////////////////////////////////////////////////////
///   This function sets enables/disables drawing of channel grid and sets its color:
///         -enable - decides whether the channel grid is shown or not
///         -color - color of channel grid

void TSpectrum2Painter::SetChanGrid(Int_t enable,Int_t color)
{
   if (enable==kChannelGridNotDrawn||enable==kChannelGridDrawn) {
      if (enable==kChannelGridDrawn) {
         fChanlineColor=color;
      }
      fChanlineEnDis=enable;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Gets display group mode and display mode:
///    -modeGroup - the following group modes might have been set: simple modes-kPicture2ModeGroupSimple, modes with shading according to light-kPicture2ModeGroupLight, modes with shading according to channels counts-kPicture2ModeGroupHeight, modes of combination of shading according to light and to channels counts-kPicture2ModeGroupLightHeight
///    -displayMode - display modes that might have been set: points, grid, contours, bars, x_lines, y_lines, bars_x, bars_y, needles, surface, triangles

void TSpectrum2Painter::GetDisplayMode(Int_t &modeGroup,Int_t &displayMode)
{
   modeGroup   = fModeGroup;
   displayMode = fDisplayMode;
}


////////////////////////////////////////////////////////////////////////////////
/// Gets pen attributes:
///    -color - color of spectrum
///    -style - style of pen (solid, dash, dot, dash-dot)
///    -width - width of pen in pixels

void TSpectrum2Painter::GetPenAttr(Int_t &color, Int_t &style, Int_t &width)
{
   color = fPenColor;
   style = fPenDash;
   width = fPenWidth;
}


////////////////////////////////////////////////////////////////////////////////
/// Gets nodes in both directions:
///    - nodesx, nodesy, only the bins at the nodes points are displayed

void TSpectrum2Painter::GetNodes(Int_t &nodesx,Int_t &nodesy)
{
   nodesx = fNodesx;
   nodesy = fNodesy;
}


////////////////////////////////////////////////////////////////////////////////
/// Gets angles of the view:
///    - alpha - angle between base line of Canvas and right lower edge of picture base plane
///    - beta - angle between base line of Canvas and left lower edge of picture base plane
///    - view - rotation angle of the view, it can be 0, 90, 180, 270 degrees

void TSpectrum2Painter::GetAngles(Int_t &alpha,Int_t &beta,Int_t &view)
{
   alpha = fAlpha;
   beta  = fBeta;
   view  = fViewAngle;
}


////////////////////////////////////////////////////////////////////////////////
/// Gets z-axis scale:
///    - scale - it can be linear, sqrt or log

void TSpectrum2Painter::GetZScale(Int_t &scale)
{
   scale = fZscale;
}


////////////////////////////////////////////////////////////////////////////////
/// Gets shading color algorithm:
///    - colorAlgorithm - rgb smooth algorithm, rgb modulo color component, cmy smooth algorithm, cmy modulo color component, cie smooth algorithm, cie modulo color component, yiq smooth algorithm, yiq modulo color component, hsv smooth algorithm, hsv modulo algorithm

void TSpectrum2Painter::GetColorAlgorithm(Int_t &colorAlgorithm)
{
   colorAlgorithm = fColorAlg;
}


////////////////////////////////////////////////////////////////////////////////
/// Gets color increments between two color levels for r, g, b components:
///    - r, g, b - color increments between two color levels

void TSpectrum2Painter::GetColorIncrements(Double_t &r,Double_t &g,Double_t &b)
{
   r = fRainbow1Step;
   g = fRainbow2Step;
   b = fRainbow3Step;
}


////////////////////////////////////////////////////////////////////////////////
/// Gets position of fictive light source in 3D space:
///    - x, y, z

void TSpectrum2Painter::GetLightPosition(Int_t &x,Int_t &y,Int_t &z)
{
   x = fXlight;
   y = fYlight;
   z = fZlight;
}


////////////////////////////////////////////////////////////////////////////////
/// Gets shading and shadow switches:
///    - shading - determines whether the picture will shaded, smoothed (no shading, shading), for rainbowed display modes only
///    - shadow - determines whether shadow will be drawn, for rainbowed display modes with shading according to light

void TSpectrum2Painter::GetShading(Int_t &shading,Int_t &shadow)
{
   shading = fShading;
   shadow  = fShadow;
}


////////////////////////////////////////////////////////////////////////////////
/// Gets Bezier smoothing switch:
///    - bezier - determines Bezier interpolation (applies only for simple display modes group for grid, x_lines, y_lines display modes)

void TSpectrum2Painter::GetBezier(Int_t &bezier)
{
   bezier = fBezier;
}


////////////////////////////////////////////////////////////////////////////////
/// Gets width between horizontal slices:
///    - width - width between contours, applies only for contours display mode

void TSpectrum2Painter::GetContourWidth(Int_t &width)
{
   width = fContWidth;
}


////////////////////////////////////////////////////////////////////////////////
/// Gets weight between shading according to fictive light source and according to channels counts:
///    - weight - weight between shading according to fictive light source and according to channels counts, applies only for kPicture2ModeGroupLightHeight modes group

void TSpectrum2Painter::GetLightHeightWeight(Double_t &weight)
{
   weight = fLHweight;
}


////////////////////////////////////////////////////////////////////////////////
/// Gets drawing attributes for channel marks:
///    - enable - decides whether the channel marks are shown or not
///    - color - color of channel marks
///    - width - width of channel marks in pixels
///    - height - height of channel marks in pixels
///    - style - style of channel marks (dot, cross, star, rectangle, X, diamond, triangle)

void TSpectrum2Painter::GetChanMarks(Int_t &enable,Int_t &color,Int_t &width,Int_t &height,Int_t &style)
{
   style  = fChanmarkStyle,width=fChanmarkWidth,height=fChanmarkHeight,color=fChanmarkColor;
   enable = fChanmarkEnDis;
}


////////////////////////////////////////////////////////////////////////////////
/// This function gets attributes for drawing channel:
///    - enable - decides whether the channel grid is shown or not
///    - color - color of channel grid

void TSpectrum2Painter::GetChanGrid(Int_t &enable,Int_t &color)
{
   color  = fChanlineColor;
   enable = fChanlineEnDis;
}


////////////////////////////////////////////////////////////////////////////////
/// This function allows to set all the possible options available in
/// TSpectrum2Painter and paint "h2".
///
/// TSpectrum2Painter offers a large set of options/attributes. In the
/// "option" parameter each of them can be set to specific values via
/// "operators" in the option itself. for instance on can do:
///
///     h2->Draw("SPEC a(30,30,0)");
///
/// to draw h2 with TSpectrum2Painter using all the default attributes except
/// the viewing angles. The operators' names are case insensitive (one can use
/// "a" or "A"). Operators parameters are separated by ",". The operators can
/// be put is any order in the option. Operators must be separated by " ".
/// No " " characters should be put in an operator. This help page describes
/// all the available operators.
///
/// The way "h2" will be painted is controlled with 2 parameters the "Display
/// modes groups" and the "Display Modes".
///
///   "Display modes groups" can take the following values:
///
///  - 0 = Simple      - it covers simple display modes using one color only
///  - 1 = Light       - in this group the shading is carried out according to
///                      the position of the fictive light source
///  - 2 = Height      - in this group the shading is carried out according to
///                      the channel contents
///  - 3 = LightHeight - combination of two previous shading algorithms. One
///                      can control the weight between both algorithms.
///
/// "Display modes" can take the following values:
///
///  -  1 = Points.
///  -  2 = Grid.
///  -  3 = Contours.
///  -  4 = Bars.
///  -  5 = LinesX.
///  -  6 = LinesY.
///  -  7 = BarsX.
///  -  8 = BarsY.
///  -  9 = Needles.
///  - 10 = Surface.
///  - 11 = Triangles.
///
/// Using this function these parameters can be set using the "dm" operator
/// in the option. Example:
///
///     h2->Draw("SPEC dm(1,2)");
///
/// will paint the 2D histogram h2 using the "Light Display mode group" and
/// the "Grid Display mode". The following table summarises all the possible
/// combinations of the "Display modes groups" and the "Display modes".
///
/// |           | Simple | Light | Height | Light-Height |
/// |-----------|--------|-------|--------|--------------|
/// | Points    |   X    |   X   |   X    |      X       |
/// | Grid      |   X    |   X   |   X    |      X       |
/// | Contours  |   X    |   -   |   X    |      -       |
/// | Bars      |   X    |   -   |   X    |      -       |
/// | LinesX    |   X    |   X   |   X    |      X       |
/// | LinesY    |   X    |   X   |   X    |      X       |
/// | BarsX     |   X    |   -   |   X    |      -       |
/// | BarsY     |   X    |   -   |   X    |      -       |
/// | Needles   |   X    |   -   |   -    |      -       |
/// | Surface   |   -    |   X   |   X    |      X       |
/// | Triangles |   X    |   X   |   X    |      X       |
///
/// The Pen Attributes can be changed using pa(color, style, width). Example:
///
///     h2->Draw("SPEC dm(1,2) pa(2,1,2)");
///
/// sets the line color to 2, line type to 1 and line width to2. Note that if
/// pa() is not specified, the histogram "h2" line attributes are used. Note
/// also that operators for SPEC option can be cumulated and specified in
/// any order.
///
/// The buffer size can be change with bf(size). Example:
///
///     h2->Draw("SPEC bf(8000)");
///
/// The spectrum painter needs a buffer to paint the spectrum. By default the
/// buffer size is set to 1600. In most cases this buffer size is enough. But
/// if the canvas size is very big, for instance 8000x5000 this buffer size is
/// too small. An error message is issued telling to use the option bf().
///
/// The number of nodes can be change with n(nodesx, nodesy). Example:
///
///     h2->Draw("SPEC n(40,40)");
///
/// Sometimes the displayed region is rather large. When displaying all
/// channels pictures become very dense and complicated. It is very difficult
/// to understand the overall shape of data. "n(nx,ny)" allows to change the
/// density of displayed channels. Only the channels coinciding with given
/// nodes are displayed.
///
/// The visualization angles can be changed with a(alpha, beta, view).
/// Example:
///
///     h2->Draw("SPEC n(40,40) dm(0,1) a(30,30,0)");
///
/// Alpha defines the angle between bottom horizontal screen line and the
/// displayed space on the right side of the picture and beta on the left
/// side, respectively. One can rotate the 3-d space around vertical axis
/// going through the center of it employing the view parameter. Allowed
/// values are 0, 90, 180 and 270 degrees.
///
/// zs(scale) changes the scale of the Z-axis Possible values are:
///
///  - 0 = Linear (default)
///  - 1 = Log
///  - 2 = Sqrt
///
/// If gPad->SetLogz() has been set, log scale on Z-axis is set automatically.
/// No need to use the zs() operator. Note that the X and Y axis are always
/// linear.
///
/// ci(r,g,b), were r,g and b are floats defines the colors increments.
/// For sophisticated shading (Light, Height and LightHeight Display Modes
/// Groups) the color palette starts from the basic pen color (see pa()
/// function). There is a predefined number of color levels (256). Color in
/// every level is calculated by adding the increments of the r, g, b
/// components to the previous level. Using this function one can change the
/// color increments between two neighbouring color levels. The function does
/// not apply dor the Simple Display Modes Group. The default values are:
/// (1,1,1).
///
/// ca(color_algorithm) allows to choose the Color Algorithm.
/// To define the colors one can employ one of the following color algorithms
/// (RGB, CMY, CIE, YIQ, HVS models). When the level of a component reaches
/// the limit value one can choose either smooth transition (by decreasing
/// the limit value) or a sharp modulo transition (continuing with 0 value).
/// This makes possible to realize various visual effects. One can choose from
/// the following set of the algorithms:
///
///  - 0 = RGB Smooth
///  - 1 = RGB Modulo
///  - 2 = CMY Smooth
///  - 3 = CMY Modulo
///  - 4 = CIE Smooth
///  - 5 = CIE Modulo
///  - 6 = YIQ Smooth
///  - 7 = YIQ Modulo
///  - 8 = HVS Smooth
///  - 9 = HVS Modulo
///
/// This function does not apply on Simple display modes group. Default
/// value is 0. Example:
///
///     h2->Draw("SPEC c1(4) dm(0,1) a(30,30,0)");
///
/// choose CMY Modulo to paint the "h2" histogram.
///
/// lp(x,y,z) set the light position.
/// In Light and LightHeight display modes groups the color palette is
/// calculated according to the fictive light source position in 3-d space.
/// Using this function one can change the position of the source and thus
/// to achieve various graphical effects. This function does not apply for
/// Simple and Height display modes groups. Default is:
/// lp(1000,1000,100).
///
/// s(shading,shadow) allows to set the shading.
/// The picture's surface is composed of triangles. If desired the edges of
/// the neighbouring triangles can be smoothed (shaded). If desired the
/// display of the shadow can be painted as well. The function does not apply
/// for Simple display modes group. The possible values for shading are:
///
///  - 0 = Not Shaded
///  - 1 = Shaded
///
/// The possible values for shadow are:
///
///  - 0 = Shadows are not painted
///  - 1 = Shadows are painted
///
/// Default values: s(1,0)
///
/// b(bezier) set the Bezier smoothing.
/// For Simple display modes group and for Grid, LinesX and LinesY display
/// modes one can smooth data using Bezier smoothing algorithm. The function
/// does not apply for other display modes groups and display modes. Possible
/// values are:
///
///  - 0 = No bezier smoothing
///  - 1 = Bezier smoothing
///
/// Default value is: b(0)
///
/// cw(width) set the contour width.
/// This function applies only for the Contours display mode. One can change
/// the width between horizontal slices and thus their density.
/// Default value: cw(50)
///
/// lhw(weight) set the light height weight.
/// For LightHeight display modes group one can change the weight between
/// both shading algorithms. The function does not apply for other display
/// modes groups. Default value is lhw(0.5).
///
/// cm(enable,color,width,height,style) allows to draw a marker on each node.
/// In addition to the surface drawn using any above given algorithm one can
/// display channel marks. One can control the color as well as the width,
/// height (in pixels) and the style of the marks. The parameter enable can
/// be set to
///
///  - 0 =  Channel marks are not drawn
///  - 1 =  Channel marks drawn
///
/// The possible styles can be chosen from the set:
///
///  - 1 = Dot
///  - 2 = Cross
///  - 3 = Star
///  - 4 = Rectangle
///  - 5 = X
///  - 6 = Diamond
///  - 7 = Triangle.
///
/// cg(enable,color) channel grid.
/// In addition to the surface drawn using any above given algorithm one can
/// display grid using the color parameter. The parameter enable can be set
/// to:
///
///  - 0 = Grid not drawn
///  - 1 = Grid drawn
///
/// See example spectrumpainter.C
///
/// \image html spectrumpainter.gif

void TSpectrum2Painter::PaintSpectrum(TH2* h2, Option_t *option, Int_t bs)
{
   TString opt = option;

   TSpectrum2Painter sp(h2, bs);

   if (gPad->GetLogz()) sp.SetZScale(kZScaleLog);
   sp.SetPenAttr(h2->GetLineColor(), h2->GetLineStyle(), h2->GetLineWidth());

   TString token;
   Int_t i1, i2, i3, i4, i5;
   Double_t f1, f2, f3;
   Ssiz_t from = 4;

   // Decode the paint options.
   while (opt.Tokenize(token, from, "[ (]")) {

      // Display Mode
      if (token=="dm") {
         opt.Tokenize(token, from, ","); i1 = token.Atoi();
         if (i1<0 || i1>3) {
            printf("PaintSpectrum: Display modes groups should be in the [0,3] range\n");
            i1 = 0;
         }
         opt.Tokenize(token, from, ")"); i2 = token.Atoi();
         if (i2<1 || i2>11) {
            printf("PaintSpectrum: Display modes should be in the [1,11] range\n");
            i2 = 1;
         }
         sp.SetDisplayMode(i1, i2);

      // Pen Attributes
      } else if (token=="pa") {
         opt.Tokenize(token, from, ","); i1 = token.Atoi();
         opt.Tokenize(token, from, ","); i2 = token.Atoi();
         opt.Tokenize(token, from, ")"); i3 = token.Atoi();
         sp.SetPenAttr(i1, i2, i3);

      // Nodes
      } else if (token=="n") {
         opt.Tokenize(token, from, ","); i1 = token.Atoi();
         opt.Tokenize(token, from, ")"); i2 = token.Atoi();
         sp.SetNodes(i1, i2);

      // Color Algorithm
      } else if (token=="ca") {
         opt.Tokenize(token, from, ")"); i1 = token.Atoi();
         if (i1<0 || i1>9) {
            printf("PaintSpectrum: Color Algorithm should be in the [0,9] range\n");
            i1 = 1;
         }
         sp.SetColorAlgorithm(i1);

      // Z Scale
      } else if (token=="zs") {
         opt.Tokenize(token, from, ")"); i1 = token.Atoi();
         if (i1<0 || i1>2) {
            printf("PaintSpectrum: Z-Scale should be in the [0,2] range\n");
            i1 = 0;
         }
         sp.SetZScale(i1);

      // Color Increment
      } else if (token=="ci") {
         opt.Tokenize(token, from, ","); f1 = token.Atof();
         opt.Tokenize(token, from, ","); f2 = token.Atof();
         opt.Tokenize(token, from, ")"); f3 = token.Atof();
         sp.SetColorIncrements(f1, f2, f3);

      // Light Height Weight
      } else if (token=="lhw") {
         opt.Tokenize(token, from, ")"); f1 = token.Atof();
         sp.SetLightHeightWeight(f1);

      // Light Position
      } else if (token=="lp") {
         opt.Tokenize(token, from, ","); i1 = token.Atoi();
         opt.Tokenize(token, from, ","); i2 = token.Atoi();
         opt.Tokenize(token, from, ")"); i3 = token.Atoi();
         sp.SetLightPosition(i1, i2, i3);

      // Contour Width
      } else if (token=="cw") {
         opt.Tokenize(token, from, ")"); i1 = token.Atoi();
         sp.SetContourWidth(i1);

      // Bezier
      } else if (token=="b") {
         opt.Tokenize(token, from, ")"); i1 = token.Atoi();
         if (i1<0 || i1>1) {
            printf("PaintSpectrum: Bezier should be in the [0,1] range\n");
            i1 = 0;
         }
         sp.SetBezier(i1);

      // Shading
      } else if (token=="s") {
         opt.Tokenize(token, from, ","); i1 = token.Atoi();
         if (i1<0 || i1>1) {
            printf("PaintSpectrum: Shading should be in the [0,1] range\n");
            i1 = 0;
         }
         opt.Tokenize(token, from, ")"); i2 = token.Atoi();
         if (i2<0 || i2>1) {
            printf("PaintSpectrum: Shadow should be in the [0,1] range\n");
            i2 = 0;
         }
         sp.SetShading(i1, i2);

      // Channel Marks
      } else if (token=="cm") {
         opt.Tokenize(token, from, ","); i1 = token.Atoi();
         opt.Tokenize(token, from, ","); i2 = token.Atoi();
         opt.Tokenize(token, from, ","); i3 = token.Atoi();
         opt.Tokenize(token, from, ","); i4 = token.Atoi();
         opt.Tokenize(token, from, ")"); i5 = token.Atoi();
         sp.SetChanMarks(i1, i2, i3, i4, i5);

      // Channel Grid
      } else if (token=="cg") {
         opt.Tokenize(token, from, ","); i1 = token.Atoi();
         opt.Tokenize(token, from, ")"); i2 = token.Atoi();
         sp.SetChanGrid(i1, i2);

      // Angles
      } else if (token=="a" || token=="a=") {
         opt.Tokenize(token, from, ","); i1 = token.Atoi();
         opt.Tokenize(token, from, ","); i2 = token.Atoi();
         opt.Tokenize(token, from, ")"); i3 = token.Atoi();
         sp.SetAngles(i1, i2, i3);

      // Buffer size
      } else if (token=="bf") {
         // Nothing to do here, The option "bf" has been handle before.
         // But it is a valid option.
         opt.Tokenize(token, from, ")");

      // Unknown option
      } else {
         if (!token.IsNull()) {
            printf("Unknown option \"%s\"\n",token.Data());
            return;
         }
      }
   }

   sp.Paint("");
}
