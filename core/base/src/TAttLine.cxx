// @(#)root/base:$Id$
// Author: Rene Brun   28/11/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TAttLine.h"
#include "TVirtualPad.h"
#include "TStyle.h"
#include "TVirtualX.h"
#include "TVirtualPadEditor.h"
#include "TColor.h"
#include <cmath>
#include <iostream>

ClassImp(TAttLine);
using std::sqrt;

/** \class TAttLine
\ingroup Base
\ingroup GraphicsAtt

Line Attributes class.

This class is used (in general by secondary inheritance)
by many other classes (graphics, histograms). It holds all the line attributes.

## Line attributes
Line attributes are:

  - [Line Color](#L1)
  - [Line Width](#L2)
  - [Line Style](#L3)

## <a name="L1"></a> Line Color
The line color is a color index (integer) pointing in the ROOT
color table.
The line color of any class inheriting from `TAttLine` can
be changed using the method `SetLineColor` and retrieved using the
method `GetLineColor`.
The following table shows the first 50 default colors.

Begin_Macro
{
   TCanvas *c = new TCanvas("c","Fill Area colors",0,0,500,200);
   c->DrawColorTable();
   return c;
}
End_Macro

### Color transparency
`SetLineColorAlpha()`, allows to set a transparent color.
In the following example the line color of the histogram `histo`
is set to blue with a transparency of 35%. The color `kBlue`
itself remains fully opaque.

~~~ {.cpp}
histo->SetLineColorAlpha(kBlue, 0.35);
~~~

The transparency is available on all platforms when the flag `OpenGL.CanvasPreferGL` is set to `1`
in `$ROOTSYS/etc/system.rootrc`, or on Mac with the Cocoa backend. On the file output
it is visible with PDF, PNG, Gif, JPEG, SVG, TeX ... but not PostScript.


## <a name="L2"></a> Line Width
The line width is expressed in pixel units.
The line width of any class inheriting from `TAttLine` can
be changed using the method `SetLineWidth` and retrieved using the
method `GetLineWidth`.
The following picture shows the line widths from 1 to 10 pixels.

Begin_Macro
{
   TCanvas *Lw = new TCanvas("Lw","test",500,200);
   TText  t;
   t.SetTextAlign(32);
   t.SetTextSize(0.08);
   Int_t i=1;
   for (float s=0.1; s<1.0 ; s+=0.092) {
      TLine *lh = new TLine(0.15,s,.85,s);
      lh->SetLineWidth(i);
      t.DrawText(0.1,s,Form("%d",i++));
      lh->Draw();
   }
}
End_Macro

## <a name="L3"></a> Line Style
Line styles are identified via integer numbers. The line style of any class
inheriting from `TAttLine` can be changed using the method
`SetLineStyle` and retrieved using the method `GetLineStyle`.

The first 10 line styles are predefined as shown on the following picture:

Begin_Macro
{
   TCanvas *Ls = new TCanvas("Ls","test",500,200);
   TText  t;
   t.SetTextAlign(32);
   t.SetTextSize(0.08);
   Int_t i=1;
   for (float s=0.1; s<1.0 ; s+=0.092) {
      TLine *lh = new TLine(0.15,s,.85,s);
      lh->SetLineStyle(i);
      lh->SetLineWidth(3);
      t.DrawText(0.1,s,Form("%d",i++));
      lh->Draw();
   }
}
End_Macro

Some line styles can be accessed via the following enum:

~~~ {.cpp}
   kSolid      = 1
   kDashed     = 2
   kDotted     = 3
   kDashDotted = 4
~~~

Additional line styles can be defined using `TStyle::SetLineStyleString`.
For example the line style number 11 can be defined as follow:
~~~ {.cpp}
   gStyle->SetLineStyleString(11,"400 200");
~~~
Existing line styles (1 to 10) can be redefined using the same method.
 */

////////////////////////////////////////////////////////////////////////////////
/// AttLine default constructor.

TAttLine::TAttLine()
{
   if (!gStyle) {fLineColor=1; fLineWidth=1; fLineStyle=1; return;}
   fLineColor = gStyle->GetLineColor();
   fLineWidth = gStyle->GetLineWidth();
   fLineStyle = gStyle->GetLineStyle();
}

////////////////////////////////////////////////////////////////////////////////
/// AttLine normal constructor.
/// Line attributes are taking from the argument list
///
///  - color : must be one of the valid color index
///  - style : 1=solid, 2=dash, 3=dash-dot, 4=dot-dot. New styles can be
///            defined using TStyle::SetLineStyleString.
///  - width : expressed in pixel units

TAttLine::TAttLine(Color_t color, Style_t style, Width_t width)
{
   fLineColor = color;
   fLineWidth = width;
   fLineStyle = style;
}

////////////////////////////////////////////////////////////////////////////////
/// AttLine destructor.

TAttLine::~TAttLine()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this line attributes to a new TAttLine.

void TAttLine::Copy(TAttLine &attline) const
{
   attline.fLineColor  = fLineColor;
   attline.fLineStyle  = fLineStyle;
   attline.fLineWidth  = fLineWidth;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from point px,py to a line.
/// Compute the closest distance of approach from point px,py to this line.
/// The distance is computed in pixels units.
///
/// Algorithm:
///~~~ {.cpp}
///   A(x1,y1)         P                             B(x2,y2)
///   -----------------+------------------------------
///                    |
///                    |
///                    |
///                    |
///                   M(x,y)
///
/// Let us call  a = distance AM     A=a**2
///              b = distance BM     B=b**2
///              c = distance AB     C=c**2
///              d = distance PM     D=d**2
///              u = distance AP     U=u**2
///              v = distance BP     V=v**2     c = u + v
///
/// D = A - U
/// D = B - V  = B -(c-u)**2
///    ==> u = (A -B +C)/2c
///~~~

Int_t TAttLine::DistancetoLine(Int_t px, Int_t py, Double_t xp1, Double_t yp1, Double_t xp2, Double_t yp2 )
{
   Double_t xl, xt, yl, yt;
   Double_t x     = px;
   Double_t y     = py;
   Double_t x1    = gPad->XtoAbsPixel(xp1);
   Double_t y1    = gPad->YtoAbsPixel(yp1);
   Double_t x2    = gPad->XtoAbsPixel(xp2);
   Double_t y2    = gPad->YtoAbsPixel(yp2);
   if (x1 < x2) {xl = x1; xt = x2;}
   else         {xl = x2; xt = x1;}
   if (y1 < y2) {yl = y1; yt = y2;}
   else         {yl = y2; yt = y1;}
   if (x < xl-2 || x> xt+2) return 9999;  //following algorithm only valid in the box
   if (y < yl-2 || y> yt+2) return 9999;  //surrounding the line
   Double_t xx1   = x  - x1;
   Double_t xx2   = x  - x2;
   Double_t x1x2  = x1 - x2;
   Double_t yy1   = y  - y1;
   Double_t yy2   = y  - y2;
   Double_t y1y2  = y1 - y2;
   Double_t a     = xx1*xx1   + yy1*yy1;
   Double_t b     = xx2*xx2   + yy2*yy2;
   Double_t c     = x1x2*x1x2 + y1y2*y1y2;
   if (c <= 0)  return 9999;
   Double_t v     = sqrt(c);
   Double_t u     = (a - b + c)/(2*v);
   Double_t d     = TMath::Abs(a - u*u);
   if (d < 0)   return 9999;

   return Int_t(sqrt(d) - 0.5*Double_t(fLineWidth));
}

////////////////////////////////////////////////////////////////////////////////
/// Change current line attributes if necessary.

void TAttLine::Modify()
{
   if (!gPad) return;
   Int_t lineWidth = TMath::Abs(fLineWidth%100);
   if (!gPad->IsBatch()) {
      gVirtualX->SetLineColor(fLineColor);
      if (fLineStyle > 0 && fLineStyle < 30) gVirtualX->SetLineStyle(fLineStyle);
      else                                   gVirtualX->SetLineStyle(1);
      gVirtualX->SetLineWidth(lineWidth);
   }

   if (fLineStyle > 0 && fLineStyle < 30) gPad->SetAttLinePS(fLineColor,fLineStyle,lineWidth);
   else                                   gPad->SetAttLinePS(fLineColor,1,lineWidth);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset this line attributes to default values.

void TAttLine::ResetAttLine(Option_t *)
{
   fLineColor  = 1;
   fLineStyle  = 1;
   fLineWidth  = 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Save line attributes as C++ statement(s) on output stream out.

void TAttLine::SaveLineAttributes(std::ostream &out, const char *name, Int_t coldef, Int_t stydef, Int_t widdef)
{
   if (fLineColor != coldef) {
      if (fLineColor > 228) {
         TColor::SaveColor(out, fLineColor);
         out<<"   "<<name<<"->SetLineColor(ci);" << std::endl;
      } else
         out<<"   "<<name<<"->SetLineColor("<<fLineColor<<");"<<std::endl;
   }
   if (fLineStyle != stydef) {
      out<<"   "<<name<<"->SetLineStyle("<<fLineStyle<<");"<<std::endl;
   }
   if (fLineWidth != widdef) {
      out<<"   "<<name<<"->SetLineWidth("<<fLineWidth<<");"<<std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Invoke the DialogCanvas Line attributes.

void TAttLine::SetLineAttributes()
{
   TVirtualPadEditor::UpdateLineAttributes(fLineColor,fLineStyle,fLineWidth);
}

////////////////////////////////////////////////////////////////////////////////
/// Set a transparent line color. lalpha defines the percentage of
/// the color opacity from 0. (fully transparent) to 1. (fully opaque).

void TAttLine::SetLineColorAlpha(Color_t lcolor, Float_t lalpha)
{
   fLineColor = TColor::GetColorTransparent(lcolor, lalpha);
}
