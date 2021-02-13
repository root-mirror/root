// @(#)root/base:$Id$
// Author: Rene Brun   12/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <iostream>
#include "TAttFill.h"
#include "TVirtualPad.h"
#include "TStyle.h"
#include "TVirtualX.h"
#include "TVirtualPadEditor.h"
#include "TColor.h"

ClassImp(TAttFill);

/** \class TAttFill
\ingroup Base
\ingroup GraphicsAtt

Fill Area Attributes class.

This class is used (in general by secondary inheritance)
by many other classes (graphics, histograms). It holds all the fill area
attributes.

## Fill Area attributes
Fill Area attributes are:

  - [Fill Area color](#F1)</a>
  - [Fill Area style](#F2)</a>

## <a name="F1"></a> Fill Area color
The fill area color is a color index (integer) pointing in the ROOT
color table.
The fill area color of any class inheriting from `TAttFill` can
be changed using the method `SetFillColor` and retrieved using the
method `GetFillColor`.
The following table shows the first 50 default colors.

Begin_Macro
{
   TCanvas *c = new TCanvas("c","Fill Area colors",0,0,500,200);
   c->DrawColorTable();
   return c;
}
End_Macro

### Color transparency
`SetFillColorAlpha()`, allows to set a transparent color.
In the following example the fill color of the histogram `histo`
is set to blue with a transparency of 35%. The color `kBlue`
itself remains fully opaque.

~~~ {.cpp}
histo->SetFillColorAlpha(kBlue, 0.35);
~~~

The transparency is available on all platforms when the flag
`OpenGL.CanvasPreferGL` is set to `1` in `$ROOTSYS/etc/system.rootrc`, or on Mac
with the Cocoa backend.
On the file output it is visible with PDF, PNG, Gif, JPEG, SVG, TeX... but not PostScript.

### The ROOT Color Wheel.
The wheel contains the recommended 216 colors to be used in web applications.
The colors in the Color Wheel are created by TColor::CreateColorWheel.

Using this color set for your text, background or graphics will give your
application a consistent appearance across different platforms and browsers.

Colors are grouped by hue, the aspect most important in human perception
Touching color chips have the same hue, but with different brightness and vividness.

Colors of slightly different hues _clash_. If you intend to display
colors of the same hue together, you should pick them from the same group.

Each color chip is identified by a mnemonic (eg kYellow) and a number.
The keywords, kRed, kBlue, kYellow, kPink, etc are defined in the header file __Rtypes.h__
that is included in all ROOT other header files. We strongly recommend to use these keywords
in your code instead of hardcoded color numbers, eg:
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

### Special case forcing black&white output.
If the current style fill area color is set to 0, then ROOT will force
a black&white output for all objects with a fill area defined and independently
of the object fill style.

## <a name="F2"></a> Fill Area style
The fill area style defines the pattern used to fill a polygon.
The fill area style of any class inheriting from `TAttFill` can
be changed using the method `SetFillStyle` and retrieved using the
method `GetFillStyle`.
### Conventions for fill styles:

  -   0    : hollow
  -   1001 : Solid
  -   3000+pattern_number (see below)
  -   For TPad only:

     -   4000 :the window is transparent.
     -   4000 to 4100 the window is 100% transparent to 100% opaque.

      The pad transparency is visible in binary outputs files like gif, jpg, png etc ..
      but not in vector graphics output files like PS, PDF and SVG. This convention
      (fill style > 4000) is kept for backward compatibility. It is better to use
      the color transparency instead.

pattern_number can have any value from 1 to 25 (see table), or any
value from 100 to 999. For the latest the numbering convention is the following:
~~~ {.cpp}
      pattern_number = ijk      (FillStyle = 3ijk)

      i (1-9) : specify the space between each hatch
                1 = 1/2mm  9 = 6mm

      j (0-9) : specify angle between 0 and 90 degrees
                0 = 0
                1 = 10
                2 = 20
                3 = 30
                4 = 45
                5 = Not drawn
                6 = 60
                7 = 70
                8 = 80
                9 = 90

      k (0-9) : specify angle between 90 and 180 degrees
                0 = 180
                1 = 170
                2 = 160
                3 = 150
                4 = 135
                5 = Not drawn
                6 = 120
                7 = 110
                8 = 100
                9 = 90
~~~
The following table shows the list of pattern styles.
The first table displays the 25 fixed patterns. They cannot be
customized unlike the hatches displayed in the second table which be
customized using:

  -  `gStyle->SetHatchesSpacing()` to define the spacing between hatches.
  -  `gStyle->SetHatchesLineWidth()` to define the hatches line width.

Begin_Macro
fillpatterns.C(500,700)
End_Macro
*/

////////////////////////////////////////////////////////////////////////////////
/// AttFill default constructor.
/// Default fill attributes are taking from the current style

TAttFill::TAttFill()
{
   if (!gStyle) {fFillColor=1; fFillStyle=0; return;}
   fFillColor = gStyle->GetFillColor();
   fFillStyle = gStyle->GetFillStyle();
}

////////////////////////////////////////////////////////////////////////////////
/// AttFill normal constructor.
///   - color Fill Color
///   - style Fill Style

TAttFill::TAttFill(Color_t color, Style_t style)
{
   fFillColor = color;
   fFillStyle = style;
}

////////////////////////////////////////////////////////////////////////////////
/// AttFill destructor.

TAttFill::~TAttFill()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this fill attributes to a new TAttFill.

void TAttFill::Copy(TAttFill &attfill) const
{
   attfill.fFillColor  = fFillColor;
   attfill.fFillStyle  = fFillStyle;
}

////////////////////////////////////////////////////////////////////////////////
/// Change current fill area attributes if necessary.

void TAttFill::Modify()
{
   if (!gPad) return;
   if (!gPad->IsBatch()) {
      gVirtualX->SetFillColor(fFillColor);
      gVirtualX->SetFillStyle(fFillStyle);
   }

   gPad->SetAttFillPS(fFillColor,fFillStyle);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset this fill attributes to default values.

void TAttFill::ResetAttFill(Option_t *)
{
   fFillColor = 1;
   fFillStyle = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Save fill attributes as C++ statement(s) on output stream out

void TAttFill::SaveFillAttributes(std::ostream &out, const char *name, Int_t coldef, Int_t stydef)
{
   if (fFillColor != coldef) {
      if (fFillColor > 228) {
         TColor::SaveColor(out, fFillColor);
         out<<"   "<<name<<"->SetFillColor(ci);" << std::endl;
      } else
         out<<"   "<<name<<"->SetFillColor("<<fFillColor<<");"<<std::endl;
   }
   if (fFillStyle != stydef) {
      out<<"   "<<name<<"->SetFillStyle("<<fFillStyle<<");"<<std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Invoke the DialogCanvas Fill attributes.

void TAttFill::SetFillAttributes()
{
   TVirtualPadEditor::UpdateFillAttributes(fFillColor,fFillStyle);
}

////////////////////////////////////////////////////////////////////////////////
/// Set a transparent fill color. falpha defines the percentage of
/// the color opacity from 0. (fully transparent) to 1. (fully opaque).

void TAttFill::SetFillColorAlpha(Color_t fcolor, Float_t falpha)
{
   fFillColor = TColor::GetColorTransparent(fcolor, falpha);
}
