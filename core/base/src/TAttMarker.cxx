// @(#)root/base:$Id$
// Author: Rene Brun   12/05/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <iostream>
#include "Strlen.h"
#include "TAttMarker.h"
#include "TVirtualPad.h"
#include "TStyle.h"
#include "TVirtualX.h"
#include "TVirtualPadEditor.h"
#include "TColor.h"

ClassImp(TAttMarker);

/** \class TAttMarker
\ingroup Base
\ingroup GraphicsAtt

Marker Attributes class.

This class is used (in general by secondary inheritance)
by many other classes (graphics, histograms). It holds all the markers
attributes.

## Marker attributes
The marker attributes are:

  - [Marker color](#M1)
  - [Marker style](#M2)
    - [Marker line width](#M21)
  - [Marker size](#M3)

## <a name="M1"></a> Marker color
The marker color is a color index (integer) pointing in the ROOT color
table.
The marker color of any class inheriting from `TAttMarker` can
be changed using the method `SetMarkerColor` and retrieved using the
method `GetMarkerColor`.
The following table shows the first 50 default colors.

Begin_Macro
{
   TCanvas *c = new TCanvas("c","Fill Area colors",0,0,500,200);
   c->DrawColorTable();
   return c;
}
End_Macro

### Color transparency

`SetMarkerColorAlpha()`, allows to set a transparent color.
In the following example the marker color of the histogram `histo`
is set to blue with a transparency of 35%. The color `kBlue`
itself remains fully opaque.

~~~ {.cpp}
histo->SetMarkerColorAlpha(kBlue, 0.35);
~~~

The transparency is available on all platforms when the flag `OpenGL.CanvasPreferGL` is set to `1`
in `$ROOTSYS/etc/system.rootrc`, or on Mac with the Cocoa backend. On the file output
it is visible with PDF, PNG, Gif, JPEG, SVG, TeX ... but not PostScript.

## <a name="M2"></a> Marker style

The Marker style defines the markers' shape.
The marker style of any class inheriting from `TAttMarker` can
be changed using the method `SetMarkerStyle` and retrieved using the
method `GetMarkerStyle`.

The following list gives the currently supported markers (screen
and PostScript) style. Each marker style is identified by an integer number
(first column) corresponding to a marker shape (second column) and can be also
accessed via a global name (third column).

~~~ {.cpp}
   Marker number         Marker shape          Marker name
        1                    dot                  kDot
        2                    +                    kPlus
        3                    *                    kStar
        4                    o                    kCircle
        5                    x                    kMultiply
        6                    small dot            kFullDotSmall
        7                    medium dot           kFullDotMedium
        8                    large scalable dot   kFullDotLarge
        9 -->19              large scalable dot
       20                    full circle          kFullCircle
       21                    full square          kFullSquare
       22                    full triangle up     kFullTriangleUp
       23                    full triangle down   kFullTriangleDown
       24                    open circle          kOpenCircle
       25                    open square          kOpenSquare
       26                    open triangle up     kOpenTriangleUp
       27                    open diamond         kOpenDiamond
       28                    open cross           kOpenCross
       29                    full star            kFullStar
       30                    open star            kOpenStar
       31                    *
       32                    open triangle down   kOpenTriangleDown
       33                    full diamond         kFullDiamond
       34                    full cross           kFullCross
       35                    open diamond cross   kOpenDiamondCross
       36                    open square diagonal kOpenSquareDiagonal
       37                    open three triangle  kOpenThreeTriangles
       38                    octagon with cross   kOctagonCross
       39                    full three trangles  kFullThreeTriangles
       40                    open four triangleX  kOpenFourTrianglesX
       41                    full four triangleX  kFullFourTrianglesX
       42                    open double diamond  kOpenDoubleDiamond
       43                    full double diamond  kFullDoubleDiamond
       44                    open four triangle+  kOpenFourTrianglesPlus
       45                    full four triangle+  kFullFourTrianglesPlus
       46                    open cross X         kOpenCrossX
       47                    full cross X         kFullCrossX
       48                    four squares X       kFourSquaresX
       49                    four squares+        kFourSquaresPlus
~~~

Begin_Macro
{
   TCanvas *c = new TCanvas("c","Marker types",0,0,500,200);
   TMarker marker;
   marker.DisplayMarkerTypes();
}
End_Macro

### <a name="M21"></a> Marker line width

The line width of a marker is not actually a marker attribute since it does
only apply to open marker symbols and marker symbols consisting of lines. All
of these marker symbols are redefined with thicker lines by style numbers
starting from 50:

~~~ {.cpp}
   Marker numbers   Line width
      50 -  67         2
      68 -  85         3
      86 - 103         4
     104 - 121         5
   ...
~~~

Begin_Macro
{
   TCanvas *c = new TCanvas("c","Marker line widths",0,0,600,266);
   TMarker marker;
   marker.DisplayMarkerLineWidths();
}
End_Macro

## <a name="M3"></a> Marker size

Various marker sizes are shown in the figure below. The default marker size=1
is shown in the top left corner. Marker sizes smaller than 1 can be
specified. The marker size does not refer to any coordinate systems, it is an
absolute value. Therefore the marker size is not affected by any change
in TPad's scale. A marker size equal to 1 correspond to 8 pixels.
That is, a square marker with size 1 will be drawn with a side equal to 8
pixels on the screen.

The marker size of any class inheriting from `TAttMarker` can
be changed using the method `SetMarkerSize` and retrieved using the
method `GetMarkerSize`.

Begin_Macro
{
   c = new TCanvas("c","Marker sizes",0,0,500,200);
   TMarker marker;
   marker.SetMarkerStyle(3);
   Double_t x = 0;
   Double_t dx = 1/6.0;
   for (Int_t i=1; i<6; i++) {
      x += dx;
      marker.SetMarkerSize(i*0.2); marker.DrawMarker(x,.165);
      marker.SetMarkerSize(i*0.8); marker.DrawMarker(x,.495);
      marker.SetMarkerSize(i*1.0); marker.DrawMarker(x,.835);
   }
}
End_Macro

Note that the marker styles number 1 6 and 7 (the dots), cannot be scaled. They
are meant to be very fast to draw and are always drawn with the same number of
pixels; therefore `SetMarkerSize` does not apply on them. To have a
"scalable dot" a filled circle should be used instead, i.e. the marker style
number 20. By default (if `SetMarkerStyle` is not specified), the marker
style used is 1. That's the most common one to draw scatter plots.
*/

////////////////////////////////////////////////////////////////////////////////
/// TAttMarker default constructor.
///
/// Default text attributes are taking from the current style.

TAttMarker::TAttMarker()
{
   if (!gStyle) {fMarkerColor=1; fMarkerStyle=1; fMarkerSize=1; return;}
   fMarkerColor = gStyle->GetMarkerColor();
   fMarkerStyle = gStyle->GetMarkerStyle();
   fMarkerSize  = gStyle->GetMarkerSize();
}

////////////////////////////////////////////////////////////////////////////////
/// TAttMarker normal constructor.
///
/// Text attributes are taking from the argument list
///  - color : Marker Color Index
///  - style : Marker style (from 1 to 30)
///  - size  : marker size (float)

TAttMarker::TAttMarker(Color_t color, Style_t style, Size_t msize)
{
   fMarkerColor = color;
   fMarkerSize  = msize;
   fMarkerStyle = style;
}

////////////////////////////////////////////////////////////////////////////////
/// TAttMarker destructor.

TAttMarker::~TAttMarker()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this marker attributes to a new TAttMarker.

void TAttMarker::Copy(TAttMarker &attmarker) const
{
   attmarker.fMarkerColor  = fMarkerColor;
   attmarker.fMarkerStyle  = fMarkerStyle;
   attmarker.fMarkerSize   = fMarkerSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Internal helper function that returns the corresponding marker style with
/// line width 1 for the given style.

Style_t TAttMarker::GetMarkerStyleBase(Style_t style)
{
   if (style <= 49)
      return style;

   switch ((style - 50) % 18) {
   case 0:
      return 2;
   case 1:
      return 3;
   case 2:
      return 5;
   case 3:
      return 24;
   case 4:
      return 25;
   case 5:
      return 26;
   case 6:
      return 27;
   case 7:
      return 28;
   case 8:
      return 30;
   case 9:
      return 32;
   case 10:
      return 35;
   case 11:
      return 36;
   case 12:
      return 37;
   case 13:
      return 38;
   case 14:
      return 40;
   case 15:
      return 42;
   case 16:
      return 44;
   case 17:
      return 46;
   default:
      return style;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Internal helper function that returns the line width of the given marker
/// style (0 = filled marker)

Width_t TAttMarker::GetMarkerLineWidth(Style_t style)
{
   if (style >= 50)
      return ((style - 50) / 18) + 2;
   else if (style == 2 || style == 3 || style == 4 || style == 5
	    || style == 24 || style == 25 || style == 26 || style == 27
	    || style == 28 || style == 30 || style == 31 || style == 32
	    || style == 35 || style == 36 || style == 37 || style == 38
	    || style == 40 || style == 42 || style == 44 || style == 46)
      return 1;
   else
      return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Change current marker attributes if necessary.

void TAttMarker::Modify()
{
   if (!gPad) return;
   if (!gPad->IsBatch()) {
      gVirtualX->SetMarkerColor(fMarkerColor);
      gVirtualX->SetMarkerSize (fMarkerSize);
      gVirtualX->SetMarkerStyle(fMarkerStyle);
   }

   gPad->SetAttMarkerPS(fMarkerColor,fMarkerStyle,fMarkerSize);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset this marker attributes to the default values.

void TAttMarker::ResetAttMarker(Option_t *)
{
   fMarkerColor  = 1;
   fMarkerStyle  = 1;
   fMarkerSize   = 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Save line attributes as C++ statement(s) on output stream out.

void TAttMarker::SaveMarkerAttributes(std::ostream &out, const char *name, Int_t coldef, Int_t stydef, Int_t sizdef)
{
   if (fMarkerColor != coldef) {
      if (fMarkerColor > 228) {
         TColor::SaveColor(out, fMarkerColor);
         out<<"   "<<name<<"->SetMarkerColor(ci);" << std::endl;
      } else
         out<<"   "<<name<<"->SetMarkerColor("<<fMarkerColor<<");"<<std::endl;
   }
   if (fMarkerStyle != stydef) {
      out<<"   "<<name<<"->SetMarkerStyle("<<fMarkerStyle<<");"<<std::endl;
   }
   if (fMarkerSize != sizdef) {
      out<<"   "<<name<<"->SetMarkerSize("<<fMarkerSize<<");"<<std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Invoke the DialogCanvas Marker attributes.

void TAttMarker::SetMarkerAttributes()
{
   TVirtualPadEditor::UpdateMarkerAttributes(fMarkerColor,fMarkerStyle,fMarkerSize);
}

////////////////////////////////////////////////////////////////////////////////
/// Set a transparent marker color. malpha defines the percentage of
/// the color opacity from 0. (fully transparent) to 1. (fully opaque).

void TAttMarker::SetMarkerColorAlpha(Color_t mcolor, Float_t malpha)
{
   fMarkerColor = TColor::GetColorTransparent(mcolor, malpha);
}
