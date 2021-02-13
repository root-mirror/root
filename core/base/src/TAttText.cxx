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
#include "Strlen.h"
#include "TROOT.h"
#include "TAttText.h"
#include "TVirtualPad.h"
#include "TStyle.h"
#include "TVirtualX.h"
#include "TError.h"
#include "TVirtualPadEditor.h"
#include "TColor.h"

ClassImp(TAttText);

/** \class TAttText
\ingroup Base
\ingroup GraphicsAtt

Text Attributes class.

This class is used (in general by secondary inheritance)
by many other classes (graphics, histograms). It holds all the text attributes.

## Text attributes
Text attributes are:

  - [Text Alignment](#T1)
  - [Text Angle](#T2)
  - [Text Color](#T3)
  - [Text Size](#T4)
  - [Text Font and Precision](#T5)
     - [Font quality and speed](#T51)
     - [How to use True Type Fonts](#T52)
     - [List of the currently supported fonts](#T53)

## <a name="T1"></a> Text Alignment

The text alignment is an integer number (`align`) allowing to control
the horizontal and vertical position of the text string with respect
to the text position.
The text alignment of any class inheriting from `TAttText` can
be changed using the method `SetTextAlign` and retrieved using the
method `GetTextAlign`.

~~~ {.cpp}
   align = 10*HorizontalAlign + VerticalAlign
~~~

For horizontal alignment the following convention applies:

~~~ {.cpp}
   1=left adjusted, 2=centered, 3=right adjusted
~~~

For vertical alignment the following convention applies:

~~~ {.cpp}
   1=bottom adjusted, 2=centered, 3=top adjusted
~~~

For example:

~~~ {.cpp}
   align = 11 = left adjusted and bottom adjusted
   align = 32 = right adjusted and vertically centered
~~~

Begin_Macro(source)
textalign.C
End_Macro

Mnemonic constants are available:

~~~ {.cpp}
kHAlignLeft   = 10, kHAlignCenter = 20, kHAlignRight = 30,
kVAlignBottom = 1,  kVAlignCenter = 2,  kVAlignTop   = 3
~~~

They allow to write:

~~~ {.cpp}
object->SetTextAlign(kHAlignLeft+kVAlignTop);
~~~

## <a name="T2"></a> Text Angle

Text angle in degrees.
The text angle of any class inheriting from `TAttText` can
be changed using the method `SetTextAngle` and retrieved using the
method `GetTextAngle`.
The following picture shows the text angle:

Begin_Macro(source)
textangle.C
End_Macro

## <a name="T3"></a> Text Color

The text color is a color index (integer) pointing in the ROOT
color table.
The text color of any class inheriting from `TAttText` can
be changed using the method `SetTextColor` and retrieved using the
method `GetTextColor`.
The following table shows the first 50 default colors.

Begin_Macro
{
   TCanvas *c = new TCanvas("c","Fill Area colors",0,0,500,200);
   c->DrawColorTable();
   return c;
}
End_Macro

### Color transparency
`SetTextColorAlpha()`, allows to set a transparent color.
In the following example the text color of the text `text`
is set to blue with a transparency of 35%. The color `kBlue`
itself remains fully opaque.

~~~ {.cpp}
text->SetTextColorAlpha(kBlue, 0.35);
~~~

The transparency is available on all platforms when the flag `OpenGL.CanvasPreferGL` is set to `1`
in `$ROOTSYS/etc/system.rootrc`, or on Mac with the Cocoa backend. On the file output
it is visible with PDF, PNG, Gif, JPEG, SVG, TeX ... but not PostScript.

## <a name="T4"></a> Text Size

If the text precision (see next paragraph) is smaller than 3, the text
size (`textsize`) is a fraction of the current pad size. Therefore the
same `textsize` value can generate text outputs with different absolute
sizes in two different pads.
The text size in pixels (`charheight`) is computed the following way:

~~~ {.cpp}
   pad_width  = gPad->XtoPixel(gPad->GetX2());
   pad_height = gPad->YtoPixel(gPad->GetY1());
   if (pad_width < pad_height)  charheight = textsize*pad_width;
   else                         charheight = textsize*pad_height;
~~~

If the text precision is equal to 3, the text size doesn't depend on the pad's
dimensions. A given `textsize` value always generates the same absolute
size. The text size (`charheight`) is given in pixels:

~~~ {.cpp}
   charheight = textsize;
~~~

Note that to scale fonts to the same size as the old True Type package a
scale factor of `0.93376068` is apply to the text size before drawing.

The text size of any class inheriting from `TAttText` can
be changed using the method `SetTextSize` and retrieved using the
method `GetTextSize`.

## <a name="T5"></a> Text Font and Precision

The text font code is combination of the font number and the precision.
~~~ {.cpp}
   Text font code = 10*fontnumber + precision
~~~
Font numbers must be between 1 and 14.

The precision can be:

  - `precision = 0` fast hardware fonts (steps in the size)
  - `precision = 1` scalable and rotatable hardware fonts (see below)
  - `precision = 2` scalable and rotatable hardware fonts
  - `precision = 3` scalable and rotatable hardware fonts. Text size
                           is given in pixels.

The text font and precision of any class inheriting from `TAttText` can
be changed using the method `SetTextFont` and retrieved using the
method `GetTextFont`.

### <a name="T51"></a> Font quality and speed

When precision 0 is used, only the original non-scaled system fonts are
used. The fonts have a minimum (4) and maximum (37) size in pixels. These
fonts are fast and are of good quality. Their size varies with large steps
and they cannot be rotated.
Precision 1 and 2 fonts have a different behaviour depending if the
True Type Fonts (TTF) are used or not. If TTF are used, you always get very good
quality scalable and rotatable fonts. However TTF are slow.

### <a name="T52"></a> How to use True Type Fonts

One can activate the TTF by adding (or activating) the following line
in the `.rootrc` file:

~~~ {.cpp}
   Unix.*.Root.UseTTFonts:     true
~~~

It is possible to check the TTF are in use in a Root session
with the command:

~~~ {.cpp}
   gEnv->Print();
~~~

If the TTF are in use the following line will appear at the beginning of the
printout given by this command:

~~~ {.cpp}
   Unix.*.Root.UseTTFonts:   true                           [Global]
~~~

### <a name="T53"></a> List of the currently supported fonts

~~~ {.cpp}
   Font number         X11 Names             Win32/TTF Names
       1 :       times-medium-i-normal      "Times New Roman"
       2 :       times-bold-r-normal        "Times New Roman"
       3 :       times-bold-i-normal        "Times New Roman"
       4 :       helvetica-medium-r-normal  "Arial"
       5 :       helvetica-medium-o-normal  "Arial"
       6 :       helvetica-bold-r-normal    "Arial"
       7 :       helvetica-bold-o-normal    "Arial"
       8 :       courier-medium-r-normal    "Courier New"
       9 :       courier-medium-o-normal    "Courier New"
      10 :       courier-bold-r-normal      "Courier New"
      11 :       courier-bold-o-normal      "Courier New"
      12 :       symbol-medium-r-normal     "Symbol"
      13 :       times-medium-r-normal      "Times New Roman"
      14 :                                  "Wingdings"
      15 :       Symbol italic (derived from Symbol)
~~~

The following picture shows how each font looks. The number on the left
is the "text font code". In this picture precision 2 was selected.

Begin_Macro
fonts.C
End_Macro
*/

////////////////////////////////////////////////////////////////////////////////
/// AttText default constructor.
///
/// Default text attributes are taken from the current style.

TAttText::TAttText()
{
   if (!gStyle) {
      ResetAttText();
      return;
   }
   fTextAlign = gStyle->GetTextAlign();
   fTextAngle = gStyle->GetTextAngle();
   fTextColor = gStyle->GetTextColor();
   fTextFont  = gStyle->GetTextFont();
   fTextSize  = gStyle->GetTextSize();
}

////////////////////////////////////////////////////////////////////////////////
/// AttText normal constructor.
///
/// Text attributes are taken from the argument list.

TAttText::TAttText(Int_t align, Float_t angle, Color_t color, Style_t font, Float_t tsize)
{
   fTextAlign = align;
   fTextAngle = angle;
   fTextColor = color;
   fTextFont  = font;
   fTextSize  = tsize;
}

////////////////////////////////////////////////////////////////////////////////
/// AttText destructor.

TAttText::~TAttText()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this text attributes to a new TAttText.

void TAttText::Copy(TAttText &atttext) const
{
   atttext.fTextAlign  = fTextAlign;
   atttext.fTextAngle  = fTextAngle;
   atttext.fTextColor  = fTextColor;
   atttext.fTextFont   = fTextFont;
   atttext.fTextSize   = fTextSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Change current text attributes if necessary.

void TAttText::Modify()
{
   if (!gPad) return;

   // Do we need to change font?
   if (!gPad->IsBatch()) {
      gVirtualX->SetTextAngle(fTextAngle);
      Float_t wh = (Float_t)gPad->XtoPixel(gPad->GetX2());
      Float_t hh = (Float_t)gPad->YtoPixel(gPad->GetY1());
      Float_t tsize;
      if (wh < hh)  tsize = fTextSize*wh;
      else          tsize = fTextSize*hh;
      if (fTextFont%10 > 2) tsize = fTextSize;

      if (gVirtualX->GetTextFont() != fTextFont) {
            gVirtualX->SetTextFont(fTextFont);
            gVirtualX->SetTextSize(tsize);
      }
      if (gVirtualX->GetTextSize() != tsize)
            gVirtualX->SetTextSize(tsize);
      gVirtualX->SetTextAlign(fTextAlign);
      gVirtualX->SetTextColor(fTextColor);
   }
   gPad->SetAttTextPS(fTextAlign,fTextAngle,fTextColor,fTextFont,fTextSize);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset this text attributes to default values.

void TAttText::ResetAttText(Option_t *)
{
   fTextAlign  = 11;
   fTextAngle  = 0;
   fTextColor  = 1;
   fTextFont   = 62;
   fTextSize   = 0.05;
}

////////////////////////////////////////////////////////////////////////////////
/// Save text attributes as C++ statement(s) on output stream out.

void TAttText::SaveTextAttributes(std::ostream &out, const char *name, Int_t alidef,
                                  Float_t angdef, Int_t coldef, Int_t fondef,
                                  Float_t sizdef)
{
   if (fTextAlign != alidef) {
      out<<"   "<<name<<"->SetTextAlign("<<fTextAlign<<");"<<std::endl;
   }
   if (fTextColor != coldef) {
      if (fTextColor > 228) {
         TColor::SaveColor(out, fTextColor);
         out<<"   "<<name<<"->SetTextColor(ci);" << std::endl;
      } else
         out<<"   "<<name<<"->SetTextColor("<<fTextColor<<");"<<std::endl;
   }
   if (fTextFont != fondef) {
      out<<"   "<<name<<"->SetTextFont("<<fTextFont<<");"<<std::endl;
   }
   if (fTextSize != sizdef) {
      out<<"   "<<name<<"->SetTextSize("<<fTextSize<<");"<<std::endl;
   }
   if (fTextAngle != angdef) {
      out<<"   "<<name<<"->SetTextAngle("<<fTextAngle<<");"<<std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Invoke the DialogCanvas Text attributes.

void TAttText::SetTextAttributes()
{
   TVirtualPadEditor::UpdateTextAttributes(fTextAlign,fTextAngle,fTextColor,
                                           fTextFont,fTextSize);
}

////////////////////////////////////////////////////////////////////////////////
/// Set a transparent marker color. talpha defines the percentage of
/// the color opacity from 0. (fully transparent) to 1. (fully opaque).

void TAttText::SetTextColorAlpha(Color_t tcolor, Float_t talpha)
{
   fTextColor = TColor::GetColorTransparent(tcolor, talpha);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the text size in pixels.
///
/// If the font precision is greater than 2, the text size is set to npixels,
/// otherwise the text size is computed as a percent of the pad size.

void TAttText::SetTextSizePixels(Int_t npixels)
{
   if (fTextFont%10 > 2) {
      fTextSize = Float_t(npixels);
   } else {
      TVirtualPad *pad = gROOT->GetSelectedPad();
      if (!pad) return;
      Float_t dy = pad->AbsPixeltoY(0) - pad->AbsPixeltoY(npixels);
      fTextSize = dy/(pad->GetY2() - pad->GetY1());
   }
}
