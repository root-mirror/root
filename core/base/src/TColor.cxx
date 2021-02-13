// @(#)root/base:$Id$
// Author: Rene Brun   12/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"
#include "TColor.h"
#include "TObjArray.h"
#include "TArrayI.h"
#include "TArrayD.h"
#include "TVirtualX.h"
#include "TError.h"
#include "TMathBase.h"
#include "TApplication.h"
#include "snprintf.h"
#include <algorithm>
#include <cmath>
#include <iostream>

ClassImp(TColor);

namespace {
   static Bool_t& TColor__GrayScaleMode() {
      static Bool_t grayScaleMode;
      return grayScaleMode;
   }
   static TArrayI& TColor__Palette() {
      static TArrayI globalPalette(0);
      return globalPalette;
   }
   static TArrayD& TColor__PalettesList() {
      static TArrayD globalPalettesList(0);
      return globalPalettesList;
   }
}

static Int_t   gHighestColorIndex = 0;   ///< Highest color index defined
static Float_t gColorThreshold    = -1.; ///< Color threshold used by GetColor
static Int_t   gDefinedColors     = 0;   ///< Number of defined colors.
static Int_t   gLastDefinedColors = 649; ///< Previous number of defined colors

#define fgGrayscaleMode TColor__GrayScaleMode()
#define fgPalette TColor__Palette()
#define fgPalettesList TColor__PalettesList()

using std::floor;

/** \class TColor
\ingroup Base
\ingroup GraphicsAtt

The color creation and management class.

  - [Introduction](#C00)
  - [Basic colors](#C01)
  - [The color wheel](#C02)
  - [Bright and dark colors](#C03)
  - [Gray scale view of of canvas with colors](#C04)
  - [Color palettes](#C05)
  - [High quality predefined palettes](#C06)
  - [Palette inversion](#C061)
  - [Color transparency](#C07)

## <a name="C00"></a> Introduction

Colors are defined by their red, green and blue components, simply called the
RGB components. The colors are also known by the hue, light and saturation
components also known as the HLS components. When a new color is created the
components of both color systems are computed.

At initialization time, a table of colors is generated. An existing color can
be retrieved by its index:

~~~ {.cpp}
   TColor *color = gROOT->GetColor(10);
~~~

Then it can be manipulated. For example its RGB components can be modified:

~~~ {.cpp}
   color->SetRGB(0.1, 0.2, 0.3);
~~~

A new color can be created the following way:

~~~ {.cpp}
   Int_t ci = 1756; // color index
   TColor *color = new TColor(ci, 0.1, 0.2, 0.3);
~~~

\since **6.07/07:**
TColor::GetFreeColorIndex() allows to make sure the new color is created with an
unused color index:

~~~ {.cpp}
   Int_t ci = TColor::GetFreeColorIndex();
   TColor *color = new TColor(ci, 0.1, 0.2, 0.3);
~~~

Two sets of colors are initialized;

  -  The basic colors: colors with index from 0 to 50.
  -  The color wheel: colors with indices from 300 to 1000.

## <a name="C01"></a> Basic colors
The following image displays the 50 basic colors.

Begin_Macro(source)
{
   TCanvas *c = new TCanvas("c","Fill Area colors",0,0,500,200);
   c->DrawColorTable();
   return c;
}
End_Macro

## <a name="C02"></a> The color wheel
The wheel contains the recommended 216 colors to be used in web applications.

The colors in the color wheel are created by `TColor::CreateColorWheel`.

Using this color set for your text, background or graphics will give your
application a consistent appearance across different platforms and browsers.

Colors are grouped by hue, the aspect most important in human perception.
Touching color chips have the same hue, but with different brightness and
vividness.

Colors of slightly different hues clash. If you intend to display
colors of the same hue together, you should pick them from the same group.

Each color chip is identified by a mnemonic (e.g. kYellow) and a number.
The keywords, kRed, kBlue, kYellow, kPink, etc are defined in the header file
Rtypes.h that is included in all ROOT other header files. It is better
to use these keywords in user code instead of hardcoded color numbers, e.g.:

~~~ {.cpp}
   myObject.SetFillColor(kRed);
   myObject.SetFillColor(kYellow-10);
   myLine.SetLineColor(kMagenta+2);
~~~

Begin_Macro(source)
{
   TColorWheel *w = new TColorWheel();
   cw = new TCanvas("cw","cw",0,0,400,400);
   w->SetCanvas(cw);
   w->Draw();
}
End_Macro

The complete list of predefined color names is the following:

~~~ {.cpp}
kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900
~~~

Note the special role of color `kWhite` (color number 0). It is the default
background color also. For instance in a PDF or PS files (as paper is usually white)
it is simply not painted. To have a white color behaving like the other color the
simplest is to define an other white color not attached to the color index 0:

~~~ {.cpp}
   Int_t ci = TColor::GetFreeColorIndex();
   TColor *color = new TColor(ci, 1., 1., 1.);
~~~

## <a name="C03"></a> Bright and dark colors
The dark and bright color are used to give 3-D effects when drawing various
boxes (see TWbox, TPave, TPaveText, TPaveLabel, etc).

  - The dark colors have an index = color_index+100
  - The bright colors have an index = color_index+150
  - Two static functions return the bright and dark color number
    corresponding to a color index. If the bright or dark color does not
    exist, they are created:
   ~~~ {.cpp}
      Int_t dark   = TColor::GetColorDark(color_index);
      Int_t bright = TColor::GetColorBright(color_index);
   ~~~

## <a name="C04"></a> Grayscale view of of canvas with colors
One can toggle between a grayscale preview and the regular colored mode using
`TCanvas::SetGrayscale()`. Note that in grayscale mode, access via RGB
will return grayscale values according to ITU standards (and close to b&w
printer gray-scales), while access via HLS returns de-saturated gray-scales. The
image below shows the ROOT color wheel in grayscale mode.

Begin_Macro(source)
{
   TColorWheel *w = new TColorWheel();
   cw = new TCanvas("cw","cw",0,0,400,400);
   cw->GetCanvas()->SetGrayscale();
   w->SetCanvas(cw);
   w->Draw();
}
End_Macro

## <a name="C05"></a> Color palettes
It is often very useful to represent a variable with a color map. The concept
of "color palette" allows to do that. One color palette is active at any time.
This "current palette" is set using:

~~~ {.cpp}
gStyle->SetPalette(...);
~~~

This function has two parameters: the number of colors in the palette and an
array of containing the indices of colors in the palette. The following small
example demonstrates how to define and use the color palette:

Begin_Macro(source)
{
   TCanvas *c1  = new TCanvas("c1","c1",0,0,600,400);
   TF2 *f1 = new TF2("f1","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",1,3,1,3);
   Int_t palette[5];
   palette[0] = 15;
   palette[1] = 20;
   palette[2] = 23;
   palette[3] = 30;
   palette[4] = 32;
   gStyle->SetPalette(5,palette);
   f1->Draw("colz");
   return c1;
}
End_Macro

 To define more a complex palette with a continuous gradient of color, one
should use the static function `TColor::CreateGradientColorTable()`.
The following example demonstrates how to proceed:

Begin_Macro(source)
{
   TCanvas *c2  = new TCanvas("c2","c2",0,0,600,400);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",1,3,1,3);
   const Int_t Number = 3;
   Double_t Red[Number]    = { 1.00, 0.00, 0.00};
   Double_t Green[Number]  = { 0.00, 1.00, 0.00};
   Double_t Blue[Number]   = { 1.00, 0.00, 1.00};
   Double_t Length[Number] = { 0.00, 0.50, 1.00 };
   Int_t nb=50;
   TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
   f2->SetContour(nb);
   f2->SetLineWidth(1);
   f2->SetLineColor(kBlack);
   f2->Draw("surf1z");
   return c2;
}
End_Macro

The function `TColor::CreateGradientColorTable()` automatically
calls `gStyle->SetPalette()`, so there is not need to add one.

After a call to `TColor::CreateGradientColorTable()` it is sometimes
useful to store the newly create palette for further use. In particular, it is
recommended to do if one wants to switch between several user define palettes.
To store a palette in an array it is enough to do:

~~~ {.cpp}
   Int_t MyPalette[100];
   Double_t Red[]    = {0., 0.0, 1.0, 1.0, 1.0};
   Double_t Green[]  = {0., 0.0, 0.0, 1.0, 1.0};
   Double_t Blue[]   = {0., 1.0, 0.0, 0.0, 1.0};
   Double_t Length[] = {0., .25, .50, .75, 1.0};
   Int_t FI = TColor::CreateGradientColorTable(5, Length, Red, Green, Blue, 100);
   for (int i=0;i<100;i++) MyPalette[i] = FI+i;
~~~

Later on to reuse the palette `MyPalette` it will be enough to do

~~~ {.cpp}
   gStyle->SetPalette(100, MyPalette);
~~~

As only one palette is active, one need to use `TExec` to be able to
display plots using different palettes on the same pad.
The tutorial multipalette.C illustrates this feature.

Begin_Macro(source)
../../../tutorials/graphs/multipalette.C
End_Macro

## <a name="C06"></a> High quality predefined palettes
\since **6.04:**
62 high quality palettes are predefined with 255 colors each.
Despite the [disadvantages of the Rainbow color map](https://root.cern.ch/rainbow-color-map),
it was kept in the list of predefined color maps.
These palettes can be accessed "by name" with `gStyle->SetPalette(num)`.
`num` can be taken within the following enum:

~~~ {.cpp}
kDeepSea=51,          kGreyScale=52,    kDarkBodyRadiator=53,
kBlueYellow= 54,      kRainBow=55,      kInvertedDarkBodyRadiator=56,
kBird=57,             kCubehelix=58,    kGreenRedViolet=59,
kBlueRedYellow=60,    kOcean=61,        kColorPrintableOnGrey=62,
kAlpine=63,           kAquamarine=64,   kArmy=65,
kAtlantic=66,         kAurora=67,       kAvocado=68,
kBeach=69,            kBlackBody=70,    kBlueGreenYellow=71,
kBrownCyan=72,        kCMYK=73,         kCandy=74,
kCherry=75,           kCoffee=76,       kDarkRainBow=77,
kDarkTerrain=78,      kFall=79,         kFruitPunch=80,
kFuchsia=81,          kGreyYellow=82,   kGreenBrownTerrain=83,
kGreenPink=84,        kIsland=85,       kLake=86,
kLightTemperature=87, kLightTerrain=88, kMint=89,
kNeon=90,             kPastel=91,       kPearl=92,
kPigeon=93,           kPlum=94,         kRedBlue=95,
kRose=96,             kRust=97,         kSandyTerrain=98,
kSienna=99,           kSolar=100,       kSouthWest=101,
kStarryNight=102,     kSunset=103,      kTemperatureMap=104,
kThermometer=105,     kValentine=106,   kVisibleSpectrum=107,
kWaterMelon=108,      kCool=109,        kCopper=110,
kGistEarth=111,       kViridis=112,     kCividis=113
~~~

<table border=0>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kDeepSea);
   f2->Draw("surf2Z"); f2->SetTitle("kDeepSea");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kGreyScale);
   f2->Draw("surf2Z"); f2->SetTitle("kGreyScale");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kDarkBodyRadiator);
   f2->Draw("surf2Z"); f2->SetTitle("kDarkBodyRadiator");
}
End_Macro
</td></tr>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kBlueYellow);
   f2->Draw("surf2Z"); f2->SetTitle("kBlueYellow");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kRainBow);
   f2->Draw("surf2Z"); f2->SetTitle("kRainBow");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kInvertedDarkBodyRadiator);
   f2->Draw("surf2Z"); f2->SetTitle("kInvertedDarkBodyRadiator");
}
End_Macro
</td></tr>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kBird);
   f2->Draw("surf2Z"); f2->SetTitle("kBird (default)");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kCubehelix);
   f2->Draw("surf2Z"); f2->SetTitle("kCubehelix");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kGreenRedViolet);
   f2->Draw("surf2Z"); f2->SetTitle("kGreenRedViolet");
}
End_Macro
</td></tr>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kBlueRedYellow);
   f2->Draw("surf2Z"); f2->SetTitle("kBlueRedYellow");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kOcean);
   f2->Draw("surf2Z"); f2->SetTitle("kOcean");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kColorPrintableOnGrey);
   f2->Draw("surf2Z"); f2->SetTitle("kColorPrintableOnGrey");
}
End_Macro
</td></tr>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kAlpine);
   f2->Draw("surf2Z"); f2->SetTitle("kAlpine");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kAquamarine);
   f2->Draw("surf2Z"); f2->SetTitle("kAquamarine");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kArmy);
   f2->Draw("surf2Z"); f2->SetTitle("kArmy");
}
End_Macro
</td></tr>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kAtlantic);
   f2->Draw("surf2Z"); f2->SetTitle("kAtlantic");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kAurora);
   f2->Draw("surf2Z"); f2->SetTitle("kAurora");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kAvocado);
   f2->Draw("surf2Z"); f2->SetTitle("kAvocado");
}
End_Macro
</td></tr>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kBeach);
   f2->Draw("surf2Z"); f2->SetTitle("kBeach");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kBlackBody);
   f2->Draw("surf2Z"); f2->SetTitle("kBlackBody");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kBlueGreenYellow);
   f2->Draw("surf2Z"); f2->SetTitle("kBlueGreenYellow");
}
End_Macro
</td></tr>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kBrownCyan);
   f2->Draw("surf2Z"); f2->SetTitle("kBrownCyan");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kCMYK);
   f2->Draw("surf2Z"); f2->SetTitle("kCMYK");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kCandy);
   f2->Draw("surf2Z"); f2->SetTitle("kCandy");
}
End_Macro
</td></tr>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kCherry);
   f2->Draw("surf2Z"); f2->SetTitle("kCherry");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kCoffee);
   f2->Draw("surf2Z"); f2->SetTitle("kCoffee");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kDarkRainBow);
   f2->Draw("surf2Z"); f2->SetTitle("kDarkRainBow");
}
End_Macro
</td></tr>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kDarkTerrain);
   f2->Draw("surf2Z"); f2->SetTitle("kDarkTerrain");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kFall);
   f2->Draw("surf2Z"); f2->SetTitle("kFall");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kFruitPunch);
   f2->Draw("surf2Z"); f2->SetTitle("kFruitPunch");
}
End_Macro
</td></tr>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kFuchsia);
   f2->Draw("surf2Z"); f2->SetTitle("kFuchsia");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kGreyYellow);
   f2->Draw("surf2Z"); f2->SetTitle("kGreyYellow");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kGreenBrownTerrain);
   f2->Draw("surf2Z"); f2->SetTitle("kGreenBrownTerrain");
}
End_Macro
</td></tr>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kGreenPink);
   f2->Draw("surf2Z"); f2->SetTitle("kGreenPink");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kIsland);
   f2->Draw("surf2Z"); f2->SetTitle("kIsland");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kLake);
   f2->Draw("surf2Z"); f2->SetTitle("kLake");
}
End_Macro
</td></tr>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kLightTemperature);
   f2->Draw("surf2Z"); f2->SetTitle("kLightTemperature");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kLightTerrain);
   f2->Draw("surf2Z"); f2->SetTitle("kLightTerrain");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kMint);
   f2->Draw("surf2Z"); f2->SetTitle("kMint");
}
End_Macro
</td></tr>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kNeon);
   f2->Draw("surf2Z"); f2->SetTitle("kNeon");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kPastel);
   f2->Draw("surf2Z"); f2->SetTitle("kPastel");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kPearl);
   f2->Draw("surf2Z"); f2->SetTitle("kPearl");
}
End_Macro
</td></tr>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kPigeon);
   f2->Draw("surf2Z"); f2->SetTitle("kPigeon");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kPlum);
   f2->Draw("surf2Z"); f2->SetTitle("kPlum");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kRedBlue);
   f2->Draw("surf2Z"); f2->SetTitle("kRedBlue");
}
End_Macro
</td></tr>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kRose);
   f2->Draw("surf2Z"); f2->SetTitle("kRose");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kRust);
   f2->Draw("surf2Z"); f2->SetTitle("kRust");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kSandyTerrain);
   f2->Draw("surf2Z"); f2->SetTitle("kSandyTerrain");
}
End_Macro
</td></tr>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kSienna);
   f2->Draw("surf2Z"); f2->SetTitle("kSienna");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kSolar);
   f2->Draw("surf2Z"); f2->SetTitle("kSolar");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kSouthWest);
   f2->Draw("surf2Z"); f2->SetTitle("kSouthWest");
}
End_Macro
</td></tr>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kStarryNight);
   f2->Draw("surf2Z"); f2->SetTitle("kStarryNight");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kSunset);
   f2->Draw("surf2Z"); f2->SetTitle("kSunset");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kTemperatureMap);
   f2->Draw("surf2Z"); f2->SetTitle("kTemperatureMap");
}
End_Macro
</td></tr>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kThermometer);
   f2->Draw("surf2Z"); f2->SetTitle("kThermometer");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kValentine);
   f2->Draw("surf2Z"); f2->SetTitle("kValentine");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kVisibleSpectrum);
   f2->Draw("surf2Z"); f2->SetTitle("kVisibleSpectrum");
}
End_Macro
</td></tr>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kWaterMelon);
   f2->Draw("surf2Z"); f2->SetTitle("kWaterMelon");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kCool);
   f2->Draw("surf2Z"); f2->SetTitle("kCool");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kCopper);
   f2->Draw("surf2Z"); f2->SetTitle("kCopper");
}
End_Macro
</td></tr>
<tr><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kGistEarth);
   f2->Draw("surf2Z"); f2->SetTitle("kGistEarth");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kViridis);
   f2->Draw("surf2Z"); f2->SetTitle("kViridis");
}
End_Macro
</td><td>
Begin_Macro
{
   c  = new TCanvas("c","c",0,0,300,300);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kCividis);
   f2->Draw("surf2Z"); f2->SetTitle("kCividis");
}
End_Macro
</td></tr>
</table>

## <a name="C061"></a> Palette inversion
Once a palette is defined, it is possible to invert the color order thanks to the
method TColor::InvertPalette. The top of the palette becomes the bottom and vice versa.

Begin_Macro(source)
{
   auto c  = new TCanvas("c","c",0,0,600,400);
   TF2 *f2 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",0.999,3.002,0.999,3.002);
   f2->SetContour(99); gStyle->SetPalette(kCherry);
   TColor::InvertPalette();
   f2->Draw("surf2Z"); f2->SetTitle("kCherry inverted");
}
End_Macro

## <a name="C07"></a> Color transparency
To make a graphics object transparent it is enough to set its color to a
transparent one. The color transparency is defined via its alpha component. The
alpha value varies from `0.` (fully transparent) to `1.` (fully
opaque). To set the alpha value of an existing color it is enough to do:

~~~ {.cpp}
   TColor *col26 = gROOT->GetColor(26);
   col26->SetAlpha(0.01);
~~~

A new color can be created transparent the following way:

~~~ {.cpp}
   Int_t ci = 1756;
   TColor *color = new TColor(ci, 0.1, 0.2, 0.3, "", 0.5); // alpha = 0.5
~~~

An example of transparency usage with parallel coordinates can be found
in parallelcoordtrans.C.

To ease the creation of a transparent color the static method
`GetColorTransparent(Int_t color, Float_t a)` is provided.
In the following example the `trans_red` color index point to
a red color 30% transparent. The alpha value of the color index
`kRed` is not modified.

~~~ {.cpp}
   Int_t trans_red = GetColorTransparent(kRed, 0.3);
~~~

This function is also used in the methods
`SetFillColorAlpha()`, `SetLineColorAlpha()`,
`SetMarkerColorAlpha()` and `SetTextColorAlpha()`.
In the following example the fill color of the histogram `histo`
is set to blue with a transparency of 35%. The color `kBlue`
itself remains fully opaque.

~~~ {.cpp}
   histo->SetFillColorAlpha(kBlue, 0.35);
~~~

The transparency is available on all platforms when the flag `OpenGL.CanvasPreferGL` is set to `1`
in `$ROOTSYS/etc/system.rootrc`, or on Mac with the Cocoa backend. On the file output
it is visible with PDF, PNG, Gif, JPEG, SVG, TeX ... but not PostScript.
The following macro gives an example of transparency usage:

Begin_Macro(source)
../../../tutorials/graphics/transparency.C
End_Macro

*/

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TColor::TColor(): TNamed()
{
   fNumber = -1;
   fRed = fGreen = fBlue = fHue = fLight = fSaturation = -1;
   fAlpha = 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Normal color constructor. Initialize a color structure.
/// Compute the RGB and HLS color components.

TColor::TColor(Int_t color, Float_t r, Float_t g, Float_t b, const char *name,
               Float_t a)
      : TNamed(name,"")
{
   TColor::InitializeColors();
   // do not enter if color number already exist
   TColor *col = gROOT->GetColor(color);
   if (col) {
      Warning("TColor", "color %d already defined", color);
      fNumber = col->GetNumber();
      fRed    = col->GetRed();
      fGreen  = col->GetGreen();
      fBlue   = col->GetBlue();
      fHue    = col->GetHue();
      fLight  = col->GetLight();
      fAlpha  = col->GetAlpha();
      fSaturation = col->GetSaturation();
      return;
   }

   fNumber = color;

   if (fNumber > gHighestColorIndex) gHighestColorIndex = fNumber;

   char aname[32];
   if (!name || !*name) {
      snprintf(aname,32, "Color%d", color);
      SetName(aname);
   }

   // enter in the list of colors
   TObjArray *lcolors = (TObjArray*)gROOT->GetListOfColors();
   lcolors->AddAtAndExpand(this, color);

   // fill color structure
   SetRGB(r, g, b);
   fAlpha = a;
   gDefinedColors++;
}

////////////////////////////////////////////////////////////////////////////////
/// Fast TColor constructor. It creates a color with an index just above the
/// current highest one. It does not name the color.
/// This is useful to create palettes.

TColor::TColor(Float_t r, Float_t g, Float_t b, Float_t a): TNamed("","")
{
   gHighestColorIndex++;
   fNumber = gHighestColorIndex;
   fRed    = r;
   fGreen  = g;
   fBlue   = b;
   fAlpha  = a;
   RGBtoHLS(r, g, b, fHue, fLight, fSaturation);

   // enter in the list of colors
   TObjArray *lcolors = (TObjArray*)gROOT->GetListOfColors();
   lcolors->AddAtAndExpand(this, fNumber);
   gDefinedColors++;
}

////////////////////////////////////////////////////////////////////////////////
/// Color destructor.

TColor::~TColor()
{
   gROOT->GetListOfColors()->Remove(this);
   if (gROOT->GetListOfColors()->GetEntries() == 0) {fgPalette.Set(0); fgPalette=0;}
}

////////////////////////////////////////////////////////////////////////////////
/// Color copy constructor.

TColor::TColor(const TColor &color) : TNamed(color)
{
   ((TColor&)color).Copy(*this);
}

TColor &TColor::operator=(const TColor &color)
{
   ((TColor &)color).Copy(*this);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize colors used by the TCanvas based graphics (via TColor objects).
/// This method should be called before the ApplicationImp is created (which
/// initializes the GUI colors).

void TColor::InitializeColors()
{
   static Bool_t initDone = kFALSE;

   if (initDone) return;
   initDone = kTRUE;

   if (gROOT->GetListOfColors()->First() == 0) {

      new TColor(kWhite,1,1,1,"background");
      new TColor(kBlack,0,0,0,"black");
      new TColor(2,1,0,0,"red");
      new TColor(3,0,1,0,"green");
      new TColor(4,0,0,1,"blue");
      new TColor(5,1,1,0,"yellow");
      new TColor(6,1,0,1,"magenta");
      new TColor(7,0,1,1,"cyan");
      new TColor(10,0.999,0.999,0.999,"white");
      new TColor(11,0.754,0.715,0.676,"editcol");

      // The color white above is defined as being nearly white.
      // Sets the associated dark color also to white.
      TColor::GetColorDark(10);
      TColor *c110 = gROOT->GetColor(110);
      if (c110) c110->SetRGB(0.999,0.999,.999);

      // Initialize Custom colors
      new TColor(20,0.8,0.78,0.67);
      new TColor(31,0.54,0.66,0.63);
      new TColor(41,0.83,0.81,0.53);
      new TColor(30,0.52,0.76,0.64);
      new TColor(32,0.51,0.62,0.55);
      new TColor(24,0.70,0.65,0.59);
      new TColor(21,0.8,0.78,0.67);
      new TColor(47,0.67,0.56,0.58);
      new TColor(35,0.46,0.54,0.57);
      new TColor(33,0.68,0.74,0.78);
      new TColor(39,0.5,0.5,0.61);
      new TColor(37,0.43,0.48,0.52);
      new TColor(38,0.49,0.6,0.82);
      new TColor(36,0.41,0.51,0.59);
      new TColor(49,0.58,0.41,0.44);
      new TColor(43,0.74,0.62,0.51);
      new TColor(22,0.76,0.75,0.66);
      new TColor(45,0.75,0.51,0.47);
      new TColor(44,0.78,0.6,0.49);
      new TColor(26,0.68,0.6,0.55);
      new TColor(28,0.53,0.4,0.34);
      new TColor(25,0.72,0.64,0.61);
      new TColor(27,0.61,0.56,0.51);
      new TColor(23,0.73,0.71,0.64);
      new TColor(42,0.87,0.73,0.53);
      new TColor(46,0.81,0.37,0.38);
      new TColor(48,0.65,0.47,0.48);
      new TColor(34,0.48,0.56,0.6);
      new TColor(40,0.67,0.65,0.75);
      new TColor(29,0.69,0.81,0.78);

      // Initialize some additional greyish non saturated colors
      new TColor(8, 0.35,0.83,0.33);
      new TColor(9, 0.35,0.33,0.85);
      new TColor(12,.3,.3,.3,"grey12");
      new TColor(13,.4,.4,.4,"grey13");
      new TColor(14,.5,.5,.5,"grey14");
      new TColor(15,.6,.6,.6,"grey15");
      new TColor(16,.7,.7,.7,"grey16");
      new TColor(17,.8,.8,.8,"grey17");
      new TColor(18,.9,.9,.9,"grey18");
      new TColor(19,.95,.95,.95,"grey19");
      new TColor(50, 0.83,0.35,0.33);

      // Initialize the Pretty Palette Spectrum Violet->Red
      //   The color model used here is based on the HLS model which
      //   is much more suitable for creating palettes than RGB.
      //   Fixing the saturation and lightness we can scan through the
      //   spectrum of visible light by using "hue" alone.
      //   In Root hue takes values from 0 to 360.
      Int_t   i;
      Float_t  saturation = 1;
      Float_t  lightness = 0.5;
      Float_t  maxHue = 280;
      Float_t  minHue = 0;
      Int_t    maxPretty = 50;
      Float_t  hue;
      Float_t r=0., g=0., b=0., h, l, s;

      for (i=0 ; i<maxPretty-1 ; i++) {
         hue = maxHue-(i+1)*((maxHue-minHue)/maxPretty);
         TColor::HLStoRGB(hue, lightness, saturation, r, g, b);
         new TColor(i+51, r, g, b);
      }

      // Initialize special colors for x3d
      TColor *s0;
      for (i = 1; i < 8; i++) {
         s0 = gROOT->GetColor(i);
         if (s0) s0->GetRGB(r,g,b);
         if (i == 1) { r = 0.6; g = 0.6; b = 0.6; }
         if (r == 1) r = 0.9; else if (r == 0) r = 0.1;
         if (g == 1) g = 0.9; else if (g == 0) g = 0.1;
         if (b == 1) b = 0.9; else if (b == 0) b = 0.1;
         TColor::RGBtoHLS(r,g,b,h,l,s);
         TColor::HLStoRGB(h,0.6*l,s,r,g,b);
         new TColor(200+4*i-3,r,g,b);
         TColor::HLStoRGB(h,0.8*l,s,r,g,b);
         new TColor(200+4*i-2,r,g,b);
         TColor::HLStoRGB(h,1.2*l,s,r,g,b);
         new TColor(200+4*i-1,r,g,b);
         TColor::HLStoRGB(h,1.4*l,s,r,g,b);
         new TColor(200+4*i  ,r,g,b);
      }

      // Create the ROOT Color Wheel
      TColor::CreateColorWheel();
   }
   // If fgPalette.fN !=0 SetPalette has been called already
   // (from rootlogon.C for instance)

   if (!fgPalette.fN) SetPalette(1,0);
}

////////////////////////////////////////////////////////////////////////////////
/// Return color as hexadecimal string. This string can be directly passed
/// to, for example, TGClient::GetColorByName(). String will be reused so
/// copy immediately if needed.

const char *TColor::AsHexString() const
{
   static TString tempbuf;

   Int_t r, g, b, a;
   r = Int_t(GetRed()   * 255);
   g = Int_t(GetGreen() * 255);
   b = Int_t(GetBlue()  * 255);
   a = Int_t(fAlpha     * 255);

   if (a != 255) {
      tempbuf.Form("#%02x%02x%02x%02x", a, r, g, b);
   } else {
      tempbuf.Form("#%02x%02x%02x", r, g, b);
   }
   return tempbuf;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this color to obj.

void TColor::Copy(TObject &obj) const
{
   TNamed::Copy((TNamed&)obj);
   ((TColor&)obj).fRed   = fRed;
   ((TColor&)obj).fGreen = fGreen;
   ((TColor&)obj).fBlue  = fBlue;
   ((TColor&)obj).fHue   = fHue;
   ((TColor&)obj).fLight = fLight;
   ((TColor&)obj).fAlpha = fAlpha;
   ((TColor&)obj).fSaturation = fSaturation;
   ((TColor&)obj).fNumber = fNumber;
}

////////////////////////////////////////////////////////////////////////////////
/// Create the Gray scale colors in the Color Wheel

void TColor::CreateColorsGray()
{
   if (gROOT->GetColor(kGray)) return;
   TColor *gray  = new TColor(kGray,204./255.,204./255.,204./255.);
   TColor *gray1 = new TColor(kGray+1,153./255.,153./255.,153./255.);
   TColor *gray2 = new TColor(kGray+2,102./255.,102./255.,102./255.);
   TColor *gray3 = new TColor(kGray+3, 51./255., 51./255., 51./255.);
   gray ->SetName("kGray");
   gray1->SetName("kGray+1");
   gray2->SetName("kGray+2");
   gray3->SetName("kGray+3");
}

////////////////////////////////////////////////////////////////////////////////
/// Create the "circle" colors in the color wheel.

void TColor::CreateColorsCircle(Int_t offset, const char *name, UChar_t *rgb)
{
   TString colorname;
   for (Int_t n=0;n<15;n++) {
      Int_t colorn = offset+n-10;
      TColor *color = gROOT->GetColor(colorn);
      if (!color) {
         color = new TColor(colorn,rgb[3*n]/255.,rgb[3*n+1]/255.,rgb[3*n+2]/255.);
         color->SetTitle(color->AsHexString());
         if      (n>10) colorname.Form("%s+%d",name,n-10);
         else if (n<10) colorname.Form("%s-%d",name,10-n);
         else           colorname.Form("%s",name);
         color->SetName(colorname);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create the "rectangular" colors in the color wheel.

void TColor::CreateColorsRectangle(Int_t offset, const char *name, UChar_t *rgb)
{
   TString colorname;
   for (Int_t n=0;n<20;n++) {
      Int_t colorn = offset+n-9;
      TColor *color = gROOT->GetColor(colorn);
      if (!color) {
         color = new TColor(colorn,rgb[3*n]/255.,rgb[3*n+1]/255.,rgb[3*n+2]/255.);
         color->SetTitle(color->AsHexString());
         if      (n>9) colorname.Form("%s+%d",name,n-9);
         else if (n<9) colorname.Form("%s-%d",name,9-n);
         else          colorname.Form("%s",name);
         color->SetName(colorname);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Static function steering the creation of all colors in the color wheel.

void TColor::CreateColorWheel()
{
   UChar_t magenta[46]= {255,204,255
                        ,255,153,255, 204,153,204
                        ,255,102,255, 204,102,204, 153,102,153
                        ,255, 51,255, 204, 51,204, 153, 51,153, 102, 51,102
                        ,255,  0,255, 204,  0,204, 153,  0,153, 102,  0,102,  51,  0, 51};

   UChar_t red[46]    = {255,204,204
                        ,255,153,153, 204,153,153
                        ,255,102,102, 204,102,102, 153,102,102
                        ,255, 51, 51, 204, 51, 51, 153, 51, 51, 102, 51, 51
                        ,255,  0,  0, 204,  0,  0, 153,  0,  0, 102,  0,  0,  51,  0,  0};

   UChar_t yellow[46] = {255,255,204
                        ,255,255,153, 204,204,153
                        ,255,255,102, 204,204,102, 153,153,102
                        ,255,255, 51, 204,204, 51, 153,153, 51, 102,102, 51
                        ,255,255,  0, 204,204,  0, 153,153,  0, 102,102,  0,  51, 51,  0};

   UChar_t green[46]  = {204,255,204
                        ,153,255,153, 153,204,153
                        ,102,255,102, 102,204,102, 102,153,102
                        , 51,255, 51,  51,204, 51,  51,153, 51,  51,102, 51
                        ,  0,255,  0,   0,204,  0,   0,153,  0,   0,102,  0,  0, 51,  0};

   UChar_t cyan[46]   = {204,255,255
                        ,153,255,255, 153,204,204
                        ,102,255,255, 102,204,204, 102,153,153
                        , 51,255,255,  51,204,204,  51,153,153,  51,102,102
                        ,  0,255,255,   0,204,204,   0,153,153,   0,102,102,   0, 51,  51};

   UChar_t blue[46]   = {204,204,255
                        ,153,153,255, 153,153,204
                        ,102,102,255, 102,102,204, 102,102,153
                        , 51, 51,255,  51, 51,204,  51, 51,153,  51, 51,102
                        ,  0,  0,255,   0,  0,204,   0,  0,153,   0,  0,102,   0,  0,  51};

   UChar_t pink[60] = {255, 51,153,  204,  0,102,  102,  0, 51,  153,  0, 51,  204, 51,102
                      ,255,102,153,  255,  0,102,  255, 51,102,  204,  0, 51,  255,  0, 51
                      ,255,153,204,  204,102,153,  153, 51,102,  153,  0,102,  204, 51,153
                      ,255,102,204,  255,  0,153,  204,  0,153,  255, 51,204,  255,  0,153};

   UChar_t orange[60]={255,204,153,  204,153,102,  153,102, 51,  153,102,  0,  204,153, 51
                      ,255,204,102,  255,153,  0,  255,204, 51,  204,153,  0,  255,204,  0
                      ,255,153, 51,  204,102,  0,  102, 51,  0,  153, 51,  0,  204,102, 51
                      ,255,153,102,  255,102,  0,  255,102, 51,  204, 51,  0,  255, 51,  0};

   UChar_t spring[60]={153,255, 51,  102,204,  0,   51,102,  0,   51,153,  0,  102,204, 51
                      ,153,255,102,  102,255,  0,  102,255, 51,   51,204,  0,   51,255, 0
                      ,204,255,153,  153,204,102,  102,153, 51,  102,153,  0,  153,204, 51
                      ,204,255,102,  153,255,  0,  204,255, 51,  153,204,  0,  204,255,  0};

   UChar_t teal[60] = {153,255,204,  102,204,153,   51,153,102,    0,153,102,   51,204,153
                      ,102,255,204,    0,255,102,   51,255,204,    0,204,153,    0,255,204
                      , 51,255,153,    0,204,102,    0,102, 51,    0,153, 51,   51,204,102
                      ,102,255,153,    0,255,153,   51,255,102,    0,204, 51,    0,255, 51};

   UChar_t azure[60] ={153,204,255,  102,153,204,   51,102,153,    0, 51,153,   51,102,204
                      ,102,153,255,    0,102,255,   51,102,255,    0, 51,204,    0, 51,255
                      , 51,153,255,    0,102,204,    0, 51,102,    0,102,153,   51,153,204
                      ,102,204,255,    0,153,255,   51,204,255,    0,153,204,    0,204,255};

   UChar_t violet[60]={204,153,255,  153,102,204,  102, 51,153,  102,  0,153,  153, 51,204
                      ,204,102,255,  153,  0,255,  204, 51,255,  153,  0,204,  204,  0,255
                      ,153, 51,255,  102,  0,204,   51,  0,102,   51,  0,153,  102, 51,204
                      ,153,102,255,  102,  0,255,  102, 51,255,   51,  0,204,   51,  0,255};

   TColor::CreateColorsCircle(kMagenta,"kMagenta",magenta);
   TColor::CreateColorsCircle(kRed,    "kRed",    red);
   TColor::CreateColorsCircle(kYellow, "kYellow", yellow);
   TColor::CreateColorsCircle(kGreen,  "kGreen",  green);
   TColor::CreateColorsCircle(kCyan,   "kCyan",   cyan);
   TColor::CreateColorsCircle(kBlue,   "kBlue",   blue);

   TColor::CreateColorsRectangle(kPink,  "kPink",  pink);
   TColor::CreateColorsRectangle(kOrange,"kOrange",orange);
   TColor::CreateColorsRectangle(kSpring,"kSpring",spring);
   TColor::CreateColorsRectangle(kTeal,  "kTeal",  teal);
   TColor::CreateColorsRectangle(kAzure, "kAzure", azure);
   TColor::CreateColorsRectangle(kViolet,"kViolet",violet);

   TColor::CreateColorsGray();
}

////////////////////////////////////////////////////////////////////////////////
/// Static function returning the color number i in current palette.

Int_t TColor::GetColorPalette(Int_t i)
{
   Int_t ncolors = fgPalette.fN;
   if (ncolors == 0) return 0;
   Int_t icol    = i%ncolors;
   if (icol < 0) icol = 0;
   return fgPalette.fArray[icol];
}

////////////////////////////////////////////////////////////////////////////////
/// Static function returning the current active palette.

const TArrayI& TColor::GetPalette()
{
   return fgPalette;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function returning number of colors in the color palette.

Int_t TColor::GetNumberOfColors()
{
   return fgPalette.fN;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function returning kTRUE if some new colors have been defined after
/// initialisation or since the last call to this method. This allows to avoid
/// the colors and palette streaming in TCanvas::Streamer if not needed.

Bool_t TColor::DefinedColors()
{
   // After initialization gDefinedColors == 649. If it is bigger it means some new
   // colors have been defined
   Bool_t hasChanged = (gDefinedColors - gLastDefinedColors) > 50;
   gLastDefinedColors = gDefinedColors;
   return hasChanged;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pixel value corresponding to this color. This pixel value can
/// be used in the GUI classes. This call does not work in batch mode since
/// it needs to communicate with the graphics system.

ULong_t TColor::GetPixel() const
{
   if (gVirtualX && !gROOT->IsBatch()) {
      if (gApplication) {
         TApplication::NeedGraphicsLibs();
         gApplication->InitializeGraphics();
      }
      return gVirtualX->GetPixel(fNumber);
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method to compute RGB from HLS. The l and s are between [0,1]
/// and h is between [0,360]. The returned r,g,b triplet is between [0,1].

void TColor::HLS2RGB(Float_t hue, Float_t light, Float_t satur,
                     Float_t &r, Float_t &g, Float_t &b)
{

   Float_t rh, rl, rs, rm1, rm2;
   rh = rl = rs = 0;
   if (hue   > 0) { rh = hue;   if (rh > 360) rh = 360; }
   if (light > 0) { rl = light; if (rl > 1)   rl = 1; }
   if (satur > 0) { rs = satur; if (rs > 1)   rs = 1; }

   if (rl <= 0.5)
      rm2 = rl*(1.0f + rs);
   else
      rm2 = rl + rs - rl*rs;
   rm1 = 2.0f*rl - rm2;

   if (!rs) { r = rl; g = rl; b = rl; return; }
   r = HLStoRGB1(rm1, rm2, rh+120.0f);
   g = HLStoRGB1(rm1, rm2, rh);
   b = HLStoRGB1(rm1, rm2, rh-120.0f);
}

////////////////////////////////////////////////////////////////////////////////
/// Static method. Auxiliary to HLS2RGB().

Float_t TColor::HLStoRGB1(Float_t rn1, Float_t rn2, Float_t huei)
{
   Float_t hue = huei;
   if (hue > 360) hue = hue - 360.0f;
   if (hue < 0)   hue = hue + 360.0f;
   if (hue < 60 ) return rn1 + (rn2-rn1)*hue/60.0f;
   if (hue < 180) return rn2;
   if (hue < 240) return rn1 + (rn2-rn1)*(240.0f-hue)/60.0f;
   return rn1;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method to compute RGB from HLS. The h,l,s are between [0,255].
/// The returned r,g,b triplet is between [0,255].

void TColor::HLS2RGB(Int_t h, Int_t l, Int_t s, Int_t &r, Int_t &g, Int_t &b)
{
   Float_t hh, ll, ss, rr, gg, bb;

   hh = Float_t(h) * 360.0f / 255.0f;
   ll = Float_t(l) / 255.0f;
   ss = Float_t(s) / 255.0f;

   TColor::HLStoRGB(hh, ll, ss, rr, gg, bb);

   r = (Int_t) (rr * 255.0f);
   g = (Int_t) (gg * 255.0f);
   b = (Int_t) (bb * 255.0f);
}

////////////////////////////////////////////////////////////////////////////////
/// Static method to compute RGB from HSV:
///
///  -  The hue value runs from 0 to 360.
///  -  The saturation is the degree of strength or purity and is from 0 to 1.
///     Purity is how much white is added to the color, so S=1 makes the purest
///     color (no white).
///  -  Brightness value also ranges from 0 to 1, where 0 is the black.
///
/// The returned r,g,b triplet is between [0,1].

void TColor::HSV2RGB(Float_t hue, Float_t satur, Float_t value,
                     Float_t &r, Float_t &g, Float_t &b)
{
   Int_t i;
   Float_t f, p, q, t;

   if (satur==0) {
      // Achromatic (grey)
      r = g = b = value;
      return;
   }

   hue /= 60.0f;   // sector 0 to 5
   i = (Int_t)floor(hue);
   f = hue-i;   // factorial part of hue
   p = value*(1-satur);
   q = value*(1-satur*f );
   t = value*(1-satur*(1-f));

   switch (i) {
      case 0:
         r = value;
         g = t;
         b = p;
         break;
      case 1:
         r = q;
         g = value;
         b = p;
         break;
      case 2:
         r = p;
         g = value;
         b = t;
         break;
      case 3:
         r = p;
         g = q;
         b = value;
         break;
      case 4:
         r = t;
         g = p;
         b = value;
         break;
      default:
         r = value;
         g = p;
         b = q;
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// List this color with its attributes.

void TColor::ls(Option_t *) const
{
   printf("Color:%d  Red=%f Green=%f Blue=%f Alpha=%f Name=%s\n",
          fNumber, fRed, fGreen, fBlue, fAlpha, GetName());
}

////////////////////////////////////////////////////////////////////////////////
/// Dump this color with its attributes.

void TColor::Print(Option_t *) const
{
   ls();
}

////////////////////////////////////////////////////////////////////////////////
/// Static method to compute HLS from RGB. The r,g,b triplet is between
/// [0,1], hue is between [0,360], light and satur are [0,1].

void TColor::RGB2HLS(Float_t rr, Float_t gg, Float_t bb,
                     Float_t &hue, Float_t &light, Float_t &satur)
{
   Float_t rnorm, gnorm, bnorm, minval, maxval, msum, mdiff, r, g, b;
   minval = maxval =0 ;
   r = g = b = 0;
   if (rr > 0) { r = rr; if (r > 1) r = 1; }
   if (gg > 0) { g = gg; if (g > 1) g = 1; }
   if (bb > 0) { b = bb; if (b > 1) b = 1; }

   minval = r;
   if (g < minval) minval = g;
   if (b < minval) minval = b;
   maxval = r;
   if (g > maxval) maxval = g;
   if (b > maxval) maxval = b;

   rnorm = gnorm = bnorm = 0;
   mdiff = maxval - minval;
   msum  = maxval + minval;
   light = 0.5f * msum;
   if (maxval != minval) {
      rnorm = (maxval - r)/mdiff;
      gnorm = (maxval - g)/mdiff;
      bnorm = (maxval - b)/mdiff;
   } else {
      satur = hue = 0;
      return;
   }

   if (light < 0.5)
      satur = mdiff/msum;
   else
      satur = mdiff/(2.0f - msum);

   if (r == maxval)
      hue = 60.0f * (6.0f + bnorm - gnorm);
   else if (g == maxval)
      hue = 60.0f * (2.0f + rnorm - bnorm);
   else
      hue = 60.0f * (4.0f + gnorm - rnorm);

   if (hue > 360)
      hue = hue - 360.0f;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method to compute HSV from RGB.
///
///  -  The input values:
///    -  r,g,b triplet is between [0,1].
///  -  The returned values:
///    -  The hue value runs from 0 to 360.
///    -  The saturation is the degree of strength or purity and is from 0 to 1.
///       Purity is how much white is added to the color, so S=1 makes the purest
///       color (no white).
///    -  Brightness value also ranges from 0 to 1, where 0 is the black.

void TColor::RGB2HSV(Float_t r, Float_t g, Float_t b,
                     Float_t &hue, Float_t &satur, Float_t &value)
{
   Float_t min, max, delta;

   min   = TMath::Min(TMath::Min(r, g), b);
   max   = TMath::Max(TMath::Max(r, g), b);
   value = max;

   delta = max - min;

   if (max != 0) {
      satur = delta/max;
   } else {
      satur = 0;
      hue   = -1;
      return;
   }

   if (r == max) {
      hue = (g-b)/delta;
   } else if (g == max) {
      hue = 2.0f+(b-r)/delta;
   } else {
      hue = 4.0f+(r-g)/delta;
   }

   hue *= 60.0f;
   if (hue < 0.0f) hue += 360.0f;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method to compute HLS from RGB. The r,g,b triplet is between
/// [0,255], hue, light and satur are between [0,255].

void TColor::RGB2HLS(Int_t r, Int_t g, Int_t b, Int_t &h, Int_t &l, Int_t &s)
{
   Float_t rr, gg, bb, hue, light, satur;

   rr = Float_t(r) / 255.0f;
   gg = Float_t(g) / 255.0f;
   bb = Float_t(b) / 255.0f;

   TColor::RGBtoHLS(rr, gg, bb, hue, light, satur);

   h = (Int_t) (hue/360.0f * 255.0f);
   l = (Int_t) (light * 255.0f);
   s = (Int_t) (satur * 255.0f);
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize this color and its associated colors.

void TColor::SetRGB(Float_t r, Float_t g, Float_t b)
{
   TColor::InitializeColors();
   fRed   = r;
   fGreen = g;
   fBlue  = b;

   if (fRed < 0) return;

   RGBtoHLS(r, g, b, fHue, fLight, fSaturation);

   Int_t nplanes = 16;
   if (gVirtualX) gVirtualX->GetPlanes(nplanes);
   if (nplanes == 0) nplanes = 16;

   // allocate color now (can be delayed when we have a large colormap)
#ifndef R__WIN32
   if (nplanes < 15)
#endif
      Allocate();

   if (fNumber > 50) return;

   // now define associated colors for WBOX shading
   Float_t dr, dg, db, lr, lg, lb;

   // set dark color
   HLStoRGB(fHue, 0.7f*fLight, fSaturation, dr, dg, db);
   TColor *dark = gROOT->GetColor(100+fNumber);
   if (dark) {
      if (nplanes > 8) dark->SetRGB(dr, dg, db);
      else             dark->SetRGB(0.3f,0.3f,0.3f);
   }

   // set light color
   HLStoRGB(fHue, 1.2f*fLight, fSaturation, lr, lg, lb);
   TColor *light = gROOT->GetColor(150+fNumber);
   if (light) {
      if (nplanes > 8) light->SetRGB(lr, lg, lb);
      else             light->SetRGB(0.8f,0.8f,0.8f);
   }
   gDefinedColors++;
}

////////////////////////////////////////////////////////////////////////////////
/// Make this color known to the graphics system.

void TColor::Allocate()
{
   if (gVirtualX && !gROOT->IsBatch())

      gVirtualX->SetRGB(fNumber, GetRed(), GetGreen(), GetBlue());
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning color number for color specified by
/// hex color string of form: "#rrggbb", where rr, gg and bb are in
/// hex between [0,FF], e.g. "#c0c0c0".
///
/// The color retrieval is done using a threshold defined by SetColorThreshold.
///
/// If specified color does not exist it will be created with as
/// name "#rrggbb" with rr, gg and bb in hex between [0,FF].

Int_t TColor::GetColor(const char *hexcolor)
{
   if (hexcolor && *hexcolor == '#') {
      Int_t r, g, b;
      if (sscanf(hexcolor+1, "%02x%02x%02x", &r, &g, &b) == 3)
         return GetColor(r, g, b);
   }
   ::Error("TColor::GetColor(const char*)", "incorrect color string");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning color number for color specified by
/// r, g and b. The r,g,b should be in the range [0,1].
///
/// The color retrieval is done using a threshold defined by SetColorThreshold.
///
/// If specified color does not exist it will be created
/// with as name "#rrggbb" with rr, gg and bb in hex between
/// [0,FF].

Int_t TColor::GetColor(Float_t r, Float_t g, Float_t b)
{
   Int_t rr, gg, bb;
   rr = Int_t(r * 255);
   gg = Int_t(g * 255);
   bb = Int_t(b * 255);

   return GetColor(rr, gg, bb);
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning color number for color specified by
/// system dependent pixel value. Pixel values can be obtained, e.g.,
/// from the GUI color picker.
///
/// The color retrieval is done using a threshold defined by SetColorThreshold.


Int_t TColor::GetColor(ULong_t pixel)
{
   Int_t r, g, b;

   Pixel2RGB(pixel, r, g, b);

   return GetColor(r, g, b);
}

////////////////////////////////////////////////////////////////////////////////
/// This method specifies the color threshold used by GetColor to retrieve a color.
///
/// \param[in] t   Color threshold. By default is equal to 1./31. or 1./255.
///                depending on the number of available color planes.
///
/// When GetColor is called, it scans the defined colors and compare them to the
/// requested color.
/// If the Red Green and Blue values passed to GetColor are Rr Gr Br
/// and Rd Gd Bd the values of a defined color. These two colors are considered equal
/// if (abs(Rr-Rd) < t  & abs(Br-Bd) < t & abs(Br-Bd) < t). If this test passes,
/// the color defined by Rd Gd Bd is returned by GetColor.
///
/// To make sure GetColor will return a color having exactly the requested
/// R G B values it is enough to specify a nul :
/// ~~~ {.cpp}
///   TColor::SetColorThreshold(0.);
/// ~~~
///
/// To reset the color threshold to its default value it is enough to do:
/// ~~~ {.cpp}
///   TColor::SetColorThreshold(-1.);
/// ~~~

void TColor::SetColorThreshold(Float_t t)
{
   gColorThreshold = t;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning color number for color specified by
/// r, g and b. The r,g,b should be in the range [0,255].
/// If the specified color does not exist it will be created
/// with as name "#rrggbb" with rr, gg and bb in hex between
/// [0,FF].
///
/// The color retrieval is done using a threshold defined by SetColorThreshold.


Int_t TColor::GetColor(Int_t r, Int_t g, Int_t b)
{
   TColor::InitializeColors();
   if (r < 0) r = 0;
   if (g < 0) g = 0;
   if (b < 0) b = 0;
   if (r > 255) r = 255;
   if (g > 255) g = 255;
   if (b > 255) b = 255;

   // Get list of all defined colors
   TObjArray *colors = (TObjArray*) gROOT->GetListOfColors();

   TColor *color = 0;

   // Look for color by name
   if ((color = (TColor*) colors->FindObject(Form("#%02x%02x%02x", r, g, b))))
      // We found the color by name, so we use that right away
      return color->GetNumber();

   Float_t rr, gg, bb;
   rr = Float_t(r)/255.0f;
   gg = Float_t(g)/255.0f;
   bb = Float_t(b)/255.0f;

   TIter next(colors);

   Float_t thres;
   if (gColorThreshold >= 0) {
      thres = gColorThreshold;
   } else {
      Int_t nplanes = 16;
      thres = 1.0f/31.0f;   // 5 bits per color : 0 - 0x1F !
      if (gVirtualX) gVirtualX->GetPlanes(nplanes);
      if (nplanes >= 24) thres = 1.0f/255.0f;       // 8 bits per color : 0 - 0xFF !
   }

   // Loop over all defined colors
   while ((color = (TColor*)next())) {
      if (TMath::Abs(color->GetRed() - rr) > thres)   continue;
      if (TMath::Abs(color->GetGreen() - gg) > thres) continue;
      if (TMath::Abs(color->GetBlue() - bb) > thres)  continue;
      // We found a matching color in the color table
      return color->GetNumber();
   }

   // We didn't find a matching color in the color table, so we
   // add it. Note name is of the form "#rrggbb" where rr, etc. are
   // hexadecimal numbers.
   color = new TColor(colors->GetLast()+1, rr, gg, bb,
                      Form("#%02x%02x%02x", r, g, b));

   return color->GetNumber();
}

////////////////////////////////////////////////////////////////////////////////
/// Static function: Returns the bright color number corresponding to n
/// If the TColor object does not exist, it is created.
/// The convention is that the bright color nb = n+150

Int_t TColor::GetColorBright(Int_t n)
{
   if (n < 0) return -1;

   // Get list of all defined colors
   TObjArray *colors = (TObjArray*) gROOT->GetListOfColors();
   Int_t ncolors = colors->GetSize();
   // Get existing color at index n
   TColor *color = 0;
   if (n < ncolors) color = (TColor*)colors->At(n);
   if (!color) return -1;

   //Get the rgb of the the new bright color corresponding to color n
   Float_t r,g,b;
   HLStoRGB(color->GetHue(), 1.2f*color->GetLight(), color->GetSaturation(), r, g, b);

   //Build the bright color (unless the slot nb is already used)
   Int_t nb = n+150;
   TColor *colorb = 0;
   if (nb < ncolors) colorb = (TColor*)colors->At(nb);
   if (colorb) return nb;
   colorb = new TColor(nb,r,g,b);
   colorb->SetName(Form("%s_bright",color->GetName()));
   colors->AddAtAndExpand(colorb,nb);
   return nb;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function: Returns the dark color number corresponding to n
/// If the TColor object does not exist, it is created.
/// The convention is that the dark color nd = n+100

Int_t TColor::GetColorDark(Int_t n)
{
   if (n < 0) return -1;

   // Get list of all defined colors
   TObjArray *colors = (TObjArray*) gROOT->GetListOfColors();
   Int_t ncolors = colors->GetSize();
   // Get existing color at index n
   TColor *color = 0;
   if (n < ncolors) color = (TColor*)colors->At(n);
   if (!color) return -1;

   //Get the rgb of the the new dark color corresponding to color n
   Float_t r,g,b;
   HLStoRGB(color->GetHue(), 0.7f*color->GetLight(), color->GetSaturation(), r, g, b);

   //Build the dark color (unless the slot nd is already used)
   Int_t nd = n+100;
   TColor *colord = 0;
   if (nd < ncolors) colord = (TColor*)colors->At(nd);
   if (colord) return nd;
   colord = new TColor(nd,r,g,b);
   colord->SetName(Form("%s_dark",color->GetName()));
   colors->AddAtAndExpand(colord,nd);
   return nd;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function: Returns the transparent color number corresponding to n.
/// The transparency level is given by the alpha value a.

Int_t TColor::GetColorTransparent(Int_t n, Float_t a)
{
   if (n < 0) return -1;

   TColor *color = gROOT->GetColor(n);
   if (color) {
      TColor *colort = new TColor(gROOT->GetListOfColors()->GetLast()+1,
                                  color->GetRed(), color->GetGreen(), color->GetBlue());
      colort->SetAlpha(a);
      colort->SetName(Form("%s_transparent",color->GetName()));
      return colort->GetNumber();
   } else {
      ::Error("TColor::GetColorTransparent", "color with index %d not defined", n);
      return -1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Static function: Returns a free color index which can be used to define
/// a user custom color.
///
/// ~~~ {.cpp}
///   Int_t ci = TColor::GetFreeColorIndex();
///   TColor *color = new TColor(ci, 0.1, 0.2, 0.3);
/// ~~~

Int_t TColor::GetFreeColorIndex()
{
   return gHighestColorIndex+1;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method that given a color index number, returns the corresponding
/// pixel value. This pixel value can be used in the GUI classes. This call
/// does not work in batch mode since it needs to communicate with the
/// graphics system.

ULong_t TColor::Number2Pixel(Int_t ci)
{
   TColor::InitializeColors();
   TColor *color = gROOT->GetColor(ci);
   if (color)
      return color->GetPixel();
   else
      ::Warning("TColor::Number2Pixel", "color with index %d not defined", ci);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert r,g,b to graphics system dependent pixel value.
/// The r,g,b triplet must be [0,1].

ULong_t TColor::RGB2Pixel(Float_t r, Float_t g, Float_t b)
{
   if (r < 0) r = 0;
   if (g < 0) g = 0;
   if (b < 0) b = 0;
   if (r > 1) r = 1;
   if (g > 1) g = 1;
   if (b > 1) b = 1;

   ColorStruct_t color;
   color.fRed   = UShort_t(r * 65535);
   color.fGreen = UShort_t(g * 65535);
   color.fBlue  = UShort_t(b * 65535);
   color.fMask  = kDoRed | kDoGreen | kDoBlue;
   gVirtualX->AllocColor(gVirtualX->GetColormap(), color);
   return color.fPixel;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert r,g,b to graphics system dependent pixel value.
/// The r,g,b triplet must be [0,255].

ULong_t TColor::RGB2Pixel(Int_t r, Int_t g, Int_t b)
{
   if (r < 0) r = 0;
   if (g < 0) g = 0;
   if (b < 0) b = 0;
   if (r > 255) r = 255;
   if (g > 255) g = 255;
   if (b > 255) b = 255;

   ColorStruct_t color;
   color.fRed   = UShort_t(r * 257);  // 65535/255
   color.fGreen = UShort_t(g * 257);
   color.fBlue  = UShort_t(b * 257);
   color.fMask  = kDoRed | kDoGreen | kDoBlue;
   gVirtualX->AllocColor(gVirtualX->GetColormap(), color);
   return color.fPixel;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert machine dependent pixel value (obtained via RGB2Pixel or
/// via Number2Pixel() or via TColor::GetPixel()) to r,g,b triplet.
/// The r,g,b triplet will be [0,1].

void TColor::Pixel2RGB(ULong_t pixel, Float_t &r, Float_t &g, Float_t &b)
{
   ColorStruct_t color;
   color.fPixel = pixel;
   gVirtualX->QueryColor(gVirtualX->GetColormap(), color);
   r = (Float_t)color.fRed / 65535.0f;
   g = (Float_t)color.fGreen / 65535.0f;
   b = (Float_t)color.fBlue / 65535.0f;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert machine dependent pixel value (obtained via RGB2Pixel or
/// via Number2Pixel() or via TColor::GetPixel()) to r,g,b triplet.
/// The r,g,b triplet will be [0,255].

void TColor::Pixel2RGB(ULong_t pixel, Int_t &r, Int_t &g, Int_t &b)
{
   ColorStruct_t color;
   color.fPixel = pixel;
   gVirtualX->QueryColor(gVirtualX->GetColormap(), color);
   r = color.fRed / 257;
   g = color.fGreen / 257;
   b = color.fBlue / 257;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert machine dependent pixel value (obtained via RGB2Pixel or
/// via Number2Pixel() or via TColor::GetPixel()) to a hexadecimal string.
/// This string can be directly passed to, for example,
/// TGClient::GetColorByName(). String will be reused so copy immediately
/// if needed.

const char *TColor::PixelAsHexString(ULong_t pixel)
{
   static TString tempbuf;
   Int_t r, g, b;
   Pixel2RGB(pixel, r, g, b);
   tempbuf.Form("#%02x%02x%02x", r, g, b);
   return tempbuf;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a color with index > 228 as a C++ statement(s) on output stream out.

void TColor::SaveColor(std::ostream &out, Int_t ci)
{
   char quote = '"';
   Float_t r,g,b,a;
   Int_t ri, gi, bi;
   TString cname;

   TColor *c = gROOT->GetColor(ci);
   if (c) {
      c->GetRGB(r, g, b);
      a = c->GetAlpha();
   } else {
      return;
   }

   if (gROOT->ClassSaved(TColor::Class())) {
      out << std::endl;
   } else {
      out << std::endl;
      out << "   Int_t ci;      // for color index setting" << std::endl;
      out << "   TColor *color; // for color definition with alpha" << std::endl;
   }

   if (a<1) {
      out<<"   ci = "<<ci<<";"<<std::endl;
      out<<"   color = new TColor(ci, "<<r<<", "<<g<<", "<<b<<", "
      <<"\" \", "<<a<<");"<<std::endl;
   } else {
      ri = (Int_t)(255*r);
      gi = (Int_t)(255*g);
      bi = (Int_t)(255*b);
      cname.Form("#%02x%02x%02x", ri, gi, bi);
      out<<"   ci = TColor::GetColor("<<quote<<cname.Data()<<quote<<");"<<std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return whether all colors return grayscale values.
Bool_t TColor::IsGrayscale()
{
   return fgGrayscaleMode;
}

////////////////////////////////////////////////////////////////////////////////
/// Set whether all colors should return grayscale values.

void TColor::SetGrayscale(Bool_t set /*= kTRUE*/)
{
   if (fgGrayscaleMode == set) return;

   fgGrayscaleMode = set;

   if (!gVirtualX || gROOT->IsBatch()) return;

   TColor::InitializeColors();
   TIter iColor(gROOT->GetListOfColors());
   TColor* color = 0;
   while ((color = (TColor*) iColor()))
      color->Allocate();
}

////////////////////////////////////////////////////////////////////////////////
/// Static function creating a color table with several connected linear gradients.
///
///   - Number: The number of end point colors that will form the gradients.
///             Must be at least 2.
///   - Stops: Where in the whole table the end point colors should lie.
///            Each entry must be on [0, 1], each entry must be greater than
///            the previous entry.
///   - Red, Green, Blue: The end point color values.
///                       Each entry must be on [0, 1]
///   - NColors: Total number of colors in the table. Must be at least 1.
///
/// Returns a positive value on success and -1 on error.
///
/// The table is constructed by tracing lines between the given points in
/// RGB space.  Each color value may have a value between 0 and 1.  The
/// difference between consecutive "Stops" values gives the fraction of
/// space in the whole table that should be used for the interval between
/// the corresponding color values.
///
/// Normally the first element of Stops should be 0 and the last should be 1.
/// If this is not true, fewer than NColors will be used in proportion with
/// the total interval between the first and last elements of Stops.
///
/// This definition is similar to the povray-definition of gradient
/// color tables.
///
/// For instance:
/// ~~~ {.cpp}
/// UInt_t Number = 3;
/// Double_t Red[3]   = { 0.0, 1.0, 1.0 };
/// Double_t Green[3] = { 0.0, 0.0, 1.0 };
/// Double_t Blue[3]  = { 1.0, 0.0, 1.0 };
/// Double_t Stops[3] = { 0.0, 0.4, 1.0 };
/// ~~~
/// This defines a table in which there are three color end points:
/// RGB = {0, 0, 1}, {1, 0, 0}, and {1, 1, 1} = blue, red, white
/// The first 40% of the table is used to go linearly from blue to red.
/// The remaining 60% of the table is used to go linearly from red to white.
///
/// If you define a very short interval such that less than one color fits
/// in it, no colors at all will be allocated.  If this occurs for all
/// intervals, ROOT will revert to the default palette.
///
/// Original code by Andreas Zoglauer (zog@mpe.mpg.de)

Int_t TColor::CreateGradientColorTable(UInt_t Number, Double_t* Stops,
                              Double_t* Red, Double_t* Green,
                              Double_t* Blue, UInt_t NColors, Float_t alpha)
{
   TColor::InitializeColors();

   UInt_t g, c;
   UInt_t nPalette = 0;
   Int_t *palette = new Int_t[NColors+1];
   UInt_t nColorsGradient;

   if(Number < 2 || NColors < 1){
      delete [] palette;
      return -1;
   }

   // Check if all RGB values are between 0.0 and 1.0 and
   // Stops goes from 0.0 to 1.0 in increasing order.
   for (c = 0; c < Number; c++) {
      if (Red[c] < 0 || Red[c] > 1.0 ||
          Green[c] < 0 || Green[c] > 1.0 ||
          Blue[c] < 0 || Blue[c] > 1.0 ||
          Stops[c] < 0 || Stops[c] > 1.0) {
         delete [] palette;
         return -1;
      }
      if (c >= 1) {
         if (Stops[c-1] > Stops[c]) {
            delete [] palette;
            return -1;
         }
      }
   }

   // Now create the colors and add them to the default palette:

   // For each defined gradient...
   for (g = 1; g < Number; g++) {
      // create the colors...
      nColorsGradient = (Int_t) (floor(NColors*Stops[g]) - floor(NColors*Stops[g-1]));
      for (c = 0; c < nColorsGradient; c++) {
         new TColor( Float_t(Red[g-1]   + c * (Red[g]   - Red[g-1])  / nColorsGradient),
                     Float_t(Green[g-1] + c * (Green[g] - Green[g-1])/ nColorsGradient),
                     Float_t(Blue[g-1]  + c * (Blue[g]  - Blue[g-1]) / nColorsGradient),
                     alpha);
         palette[nPalette] = gHighestColorIndex;
         nPalette++;
      }
   }

   TColor::SetPalette(nPalette, palette);
   delete [] palette;
   return gHighestColorIndex + 1 - NColors;
}


////////////////////////////////////////////////////////////////////////////////
/// Static function.
/// The color palette is used by the histogram classes
///  (see TH1::Draw options).
/// For example TH1::Draw("col") draws a 2-D histogram with cells
/// represented by a box filled with a color CI function of the cell content.
/// if the cell content is N, the color CI used will be the color number
/// in colors[N],etc. If the maximum cell content is > ncolors, all
/// cell contents are scaled to ncolors.
///
/// `if ncolors <= 0` a default palette (see below) of 50 colors is
/// defined. The colors defined in this palette are OK for coloring pads, labels.
///
/// ~~~ {.cpp}
/// index 0->9   : grey colors from light to dark grey
/// index 10->19 : "brown" colors
/// index 20->29 : "blueish" colors
/// index 30->39 : "redish" colors
/// index 40->49 : basic colors
/// ~~~
///
/// `if ncolors == 1 && colors == 0`, a Rainbow Color map is created
/// with 50 colors. It is kept for backward compatibility. Better palettes like
/// kBird are recommended.
///
/// High quality predefined palettes with 255 colors are available when `colors == 0`.
/// The following value of `ncolors` give access to:
///
/// ~~~ {.cpp}
/// if ncolors = 51 and colors=0, a Deep Sea palette is used.
/// if ncolors = 52 and colors=0, a Grey Scale palette is used.
/// if ncolors = 53 and colors=0, a Dark Body Radiator palette is used.
/// if ncolors = 54 and colors=0, a Two-Color Hue palette is used.(dark blue through neutral gray to bright yellow)
/// if ncolors = 55 and colors=0, a Rain Bow palette is used.
/// if ncolors = 56 and colors=0, an Inverted Dark Body Radiator palette is used.
/// if ncolors = 57 and colors=0, a monotonically increasing L value palette is used.
/// if ncolors = 58 and colors=0, a Cubehelix palette is used
///                                 (Cf. Dave Green's "cubehelix" colour scheme at http://www.mrao.cam.ac.uk/~dag/CUBEHELIX/)
/// if ncolors = 59 and colors=0, a Green Red Violet palette is used.
/// if ncolors = 60 and colors=0, a Blue Red Yellow palette is used.
/// if ncolors = 61 and colors=0, an Ocean palette is used.
/// if ncolors = 62 and colors=0, a Color Printable On Grey palette is used.
/// if ncolors = 63 and colors=0, an Alpine palette is used.
/// if ncolors = 64 and colors=0, an Aquamarine palette is used.
/// if ncolors = 65 and colors=0, an Army palette is used.
/// if ncolors = 66 and colors=0, an Atlantic palette is used.
/// if ncolors = 67 and colors=0, an Aurora palette is used.
/// if ncolors = 68 and colors=0, an Avocado palette is used.
/// if ncolors = 69 and colors=0, a Beach palette is used.
/// if ncolors = 70 and colors=0, a Black Body palette is used.
/// if ncolors = 71 and colors=0, a Blue Green Yellow palette is used.
/// if ncolors = 72 and colors=0, a Brown Cyan palette is used.
/// if ncolors = 73 and colors=0, a CMYK palette is used.
/// if ncolors = 74 and colors=0, a Candy palette is used.
/// if ncolors = 75 and colors=0, a Cherry palette is used.
/// if ncolors = 76 and colors=0, a Coffee palette is used.
/// if ncolors = 77 and colors=0, a Dark Rain Bow palette is used.
/// if ncolors = 78 and colors=0, a Dark Terrain palette is used.
/// if ncolors = 79 and colors=0, a Fall palette is used.
/// if ncolors = 80 and colors=0, a Fruit Punch palette is used.
/// if ncolors = 81 and colors=0, a Fuchsia palette is used.
/// if ncolors = 82 and colors=0, a Grey Yellow palette is used.
/// if ncolors = 83 and colors=0, a Green Brown Terrain palette is used.
/// if ncolors = 84 and colors=0, a Green Pink palette is used.
/// if ncolors = 85 and colors=0, an Island palette is used.
/// if ncolors = 86 and colors=0, a Lake palette is used.
/// if ncolors = 87 and colors=0, a Light Temperature palette is used.
/// if ncolors = 88 and colors=0, a Light Terrain palette is used.
/// if ncolors = 89 and colors=0, a Mint palette is used.
/// if ncolors = 90 and colors=0, a Neon palette is used.
/// if ncolors = 91 and colors=0, a Pastel palette is used.
/// if ncolors = 92 and colors=0, a Pearl palette is used.
/// if ncolors = 93 and colors=0, a Pigeon palette is used.
/// if ncolors = 94 and colors=0, a Plum palette is used.
/// if ncolors = 95 and colors=0, a Red Blue palette is used.
/// if ncolors = 96 and colors=0, a Rose palette is used.
/// if ncolors = 97 and colors=0, a Rust palette is used.
/// if ncolors = 98 and colors=0, a Sandy Terrain palette is used.
/// if ncolors = 99 and colors=0, a Sienna palette is used.
/// if ncolors = 100 and colors=0, a Solar palette is used.
/// if ncolors = 101 and colors=0, a South West palette is used.
/// if ncolors = 102 and colors=0, a Starry Night palette is used.
/// if ncolors = 103 and colors=0, a Sunset palette is used.
/// if ncolors = 104 and colors=0, a Temperature Map palette is used.
/// if ncolors = 105 and colors=0, a Thermometer palette is used.
/// if ncolors = 106 and colors=0, a Valentine palette is used.
/// if ncolors = 107 and colors=0, a Visible Spectrum palette is used.
/// if ncolors = 108 and colors=0, a Water Melon palette is used.
/// if ncolors = 109 and colors=0, a Cool palette is used.
/// if ncolors = 110 and colors=0, a Copper palette is used.
/// if ncolors = 111 and colors=0, a Gist Earth palette is used.
/// if ncolors = 112 and colors=0, a Viridis palette is used.
/// if ncolors = 113 and colors=0, a Cividis palette is used.
/// ~~~
/// These palettes can also be accessed by names:
/// ~~~ {.cpp}
/// kDeepSea=51,          kGreyScale=52,    kDarkBodyRadiator=53,
/// kBlueYellow= 54,      kRainBow=55,      kInvertedDarkBodyRadiator=56,
/// kBird=57,             kCubehelix=58,    kGreenRedViolet=59,
/// kBlueRedYellow=60,    kOcean=61,        kColorPrintableOnGrey=62,
/// kAlpine=63,           kAquamarine=64,   kArmy=65,
/// kAtlantic=66,         kAurora=67,       kAvocado=68,
/// kBeach=69,            kBlackBody=70,    kBlueGreenYellow=71,
/// kBrownCyan=72,        kCMYK=73,         kCandy=74,
/// kCherry=75,           kCoffee=76,       kDarkRainBow=77,
/// kDarkTerrain=78,      kFall=79,         kFruitPunch=80,
/// kFuchsia=81,          kGreyYellow=82,   kGreenBrownTerrain=83,
/// kGreenPink=84,        kIsland=85,       kLake=86,
/// kLightTemperature=87, kLightTerrain=88, kMint=89,
/// kNeon=90,             kPastel=91,       kPearl=92,
/// kPigeon=93,           kPlum=94,         kRedBlue=95,
/// kRose=96,             kRust=97,         kSandyTerrain=98,
/// kSienna=99,           kSolar=100,       kSouthWest=101,
/// kStarryNight=102,     kSunset=103,      kTemperatureMap=104,
/// kThermometer=105,     kValentine=106,   kVisibleSpectrum=107,
/// kWaterMelon=108,      kCool=109,        kCopper=110,
/// kGistEarth=111        kViridis=112,     kCividis=113
/// ~~~
/// For example:
/// ~~~ {.cpp}
/// gStyle->SetPalette(kBird);
/// ~~~
/// Set the current palette as "Bird" (number 57).
///
/// The color numbers specified in the palette can be viewed by selecting
/// the item "colors" in the "VIEW" menu of the canvas toolbar.
/// The color parameters can be changed via TColor::SetRGB.
///
/// Note that when drawing a 2D histogram `h2` with the option "COL" or
/// "COLZ" or with any "CONT" options using the color map, the number of colors
/// used is defined by the number of contours `n` specified with:
/// `h2->SetContour(n)`

void TColor::SetPalette(Int_t ncolors, Int_t *colors, Float_t alpha)
{
   Int_t i;

   static Int_t paletteType = 0;

   Int_t palette[50] = {19,18,17,16,15,14,13,12,11,20,
                        21,22,23,24,25,26,27,28,29,30, 8,
                        31,32,33,34,35,36,37,38,39,40, 9,
                        41,42,43,44,45,47,48,49,46,50, 2,
                         7, 6, 5, 4, 3, 2,1};

   // set default palette (pad type)
   if (ncolors <= 0) {
      ncolors = 50;
      fgPalette.Set(ncolors);
      for (i=0;i<ncolors;i++) fgPalette.fArray[i] = palette[i];
      paletteType = 1;
      return;
   }

   // set Rainbow Color map. Kept for backward compatibility.
   if (ncolors == 1 && colors == 0) {
      ncolors = 50;
      fgPalette.Set(ncolors);
      for (i=0;i<ncolors-1;i++) fgPalette.fArray[i] = 51+i;
      fgPalette.fArray[ncolors-1] = kRed; // the last color of this palette is red
      paletteType = 2;
      return;
   }

   // High quality palettes (255 levels)
   if (colors == 0 && ncolors>50) {

      if (!fgPalettesList.fN) fgPalettesList.Set(63);        // Right now 63 high quality palettes
      Int_t Idx = (Int_t)fgPalettesList.fArray[ncolors-51];  // High quality palettes indices start at 51

      // This high quality palette has already been created. Reuse it.
      if (Idx > 0) {
         Double_t alphas = 10*(fgPalettesList.fArray[ncolors-51]-Idx);
         Bool_t same_alpha = TMath::Abs(alpha-alphas) < 0.0001;
         if (paletteType == ncolors && same_alpha) return; // The current palette is already this one.
         fgPalette.Set(255); // High quality palettes have 255 entries
         for (i=0;i<255;i++) fgPalette.fArray[i] = Idx+i;
         paletteType = ncolors;

         // restore the palette transparency if needed
          if (alphas>0 && !same_alpha) {
             TColor *ca;
             for (i=0;i<255;i++) {
                ca = gROOT->GetColor(Idx+i);
                ca->SetAlpha(alpha);
             }
             fgPalettesList.fArray[paletteType-51] = (Double_t)Idx+alpha/10.;
          }
         return;
      }

      TColor::InitializeColors();
      Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};

      switch (ncolors) {
      // Deep Sea
      case 51:
         {
            Double_t red[9]   = {  0./255.,  9./255., 13./255., 17./255., 24./255.,  32./255.,  27./255.,  25./255.,  29./255.};
            Double_t green[9] = {  0./255.,  0./255.,  0./255.,  2./255., 37./255.,  74./255., 113./255., 160./255., 221./255.};
            Double_t blue[9]  = { 28./255., 42./255., 59./255., 78./255., 98./255., 129./255., 154./255., 184./255., 221./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Grey Scale
      case 52:
         {
            Double_t red[9]   = { 0./255., 32./255., 64./255., 96./255., 128./255., 160./255., 192./255., 224./255., 255./255.};
            Double_t green[9] = { 0./255., 32./255., 64./255., 96./255., 128./255., 160./255., 192./255., 224./255., 255./255.};
            Double_t blue[9]  = { 0./255., 32./255., 64./255., 96./255., 128./255., 160./255., 192./255., 224./255., 255./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Dark Body Radiator
      case 53:
         {
            Double_t red[9]   = { 0./255., 45./255., 99./255., 156./255., 212./255., 230./255., 237./255., 234./255., 242./255.};
            Double_t green[9] = { 0./255.,  0./255.,  0./255.,  45./255., 101./255., 168./255., 238./255., 238./255., 243./255.};
            Double_t blue[9]  = { 0./255.,  1./255.,  1./255.,   3./255.,   9./255.,   8./255.,  11./255.,  95./255., 230./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Two-color hue (dark blue through neutral gray to bright yellow)
      case 54:
         {
            Double_t red[9]   = {  0./255.,  22./255., 44./255., 68./255., 93./255., 124./255., 160./255., 192./255., 237./255.};
            Double_t green[9] = {  0./255.,  16./255., 41./255., 67./255., 93./255., 125./255., 162./255., 194./255., 241./255.};
            Double_t blue[9]  = { 97./255., 100./255., 99./255., 99./255., 93./255.,  68./255.,  44./255.,  26./255.,  74./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Rain Bow
      case 55:
         {
            Double_t red[9]   = {  0./255.,   5./255.,  15./255.,  35./255., 102./255., 196./255., 208./255., 199./255., 110./255.};
            Double_t green[9] = {  0./255.,  48./255., 124./255., 192./255., 206./255., 226./255.,  97./255.,  16./255.,   0./255.};
            Double_t blue[9]  = { 99./255., 142./255., 198./255., 201./255.,  90./255.,  22./255.,  13./255.,   8./255.,   2./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Inverted Dark Body Radiator
      case 56:
         {
            Double_t red[9]   = { 242./255., 234./255., 237./255., 230./255., 212./255., 156./255., 99./255., 45./255., 0./255.};
            Double_t green[9] = { 243./255., 238./255., 238./255., 168./255., 101./255.,  45./255.,  0./255.,  0./255., 0./255.};
            Double_t blue[9]  = { 230./255.,  95./255.,  11./255.,   8./255.,   9./255.,   3./255.,  1./255.,  1./255., 0./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Bird
      case 57:
         {
            Double_t red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
            Double_t green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
            Double_t blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Cubehelix
      case 58:
         {
            Double_t red[9]   = { 0.0000, 0.0956, 0.0098, 0.2124, 0.6905, 0.9242, 0.7914, 0.7596, 1.0000};
            Double_t green[9] = { 0.0000, 0.1147, 0.3616, 0.5041, 0.4577, 0.4691, 0.6905, 0.9237, 1.0000};
            Double_t blue[9]  = { 0.0000, 0.2669, 0.3121, 0.1318, 0.2236, 0.6741, 0.9882, 0.9593, 1.0000};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Green Red Violet
      case 59:
         {
            Double_t red[9]   = {13./255., 23./255., 25./255., 63./255., 76./255., 104./255., 137./255., 161./255., 206./255.};
            Double_t green[9] = {95./255., 67./255., 37./255., 21./255.,  0./255.,  12./255.,  35./255.,  52./255.,  79./255.};
            Double_t blue[9]  = { 4./255.,  3./255.,  2./255.,  6./255., 11./255.,  22./255.,  49./255.,  98./255., 208./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Blue Red Yellow
      case 60:
         {
            Double_t red[9]   = {0./255.,  61./255.,  89./255., 122./255., 143./255., 160./255., 185./255., 204./255., 231./255.};
            Double_t green[9] = {0./255.,   0./255.,   0./255.,   0./255.,  14./255.,  37./255.,  72./255., 132./255., 235./255.};
            Double_t blue[9]  = {0./255., 140./255., 224./255., 144./255.,   4./255.,   5./255.,   6./255.,   9./255.,  13./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Ocean
      case 61:
         {
            Double_t red[9]   = { 14./255.,  7./255.,  2./255.,  0./255.,  5./255.,  11./255.,  55./255., 131./255., 229./255.};
            Double_t green[9] = {105./255., 56./255., 26./255.,  1./255., 42./255.,  74./255., 131./255., 171./255., 229./255.};
            Double_t blue[9]  = {  2./255., 21./255., 35./255., 60./255., 92./255., 113./255., 160./255., 185./255., 229./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Color Printable On Grey
      case 62:
         {
            Double_t red[9]   = { 0./255.,   0./255.,   0./255.,  70./255., 148./255., 231./255., 235./255., 237./255., 244./255.};
            Double_t green[9] = { 0./255.,   0./255.,   0./255.,   0./255.,   0./255.,  69./255.,  67./255., 216./255., 244./255.};
            Double_t blue[9]  = { 0./255., 102./255., 228./255., 231./255., 177./255., 124./255., 137./255.,  20./255., 244./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Alpine
      case 63:
         {
            Double_t red[9]   = { 50./255., 56./255., 63./255., 68./255.,  93./255., 121./255., 165./255., 192./255., 241./255.};
            Double_t green[9] = { 66./255., 81./255., 91./255., 96./255., 111./255., 128./255., 155./255., 189./255., 241./255.};
            Double_t blue[9]  = { 97./255., 91./255., 75./255., 65./255.,  77./255., 103./255., 143./255., 167./255., 217./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Aquamarine
      case 64:
         {
            Double_t red[9]   = { 145./255., 166./255., 167./255., 156./255., 131./255., 114./255., 101./255., 112./255., 132./255.};
            Double_t green[9] = { 158./255., 178./255., 179./255., 181./255., 163./255., 154./255., 144./255., 152./255., 159./255.};
            Double_t blue[9]  = { 190./255., 199./255., 201./255., 192./255., 176./255., 169./255., 160./255., 166./255., 190./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Army
      case 65:
         {
            Double_t red[9]   = { 93./255.,   91./255.,  99./255., 108./255., 130./255., 125./255., 132./255., 155./255., 174./255.};
            Double_t green[9] = { 126./255., 124./255., 128./255., 129./255., 131./255., 121./255., 119./255., 153./255., 173./255.};
            Double_t blue[9]  = { 103./255.,  94./255.,  87./255.,  85./255.,  80./255.,  85./255., 107./255., 120./255., 146./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Atlantic
      case 66:
         {
            Double_t red[9]   = { 24./255., 40./255., 69./255.,  90./255., 104./255., 114./255., 120./255., 132./255., 103./255.};
            Double_t green[9] = { 29./255., 52./255., 94./255., 127./255., 150./255., 162./255., 159./255., 151./255., 101./255.};
            Double_t blue[9]  = { 29./255., 52./255., 96./255., 132./255., 162./255., 181./255., 184./255., 186./255., 131./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Aurora
      case 67:
         {
            Double_t red[9]   = { 46./255., 38./255., 61./255., 92./255., 113./255., 121./255., 132./255., 150./255., 191./255.};
            Double_t green[9] = { 46./255., 36./255., 40./255., 69./255., 110./255., 135./255., 131./255.,  92./255.,  34./255.};
            Double_t blue[9]  = { 46./255., 80./255., 74./255., 70./255.,  81./255., 105./255., 165./255., 211./255., 225./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Avocado
      case 68:
         {
            Double_t red[9]   = { 0./255.,  4./255., 12./255.,  30./255.,  52./255., 101./255., 142./255., 190./255., 237./255.};
            Double_t green[9] = { 0./255., 40./255., 86./255., 121./255., 140./255., 172./255., 187./255., 213./255., 240./255.};
            Double_t blue[9]  = { 0./255.,  9./255., 14./255.,  18./255.,  21./255.,  23./255.,  27./255.,  35./255., 101./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Beach
      case 69:
         {
            Double_t red[9]   = { 198./255., 206./255., 206./255., 211./255., 198./255., 181./255., 161./255., 171./255., 244./255.};
            Double_t green[9] = { 103./255., 133./255., 150./255., 172./255., 178./255., 174./255., 163./255., 175./255., 244./255.};
            Double_t blue[9]  = {  49./255.,  54./255.,  55./255.,  66./255.,  91./255., 130./255., 184./255., 224./255., 244./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Black Body
      case 70:
         {
            Double_t red[9]   = { 243./255., 243./255., 240./255., 240./255., 241./255., 239./255., 186./255., 151./255., 129./255.};
            Double_t green[9] = {   0./255.,  46./255.,  99./255., 149./255., 194./255., 220./255., 183./255., 166./255., 147./255.};
            Double_t blue[9]  = {   6./255.,   8./255.,  36./255.,  91./255., 169./255., 235./255., 246./255., 240./255., 233./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Blue Green Yellow
      case 71:
         {
            Double_t red[9]   = { 22./255., 19./255.,  19./255.,  25./255.,  35./255.,  53./255.,  88./255., 139./255., 210./255.};
            Double_t green[9] = {  0./255., 32./255.,  69./255., 108./255., 135./255., 159./255., 183./255., 198./255., 215./255.};
            Double_t blue[9]  = { 77./255., 96./255., 110./255., 116./255., 110./255., 100./255.,  90./255.,  78./255.,  70./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Brown Cyan
      case 72:
         {
            Double_t red[9]   = { 68./255., 116./255., 165./255., 182./255., 189./255., 180./255., 145./255., 111./255.,  71./255.};
            Double_t green[9] = { 37./255.,  82./255., 135./255., 178./255., 204./255., 225./255., 221./255., 202./255., 147./255.};
            Double_t blue[9]  = { 16./255.,  55./255., 105./255., 147./255., 196./255., 226./255., 232./255., 224./255., 178./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // CMYK
      case 73:
         {
            Double_t red[9]   = {  61./255.,  99./255., 136./255., 181./255., 213./255., 225./255., 198./255., 136./255., 24./255.};
            Double_t green[9] = { 149./255., 140./255.,  96./255.,  83./255., 132./255., 178./255., 190./255., 135./255., 22./255.};
            Double_t blue[9]  = { 214./255., 203./255., 168./255., 135./255., 110./255., 100./255., 111./255., 113./255., 22./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Candy
      case 74:
         {
            Double_t red[9]   = { 76./255., 120./255., 156./255., 183./255., 197./255., 180./255., 162./255., 154./255., 140./255.};
            Double_t green[9] = { 34./255.,  35./255.,  42./255.,  69./255., 102./255., 137./255., 164./255., 188./255., 197./255.};
            Double_t blue[9]  = { 64./255.,  69./255.,  78./255., 105./255., 142./255., 177./255., 205./255., 217./255., 198./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Cherry
      case 75:
         {
            Double_t red[9]   = { 37./255., 102./255., 157./255., 188./255., 196./255., 214./255., 223./255., 235./255., 251./255.};
            Double_t green[9] = { 37./255.,  29./255.,  25./255.,  37./255.,  67./255.,  91./255., 132./255., 185./255., 251./255.};
            Double_t blue[9]  = { 37./255.,  32./255.,  33./255.,  45./255.,  66./255.,  98./255., 137./255., 187./255., 251./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Coffee
      case 76:
         {
            Double_t red[9]   = { 79./255., 100./255., 119./255., 137./255., 153./255., 172./255., 192./255., 205./255., 250./255.};
            Double_t green[9] = { 63./255.,  79./255.,  93./255., 103./255., 115./255., 135./255., 167./255., 196./255., 250./255.};
            Double_t blue[9]  = { 51./255.,  59./255.,  66./255.,  61./255.,  62./255.,  70./255., 110./255., 160./255., 250./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Dark Rain Bow
      case 77:
         {
            Double_t red[9]   = {  43./255.,  44./255., 50./255.,  66./255., 125./255., 172./255., 178./255., 155./255., 157./255.};
            Double_t green[9] = {  63./255.,  63./255., 85./255., 101./255., 138./255., 163./255., 122./255.,  51./255.,  39./255.};
            Double_t blue[9]  = { 121./255., 101./255., 58./255.,  44./255.,  47./255.,  55./255.,  57./255.,  44./255.,  43./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Dark Terrain
      case 78:
         {
            Double_t red[9]   = {  0./255., 41./255., 62./255., 79./255., 90./255., 87./255., 99./255., 140./255., 228./255.};
            Double_t green[9] = {  0./255., 57./255., 81./255., 93./255., 85./255., 70./255., 71./255., 125./255., 228./255.};
            Double_t blue[9]  = { 95./255., 91./255., 91./255., 82./255., 60./255., 43./255., 44./255., 112./255., 228./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Fall
      case 79:
         {
            Double_t red[9]   = { 49./255., 59./255., 72./255., 88./255., 114./255., 141./255., 176./255., 205./255., 222./255.};
            Double_t green[9] = { 78./255., 72./255., 66./255., 57./255.,  59./255.,  75./255., 106./255., 142./255., 173./255.};
            Double_t blue[9]  = { 78./255., 55./255., 46./255., 40./255.,  39./255.,  39./255.,  40./255.,  41./255.,  47./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Fruit Punch
      case 80:
         {
            Double_t red[9]   = { 243./255., 222./255., 201./255., 185./255., 165./255., 158./255., 166./255., 187./255., 219./255.};
            Double_t green[9] = {  94./255., 108./255., 132./255., 135./255., 125./255.,  96./255.,  68./255.,  51./255.,  61./255.};
            Double_t blue[9]  = {   7./255.,  9./255.,   12./255.,  19./255.,  45./255.,  89./255., 118./255., 146./255., 118./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Fuchsia
      case 81:
         {
            Double_t red[9]   = { 19./255., 44./255., 74./255., 105./255., 137./255., 166./255., 194./255., 206./255., 220./255.};
            Double_t green[9] = { 19./255., 28./255., 40./255.,  55./255.,  82./255., 110./255., 159./255., 181./255., 220./255.};
            Double_t blue[9]  = { 19./255., 42./255., 68./255.,  96./255., 129./255., 157./255., 188./255., 203./255., 220./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Grey Yellow
      case 82:
         {
            Double_t red[9]   = { 33./255., 44./255., 70./255.,  99./255., 140./255., 165./255., 199./255., 211./255., 216./255.};
            Double_t green[9] = { 38./255., 50./255., 76./255., 105./255., 140./255., 165./255., 191./255., 189./255., 167./255.};
            Double_t blue[9]  = { 55./255., 67./255., 97./255., 124./255., 140./255., 166./255., 163./255., 129./255.,  52./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Green Brown Terrain
      case 83:
         {
            Double_t red[9]   = { 0./255., 33./255., 73./255., 124./255., 136./255., 152./255., 159./255., 171./255., 223./255.};
            Double_t green[9] = { 0./255., 43./255., 92./255., 124./255., 134./255., 126./255., 121./255., 144./255., 223./255.};
            Double_t blue[9]  = { 0./255., 43./255., 68./255.,  76./255.,  73./255.,  64./255.,  72./255., 114./255., 223./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Green Pink
      case 84:
         {
            Double_t red[9]   = {  5./255.,  18./255.,  45./255., 124./255., 193./255., 223./255., 205./255., 128./255., 49./255.};
            Double_t green[9] = { 48./255., 134./255., 207./255., 230./255., 193./255., 113./255.,  28./255.,   0./255.,  7./255.};
            Double_t blue[9]  = {  6./255.,  15./255.,  41./255., 121./255., 193./255., 226./255., 208./255., 130./255., 49./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Island
      case 85:
         {
            Double_t red[9]   = { 180./255., 106./255., 104./255., 135./255., 164./255., 188./255., 189./255., 165./255., 144./255.};
            Double_t green[9] = {  72./255., 126./255., 154./255., 184./255., 198./255., 207./255., 205./255., 190./255., 179./255.};
            Double_t blue[9]  = {  41./255., 120./255., 158./255., 188./255., 194./255., 181./255., 145./255., 100./255.,  62./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Lake
      case 86:
         {
            Double_t red[9]   = {  57./255.,  72./255.,  94./255., 117./255., 136./255., 154./255., 174./255., 192./255., 215./255.};
            Double_t green[9] = {   0./255.,  33./255.,  68./255., 109./255., 140./255., 171./255., 192./255., 196./255., 209./255.};
            Double_t blue[9]  = { 116./255., 137./255., 173./255., 201./255., 200./255., 201./255., 203./255., 190./255., 187./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Light Temperature
      case 87:
         {
            Double_t red[9]   = {  31./255.,  71./255., 123./255., 160./255., 210./255., 222./255., 214./255., 199./255., 183./255.};
            Double_t green[9] = {  40./255., 117./255., 171./255., 211./255., 231./255., 220./255., 190./255., 132./255.,  65./255.};
            Double_t blue[9]  = { 234./255., 214./255., 228./255., 222./255., 210./255., 160./255., 105./255.,  60./255.,  34./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Light Terrain
      case 88:
         {
            Double_t red[9]   = { 123./255., 108./255., 109./255., 126./255., 154./255., 172./255., 188./255., 196./255., 218./255.};
            Double_t green[9] = { 184./255., 138./255., 130./255., 133./255., 154./255., 175./255., 188./255., 196./255., 218./255.};
            Double_t blue[9]  = { 208./255., 130./255., 109./255.,  99./255., 110./255., 122./255., 150./255., 171./255., 218./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Mint
      case 89:
         {
            Double_t red[9]   = { 105./255., 106./255., 122./255., 143./255., 159./255., 172./255., 176./255., 181./255., 207./255.};
            Double_t green[9] = { 252./255., 197./255., 194./255., 187./255., 174./255., 162./255., 153./255., 136./255., 125./255.};
            Double_t blue[9]  = { 146./255., 133./255., 144./255., 155./255., 163./255., 167./255., 166./255., 162./255., 174./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Neon
      case 90:
         {
            Double_t red[9]   = { 171./255., 141./255., 145./255., 152./255., 154./255., 159./255., 163./255., 158./255., 177./255.};
            Double_t green[9] = { 236./255., 143./255., 100./255.,  63./255.,  53./255.,  55./255.,  44./255.,  31./255.,   6./255.};
            Double_t blue[9]  = {  59./255.,  48./255.,  46./255.,  44./255.,  42./255.,  54./255.,  82./255., 112./255., 179./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Pastel
      case 91:
         {
            Double_t red[9]   = { 180./255., 190./255., 209./255., 223./255., 204./255., 228./255., 205./255., 152./255.,  91./255.};
            Double_t green[9] = {  93./255., 125./255., 147./255., 172./255., 181./255., 224./255., 233./255., 198./255., 158./255.};
            Double_t blue[9]  = { 236./255., 218./255., 160./255., 133./255., 114./255., 132./255., 162./255., 220./255., 218./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Pearl
      case 92:
         {
            Double_t red[9]   = { 225./255., 183./255., 162./255., 135./255., 115./255., 111./255., 119./255., 145./255., 211./255.};
            Double_t green[9] = { 205./255., 177./255., 166./255., 135./255., 124./255., 117./255., 117./255., 132./255., 172./255.};
            Double_t blue[9]  = { 186./255., 165./255., 155./255., 135./255., 126./255., 130./255., 150./255., 178./255., 226./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Pigeon
      case 93:
         {
            Double_t red[9]   = { 39./255., 43./255., 59./255., 63./255., 80./255., 116./255., 153./255., 177./255., 223./255.};
            Double_t green[9] = { 39./255., 43./255., 59./255., 74./255., 91./255., 114./255., 139./255., 165./255., 223./255.};
            Double_t blue[9]  = { 39./255., 50./255., 59./255., 70./255., 85./255., 115./255., 151./255., 176./255., 223./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Plum
      case 94:
         {
            Double_t red[9]   = { 0./255., 38./255., 60./255., 76./255., 84./255., 89./255., 101./255., 128./255., 204./255.};
            Double_t green[9] = { 0./255., 10./255., 15./255., 23./255., 35./255., 57./255.,  83./255., 123./255., 199./255.};
            Double_t blue[9]  = { 0./255., 11./255., 22./255., 40./255., 63./255., 86./255.,  97./255.,  94./255.,  85./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Red Blue
      case 95:
         {
            Double_t red[9]   = { 94./255., 112./255., 141./255., 165./255., 167./255., 140./255.,  91./255.,  49./255.,  27./255.};
            Double_t green[9] = { 27./255.,  46./255.,  88./255., 135./255., 166./255., 161./255., 135./255.,  97./255.,  58./255.};
            Double_t blue[9]  = { 42./255.,  52./255.,  81./255., 106./255., 139./255., 158./255., 155./255., 137./255., 116./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Rose
      case 96:
         {
            Double_t red[9]   = { 30./255., 49./255., 79./255., 117./255., 135./255., 151./255., 146./255., 138./255., 147./255.};
            Double_t green[9] = { 63./255., 60./255., 72./255.,  90./255.,  94./255.,  94./255.,  68./255.,  46./255.,  16./255.};
            Double_t blue[9]  = { 18./255., 28./255., 41./255.,  56./255.,  62./255.,  63./255.,  50./255.,  36./255.,  21./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Rust
      case 97:
         {
            Double_t red[9]   = {  0./255., 30./255., 63./255., 101./255., 143./255., 152./255., 169./255., 187./255., 230./255.};
            Double_t green[9] = {  0./255., 14./255., 28./255.,  42./255.,  58./255.,  61./255.,  67./255.,  74./255.,  91./255.};
            Double_t blue[9]  = { 39./255., 26./255., 21./255.,  18./255.,  15./255.,  14./255.,  14./255.,  13./255.,  13./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Sandy Terrain
      case 98:
         {
            Double_t red[9]   = { 149./255., 140./255., 164./255., 179./255., 182./255., 181./255., 131./255., 87./255., 61./255.};
            Double_t green[9] = {  62./255.,  70./255., 107./255., 136./255., 144./255., 138./255., 117./255., 87./255., 74./255.};
            Double_t blue[9]  = {  40./255.,  38./255.,  45./255.,  49./255.,  49./255.,  49./255.,  38./255., 32./255., 34./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Sienna
      case 99:
         {
            Double_t red[9]   = { 99./255., 112./255., 148./255., 165./255., 179./255., 182./255., 183./255., 183./255., 208./255.};
            Double_t green[9] = { 39./255.,  40./255.,  57./255.,  79./255., 104./255., 127./255., 148./255., 161./255., 198./255.};
            Double_t blue[9]  = { 15./255.,  16./255.,  18./255.,  33./255.,  51./255.,  79./255., 103./255., 129./255., 177./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Solar
      case 100:
         {
            Double_t red[9]   = { 99./255., 116./255., 154./255., 174./255., 200./255., 196./255., 201./255., 201./255., 230./255.};
            Double_t green[9] = {  0./255.,   0./255.,   8./255.,  32./255.,  58./255.,  83./255., 119./255., 136./255., 173./255.};
            Double_t blue[9]  = {  5./255.,   6./255.,   7./255.,   9./255.,   9./255.,  14./255.,  17./255.,  19./255.,  24./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // South West
      case 101:
         {
            Double_t red[9]   = { 82./255., 106./255., 126./255., 141./255., 155./255., 163./255., 142./255., 107./255.,  66./255.};
            Double_t green[9] = { 62./255.,  44./255.,  69./255., 107./255., 135./255., 152./255., 149./255., 132./255., 119./255.};
            Double_t blue[9]  = { 39./255.,  25./255.,  31./255.,  60./255.,  73./255.,  68./255.,  49./255.,  72./255., 188./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Starry Night
      case 102:
         {
            Double_t red[9]   = { 18./255., 29./255., 44./255.,  72./255., 116./255., 158./255., 184./255., 208./255., 221./255.};
            Double_t green[9] = { 27./255., 46./255., 71./255., 105./255., 146./255., 177./255., 189./255., 190./255., 183./255.};
            Double_t blue[9]  = { 39./255., 55./255., 80./255., 108./255., 130./255., 133./255., 124./255., 100./255.,  76./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Sunset
      case 103:
         {
            Double_t red[9]   = { 0./255., 48./255., 119./255., 173./255., 212./255., 224./255., 228./255., 228./255., 245./255.};
            Double_t green[9] = { 0./255., 13./255.,  30./255.,  47./255.,  79./255., 127./255., 167./255., 205./255., 245./255.};
            Double_t blue[9]  = { 0./255., 68./255.,  75./255.,  43./255.,  16./255.,  22./255.,  55./255., 128./255., 245./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Temperature Map
      case 104:
         {
            Double_t red[9]   = {  34./255.,  70./255., 129./255., 187./255., 225./255., 226./255., 216./255., 193./255., 179./255.};
            Double_t green[9] = {  48./255.,  91./255., 147./255., 194./255., 226./255., 229./255., 196./255., 110./255.,  12./255.};
            Double_t blue[9]  = { 234./255., 212./255., 216./255., 224./255., 206./255., 110./255.,  53./255.,  40./255.,  29./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Thermometer
      case 105:
         {
            Double_t red[9]   = {  30./255.,  55./255., 103./255., 147./255., 174./255., 203./255., 188./255., 151./255., 105./255.};
            Double_t green[9] = {   0./255.,  65./255., 138./255., 182./255., 187./255., 175./255., 121./255.,  53./255.,   9./255.};
            Double_t blue[9]  = { 191./255., 202./255., 212./255., 208./255., 171./255., 140./255.,  97./255.,  57./255.,  30./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Valentine
      case 106:
         {
            Double_t red[9]   = { 112./255., 97./255., 113./255., 125./255., 138./255., 159./255., 178./255., 188./255., 225./255.};
            Double_t green[9] = {  16./255., 17./255.,  24./255.,  37./255.,  56./255.,  81./255., 110./255., 136./255., 189./255.};
            Double_t blue[9]  = {  38./255., 35./255.,  46./255.,  59./255.,  78./255., 103./255., 130./255., 152./255., 201./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Visible Spectrum
      case 107:
         {
            Double_t red[9]   = { 18./255.,  72./255.,   5./255.,  23./255.,  29./255., 201./255., 200./255., 98./255., 29./255.};
            Double_t green[9] = {  0./255.,   0./255.,  43./255., 167./255., 211./255., 117./255.,   0./255.,  0./255.,  0./255.};
            Double_t blue[9]  = { 51./255., 203./255., 177./255.,  26./255.,  10./255.,   9./255.,   8./255.,  3./255.,  0./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Water Melon
      case 108:
         {
            Double_t red[9]   = { 19./255., 42./255., 64./255.,  88./255., 118./255., 147./255., 175./255., 187./255., 205./255.};
            Double_t green[9] = { 19./255., 55./255., 89./255., 125./255., 154./255., 169./255., 161./255., 129./255.,  70./255.};
            Double_t blue[9]  = { 19./255., 32./255., 47./255.,  70./255., 100./255., 128./255., 145./255., 130./255.,  75./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Cool
      case 109:
         {
            Double_t red[9]   = {  33./255.,  31./255.,  42./255.,  68./255.,  86./255., 111./255., 141./255., 172./255., 227./255.};
            Double_t green[9] = { 255./255., 175./255., 145./255., 106./255.,  88./255.,  55./255.,  15./255.,   0./255.,   0./255.};
            Double_t blue[9]  = { 255./255., 205./255., 202./255., 203./255., 208./255., 205./255., 203./255., 206./255., 231./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Copper
      case 110:
         {
            Double_t red[9]   = { 0./255., 25./255., 50./255., 79./255., 110./255., 145./255., 181./255., 201./255., 254./255.};
            Double_t green[9] = { 0./255., 16./255., 30./255., 46./255.,  63./255.,  82./255., 101./255., 124./255., 179./255.};
            Double_t blue[9]  = { 0./255., 12./255., 21./255., 29./255.,  39./255.,  49./255.,  61./255.,  74./255., 103./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Gist Earth
      case 111:
         {
            Double_t red[9]   = { 0./255., 13./255.,  30./255.,  44./255.,  72./255., 120./255., 156./255., 200./255., 247./255.};
            Double_t green[9] = { 0./255., 36./255.,  84./255., 117./255., 141./255., 153./255., 151./255., 158./255., 247./255.};
            Double_t blue[9]  = { 0./255., 94./255., 100./255.,  82./255.,  56./255.,  66./255.,  76./255., 131./255., 247./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Viridis
      case 112:
         {
            Double_t red[9]   = { 26./255., 51./255.,  43./255.,  33./255.,  28./255.,  35./255.,  74./255., 144./255., 246./255.};
            Double_t green[9] = {  9./255., 24./255.,  55./255.,  87./255., 118./255., 150./255., 180./255., 200./255., 222./255.};
            Double_t blue[9]  = { 30./255., 96./255., 112./255., 114./255., 112./255., 101./255.,  72./255.,  35./255.,   0./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Cividis
      case 113:
         {
            Double_t red[9]   = {  0./255.,   5./255.,  65./255.,  97./255., 124./255., 156./255., 189./255., 224./255., 255./255.};
            Double_t green[9] = { 32./255.,  54./255.,  77./255., 100./255., 123./255., 148./255., 175./255., 203./255., 234./255.};
            Double_t blue[9]  = { 77./255., 110./255., 107./255., 111./255., 120./255., 119./255., 111./255.,  94./255.,  70./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      default:
         ::Error("SetPalette", "Unknown palette number %d", ncolors);
         return;
      }
      paletteType = ncolors;
      if (Idx>0) fgPalettesList.fArray[paletteType-51] = (Double_t)Idx;
      else       fgPalettesList.fArray[paletteType-51] = 0.;
      if (alpha > 0.) fgPalettesList.fArray[paletteType-51] += alpha/10.0f;
      return;
   }

   // set user defined palette
   if (colors)  {
      fgPalette.Set(ncolors);
      for (i=0;i<ncolors;i++) fgPalette.fArray[i] = colors[i];
   } else {
      fgPalette.Set(TMath::Min(50,ncolors));
      for (i=0;i<TMath::Min(50,ncolors);i++) fgPalette.fArray[i] = palette[i];
   }
   paletteType = 3;
}


////////////////////////////////////////////////////////////////////////////////
/// Invert the current color palette.
/// The top of the palette becomes the bottom and vice versa.

void TColor::InvertPalette()
{
   std::reverse(fgPalette.fArray, fgPalette.fArray + fgPalette.GetSize());
}
