// @(#)root/histpainter:$Id$
// Author: Rene Brun, Olivier Couet

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cctype>
#include <iostream>

#include "TROOT.h"
#include "TSystem.h"
#include "THistPainter.h"
#include "TH2.h"
#include "TH2Poly.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "THStack.h"
#include "TF2.h"
#include "TF3.h"
#include "TCutG.h"
#include "TMatrixDBase.h"
#include "TMatrixFBase.h"
#include "TVectorD.h"
#include "TVectorF.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TPaveStats.h"
#include "TFrame.h"
#include "TLatex.h"
#include "TPolyLine.h"
#include "TPoints.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TPie.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TPainter3dAlgorithms.h"
#include "TGraph2DPainter.h"
#include "TGraphDelaunay2D.h"
#include "TView.h"
#include "TMath.h"
#include "TRandom2.h"
#include "TObjArray.h"
#include "Hoption.h"
#include "Hparam.h"
#include "TPluginManager.h"
#include "TPaletteAxis.h"
#include "TCrown.h"
#include "TArrow.h"
#include "TVirtualPadEditor.h"
#include "TVirtualX.h"
#include "TEnv.h"
#include "TPoint.h"
#include "TImage.h"
#include "TCandle.h"
#include "strlcpy.h"

/*! \class THistPainter
\ingroup Histpainter
\brief The histogram painter class. Implements all histograms' drawing's options.

- [Introduction](#HP00)
- [Histograms' plotting options](#HP01)
   - [Options supported for 1D and 2D histograms](#HP01a)
   - [Options supported for 1D histograms](#HP01b)
   - [Options supported for 2D histograms](#HP01c)
   - [Options supported for 3D histograms](#HP01d)
   - [Options supported for histograms' stacks (THStack)](#HP01e)
- [Setting the Style](#HP02)
- [Setting line, fill, marker, and text attributes](#HP03)
- [Setting Tick marks on the histogram axis](#HP04)
- [Giving titles to the X, Y and Z axis](#HP05)
- [The option "SAME"](#HP060)
   - [Limitations](#HP060a)
- [Colors automatically picked in palette](#HP061)
- [Superimposing two histograms with different scales in the same pad](#HP06)
- [Statistics Display](#HP07)
- [Fit Statistics](#HP08)
- [The error bars options](#HP09)
- [The bar chart option](#HP100)
- [The "BAR" and "HBAR" options](#HP10)
- [The SCATter plot option (default for 2D histograms)](#HP11)
- [The ARRow option](#HP12)
- [The BOX option](#HP13)
- [The COLor option](#HP14)
- [The CANDLE and VIOLIN options](#HP140)
   - [The CANDLE option](#HP140a)
   - [The VIOLIN option](#HP140b)
- [The TEXT and TEXTnn Option](#HP15)
- [The CONTour options](#HP16)
   - [The LIST option](#HP16a)
   - [The AITOFF, MERCATOR, SINUSOIDAL and PARABOLIC options](#HP16b)
- [The LEGO options](#HP17)
- [The "SURFace" options](#HP18)
- [Cylindrical, Polar, Spherical and PseudoRapidity/Phi options](#HP19)
- [Base line for bar-charts and lego plots](#HP20)
- [TH2Poly Drawing](#HP20a)
- [The SPEC option](#HP21)
- [Option "Z" : Adding the color palette on the right side of the pad](#HP22)
- [Setting the color palette](#HP23)
- [Drawing a sub-range of a 2-D histogram; the [cutg] option](#HP24)
- [Drawing options for 3D histograms](#HP25)
- [Drawing option for histograms' stacks](#HP26)
- [Drawing of 3D implicit functions](#HP27)
- [Associated functions drawing](#HP28)
- [Drawing using OpenGL](#HP29)
   - [General information: plot types and supported options](#HP29a)
   - [TH3 as color boxes](#HP290)
   - [TH3 as boxes (spheres)](#HP29b)
   - [TH3 as iso-surface(s)](#HP29c)
   - [TF3 (implicit function)](#HP29d)
   - [Parametric surfaces](#HP29e)
   - [Interaction with the plots](#HP29f)
   - [Selectable parts](#HP29g)
   - [Rotation and zooming](#HP29h)
   - [Panning](#HP29i)
   - [Box cut](#HP29j)
   - [Plot specific interactions (dynamic slicing etc.)](#HP29k)
   - [Surface with option "GLSURF"](#HP29l)
   - [TF3](#HP29m)
   - [Box](#HP29n)
   - [Iso](#HP29o)
   - [Parametric plot](#HP29p)
- [Highlight mode for histogram](#HP30)
   - [Highlight mode and user function](#HP30a)


## <a name="HP00"></a> Introduction


Histograms are drawn via the `THistPainter` class. Each histogram has a
pointer to its own painter (to be usable in a multithreaded program). When the
canvas has to be redrawn, the `Paint` function of each objects in the
pad is called. In case of histograms, `TH1::Paint` invokes directly
`THistPainter::Paint`.

To draw a histogram `h` it is enough to do:

    h->Draw();

`h` can be of any kind: 1D, 2D or 3D. To choose how the histogram will
be drawn, the `Draw()` method can be invoked with an option. For instance
to draw a 2D histogram as a lego plot it is enough to do:

    h->Draw("lego");

`THistPainter` offers many options to paint 1D, 2D and 3D histograms.

When the `Draw()` method of a histogram is called for the first time
(`TH1::Draw`), it creates a `THistPainter` object and saves a
pointer to this "painter" as a data member of the histogram. The
`THistPainter` class specializes in the drawing of histograms. It is
separated from the histogram so that one can have histograms without the
graphics overhead, for example in a batch program. Each histogram having its own
painter (rather than a central singleton painter painting all histograms), allows
two histograms to be drawn in two threads without overwriting the painter's
values.

When a displayed histogram is filled again, there is no need to call the
`Draw()` method again; the image will be refreshed the next time the
pad will be updated.

A pad is updated after one of these three actions:

1. a carriage control on the ROOT command line,
2. a click inside the pad,
3. a call to `TPad::Update`.


By default a call to `TH1::Draw()` clears the pad of all objects
before drawing the new image of the histogram. One can use the `SAME`
option to leave the previous display intact and superimpose the new histogram.
The same histogram can be drawn with different graphics options in different
pads.

When a displayed histogram is deleted, its image is automatically removed
from the pad.

To create a copy of the histogram when drawing it, one can use
`TH1::DrawClone()`. This will clone the histogram and allow to change
and delete the original one without affecting the clone.


### <a name="HP01"></a> Histograms' plotting options


Most options can be concatenated with or without spaces or commas, for example:

    h->Draw("E1 SAME");

The options are not case sensitive:

    h->Draw("e1 same");


The default drawing option can be set with `TH1::SetOption` and retrieve
using `TH1::GetOption`:

    root [0] h->Draw();          // Draw "h" using the standard histogram representation.
    root [1] h->Draw("E");       // Draw "h" using error bars
    root [3] h->SetOption("E");  // Change the default drawing option for "h"
    root [4] h->Draw();          // Draw "h" using error bars
    root [5] h->GetOption();     // Retrieve the default drawing option for "h"
    (const Option_t* 0xa3ff948)"E"


#### <a name="HP01a"></a> Options supported for 1D and 2D histograms

| Option   | Description                                                       |
|----------|-------------------------------------------------------------------|
| "E"      | Draw error bars. |
| "AXIS"   | Draw only axis. |
| "AXIG"   | Draw only grid (if the grid is requested). |
| <a name="OPTHIST">"HIST"</a>   | When an histogram has errors it is visualized by default with error bars. To visualize it without errors use the option "HIST" together with the required option (eg "hist same c").  The "HIST" option can also be used to plot only the histogram and not the associated function(s). |
| "FUNC"   | When an histogram has a fitted function, this option allows to draw the fit result only. |
| "SAME"   | Superimpose on previous picture in the same pad. |
| "SAMES"  | Same as "SAME" and draw the statistics box|
| "PFC"    | Palette Fill Color: histogram's fill color is taken in the current palette. |
| "PLC"    | Palette Line Color: histogram's line color is taken in the current palette. |
| "PMC"    | Palette Marker Color: histogram's marker color is taken in the current palette. |
| "LEGO"   | Draw a lego plot with hidden line removal. |
| "LEGO1"  | Draw a lego plot with hidden surface removal. |
| "LEGO2"  | Draw a lego plot using colors to show the cell contents When the option "0" is used with any LEGO option, the empty bins are not drawn.|
| "LEGO3"  | Draw a lego plot with hidden surface removal, like LEGO1 but the border lines of each lego-bar are not drawn.|
| "LEGO4"  | Draw a lego plot with hidden surface removal, like LEGO1 but without the shadow effect on each lego-bar.|
| "TEXT"   | Draw bin contents as text (format set via `gStyle->SetPaintTextFormat`).|
| "TEXTnn" | Draw bin contents as text at angle nn (0 < nn < 90). |
| "X+"     | The X-axis is drawn on the top side of the plot. |
| "Y+"     | The Y-axis is drawn on the right side of the plot. |
| "MIN0"   | Set minimum value for the Y axis to 0, equivalent to gStyle->SetHistMinimumZero(). |

#### <a name="HP01b"></a> Options supported for 1D histograms

| Option   | Description                                                       |
|----------|-------------------------------------------------------------------|
| " "      | Default.                                                          |
| "AH"     | Draw histogram without axis. "A" can be combined with any drawing option. For instance, "AC" draws the histogram as a smooth Curve without axis.|
| "]["     | When this option is selected the first and last vertical lines of the histogram are not drawn.|
| "B"      | Bar chart option.|
| "BAR"    | Like option "B", but bars can be drawn with a 3D effect.|
| "HBAR"   | Like option "BAR", but bars are drawn horizontally.|
| "C"      | Draw a smooth Curve through the histogram bins.|
| "E0"     | Draw error bars. Markers are drawn for bins with 0 contents. Combined with E1 or E2 it avoids error bars clipping|
| "E1"     | Draw error bars with perpendicular lines at the edges.|
| "E2"     | Draw error bars with rectangles.|
| "E3"     | Draw a fill area through the end points of the vertical error bars.|
| "E4"     | Draw a smoothed filled area through the end points of the error bars.|
| "E5"     | Like E3 but ignore the bins with 0 contents.|
| "E6"     | Like E4 but ignore the bins with 0 contents.|
| "X0"     | When used with one of the "E" option, it suppress the error bar along X as `gStyle->SetErrorX(0)` would do.|
| "L"      | Draw a line through the bin contents.|
| "P"      | Draw current marker at each bin except empty bins.|
| "P0"     | Draw current marker at each bin including empty bins.|
| "PIE"    | Draw histogram as a Pie Chart.|
| "*H"     | Draw histogram with a * at each bin.|
| "LF2"    | Draw histogram like with option "L" but with a fill area. Note that "L" draws also a fill area if the hist fill color is set but the fill area corresponds to the histogram contour.|


#### <a name="HP01c"></a> Options supported for 2D histograms

| Option    | Description                                                      |
|-----------|------------------------------------------------------------------|
| " "       | Default (scatter plot).|
| "ARR"     | Arrow mode. Shows gradient between adjacent cells.|
| "BOX"     | A box is drawn for each cell with surface proportional to the content's absolute value. A negative content is marked with a X.|
| "BOX1"    | A button is drawn for each cell with surface proportional to content's absolute value. A sunken button is drawn for negative values a raised one for positive.|
| "COL"     | A box is drawn for each cell with a color scale varying with contents. All the none empty bins are painted. Empty bins are not painted unless some bins have a negative content because in that case the null bins might be not empty.  `TProfile2D` histograms are handled differently because, for this type of 2D histograms, it is possible to know if an empty bin has been filled or not. So even if all the bins' contents are positive some empty bins might be painted. And vice versa, if some bins have a negative content some empty bins might be not painted.|
| "COLZ"    | Same as "COL". In addition the color palette is also drawn.|
| "COL2"    | Alternative rendering algorithm to "COL". Can significantly improve rendering performance for large, non-sparse 2-D histograms.|
| "COLZ2"   | Same as "COL2". In addition the color palette is also drawn.|
| "Z CJUST" | In combination with colored options "COL","CONT0" etc: Justify labels in the color palette at color boudaries. For more details see `TPaletteAxis`|
| "CANDLE"  | Draw a candle plot along X axis.|
| "CANDLEX" | Same as "CANDLE".|
| "CANDLEY" | Draw a candle plot along Y axis.|
| "CANDLEXn"| Draw a candle plot along X axis. Different candle-styles with n from 1 to 6.|
| "CANDLEYn"| Draw a candle plot along Y axis. Different candle-styles with n from 1 to 6.|
| "VIOLIN"  | Draw a violin plot along X axis.|
| "VIOLINX" | Same as "VIOLIN".|
| "VIOLINY" | Draw a violin plot along Y axis.|
| "VIOLINXn"| Draw a violin plot along X axis. Different violin-styles with n being 1 or 2.|
| "VIOLINYn"| Draw a violin plot along Y axis. Different violin-styles with n being 1 or 2.|
| "CONT"    | Draw a contour plot (same as CONT0).|
| "CONT0"   | Draw a contour plot using surface colors to distinguish contours.|
| "CONT1"   | Draw a contour plot using line styles to distinguish contours.|
| "CONT2"   | Draw a contour plot using the same line style for all contours.|
| "CONT3"   | Draw a contour plot using fill area colors.|
| "CONT4"   | Draw a contour plot using surface colors (SURF option at theta = 0).|
| "CONT5"   | (TGraph2D only) Draw a contour plot using Delaunay triangles.|
| "LIST"    | Generate a list of TGraph objects for each contour.|
| "CYL"     | Use Cylindrical coordinates. The X coordinate is mapped on the angle and the Y coordinate on the cylinder length.|
| "POL"     | Use Polar coordinates. The X coordinate is mapped on the angle and the Y coordinate on the radius.|
| "SAME0"   | Same as "SAME" but do not use the z-axis range of the first plot. |
| "SAMES0"  | Same as "SAMES" but do not use the z-axis range of the first plot. |
| "SPH"     | Use Spherical coordinates. The X coordinate is mapped on the latitude and the Y coordinate on the longitude.|
| "PSR"     | Use PseudoRapidity/Phi coordinates. The X coordinate is mapped on Phi.|
| "SURF"    | Draw a surface plot with hidden line removal.|
| "SURF1"   | Draw a surface plot with hidden surface removal.|
| "SURF2"   | Draw a surface plot using colors to show the cell contents.|
| "SURF3"   | Same as SURF with in addition a contour view drawn on the top.|
| "SURF4"   | Draw a surface using Gouraud shading.|
| "SURF5"   | Same as SURF3 but only the colored contour is drawn. Used with option CYL, SPH or PSR it allows to draw colored contours on a sphere, a cylinder or a in pseudo rapidity space. In cartesian or polar coordinates, option SURF3 is used.|
| "LEGO9"   | Draw the 3D axis only. Mainly needed for internal use |
| "FB"      | With LEGO or SURFACE, suppress the Front-Box.|
| "BB"      | With LEGO or SURFACE, suppress the Back-Box.|
| "A"       | With LEGO or SURFACE, suppress the axis.|
| "SCAT"    | Draw a scatter-plot (default).|
| "[cutg]"  | Draw only the sub-range selected by the TCutG named "cutg".|


#### <a name="HP01d"></a> Options supported for 3D histograms

| Option   | Description                                                       |
|----------|-------------------------------------------------------------------|
| " "      | Default (scatter plot).|
| "ISO"    | Draw a Gouraud shaded 3d iso surface through a 3d histogram. It paints one surface at the value computed as follow: `SumOfWeights/(NbinsX*NbinsY*NbinsZ)`.|
| "BOX"    | Draw a for each cell with volume proportional to the content's absolute value. An hidden line removal algorithm is used|
| "BOX1"   | Same as BOX but an hidden surface removal algorithm is used|
| "BOX2"   | The boxes' colors are picked in the current palette according to the bins' contents|
| "BOX2Z"  | Same as "BOX2". In addition the color palette is also drawn.|
| "BOX3"   | Same as BOX1, but the border lines of each lego-bar are not drawn.|
| "LEGO"   | Same as `BOX`.|


#### <a name="HP01e"></a> Options supported for histograms' stacks (`THStack`)

| Option     | Description                                                     |
|------------|-----------------------------------------------------------------|
| " "        | Default, the histograms are drawn on top of each other (as lego plots for 2D histograms).|
| "NOSTACK"  | Histograms in the stack are all paint in the same pad as if the option `SAME` had been specified.|
| "NOSTACKB" | Histograms are  drawn next to each other as bar charts.|
| "PADS"     | The current pad/canvas is subdivided into a number of pads equal to the number of histograms in the stack and each histogram is paint into a separate pad.|
| "PFC"      | Palette Fill Color: stack's fill color is taken in the current palette. |
| "PLC"      | Palette Line Color: stack's line color is taken in the current palette. |
| "PMC"      | Palette Marker Color: stack's marker color is taken in the current palette. |



### <a name="HP02"></a> Setting the Style


Histograms use the current style (`gStyle`). When one changes the current
style and would like to propagate the changes to the histogram,
`TH1::UseCurrentStyle` should be called. Call `UseCurrentStyle` on
each histogram is needed.

To force all the histogram to use the current style use:

    gROOT->ForceStyle();

All the histograms read after this call will use the current style.


### <a name="HP03"></a> Setting line, fill, marker, and text attributes


The histogram classes inherit from the attribute classes:
`TAttLine`, `TAttFill` and `TAttMarker`.
See the description of these classes for the list of options.


### <a name="HP04"></a> Setting Tick marks on the histogram axis


The `TPad::SetTicks` method specifies the type of tick marks on the axis.
If ` tx = gPad->GetTickx()` and `ty = gPad->GetTicky()` then:

    tx = 1;   tick marks on top side are drawn (inside)
    tx = 2;   tick marks and labels on top side are drawn
    ty = 1;   tick marks on right side are drawn (inside)
    ty = 2;   tick marks and labels on right side are drawn

By default only the left Y axis and X bottom axis are drawn
(`tx = ty = 0`)

`TPad::SetTicks(tx,ty)` allows to set these options.
See also The `TAxis` functions to set specific axis attributes.

In case multiple color filled histograms are drawn on the same pad, the fill
area may hide the axis tick marks. One can force a redraw of the axis over all
the histograms by calling:

    gPad->RedrawAxis();


### <a name="HP05"></a> Giving titles to the X, Y and Z axis


    h->GetXaxis()->SetTitle("X axis title");
    h->GetYaxis()->SetTitle("Y axis title");

The histogram title and the axis titles can be any `TLatex` string.
The titles are part of the persistent histogram.


### <a name="HP060"></a> The option "SAME"


By default, when an histogram is drawn, the current pad is cleared before
drawing. In order to keep the previous drawing and draw on top of it the
option `SAME` should be use. The histogram drawn with the option
`SAME` uses the coordinates system available in the current pad.

This option can be used alone or combined with any valid drawing option but
some combinations must be use with care.

#### <a name="HP060a"></a> Limitations

- It does not work when combined with the `LEGO` and `SURF` options unless the
  histogram plotted with the option `SAME` has exactly the same
  ranges on the X, Y and Z axis as the currently drawn histogram. To superimpose
  lego plots [histograms' stacks](#HP26) should be used.


### <a name="HP061"></a> Colors automatically picked in palette

\since **ROOT version 6.09/01**

When several histograms are painted in the same canvas thanks to the option "SAME"
or via a `THStack` it might be useful to have an easy and automatic way to choose
their color. The simplest way is to pick colors in the current active color
palette. Palette coloring for histogram is activated thanks to the options `PFC`
(Palette Fill Color), `PLC` (Palette Line Color) and `PMC` (Palette Marker Color).
When one of these options is given to `TH1::Draw` the histogram get its color
from the current color palette defined by `gStyle->SetPalette(…)`. The color
is determined according to the number of objects having palette coloring in
the current pad.

Begin_Macro(source)
../../../tutorials/hist/histpalettecolor.C
End_Macro

Begin_Macro(source)
../../../tutorials/hist/thstackpalettecolor.C
End_Macro

Begin_Macro(source)
../../../tutorials/hist/thstack2palettecolor.C
End_Macro

### <a name="HP06"></a> Superimposing two histograms with different scales in the same pad


The following example creates two histograms, the second histogram is the bins
integral of the first one. It shows a procedure to draw the two histograms in
the same pad and it draws the scale of the second histogram using a new vertical
axis on the right side. See also the tutorial `transpad.C` for a variant
of this example.

Begin_Macro(source)
{
   auto c1 = new TCanvas("c1","c1",600,400);
   // create/fill draw h1
   gStyle->SetOptStat(kFALSE);
   auto h1 = new TH1F("h1","Superimposing two histograms with different scales",100,-3,3);
   Int_t i;
   for (i=0;i<10000;i++) h1->Fill(gRandom->Gaus(0,1));
   h1->Draw();
   c1->Update();

   // create hint1 filled with the bins integral of h1
   auto hint1 = new TH1F("hint1","h1 bins integral",100,-3,3);
   float sum = 0.f;
   for (i=1;i<=100;i++) {
      sum += h1->GetBinContent(i);
      hint1->SetBinContent(i,sum);
   }

   // scale hint1 to the pad coordinates
   float rightmax = 1.1*hint1->GetMaximum();
   float scale = gPad->GetUymax()/rightmax;
   hint1->SetLineColor(kRed);
   hint1->Scale(scale);
   hint1->Draw("same");

   // draw an axis on the right side
   auto axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
   gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
   axis->SetLineColor(kRed);
   axis->SetTextColor(kRed);
   axis->Draw();
}
End_Macro


### <a name="HP07"></a> Statistics Display


The type of information shown in the histogram statistics box can be selected
with:

    gStyle->SetOptStat(mode);

The `mode` has up to nine digits that can be set to on (1 or 2), off (0).

    mode = ksiourmen  (default = 000001111)
    k = 1;  kurtosis printed
    k = 2;  kurtosis and kurtosis error printed
    s = 1;  skewness printed
    s = 2;  skewness and skewness error printed
    i = 1;  integral of bins printed
    i = 2;  integral of bins with option "width" printed
    o = 1;  number of overflows printed
    u = 1;  number of underflows printed
    r = 1;  standard deviation printed
    r = 2;  standard deviation and standard deviation error printed
    m = 1;  mean value printed
    m = 2;  mean and mean error values printed
    e = 1;  number of entries printed
    n = 1;  name of histogram is printed

For example:

    gStyle->SetOptStat(11);

displays only the name of histogram and the number of entries, whereas:

    gStyle->SetOptStat(1101);

displays the name of histogram, mean value and standard deviation.

<b>WARNING 1:</b> never do:

    gStyle->SetOptStat(0001111);

but instead do:

    gStyle->SetOptStat(1111);

because `0001111` will be taken as an octal number!

<b>WARNING 2:</b> for backward compatibility with older versions

    gStyle->SetOptStat(1);

is taken as:

    gStyle->SetOptStat(1111)

To print only the name of the histogram do:

    gStyle->SetOptStat(1000000001);

<b>NOTE</b> that in case of 2D histograms, when selecting only underflow
(10000) or overflow (100000), the statistics box will show all combinations
of underflow/overflows and not just one single number.

The parameter mode can be any combination of the letters `kKsSiIourRmMen`

    k :  kurtosis printed
    K :  kurtosis and kurtosis error printed
    s :  skewness printed
    S :  skewness and skewness error printed
    i :  integral of bins printed
    I :  integral of bins with option "width" printed
    o :  number of overflows printed
    u :  number of underflows printed
    r :  standard deviation printed
    R :  standard deviation and standard deviation error printed
    m :  mean value printed
    M :  mean value mean error values printed
    e :  number of entries printed
    n :  name of histogram is printed

For example, to print only name of histogram and number of entries do:

    gStyle->SetOptStat("ne");

To print only the name of the histogram do:

    gStyle->SetOptStat("n");

The default value is:

    gStyle->SetOptStat("nemr");

When a histogram is painted, a `TPaveStats` object is created and added
to the list of functions of the histogram. If a `TPaveStats` object
already exists in the histogram list of functions, the existing object is just
updated with the current histogram parameters.

Once a histogram is painted, the statistics box can be accessed using
`h->FindObject("stats")`. In the command line it is enough to do:

    Root > h->Draw()
    Root > TPaveStats *st = (TPaveStats*)h->FindObject("stats")

because after `h->Draw()` the histogram is automatically painted. But
in a script file the painting should be forced using `gPad->Update()`
in order to make sure the statistics box is created:

    h->Draw();
    gPad->Update();
    TPaveStats *st = (TPaveStats*)h->FindObject("stats");

Without `gPad->Update()` the line `h->FindObject("stats")` returns a null pointer.

When a histogram is drawn with the option `SAME`, the statistics box
is not drawn. To force the statistics box drawing with the option
`SAME`, the option `SAMES` must be used.
If the new statistics box hides the previous statistics box, one can change
its position with these lines (`h` being the pointer to the histogram):

    Root > TPaveStats *st = (TPaveStats*)h->FindObject("stats")
    Root > st->SetX1NDC(newx1); //new x start position
    Root > st->SetX2NDC(newx2); //new x end position

To change the type of information for an histogram with an existing
`TPaveStats` one should do:

    st->SetOptStat(mode);

Where `mode` has the same meaning than when calling `gStyle->SetOptStat(mode)`
(see above).

One can delete the statistics box for a histogram `TH1* h` with:

    h->SetStats(0)

and activate it again with:

    h->SetStats(1).

Labels used in the statistics box ("Mean", "Std Dev", ...) can be changed from
`$ROOTSYS/etc/system.rootrc` or `.rootrc` (look for the string `Hist.Stats.`).


### <a name="HP08"></a> Fit Statistics


The type of information about fit parameters printed in the histogram statistics
box can be selected via the parameter mode. The parameter mode can be
`= pcev`  (default `= 0111`)

    p = 1;  print Probability
    c = 1;  print Chisquare/Number of degrees of freedom
    e = 1;  print errors (if e=1, v must be 1)
    v = 1;  print name/values of parameters

Example:

    gStyle->SetOptFit(1011);

print fit probability, parameter names/values and errors.

1. When `v = 1` is specified, only the non-fixed parameters are shown.
2. When `v = 2` all parameters are shown.

Note: `gStyle->SetOptFit(1)` means "default value", so it is equivalent
to `gStyle->SetOptFit(111)`


### <a name="HP09"></a> The error bars options


| Option   | Description                                                       |
|----------|-------------------------------------------------------------------|
| "E"      | Default. Shows only the error bars, not a marker.|
| "E1"     | Small lines are drawn at the end of the error bars.|
| "E2"     | Error rectangles are drawn.|
| "E3"     | A filled area is drawn through the end points of the vertical error bars.|
| "E4"     | A smoothed filled area is drawn through the end points of the vertical error bars.|
| "E0"     | Draw error bars. Markers are drawn for bins with 0 contents. Combined with E1 or E2 it avoids error bars clipping|
| "E5"     | Like E3 but ignore the bins with 0 contents.|
| "E6"     | Like E4 but ignore the bins with 0 contents.|
| "X0"     | When used with one of the "E" option, it suppress the error bar along X as `gStyle->SetErrorX(0)` would do.|

Begin_Macro(source)
{
   auto c1 = new TCanvas("c1","c1",600,400);
   auto he = new TH1F("he","Distribution drawn with error bars (option E1)  ",100,-3,3);
   for (int i=0; i<10000; i++) he->Fill(gRandom->Gaus(0,1));
   gStyle->SetEndErrorSize(3);
   gStyle->SetErrorX(1.);
   he->SetMarkerStyle(20);
   he->Draw("E1");
}
End_Macro

The options "E3" and "E4" draw an error band through the end points of the
vertical error bars. With "E4" the error band is smoothed. Because of the
smoothing algorithm used some artefacts may appear at the end of the band
like in the following example. In such cases "E3" should be used instead
of "E4".

Begin_Macro(source)
{
   auto ce4 = new TCanvas("ce4","ce4",600,400);
   ce4->Divide(2,1);
   auto he4 = new TH1F("he4","Distribution drawn with option E4",100,-3,3);
   Int_t i;
   for (i=0;i<10000;i++) he4->Fill(gRandom->Gaus(0,1));
   he4->SetFillColor(kRed);
   he4->GetXaxis()->SetRange(40,48);
   ce4->cd(1);
   he4->Draw("E4");
   ce4->cd(2);
   auto he3 = (TH1F*)he4->DrawClone("E3");
   he3->SetTitle("Distribution drawn option E3");
}
End_Macro

2D histograms can be drawn with error bars as shown is the following example:

Begin_Macro(source)
{
   auto c2e = new TCanvas("c2e","c2e",600,400);
   auto h2e = new TH2F("h2e","TH2 drawn with option E",40,-4,4,40,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      h2e->Fill(px,5*py);
   }
   h2e->Draw("E");
}
End_Macro


### <a name="HP100"></a> The bar chart option


The option "B" allows to draw simple vertical bar charts.
The bar width is controlled with `TH1::SetBarWidth()`,
and the bar offset within the bin, with `TH1::SetBarOffset()`.
These two settings are useful to draw several histograms on the
same plot as shown in the following example:

Begin_Macro(source)
{
   int i;
   const Int_t nx = 8;
   string os_X[nx]   = {"8","32","128","512","2048","8192","32768","131072"};
   float d_35_0[nx] = {0.75, -3.30, -0.92, 0.10, 0.08, -1.69, -1.29, -2.37};
   float d_35_1[nx] = {1.01, -3.02, -0.65, 0.37, 0.34, -1.42, -1.02, -2.10};

   auto cb = new TCanvas("cb","cb",600,400);
   cb->SetGrid();

   gStyle->SetHistMinimumZero();

   auto h1b = new TH1F("h1b","Option B example",nx,0,nx);
   h1b->SetFillColor(4);
   h1b->SetBarWidth(0.4);
   h1b->SetBarOffset(0.1);
   h1b->SetStats(0);
   h1b->SetMinimum(-5);
   h1b->SetMaximum(5);

   for (i=1; i<=nx; i++) {
      h1b->SetBinContent(i, d_35_0[i-1]);
      h1b->GetXaxis()->SetBinLabel(i,os_X[i-1].c_str());
   }

   h1b->Draw("b");

   auto h2b = new TH1F("h2b","h2b",nx,0,nx);
   h2b->SetFillColor(38);
   h2b->SetBarWidth(0.4);
   h2b->SetBarOffset(0.5);
   h2b->SetStats(0);
   for (i=1;i<=nx;i++) h2b->SetBinContent(i, d_35_1[i-1]);

   h2b->Draw("b same");
}
End_Macro


### <a name="HP10"></a> The "BAR" and "HBAR" options


When the option `bar` or `hbar` is specified, a bar chart is drawn. A vertical
bar-chart is drawn with the options `bar`, `bar0`, `bar1`, `bar2`, `bar3`, `bar4`.
An horizontal bar-chart is drawn with the options `hbar`, `hbar0`, `hbar1`,
`hbar2`, `hbar3`, `hbar4` (hbars.C).

- The bar is filled with the histogram fill color.
- The left side of the bar is drawn with a light fill color.
- The right side of the bar is drawn with a dark fill color.
- The percentage of the bar drawn with either the light or dark color is:
   - 0% for option "(h)bar" or "(h)bar0"
   - 10% for option "(h)bar1"
   - 20% for option "(h)bar2"
   - 30% for option "(h)bar3"
   - 40% for option "(h)bar4"

When an histogram has errors the option ["HIST"](#OPTHIST) together with the `(h)bar` option.

Begin_Macro(source)
../../../tutorials/hist/hbars.C
End_Macro

To control the bar width (default is the bin width) `TH1::SetBarWidth()`
should be used.

To control the bar offset (default is 0) `TH1::SetBarOffset()` should
be used.

These two parameters are useful when several histograms are plotted using
the option `SAME`. They allow to plot the histograms next to each other.


### <a name="HP11"></a> The SCATter plot option (default for 2D histograms)


For each cell (i,j) a number of points proportional to the cell content is
drawn. A maximum of `kNMAX` points per cell is drawn. If the maximum is above
`kNMAX` contents are normalized to `kNMAX` (`kNMAX=2000`).
If option is of the form `scat=ff`, (eg `scat=1.8`,
`scat=1e-3`), then `ff` is used as a scale factor to compute the
number of dots. `scat=1` is the default.

By default the scatter plot is painted with a "dot marker" which not scalable
(see the `TAttMarker` documentation). To change the marker size, a scalable marker
type should be used. For instance a circle (marker style 20).

Begin_Macro(source)
{
   auto c1 = new TCanvas("c1","c1",600,400);
   auto hscat = new TH2F("hscat","Option SCATter example (default for 2D histograms)  ",40,-4,4,40,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hscat->Fill(px,5*py);
      hscat->Fill(3+0.5*px,2*py-10.);
   }
   hscat->Draw("scat=0.5");
}
End_Macro


### <a name="HP12"></a> The ARRow option


Shows gradient between adjacent cells. For each cell (i,j) an arrow is drawn
The orientation of the arrow follows the cell gradient.

Begin_Macro(source)
{
   auto c1   = new TCanvas("c1","c1",600,400);
   auto harr = new TH2F("harr","Option ARRow example",20,-4,4,20,-20,20);
   harr->SetLineColor(kRed);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      harr->Fill(px,5*py);
      harr->Fill(3+0.5*px,2*py-10.,0.1);
   }
   harr->Draw("ARR");
}
End_Macro

\since **ROOT version 6.17/01**

The option `ARR` can be combined with the option `COL` or `COLZ`.

Begin_Macro(source)
{
   auto c1   = new TCanvas("c1","c1",600,400);
   auto harr = new TH2F("harr","Option ARR + COLZ example",20,-4,4,20,-20,20);
   harr->SetStats(0);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      harr->Fill(px,5*py);
      harr->Fill(3+0.5*px,2*py-10.,0.1);
   }
   harr->Draw("ARR COLZ");
}
End_Macro


### <a name="HP13"></a> The BOX option


For each cell (i,j) a box is drawn. The size (surface) of the box is
proportional to the absolute value of the cell content.
The cells with a negative content are drawn with a `X` on top of the box.

Begin_Macro(source)
{
   auto c1   = new TCanvas("c1","c1",600,400);
   auto hbox = new TH2F("hbox","Option BOX example",3,0,3,3,0,3);
   hbox->SetFillColor(42);
   hbox->Fill(0.5, 0.5,  1.);
   hbox->Fill(0.5, 1.5,  4.);
   hbox->Fill(0.5, 2.5,  3.);
   hbox->Fill(1.5, 0.5,  2.);
   hbox->Fill(1.5, 1.5, 12.);
   hbox->Fill(1.5, 2.5, -6.);
   hbox->Fill(2.5, 0.5, -4.);
   hbox->Fill(2.5, 1.5,  6.);
   hbox->Fill(2.5, 2.5,  0.5);
   hbox->Draw("BOX");
}
End_Macro

With option `BOX1` a button is drawn for each cell with surface
proportional to content's absolute value. A sunken button is drawn for
negative values a raised one for positive.

Begin_Macro(source)
{
   auto c1    = new TCanvas("c1","c1",600,400);
   auto hbox1 = new TH2F("hbox1","Option BOX1 example",3,0,3,3,0,3);
   hbox1->SetFillColor(42);
   hbox1->Fill(0.5, 0.5,  1.);
   hbox1->Fill(0.5, 1.5,  4.);
   hbox1->Fill(0.5, 2.5,  3.);
   hbox1->Fill(1.5, 0.5,  2.);
   hbox1->Fill(1.5, 1.5, 12.);
   hbox1->Fill(1.5, 2.5, -6.);
   hbox1->Fill(2.5, 0.5, -4.);
   hbox1->Fill(2.5, 1.5,  6.);
   hbox1->Fill(2.5, 2.5,  0.5);
   hbox1->Draw("BOX1");
}
End_Macro

When the option `SAME` (or "SAMES") is used with the option `BOX`,
the boxes' sizes are computing taking the previous plots into account. The range
along the Z axis is imposed by the first plot (the one without option
`SAME`); therefore the order in which the plots are done is relevant.

Begin_Macro(source)
{
   auto c1  = new TCanvas("c1","c1",600,400);
   auto hb1 = new TH2F("hb1","Example of BOX plots with option SAME ",40,-3,3,40,-3,3);
   auto hb2 = new TH2F("hb2","hb2",40,-3,3,40,-3,3);
   auto hb3 = new TH2F("hb3","hb3",40,-3,3,40,-3,3);
   auto hb4 = new TH2F("hb4","hb4",40,-3,3,40,-3,3);
   for (Int_t i=0;i<1000;i++) {
      double x,y;
      gRandom->Rannor(x,y);
      if (x>0 && y>0) hb1->Fill(x,y,4);
      if (x<0 && y<0) hb2->Fill(x,y,3);
      if (x>0 && y<0) hb3->Fill(x,y,2);
      if (x<0 && y>0) hb4->Fill(x,y,1);
   }
   hb1->SetFillColor(1);
   hb2->SetFillColor(2);
   hb3->SetFillColor(3);
   hb4->SetFillColor(4);
   hb1->Draw("box");
   hb2->Draw("box same");
   hb3->Draw("box same");
   hb4->Draw("box same");
}
End_Macro

\since **ROOT version 6.17/01:**

Sometimes the change of the range of the Z axis is unwanted, in which case, one
can use `SAME0` (or `SAMES0`) option to opt out of this change.

Begin_Macro(source)
{
   auto h2 = new TH2F("h2"," ",10,0,10,10,20,30);
   auto hf = (TH2F*)h2->Clone("hf");
   h2->SetBit(TH1::kNoStats);
   hf->SetBit(TH1::kNoStats);
   h2->Fill(5,22);
   h2->Fill(5,23);
   h2->Fill(6,22);
   h2->Fill(6,23);
   hf->Fill(6,23);
   hf->Fill(6,23);
   hf->Fill(6,23);
   hf->Fill(6,23);
   hf->Fill(5,23);

   auto hf_copy1 = hf->Clone("hf_copy1");
   auto lt = new TLatex();

   auto cx = new TCanvas(); cx->Divide(2,1);

   cx->cd(1);
   h2->Draw("box");
   hf->Draw("text colz same");
   lt->DrawLatexNDC(0.3,0.5,"SAME");

   cx->cd(2);
   h2->Draw("box");
   hf_copy1->Draw("text colz same0");
   lt->DrawLatexNDC(0.3,0.5,"SAME0");
}
End_Macro


### <a name="HP14"></a> The COLor option


For each cell (i,j) a box is drawn with a color proportional to the cell
content.

The color table used is defined in the current style.

If the histogram's minimum and maximum are the same (flat histogram), the
mapping on colors is not possible, therefore nothing is painted. To paint a
flat histogram it is enough to set the histogram minimum
(`TH1::SetMinimum()`) different from the bins' content.

The default number of color levels used to paint the cells is 20.
It can be changed with `TH1::SetContour()` or
`TStyle::SetNumberContours()`. The higher this number is, the smoother
is the color change between cells.

The color palette in TStyle can be modified via `gStyle->SetPalette()`.

All the non-empty bins are painted. Empty bins are not painted unless
some bins have a negative content because in that case the null bins
might be not empty.

`TProfile2D` histograms are handled differently because, for this type of 2D
histograms, it is possible to know if an empty bin has been filled or not. So even
if all the bins' contents are positive some empty bins might be painted. And vice versa,
if some bins have a negative content some empty bins might be not painted.

Combined with the option `COL`, the option `Z` allows to
display the color palette defined by `gStyle->SetPalette()`.

In the following example, the histogram has only positive bins; the empty
bins (containing 0) are not drawn.

Begin_Macro(source)
{
   auto c1    = new TCanvas("c1","c1",600,400);
   auto hcol1 = new TH2F("hcol1","Option COLor example ",40,-4,4,40,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hcol1->Fill(px,5*py);
   }
   hcol1->Draw("COLZ");
}
End_Macro

In the first plot of following example, the histogram has some negative bins;
the empty bins (containing 0) are drawn. In some cases one wants to not draw
empty bins (containing 0) of histograms having a negative minimum. The option
`1`, used to produce the second plot in the following picture, allows to do that.

Begin_Macro(source)
{
   auto c1 = new TCanvas("c1","c1",600,600);
   c1->Divide(1,2);
   auto hcol23 = new TH2F("hcol23","Option COLZ example ",40,-4,4,40,-20,20);
   auto hcol24 = new TH2F("hcol24","Option COLZ1 example ",40,-4,4,40,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hcol23->Fill(px,5*py);
      hcol24->Fill(px,5*py);
   }
   hcol23->Fill(0.,0.,-200.);
   hcol24->Fill(0.,0.,-200.);
   c1->cd(1); hcol23->Draw("COLZ");
   c1->cd(2); hcol24->Draw("COLZ1");
}
End_Macro

When the maximum of the histogram is set to a smaller value than the real maximum,
 the bins having a content between the new maximum and the real maximum are
painted with the color corresponding to the new maximum.

When the minimum of the histogram is set to a greater value than the real minimum,
 the bins having a value between the real minimum and the new minimum are not drawn
 unless the option `0` is set.

The following example illustrates the option `0` combined with the option `COL`.

Begin_Macro(source)
{
   auto c1 = new TCanvas("c1","c1",600,600);
   c1->Divide(1,2);
   auto hcol21 = new TH2F("hcol21","Option COLZ",40,-4,4,40,-20,20);
   auto hcol22 = new TH2F("hcol22","Option COLZ0",40,-4,4,40,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hcol21->Fill(px,5*py);
      hcol22->Fill(px,5*py);
   }
   hcol21->SetBit(TH1::kNoStats);
   hcol22->SetBit(TH1::kNoStats);
   c1->cd(1); hcol21->Draw("COLZ");
   c1->cd(2); hcol22->Draw("COLZ0");
   hcol22->SetMaximum(100);
   hcol22->SetMinimum(40);
}
End_Macro

\since **ROOT version 6.09/01:**

When the option SAME (or "SAMES") is used with the option COL, the boxes' color
are computing taking the previous plots into account. The range along the Z axis
is imposed by the first plot (the one without option SAME); therefore the order
in which the plots are done is relevant. Same as [in the `BOX` option](#HP13), one can use
`SAME0` (or `SAMES0`) to opt out of this imposition.

Begin_Macro(source)
{
   auto c = new TCanvas("c","Example of col plots with option SAME",200,10,700,500);
   auto h1 = new TH2F("h1","h1",40,-3,3,40,-3,3);
   auto h2 = new TH2F("h2","h2",40,-3,3,40,-3,3);
   auto h3 = new TH2F("h3","h3",40,-3,3,40,-3,3);
   auto h4 = new TH2F("h4","h4",40,-3,3,40,-3,3);
   h1->SetBit(TH1::kNoStats);
   for (Int_t i=0;i<5000;i++) {
      double x,y;
      gRandom->Rannor(x,y);
      if(x>0 && y>0) h1->Fill(x,y,4);
      if(x<0 && y<0) h2->Fill(x,y,3);
      if(x>0 && y<0) h3->Fill(x,y,2);
      if(x<0 && y>0) h4->Fill(x,y,1);
   }
   h1->Draw("colz");
   h2->Draw("col same");
   h3->Draw("col same");
   h4->Draw("col same");
}
End_Macro

The option `COL` can be combined with the option `POL`:

Begin_Macro(source)
{
   auto c1 = new TCanvas("c1","c1",600,400);
   auto hcol1 = new TH2F("hcol1","Option COLor combined with POL",40,-4,4,40,-4,4);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hcol1->Fill(px,py);
   }
   hcol1->Draw("COLZPOL");
}
End_Macro

\since **ROOT version 6.07/03:**

A second rendering technique is also available with the COL2 and COLZ2 options.

These options provide potential performance improvements compared to the standard
COL option. The performance comparison of the COL2 to the COL option depends on
the histogram and the size of the rendering region in the current pad. In general,
a small (approx. less than 100 bins per axis), sparsely populated TH2 will render
faster with the COL option.

However, for larger histograms (approx. more than 100 bins per axis)
that are not sparse, the COL2 option will provide up to 20 times performance improvements.
For example, a 1000x1000 bin TH2 that is not sparse will render an order of magnitude
faster with the COL2 option.

The COL2 option will also scale its performance based on the size of the
pixmap the histogram image is being rendered into. It also is much better optimized for
sessions where the user is forwarding X11 windows through an `ssh` connection.

For the most part, the COL2 and COLZ2 options are a drop in replacement to the COL
and COLZ options. There is one major difference and that concerns the treatment of
bins with zero content. The COL2 and COLZ2 options color these bins the color of zero.

COL2 option renders the histogram as a bitmap. Therefore it cannot be saved in vector
graphics file format like PostScript or PDF (an empty image will be generated). It can
be saved only in bitmap files like PNG format for instance.


### <a name="HP140"></a> The CANDLE and VIOLIN options

The mechanism behind Candle plots and Violin plots is very similar. Because of this they are
implemented in the same class TCandle. The keywords CANDLE or VIOLIN will initiate the drawing of
the corresponding plots. Followed by the keyword the user can select a plot direction (X or V for
vertical projections, or Y or H for horizontal projections) and/or predefined definitions
(1-6 for candles, 1-2 for violins). The order doesn't matter. Default is X and 1.

Instead of using the predefined representations, the candle and violin parameters can be
changed individually. In that case the option have the following form:

    CANDLEX(<option-string>)
    CANDLEY(<option-string>)
    VIOLINX(<option-string>)
    VIOLINY(<option-string>).

All zeros at the beginning of `option-string` can be omitted.

`option-string` consists eight values, defined as follow:

    "CANDLEX(zhpawMmb)"

Where:

  -  `b = 0`;  no box drawn
  -  `b = 1`;  the box is drawn. As the candle-plot is also called a box-plot it
               makes sense in the very most cases to always draw the box
  -  `b = 2`;  draw a filled box with border

  -  `m = 0`;  no median drawn
  -  `m = 1`;  median is drawn as a line
  -  `m = 2`;  median is drawn with errors (notches)
  -  `m = 3`;  median is drawn as a circle

  -  `M = 0`;  no mean drawn
  -  `M = 1`;  mean is drawn as a dashed line
  -  `M = 3`;  mean is drawn as a circle

  -  `w = 0`;  no whisker drawn
  -  `w = 1`;  whisker is drawn to end of distribution.
  -  `w = 2`;  whisker is drawn to max 1.5*iqr

  -  `a = 0`;  no anchor drawn
  -  `a = 1`;  the anchors are drawn

  -  `p = 0`;  no points drawn
  -  `p = 1`;  only outliers are drawn
  -  `p = 2`;  all datapoints are drawn
  -  `p = 3`:  all datapoints are drawn scattered

  -  `h = 0`;  no histogram is drawn
  -  `h = 1`;  histogram at the left or bottom side is drawn
  -  `h = 2`;  histogram at the right or top side is drawn
  -  `h = 3`;  histogram at left and right or top and bottom (violin-style) is drawn

  -  `z = 0`;  no zero indicator line is drawn
  -  `z = 1`;  zero indicator line is drawn.

As one can see all individual options for both candle and violin plots can be accessed by this
mechanism. In deed the keywords CANDLE(<option-string>) and VIOLIN(<option-string>) have the same
meaning. So you can parametrise an option-string for a candle plot and use the keywords VIOLIN and
vice versa, if you wish.

Using a logarithmic x- or y-axis is possible for candle and violin charts.

\since **ROOT version 6.11/01**

a logarithmic z-axis is possible, too but will only affect violin charts of course.

#### <a name="HP140a"></a> The CANDLE option

<a href="http://en.wikipedia.org/wiki/Box_plot">A Candle plot</a> (also known as
a "box plot" or "whisker plot") was invented in 1977 by John Tukey. It is a convenient
way to describe graphically a data distribution (D) with only five numbers:

  1. The minimum value of the distribution D (bottom or left whisker).
  2. The lower quartile (Q1): 25% of the data points in D are less than Q1 (bottom of the box).
  3. The median (M): 50% of the data points in D are less than M.
  4. The upper quartile (Q3): 75% of the data points in D are less than Q3 (top of the box).
  5. The maximum value of the distribution D (top or right whisker).

In this implementation a TH2 is considered as a collection of TH1 along
X (option `CANDLE` or `CANDLEX`) or Y (option `CANDLEY`).
Each TH1 is represented as one candle.

Begin_Macro(source)
../../../tutorials/hist/candleplotwhiskers.C
End_Macro

The candle reduces the information coming from a whole distribution into few values.
Independently from the number of entries or the significance of the underlying distribution
a candle will always look like a candle. So candle plots should be used carefully in
particular with unknown distributions. The definition of a candle is based on
__unbinned data__. Here, candles are created from binned data. Because of this, the
deviation is connected to the bin width used. The calculation of the quantiles
normally done on unbinned data also. Because data are binned, this will
only work the best possible way within the resolution of one bin

Because of all these facts one should take care that:

  - there are enough points per candle
  - the bin width is small enough (more bins will increase the maximum
    available resolution of the quantiles although there will be some
    bins with no entries)
  - never make a candle-plot if the underlying distribution is double-distributed
  - only create candles of distributions that are more-or-less gaussian (the
    MPV should be not too far away from the mean).

#### What a candle is made of

\since **ROOT version 6.07/05**

##### The box
The box displays the position of the inter-quantile-range of the underlying
distribution. The box contains 25% of the distribution below the median
and 25% of the distribution above the median. If the underlying distribution is large
enough and gaussian shaped the end-points of the box represent \f$ 0.6745\times\sigma \f$
(Where \f$ \sigma \f$ is the standard deviation of the gaussian). The width and
the position of the box can be modified by SetBarWidth() and SetBarOffset().
The +-25% quantiles are calculated by the GetQuantiles() methods.

\since **ROOT version 6.11/01**

Using the static function TCandle::SetBoxRange(double) the box definition will be
overwritten. E.g. using a box range of 0.68 will redefine the area of the lower box edge
to the upper box edge in order to cover 68% of the distribution illustrated by that candle.
The static function will affect all candle-charts in the running program.
Default is 0.5.

Using the static function TCandle::SetScaledCandle(bool) the width of the box (and the
whole candle) can be influenced. Deactivated, the width is constant (to be set by
SetBarWidth() ). Activated, the width of the boxes will be scaled to each other based on the
amount of data in the corresponding candle, the maximum width can be influenced by
SetBarWidth(). The static function will affect all candle-charts in the running program.
Default is false. Scaling between multiple candle-charts (using "same" or THStack) is not
supported, yet

##### The Median
For a sorted list of numbers, the median is the value in the middle of the list.
E.g. if a sorted list is made of five numbers "1,2,3,6,7" 3 will be the median
because it is in the middle of the list. If the number of entries is even the
average of the two values in the middle will be used. As histograms are binned
data, the situation is a bit more complex. The following example shows this:

~~~ {.cpp}
void quantiles() {
   auto h = new TH1I("h","h",10,0,10);
   //h->Fill(3);
   //h->Fill(3);
   h->Fill(4);
   h->Draw();
   double p = 0.;
   double q = 0.;
   h->GetQuantiles(1,&q,&p);

   cout << "Median is: " << q << std::endl;
}
~~~

Here the bin-width is 1.0. If the two Fill(3) are commented out, as there are currently,
the example will return a calculated median of 4.5, because that's the bin center
of the bin in which the value 4.0 has been dropped. If the two Fill(3) are not
commented out, it will return 3.75, because the algorithm tries to evenly distribute
the individual values of a bin with bin content > 0. It means the sorted list
would be "3.25, 3.75, 4.5".

The consequence is a median of 3.75. This shows how important it is to use a
small enough bin-width when using candle-plots on binned data.
If the distribution is large enough and gaussian shaped the median will be exactly
equal to the mean.
The median can be shown as a line or as a circle or not shown at all.

In order to show the significance of the median notched candle plots apply a "notch" or
narrowing of the box around the median. The significance is defined by
\f$ 1.57\times\frac{iqr}{N} \f$ and will be represented as the size of the notch
(where iqr is the size of the box and N is the number of entries of the whole
distribution). Candle plots like these are usually called "notched candle plots".

In case the significance of the median is greater that the size of the box, the
box will have an unnatural shape. Usually it means the chart has not enough data,
or that representing this uncertainty is not useful

##### The Mean
The mean can be drawn as a dashed line or as a circle or not drawn at all.
The mean is the arithmetic average of the values in the distribution.
It is calculated using GetMean(). Because histograms are
binned data, the mean value can differ from a calculation on the raw-data.
If the distribution is large enough and gaussian shaped the mean will be
exactly the median.

##### The Whiskers
The whiskers represent the part of the distribution not covered by the box.
The upper 25% and the lower 25% of the distribution are located within the whiskers.
Two representations are available.

  - A simple one (using w=1) defining the lower whisker from the lowest data value
    to the bottom of the box, and the upper whisker from the top of the box to the
    highest data value. In this representation the whisker-lines are dashed.
  - A more complex one having a further restriction. The whiskers are still connected
    to the box but their length cannot exceed \f$ 1.5\times iqr \f$. So it might
    be that the outermost part of the underlying distribution will not be covered
    by the whiskers. Usually these missing parts will be represented by the outliers
    (see points). Of course the upper and the lower whisker may differ in length.
    In this representation the whiskers are drawn as solid lines.

\since **ROOT version 6.11/01**

Using the static function TCandle::SetWhiskerRange(double) the whisker definition w=1
will be overwritten. E.g. using a whisker-range of 0.95 and w=1 will redefine the area of
the lower whisker to the upper whisker in order to cover 95% of the distribution inside
that candle. The static function will affect all candle-charts in the running program.
Default is 1.

If the distribution is large enough and gaussian shaped, the maximum length of
the whisker will be located at \f$ \pm 2.698 \sigma \f$ (when using the
1.5*iqr-definition (w=2), where \f$ \sigma \f$ is the standard deviation
(see picture above). In that case 99.3% of the total distribution will be covered
by the box and the whiskers, whereas 0.7% are represented by the outliers.

##### The Anchors
The anchors have no special meaning in terms of statistical calculation. They mark
the end of the whiskers and they have the width of the box. Both representation
with and without anchors are common.

##### The Points
Depending on the configuration the points can have different meanings:
  - If p=1 the points represent the outliers. If they are shown, it means
    some parts of the underlying distribution are not covered by the whiskers.
    This can only occur when the whiskers are set to option w=2. Here the whiskers
    can have a maximum length of \f$ 1.5 \times iqr \f$. So any points outside the
    whiskers will be drawn as outliers. The outliers will be represented by crosses.
  - If p=2 all points in the distribution will be painted as crosses. This is
    useful for small datasets only (up to 10 or 20 points per candle).
    The outliers are shown along the candle. Because the underlying distribution
    is binned, is frequently occurs that a bin contains more than one value.
    Because of this the points will be randomly scattered within their bin along
    the candle axis. If the bin content for a bin is exactly 1 (usually
    this happens for the outliers) if will be drawn in the middle of the bin along
    the candle axis. As the maximum number of points per candle is limited by kNMax/2
    on very large datasets scaling will be performed automatically. In that case one
    would loose all outliers because they have usually a bin content of 1 (and a
    bin content between 0 and 1 after the scaling). Because of this all bin contents
    between 0 and 1 - after the scaling - will be forced to be 1.
  - As the drawing of all values on large datasets can lead to big amounts of crosses,
    one can show all values as a scatter plot instead by choosing p=3. The points will be
    drawn as dots and will be scattered within the width of the candle. The color
    of the points will be the color of the candle-chart.

##### Other Options
Is is possible to combine all options of candle and violin plots with each other. E.g. a box-plot
with a histogram.

#### How to use the candle-plots drawing option

There are six predefined candle-plot representations:

  - "CANDLEX1": Standard candle (whiskers cover the whole distribution)
  - "CANDLEX2": Standard candle with better whisker definition + outliers.
                It is a good compromise
  - "CANDLEX3": Like candle2 but with a mean as a circle.
                It is easier to distinguish mean and median
  - "CANDLEX4": Like candle3 but showing the uncertainty of the median as well
                (notched candle plots).
                For bigger datasets per candle
  - "CANDLEX5": Like candle2 but showing all data points.
                For very small datasets
  - "CANDLEX6": Like candle2 but showing all datapoints scattered.
                For huge datasets


The following picture shows how the six predefined representations look.

Begin_Macro
{
   auto c1 = new TCanvas("c1","c1",700,800);
   c1->Divide(2,3);
   gStyle->SetOptStat(kFALSE);

   auto hcandle = new TH2F("hcandle"," ",10,-4,4,40,-20,20);
   float px, py;
   for (Int_t i = 0; i < 15000; i++) {
      gRandom->Rannor(px,py);
      hcandle->Fill(px,5*py);
   }
   hcandle->SetMarkerSize(0.5);

   TH2F *h2;
   for (Int_t i=1; i<7; i++) {
      c1->cd(i);
      h2 = (TH2F*)hcandle->DrawClone(Form("CANDLE%d",i));
      h2->SetTitle(Form("CANDLE%d",i));
   }
}
End_Macro


#### Example 1
Box and improved whisker, no mean, no median, no anchor no outliers

    h1->Draw("CANDLEX(2001)");

#### Example 2
A Candle-definition like "CANDLEX2" (New standard candle with better whisker definition + outliers)

    h1->Draw("CANDLEX(112111)");

#### Example 3
The following example shows how several candle plots can be super-imposed using
the option SAME. Note that the bar-width and bar-offset are active on candle plots.
Also the color, the line width, the size of the points and so on can be changed by the
standard attribute setting methods such as SetLineColor() SetLineWidth().

Begin_Macro(source)
../../../tutorials/hist/candleplot.C
End_Macro

#### <a name="HP140b"></a> The VIOLIN option

<a href="http://en.wikipedia.org/wiki/Violin_plot">A violin plot</a> is a candle plot
that also encodes the pdf information at each point.


Quartiles and mean are also represented at each point, with a marker
and two lines.

In this implementation a TH2 is considered as a collection of TH1 along
X (option `VIOLIN` or `VIOLINX`) or Y (option `VIOLINY`).

#### What a violin is made of

\since **ROOT version 6.09/02**

##### The histogram
The histogram is typically drawn to both directions with respect to the middle-line of the
corresponding bin. This can be achieved by using h=3. It is possible to draw a histogram only to
one side (h=1, or h=2).
The maximum number of bins in the histogram is limited to 500, if the number of bins in the used
histogram is higher it will be rebinned automatically. The maximum height of the histogram can
be modified by using SetBarWidth() and the position can be changed with SetBarOffset().
A solid fill style is recommended.

\since **ROOT version 6.11/01**

Using the static function TCandle::SetScaledViolin(bool) the height of the histogram or the
violin can be influenced. Activated, the height of the bins of the individual violins will be
scaled with respect to each other, the maximum height can be influenced by SetBarWidth().
Deactivated, the height of the bin with the maximum content of each individual violin is
set to a constant value using SetBarWidth(). The static function will affect all violin-charts
in the running program. Default is true. Scaling between multiple violin-charts
(using "same" or THStack) is not supported, yet.

##### The zero indicator line
Typical for violin charts is a line in the background over the whole histogram indicating
the bins with zero entries. The zero indicator line can be activated with z=1. The line color
will always be the same as the fill-color of the histogram.

##### The Mean
The Mean is illustrated with the same mechanism as used for candle plots. Usually a circle is used.

##### Whiskers
The whiskers are illustrated by the same mechanism as used for candle plots. There is only one
difference. When using the simple whisker definition (w=1) and the zero indicator line (z=1), then
the whiskers will be forced to be solid (usually hashed)

##### Points
The points are illustrated by the same mechanism as used for candle plots. E.g. VIOLIN2 uses
better whisker definition (w=2) and outliers (p=1).

##### Other options
It is possible to combine all options of candle or violin plots with each other. E.g. a violin plot
including a box-plot.

#### How to use the violin-plots drawing option

There are two predefined violin-plot representations:
  - "VIOLINX1": Standard violin (histogram, mean, whisker over full distribution,
                zero indicator line)
  - "VIOLINX2": Line VIOLINX1 both with better whisker definition + outliers.

A solid fill style is recommended for this plot (as opposed to a hollow or
hashed style).

Begin_Macro(source)
{
    auto c1 = new TCanvas("c1","c1",600,400);
    Int_t nx(6), ny(40);
    double xmin(0.0), xmax(+6.0), ymin(0.0), ymax(+4.0);
    auto hviolin = new TH2F("hviolin", "Option VIOLIN example", nx, xmin, xmax, ny, ymin, ymax);
    TF1 f1("f1", "gaus", +0,0 +4.0);
    double x,y;
    for (Int_t iBin=1; iBin<hviolin->GetNbinsX(); ++iBin) {
        double xc = hviolin->GetXaxis()->GetBinCenter(iBin);
        f1.SetParameters(1, 2.0+TMath::Sin(1.0+xc), 0.2+0.1*(xc-xmin)/xmax);
        for(Int_t i=0; i<10000; ++i){
            x = xc;
            y = f1.GetRandom();
            hviolin->Fill(x, y);
        }
    }
    hviolin->SetFillColor(kGray);
    hviolin->SetMarkerStyle(20);
    hviolin->SetMarkerSize(0.5);
    hviolin->Draw("VIOLIN");
    c1->Update();
}
End_Macro

The next example illustrates a time development of a certain value:

Begin_Macro(source)
../../../tutorials/hist/candledecay.C
End_Macro


### <a name="HP15"></a> The TEXT and TEXTnn Option


For each bin the content is printed. The text attributes are:

-  text font = current TStyle font (`gStyle->SetTextFont()`).
-  text size = 0.02*padheight*markersize (if `h` is the histogram drawn
   with the option `TEXT` the marker size can be changed with
   `h->SetMarkerSize(markersize)`).
-  text color = marker color.

By default the format `g` is used. This format can be redefined
by calling `gStyle->SetPaintTextFormat()`.

It is also possible to use `TEXTnn` in order to draw the text with
the angle `nn` (`0 < nn < 90`).

For 2D histograms the text is plotted in the center of each non empty cells.
It is possible to plot empty cells by calling `gStyle->SetHistMinimumZero()`
or providing MIN0 draw option. For 1D histogram the text is plotted at a y
position equal to the bin content.

For 2D histograms when the option "E" (errors) is combined with the option
text ("TEXTE"), the error for each bin is also printed.

Begin_Macro(source)
{
   auto c01 = new TCanvas("c01","c01",700,400);
   c01->Divide(2,1);
   auto htext1 = new TH1F("htext1","Option TEXT on 1D histograms ",10,-4,4);
   auto htext2 = new TH2F("htext2","Option TEXT on 2D histograms ",10,-4,4,10,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      htext1->Fill(px,0.1);
      htext2->Fill(px,5*py,0.1);
   }
   gStyle->SetPaintTextFormat("4.1f m");
   htext2->SetMarkerSize(1.8);
   c01->cd(1);
   htext2->Draw("TEXT45");
   c01->cd(2);
   htext1->Draw();
   htext1->Draw("HIST TEXT0 SAME");
}
End_Macro

\since **ROOT version 6.07/07:**

In case several histograms are drawn on top ot each other (using option `SAME`),
the text can be shifted using `SetBarOffset()`. It specifies an offset for the
text position in each cell, in percentage of the bin width.

Begin_Macro(source)
{
   auto c03 = new TCanvas("c03","c03",700,400);
   gStyle->SetOptStat(0);
   auto htext3 = new TH2F("htext3","Several 2D histograms drawn with option TEXT",10,-4,4,10,-20,20);
   auto htext4 = new TH2F("htext4","htext4",10,-4,4,10,-20,20);
   auto htext5 = new TH2F("htext5","htext5",10,-4,4,10,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      htext3->Fill(4*px,20*py,0.1);
      htext4->Fill(4*px,20*py,0.5);
      htext5->Fill(4*px,20*py,1.0);
   }
   htext4->SetMarkerSize(1.8);
   htext5->SetMarkerSize(1.8);
   htext5->SetMarkerColor(kRed);
   htext3->Draw("COL");
   htext4->SetBarOffset(0.2);
   htext4->Draw("TEXT SAME");
   htext5->SetBarOffset(-0.2);
   htext5->Draw("TEXT SAME");
}
End_Macro

In the case of profile histograms it is possible to print the number
of entries instead of the bin content. It is enough to combine the
option "E" (for entries) with the option "TEXT".

Begin_Macro(source)
{
   auto c02 = new TCanvas("c02","c02",700,400);
   c02->Divide(2,1);
   gStyle->SetPaintTextFormat("g");

   auto profile = new TProfile("profile","profile",10,0,10);
   profile->SetMarkerSize(2.2);
   profile->Fill(0.5,1);
   profile->Fill(1.5,2);
   profile->Fill(2.5,3);
   profile->Fill(3.5,4);
   profile->Fill(4.5,5);
   profile->Fill(5.5,5);
   profile->Fill(6.5,4);
   profile->Fill(7.5,3);
   profile->Fill(8.5,2);
   profile->Fill(9.5,1);
   c02->cd(1); profile->Draw("HIST TEXT0");
   c02->cd(2); profile->Draw("HIST TEXT0E");
}
End_Macro

### <a name="HP16"></a> The CONTour options


The following contour options are supported:

| Option   | Description                                                                 |
|----------|-----------------------------------------------------------------------------|
| "CONT"   | Draw a contour plot (same as CONT0).                                        |
| "CONT0"  | Draw a contour plot using surface colors to distinguish contours.           |
| "CONT1"  | Draw a contour plot using the line colors to distinguish contours.          |
| "CONT2"  | Draw a contour plot using the line styles (1 to 5) to distinguish contours. |
| "CONT3"  | Draw a contour plot using the same line style for all contours.             |
| "CONT4"  | Draw a contour plot using surface colors (`SURF` option at theta = 0).      |
| "CONT5"  | Draw a contour plot using Delaunay triangles.                               |



The following example shows a 2D histogram plotted with the option
`CONTZ`. The option `CONT` draws a contour plot using surface
colors to distinguish contours.  Combined with the option `CONT` (or
`CONT0`), the option `Z` allows to display the color palette
defined by `gStyle->SetPalette()`.

Begin_Macro(source)
{
   auto c1 = new TCanvas("c1","c1",600,400);
   auto hcontz = new TH2F("hcontz","Option CONTZ example ",40,-4,4,40,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hcontz->Fill(px-1,5*py);
      hcontz->Fill(2+0.5*px,2*py-10.,0.1);
   }
   hcontz->Draw("CONTZ");
}
End_Macro

The following example shows a 2D histogram plotted with the option
`CONT1Z`. The option `CONT1` draws a contour plot using the
line colors to distinguish contours. Combined with the option `CONT1`,
the option `Z` allows to display the color palette defined by
`gStyle->SetPalette()`.

Begin_Macro(source)
{
   auto c1 = new TCanvas("c1","c1",600,400);
   auto hcont1 = new TH2F("hcont1","Option CONT1Z example ",40,-4,4,40,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hcont1->Fill(px-1,5*py);
      hcont1->Fill(2+0.5*px,2*py-10.,0.1);
   }
   hcont1->Draw("CONT1Z");
}
End_Macro

The following example shows a 2D histogram plotted with the option
`CONT2`. The option `CONT2` draws a contour plot using the
line styles (1 to 5) to distinguish contours.

Begin_Macro(source)
{
   auto c1 = new TCanvas("c1","c1",600,400);
   auto hcont2 = new TH2F("hcont2","Option CONT2 example ",40,-4,4,40,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hcont2->Fill(px-1,5*py);
      hcont2->Fill(2+0.5*px,2*py-10.,0.1);
   }
   hcont2->Draw("CONT2");
}
End_Macro

The following example shows a 2D histogram plotted with the option
`CONT3`. The option `CONT3` draws contour plot using the same line style for
all contours.

Begin_Macro(source)
{
   auto c1 = new TCanvas("c1","c1",600,400);
   auto hcont3 = new TH2F("hcont3","Option CONT3 example ",40,-4,4,40,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hcont3->Fill(px-1,5*py);
      hcont3->Fill(2+0.5*px,2*py-10.,0.1);
   }
   hcont3->SetLineStyle(kDotted);
   hcont3->Draw("CONT3");
}
End_Macro

The following example shows a 2D histogram plotted with the option
`CONT4`. The option `CONT4` draws a contour plot using surface
colors to distinguish contours (`SURF` option at theta = 0). Combined
with the option `CONT` (or `CONT0`), the option `Z`
allows to display the color palette defined by `gStyle->SetPalette()`.

Begin_Macro(source)
{
   auto c1 = new TCanvas("c1","c1",600,400);
   auto hcont4 = new TH2F("hcont4","Option CONT4Z example ",40,-4,4,40,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hcont4->Fill(px-1,5*py);
      hcont4->Fill(2+0.5*px,2*py-10.,0.1);
   }
   hcont4->Draw("CONT4Z");
}
End_Macro

The default number of contour levels is 20 equidistant levels and can be changed
with `TH1::SetContour()` or `TStyle::SetNumberContours()`.

#### <a name="HP16a"></a> The LIST option

When option `LIST` is specified together with option
`CONT`, the points used to draw the contours are saved in
`TGraph` objects:

    h->Draw("CONT LIST");
    gPad->Update();

The contour are saved in `TGraph` objects once the pad is painted.
Therefore to use this functionality in a macro, `gPad->Update()`
should be performed after the histogram drawing. Once the list is
built, the contours are accessible in the following way:

    TObjArray *contours = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    Int_t ncontours     = contours->GetSize();
    TList *list         = (TList*)contours->At(i);

Where `i` is a contour number, and list contains a list of
`TGraph` objects.
For one given contour, more than one disjoint polyline may be generated.
The number of TGraphs per contour is given by:

    list->GetSize();

To access the first graph in the list one should do:

    TGraph *gr1 = (TGraph*)list->First();


The following example (ContourList.C) shows how to use this functionality.

Begin_Macro(source)
../../../tutorials/hist/ContourList.C
End_Macro

#### <a name="HP16b"></a> The AITOFF, MERCATOR, SINUSOIDAL and PARABOLIC options

The following options select the `CONT4` option and are useful for
sky maps or exposure maps (earth.C).

| Option       | Description                                                   |
|--------------|---------------------------------------------------------------|
| "AITOFF"     | Draw a contour via an AITOFF projection.|
| "MERCATOR"   | Draw a contour via an Mercator projection.|
| "SINUSOIDAL" | Draw a contour via an Sinusoidal projection.|
| "PARABOLIC"  | Draw a contour via an Parabolic projection.|

Begin_Macro(source)
../../../tutorials/graphics/earth.C
End_Macro


### <a name="HP17"></a> The LEGO options


In a lego plot the cell contents are drawn as 3-d boxes. The height of each box
is proportional to the cell content. The lego aspect is control with the
following options:

| Option   | Description                                                       |
|----------|-------------------------------------------------------------------|
| "LEGO"   | Draw a lego plot using the hidden lines removal technique.|
| "LEGO1"  | Draw a lego plot using the hidden surface removal technique.|
| "LEGO2"  | Draw a lego plot using colors to show the cell contents.|
| "LEGO3"  | Draw a lego plot with hidden surface removal, like LEGO1 but the border lines of each lego-bar are not drawn.|
| "LEGO4"  | Draw a lego plot with hidden surface removal, like LEGO1 but without the shadow effect on each lego-bar.|
| "0"      | When used with any LEGO option, the empty bins are not drawn.|


See the limitations with [the option "SAME"](#HP060a).

Line attributes can be used in lego plots to change the edges' style.

The following example shows a 2D histogram plotted with the option
`LEGO`. The option `LEGO` draws a lego plot using the hidden
lines removal technique.

Begin_Macro(source)
{
   auto c2 = new TCanvas("c2","c2",600,400);
   auto hlego = new TH2F("hlego","Option LEGO example ",40,-4,4,40,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hlego->Fill(px-1,5*py);
      hlego->Fill(2+0.5*px,2*py-10.,0.1);
   }
   hlego->Draw("LEGO");
}
End_Macro

The following example shows a 2D histogram plotted with the option
`LEGO1`. The option `LEGO1` draws a lego plot using the
hidden surface removal technique. Combined with any `LEGOn` option, the
option `0` allows to not drawn the empty bins.

Begin_Macro(source)
{
   auto c2 = new TCanvas("c2","c2",600,400);
   auto hlego1 = new TH2F("hlego1","Option LEGO1 example (with option 0)  ",40,-4,4,40,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hlego1->Fill(px-1,5*py);
      hlego1->Fill(2+0.5*px,2*py-10.,0.1);
   }
   hlego1->SetFillColor(kYellow);
   hlego1->Draw("LEGO1 0");
}
End_Macro

The following example shows a 2D histogram plotted with the option
`LEGO3`. Like the option `LEGO1`, the option `LEGO3`
draws a lego plot using the hidden surface removal technique but doesn't draw
the border lines of each individual lego-bar. This is very useful for histograms
having many bins. With such histograms the option `LEGO1` gives a black
image because of the border lines. This option also works with stacked legos.

Begin_Macro(source)
{
   auto c2 = new TCanvas("c2","c2",600,400);
   auto hlego3 = new TH2F("hlego3","Option LEGO3 example",40,-4,4,40,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hlego3->Fill(px-1,5*py);
      hlego3->Fill(2+0.5*px,2*py-10.,0.1);
   }
   hlego3->SetFillColor(kRed);
   hlego3->Draw("LEGO3");
}
End_Macro

The following example shows a 2D histogram plotted with the option
`LEGO2`. The option `LEGO2` draws a lego plot using colors to
show the cell contents.  Combined with the option `LEGO2`, the option
`Z` allows to display the color palette defined by
`gStyle->SetPalette()`.

Begin_Macro(source)
{
   auto c2 = new TCanvas("c2","c2",600,400);
   auto hlego2 = new TH2F("hlego2","Option LEGO2Z example ",40,-4,4,40,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hlego2->Fill(px-1,5*py);
      hlego2->Fill(2+0.5*px,2*py-10.,0.1);
   }
   hlego2->Draw("LEGO2Z");
}
End_Macro



### <a name="HP18"></a> The "SURFace" options


In a surface plot, cell contents are represented as a mesh.
The height of the mesh is proportional to the cell content.

| Option   | Description                                                       |
|----------|-------------------------------------------------------------------|
| "SURF"   | Draw a surface plot using the hidden line removal technique.|
| "SURF1"  | Draw a surface plot using the hidden surface removal technique.|
| "SURF2"  | Draw a surface plot using colors to show the cell contents.|
| "SURF3"  | Same as `SURF` with an additional filled contour plot on top.|
| "SURF4"  | Draw a surface using the Gouraud shading technique.|
| "SURF5"  | Used with one of the options CYL, PSR and CYL this option allows to draw a a filled contour plot.|
| "SURF6"  | This option should not be used directly. It is used internally when the CONT is used with option the option SAME on a 3D plot.|
| "SURF7"  | Same as `SURF2` with an additional line contour plot on top.|



See the limitations with [the option "SAME"](#HP060a).

The following example shows a 2D histogram plotted with the option
`SURF`. The option `SURF` draws a lego plot using the hidden
lines removal technique.

Begin_Macro(source)
{
   auto c2 = new TCanvas("c2","c2",600,400);
   auto hsurf = new TH2F("hsurf","Option SURF example ",30,-4,4,30,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hsurf->Fill(px-1,5*py);
      hsurf->Fill(2+0.5*px,2*py-10.,0.1);
   }
   hsurf->Draw("SURF");
}
End_Macro

The following example shows a 2D histogram plotted with the option
`SURF1`. The option `SURF1` draws a surface plot using the
hidden surface removal technique.  Combined with the option `SURF1`,
the option `Z` allows to display the color palette defined by
`gStyle->SetPalette()`.

Begin_Macro(source)
{
   auto c2 = new TCanvas("c2","c2",600,400);
   auto hsurf1 = new TH2F("hsurf1","Option SURF1 example ",30,-4,4,30,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hsurf1->Fill(px-1,5*py);
      hsurf1->Fill(2+0.5*px,2*py-10.,0.1);
   }
   hsurf1->Draw("SURF1");
}
End_Macro

The following example shows a 2D histogram plotted with the option
`SURF2`. The option `SURF2` draws a surface plot using colors
to show the cell contents. Combined with the option `SURF2`, the option
`Z` allows to display the color palette defined by
`gStyle->SetPalette()`.

Begin_Macro(source)
{
   auto c2 = new TCanvas("c2","c2",600,400);
   auto hsurf2 = new TH2F("hsurf2","Option SURF2 example ",30,-4,4,30,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hsurf2->Fill(px-1,5*py);
      hsurf2->Fill(2+0.5*px,2*py-10.,0.1);
   }
   hsurf2->Draw("SURF2");
}
End_Macro

The following example shows a 2D histogram plotted with the option
`SURF3`. The option `SURF3` draws a surface plot using the
hidden line removal technique with, in addition, a filled contour view drawn on the
top.  Combined with the option `SURF3`, the option `Z` allows
to display the color palette defined by `gStyle->SetPalette()`.

Begin_Macro(source)
{
   auto c2 = new TCanvas("c2","c2",600,400);
   auto hsurf3 = new TH2F("hsurf3","Option SURF3 example ",30,-4,4,30,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hsurf3->Fill(px-1,5*py);
      hsurf3->Fill(2+0.5*px,2*py-10.,0.1);
   }
   hsurf3->Draw("SURF3");
}
End_Macro

The following example shows a 2D histogram plotted with the option
`SURF4`. The option `SURF4` draws a surface using the Gouraud
shading technique.

Begin_Macro(source)
{
   auto c2 = new TCanvas("c2","c2",600,400);
   auto hsurf4 = new TH2F("hsurf4","Option SURF4 example ",30,-4,4,30,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hsurf4->Fill(px-1,5*py);
      hsurf4->Fill(2+0.5*px,2*py-10.,0.1);
   }
   hsurf4->SetFillColor(kOrange);
   hsurf4->Draw("SURF4");
}
End_Macro

The following example shows a 2D histogram plotted with the option
`SURF5 CYL`.  Combined with the option `SURF5`, the option
`Z` allows to display the color palette defined by `gStyle->SetPalette()`.

Begin_Macro(source)
{
   auto c2 = new TCanvas("c2","c2",600,400);
   auto hsurf5 = new TH2F("hsurf4","Option SURF5 example ",30,-4,4,30,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hsurf5->Fill(px-1,5*py);
      hsurf5->Fill(2+0.5*px,2*py-10.,0.1);
   }
   hsurf5->Draw("SURF5 CYL");
}
End_Macro

The following example shows a 2D histogram plotted with the option
`SURF7`. The option `SURF7` draws a surface plot using the
hidden surfaces removal technique with, in addition, a line contour view drawn on the
top.  Combined with the option `SURF7`, the option `Z` allows
to display the color palette defined by `gStyle->SetPalette()`.

Begin_Macro(source)
{
   auto c2 = new TCanvas("c2","c2",600,400);
   auto hsurf7 = new TH2F("hsurf3","Option SURF7 example ",30,-4,4,30,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hsurf7->Fill(px-1,5*py);
      hsurf7->Fill(2+0.5*px,2*py-10.,0.1);
   }
   hsurf7->Draw("SURF7");
}
End_Macro

As shown in the following example, when a contour plot is painted on top of a
surface plot using the option `SAME`, the contours appear in 3D on the
surface.

Begin_Macro(source)
{
   auto c20=new TCanvas("c20","c20",600,400);
   int NBins = 50;
   double d = 2;
   auto hsc = new TH2F("hsc", "Surface and contour with option SAME ", NBins, -d, d, NBins, -d, d);
   for (int bx = 1;  bx <= NBins; ++bx) {
      for (int by = 1;  by <= NBins; ++by) {
         double x = hsc->GetXaxis()->GetBinCenter(bx);
         double y = hsc->GetYaxis()->GetBinCenter(by);
         hsc->SetBinContent(bx, by, exp(-x*x)*exp(-y*y));
      }
   }
   hsc->Draw("surf2");
   hsc->Draw("CONT1 SAME");
}
End_Macro


### <a name="HP19"></a> Cylindrical, Polar, Spherical and PseudoRapidity/Phi options


Legos and surfaces plots are represented by default in Cartesian coordinates.
Combined with any `LEGOn` or `SURFn` options the following
options allow to draw a lego or a surface in other coordinates systems.

| Option   | Description                                                       |
|----------|-------------------------------------------------------------------|
| "CYL"    | Use Cylindrical coordinates. The X coordinate is mapped on the angle and the Y coordinate on the cylinder length.|
| "POL"    | Use Polar coordinates. The X coordinate is mapped on the angle and the Y coordinate on the radius.|
| "SPH"    | Use Spherical coordinates. The X coordinate is mapped on the latitude and the Y coordinate on the longitude.|
| "PSR"    | Use PseudoRapidity/Phi coordinates. The X coordinate is mapped on Phi.|



<b>WARNING:</b> Axis are not drawn with these options.

The following example shows the same histogram as a lego plot is the four
different coordinates systems.

Begin_Macro(source)
{
   auto c3 = new TCanvas("c3","c3",600,400);
   c3->Divide(2,2);
   auto hlcc = new TH2F("hlcc","Cylindrical coordinates",20,-4,4,20,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hlcc->Fill(px-1,5*py);
      hlcc->Fill(2+0.5*px,2*py-10.,0.1);
   }
   hlcc->SetFillColor(kYellow);
   c3->cd(1); hlcc->Draw("LEGO1 CYL");
   c3->cd(2); auto hlpc = (TH2F*) hlcc->DrawClone("LEGO1 POL");
   hlpc->SetTitle("Polar coordinates");
   c3->cd(3); auto hlsc = (TH2F*) hlcc->DrawClone("LEGO1 SPH");
   hlsc->SetTitle("Spherical coordinates");
   c3->cd(4); auto hlprpc = (TH2F*) hlcc->DrawClone("LEGO1 PSR");
   hlprpc->SetTitle("PseudoRapidity/Phi coordinates");
}
End_Macro

The following example shows the same histogram as a surface plot is the four different coordinates systems.

Begin_Macro(source)
{
   auto c4 = new TCanvas("c4","c4",600,400);
   c4->Divide(2,2);
   auto hscc = new TH2F("hscc","Cylindrical coordinates",20,-4,4,20,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hscc->Fill(px-1,5*py);
      hscc->Fill(2+0.5*px,2*py-10.,0.1);
   }
   c4->cd(1); hscc->Draw("SURF1 CYL");
   c4->cd(2); auto hspc = (TH2F*) hscc->DrawClone("SURF1 POL");
   hspc->SetTitle("Polar coordinates");
   c4->cd(3); auto hssc = (TH2F*) hscc->DrawClone("SURF1 SPH");
   hssc->SetTitle("Spherical coordinates");
   c4->cd(4); auto hsprpc = (TH2F*) hscc->DrawClone("SURF1 PSR");
   hsprpc->SetTitle("PseudoRapidity/Phi coordinates");
}
End_Macro


### <a name="HP20"></a> Base line for bar-charts and lego plots


By default the base line used to draw the boxes for bar-charts and lego plots is
the histogram minimum. It is possible to force this base line to be 0, using MIN0 draw
option or with the command:

    gStyle->SetHistMinimumZero();

Begin_Macro(source)
{
   auto c5 = new TCanvas("c5","c5",700,400);
   c5->Divide(2,1);
   auto hz1 = new TH1F("hz1","Bar-chart drawn from 0",20,-3,3);
   auto hz2 = new TH2F("hz2","Lego plot drawn from 0",20,-3,3,20,-3,3);
   Int_t i;
   double x,y;
   hz1->SetFillColor(kBlue);
   hz2->SetFillColor(kBlue);
   for (i=0;i<10000;i++) {
      x = gRandom->Gaus(0,1);
      y = gRandom->Gaus(0,1);
      if (x>0) {
         hz1->Fill(x,1);
         hz2->Fill(x,y,1);
      } else {
         hz1->Fill(x,-1);
         hz2->Fill(x,y,-2);
      }
   }
   c5->cd(1); hz1->Draw("bar2 min0");
   c5->cd(2); hz2->Draw("lego1 min0");
}
End_Macro

This option also works for horizontal plots. The example given in the section
["The bar chart option"](#HP100) appears as follow:

Begin_Macro(source)
{
   int i;
   const Int_t nx = 8;
   string os_X[nx]   = {"8","32","128","512","2048","8192","32768","131072"};
   float d_35_0[nx] = {0.75, -3.30, -0.92, 0.10, 0.08, -1.69, -1.29, -2.37};
   float d_35_1[nx] = {1.01, -3.02, -0.65, 0.37, 0.34, -1.42, -1.02, -2.10};

   auto cbh = new TCanvas("cbh","cbh",400,600);
   cbh->SetGrid();

   auto h1bh = new TH1F("h1bh","Option HBAR centered on 0",nx,0,nx);
   h1bh->SetFillColor(4);
   h1bh->SetBarWidth(0.4);
   h1bh->SetBarOffset(0.1);
   h1bh->SetStats(0);
   h1bh->SetMinimum(-5);
   h1bh->SetMaximum(5);

   for (i=1; i<=nx; i++) {
      h1bh->Fill(os_X[i-1].c_str(), d_35_0[i-1]);
      h1bh->GetXaxis()->SetBinLabel(i,os_X[i-1].c_str());
   }

   h1bh->Draw("hbar min0");

   auto h2bh = new TH1F("h2bh","h2bh",nx,0,nx);
   h2bh->SetFillColor(38);
   h2bh->SetBarWidth(0.4);
   h2bh->SetBarOffset(0.5);
   h2bh->SetStats(0);
   for (i=1;i<=nx;i++) h2bh->Fill(os_X[i-1].c_str(), d_35_1[i-1]);

   h2bh->Draw("hbar min0 same");
}
End_Macro


### <a name="HP20a"></a> TH2Poly Drawing


The following options are supported:

| Option   | Description                                                       |
|----------|-------------------------------------------------------------------|
| "SCAT"   | Draw a scatter plot (default).|
| "COL"    | Draw a color plot. All the bins are painted even the empty bins.|
| "COLZ"   | Same as "COL". In addition the color palette is also drawn.|
| "0"      | When used with any COL options, the empty bins are not drawn.|
| "TEXT"   | Draw bin contents as text (format set via `gStyle->SetPaintTextFormat`).|
| "TEXTN"  | Draw bin names as text.|
| "TEXTnn" | Draw bin contents as text at angle nn (0 < nn < 90).|
| "L"      | Draw the bins boundaries as lines. The lines attributes are the TGraphs ones.|
| "P"      | Draw the bins boundaries as markers. The markers attributes are the TGraphs ones.|
| "F"      | Draw the bins boundaries as filled polygons.  The filled polygons attributes are the TGraphs ones.|



`TH2Poly` can be drawn as a color plot (option COL). `TH2Poly` bins can have any
shapes. The bins are defined as graphs. The following macro is a very simple
example showing how to book a TH2Poly and draw it.

Begin_Macro(source)
{
   auto ch2p1 = new TCanvas("ch2p1","ch2p1",600,400);
   auto h2p = new TH2Poly();
   h2p->SetName("h2poly_name");
   h2p->SetTitle("h2poly_title");
   double px1[] = {0, 5, 6};
   double py1[] = {0, 0, 5};
   double px2[] = {0, -1, -1, 0};
   double py2[] = {0, 0, -1, 3};
   double px3[] = {4, 3, 0, 1, 2.4};
   double py3[] = {4, 3.7, 1, 3.7, 2.5};
   h2p->AddBin(3, px1, py1);
   h2p->AddBin(4, px2, py2);
   h2p->AddBin(5, px3, py3);
   h2p->Fill(0.1, 0.01, 3);
   h2p->Fill(-0.5, -0.5, 7);
   h2p->Fill(-0.7, -0.5, 1);
   h2p->Fill(1, 3, 1.5);
   double fx[] = {0.1, -0.5, -0.7, 1};
   double fy[] = {0.01, -0.5, -0.5, 3};
   double fw[] = {3, 1, 1, 1.5};
   h2p->FillN(4, fx, fy, fw);
   h2p->Draw("col");
}
End_Macro

Rectangular bins are a frequent case. The special version of
the `AddBin` method allows to define them more easily like
shown in the following example (th2polyBoxes.C).

Begin_Macro(source)
../../../tutorials/hist/th2polyBoxes.C
End_Macro

One `TH2Poly` bin can be a list of polygons. Such bins are defined
by calling `AddBin` with a `TMultiGraph`. The following example
shows a such case:

Begin_Macro(source)
{
   auto ch2p2 = new TCanvas("ch2p2","ch2p2",600,400);

   Int_t i, bin;
   const Int_t nx = 48;
   const char *states [nx] = {
      "alabama",      "arizona",        "arkansas",       "california",
      "colorado",     "connecticut",    "delaware",       "florida",
      "georgia",      "idaho",          "illinois",       "indiana",
      "iowa",         "kansas",         "kentucky",       "louisiana",
      "maine",        "maryland",       "massachusetts",  "michigan",
      "minnesota",    "mississippi",    "missouri",       "montana",
      "nebraska",     "nevada",         "new_hampshire",  "new_jersey",
      "new_mexico",   "new_york",       "north_carolina", "north_dakota",
      "ohio",         "oklahoma",       "oregon",         "pennsylvania",
      "rhode_island", "south_carolina", "south_dakota",   "tennessee",
      "texas",        "utah",           "vermont",        "virginia",
      "washington",   "west_virginia",  "wisconsin",      "wyoming"
   };
   Double_t pop[nx] = {
    4708708, 6595778,  2889450, 36961664, 5024748,  3518288,  885122, 18537969,
    9829211, 1545801, 12910409,  6423113, 3007856,  2818747, 4314113,  4492076,
    1318301, 5699478,  6593587,  9969727, 5266214,  2951996, 5987580,   974989,
    1796619, 2643085,  1324575,  8707739, 2009671, 19541453, 9380884,   646844,
   11542645, 3687050,  3825657, 12604767, 1053209,  4561242,  812383,  6296254,
   24782302, 2784572,   621760,  7882590, 6664195,  1819777, 5654774,   544270
   };

   Double_t lon1 = -130;
   Double_t lon2 = -65;
   Double_t lat1 = 24;
   Double_t lat2 = 50;
   auto p = new TH2Poly("USA","USA Population",lon1,lon2,lat1,lat2);

   TFile::SetCacheFileDir(".");
   auto f = TFile::Open("http://root.cern.ch/files/usa.root", "CACHEREAD");

   TMultiGraph *mg;
   TKey *key;
   TIter nextkey(gDirectory->GetListOfKeys());
   while ((key = (TKey*)nextkey())) {
      TObject *obj = key->ReadObj();
      if (obj->InheritsFrom("TMultiGraph")) {
         mg = (TMultiGraph*)obj;
         bin = p->AddBin(mg);
      }
   }

   for (i=0; i<nx; i++) p->Fill(states[i], pop[i]);

   gStyle->SetOptStat(11);
   p->Draw("COLZ L");
}
End_Macro

`TH2Poly` histograms can also be plotted using the GL interface using
the option "GLLEGO".

\since **ROOT version 6.09/01**

In some cases it can be useful to not draw the empty bins. the option "0"
combined with the option "COL" et COLZ allows to do that.

Begin_Macro(source)
{
   auto chc = new TCanvas("chc","chc",600,400);

   auto hc = new TH2Poly();
   hc->Honeycomb(0,0,.1,25,25);
   hc->SetName("hc");
   hc->SetTitle("Option COLZ 0");
   TRandom ran;
   for (int i = 0; i<300; i++) hc->Fill(ran.Gaus(2.,1), ran.Gaus(2.,1));
   hc->Draw("colz 0");
}
End_Macro

### <a name="HP21"></a> The SPEC option


This option allows to use the `TSpectrum2Painter` tools. See the full
documentation in `TSpectrum2Painter::PaintSpectrum`.


### <a name="HP22"></a> Option "Z" : Adding the color palette on the right side of the pad


When this option is specified, a color palette with an axis indicating the value
of the corresponding color is drawn on the right side of the picture. In case,
not enough space is left, one can increase the size of the right margin by
calling `TPad::SetRightMargin()`. The attributes used to display the
palette axis values are taken from the Z axis of the object. For example, to
set the labels size on the palette axis do:

    hist->GetZaxis()->SetLabelSize().

<b>WARNING:</b> The palette axis is always drawn vertically.


### <a name="HP23"></a> Setting the color palette


To change the color palette `TStyle::SetPalette` should be used, eg:

    gStyle->SetPalette(ncolors,colors);

For example the option `COL` draws a 2D histogram with cells
represented by a box filled with a color index which is a function
of the cell content.
If the cell content is N, the color index used will be the color number
in `colors[N]`, etc. If the maximum cell content is greater than
`ncolors`, all cell contents are scaled to `ncolors`.

If ` ncolors <= 0`, a default palette (see below) of 50 colors is
defined. This palette is recommended for pads, labels ...

`if ncolors == 1 && colors == 0`, then a Pretty Palette with a
Spectrum Violet->Red is created with 50 colors. That's the default rain bow
palette.

Other pre-defined palettes with 255 colors are available when `colors == 0`.
The following value of `ncolors` give access to:


    if ncolors = 51 and colors=0, a Deep Sea palette is used.
    if ncolors = 52 and colors=0, a Grey Scale palette is used.
    if ncolors = 53 and colors=0, a Dark Body Radiator palette is used.
    if ncolors = 54 and colors=0, a two-color hue palette palette is used.(dark blue through neutral gray to bright yellow)
    if ncolors = 55 and colors=0, a Rain Bow palette is used.
    if ncolors = 56 and colors=0, an inverted Dark Body Radiator palette is used.


If `ncolors > 0 && colors == 0`, the default palette is used with a maximum of ncolors.

The default palette defines:

-  index  0  to  9 : shades of grey
-  index 10  to 19 : shades of brown
-  index 20  to 29 : shades of blue
-  index 30  to 39 : shades of red
-  index 40  to 49 : basic colors

The color numbers specified in the palette can be viewed by selecting
the item `colors` in the `VIEW` menu of the canvas tool bar.
The red, green, and blue components of a color can be changed thanks to
`TColor::SetRGB()`.

\since **ROOT version 6.19/01**

As default labels and ticks are drawn by `TGAxis` at equidistant (lin or log)
points as controlled by SetNdivisions.
If option "CJUST" is given labels and ticks are justified at the
color boundaries defined by the contour levels.
For more details see `TPaletteAxis`

### <a name="HP24"></a> Drawing a sub-range of a 2D histogram; the [cutg] option


Using a `TCutG` object, it is possible to draw a sub-range of a 2D
histogram. One must create a graphical cut (mouse or C++) and specify the name
of the cut between `[]` in the `Draw()` option.
For example (fit2a.C), with a `TCutG` named `cutg`, one can call:

    myhist->Draw("surf1 [cutg]");

To invert the cut, it is enough to put a `-` in front of its name:

    myhist->Draw("surf1 [-cutg]");

It is possible to apply several cuts (`,` means logical AND):

    myhist->Draw("surf1 [cutg1,cutg2]");

Begin_Macro(source)
../../../tutorials/fit/fit2a.C
End_Macro

### <a name="HP25"></a> Drawing options for 3D histograms


| Option   | Description                                                       |
|----------|-------------------------------------------------------------------|
| "ISO"    | Draw a Gouraud shaded 3d iso surface through a 3d histogram. It paints one surface at the value computed as follow: `SumOfWeights/(NbinsX*NbinsY*NbinsZ)`|
| "BOX"    | Draw a for each cell with volume proportional to the content's absolute value. An hidden line removal algorithm is used|
| "BOX1"   | Same as BOX but an hidden surface removal algorithm is used|
| "BOX2"   | The boxes' colors are picked in the current palette according to the bins' contents|
| "BOX2Z"  | Same as "BOX2". In addition the color palette is also drawn.|
| "BOX3"   | Same as BOX1, but the border lines of each lego-bar are not drawn.|

Note that instead of `BOX` one can also use `LEGO`.

By default, like 2D histograms, 3D histograms are drawn as scatter plots.

The following example shows a 3D histogram plotted as a scatter plot.

Begin_Macro(source)
{
   auto c06 = new TCanvas("c06","c06",600,400);
   gStyle->SetOptStat(kFALSE);
   auto h3scat = new TH3F("h3scat","Option SCAT (default) ",15,-2,2,15,-2,2,15,0,4);
   double x, y, z;
   for (Int_t i=0;i<10000;i++) {
      gRandom->Rannor(x, y);
      z = x*x + y*y;
      h3scat->Fill(x,y,z);
   }
   h3scat->Draw();
}
End_Macro

The following example shows a 3D histogram plotted with the option `BOX`.

Begin_Macro(source)
{
   auto c16 = new TCanvas("c16","c16",600,400);
   gStyle->SetOptStat(kFALSE);
   auto h3box = new TH3F("h3box","Option BOX",15,-2,2,15,-2,2,15,0,4);
   double x, y, z;
   for (Int_t i=0;i<10000;i++) {
      gRandom->Rannor(x, y);
      z = x*x + y*y;
      h3box->Fill(x,y,z);
   }
   h3box->Draw("BOX");
}
End_Macro

The following example shows a 3D histogram plotted with the option `BOX1`.

Begin_Macro(source)
{
   auto c36 = new TCanvas("c36","c36",600,400);
   gStyle->SetOptStat(kFALSE);
   auto h3box = new TH3F("h3box","Option BOX1",10,-2.,2.,10,-2.,2.,10,-0.5,2.);
   double x, y, z;
   for (Int_t i=0;i<10000;i++) {
      gRandom->Rannor(x, y);
      z = abs(sin(x)/x + cos(y)*y);
      h3box->Fill(x,y,z);
   }
   h3box->SetFillColor(9);
   h3box->Draw("BOX1");
}
End_Macro

The following example shows a 3D histogram plotted with the option `BOX2`.

Begin_Macro(source)
{
   auto c56 = new TCanvas("c56","c56",600,400);
   gStyle->SetOptStat(kFALSE);
   auto h3box = new TH3F("h3box","Option BOX2",10,-2.,2.,10,-2.,2.,10,-0.5,2.);
   double x, y, z;
   for (Int_t i=0;i<10000;i++) {
      gRandom->Rannor(x, y);
      z = abs(sin(x)/x + cos(y)*y);
      h3box->Fill(x,y,z);
   }
   h3box->Draw("BOX2 Z");
}
End_Macro

The following example shows a 3D histogram plotted with the option `BOX3`.

Begin_Macro(source)
{
   auto c46 = new TCanvas("c46","c46",600,400);
   c46->SetFillColor(38);
   gStyle->SetOptStat(kFALSE);
   auto h3box = new TH3F("h3box","Option BOX3",15,-2,2,15,-2,2,15,0,4);
   double x, y, z;
   for (Int_t i=0;i<10000;i++) {
      gRandom->Rannor(x, y);
      z = x*x + y*y;
      h3box->Fill(x,y,z);
   }
   h3box->Draw("BOX3");
}
End_Macro

For all the `BOX` options each bin is drawn as a 3D box with a volume proportional
to the absolute value of the bin content. The bins with a negative content are
drawn with a X on each face of the box as shown in the following example:

Begin_Macro(source)
{
   auto c = new TCanvas("c","c",600,400);
   gStyle->SetOptStat(kFALSE);
   auto h3box = new TH3F("h3box","Option BOX1 with negative bins",3, 0., 4., 3, 0.,4., 3, 0., 4.);
   h3box->Fill(0., 2., 2.,  10.);
   h3box->Fill(2., 2., 2.,   5.);
   h3box->Fill(2., 2., .5,   2.);
   h3box->Fill(2., 2., 3.,   -1.);
   h3box->Fill(3., 2., 2., -10.);
   h3box->SetFillColor(8);
   h3box->Draw("box1");
}
End_Macro

The following example shows a 3D histogram plotted with the option `ISO`.

Begin_Macro(source)
{
   auto c26 = new TCanvas("c26","c26",600,400);
   gStyle->SetOptStat(kFALSE);
   auto h3iso = new TH3F("h3iso","Option ISO",15,-2,2,15,-2,2,15,0,4);
   double x, y, z;
   for (Int_t i=0;i<10000;i++) {
      gRandom->Rannor(x, y);
      z = x*x + y*y;
      h3iso->Fill(x,y,z);
   }
   h3iso->SetFillColor(kCyan);
   h3iso->Draw("ISO");
}
End_Macro


### <a name="HP26"></a> Drawing option for histograms' stacks


Stacks of histograms are managed with the `THStack`. A `THStack`
is a collection of `TH1` (or derived) objects. For painting only the
`THStack` containing `TH1` only or
`THStack` containing `TH2` only will be considered.

By default, histograms are shown stacked:

1. The first histogram is paint.
2. The the sum of the first and second, etc...

If the option `NOSTACK` is specified, the histograms are all paint in
the same pad as if the option `SAME` had been specified. This allows to
compute X and Y scales common to all the histograms, like
`TMultiGraph` does for graphs.

If the option `PADS` is specified, the current pad/canvas is
subdivided into a number of pads equal to the number of histograms and each
histogram is paint into a separate pad.

The following example shows various types of stacks (hstack.C).

Begin_Macro(source)
../../../tutorials/hist/hstack.C
End_Macro

The option `nostackb` allows to draw the histograms next to each
other as bar charts:

Begin_Macro(source)
{
   auto cst0 = new TCanvas("cst0","cst0",600,400);
   auto hs = new THStack("hs","Stacked 1D histograms: option #font[82]{\"nostackb\"}");

   auto h1 = new TH1F("h1","h1",10,-4,4);
   h1->FillRandom("gaus",20000);
   h1->SetFillColor(kRed);
   hs->Add(h1);

   auto h2 = new TH1F("h2","h2",10,-4,4);
   h2->FillRandom("gaus",15000);
   h2->SetFillColor(kBlue);
   hs->Add(h2);

   auto h3 = new TH1F("h3","h3",10,-4,4);
   h3->FillRandom("gaus",10000);
   h3->SetFillColor(kGreen);
   hs->Add(h3);

   hs->Draw("nostackb");
   hs->GetXaxis()->SetNdivisions(-10);
   cst0->SetGridx();
}
End_Macro

If at least one of the histograms in the stack has errors, the whole stack is
visualized by default with error bars. To visualize it without errors the
option `HIST` should be used.

Begin_Macro(source)
{
   auto cst1 = new TCanvas("cst1","cst1",700,400);
   cst1->Divide(2,1);

   auto hst11 = new TH1F("hst11", "", 20, -10, 10);
   hst11->Sumw2();
   hst11->FillRandom("gaus", 1000);
   hst11->SetFillColor(kViolet);
   hst11->SetLineColor(kViolet);

   auto hst12 = new TH1F("hst12", "", 20, -10, 10);
   hst12->FillRandom("gaus", 500);
   hst12->SetFillColor(kBlue);
   hst12->SetLineColor(kBlue);

   THStack st1("st1", "st1");
   st1.Add(hst11);
   st1.Add(hst12);

   cst1->cd(1); st1.Draw();
   cst1->cd(2); st1.Draw("hist");
}
End_Macro

### <a name="HP27"></a> Drawing of 3D implicit functions


3D implicit functions (`TF3`) can be drawn as iso-surfaces.
The implicit function f(x,y,z) = 0 is drawn in cartesian coordinates.
In the following example the options "FB" and "BB" suppress the
"Front Box" and "Back Box" around the plot.

Begin_Macro(source)
{
   auto c2 = new TCanvas("c2","c2",600,400);
   auto f3 = new TF3("f3","sin(x*x+y*y+z*z-36)",-2,2,-2,2,-2,2);
   f3->SetClippingBoxOn(0,0,0);
   f3->SetFillColor(30);
   f3->SetLineColor(15);
   f3->Draw("FBBB");
}
End_Macro


### <a name="HP28"></a> Associated functions drawing


An associated function is created by `TH1::Fit`. More than on fitted
function can be associated with one histogram (see `TH1::Fit`).

A `TF1` object `f1` can be added to the list of associated
functions of an histogram `h` without calling `TH1::Fit`
simply doing:

    h->GetListOfFunctions()->Add(f1);

or

    h->GetListOfFunctions()->Add(f1,someoption);

To retrieve a function by name from this list, do:

    TF1 *f1 = (TF1*)h->GetListOfFunctions()->FindObject(name);

or

    TF1 *f1 = h->GetFunction(name);

Associated functions are automatically painted when an histogram is drawn.
To avoid the painting of the associated functions the option `HIST`
should be added to the list of the options used to paint the histogram.


### <a name="HP29"></a> Drawing using OpenGL


The class `TGLHistPainter` allows to paint data set using the OpenGL 3D
graphics library. The plotting options start with `GL` keyword.
In addition, in order to inform canvases that OpenGL should be used to render
3D representations, the following option should be set:

    gStyle->SetCanvasPreferGL(true);


#### <a name="HP29a"></a> General information: plot types and supported options

The following types of plots are provided:

For lego plots the supported options are:

| Option   | Description                                                       |
|----------|-------------------------------------------------------------------|
| "GLLEGO" | Draw a lego plot. It works also for `TH2Poly`.|
| "GLLEGO2"| Bins with color levels.|
| "GLLEGO3"| Cylindrical bars.|



Lego painter in cartesian supports logarithmic scales for X, Y, Z.
In polar only Z axis can be logarithmic, in cylindrical only Y.

For surface plots (`TF2` and `TH2`) the supported options are:

| Option    | Description                                                      |
|-----------|------------------------------------------------------------------|
| "GLSURF"  | Draw a surface.|
| "GLSURF1" | Surface with color levels|
| "GLSURF2" | The same as "GLSURF1" but without polygon outlines.|
| "GLSURF3" | Color level projection on top of plot (works only in cartesian coordinate system).|
| "GLSURF4" | Same as "GLSURF" but without polygon outlines.|



The surface painting in cartesian coordinates supports logarithmic scales along
X, Y, Z axis. In polar coordinates only the Z axis can be logarithmic,
in cylindrical coordinates only the Y axis.

Additional options to SURF and LEGO - Coordinate systems:

| Option   | Description                                                       |
|----------|-------------------------------------------------------------------|
| " "      | Default, cartesian coordinates system.|
| "POL"    | Polar coordinates system.|
| "CYL"    | Cylindrical coordinates system.|
| "SPH"    | Spherical coordinates system.|



#### <a name="HP290"></a> TH3 as color boxes

The supported option is:

| Option   | Description                                                       |
|----------|-------------------------------------------------------------------|
| "GLCOL"  | H3 is drawn using semi-transparent colored boxes.  See `$ROOTSYS/tutorials/gl/glvox1.C`.|



#### <a name="HP29b"></a> TH3 as boxes (spheres)

The supported options are:

| Option   | Description                                                       |
|----------|-------------------------------------------------------------------|
| "GLBOX"  | TH3 as a set of boxes, size of box is proportional to bin content.|
| "GLBOX1" | The same as "glbox", but spheres are drawn instead of boxes.|



#### <a name="HP29c"></a> TH3 as iso-surface(s)

The supported option is:

| Option   | Description                                                       |
|----------|-------------------------------------------------------------------|
| "GLISO"  | TH3 is drawn using iso-surfaces.|



#### <a name="HP29d"></a> TF3 (implicit function)

The supported option is:

| Option   | Description                                                       |
|----------|-------------------------------------------------------------------|
| "GL"     | Draw a TF3.|



#### <a name="HP29e"></a> Parametric surfaces

`$ROOTSYS/tutorials/gl/glparametric.C` shows how to create parametric
equations and visualize the surface.

#### <a name="HP29f"></a> Interaction with the plots

All the interactions are implemented via standard methods
`DistancetoPrimitive()` and `ExecuteEvent()`. That's why all the
interactions with the OpenGL plots are possible only when the mouse cursor is
in the plot's area (the plot's area is the part of a the pad occupied by
gl-produced picture). If the mouse cursor is not above gl-picture, the standard
pad interaction is performed.

#### <a name="HP29g"></a> Selectable parts

Different parts of the plot can be selected:

-  xoz, yoz, xoy back planes: When such a plane selected, it's highlighted in green
   if the dynamic slicing by this plane is supported, and it's highlighted in red,
   if the dynamic slicing is not supported.
-  The plot itself:
   On surfaces, the selected surface is outlined in red. (TF3 and
   ISO are not outlined). On lego plots, the selected bin is
   highlighted. The bin number and content are displayed in pad's
   status bar. In box plots, the box or sphere is highlighted and
   the bin info is displayed in pad's status bar.


#### <a name="HP29h"></a> Rotation and zooming


- Rotation:
  When the plot is selected, it can be rotated by pressing and
  holding the left mouse button and move the cursor.
- Zoom/Unzoom:
  Mouse wheel or 'j', 'J', 'k', 'K' keys.


#### <a name="HP29i"></a> Panning

The selected plot can be moved in a pad's area by pressing and
holding the left mouse button and the shift key.

#### <a name="HP29j"></a> Box cut

Surface, iso, box, TF3 and parametric painters support box cut by
pressing the 'c' or 'C' key when the mouse cursor is in a plot's
area. That will display a transparent box, cutting away part of the
surface (or boxes) in order to show internal part of plot. This box
can be moved inside the plot's area (the full size of the box is
equal to the plot's surrounding box) by selecting one of the box
cut axes and pressing the left mouse button to move it.

#### <a name="HP29k"></a> Plot specific interactions (dynamic slicing etc.)

Currently, all gl-plots support some form of slicing. When back plane
is selected (and if it's highlighted in green) you can press and hold
left mouse button and shift key and move this back plane inside
plot's area, creating the slice. During this "slicing" plot becomes
semi-transparent. To remove all slices (and projected curves for
surfaces) double click with left mouse button in a plot's area.

#### <a name="HP29l"></a> Surface with option "GLSURF"

The surface profile is displayed on the slicing plane.
The profile projection is drawn on the back plane
by pressing `'p'` or `'P'` key.

#### <a name="HP29m"></a> TF3

The contour plot is drawn on the slicing plane. For TF3 the color
scheme can be changed by pressing 's' or 'S'.

#### <a name="HP29n"></a> Box

The contour plot corresponding to slice plane position is drawn in real time.

#### <a name="HP29o"></a> Iso

Slicing is similar to "GLBOX" option.

#### <a name="HP29p"></a> Parametric plot

No slicing. Additional keys: 's' or 'S' to change color scheme -
about 20 color schemes supported ('s' for "scheme"); 'l' or 'L' to
increase number of polygons ('l' for "level" of details), 'w' or 'W'
to show outlines ('w' for "wireframe").

#### <a name="HP30"></a> Highlight mode for histogram

\since **ROOT version 6.15/01**

\image html hlHisto3_top.gif "Highlight mode"

Highlight mode is implemented for `TH1` (and for `TGraph`) class. When
highlight mode is on, mouse movement over the bin will be represented
graphically. Bin will be highlighted as "bin box" (presented by box
object). Moreover, any highlight (change of bin) emits signal
`TCanvas::Highlighted()` which allows the user to react and call their own
function. For a better understanding see also the tutorials
`$ROOTSYS/tutorials/hist/hlHisto*.C` files.

Highlight mode is switched on/off by `TH1::SetHighlight()` function
or interactively from `TH1` context menu. `TH1::IsHighlight()` to verify
whether the highlight mode enabled or disabled, default it is disabled.

~~~ {.cpp}
    root [0] .x $ROOTSYS/tutorials/hsimple.C
    root [1] hpx->SetHighlight(kTRUE)   // or interactively from TH1 context menu
    root [2] hpx->IsHighlight()
    (bool) true
~~~

\image html hlsimple_nofun.gif "Highlight mode for histogram"

#### <a name="HP30a"></a> Highlight mode and user function

The user can use (connect) `TCanvas::Highlighted()` signal, which is always
emitted if there is a highlight bin and call user function via signal
and slot communication mechanism. `TCanvas::Highlighted()` is similar
`TCanvas::Picked()`

-  when selected object (histogram as a whole) is different from previous
then emit `Picked()` signal
-  when selected (highlighted) bin from histogram is different from previous
then emit `Highlighted()` signal

Any user function (or functions) has to be defined
`UserFunction(TVirtualPad *pad, TObject *obj, Int_t x, Int_t y)`.
In example (see below) has name `PrintInfo()`. All parameters of user
function are taken from

    void TCanvas::Highlighted(TVirtualPad *pad, TObject *obj, Int_t x, Int_t y)

-  `pad` is pointer to pad with highlighted histogram
-  `obj` is pointer to highlighted histogram
-  `x` is highlighted x bin for 1D histogram
-  `y` is highlighted y bin for 2D histogram (for 1D histogram not in use)

Example how to create a connection from any `TCanvas` object to a user
`UserFunction()` slot (see also `TQObject::Connect()` for additional info)

    TQObject::Connect("TCanvas", "Highlighted(TVirtualPad*,TObject*,Int_t,Int_t)",
                          0, 0, "UserFunction(TVirtualPad*,TObject*,Int_t,Int_t)");

or use non-static "simplified" function
`TCanvas::HighlightConnect(const char *slot)`

    c1->HighlightConnect("UserFunction(TVirtualPad*,TObject*,Int_t,Int_t)");

NOTE the signal and slot string must have a form
"(TVirtualPad*,TObject*,Int_t,Int_t)"

    root [0] .x $ROOTSYS/tutorials/hsimple.C
    root [1] hpx->SetHighlight(kTRUE)
    root [2] .x hlprint.C

file `hlprint.C`
~~~ {.cpp}
void PrintInfo(TVirtualPad *pad, TObject *obj, Int_t x, Int_t y)
{
   auto h = (TH1F *)obj;
   if (!h->IsHighlight()) // after highlight disabled
      h->SetTitle("highlight disable");
   else
      h->SetTitle(TString::Format("bin[%03d] (%5.2f) content %g", x,
                                  h->GetBinCenter(x), h->GetBinContent(x)));
   pad->Update();
}

void hlprint()
{
   if (!gPad) return;
   gPad->GetCanvas()->HighlightConnect("PrintInfo(TVirtualPad*,TObject*,Int_t,Int_t)");
}
~~~

\image html hlsimple.gif "Highlight mode and simple user function"

For more complex demo please see for example `$ROOTSYS/tutorials/tree/temperature.C` file.

*/

TH1 *gCurrentHist = 0;

Hoption_t Hoption;
Hparam_t  Hparam;

const Int_t kNMAX = 2000;

const Int_t kMAXCONTOUR  = 104;
const UInt_t kCannotRotate = BIT(11);

static TBox *gXHighlightBox = 0;   // highlight X box
static TBox *gYHighlightBox = 0;   // highlight Y box

static TString gStringEntries;
static TString gStringMean;
static TString gStringMeanX;
static TString gStringMeanY;
static TString gStringMeanZ;
static TString gStringStdDev;
static TString gStringStdDevX;
static TString gStringStdDevY;
static TString gStringStdDevZ;
static TString gStringUnderflow;
static TString gStringOverflow;
static TString gStringIntegral;
static TString gStringIntegralBinWidth;
static TString gStringSkewness;
static TString gStringSkewnessX;
static TString gStringSkewnessY;
static TString gStringSkewnessZ;
static TString gStringKurtosis;
static TString gStringKurtosisX;
static TString gStringKurtosisY;
static TString gStringKurtosisZ;

ClassImp(THistPainter);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

THistPainter::THistPainter()
{

   fH = 0;
   fXaxis = 0;
   fYaxis = 0;
   fZaxis = 0;
   fFunctions = 0;
   fXbuf  = 0;
   fYbuf  = 0;
   fNcuts = 0;
   fStack = 0;
   fLego  = 0;
   fPie   = 0;
   fGraph2DPainter = 0;
   fShowProjection = 0;
   fShowOption = "";
   for (int i=0; i<kMaxCuts; i++) {
      fCuts[i] = 0;
      fCutsOpt[i] = 0;
   }
   fXHighlightBin = -1;
   fYHighlightBin = -1;
   fCurrentF3 = nullptr;

   gStringEntries          = gEnv->GetValue("Hist.Stats.Entries",          "Entries");
   gStringMean             = gEnv->GetValue("Hist.Stats.Mean",             "Mean");
   gStringMeanX            = gEnv->GetValue("Hist.Stats.MeanX",            "Mean x");
   gStringMeanY            = gEnv->GetValue("Hist.Stats.MeanY",            "Mean y");
   gStringMeanZ            = gEnv->GetValue("Hist.Stats.MeanZ",            "Mean z");
   gStringStdDev           = gEnv->GetValue("Hist.Stats.StdDev",           "Std Dev");
   gStringStdDevX          = gEnv->GetValue("Hist.Stats.StdDevX",          "Std Dev x");
   gStringStdDevY          = gEnv->GetValue("Hist.Stats.StdDevY",          "Std Dev y");
   gStringStdDevZ          = gEnv->GetValue("Hist.Stats.StdDevZ",          "Std Dev z");
   gStringUnderflow        = gEnv->GetValue("Hist.Stats.Underflow",        "Underflow");
   gStringOverflow         = gEnv->GetValue("Hist.Stats.Overflow",         "Overflow");
   gStringIntegral         = gEnv->GetValue("Hist.Stats.Integral",         "Integral");
   gStringIntegralBinWidth = gEnv->GetValue("Hist.Stats.IntegralBinWidth", "Integral(w)");
   gStringSkewness         = gEnv->GetValue("Hist.Stats.Skewness",         "Skewness");
   gStringSkewnessX        = gEnv->GetValue("Hist.Stats.SkewnessX",        "Skewness x");
   gStringSkewnessY        = gEnv->GetValue("Hist.Stats.SkewnessY",        "Skewness y");
   gStringSkewnessZ        = gEnv->GetValue("Hist.Stats.SkewnessZ",        "Skewness z");
   gStringKurtosis         = gEnv->GetValue("Hist.Stats.Kurtosis",         "Kurtosis");
   gStringKurtosisX        = gEnv->GetValue("Hist.Stats.KurtosisX",        "Kurtosis x");
   gStringKurtosisY        = gEnv->GetValue("Hist.Stats.KurtosisY",        "Kurtosis y");
   gStringKurtosisZ        = gEnv->GetValue("Hist.Stats.KurtosisZ",        "Kurtosis z");
}

////////////////////////////////////////////////////////////////////////////////
/// Default destructor.

THistPainter::~THistPainter()
{
   if (fPie) delete fPie;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the distance from the point px,py to a line.
///
/// Compute the closest distance of approach from point px,py to elements of
/// an histogram. The distance is computed in pixels units.
///
/// Algorithm: Currently, this simple model computes the distance from the mouse
/// to the histogram contour only.

Int_t THistPainter::DistancetoPrimitive(Int_t px, Int_t py)
{

   Double_t defaultLabelSize = 0.04; // See TAttAxis.h for source of this value

   const Int_t big = 9999;
   const Int_t kMaxDiff = 7;

   if (fPie) return fPie->DistancetoPrimitive(px, py);

   Double_t x  = gPad->AbsPixeltoX(px);
   Double_t x1 = gPad->AbsPixeltoX(px+1);

   Int_t puxmin = gPad->XtoAbsPixel(gPad->GetUxmin());
   Int_t puymin = gPad->YtoAbsPixel(gPad->GetUymin());
   Int_t puxmax = gPad->XtoAbsPixel(gPad->GetUxmax());
   Int_t puymax = gPad->YtoAbsPixel(gPad->GetUymax());
   Int_t curdist = big;
   Int_t yxaxis, dyaxis,xyaxis, dxaxis;
   Bool_t dsame;
   TObject *PadPointer = gPad->GetPadPointer();
   if (!PadPointer) return 0;
   TString doption = PadPointer->GetDrawOption();
   Double_t factor = 1;
   if (fH->GetNormFactor() != 0) {
      factor = fH->GetNormFactor()/fH->GetSumOfWeights();
   }
   //     return if point is not in the histogram area

   //     If a 3D view exists, check distance to axis
   TView *view = gPad->GetView();
   Int_t d1,d2,d3;
   if (view && Hoption.Contour != 14) {
      Double_t ratio;
      d3 = view->GetDistancetoAxis(3, px, py, ratio);
      if (d3 <= kMaxDiff) {gPad->SetSelected(fZaxis); return 0;}
      d1 = view->GetDistancetoAxis(1, px, py, ratio);
      if (d1 <= kMaxDiff) {gPad->SetSelected(fXaxis); return 0;}
      d2 = view->GetDistancetoAxis(2, px, py, ratio);
      if (d2 <= kMaxDiff) {gPad->SetSelected(fYaxis); return 0;}
      if ( px > puxmin && px < puxmax && py > puymax && py < puymin) curdist = 1;
      goto FUNCTIONS;
   }
   //     check if point is close to an axis
   doption.ToLower();
   dsame = kFALSE;
   if (doption.Contains("same")) dsame = kTRUE;

   dyaxis = Int_t(2*(puymin-puymax)*TMath::Max(Double_t(fYaxis->GetLabelSize()), defaultLabelSize));
   if (doption.Contains("y+")) {
      xyaxis = puxmax + Int_t((puxmax-puxmin)*fYaxis->GetLabelOffset());
      if (px <= xyaxis+dyaxis && px >= xyaxis && py >puymax && py < puymin) {
         if (!dsame) {
            if (gPad->IsVertical()) gPad->SetSelected(fYaxis);
            else                    gPad->SetSelected(fXaxis);
            return 0;
         }
      }
   } else {
      xyaxis = puxmin - Int_t((puxmax-puxmin)*fYaxis->GetLabelOffset());
      if (px >= xyaxis-dyaxis && px <= xyaxis && py >puymax && py < puymin) {
         if (!dsame) {
            if (gPad->IsVertical()) gPad->SetSelected(fYaxis);
            else                    gPad->SetSelected(fXaxis);
            return 0;
         }
      }
   }

   dxaxis = Int_t((puymin-puymax)*TMath::Max(Double_t(fXaxis->GetLabelSize()), defaultLabelSize));
   if (doption.Contains("x+")) {
      yxaxis = puymax - Int_t((puymin-puymax)*fXaxis->GetLabelOffset());
      if (py >= yxaxis-dxaxis && py <= yxaxis && px <puxmax && px > puxmin) {
         if (!dsame) {
            if (gPad->IsVertical()) gPad->SetSelected(fXaxis);
            else                    gPad->SetSelected(fYaxis);
            return 0;
         }
      }
   } else {
      yxaxis = puymin + Int_t((puymin-puymax)*fXaxis->GetLabelOffset());
      if (yxaxis < puymin) yxaxis = puymin;
      if (py <= yxaxis+dxaxis && py >= yxaxis && px <puxmax && px > puxmin) {
         if (!dsame) {
            if (gPad->IsVertical()) gPad->SetSelected(fXaxis);
            else                    gPad->SetSelected(fYaxis);
            return 0;
         }
      }
   }

   if (fH->IsHighlight()) { // only if highlight is enable
      if ((px > puxmin) && (py < puymin) && (px < puxmax) && (py > puymax))
         HighlightBin(px, py);
   }

   //     if object is 2D or 3D return this object
   if (fH->GetDimension() == 2) {
      if (fH->InheritsFrom(TH2Poly::Class())) {
         TH2Poly *th2 = (TH2Poly*)fH;
         Double_t xmin, ymin, xmax, ymax;
         gPad->GetRangeAxis(xmin, ymin, xmax, ymax);
         Double_t pxu = gPad->AbsPixeltoX(px);
         Double_t pyu = gPad->AbsPixeltoY(py);
         if ((pxu>xmax) || (pxu < xmin) || (pyu>ymax) || (pyu < ymin)) {
            curdist = big;
            goto FUNCTIONS;
         } else {
            Int_t bin = th2->FindBin(pxu, pyu);
            if (bin>0) curdist = 1;
            else       curdist = big;
            goto FUNCTIONS;
         }
      }
      Int_t delta2 = 5; //Give a margin of delta2 pixels to be in the 2-d area
      if ( px > puxmin + delta2
        && px < puxmax - delta2
        && py > puymax + delta2
        && py < puymin - delta2) {curdist =1; goto FUNCTIONS;}
   }

   //     point is inside histogram area. Find channel number
   if (gPad->IsVertical()) {
      Int_t bin      = fXaxis->FindFixBin(gPad->PadtoX(x));
      Int_t binsup   = fXaxis->FindFixBin(gPad->PadtoX(x1));
      Double_t binval = factor*fH->GetBinContent(bin);
      Int_t pybin    = gPad->YtoAbsPixel(gPad->YtoPad(binval));
      if (binval == 0 && pybin < puymin) pybin = 10000;
      // special case if more than one bin for the pixel
      if (binsup-bin>1) {
         Double_t binvalmin, binvalmax;
         binvalmin=binval;
         binvalmax=binval;
         for (Int_t ibin=bin+1; ibin<binsup; ibin++) {
            Double_t binvaltmp = factor*fH->GetBinContent(ibin);
            if (binvalmin>binvaltmp) binvalmin=binvaltmp;
            if (binvalmax<binvaltmp) binvalmax=binvaltmp;
         }
         Int_t pybinmin = gPad->YtoAbsPixel(gPad->YtoPad(binvalmax));
         Int_t pybinmax = gPad->YtoAbsPixel(gPad->YtoPad(binvalmin));
         if (py<pybinmax+kMaxDiff/2 && py>pybinmin-kMaxDiff/2) pybin = py;
      }
      if (bin != binsup) { // Mouse on bin border
        Double_t binsupval = factor*fH->GetBinContent(binsup);
        Int_t pybinsub     = gPad->YtoAbsPixel(gPad->YtoPad(binsupval));
        if (py <= TMath::Max(pybinsub,pybin) && py >= TMath::Min(pybinsub,pybin) && pybin != 10000) return 0;
      }
      if (TMath::Abs(py - pybin) <= kMaxDiff) return TMath::Abs(py - pybin);
   } else {
      Double_t y  = gPad->AbsPixeltoY(py);
      Double_t y1 = gPad->AbsPixeltoY(py+1);
      Int_t bin      = fXaxis->FindFixBin(gPad->PadtoY(y));
      Int_t binsup   = fXaxis->FindFixBin(gPad->PadtoY(y1));
      Double_t binval = factor*fH->GetBinContent(bin);
      Int_t pxbin    = gPad->XtoAbsPixel(gPad->XtoPad(binval));
      if (binval == 0 && pxbin > puxmin) pxbin = 10000;
      // special case if more than one bin for the pixel
      if (binsup-bin>1) {
         Double_t binvalmin, binvalmax;
         binvalmin=binval;
         binvalmax=binval;
         for (Int_t ibin=bin+1; ibin<binsup; ibin++) {
            Double_t binvaltmp = factor*fH->GetBinContent(ibin);
            if (binvalmin>binvaltmp) binvalmin=binvaltmp;
            if (binvalmax<binvaltmp) binvalmax=binvaltmp;
         }
         Int_t pxbinmin = gPad->XtoAbsPixel(gPad->XtoPad(binvalmax));
         Int_t pxbinmax = gPad->XtoAbsPixel(gPad->XtoPad(binvalmin));
         if (px<pxbinmax+kMaxDiff/2 && px>pxbinmin-kMaxDiff/2) pxbin = px;
      }
      if (TMath::Abs(px - pxbin) <= kMaxDiff) return TMath::Abs(px - pxbin);
   }
   //     Loop on the list of associated functions and user objects
FUNCTIONS:
   TObject *f;
   TIter   next(fFunctions);
   while ((f = (TObject*) next())) {
      Int_t dist;
      if (f->InheritsFrom(TF1::Class())) dist = f->DistancetoPrimitive(-px,py);
      else                               dist = f->DistancetoPrimitive(px,py);
      if (dist < kMaxDiff) {gPad->SetSelected(f); return dist;}
   }
   return curdist;
}

////////////////////////////////////////////////////////////////////////////////
/// Display a panel with all histogram drawing options.

void THistPainter::DrawPanel()
{

   gCurrentHist = fH;
   if (!gPad) {
      Error("DrawPanel", "need to draw histogram first");
      return;
   }
   TVirtualPadEditor *editor = TVirtualPadEditor::GetPadEditor();
   editor->Show();
   gROOT->ProcessLine(Form("((TCanvas*)0x%lx)->Selected((TVirtualPad*)0x%lx,(TObject*)0x%lx,1)",
                           (ULong_t)gPad->GetCanvas(), (ULong_t)gPad, (ULong_t)fH));
}

////////////////////////////////////////////////////////////////////////////////
/// Execute the actions corresponding to `event`.
///
/// This function is called when a histogram is clicked with the locator at
/// the pixel position px,py.

void THistPainter::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{

   if (!gPad) return;

   static Int_t bin, px1, py1, px2, py2, pyold;
   static TBox *zoombox;
   Double_t zbx1,zbx2,zby1,zby2;

   Int_t bin1, bin2;
   Double_t xlow, xup, ylow, binval, x, baroffset, barwidth, binwidth;
   Bool_t opaque  = gPad->OpaqueMoving();

   if (!gPad->IsEditable()) return;

   if (fPie) {
      fPie->ExecuteEvent(event, px, py);
      return;
   }
   //     come here if we have a lego/surface in the pad
   TView *view = gPad->GetView();

   if (!fShowProjection && view && view->TestBit(kCannotRotate) == 0) {
      view->ExecuteRotateView(event, px, py);
      return;
   }

   TAxis *xaxis    = fH->GetXaxis();
   TAxis *yaxis    = fH->GetYaxis();
   Int_t dimension = fH->GetDimension();

   // In case of option SAME the axis must be the ones of the first drawn histogram
   TString IsSame = fH->GetDrawOption();
   IsSame.ToLower();
   if (IsSame.Index("same")>=0) {
      TH1 *h1;
      TIter next(gPad->GetListOfPrimitives());
      while ((h1 = (TH1 *)next())) {
         if (!h1->InheritsFrom(TH1::Class())) continue;
         xaxis    = h1->GetXaxis();
         yaxis    = h1->GetYaxis();
         break;
      }
   }

   Double_t factor = 1;
   if (fH->GetNormFactor() != 0) {
      factor = fH->GetNormFactor()/fH->GetSumOfWeights();
   }

   switch (event) {

   case kButton1Down:

      if (!opaque) gVirtualX->SetLineColor(-1);
      fH->TAttLine::Modify();

      if (opaque && dimension ==2) {
         zbx1 = gPad->AbsPixeltoX(px);
         zbx2 = gPad->AbsPixeltoX(px);
         zby1 = gPad->AbsPixeltoY(py);
         zby2 = gPad->AbsPixeltoY(py);
         px1 = px;
         py1 = py;
         if (gPad->GetLogx()) {
            zbx1 = TMath::Power(10,zbx1);
            zbx2 = TMath::Power(10,zbx2);
         }
         if (gPad->GetLogy()) {
            zby1 = TMath::Power(10,zby1);
            zby2 = TMath::Power(10,zby2);
         }
         zoombox = new TBox(zbx1, zby1, zbx2, zby2);
         Int_t ci = TColor::GetColor("#7d7dff");
         TColor *zoomcolor = gROOT->GetColor(ci);
         if (!TCanvas::SupportAlpha() || !zoomcolor) zoombox->SetFillStyle(3002);
         else                                        zoomcolor->SetAlpha(0.5);
         zoombox->SetFillColor(ci);
         zoombox->Draw();
         gPad->Modified();
         gPad->Update();
      }
      // No break !!!

   case kMouseMotion:

      if (fShowProjection) {ShowProjection3(px,py); break;}

      gPad->SetCursor(kPointer);
      if (dimension ==1) {
         if (Hoption.Bar) {
            baroffset = fH->GetBarOffset();
            barwidth  = fH->GetBarWidth();
         } else {
            baroffset = 0;
            barwidth  = 1;
         }
         x        = gPad->AbsPixeltoX(px);
         bin      = fXaxis->FindFixBin(gPad->PadtoX(x));
         binwidth = fXaxis->GetBinWidth(bin);
         xlow     = gPad->XtoPad(fXaxis->GetBinLowEdge(bin) + baroffset*binwidth);
         xup      = gPad->XtoPad(xlow + barwidth*binwidth);
         ylow     = gPad->GetUymin();
         px1      = gPad->XtoAbsPixel(xlow);
         px2      = gPad->XtoAbsPixel(xup);
         py1      = gPad->YtoAbsPixel(ylow);
         py2      = py;
         pyold    = py;
         if (gROOT->GetEditHistograms()) gPad->SetCursor(kArrowVer);
      }

      break;

   case kButton1Motion:

      if (dimension ==1) {
         if (gROOT->GetEditHistograms()) {
            if (!opaque) {
               gVirtualX->DrawBox(px1, py1, px2, py2,TVirtualX::kHollow);  // Draw the old box
               py2 += py - pyold;
               gVirtualX->DrawBox(px1, py1, px2, py2,TVirtualX::kHollow);  // Draw the new box
               pyold = py;
            } else {
               py2 += py - pyold;
               pyold = py;
               binval = gPad->PadtoY(gPad->AbsPixeltoY(py2))/factor;
               fH->SetBinContent(bin,binval);
               gPad->Modified(kTRUE);
            }
         }
      }

      if (opaque && dimension ==2) {
         if (TMath::Abs(px1-px)>5 && TMath::Abs(py1-py)>5) {
            zbx2 = gPad->AbsPixeltoX(px);
            zby2 = gPad->AbsPixeltoY(py);
            if (gPad->GetLogx()) zbx2 = TMath::Power(10,zbx2);
            if (gPad->GetLogy()) zby2 = TMath::Power(10,zby2);
            zoombox->SetX2(zbx2);
            zoombox->SetY2(zby2);
            gPad->Modified();
            gPad->Update();
         }
      }

      break;

   case kWheelUp:

      if (dimension ==2) {
         bin1 = xaxis->GetFirst()+1;
         bin2 = xaxis->GetLast()-1;
         bin1 = TMath::Max(bin1, 1);
         bin2 = TMath::Min(bin2, xaxis->GetNbins());
         if (bin2>bin1) xaxis->SetRange(bin1,bin2);
         bin1 = yaxis->GetFirst()+1;
         bin2 = yaxis->GetLast()-1;
         bin1 = TMath::Max(bin1, 1);
         bin2 = TMath::Min(bin2, yaxis->GetNbins());
         if (bin2>bin1) yaxis->SetRange(bin1,bin2);
      }
      gPad->Modified();
      gPad->Update();

      break;

   case kWheelDown:

      if (dimension == 2) {
         bin1 = xaxis->GetFirst()-1;
         bin2 = xaxis->GetLast()+1;
         bin1 = TMath::Max(bin1, 1);
         bin2 = TMath::Min(bin2, xaxis->GetNbins());
         if (bin2>bin1) xaxis->SetRange(bin1,bin2);
         bin1 = yaxis->GetFirst()-1;
         bin2 = yaxis->GetLast()+1;
         bin1 = TMath::Max(bin1, 1);
         bin2 = TMath::Min(bin2, yaxis->GetNbins());
         if (bin2>bin1) yaxis->SetRange(bin1,bin2);
      }
      gPad->Modified();
      gPad->Update();

      break;

   case kButton1Up:
      if (dimension ==1) {
         if (gROOT->GetEditHistograms()) {
            binval = gPad->PadtoY(gPad->AbsPixeltoY(py2))/factor;
            fH->SetBinContent(bin,binval);
            PaintInit();   // recalculate Hparam structure and recalculate range
         }

         // might resize pad pixmap so should be called before any paint routine
         RecalculateRange();
      }
      if (opaque && dimension ==2) {
         if (zoombox) {
            Double_t x1 = TMath::Min(zoombox->GetX1(), zoombox->GetX2());
            Double_t x2 = TMath::Max(zoombox->GetX1(), zoombox->GetX2());
            Double_t y1 = TMath::Min(zoombox->GetY1(), zoombox->GetY2());
            Double_t y2 = TMath::Max(zoombox->GetY1(), zoombox->GetY2());
            x1 = TMath::Max(x1,xaxis->GetXmin());
            x2 = TMath::Min(x2,xaxis->GetXmax());
            y1 = TMath::Max(y1,yaxis->GetXmin());
            y2 = TMath::Min(y2,yaxis->GetXmax());
            if (x1<x2 && y1<y2) {
               xaxis->SetRangeUser(x1, x2);
               yaxis->SetRangeUser(y1, y2);
            }
            zoombox->Delete();
            zoombox = 0;
         }
      }
      gPad->Modified(kTRUE);
      if (opaque) gVirtualX->SetLineColor(-1);

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
/// Get a contour (as a list of TGraphs) using the Delaunay triangulation.

TList *THistPainter::GetContourList(Double_t contour) const
{



   // Check if fH contains a TGraphDelaunay2D
   TList *hl = fH->GetListOfFunctions();
   TGraphDelaunay2D *dt = (TGraphDelaunay2D*)hl->FindObject("TGraphDelaunay2D");
   // try with the old painter
   TGraphDelaunay *dtOld = nullptr;
   if (!dt) dtOld =  (TGraphDelaunay*)hl->FindObject("TGraphDelaunay");

   if (!dt && !dtOld) return nullptr;

   gCurrentHist = fH;

   if (!fGraph2DPainter) {
      if (dt) ((THistPainter*)this)->fGraph2DPainter = new TGraph2DPainter(dt);
      else ((THistPainter*)this)->fGraph2DPainter = new TGraph2DPainter(dtOld);
   }

   return fGraph2DPainter->GetContourList(contour);
}

////////////////////////////////////////////////////////////////////////////////
/// Display the histogram info (bin number, contents, integral up to bin
/// corresponding to cursor position px,py.

char *THistPainter::GetObjectInfo(Int_t px, Int_t py) const
{

   if (!gPad) return (char*)"";

   Double_t x  = gPad->PadtoX(gPad->AbsPixeltoX(px));
   Double_t y  = gPad->PadtoY(gPad->AbsPixeltoY(py));
   Double_t x1 = gPad->PadtoX(gPad->AbsPixeltoX(px+1));
   TString drawOption = fH->GetDrawOption();
   drawOption.ToLower();
   Double_t xmin, xmax, uxmin,uxmax;
   Double_t ymin, ymax, uymin,uymax;
   if (fH->GetDimension() == 2) {
      if (gPad->GetView() || drawOption.Index("cont") >= 0) {
         uxmin=gPad->GetUxmin();
         uxmax=gPad->GetUxmax();
         xmin = fXaxis->GetBinLowEdge(fXaxis->GetFirst());
         xmax = fXaxis->GetBinUpEdge(fXaxis->GetLast());
         x = xmin +(xmax-xmin)*(x-uxmin)/(uxmax-uxmin);
         uymin=gPad->GetUymin();
         uymax=gPad->GetUymax();
         ymin = fYaxis->GetBinLowEdge(fYaxis->GetFirst());
         ymax = fYaxis->GetBinUpEdge(fYaxis->GetLast());
         y = ymin +(ymax-ymin)*(y-uymin)/(uymax-uymin);
      }
   }
   Int_t binx,biny,binmin=0,binx1;
   if (gPad->IsVertical()) {
      binx   = fXaxis->FindFixBin(x);
      if (drawOption.Index("same") >= 0) {
         TH1 *h1;
         TIter next(gPad->GetListOfPrimitives());
         while ((h1 = (TH1 *)next())) {
            if (!h1->InheritsFrom(TH1::Class())) continue;
            binmin = h1->GetXaxis()->GetFirst();
            break;
         }
      } else {
         binmin = fXaxis->GetFirst();
      }
      binx1  = fXaxis->FindFixBin(x1);
      // special case if more than 1 bin in x per pixel
      if (binx1-binx>1 && fH->GetDimension() == 1) {
         Double_t binval=fH->GetBinContent(binx);
         Int_t binnear=binx;
         for (Int_t ibin=binx+1; ibin<binx1; ibin++) {
            Double_t binvaltmp = fH->GetBinContent(ibin);
            if (TMath::Abs(y-binvaltmp) < TMath::Abs(y-binval)) {
               binval=binvaltmp;
               binnear=ibin;
            }
         }
         binx = binnear;
      }
   } else {
      x1 = gPad->PadtoY(gPad->AbsPixeltoY(py+1));
      binx   = fXaxis->FindFixBin(y);
      if (drawOption.Index("same") >= 0) {
         TH1 *h1;
         TIter next(gPad->GetListOfPrimitives());
         while ((h1 = (TH1 *)next())) {
            if (!h1->InheritsFrom(TH1::Class())) continue;
            binmin = h1->GetXaxis()->GetFirst();
            break;
         }
      } else {
         binmin = fXaxis->GetFirst();
      }
      binx1  = fXaxis->FindFixBin(x1);
      // special case if more than 1 bin in x per pixel
      if (binx1-binx>1 && fH->GetDimension() == 1) {
         Double_t binval=fH->GetBinContent(binx);
         Int_t binnear=binx;
         for (Int_t ibin=binx+1; ibin<binx1; ibin++) {
            Double_t binvaltmp = fH->GetBinContent(ibin);
            if (TMath::Abs(x-binvaltmp) < TMath::Abs(x-binval)) {
               binval=binvaltmp;
               binnear=ibin;
            }
         }
         binx = binnear;
      }
   }
   if (fH->GetDimension() == 1) {
      if (fH->InheritsFrom(TProfile::Class())) {
         TProfile *tp = (TProfile*)fH;
         fObjectInfo.Form("(x=%g, y=%g, binx=%d, binc=%g, bine=%g, binn=%d)",
            x, y, binx, fH->GetBinContent(binx), fH->GetBinError(binx),
            (Int_t) tp->GetBinEntries(binx));
      }
      else {
         Double_t integ = 0;
         for (Int_t bin=binmin;bin<=binx;bin++) {integ += fH->GetBinContent(bin);}
         fObjectInfo.Form("(x=%g, y=%g, binx=%d, binc=%g, Sum=%g)",
            x,y,binx,fH->GetBinContent(binx),integ);
      }
   } else if (fH->GetDimension() == 2) {
      if (fH->InheritsFrom(TH2Poly::Class())) {
         TH2Poly *th2 = (TH2Poly*)fH;
         biny = th2->FindBin(x,y);
         fObjectInfo.Form("%s (x=%g, y=%g, bin=%d, binc=%g)",
            th2->GetBinTitle(biny),x,y,biny,th2->GetBinContent(biny));
      }
      else if (fH->InheritsFrom(TProfile2D::Class())) {
         TProfile2D *tp = (TProfile2D*)fH;
         biny = fYaxis->FindFixBin(y);
         Int_t bin = fH->GetBin(binx,biny);
         fObjectInfo.Form("(x=%g, y=%g, binx=%d, biny=%d, binc=%g, bine=%g, binn=%d)",
            x, y, binx, biny, fH->GetBinContent(bin),
            fH->GetBinError(bin), (Int_t) tp->GetBinEntries(bin));
      } else {
         biny = fYaxis->FindFixBin(y);
         fObjectInfo.Form("(x=%g, y=%g, binx=%d, biny=%d, binc=%g bine=%g)",
            x,y,binx,biny,fH->GetBinContent(binx,biny),
            fH->GetBinError(binx,biny));
      }
   } else {
      // 3d case: retrieving the x,y,z bin is not yet implemented
      // print just the x,y info
      fObjectInfo.Form("(x=%g, y=%g)",x,y);
   }

   return (char *)fObjectInfo.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Set highlight (enable/disable) mode for fH

void THistPainter::SetHighlight()
{
   if (fH->IsHighlight()) return;

   fXHighlightBin = -1;
   fYHighlightBin = -1;
   // delete previous highlight box
   if (gXHighlightBox) { gXHighlightBox->Delete(); gXHighlightBox = 0; }
   if (gYHighlightBox) { gYHighlightBox->Delete(); gYHighlightBox = 0; }
   // emit Highlighted() signal (user can check on disabled)
   if (gPad->GetCanvas()) gPad->GetCanvas()->Highlighted(gPad, fH, fXHighlightBin, fYHighlightBin);
}

////////////////////////////////////////////////////////////////////////////////
/// Check on highlight bin

void THistPainter::HighlightBin(Int_t px, Int_t py)
{
   // call from DistancetoPrimitive (only if highlight is enable)

   Double_t x = gPad->PadtoX(gPad->AbsPixeltoX(px));
   Double_t y = gPad->PadtoY(gPad->AbsPixeltoY(py));
   Int_t binx = fXaxis->FindFixBin(x);
   Int_t biny = fYaxis->FindFixBin(y);
   if (!gPad->IsVertical()) binx = fXaxis->FindFixBin(y);

   Bool_t changedBin = kFALSE;
   if (binx != fXHighlightBin) {
      fXHighlightBin = binx;
      changedBin = kTRUE;
   } else if (fH->GetDimension() == 1) return;
   if (biny != fYHighlightBin) {
      fYHighlightBin = biny;
      changedBin = kTRUE;
   }
   if (!changedBin) return;

   //   Info("HighlightBin", "histo: %p '%s'\txbin: %d, ybin: %d",
   //        (void *)fH, fH->GetName(), fXHighlightBin, fYHighlightBin);

   // paint highlight bin as box (recursive calls PaintHighlightBin)
   gPad->Modified(kTRUE);
   gPad->Update();

   // emit Highlighted() signal
   if (gPad->GetCanvas()) gPad->GetCanvas()->Highlighted(gPad, fH, fXHighlightBin, fYHighlightBin);
}

////////////////////////////////////////////////////////////////////////////////
/// Paint highlight bin as TBox object

void THistPainter::PaintHighlightBin(Option_t * /*option*/)
{
   // call from PaintTitle

   if (!fH->IsHighlight()) return;

   Double_t uxmin = gPad->GetUxmin();
   Double_t uxmax = gPad->GetUxmax();
   Double_t uymin = gPad->GetUymin();
   Double_t uymax = gPad->GetUymax();
   if (gPad->GetLogx()) {
      uxmin = TMath::Power(10.0, uxmin);
      uxmax = TMath::Power(10.0, uxmax);
   }
   if (gPad->GetLogy()) {
      uymin = TMath::Power(10.0, uymin);
      uymax = TMath::Power(10.0, uymax);
   }

   // testing specific possibility (after zoom, draw with "same", log, etc.)
   Double_t hcenter;
   if (gPad->IsVertical()) {
      hcenter = fXaxis->GetBinCenter(fXHighlightBin);
      if ((hcenter < uxmin) || (hcenter > uxmax)) return;
   } else {
      hcenter = fYaxis->GetBinCenter(fXHighlightBin);
      if ((hcenter < uymin) || (hcenter > uymax)) return;
   }
   if (fH->GetDimension() == 2) {
      hcenter = fYaxis->GetBinCenter(fYHighlightBin);
      if ((hcenter < uymin) || (hcenter > uymax)) return;
   }

   // paint X highlight bin (for 1D or 2D)
   Double_t hbx1, hbx2, hby1, hby2;
   if (gPad->IsVertical()) {
      hbx1 = fXaxis->GetBinLowEdge(fXHighlightBin);
      hbx2 = fXaxis->GetBinUpEdge(fXHighlightBin);
      hby1 = uymin;
      hby2 = uymax;
   } else {
      hbx1 = uxmin;
      hbx2 = uxmax;
      hby1 = fYaxis->GetBinLowEdge(fXHighlightBin);
      hby2 = fYaxis->GetBinUpEdge(fXHighlightBin);
   }

   if (!gXHighlightBox) {
      gXHighlightBox = new TBox(hbx1, hby1, hbx2, hby2);
      gXHighlightBox->SetBit(kCannotPick);
      gXHighlightBox->SetFillColor(TColor::GetColor("#9797ff"));
      if (!TCanvas::SupportAlpha()) gXHighlightBox->SetFillStyle(3001);
      else gROOT->GetColor(gXHighlightBox->GetFillColor())->SetAlpha(0.5);
   }
   gXHighlightBox->SetX1(hbx1);
   gXHighlightBox->SetX2(hbx2);
   gXHighlightBox->SetY1(hby1);
   gXHighlightBox->SetY2(hby2);
   gXHighlightBox->Paint();

   //   Info("PaintHighlightBin", "histo: %p '%s'\txbin: %d, ybin: %d",
   //        (void *)fH, fH->GetName(), fXHighlightBin, fYHighlightBin);

   // paint Y highlight bin (only for 2D)
   if (fH->GetDimension() != 2) return;
   hbx1 = uxmin;
   hbx2 = uxmax;
   hby1 = fYaxis->GetBinLowEdge(fYHighlightBin);
   hby2 = fYaxis->GetBinUpEdge(fYHighlightBin);

   if (!gYHighlightBox) {
      gYHighlightBox = new TBox(hbx1, hby1, hbx2, hby2);
      gYHighlightBox->SetBit(kCannotPick);
      gYHighlightBox->SetFillColor(gXHighlightBox->GetFillColor());
      gYHighlightBox->SetFillStyle(gXHighlightBox->GetFillStyle());
   }
   gYHighlightBox->SetX1(hbx1);
   gYHighlightBox->SetX2(hbx2);
   gYHighlightBox->SetY1(hby1);
   gYHighlightBox->SetY2(hby2);
   gYHighlightBox->Paint();
}

////////////////////////////////////////////////////////////////////////////////
/// Return `kTRUE` if the cell `ix`, `iy` is inside one of the graphical cuts.

Bool_t THistPainter::IsInside(Int_t ix, Int_t iy)
{

   for (Int_t i=0;i<fNcuts;i++) {
      Double_t x = fXaxis->GetBinCenter(ix);
      Double_t y = fYaxis->GetBinCenter(iy);
      if (fCutsOpt[i] > 0) {
         if (!fCuts[i]->IsInside(x,y)) return kFALSE;
      } else {
         if (fCuts[i]->IsInside(x,y))  return kFALSE;
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return `kTRUE` if the point `x`, `y` is inside one of the graphical cuts.

Bool_t THistPainter::IsInside(Double_t x, Double_t y)
{

   for (Int_t i=0;i<fNcuts;i++) {
      if (fCutsOpt[i] > 0) {
         if (!fCuts[i]->IsInside(x,y)) return kFALSE;
      } else {
         if (fCuts[i]->IsInside(x,y))  return kFALSE;
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Decode string `choptin` and fill Hoption structure.

Int_t THistPainter::MakeChopt(Option_t *choptin)
{

   char *l;
   char chopt[128];
   Int_t nch = strlen(choptin);
   strlcpy(chopt,choptin,128);
   Int_t hdim = fH->GetDimension();

   Hoption.Axis = Hoption.Bar    = Hoption.Curve   = Hoption.Error   = 0;
   Hoption.Hist = Hoption.Line   = Hoption.Mark    = Hoption.Fill    = 0;
   Hoption.Same = Hoption.Func   = Hoption.Scat                      = 0;
   Hoption.Star = Hoption.Arrow  = Hoption.Box     = Hoption.Text    = 0;
   Hoption.Char = Hoption.Color  = Hoption.Contour = Hoption.Logx    = 0;
   Hoption.Logy = Hoption.Logz   = Hoption.Lego    = Hoption.Surf    = 0;
   Hoption.Off  = Hoption.Tri    = Hoption.Proj    = Hoption.AxisPos = 0;
   Hoption.Spec = Hoption.Pie    = Hoption.Candle  = 0;

   //    special 2D options
   Hoption.List     = 0;
   Hoption.Zscale   = 0;
   Hoption.FrontBox = 1;
   Hoption.BackBox  = 1;
   Hoption.System   = kCARTESIAN;

   Hoption.Zero     = 0;

   Hoption.MinimumZero = gStyle->GetHistMinimumZero() ? 1 : 0;

   //check for graphical cuts
   MakeCuts(chopt);

   for (Int_t i=0;i<nch;i++) chopt[i] = toupper(chopt[i]);
   if (hdim > 1) Hoption.Scat = 1;
   if (!nch) Hoption.Hist = 1;
   if (fFunctions->First()) Hoption.Func = 1;
   if (fH->GetSumw2N() && hdim == 1) Hoption.Error = 2;

   char *l1 = strstr(chopt,"PFC"); // Automatic Fill Color
   char *l2 = strstr(chopt,"PLC"); // Automatic Line Color
   char *l3 = strstr(chopt,"PMC"); // Automatic Marker Color
   if (l1 || l2 || l3) {
      Int_t i = gPad->NextPaletteColor();
      if (l1) {memcpy(l1,"   ",3); fH->SetFillColor(i);}
      if (l2) {memcpy(l2,"   ",3); fH->SetLineColor(i);}
      if (l3) {memcpy(l3,"   ",3); fH->SetMarkerColor(i);}
      Hoption.Hist = 1; // Make sure something is drawn in case there is no drawing option specified.
   }

   l = strstr(chopt,"MIN0");
   if (l) {
      Hoption.MinimumZero = 1;
      memcpy(l,"    ",4);
   }

   l = strstr(chopt,"SPEC");
   if (l) {
      Hoption.Scat = 0;
      memcpy(l,"    ",4);
      Int_t bs=0;
      l = strstr(chopt,"BF(");
      if (l) {
         if (sscanf(&l[3],"%d",&bs) > 0) {
            Int_t i=0;
            while (l[i]!=')') {
               l[i] = ' ';
               i++;
            }
            l[i] = ' ';
         }
      }
      Hoption.Spec = TMath::Max(1600,bs);
      return 1;
   }

   l = strstr(chopt,"GL");
   if (l) {
      memcpy(l,"  ",2);
   }
   l = strstr(chopt,"X+");
   if (l) {
      Hoption.AxisPos = 10;
      memcpy(l,"  ",2);
   }
   l = strstr(chopt,"Y+");
   if (l) {
      Hoption.AxisPos += 1;
      memcpy(l,"  ",2);
   }
   if ((Hoption.AxisPos == 10 || Hoption.AxisPos == 1) && (nch == 2)) Hoption.Hist = 1;
   if (Hoption.AxisPos == 11 && nch == 4) Hoption.Hist = 1;

   l = strstr(chopt,"SAMES");
   if (l) {
      if (nch == 5) Hoption.Hist = 1;
      Hoption.Same = 2;
      memcpy(l,"     ",5);
      if (l[5] == '0') { Hoption.Same += 10; l[5] = ' '; }
   }
   l = strstr(chopt,"SAME");
   if (l) {
      if (nch == 4) Hoption.Hist = 1;
      Hoption.Same = 1;
      memcpy(l,"    ",4);
      if (l[4] == '0') { Hoption.Same += 10; l[4] = ' '; }
   }

   l = strstr(chopt,"PIE");
   if (l) {
      Hoption.Pie = 1;
      memcpy(l,"   ",3);
   }


   l = strstr(chopt,"CANDLE");
   if (l) {
      TCandle candle;
      Hoption.Candle = candle.ParseOption(l);
      Hoption.Scat = 0;
   }

   l = strstr(chopt,"VIOLIN");
   if (l) {
      TCandle candle;
      Hoption.Candle = candle.ParseOption(l);
      Hoption.Scat = 0;
   }

   l = strstr(chopt,"LEGO");
   if (l) {
      Hoption.Scat = 0;
      Hoption.Lego = 1; memcpy(l,"    ",4);
      if (l[4] == '1') { Hoption.Lego = 11; l[4] = ' '; }
      if (l[4] == '2') { Hoption.Lego = 12; l[4] = ' '; }
      if (l[4] == '3') { Hoption.Lego = 13; l[4] = ' '; }
      if (l[4] == '4') { Hoption.Lego = 14; l[4] = ' '; }
      if (l[4] == '9') { Hoption.Lego = 19; l[4] = ' '; }
      l = strstr(chopt,"FB"); if (l) { Hoption.FrontBox = 0; memcpy(l,"  ",2); }
      l = strstr(chopt,"BB"); if (l) { Hoption.BackBox = 0;  memcpy(l,"  ",2); }
      l = strstr(chopt,"0");  if (l) { Hoption.Zero = 1;  memcpy(l," ",1); }
   }

   l = strstr(chopt,"SURF");
   if (l) {
      Hoption.Scat = 0;
      Hoption.Surf = 1; memcpy(l,"    ",4);
      if (l[4] == '1') { Hoption.Surf = 11; l[4] = ' '; }
      if (l[4] == '2') { Hoption.Surf = 12; l[4] = ' '; }
      if (l[4] == '3') { Hoption.Surf = 13; l[4] = ' '; }
      if (l[4] == '4') { Hoption.Surf = 14; l[4] = ' '; }
      if (l[4] == '5') { Hoption.Surf = 15; l[4] = ' '; }
      if (l[4] == '6') { Hoption.Surf = 16; l[4] = ' '; }
      if (l[4] == '7') { Hoption.Surf = 17; l[4] = ' '; }
      l = strstr(chopt,"FB");   if (l) { Hoption.FrontBox = 0; memcpy(l,"  ",2); }
      l = strstr(chopt,"BB");   if (l) { Hoption.BackBox = 0;  memcpy(l,"  ",2); }
   }

   l = strstr(chopt,"TF3");
   if (l) {
      memcpy(l,"    ",3);
      l = strstr(chopt,"FB");   if (l) { Hoption.FrontBox = 0; memcpy(l,"  ",2); }
      l = strstr(chopt,"BB");   if (l) { Hoption.BackBox = 0;  memcpy(l,"  ",2); }
   }

   l = strstr(chopt,"ISO");
   if (l) {
      memcpy(l,"    ",3);
      l = strstr(chopt,"FB");   if (l) { Hoption.FrontBox = 0; memcpy(l,"  ",2); }
      l = strstr(chopt,"BB");   if (l) { Hoption.BackBox = 0;  memcpy(l,"  ",2); }
   }

   l = strstr(chopt,"LIST");    if (l) { Hoption.List = 1;  memcpy(l,"    ",4);}

   l = strstr(chopt,"CONT");
   if (l) {
      memcpy(l,"    ",4);
      if (hdim>1) {
         Hoption.Scat = 0;
         Hoption.Contour = 1;
         if (l[4] == '1') { Hoption.Contour = 11; l[4] = ' '; }
         if (l[4] == '2') { Hoption.Contour = 12; l[4] = ' '; }
         if (l[4] == '3') { Hoption.Contour = 13; l[4] = ' '; }
         if (l[4] == '4') { Hoption.Contour = 14; l[4] = ' '; }
         if (l[4] == '5') { Hoption.Contour = 15; l[4] = ' '; }
      } else {
         Hoption.Hist = 1;
      }
   }
   l = strstr(chopt,"HBAR");
   if (l) {
      Hoption.Hist = 0;
      Hoption.Bar = 20; memcpy(l,"    ",4);
      if (l[4] == '1') { Hoption.Bar = 21; l[4] = ' '; }
      if (l[4] == '2') { Hoption.Bar = 22; l[4] = ' '; }
      if (l[4] == '3') { Hoption.Bar = 23; l[4] = ' '; }
      if (l[4] == '4') { Hoption.Bar = 24; l[4] = ' '; }
   }
   l = strstr(chopt,"BAR");
   if (l) {
      Hoption.Hist = 0;
      Hoption.Bar = 10; memcpy(l,"   ",3);
      if (l[3] == '1') { Hoption.Bar = 11; l[3] = ' '; }
      if (l[3] == '2') { Hoption.Bar = 12; l[3] = ' '; }
      if (l[3] == '3') { Hoption.Bar = 13; l[3] = ' '; }
      if (l[3] == '4') { Hoption.Bar = 14; l[3] = ' '; }
   }

   l = strstr(chopt,"ARR" );
   if (l) {
      memcpy(l,"   ", 3);
      if (hdim>1) {
         Hoption.Arrow  = 1;
         Hoption.Scat = 0;
         l = strstr(chopt,"COL"); if (l) { Hoption.Arrow  = 2;  memcpy(l,"   ",3); }
         l = strstr(chopt,"Z");   if (l) { Hoption.Zscale = 1;  memcpy(l," ",1); }
      } else {
         Hoption.Hist = 1;
      }
   }
   l = strstr(chopt,"BOX" );
   if (l) {
      memcpy(l,"   ", 3);
      if (hdim>1) {
         Hoption.Scat = 0;
         Hoption.Box  = 1;
         if (l[3] == '1') { Hoption.Box = 11; l[3] = ' '; }
         if (l[3] == '2') { Hoption.Box = 12; l[3] = ' '; }
         if (l[3] == '3') { Hoption.Box = 13; l[3] = ' '; }
      } else {
         Hoption.Hist = 1;
      }
   }
   l = strstr(chopt,"COLZ");
   if (l) {
      memcpy(l,"    ",4);
      if (hdim>1) {
         Hoption.Color  = 1;
         Hoption.Scat   = 0;
         Hoption.Zscale = 1;
         if (l[4] == '2') { Hoption.Color = 3; l[4] = ' '; }
         l = strstr(chopt,"0");  if (l) { Hoption.Zero  = 1;  memcpy(l," ",1); }
         l = strstr(chopt,"1");  if (l) { Hoption.Color = 2;  memcpy(l," ",1); }
      } else {
         Hoption.Hist = 1;
      }
   }
   l = strstr(chopt,"COL" );
   if (l) {
      memcpy(l,"   ", 3);
      if (hdim>1) {
         Hoption.Color = 1;
         Hoption.Scat  = 0;
         if (l[3] == '2') { Hoption.Color = 3; l[3] = ' '; }
         l = strstr(chopt,"0");  if (l) { Hoption.Zero  = 1;  memcpy(l," ",1); }
         l = strstr(chopt,"1");  if (l) { Hoption.Color = 2;  memcpy(l," ",1); }
      } else {
         Hoption.Hist = 1;
      }
   }
   l = strstr(chopt,"CHAR"); if (l) { Hoption.Char   = 1; memcpy(l,"    ",4); Hoption.Scat = 0; }
   l = strstr(chopt,"FUNC"); if (l) { Hoption.Func   = 2; memcpy(l,"    ",4); Hoption.Hist = 0; }
   l = strstr(chopt,"HIST"); if (l) { Hoption.Hist   = 2; memcpy(l,"    ",4); Hoption.Func = 0; Hoption.Error = 0;}
   l = strstr(chopt,"AXIS"); if (l) { Hoption.Axis   = 1; memcpy(l,"    ",4); }
   l = strstr(chopt,"AXIG"); if (l) { Hoption.Axis   = 2; memcpy(l,"    ",4); }
   l = strstr(chopt,"SCAT"); if (l) { Hoption.Scat   = 1; memcpy(l,"    ",4); }
   l = strstr(chopt,"TEXT");
   if (l) {
      Int_t angle;
      if (sscanf(&l[4],"%d",&angle) > 0) {
         if (angle < 0)  angle=0;
         if (angle > 90) angle=90;
         Hoption.Text = 1000+angle;
      } else {
         Hoption.Text = 1;
      }
      memcpy(l,"    ", 4);
      l = strstr(chopt,"N");
      if (l && fH->InheritsFrom(TH2Poly::Class())) Hoption.Text += 3000;
      Hoption.Scat = 0;
   }
   l = strstr(chopt,"POL");  if (l) { Hoption.System = kPOLAR;       memcpy(l,"   ",3); }
   l = strstr(chopt,"CYL");  if (l) { Hoption.System = kCYLINDRICAL; memcpy(l,"   ",3); }
   l = strstr(chopt,"SPH");  if (l) { Hoption.System = kSPHERICAL;   memcpy(l,"   ",3); }
   l = strstr(chopt,"PSR");  if (l) { Hoption.System = kRAPIDITY;    memcpy(l,"   ",3); }

   l = strstr(chopt,"TRI");
   if (l) {
      Hoption.Scat = 0;
      Hoption.Color  = 0;
      Hoption.Tri = 1; memcpy(l,"   ",3);
      l = strstr(chopt,"FB");   if (l) { Hoption.FrontBox = 0; memcpy(l,"  ",2); }
      l = strstr(chopt,"BB");   if (l) { Hoption.BackBox = 0;  memcpy(l,"  ",2); }
      l = strstr(chopt,"ERR");  if (l) memcpy(l,"   ",3);
   }

   l = strstr(chopt,"AITOFF");
   if (l) {
      Hoption.Proj = 1; memcpy(l,"     ",6);       //Aitoff projection
   }
   l = strstr(chopt,"MERCATOR");
   if (l) {
      Hoption.Proj = 2; memcpy(l,"       ",8);     //Mercator projection
   }
   l = strstr(chopt,"SINUSOIDAL");
   if (l) {
      Hoption.Proj = 3; memcpy(l,"         ",10);  //Sinusoidal projection
   }
   l = strstr(chopt,"PARABOLIC");
   if (l) {
      Hoption.Proj = 4; memcpy(l,"        ",9);    //Parabolic projection
   }
   if (Hoption.Proj > 0) {
      Hoption.Scat = 0;
      Hoption.Contour = 14;
   }

   if (strstr(chopt,"A"))   Hoption.Axis = -1;
   if (strstr(chopt,"B"))   Hoption.Bar  = 1;
   if (strstr(chopt,"C") && !strstr(chopt,"CJUST")) { Hoption.Curve =1; Hoption.Hist = -1;}
   if (strstr(chopt,"F"))   Hoption.Fill =1;
   if (strstr(chopt,"][")) {Hoption.Off  =1; Hoption.Hist =1;}
   if (strstr(chopt,"F2"))  Hoption.Fill =2;
   if (strstr(chopt,"L")) { Hoption.Line =1; Hoption.Hist = -1;}
   if (strstr(chopt,"P")) { Hoption.Mark =1; Hoption.Hist = -1;}
   if (strstr(chopt,"Z"))   Hoption.Zscale =1;
   if (strstr(chopt,"*"))   Hoption.Star =1;
   if (strstr(chopt,"H"))   Hoption.Hist =2;
   if (strstr(chopt,"P0"))  Hoption.Mark =10;

   if (fH->InheritsFrom(TH2Poly::Class())) {
      if (Hoption.Fill+Hoption.Line+Hoption.Mark != 0 ) Hoption.Scat = 0;
   }

   if (strstr(chopt,"E")) {
      if (hdim == 1) {
         Hoption.Error = 1;
         if (strstr(chopt,"E1"))  Hoption.Error = 11;
         if (strstr(chopt,"E2"))  Hoption.Error = 12;
         if (strstr(chopt,"E3"))  Hoption.Error = 13;
         if (strstr(chopt,"E4"))  Hoption.Error = 14;
         if (strstr(chopt,"E5"))  Hoption.Error = 15;
         if (strstr(chopt,"E6"))  Hoption.Error = 16;
         if (strstr(chopt,"E0"))  Hoption.Error += 40;
         if (strstr(chopt,"X0")) {
            if (Hoption.Error == 1)  Hoption.Error += 20;
            Hoption.Error += 10;
         }
         if (Hoption.Text && fH->InheritsFrom(TProfile::Class())) {
            Hoption.Text += 2000;
            Hoption.Error = 0;
         }
      } else {
         if (Hoption.Error == 0) {
            Hoption.Error = 100;
            Hoption.Scat  = 0;
         }
         if (Hoption.Text) {
            Hoption.Text += 2000;
            Hoption.Error = 0;
         }
      }
   }

   if (Hoption.Surf == 15) {
      if (Hoption.System == kPOLAR || Hoption.System == kCARTESIAN) {
         Hoption.Surf = 13;
         Warning("MakeChopt","option SURF5 is not supported in Cartesian and Polar modes");
      }
   }

   //      Copy options from current style
   Hoption.Logx = gPad->GetLogx();
   Hoption.Logy = gPad->GetLogy();
   Hoption.Logz = gPad->GetLogz();

   //       Check options incompatibilities
   if (Hoption.Bar  == 1) Hoption.Hist = -1;
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Decode string `choptin` and fill Graphical cuts structure.

Int_t THistPainter::MakeCuts(char *choptin)
{

   fNcuts = 0;
   char *left = (char*)strchr(choptin,'[');
   if (!left) return 0;
   char *right = (char*)strchr(choptin,']');
   if (!right) return 0;
   Int_t nch = right-left;
   if (nch < 2) return 0;
   char *cuts = left+1;
   *right = 0;
   char *comma, *minus;
   Int_t i;
   while (1) {
      comma = strchr(cuts,',');
      if (comma) *comma = 0;
      minus = strchr(cuts,'-');
      if (minus) cuts = minus+1;
      while (*cuts == ' ') cuts++;
      Int_t nc = strlen(cuts);
      while (cuts[nc-1] == ' ') {cuts[nc-1] = 0; nc--;}
      TIter next(gROOT->GetListOfSpecials());
      TCutG *cut=0;
      TObject *obj;
      while ((obj = next())) {
         if (!obj->InheritsFrom(TCutG::Class())) continue;
         if (strcmp(obj->GetName(),cuts)) continue;
         cut = (TCutG*)obj;
         break;
      }
      if (cut) {
         fCuts[fNcuts] = cut;
         fCutsOpt[fNcuts] = 1;
         if (minus) fCutsOpt[fNcuts] = -1;
         fNcuts++;
      }
      if (!comma) break;
      cuts = comma+1;
   }
   for (i=0;i<=nch;i++) left[i] = ' ';
   return fNcuts;
}

////////////////////////////////////////////////////////////////////////////////
/// [Control routine to paint any kind of histograms](#HP00)

void THistPainter::Paint(Option_t *option)
{

   if (fH->GetBuffer()) fH->BufferEmpty(-1);

   //For iOS: put the histogram on the top of stack of pickable objects.
   const TPickerStackGuard topPush(fH);

   gPad->SetVertical(kTRUE);

   TH1 *oldhist = gCurrentHist;
   gCurrentHist = fH;
   TH1 *hsave   = fH;
   Double_t minsav = fH->GetMinimumStored();

   if (!MakeChopt(option)) return; //check options and fill Hoption structure

   // Paint using TSpectrum2Painter
   if (Hoption.Spec) {
      if (!TableInit()) return;
      if (!TClass::GetClass("TSpectrum2Painter")) gSystem->Load("libSpectrumPainter");
      gROOT->ProcessLineFast(Form("TSpectrum2Painter::PaintSpectrum((TH2F*)0x%lx,\"%s\",%d)",
                                  (ULong_t)fH, option, Hoption.Spec));
      return;
   }

   if (Hoption.Pie) {
      if (fH->GetDimension() == 1) {
         if (!fPie) fPie = new TPie(fH);
         fPie->Paint(option);
      } else {
         Error("Paint", "Option PIE is for 1D histograms only");
      }
      return;
   } else {
      if (fPie) delete fPie;
      fPie = 0;
   }

   fXbuf  = new Double_t[kNMAX];
   fYbuf  = new Double_t[kNMAX];
   if (fH->GetDimension() > 2) {
      PaintH3(option);
      fH->SetMinimum(minsav);
      if (Hoption.Func) {
         Hoption_t hoptsave = Hoption;
         Hparam_t  hparsave = Hparam;
         PaintFunction(option);
         SetHistogram(hsave);
         Hoption = hoptsave;
         Hparam  = hparsave;
      }
      gCurrentHist = oldhist;
      delete [] fXbuf; delete [] fYbuf;
      return;
   }
   TView *view = gPad->GetView();
   if (view) {
      if (!Hoption.Lego && !Hoption.Surf && !Hoption.Tri) {
         delete view;
         gPad->SetView(0);
      }
   }
   if (fH->GetDimension() > 1 || Hoption.Lego || Hoption.Surf) {
      // In case of 1D histogram, Z axis becomes Y axis.
      Int_t logysav=0, logzsav=0;
      if (fH->GetDimension() == 1) {
         logysav = Hoption.Logy;
         logzsav = Hoption.Logz;
         Hoption.Logz = 0;
         if (Hoption.Logy) {
            Hoption.Logz = 1;
            Hoption.Logy = 0;
         }
      }
      PaintTable(option);
      if (Hoption.Func) {
         Hoption_t hoptsave = Hoption;
         Hparam_t  hparsave = Hparam;
         PaintFunction(option);
         SetHistogram(hsave);
         Hoption = hoptsave;
         Hparam  = hparsave;
      }
      fH->SetMinimum(minsav);
      gCurrentHist = oldhist;
      delete [] fXbuf; delete [] fYbuf;
      if (fH->GetDimension() == 1) {
         Hoption.Logy = logysav;
         Hoption.Logz = logzsav;
      }
      return;
   }

   if (Hoption.Bar >= 20) {
      PaintBarH(option);
      delete [] fXbuf; delete [] fYbuf;
      return;
   }

   // fill Hparam structure with histo parameters
   if (!PaintInit()) {
      delete [] fXbuf; delete [] fYbuf;
      return;
   }

   //          Picture surround (if new page) and page number (if requested).
   //          Histogram surround (if not option "Same").
   PaintFrame();

   //          Paint histogram axis only
   Bool_t gridx = gPad->GetGridx();
   Bool_t gridy = gPad->GetGridy();
   if (Hoption.Axis > 0) {
      if (Hoption.Axis > 1) PaintAxis(kTRUE);  //axis with grid
      else {
         if (gridx) gPad->SetGridx(0);
         if (gridy) gPad->SetGridy(0);
         PaintAxis(kFALSE);
         if (gridx) gPad->SetGridx(1);
         if (gridy) gPad->SetGridy(1);
      }
      if ((Hoption.Same%10) ==1) Hoption.Same += 1;
      goto paintstat;
   }
   if (gridx || gridy) PaintAxis(kTRUE); //    Draw the grid only

   //          test for options BAR or HBAR
   if (Hoption.Bar >= 10) {
      PaintBar(option);
   }

   //          do not draw histogram if error bars required
   if (!Hoption.Error) {
      if (Hoption.Hist && Hoption.Bar<10) PaintHist(option);
   }

   //         test for error bars or option E
   if (Hoption.Error) {
      PaintErrors(option);
      if (Hoption.Hist == 2) PaintHist(option);
   }

   if (Hoption.Text) PaintText(option);

   //         test for associated function
   if (Hoption.Func) {
      Hoption_t hoptsave = Hoption;
      Hparam_t  hparsave = Hparam;
      PaintFunction(option);
      SetHistogram(hsave);
      Hoption = hoptsave;
      Hparam  = hparsave;
   }

   if (gridx) gPad->SetGridx(0);
   if (gridy) gPad->SetGridy(0);
   PaintAxis(kFALSE);
   if (gridx) gPad->SetGridx(1);
   if (gridy) gPad->SetGridy(1);

   PaintTitle();  // Draw histogram title

   // Draw box with histogram statistics and/or fit parameters
paintstat:
   if ((Hoption.Same%10) != 1 && !fH->TestBit(TH1::kNoStats)) {  // bit set via TH1::SetStats
      TIter next(fFunctions);
      TObject *obj = 0;
      while ((obj = next())) {
         if (obj->InheritsFrom(TF1::Class())) break;
         obj = 0;
      }

      //Stat is painted twice (first, it will be in canvas' list of primitives),
      //second, it will be here, this is not required on iOS.
      //Condition is ALWAYS true on a platform different from iOS.
      if (!gPad->PadInSelectionMode() && !gPad->PadInHighlightMode())
         PaintStat(gStyle->GetOptStat(),(TF1*)obj);
   }
   fH->SetMinimum(minsav);
   gCurrentHist = oldhist;
   delete [] fXbuf; fXbuf = 0;
   delete [] fYbuf; fYbuf = 0;

}

////////////////////////////////////////////////////////////////////////////////
/// [Control function to draw a table as an arrow plot](#HP12)

void THistPainter::PaintArrows(Option_t *)
{
   Double_t xk, xstep, yk, ystep;
   Double_t dx, dy, x1, x2, y1, y2, xc, yc, dxn, dyn;
   Int_t   ncx  = Hparam.xlast - Hparam.xfirst + 1;
   Int_t   ncy  = Hparam.ylast - Hparam.yfirst + 1;
   Double_t xrg = gPad->GetUxmin();
   Double_t yrg = gPad->GetUymin();
   Double_t xln = gPad->GetUxmax() - xrg;
   Double_t yln = gPad->GetUymax() - yrg;
   Double_t cx  = (xln/Double_t(ncx) -0.03)/2;
   Double_t cy  = (yln/Double_t(ncy) -0.03)/2;
   Double_t dn  = 1.E-30;

   auto arrow = new TArrow();
   arrow->SetAngle(30);
   arrow->SetFillStyle(1001);
   arrow->SetFillColor(fH->GetLineColor());
   arrow->SetLineColor(fH->GetLineColor());
   arrow->SetLineWidth(fH->GetLineWidth());

   // Initialize the levels on the Z axis
   Int_t ncolors=0, ndivz=0;
   Double_t scale=0.;
   if (Hoption.Arrow>1) {
      ncolors = gStyle->GetNumberOfColors();
      Int_t ndiv    = fH->GetContour();
      if (ndiv == 0 ) {
         ndiv = gStyle->GetNumberContours();
         fH->SetContour(ndiv);
      }
      ndivz  = TMath::Abs(ndiv);
      if (fH->TestBit(TH1::kUserContour) == 0) fH->SetContour(ndiv);
      scale = ndivz/(fH->GetMaximum()-fH->GetMinimum());
   }

   for (Int_t id=1;id<=2;id++) {
      for (Int_t j=Hparam.yfirst; j<=Hparam.ylast;j++) {
         yk    = fYaxis->GetBinLowEdge(j);
         ystep = fYaxis->GetBinWidth(j);
         for (Int_t i=Hparam.xfirst; i<=Hparam.xlast;i++) {
            xk    = fXaxis->GetBinLowEdge(i);
            xstep = fXaxis->GetBinWidth(i);
            if (!IsInside(xk+0.5*xstep,yk+0.5*ystep)) continue;
            if (i == Hparam.xfirst) {
               dx = fH->GetBinContent(i+1, j) - fH->GetBinContent(i, j);
            } else if (i == Hparam.xlast) {
               dx = fH->GetBinContent(i, j) - fH->GetBinContent(i-1, j);
            } else {
               dx = 0.5*(fH->GetBinContent(i+1, j) - fH->GetBinContent(i-1, j));
            }
            if (j == Hparam.yfirst) {
               dy = fH->GetBinContent(i, j+1) - fH->GetBinContent(i, j);
            } else if (j == Hparam.ylast) {
               dy = fH->GetBinContent(i, j) - fH->GetBinContent(i, j-1);
            } else {
               dy = 0.5*(fH->GetBinContent(i, j+1) - fH->GetBinContent(i, j-1));
            }
            if (id == 1) {
               dn = TMath::Max(dn, TMath::Abs(dx));
               dn = TMath::Max(dn, TMath::Abs(dy));
            } else if (id == 2) {
               xc  = xrg + xln*(Double_t(i - Hparam.xfirst+1)-0.5)/Double_t(ncx);
               dxn = cx*dx/dn;
               x1  = xc - dxn;
               x2  = xc + dxn;
               yc  = yrg + yln*(Double_t(j - Hparam.yfirst+1)-0.5)/Double_t(ncy);
               dyn = cy*dy/dn;
               y1  = yc - dyn;
               y2  = yc + dyn;
               if (Hoption.Arrow>1) {
                  int color = Int_t(0.01+(fH->GetBinContent(i, j)-fH->GetMinimum())*scale);
                  Int_t theColor = Int_t((color+0.99)*Float_t(ncolors)/Float_t(ndivz));
                  if (theColor > ncolors-1) theColor = ncolors-1;
                  arrow->SetFillColor(gStyle->GetColorPalette(theColor));
                  arrow->SetLineColor(gStyle->GetColorPalette(theColor));
               }
               if (TMath::Abs(x2-x1) > 0. || TMath::Abs(y2-y1) > 0.) {
                  arrow->PaintArrow(x1, y1, x2, y2, 0.015, "|>");
               } else {
                  arrow->PaintArrow(x1, y1, x2, y2, 0.005, "|>");
               }
            }
         }
      }
   }

   if (Hoption.Zscale) PaintPalette();
}

////////////////////////////////////////////////////////////////////////////////
/// Draw axis (2D case) of an histogram.
///
/// If `drawGridOnly` is `TRUE`, only the grid is painted (if needed). This allows
/// to draw the grid and the axis separately. In `THistPainter::Paint` this
/// feature is used to make sure that the grid is drawn in the background and
/// the axis tick marks in the foreground of the pad.

void THistPainter::PaintAxis(Bool_t drawGridOnly)
{

   //On iOS, grid should not be pickable and can not be highlighted.
   //Condition is never true on a platform different from iOS.
   if (drawGridOnly && (gPad->PadInHighlightMode() || gPad->PadInSelectionMode()))
      return;

   if (Hoption.Axis == -1) return;
   if (Hoption.Same && Hoption.Axis <= 0) return;

   // Repainting alphanumeric labels axis on a plot done with
   // the option HBAR (horizontal) needs some adjustments.
   TAxis *xaxis = 0;
   TAxis *yaxis = 0;
   if (Hoption.Same && Hoption.Axis) { // Axis repainted (TPad::RedrawAxis)
      if (fXaxis->GetLabels() || fYaxis->GetLabels()) { // One axis has alphanumeric labels
         TIter next(gPad->GetListOfPrimitives());
         TObject *obj;
         // Check if the first TH1 of THStack in the pad is drawn with the option HBAR
         while ((obj = next())) {
            if (!obj->InheritsFrom(TH1::Class()) &&
                !obj->InheritsFrom(THStack::Class())) continue;
            TString opt = obj->GetDrawOption();
            opt.ToLower();
            // if drawn with HBAR, the axis should be inverted and the pad set to horizontal
            if (strstr(opt,"hbar")) {
               gPad->SetVertical(kFALSE);
               xaxis = fXaxis;
               yaxis = fYaxis;
               if (!strcmp(xaxis->GetName(),"xaxis")) {
                  fXaxis = yaxis;
                  fYaxis = xaxis;
               }
            }
            break;
         }
      }
   }

   static char chopt[10] = "";
   Double_t gridl = 0;
   Int_t ndiv, ndivx, ndivy, nx1, nx2, ndivsave;
   Int_t useHparam = 0;
   Double_t umin, umax, uminsave, umaxsave;
   Short_t xAxisPos = Hoption.AxisPos/10;
   Short_t yAxisPos = Hoption.AxisPos - 10*xAxisPos;

   Double_t axmin = gPad->GetUxmin();
   Double_t axmax = gPad->GetUxmax();
   Double_t aymin = gPad->GetUymin();
   Double_t aymax = gPad->GetUymax();
   char *cw = 0;
   TGaxis axis;

   // In case of option 'cont4' or in case of option 'same' over a 'cont4 plot'
   // Hparam must be use for the axis limits.
   if (Hoption.Contour == 14) useHparam = 1;
   if (Hoption.Same) {
      TObject *obj;
      TIter next(gPad->GetListOfPrimitives());
      while ((obj=next())) {
         if (strstr(obj->GetDrawOption(),"cont4")) {
            useHparam = 1;
            break;
         }
      }
   }

   // Paint X axis

   //To make X-axis selectable on iOS device.
   if (gPad->PadInSelectionMode())
      gPad->PushSelectableObject(fXaxis);

   //This condition is ALWAYS true, unless it works on iOS (can be false on iOS).
   if (gPad->PadInSelectionMode() || !gPad->PadInHighlightMode() || (gPad->PadInHighlightMode() && gPad->GetSelected() == fXaxis)) {
      ndivx = fXaxis->GetNdivisions();
      if (ndivx > 1000) {
         nx2   = ndivx/100;
         nx1   = TMath::Max(1, ndivx%100);
         ndivx = 100*nx2 + Int_t(Float_t(nx1)*gPad->GetAbsWNDC());
      }
      axis.SetTextAngle(0);
      axis.ImportAxisAttributes(fXaxis);

      chopt[0] = 0;
      strlcat(chopt, "SDH",10);
      if (ndivx < 0) strlcat(chopt, "N",10);
      if (gPad->GetGridx()) {
         gridl = (aymax-aymin)/(gPad->GetY2() - gPad->GetY1());
         strlcat(chopt, "W",10);
      }

      // Define X-Axis limits
      if (Hoption.Logx) {
         strlcat(chopt, "G",10);
         ndiv = TMath::Abs(ndivx);
         if (useHparam) {
            umin = TMath::Power(10,Hparam.xmin);
            umax = TMath::Power(10,Hparam.xmax);
         } else {
            umin = TMath::Power(10,axmin);
            umax = TMath::Power(10,axmax);
         }
      } else {
         ndiv = TMath::Abs(ndivx);
         if (useHparam) {
            umin = Hparam.xmin;
            umax = Hparam.xmax;
         } else {
            umin = axmin;
            umax = axmax;
         }
      }

      // Display axis as time
      if (fXaxis->GetTimeDisplay()) {
         strlcat(chopt,"t",10);
         if (strlen(fXaxis->GetTimeFormatOnly()) == 0) {
            axis.SetTimeFormat(fXaxis->ChooseTimeFormat(Hparam.xmax-Hparam.xmin));
         }
      }

      // The main X axis can be on the bottom or on the top of the pad
      Double_t xAxisYPos1, xAxisYPos2;
      if (xAxisPos == 1) {
         // Main X axis top
         xAxisYPos1 = aymax;
         xAxisYPos2 = aymin;
      } else {
         // Main X axis bottom
         xAxisYPos1 = aymin;
         xAxisYPos2 = aymax;
      }

      // Paint the main X axis (always)
      uminsave = umin;
      umaxsave = umax;
      ndivsave = ndiv;
      axis.SetOption(chopt);
      if (xAxisPos) {
         strlcat(chopt, "-",10);
         gridl = -gridl;
      }
      if (Hoption.Same && Hoption.Axis) { // Axis repainted (TPad::RedrawAxis)
         axis.SetLabelSize(0.);
         axis.SetTitle("");
      }
      axis.PaintAxis(axmin, xAxisYPos1,
                     axmax, xAxisYPos1,
                     umin, umax,  ndiv, chopt, gridl, drawGridOnly);

      // Paint additional X axis (if needed)
      // On iOS, this additional X axis is neither pickable, nor highlighted.
      // Additional checks PadInSelectionMode etc. does not effect non-iOS platform.
      if (gPad->GetTickx() && !gPad->PadInSelectionMode() && !gPad->PadInHighlightMode()) {
         if (xAxisPos) {
            cw=strstr(chopt,"-");
            *cw='z';
         } else {
            strlcat(chopt, "-",10);
         }
         if (gPad->GetTickx() < 2) strlcat(chopt, "U",10);
         if ((cw=strstr(chopt,"W"))) *cw='z';
         axis.SetTitle("");
         axis.PaintAxis(axmin, xAxisYPos2,
                        axmax, xAxisYPos2,
                        uminsave, umaxsave,  ndivsave, chopt, gridl, drawGridOnly);
      }
   }//End of "if pad in selection mode etc".

   // Paint Y axis
   //On iOS, Y axis must pushed into the stack of selectable objects.
   if (gPad->PadInSelectionMode())
      gPad->PushSelectableObject(fYaxis);

   //This conditions is ALWAYS true on a platform, different from iOS (on iOS can be true, can be false).
   if (gPad->PadInSelectionMode() || !gPad->PadInHighlightMode() || (gPad->PadInHighlightMode() && gPad->GetSelected() == fYaxis)) {
      ndivy = fYaxis->GetNdivisions();
      axis.ImportAxisAttributes(fYaxis);

      chopt[0] = 0;
      strlcat(chopt, "SDH",10);
      if (ndivy < 0) strlcat(chopt, "N",10);
      if (gPad->GetGridy()) {
         gridl = (axmax-axmin)/(gPad->GetX2() - gPad->GetX1());
         strlcat(chopt, "W",10);
      }

      // Define Y-Axis limits
      if (Hoption.Logy) {
         strlcat(chopt, "G",10);
         ndiv = TMath::Abs(ndivy);
         if (useHparam) {
            umin = TMath::Power(10,Hparam.ymin);
            umax = TMath::Power(10,Hparam.ymax);
         } else {
            umin = TMath::Power(10,aymin);
            umax = TMath::Power(10,aymax);
         }
      } else {
         ndiv = TMath::Abs(ndivy);
         if (useHparam) {
            umin = Hparam.ymin;
            umax = Hparam.ymax;
         } else {
            umin = aymin;
            umax = aymax;
         }
      }

      // Display axis as time
      if (fYaxis->GetTimeDisplay()) {
         strlcat(chopt,"t",10);
         if (strlen(fYaxis->GetTimeFormatOnly()) == 0) {
            axis.SetTimeFormat(fYaxis->ChooseTimeFormat(Hparam.ymax-Hparam.ymin));
         }
      }

      // The main Y axis can be on the left or on the right of the pad
      Double_t yAxisXPos1, yAxisXPos2;
      if (yAxisPos == 1) {
         // Main Y axis left
         yAxisXPos1 = axmax;
         yAxisXPos2 = axmin;
      } else {
         // Main Y axis right
         yAxisXPos1 = axmin;
         yAxisXPos2 = axmax;
      }

      // Paint the main Y axis (always)
      uminsave = umin;
      umaxsave = umax;
      ndivsave = ndiv;
      axis.SetOption(chopt);
      if (yAxisPos) {
         strlcat(chopt, "+L",10);
         gridl = -gridl;
      }
      if (Hoption.Same && Hoption.Axis) { // Axis repainted (TPad::RedrawAxis)
         axis.SetLabelSize(0.);
         axis.SetTitle("");
      }
      axis.PaintAxis(yAxisXPos1, aymin,
                     yAxisXPos1, aymax,
                     umin, umax,  ndiv, chopt, gridl, drawGridOnly);

      // Paint the additional Y axis (if needed)
      // Additional checks for pad mode are required on iOS: this "second" axis is
      // neither pickable, nor highlighted. Additional checks have no effect on non-iOS platform.
      if (gPad->GetTicky() && !gPad->PadInSelectionMode() && !gPad->PadInHighlightMode()) {
         if (gPad->GetTicky() < 2) {
            strlcat(chopt, "U",10);
            axis.SetTickSize(-fYaxis->GetTickLength());
         } else {
            strlcat(chopt, "+L",10);
         }
         if ((cw=strstr(chopt,"W"))) *cw='z';
         axis.SetTitle("");
         axis.PaintAxis(yAxisXPos2, aymin,
                        yAxisXPos2, aymax,
                        uminsave, umaxsave,  ndivsave, chopt, gridl, drawGridOnly);
      }
   }//End of "if pad is in selection mode etc."

   // Reset the axis if they have been inverted in case of option HBAR
   if (xaxis) {
      fXaxis = xaxis;
      fYaxis = yaxis;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// [Draw a bar-chart in a normal pad.](#HP10)

void THistPainter::PaintBar(Option_t *)
{

   Int_t bar = Hoption.Bar - 10;
   Double_t xmin,xmax,ymin,ymax,umin,umax,w,y;
   Double_t offset = fH->GetBarOffset();
   Double_t width  = fH->GetBarWidth();
   TBox box;
   Int_t hcolor = fH->GetFillColor();
   if (hcolor == gPad->GetFrameFillColor()) ++hcolor;
   Int_t hstyle = fH->GetFillStyle();
   box.SetFillColor(hcolor);
   box.SetFillStyle(hstyle);
   box.SetLineStyle(fH->GetLineStyle());
   box.SetLineColor(fH->GetLineColor());
   box.SetLineWidth(fH->GetLineWidth());
   for (Int_t bin=fXaxis->GetFirst();bin<=fXaxis->GetLast();bin++) {
      y    = fH->GetBinContent(bin);
      xmin = gPad->XtoPad(fXaxis->GetBinLowEdge(bin));
      xmax = gPad->XtoPad(fXaxis->GetBinUpEdge(bin));
      ymin = gPad->GetUymin();
      ymax = gPad->YtoPad(y);
      if (ymax < gPad->GetUymin()) continue;
      if (ymax > gPad->GetUymax()) ymax = gPad->GetUymax();
      if (ymin < gPad->GetUymin()) ymin = gPad->GetUymin();
      if (Hoption.MinimumZero && ymin < 0)
         ymin=TMath::Min(0.,gPad->GetUymax());
      w    = (xmax-xmin)*width;
      xmin += offset*(xmax-xmin);
      xmax = xmin + w;
      if (bar < 1) {
         box.PaintBox(xmin,ymin,xmax,ymax);
      } else {
         umin = xmin + bar*(xmax-xmin)/10.;
         umax = xmax - bar*(xmax-xmin)/10.;
         //box.SetFillColor(hcolor+150); //bright
         box.SetFillColor(TColor::GetColorBright(hcolor)); //bright
         box.PaintBox(xmin,ymin,umin,ymax);
         box.SetFillColor(hcolor);
         box.PaintBox(umin,ymin,umax,ymax);
         box.SetFillColor(TColor::GetColorDark(hcolor)); //dark
         box.PaintBox(umax,ymin,xmax,ymax);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// [Draw a bar char in a rotated pad (X vertical, Y horizontal)](#HP10)

void THistPainter::PaintBarH(Option_t *)
{

   gPad->SetVertical(kFALSE);

   PaintInitH();

   TAxis *xaxis = fXaxis;
   TAxis *yaxis = fYaxis;
   if (!strcmp(xaxis->GetName(),"xaxis")) {
      fXaxis = yaxis;
      fYaxis = xaxis;
   }

   PaintFrame();
   PaintAxis(kFALSE);

   Int_t bar = Hoption.Bar - 20;
   Double_t xmin,xmax,ymin,ymax,umin,umax,w;
   Double_t offset = fH->GetBarOffset();
   Double_t width  = fH->GetBarWidth();
   TBox box;
   Int_t hcolor = fH->GetFillColor();
   if (hcolor == gPad->GetFrameFillColor()) ++hcolor;
   Int_t hstyle = fH->GetFillStyle();
   box.SetFillColor(hcolor);
   box.SetFillStyle(hstyle);
   box.SetLineStyle(fH->GetLineStyle());
   box.SetLineColor(fH->GetLineColor());
   box.SetLineWidth(fH->GetLineWidth());
   for (Int_t bin=fYaxis->GetFirst();bin<=fYaxis->GetLast();bin++) {
      ymin = gPad->YtoPad(fYaxis->GetBinLowEdge(bin));
      ymax = gPad->YtoPad(fYaxis->GetBinUpEdge(bin));
      xmin = gPad->GetUxmin();
      xmax = gPad->XtoPad(fH->GetBinContent(bin));
      if (xmax < gPad->GetUxmin()) continue;
      if (xmax > gPad->GetUxmax()) xmax = gPad->GetUxmax();
      if (xmin < gPad->GetUxmin()) xmin = gPad->GetUxmin();
      if (Hoption.MinimumZero && xmin < 0)
         xmin=TMath::Min(0.,gPad->GetUxmax());
      w    = (ymax-ymin)*width;
      ymin += offset*(ymax-ymin);
      ymax = ymin + w;
      if (bar < 1) {
         box.PaintBox(xmin,ymin,xmax,ymax);
      } else {
         umin = ymin + bar*(ymax-ymin)/10.;
         umax = ymax - bar*(ymax-ymin)/10.;
         box.SetFillColor(TColor::GetColorDark(hcolor)); //dark
         box.PaintBox(xmin,ymin,xmax,umin);
         box.SetFillColor(hcolor);
         box.PaintBox(xmin,umin,xmax,umax);
         box.SetFillColor(TColor::GetColorBright(hcolor)); //bright
         box.PaintBox(xmin,umax,xmax,ymax);
      }
   }

   PaintTitle();

   //    Draw box with histogram statistics and/or fit parameters
   if ((Hoption.Same%10) != 1 && !fH->TestBit(TH1::kNoStats)) {  // bit set via TH1::SetStats
      TIter next(fFunctions);
      TObject *obj = 0;
      while ((obj = next())) {
         if (obj->InheritsFrom(TF1::Class())) break;
         obj = 0;
      }
      PaintStat(gStyle->GetOptStat(),(TF1*)obj);
   }

   fXaxis = xaxis;
   fYaxis = yaxis;
}

////////////////////////////////////////////////////////////////////////////////
/// [Control function to draw a 2D histogram as a box plot](#HP13)

void THistPainter::PaintBoxes(Option_t *)
{

   Style_t fillsav   = fH->GetFillStyle();
   Style_t colsav    = fH->GetFillColor();
   if (fH->GetFillColor() == 0)  fH->SetFillStyle(0);
   if (Hoption.Box == 11) fH->SetFillStyle(1001);
   fH->TAttLine::Modify();
   fH->TAttFill::Modify();

   Double_t z, xk,xstep, yk, ystep, xcent, ycent, xlow, xup, ylow, yup;
   Double_t ux1 = gPad->PixeltoX(1);
   Double_t ux0 = gPad->PixeltoX(0);
   Double_t uy1 = gPad->PixeltoY(1);
   Double_t uy0 = gPad->PixeltoY(0);
   Double_t dxmin = 0.51*(gPad->PadtoX(ux1)-gPad->PadtoX(ux0));
   Double_t dymin = 0.51*(gPad->PadtoY(uy0)-gPad->PadtoY(uy1));

   Double_t zmin = TMath::Max(fH->GetMinimum(),0.);
   Double_t zmax = TMath::Max(TMath::Abs(fH->GetMaximum()),
                              TMath::Abs(fH->GetMinimum()));
   Double_t zminlin = zmin, zmaxlin = zmax;

   // In case of option SAME, zmin and zmax values are taken from the
   // first plotted 2D histogram.
   if (Hoption.Same > 0 && Hoption.Same < 10) {
      TH2 *h2;
      TIter next(gPad->GetListOfPrimitives());
      while ((h2 = (TH2 *)next())) {
         if (!h2->InheritsFrom(TH2::Class())) continue;
         zmin = TMath::Max(h2->GetMinimum(), 0.);
         zmax = TMath::Max(TMath::Abs(h2->GetMaximum()),
                           TMath::Abs(h2->GetMinimum()));
         zminlin = zmin;
         zmaxlin = zmax;
         if (Hoption.Logz) {
            if (zmin <= 0) {
               zmin = TMath::Log10(zmax*0.001);
            } else {
               zmin = TMath::Log10(zmin);
            }
            zmax = TMath::Log10(zmax);
         }
         break;
      }
   } else {
      if (Hoption.Logz) {
         if (zmin > 0) {
            zmin = TMath::Log10(zmin);
            zmax = TMath::Log10(zmax);
         } else {
            return;
         }
      }
   }

   Double_t zratio, dz = zmax - zmin;
   Bool_t kZminNeg     = kFALSE;
   if (fH->GetMinimum()<0) kZminNeg = kTRUE;
   Bool_t kZNeg        = kFALSE;

   // Define the dark and light colors the "button style" boxes.
   Color_t color = fH->GetFillColor();
   Color_t light=0, dark=0;
   if (Hoption.Box == 11) {
      light = TColor::GetColorBright(color);
      dark  = TColor::GetColorDark(color);
   }

   // Loop over all the bins and draw the boxes
   for (Int_t j=Hparam.yfirst; j<=Hparam.ylast;j++) {
      yk    = fYaxis->GetBinLowEdge(j);
      ystep = fYaxis->GetBinWidth(j);
      ycent = 0.5*ystep;
      for (Int_t i=Hparam.xfirst; i<=Hparam.xlast;i++) {
         Int_t bin  = j*(fXaxis->GetNbins()+2) + i;
         xk    = fXaxis->GetBinLowEdge(i);
         xstep = fXaxis->GetBinWidth(i);
         if (!IsInside(xk+0.5*xstep,yk+0.5*ystep)) continue;
         xcent = 0.5*xstep;
         z     = Hparam.factor*fH->GetBinContent(bin);
         kZNeg = kFALSE;

         if (TMath::Abs(z) <  zminlin) continue; // Can be the case with ...
         if (TMath::Abs(z) >  zmaxlin) z = zmaxlin; // ... option Same
         if (kZminNeg && z==0) continue;      // Do not draw empty bins if case of histo with netgative bins.

         if (z < 0) {
            if (Hoption.Logz) continue;
            z = -z;
            kZNeg = kTRUE;
         }
         if (Hoption.Logz) {
            if (z != 0) z = TMath::Log10(z);
            else        z = zmin;
         }

         if (dz == 0) continue;
         zratio = TMath::Sqrt((z-zmin)/dz);
         if (zratio == 0) continue;

         xup  = xcent*zratio + xk + xcent;
         xlow = 2*(xk + xcent) - xup;
         if (xup-xlow < dxmin) xup = xlow+dxmin;
         if (Hoption.Logx) {
            if (xup > 0)  xup  = TMath::Log10(xup);
            else continue;
            if (xlow > 0) xlow = TMath::Log10(xlow);
            else continue;
         }

         yup  = ycent*zratio + yk + ycent;
         ylow = 2*(yk + ycent) - yup;
         if (yup-ylow < dymin) yup = ylow+dymin;
         if (Hoption.Logy) {
            if (yup > 0)  yup  = TMath::Log10(yup);
            else continue;
            if (ylow > 0) ylow = TMath::Log10(ylow);
            else continue;
         }

         xlow = TMath::Max(xlow, gPad->GetUxmin());
         ylow = TMath::Max(ylow, gPad->GetUymin());
         xup  = TMath::Min(xup , gPad->GetUxmax());
         yup  = TMath::Min(yup , gPad->GetUymax());

         if (xlow >= xup) continue;
         if (ylow >= yup) continue;

         if (Hoption.Box == 1) {
            fH->SetFillColor(color);
            fH->TAttFill::Modify();
            gPad->PaintBox(xlow, ylow, xup, yup);
            if (kZNeg) {
               gPad->PaintLine(xlow, ylow, xup, yup);
               gPad->PaintLine(xlow, yup, xup, ylow);
            }
         } else if (Hoption.Box == 11) {
            // Draw the center of the box
            fH->SetFillColor(color);
            fH->TAttFill::Modify();
            gPad->PaintBox(xlow, ylow, xup, yup);

            // Draw top&left part of the box
            Double_t x[7], y[7];
            Double_t bwidth = 0.1;
            x[0] = xlow;                     y[0] = ylow;
            x[1] = xlow + bwidth*(xup-xlow); y[1] = ylow + bwidth*(yup-ylow);
            x[2] = x[1];                     y[2] = yup - bwidth*(yup-ylow);
            x[3] = xup - bwidth*(xup-xlow);  y[3] = y[2];
            x[4] = xup;                      y[4] = yup;
            x[5] = xlow;                     y[5] = yup;
            x[6] = xlow;                     y[6] = ylow;
            if (kZNeg) fH->SetFillColor(dark);
            else       fH->SetFillColor(light);
            fH->TAttFill::Modify();
            gPad->PaintFillArea(7, x, y);

            // Draw bottom&right part of the box
            x[0] = xlow;                     y[0] = ylow;
            x[1] = xlow + bwidth*(xup-xlow); y[1] = ylow + bwidth*(yup-ylow);
            x[2] = xup - bwidth*(xup-xlow);  y[2] = y[1];
            x[3] = x[2];                     y[3] = yup - bwidth*(yup-ylow);
            x[4] = xup;                      y[4] = yup;
            x[5] = xup;                      y[5] = ylow;
            x[6] = xlow;                     y[6] = ylow;
            if (kZNeg) fH->SetFillColor(light);
            else       fH->SetFillColor(dark);
            fH->TAttFill::Modify();
            gPad->PaintFillArea(7, x, y);
         }
      }
   }

   if (Hoption.Zscale) PaintPalette();
   fH->SetFillStyle(fillsav);
   fH->SetFillColor(colsav);
   fH->TAttFill::Modify();
}



////////////////////////////////////////////////////////////////////////////////
/// [Control function to draw a 2D histogram as a candle (box) plot or violin plot](#HP14)

void THistPainter::PaintCandlePlot(Option_t *)
{
   TH1D *hproj;
   TH2D *h2 = (TH2D*)fH;

   TCandle myCandle;
   myCandle.SetOption((TCandle::CandleOption)Hoption.Candle);
   myCandle.SetMarkerColor(fH->GetLineColor());
   myCandle.SetLineColor(fH->GetLineColor());
   myCandle.SetLineWidth(fH->GetLineWidth());
   myCandle.SetFillColor(fH->GetFillColor());
   myCandle.SetFillStyle(fH->GetFillStyle());
   myCandle.SetMarkerSize(fH->GetMarkerSize());
   myCandle.SetMarkerStyle(fH->GetMarkerStyle());
   myCandle.SetLog(Hoption.Logx,Hoption.Logy, Hoption.Logz);

   Bool_t swapXY = myCandle.IsHorizontal();
   const Double_t standardCandleWidth = 0.66;
   const Double_t standardHistoWidth = 0.8;

   double allMaxContent = h2->GetBinContent(h2->GetMaximumBin());
   double allMaxIntegral = 0;

   if (!swapXY) { // Vertical candle
      //Determining the slice with the maximum content
      for (Int_t i=Hparam.xfirst; i<=Hparam.xlast; i++) {
         hproj = h2->ProjectionY("_px", i, i);
         if (hproj->Integral() > allMaxIntegral) allMaxIntegral = hproj->Integral();
      }
      for (Int_t i=Hparam.xfirst; i<=Hparam.xlast; i++) {
         Double_t binPosX = fXaxis->GetBinLowEdge(i);
         Double_t binWidth = fXaxis->GetBinWidth(i);
         hproj = h2->ProjectionY("_px", i, i);
         if (hproj->GetEntries() !=0) {
            Double_t candleWidth = fH->GetBarWidth();
            Double_t offset = fH->GetBarOffset()*binWidth;
            double myMaxContent = hproj->GetBinContent(hproj->GetMaximumBin());
            double myIntegral = hproj->Integral();
            Double_t histoWidth = candleWidth;
            if (candleWidth > 0.999 && candleWidth < 1.001) {
                candleWidth = standardCandleWidth;
                histoWidth = standardHistoWidth;
            }
            if (Hoption.Logz && myMaxContent > 0) {
                histoWidth *= myMaxContent/TMath::Log10(myMaxContent);
                if (myCandle.IsViolinScaled() && myMaxContent > 0 && allMaxContent > 0) histoWidth *= TMath::Log10(myMaxContent)/TMath::Log10(allMaxContent);
            } else if (myCandle.IsViolinScaled()) histoWidth *= myMaxContent/allMaxContent;
            if (myCandle.IsCandleScaled()) candleWidth *= myIntegral/allMaxIntegral;

            myCandle.SetAxisPosition(binPosX+binWidth/2. + offset);
            myCandle.SetCandleWidth(candleWidth*binWidth);
            myCandle.SetHistoWidth(histoWidth*binWidth);
            myCandle.SetHistogram(hproj);
            myCandle.Paint();
         }
      }
   } else { // Horizontal candle
      //Determining the slice with the maximum content
      for (Int_t i=Hparam.yfirst; i<=Hparam.ylast; i++) {
         hproj = h2->ProjectionX("_py", i, i);
         if (hproj->Integral() > allMaxIntegral) allMaxIntegral = hproj->Integral();
      }
      for (Int_t i=Hparam.yfirst; i<=Hparam.ylast; i++) {
         Double_t binPosY = fYaxis->GetBinLowEdge(i);
         Double_t binWidth = fYaxis->GetBinWidth(i);
         hproj = h2->ProjectionX("_py", i, i);
         if (hproj->GetEntries() !=0) {
            Double_t candleWidth = fH->GetBarWidth();
            Double_t offset = fH->GetBarOffset()*binWidth;
            double myMaxContent = hproj->GetBinContent(hproj->GetMaximumBin());
            double myIntegral = hproj->Integral();
            Double_t histoWidth = candleWidth;
            if (candleWidth > 0.999 && candleWidth < 1.001) {
                candleWidth = standardCandleWidth;
                histoWidth = standardHistoWidth;
            }
            if (Hoption.Logz && myMaxContent > 0) {
                histoWidth *= myMaxContent/TMath::Log10(myMaxContent);
                if (myCandle.IsViolinScaled() && myMaxContent > 0 && allMaxContent > 0) histoWidth *= TMath::Log10(myMaxContent)/TMath::Log10(allMaxContent);
            } else if (myCandle.IsViolinScaled()) histoWidth *= myMaxContent/allMaxContent;
            if (myCandle.IsCandleScaled()) candleWidth *= myIntegral/allMaxIntegral;

            myCandle.SetAxisPosition(binPosY+binWidth/2. + offset);
            myCandle.SetCandleWidth(candleWidth*binWidth);
            myCandle.SetHistoWidth(histoWidth*binWidth);
            myCandle.SetHistogram(hproj);
            myCandle.Paint();
         }
      }
   }
}



////////////////////////////////////////////////////////////////////////////////
/// Returns the rendering regions for an axis to use in the COL2 option
///
/// The algorithm analyses the size of the axis compared to the size of
/// the rendering region. It figures out the boundaries to use for each color
/// of the rendering region. Only one axis is computed here.
///
/// This allows for a single computation of the boundaries before iterating
/// through all of the bins.
///
/// \param pAxis     the axis to consider
/// \param nPixels   the number of pixels to render axis into
/// \param isLog     whether the axis is log scale

std::vector<THistRenderingRegion>
THistPainter::ComputeRenderingRegions(TAxis* pAxis, Int_t nPixels, Bool_t isLog)
{
   std::vector<THistRenderingRegion> regions;

   enum STRATEGY { Bins, Pixels } strategy;

   Int_t nBins = (pAxis->GetLast() - pAxis->GetFirst() + 1);

   if (nBins >= nPixels) {
      // more bins than pixels... we should loop over pixels and sample
      strategy = Pixels;
   } else {
      // fewer bins than pixels... we should loop over bins
      strategy = Bins;
   }

   if (isLog) {

      Double_t xMin = pAxis->GetBinLowEdge(pAxis->GetFirst());
      Int_t binOffset=0;
      while (xMin <= 0 && ((pAxis->GetFirst()+binOffset) != pAxis->GetLast()) ) {
         binOffset++;
         xMin = pAxis->GetBinLowEdge(pAxis->GetFirst()+binOffset);
      }
      if (xMin <= 0) {
         // this should cause an error if we have
         return regions;
      }
      Double_t xMax = pAxis->GetBinUpEdge(pAxis->GetLast());

      if (strategy == Bins) {
         // logarithmic plot. we find the pixel for the bin
         // pixel = eta * log10(V) - alpha
         //  where eta = nPixels/(log10(Vmax)-log10(Vmin))
         //  and alpha = nPixels*log10(Vmin)/(log10(Vmax)-log10(Vmin))
         // and V is axis value
         Double_t eta = (nPixels-1.0)/(TMath::Log10(xMax) - TMath::Log10(xMin));
         Double_t offset = -1.0 * eta * TMath::Log10(xMin);

         for (Int_t bin=pAxis->GetFirst()+binOffset; bin<=pAxis->GetLast(); bin++) {

            // linear plot. we simply need to find the appropriate bin
            // for the
            Double_t xLowValue  = pAxis->GetBinLowEdge(bin);
            Double_t xUpValue   = pAxis->GetBinUpEdge(bin);
            Int_t xPx0          = eta*TMath::Log10(xLowValue)+ offset;
            Int_t xPx1          = eta*TMath::Log10(xUpValue) + offset;
            THistRenderingRegion region = {std::make_pair(xPx0, xPx1),
                                           std::make_pair(bin, bin+1)};
            regions.push_back(region);
         }

      } else {

         // loop over pixels

         Double_t beta = (TMath::Log10(xMax) - TMath::Log10(xMin))/(nPixels-1.0);

         for (Int_t pixelIndex=0; pixelIndex<(nPixels-1); pixelIndex++) {
            // linear plot
            Int_t binLow  = pAxis->FindBin(xMin*TMath::Power(10.0, beta*pixelIndex));
            Int_t binHigh = pAxis->FindBin(xMin*TMath::Power(10.0, beta*(pixelIndex+1)));
            THistRenderingRegion region = { std::make_pair(pixelIndex, pixelIndex+1),
               std::make_pair(binLow, binHigh)};
            regions.push_back(region);
         }
      }
   } else {
      // standard linear plot

      if (strategy == Bins) {
         // loop over bins
         for (Int_t bin=pAxis->GetFirst(); bin<=pAxis->GetLast(); bin++) {

            // linear plot. we simply need to find the appropriate bin
            // for the
            Int_t xPx0     = ((bin - pAxis->GetFirst()) * nPixels)/nBins;
            Int_t xPx1     = xPx0 + nPixels/nBins;

            // make sure we don't compute beyond our bounds
            if (xPx1>= nPixels) xPx1 = nPixels-1;

            THistRenderingRegion region = {std::make_pair(xPx0, xPx1),
               std::make_pair(bin, bin+1)};
            regions.push_back(region);
         }
      } else {
         // loop over pixels
         for (Int_t pixelIndex=0; pixelIndex<nPixels-1; pixelIndex++) {
            // linear plot
            Int_t binLow  = (nBins*pixelIndex)/nPixels + pAxis->GetFirst();
            Int_t binHigh = binLow + nBins/nPixels;
            THistRenderingRegion region = { std::make_pair(pixelIndex, pixelIndex+1),
                                            std::make_pair(binLow, binHigh)};
            regions.push_back(region);
         }
      }
   }

   return regions;
}

////////////////////////////////////////////////////////////////////////////////
/// [Rendering scheme for the COL2 and COLZ2 options] (#HP14)

void THistPainter::PaintColorLevelsFast(Option_t*)
{

   if (Hoption.System != kCARTESIAN) {
     Error("THistPainter::PaintColorLevelsFast(Option_t*)",
           "Only cartesian coordinates supported by 'COL2' option. Using 'COL' option instead.");
     PaintColorLevels(nullptr);
     return;
   }

   Double_t z;

   // Use existing max or min values. If either is already set
   // the appropriate value to use.
   Double_t zmin = fH->GetMinimumStored();
   Double_t zmax = fH->GetMaximumStored();
   Double_t originalZMin = zmin;
   Double_t originalZMax = zmax;
   if ((zmin == -1111) && (zmax == -1111)) {
      fH->GetMinimumAndMaximum(zmin, zmax);
      fH->SetMinimum(zmin);
      fH->SetMaximum(zmax);
   } else if (zmin == -1111) {
      zmin = fH->GetMinimum();
      fH->SetMinimum(zmin);
   } else if (zmax == -1111) {
      zmax = fH->GetMaximum();
      fH->SetMaximum(zmax);
   }

   Double_t dz = zmax - zmin;
   if (dz <= 0) { // Histogram filled with a constant value
      zmax += 0.1*TMath::Abs(zmax);
      zmin -= 0.1*TMath::Abs(zmin);
      dz = zmax - zmin;
   }

   if (Hoption.Logz) {
      if (zmin > 0) {
         zmin = TMath::Log10(zmin);
         zmax = TMath::Log10(zmax);
         dz = zmax - zmin;
      } else {
         Error("THistPainter::PaintColorLevelsFast(Option_t*)",
               "Cannot plot logz because bin content is less than 0.");
         return;
      }
   }

   // Initialize the levels on the Z axis
   Int_t ndiv   = fH->GetContour();
   if (ndiv == 0 ) {
      ndiv = gStyle->GetNumberContours();
      fH->SetContour(ndiv);
   }
   std::vector<Double_t> colorBounds(ndiv);
   std::vector<Double_t> contours(ndiv, 0);
   if (fH->TestBit(TH1::kUserContour) == 0) {
      fH->SetContour(ndiv);
   } else {
      fH->GetContour(contours.data());
   }

   Double_t step = 1.0/ndiv;
   for (Int_t i=0; i<ndiv; ++i) {
      colorBounds[i] = step*i;
   }

   auto pFrame = gPad->GetFrame();
   Int_t px0 = gPad->XtoPixel(pFrame->GetX1());
   Int_t px1 = gPad->XtoPixel(pFrame->GetX2());
   Int_t py0 = gPad->YtoPixel(pFrame->GetY1());
   Int_t py1 = gPad->YtoPixel(pFrame->GetY2());
   Int_t nXPixels = px1-px0;
   Int_t nYPixels = py0-py1; // y=0 is at the top of the screen

   std::vector<Double_t> buffer(nXPixels*nYPixels, 0);

   auto xRegions = ComputeRenderingRegions(fXaxis, nXPixels, Hoption.Logx);
   auto yRegions = ComputeRenderingRegions(fYaxis, nYPixels, Hoption.Logy);
   if (xRegions.size() == 0 || yRegions.size() == 0) {
      Error("THistPainter::PaintColorLevelFast(Option_t*)",
            "Encountered error while computing rendering regions.");
      return;
   }

   Bool_t minExists = kFALSE;
   Bool_t maxExists = kFALSE;
   Double_t minValue = 1.;
   Double_t maxValue = 0.;
   for (auto& yRegion : yRegions) {
     for (auto& xRegion : xRegions ) {

       const auto& xBinRange = xRegion.fBinRange;
       const auto& yBinRange = yRegion.fBinRange;

       // sample the range
       z = fH->GetBinContent(xBinRange.second-1, yBinRange.second-1);

       if (Hoption.Logz) {
         if (z > 0) z = TMath::Log10(z);
         else       z = zmin;
       }

       // obey the user's max and min values if they were set
       if (z > zmax) z = zmax;
       if (z < zmin) z = zmin;

       if (fH->TestBit(TH1::kUserContour) == 1) {
          // contours are absolute values
          auto index  = TMath::BinarySearch(contours.size(), contours.data(), z);
          z = colorBounds[index];
       } else {
          Int_t index = 0;
          if (dz != 0) {
             index = 0.001 + ((z - zmin)/dz)*ndiv;
          }

          if (index == static_cast<Int_t>(colorBounds.size())) {
             index--;
          }

          // Do a little bookkeeping to use later for getting libAfterImage to produce
          // the correct colors
          if (index == 0) {
             minExists = kTRUE;
          } else if (index == static_cast<Int_t>(colorBounds.size()-1)) {
             maxExists = kTRUE;
          }

          z = colorBounds[index];

          if (z < minValue) {
             minValue = z;
          }
          if (z > maxValue) {
             maxValue = z;
          }
       }

       // fill in the actual pixels
       const auto& xPixelRange = xRegion.fPixelRange;
       const auto& yPixelRange = yRegion.fPixelRange;
       for (Int_t xPx = xPixelRange.first; xPx <= xPixelRange.second; ++xPx) {
         for (Int_t yPx = yPixelRange.first; yPx <= yPixelRange.second; ++yPx) {
           Int_t pixel = yPx*nXPixels + xPx;
           buffer[pixel] = z;
         }
       }
     } // end px loop
   } // end py loop

   // This is a bit of a hack to ensure that we span the entire color range and
   // don't screw up the colors for a sparse histogram. No one will notice that I set a
   // single pixel on the edge of the image to a different color. This is even more
   // true because the chosen pixels will be covered by the axis.
   if (minValue != maxValue) {
      if ( !minExists) {
         buffer.front() = 0;
      }

      if ( !maxExists) {
         buffer[buffer.size()-nXPixels] = 0.95;
      }
   }

   // Generate the TImage
   TImagePalette* pPalette = TImagePalette::CreateCOLPalette(ndiv);
   TImage* pImage = TImage::Create();
   pImage->SetImageQuality(TAttImage::kImgBest);
   pImage->SetImage(buffer.data(), nXPixels, nYPixels, pPalette);
   delete pPalette;

   Window_t wid = static_cast<Window_t>(gVirtualX->GetWindowID(gPad->GetPixmapID()));
   pImage->PaintImage(wid, px0, py1, 0, 0, nXPixels, nYPixels);
   delete pImage;

   if (Hoption.Zscale) PaintPalette();

   // Reset the maximum and minimum values to their original values
   // when this function was called. If we don't do this, an initial
   // value of -1111 will be replaced with the true max or min values.
   fH->SetMinimum(originalZMin);
   fH->SetMaximum(originalZMax);
}

////////////////////////////////////////////////////////////////////////////////
/// [Control function to draw a 2D histogram as a color plot.](#HP14)

void THistPainter::PaintColorLevels(Option_t*)
{
   Double_t z, zc, xk, xstep, yk, ystep, xlow, xup, ylow, yup;

   Double_t zmin = fH->GetMinimum();
   Double_t zmax = fH->GetMaximum();

   Double_t dz = zmax - zmin;
   if (dz <= 0) { // Histogram filled with a constant value
      zmax += 0.1*TMath::Abs(zmax);
      zmin -= 0.1*TMath::Abs(zmin);
      dz = zmax - zmin;
   }

   // In case of option SAME, zmin and zmax values are taken from the
   // first plotted 2D histogram.
   if (Hoption.Same > 0 && Hoption.Same < 10) {
      TH2 *h2;
      TIter next(gPad->GetListOfPrimitives());
      while ((h2 = (TH2 *)next())) {
         if (!h2->InheritsFrom(TH2::Class())) continue;
         zmin = h2->GetMinimum();
         zmax = h2->GetMaximum();
         fH->SetMinimum(zmin);
         fH->SetMaximum(zmax);
         if (Hoption.Logz) {
            if (zmin <= 0) {
               zmin = TMath::Log10(zmax*0.001);
            } else {
               zmin = TMath::Log10(zmin);
            }
            zmax = TMath::Log10(zmax);
         }
         dz = zmax - zmin;
         break;
      }
   } else {
      if (Hoption.Logz) {
         if (zmin > 0) {
            zmin = TMath::Log10(zmin);
            zmax = TMath::Log10(zmax);
            dz   = zmax - zmin;
         } else {
            return;
         }
      }
   }

   Style_t fillsav   = fH->GetFillStyle();
   Style_t colsav    = fH->GetFillColor();
   fH->SetFillStyle(1001);
   fH->TAttFill::Modify();

   // Initialize the levels on the Z axis
   Int_t ncolors  = gStyle->GetNumberOfColors();
   Int_t ndiv   = fH->GetContour();
   if (ndiv == 0 ) {
      ndiv = gStyle->GetNumberContours();
      fH->SetContour(ndiv);
   }
   Int_t ndivz  = TMath::Abs(ndiv);
   if (fH->TestBit(TH1::kUserContour) == 0) fH->SetContour(ndiv);
   Double_t scale = (dz ? ndivz / dz : 1.0);

   Int_t color;
   TProfile2D* prof2d = dynamic_cast<TProfile2D*>(fH);
   for (Int_t j=Hparam.yfirst; j<=Hparam.ylast;j++) {
      yk    = fYaxis->GetBinLowEdge(j);
      ystep = fYaxis->GetBinWidth(j);
      for (Int_t i=Hparam.xfirst; i<=Hparam.xlast;i++) {
         Int_t bin  = j*(fXaxis->GetNbins()+2) + i;
         xk    = fXaxis->GetBinLowEdge(i);
         xstep = fXaxis->GetBinWidth(i);
         if (Hoption.System == kPOLAR && xk<0) xk= 2*TMath::Pi()+xk;
         if (!IsInside(xk+0.5*xstep,yk+0.5*ystep)) continue;
         z     = fH->GetBinContent(bin);
         // if fH is a profile histogram do not draw empty bins
         if (prof2d) {
            const Double_t binEntries = prof2d->GetBinEntries(bin);
            if (binEntries == 0)
               continue;
         } else {
            // don't draw the empty bins for non-profile histograms
            // with positive content
            if (z == 0) {
               if (zmin >= 0 || Hoption.Logz) continue;
               if (Hoption.Color == 2) continue;
            }
         }

         if (Hoption.Logz) {
            if (z > 0) z = TMath::Log10(z);
            else       z = zmin;
         }
         if (z < zmin && !Hoption.Zero) continue;
         xup  = xk + xstep;
         xlow = xk;
         if (Hoption.Logx) {
            if (xup > 0)  xup  = TMath::Log10(xup);
            else continue;
            if (xlow > 0) xlow = TMath::Log10(xlow);
            else continue;
         }
         yup  = yk + ystep;
         ylow = yk;
         if (Hoption.System != kPOLAR) {
            if (Hoption.Logy) {
               if (yup > 0)  yup  = TMath::Log10(yup);
               else continue;
               if (ylow > 0) ylow = TMath::Log10(ylow);
               else continue;
            }
            if (xup  < gPad->GetUxmin()) continue;
            if (yup  < gPad->GetUymin()) continue;
            if (xlow > gPad->GetUxmax()) continue;
            if (ylow > gPad->GetUymax()) continue;
            if (xlow < gPad->GetUxmin()) xlow = gPad->GetUxmin();
            if (ylow < gPad->GetUymin()) ylow = gPad->GetUymin();
            if (xup  > gPad->GetUxmax()) xup  = gPad->GetUxmax();
            if (yup  > gPad->GetUymax()) yup  = gPad->GetUymax();
         }

         if (fH->TestBit(TH1::kUserContour)) {
            zc = fH->GetContourLevelPad(0);
            if (z < zc) continue;
            color = -1;
            for (Int_t k=0; k<ndiv; k++) {
               zc = fH->GetContourLevelPad(k);
               if (z < zc) {
                  continue;
               } else {
                  color++;
               }
            }
         } else {
            color = Int_t(0.01+(z-zmin)*scale);
         }

         Int_t theColor = Int_t((color+0.99)*Float_t(ncolors)/Float_t(ndivz));
         if (theColor > ncolors-1) theColor = ncolors-1;
         fH->SetFillColor(gStyle->GetColorPalette(theColor));
         fH->TAttFill::Modify();
         if (Hoption.System != kPOLAR) {
            gPad->PaintBox(xlow, ylow, xup, yup);
         } else  {
            TCrown crown(0,0,ylow,yup,xlow*TMath::RadToDeg(),xup*TMath::RadToDeg());
            crown.SetFillColor(gStyle->GetColorPalette(theColor));
            crown.Paint();
         }
      }
   }

   if (Hoption.Zscale) PaintPalette();

   fH->SetFillStyle(fillsav);
   fH->SetFillColor(colsav);
   fH->TAttFill::Modify();

}

////////////////////////////////////////////////////////////////////////////////
/// [Control function to draw a 2D histogram as a contour plot.](#HP16)

void THistPainter::PaintContour(Option_t *option)
{

   Int_t i, j, count, ncontour, icol, n, lj, m, ix, jx, ljfill;
   Int_t itars, mode, ir[4];
   Double_t xsave, ysave, thesave,phisave,x[4], y[4], zc[4];

   if (Hoption.Contour == 14) {
      Hoption.Surf = 12;
      Hoption.Axis = 1;
      thesave = gPad->GetTheta();
      phisave = gPad->GetPhi();
      gPad->SetPhi(0.);
      gPad->SetTheta(90.);
      PaintSurface(option);
      gPad->SetPhi(phisave);
      gPad->SetTheta(thesave);
      TView *view = gPad->GetView();
      if (view) view->SetBit(kCannotRotate); //tested in ExecuteEvent
      PaintAxis();
      return;
   }

   if (Hoption.Same) {
      // If the contour is painted on a 3d plot, the contour lines are
      // paint in 3d too.
      TObject *obj;
      TIter next(gPad->GetListOfPrimitives());
      while ((obj=next())) {
         if (strstr(obj->GetDrawOption(),"surf") ||
             strstr(obj->GetDrawOption(),"lego") ||
             strstr(obj->GetDrawOption(),"tri")) {
               Hoption.Surf = 16;
               PaintSurface(option);
               return;
         }
      }
   }

   if (Hoption.Contour == 15) {
      TGraphDelaunay2D *dt = nullptr;
      TGraphDelaunay *dtOld = nullptr;
      TList *hl = fH->GetListOfFunctions();
      dt = (TGraphDelaunay2D*)hl->FindObject("TGraphDelaunay2D");
      if (!dt) dtOld = (TGraphDelaunay*)hl->FindObject("TGraphDelaunay");
      if (!dt && !dtOld) return;
      if (!fGraph2DPainter) {
         if (dt) fGraph2DPainter = new TGraph2DPainter(dt);
         else fGraph2DPainter = new TGraph2DPainter(dtOld);
      }
      fGraph2DPainter->Paint(option);
      return;
   }

   gPad->SetBit(TGraph::kClipFrame);

   Double_t *levels  = new Double_t[2*kMAXCONTOUR];
   Double_t *xarr    = new Double_t[2*kMAXCONTOUR];
   Double_t *yarr    = new Double_t[2*kMAXCONTOUR];
   Int_t  *itarr     = new Int_t[2*kMAXCONTOUR];

   Int_t npmax = 0;
   for (i=0;i<2*kMAXCONTOUR;i++) itarr[i] = 0;

   ncontour = fH->GetContour();
   if (ncontour == 0) {
      ncontour = gStyle->GetNumberContours();
      fH->SetContour(ncontour);
   }
   if (ncontour > kMAXCONTOUR) {
      Warning("PaintContour", "maximum number of contours is %d, asked for %d",
              kMAXCONTOUR, ncontour);
      ncontour = kMAXCONTOUR-1;
   }
   if (fH->TestBit(TH1::kUserContour) == 0) fH->SetContour(ncontour);

   for (i=0;i<ncontour;i++) levels[i] = fH->GetContourLevelPad(i);
   Int_t linesav   = fH->GetLineStyle();
   Int_t colorsav  = fH->GetLineColor();
   Int_t fillsav  = fH->GetFillColor();
   if (Hoption.Contour == 13) {
      fH->TAttLine::Modify();
   }

   TPolyLine **polys = 0;
   TPolyLine *poly=0;
   TObjArray *contours = 0;
   TList *list = 0;
   TGraph *graph = 0;
   Int_t *np = 0;
   if (Hoption.Contour == 1) {
      np = new Int_t[ncontour];
      for (i=0;i<ncontour;i++) np[i] = 0;
      polys = new TPolyLine*[ncontour];
      for (i=0;i<ncontour;i++) {
         polys[i] = new TPolyLine(100);
      }
      if (Hoption.List == 1) {
         contours = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
         if (contours) {
            gROOT->GetListOfSpecials()->Remove(contours);
            count = contours->GetSize();
            for (i=0;i<count;i++) {
               list = (TList*)contours->At(i);
               if (list) list->Delete();
            }
         }
         contours = new TObjArray(ncontour);
         contours->SetName("contours");
         gROOT->GetListOfSpecials()->Add(contours);
         for (i=0;i<ncontour;i++) {
            list = new TList();
            contours->Add(list);
         }
      }
   }
   Int_t theColor;
   Int_t ncolors = gStyle->GetNumberOfColors();
   Int_t ndivz   = TMath::Abs(ncontour);

   Int_t k,ipoly;
   for (j=Hparam.yfirst; j<Hparam.ylast; j++) {
      y[0] = fYaxis->GetBinCenter(j);
      y[1] = y[0];
      y[2] = fYaxis->GetBinCenter(j+1);
      y[3] = y[2];
      for (i=Hparam.xfirst; i<Hparam.xlast; i++) {
         zc[0] = fH->GetBinContent(i,   j);
         zc[1] = fH->GetBinContent(i+1, j);
         zc[2] = fH->GetBinContent(i+1, j+1);
         zc[3] = fH->GetBinContent(i,   j+1);
         if (!IsInside(fXaxis->GetBinCenter(i),fYaxis->GetBinCenter(j))) continue;
         if (Hoption.Logz) {
            if (zc[0] > 0)   zc[0] = TMath::Log10(zc[0]);
            else             zc[0] = Hparam.zmin;
            if (zc[1] > 0)   zc[1] = TMath::Log10(zc[1]);
            else             zc[1] = Hparam.zmin;
            if (zc[2] > 0)   zc[2] = TMath::Log10(zc[2]);
            else             zc[2] = Hparam.zmin;
            if (zc[3] > 0)   zc[3] = TMath::Log10(zc[3]);
            else             zc[3] = Hparam.zmin;
         }
         for (k=0;k<4;k++) {
            ir[k] = TMath::BinarySearch(ncontour,levels,zc[k]);
         }
         if (ir[0] != ir[1] || ir[1] != ir[2] || ir[2] != ir[3] || ir[3] != ir[0]) {
            x[0] = fXaxis->GetBinCenter(i);
            x[3] = x[0];
            x[1] = fXaxis->GetBinCenter(i+1);
            x[2] = x[1];
            if (zc[0] <= zc[1]) n = 0; else n = 1;
            if (zc[2] <= zc[3]) m = 2; else m = 3;
            if (zc[n] > zc[m]) n = m;
            n++;
            lj=1;
            for (ix=1;ix<=4;ix++) {
               m = n%4 + 1;
               ljfill = PaintContourLine(zc[n-1],ir[n-1],x[n-1],y[n-1],zc[m-1],
                     ir[m-1],x[m-1],y[m-1],&xarr[lj-1],&yarr[lj-1],&itarr[lj-1], levels);
               lj += 2*ljfill;
               n = m;
            }

            if (zc[0] <= zc[1]) n = 0; else n = 1;
            if (zc[2] <= zc[3]) m = 2; else m = 3;
            if (zc[n] > zc[m]) n = m;
            n++;
            lj=2;
            for (ix=1;ix<=4;ix++) {
               if (n == 1) m = 4;
               else        m = n-1;
               ljfill = PaintContourLine(zc[n-1],ir[n-1],x[n-1],y[n-1],zc[m-1],
                     ir[m-1],x[m-1],y[m-1],&xarr[lj-1],&yarr[lj-1],&itarr[lj-1], levels);
               lj += 2*ljfill;
               n = m;
            }

   //     Re-order endpoints

            count = 0;
            for (ix=1; ix<=lj-5; ix +=2) {
               //count = 0;
               while (itarr[ix-1] != itarr[ix]) {
                  xsave = xarr[ix];
                  ysave = yarr[ix];
                  itars = itarr[ix];
                  for (jx=ix; jx<=lj-5; jx +=2) {
                     xarr[jx]  = xarr[jx+2];
                     yarr[jx]  = yarr[jx+2];
                     itarr[jx] = itarr[jx+2];
                  }
                  xarr[lj-3]  = xsave;
                  yarr[lj-3]  = ysave;
                  itarr[lj-3] = itars;
                  if (count > 100) break;
                  count++;
               }
            }

            if (count > 100) continue;
            for (ix=1; ix<=lj-2; ix +=2) {
               theColor = Int_t((itarr[ix-1]+0.99)*Float_t(ncolors)/Float_t(ndivz));
               icol = gStyle->GetColorPalette(theColor);
               if (Hoption.Contour == 11) {
                  fH->SetLineColor(icol);
               }
               if (Hoption.Contour == 12) {
                  mode = icol%5;
                  if (mode == 0) mode = 5;
                  fH->SetLineStyle(mode);
               }
               if (Hoption.Contour != 1) {
                  fH->TAttLine::Modify();
                  gPad->PaintPolyLine(2,&xarr[ix-1],&yarr[ix-1]);
                  continue;
               }

               ipoly = itarr[ix-1];
               if (ipoly >=0 && ipoly <ncontour) {
                  poly = polys[ipoly];
                  poly->SetPoint(np[ipoly]  ,xarr[ix-1],yarr[ix-1]);
                  poly->SetPoint(np[ipoly]+1,xarr[ix],  yarr[ix]);
                  np[ipoly] += 2;
                  if (npmax < np[ipoly]) npmax = np[ipoly];
               }
            }
         } // end of if (ir[0]
      } //end of for (i
   } //end of for (j

   Double_t xmin,ymin;
   Double_t *xp, *yp;
   Int_t nadd,iminus,iplus;
   Double_t *xx, *yy;
   Int_t istart;
   Int_t first = ncontour;
   Int_t *polysort = 0;
   Int_t contListNb;
   if (Hoption.Contour != 1) goto theEND;

   //The 2 points line generated above are now sorted/merged to generate
   //a list of consecutive points.
   // If the option "List" has been specified, the list of points is saved
   // in the form of TGraph objects in the ROOT list of special objects.
   xmin = gPad->GetUxmin();
   ymin = gPad->GetUymin();
   xp = new Double_t[2*npmax];
   yp = new Double_t[2*npmax];
   polysort = new Int_t[ncontour];
   //find first positive contour
   for (ipoly=0;ipoly<ncontour;ipoly++) {
      if (levels[ipoly] >= 0) {first = ipoly; break;}
   }
   //store negative contours from 0 to minimum, then all positive contours
   k = 0;
   for (ipoly=first-1;ipoly>=0;ipoly--) {polysort[k] = ipoly; k++;}
   for (ipoly=first;ipoly<ncontour;ipoly++) {polysort[k] = ipoly; k++;}
   // we can now draw sorted contours
   contListNb = 0;
   fH->SetFillStyle(1001);
   for (k=0;k<ncontour;k++) {
      ipoly = polysort[k];
      if (np[ipoly] == 0) continue;
      if (Hoption.List) list = (TList*)contours->At(contListNb);
      contListNb++;
      poly = polys[ipoly];
      xx = poly->GetX();
      yy = poly->GetY();
      istart = 0;
      while (1) {
         iminus = npmax;
         iplus  = iminus+1;
         xp[iminus]= xx[istart];   yp[iminus] = yy[istart];
         xp[iplus] = xx[istart+1]; yp[iplus]  = yy[istart+1];
         xx[istart]   = xmin; yy[istart]   = ymin;
         xx[istart+1] = xmin; yy[istart+1] = ymin;
         while (1) {
            nadd = 0;
            for (i=2;i<np[ipoly];i+=2) {
               if ((iplus < 2*npmax-1) && (xx[i] == xp[iplus]) && (yy[i] == yp[iplus])) {
                  iplus++;
                  xp[iplus] = xx[i+1]; yp[iplus]  = yy[i+1];
                  xx[i]   = xmin; yy[i]   = ymin;
                  xx[i+1] = xmin; yy[i+1] = ymin;
                  nadd++;
               }
               if ((iminus > 0) && (xx[i+1] == xp[iminus]) && (yy[i+1] == yp[iminus])) {
                  iminus--;
                  xp[iminus] = xx[i];   yp[iminus]  = yy[i];
                  xx[i]   = xmin; yy[i]   = ymin;
                  xx[i+1] = xmin; yy[i+1] = ymin;
                  nadd++;
               }
            }
            if (nadd == 0) break;
         }
         theColor = Int_t((ipoly+0.99)*Float_t(ncolors)/Float_t(ndivz));
         icol = gStyle->GetColorPalette(theColor);
         if (ndivz > 1) fH->SetFillColor(icol);
         fH->TAttFill::Modify();
         gPad->PaintFillArea(iplus-iminus+1,&xp[iminus],&yp[iminus]);
         if (Hoption.List) {
            graph = new TGraph(iplus-iminus+1,&xp[iminus],&yp[iminus]);
            graph->SetFillColor(icol);
            graph->SetLineWidth(fH->GetLineWidth());
            list->Add(graph);
         }
         //check if more points are left
         istart = 0;
         for (i=2;i<np[ipoly];i+=2) {
            if (xx[i] != xmin && yy[i] != ymin) {
               istart = i;
               break;
            }
         }
         if (istart == 0) break;
      }
   }

   for (i=0;i<ncontour;i++) delete polys[i];
   delete [] polys;
   delete [] xp;
   delete [] yp;
   delete [] polysort;

theEND:
   gPad->ResetBit(TGraph::kClipFrame);
   if (Hoption.Zscale) PaintPalette();
   fH->SetLineStyle(linesav);
   fH->SetLineColor(colorsav);
   fH->SetFillColor(fillsav);
   if (np) delete [] np;
   delete [] xarr;
   delete [] yarr;
   delete [] itarr;
   delete [] levels;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill the matrix `xarr` and `yarr` for Contour Plot.

Int_t THistPainter::PaintContourLine(Double_t elev1, Int_t icont1, Double_t x1, Double_t y1,
                            Double_t elev2, Int_t icont2, Double_t x2, Double_t y2,
                            Double_t *xarr, Double_t *yarr, Int_t *itarr, Double_t *levels)
{

   Bool_t vert;
   Double_t tlen, tdif, elev, diff, pdif, xlen;
   Int_t n, i, icount;

   if (x1 == x2) {
      vert = kTRUE;
      tlen = y2 - y1;
   } else {
      vert = kFALSE;
      tlen = x2 - x1;
   }

   n = icont1 +1;
   tdif = elev2 - elev1;
   i = 0;
   icount = 0;
   while (n <= icont2 && i <= kMAXCONTOUR/2 -3) {
      //elev = fH->GetContourLevel(n);
      elev = levels[n];
      diff = elev - elev1;
      pdif = diff/tdif;
      xlen = tlen*pdif;
      if (vert) {
         if (Hoption.Logx)
            xarr[i] = TMath::Log10(x1);
         else
            xarr[i] = x1;
         if (Hoption.Logy)
            yarr[i] = TMath::Log10(y1 + xlen);
         else
            yarr[i] = y1 + xlen;
      } else {
         if (Hoption.Logx)
            xarr[i] = TMath::Log10(x1 + xlen);
         else
            xarr[i] = x1 + xlen;
         if (Hoption.Logy)
            yarr[i] = TMath::Log10(y1);
         else
            yarr[i] = y1;
      }
      itarr[i] = n;
      icount++;
      i +=2;
      n++;
   }
   return icount;
}

////////////////////////////////////////////////////////////////////////////////
/// [Draw 1D histograms error bars.](#HP09)

void THistPainter::PaintErrors(Option_t *)
{

   // On iOS, we do not highlight histogram, if it's not picked at the moment
   // (but part of histogram (axis or pavestat) was picked, that's why this code
   // is called at all. This conditional statement never executes on non-iOS platform.
   if (gPad->PadInHighlightMode() && gPad->GetSelected() != fH) return;

   const Int_t kBASEMARKER=8;
   Double_t xp, yp, ex1, ex2, ey1, ey2;
   Double_t delta;
   Double_t s2x, s2y, bxsize, bysize, symbolsize, xerror, sbase;
   Double_t xi1, xi2, xi3, xi4, yi1, yi2, yi3, yi4;
   Double_t xmin, xmax, ymin, ymax;
   Double_t logxmin = 0;
   Double_t logymin = 0;
   Double_t offset = 0.;
   Double_t width  = 0.;
   Int_t i, k, npoints, first, last, fixbin;
   Int_t if1 = 0;
   Int_t if2 = 0;
   Int_t drawmarker, errormarker;
   Int_t option0, option1, option2, option3, option4, optionE, optionEX0, optionI0;

   Double_t *xline = 0;
   Double_t *yline = 0;
   option0 = option1 = option2 = option3 = option4 = optionE = optionEX0 = optionI0 = 0;
   if (Hoption.Error >= 40) {Hoption.Error -=40; option0 = 1;}
   if (Int_t(Hoption.Error/10) == 2) {optionEX0 = 1; Hoption.Error -= 10;}
   if (Hoption.Error == 31) {optionEX0 = 1; Hoption.Error = 1;}
   if (Hoption.Error == 11) option1 = 1;
   if (Hoption.Error == 12) option2 = 1;
   if (Hoption.Error == 13) option3 = 1;
   if (Hoption.Error == 14) {option4 = 1; option3 = 1;}
   if (Hoption.Error == 15) {optionI0 = 1; option3 = 1;}
   if (Hoption.Error == 16) {optionI0 = 1; option4 = 1; option3 = 1;}
   if (option2+option3 == 0) optionE = 1;
   if (Hoption.Error == 0) optionE = 0;
   if (fXaxis->GetXbins()->fN) fixbin = 0;
   else                        fixbin = 1;

   offset = fH->GetBarOffset();
   width = fH->GetBarWidth();

   errormarker = fH->GetMarkerStyle();
   if (optionEX0) {
      xerror = 0;
   } else {
      xerror = gStyle->GetErrorX();
   }
   symbolsize  = fH->GetMarkerSize();
   if (errormarker == 1) symbolsize = 0.01;
   sbase       = symbolsize*kBASEMARKER;
   // set the graphics attributes

   fH->TAttLine::Modify();
   fH->TAttFill::Modify();
   fH->TAttMarker::Modify();

   // set the first and last bin

   Double_t factor = Hparam.factor;
   first      = Hparam.xfirst;
   last       = Hparam.xlast;
   npoints    = last - first  +1;
   xmin       = gPad->GetUxmin();
   xmax       = gPad->GetUxmax();
   ymin       = gPad->GetUymin();
   ymax       = gPad->GetUymax();


   if (option3) {
      xline = new Double_t[2*npoints];
      yline = new Double_t[2*npoints];
      if (!xline || !yline) {
         Error("PaintErrors", "too many points, out of memory");
         return;
      }
      if1 = 1;
      if2 = 2*npoints;
   }

   //  compute the offset of the error bars due to the symbol size
   s2x    = gPad->PixeltoX(Int_t(0.5*sbase)) - gPad->PixeltoX(0);
   s2y    =-gPad->PixeltoY(Int_t(0.5*sbase)) + gPad->PixeltoY(0);

   // compute size of the lines at the end of the error bars
   Int_t dxend = Int_t(gStyle->GetEndErrorSize());
   bxsize    = gPad->PixeltoX(dxend) - gPad->PixeltoX(0);
   bysize    =-gPad->PixeltoY(dxend) + gPad->PixeltoY(0);


   if (fixbin) {
      if (Hoption.Logx) xp = TMath::Power(10,Hparam.xmin) + 0.5*Hparam.xbinsize;
      else              xp = Hparam.xmin + 0.5*Hparam.xbinsize;
   } else {
      delta = fH->GetBinWidth(first);
      xp    = fH->GetBinLowEdge(first) + 0.5*delta;
   }

   // if errormarker = 0 or symbolsize = 0. no symbol is drawn
   if (Hoption.Logx) logxmin = TMath::Power(10,Hparam.xmin);
   if (Hoption.Logy) logymin = TMath::Power(10,Hparam.ymin);

   //    ---------------------- Loop over the points---------------------
   for (k=first; k<=last; k++) {

      //          get the data
      //     xp    = X position of the current point
      //     yp    = Y position of the current point
      //     ex1   = Low X error
      //     ex2   = Up X error
      //     ey1   = Low Y error
      //     ey2   = Up Y error
      //     (xi,yi) = Error bars coordinates

      // apply offset on errors for bar histograms
      Double_t xminTmp = gPad->XtoPad(fXaxis->GetBinLowEdge(k));
      Double_t xmaxTmp = gPad->XtoPad(fXaxis->GetBinUpEdge(k));
      if (Hoption.Logx) {
        xminTmp = TMath::Power(10, xminTmp);
        xmaxTmp = TMath::Power(10, xmaxTmp);
      }
      Double_t w    = (xmaxTmp-xminTmp)*width;
      xminTmp += offset*(xmaxTmp-xminTmp);
      xmaxTmp = xminTmp + w;
      xp = (xminTmp+xmaxTmp)/2.;

      if (Hoption.Logx) {
         if (xp <= 0) goto L30;
         if (xp < logxmin) goto L30;
         if (xp > TMath::Power(10,xmax)) break;
      } else {
         if (xp < xmin) goto L30;
         if (xp > xmax) break;
      }
      yp = factor*fH->GetBinContent(k);
      if (optionI0 && yp==0) goto L30;
      if (fixbin) {
         ex1 = xerror*Hparam.xbinsize;
      } else {
         delta = fH->GetBinWidth(k);
         ex1 = xerror*delta;
      }
      if (fH->GetBinErrorOption() == TH1::kNormal) {
         ey1 = factor*fH->GetBinError(k);
         ey2 = ey1;
      } else {
         ey1 = factor*fH->GetBinErrorLow(k);
         ey2 = factor*fH->GetBinErrorUp(k);
      }
      ex2 = ex1;

      xi4 = xp;
      xi3 = xp;
      xi2 = xp + ex2;
      xi1 = xp - ex1;

      yi1 = yp;
      yi2 = yp;
      yi3 = yp - ey1;
      yi4 = yp + ey2;

      //          take the LOG if necessary
      if (Hoption.Logx) {
         xi1 = TMath::Log10(TMath::Max(xi1,logxmin));
         xi2 = TMath::Log10(TMath::Max(xi2,logxmin));
         xi3 = TMath::Log10(TMath::Max(xi3,logxmin));
         xi4 = TMath::Log10(TMath::Max(xi4,logxmin));
      }
      if (Hoption.Logy) {
         yi1 = TMath::Log10(TMath::Max(yi1,logymin));
         yi2 = TMath::Log10(TMath::Max(yi2,logymin));
         yi3 = TMath::Log10(TMath::Max(yi3,logymin));
         yi4 = TMath::Log10(TMath::Max(yi4,logymin));
      }

      // test if error bars are not outside the limits
      //  otherwise they are truncated

      xi1 = TMath::Max(xi1,xmin);
      xi2 = TMath::Min(xi2,xmax);
      yi3 = TMath::Max(yi3,ymin);
      yi4 = TMath::Min(yi4,ymax);

      //  test if the marker is on the frame limits. If "Yes", the
      //  marker will not be drawn and the error bars will be readjusted.

      drawmarker = kTRUE;
      if (!option0 && !option3) {
         if (Hoption.Logy && yp < logymin) goto L30;
         if (yi1 < ymin || yi1 > ymax) goto L30;
         if (Hoption.Error != 0 && yp == 0 && ey1 <= 0) drawmarker = kFALSE;
      }
      if (!symbolsize || !errormarker) drawmarker = kFALSE;

      //  draw the error rectangles
      if (option2) {
         if (yi3 >= ymax) goto L30;
         if (yi4 <= ymin) goto L30;
         gPad->PaintBox(xi1,yi3,xi2,yi4);
      }

      //  keep points for fill area drawing
      if (option3) {
         xline[if1-1] = xi3;
         xline[if2-1] = xi3;
         yline[if1-1] = yi4;
         yline[if2-1] = yi3;
         if1++;
         if2--;
      }

      //          draw the error bars
      if (Hoption.Logy && yp < logymin) drawmarker = kFALSE;
      if (optionE && drawmarker) {
         if ((yi3 < yi1 - s2y) && (yi3 < ymax)) gPad->PaintLine(xi3,yi3,xi4,TMath::Min(yi1 - s2y,ymax));
         if ((yi1 + s2y < yi4) && (yi4 > ymin)) gPad->PaintLine(xi3,TMath::Max(yi1 + s2y, ymin),xi4,yi4);
         // don't duplicate the horizontal line
         if (Hoption.Hist != 2) {
            if (yi1<ymax && yi1>ymin) {
              if (xi1 < xi3 - s2x) gPad->PaintLine(xi1,yi1,xi3 - s2x,yi2);
              if (xi3 + s2x < xi2) gPad->PaintLine(xi3 + s2x,yi1,xi2,yi2);
            }
         }
      }
      if (optionE && !drawmarker && (ey1 != 0 || ey2 !=0)) {
         if ((yi3 < yi1) && (yi3 < ymax)) gPad->PaintLine(xi3,yi3,xi4,TMath::Min(yi1,ymax));
         if ((yi1 < yi4) && (yi4 > ymin)) gPad->PaintLine(xi3,TMath::Max(yi1,ymin),xi4,yi4);
         // don't duplicate the horizontal line
         if (Hoption.Hist != 2) {
            if (yi1<ymax && yi1>ymin) {
               if (xi1 < xi3) gPad->PaintLine(xi1,yi1,xi3,yi2);
               if (xi3 < xi2) gPad->PaintLine(xi3,yi1,xi2,yi2);
            }
         }
      }

      //          draw line at the end of the error bars

      if (option1 && drawmarker) {

         if (yi3 < yi1-s2y && yi3 < ymax && yi3 > ymin) gPad->PaintLine(xi3 - bxsize, yi3         , xi3 + bxsize, yi3);
         if (yi4 > yi1+s2y && yi4 < ymax && yi4 > ymin) gPad->PaintLine(xi3 - bxsize, yi4         , xi3 + bxsize, yi4);
         if (yi1 <= ymax && yi1 >= ymin) {
            if (xi1 < xi3-s2x) gPad->PaintLine(xi1         , yi1 - bysize, xi1         , yi1 + bysize);
            if (xi2 > xi3+s2x) gPad->PaintLine(xi2         , yi1 - bysize, xi2         , yi1 + bysize);
         }
      }

      //          draw the marker

      if (drawmarker) gPad->PaintPolyMarker(1, &xi3, &yi1);

L30:
      if (fixbin) xp += Hparam.xbinsize;
      else {
         if (k < last) {
            delta = fH->GetBinWidth(k+1);
            xp    = fH->GetBinLowEdge(k+1) + 0.5*delta;
         }
      }
   }  //end of for loop

   //          draw the filled area

   if (option3) {
      TGraph graph;
      graph.SetLineStyle(fH->GetLineStyle());
      graph.SetLineColor(fH->GetLineColor());
      graph.SetLineWidth(fH->GetLineWidth());
      graph.SetFillStyle(fH->GetFillStyle());
      graph.SetFillColor(fH->GetFillColor());
      Int_t logx = gPad->GetLogx();
      Int_t logy = gPad->GetLogy();
      gPad->SetLogx(0);
      gPad->SetLogy(0);

      // In some cases the number of points in the fill area is smaller than
      // 2*npoints. In such cases the array xline and yline must be arranged
      // before being plotted. The next loop does that.
      if (if2 > npoints) {
         for (i=1; i<if1; i++) {
            xline[if1-2+i] = xline[if2-1+i];
            yline[if1-2+i] = yline[if2-1+i];
         }
         npoints = if1-1;
      }
      if (option4) graph.PaintGraph(2*npoints,xline,yline,"FC");
      else         graph.PaintGraph(2*npoints,xline,yline,"F");
      gPad->SetLogx(logx);
      gPad->SetLogy(logy);
      delete [] xline;
      delete [] yline;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw 2D histograms errors.

void THistPainter::Paint2DErrors(Option_t *)
{

   fH->TAttMarker::Modify();
   fH->TAttLine::Modify();

   // Define the 3D view
   fXbuf[0] = Hparam.xmin;
   fYbuf[0] = Hparam.xmax;
   fXbuf[1] = Hparam.ymin;
   fYbuf[1] = Hparam.ymax;
   fXbuf[2] = Hparam.zmin;
   fYbuf[2] = Hparam.zmax*(1. + gStyle->GetHistTopMargin());
   fLego = new TPainter3dAlgorithms(fXbuf, fYbuf);
   TView *view = gPad->GetView();
   if (!view) {
      Error("Paint2DErrors", "no TView in current pad");
      return;
   }
   Double_t thedeg =  90 - gPad->GetTheta();
   Double_t phideg = -90 - gPad->GetPhi();
   Double_t psideg = view->GetPsi();
   Int_t irep;
   view->SetView(phideg, thedeg, psideg, irep);

   // Set color/style for back box
   fLego->SetFillStyle(gPad->GetFrameFillStyle());
   fLego->SetFillColor(gPad->GetFrameFillColor());
   fLego->TAttFill::Modify();
   Int_t backcolor = gPad->GetFrameFillColor();
   if (Hoption.System != kCARTESIAN) backcolor = 0;
   view->PadRange(backcolor);
   fLego->SetFillStyle(fH->GetFillStyle());
   fLego->SetFillColor(fH->GetFillColor());
   fLego->TAttFill::Modify();

   // Paint the Back Box if needed
   if (Hoption.BackBox && !Hoption.Same && !Hoption.Lego && !Hoption.Surf) {
      fLego->InitMoveScreen(-1.1,1.1);
      fLego->DefineGridLevels(fZaxis->GetNdivisions()%100);
      fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMove1);
      fLego->BackBox(90);
   }

   // Paint the Errors
   Double_t x, ex, x1, x2;
   Double_t y, ey, y1, y2;
   Double_t z, ez1, ez2, z1, z2;
   Double_t temp1[3],temp2[3];
   Double_t xyerror;
   if (Hoption.Error == 110) {
      xyerror = 0;
   } else {
      xyerror = gStyle->GetErrorX();
   }

   Double_t xk, xstep, yk, ystep;
   for (Int_t j=Hparam.yfirst; j<=Hparam.ylast;j++) {
      y  = fYaxis->GetBinCenter(j);
      ey = fYaxis->GetBinWidth(j)*xyerror;
      y1 = y-ey;
      y2 = y+ey;
      if (Hoption.Logy) {
         if (y > 0)  y = TMath::Log10(y);
         else        continue;
         if (y1 > 0) y1 = TMath::Log10(y1);
         else        y1 = Hparam.ymin;
         if (y2 > 0) y2 = TMath::Log10(y2);
         else        y2 = Hparam.ymin;
      }
      yk    = fYaxis->GetBinLowEdge(j);
      ystep = fYaxis->GetBinWidth(j);
      for (Int_t i=Hparam.xfirst; i<=Hparam.xlast;i++) {
         xk    = fXaxis->GetBinLowEdge(i);
         xstep = fXaxis->GetBinWidth(i);
         if (!IsInside(xk+0.5*xstep,yk+0.5*ystep)) continue;
         Int_t bin = fH->GetBin(i,j);
         x  = fXaxis->GetBinCenter(i);
         ex = fXaxis->GetBinWidth(i)*xyerror;
         x1 = x-ex;
         x2 = x+ex;
         if (Hoption.Logx) {
            if (x > 0)  x = TMath::Log10(x);
            else        continue;
            if (x1 > 0) x1 = TMath::Log10(x1);
            else        x1 = Hparam.xmin;
            if (x2 > 0) x2 = TMath::Log10(x2);
            else        x2 = Hparam.xmin;
         }
         z  = fH->GetBinContent(bin);
         if (fH->GetBinErrorOption() == TH1::kNormal) {
            ez1 = fH->GetBinError(bin);
            ez2 = ez1;
         }
         else {
            ez1 = fH->GetBinErrorLow(bin);
            ez2 = fH->GetBinErrorUp(bin);
         }
         z1 = z - ez1;
         z2 = z + ez2;
         if (Hoption.Logz) {
            if (z > 0)   z = TMath::Log10(z);
            else         z = Hparam.zmin;
            if (z1 > 0) z1 = TMath::Log10(z1);
            else        z1 = Hparam.zmin;
            if (z2 > 0) z2 = TMath::Log10(z2);
            else        z2 = Hparam.zmin;

         }
         if (z <= Hparam.zmin) continue;
         if (z >  Hparam.zmax) z = Hparam.zmax;

         temp1[0] = x1;
         temp1[1] = y;
         temp1[2] = z;
         temp2[0] = x2;
         temp2[1] = y;
         temp2[2] = z;
         gPad->PaintLine3D(temp1, temp2);
         temp1[0] = x;
         temp1[1] = y1;
         temp1[2] = z;
         temp2[0] = x;
         temp2[1] = y2;
         temp2[2] = z;
         gPad->PaintLine3D(temp1, temp2);
         temp1[0] = x;
         temp1[1] = y;
         temp1[2] = z1;
         temp2[0] = x;
         temp2[1] = y;
         temp2[2] = z2;
         gPad->PaintLine3D(temp1, temp2);
         temp1[0] = x;
         temp1[1] = y;
         temp1[2] = z;
         view->WCtoNDC(temp1, &temp2[0]);
         gPad->PaintPolyMarker(1, &temp2[0], &temp2[1]);
      }
   }

   // Paint the Front Box if needed
   if (Hoption.FrontBox) {
      fLego->InitMoveScreen(-1.1,1.1);
      fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMove2);
      fLego->FrontBox(90);
   }

   // Paint the Axis if needed
   if (!Hoption.Axis && !Hoption.Same && !Hoption.Lego && !Hoption.Surf) {
      TGaxis *axis = new TGaxis();
      PaintLegoAxis(axis, 90);
      delete axis;
   }

   delete fLego; fLego = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate range and clear pad (canvas).

void THistPainter::PaintFrame()
{

   if (Hoption.Same) return;

   RecalculateRange();

   if (Hoption.Lego || Hoption.Surf || Hoption.Tri ||
       Hoption.Contour == 14 || Hoption.Error >= 100) {
      TObject *frame = gPad->FindObject("TFrame");
      if (frame) gPad->GetListOfPrimitives()->Remove(frame);
      return;
   }

   //The next statement is always executed on non-iOS platform,
   //on iOS depends on pad mode.
   if (!gPad->PadInSelectionMode() && !gPad->PadInHighlightMode())
      gPad->PaintPadFrame(Hparam.xmin,Hparam.ymin,Hparam.xmax,Hparam.ymax);
}

////////////////////////////////////////////////////////////////////////////////
///  [Paint functions associated to an histogram.](#HP28")

void THistPainter::PaintFunction(Option_t *)
{

   TObjOptLink *lnk = (TObjOptLink*)fFunctions->FirstLink();
   TObject *obj;

   while (lnk) {
      obj = lnk->GetObject();
      TVirtualPad *padsave = gPad;
      if (obj->InheritsFrom(TF2::Class())) {
         if (obj->TestBit(TF2::kNotDraw) == 0) {
            if (Hoption.Lego || Hoption.Surf || Hoption.Error >= 100) {
               TF2 *f2 = (TF2*)obj;
               f2->SetMinimum(fH->GetMinimum());
               f2->SetMaximum(fH->GetMaximum());
               f2->SetRange(fH->GetXaxis()->GetXmin(), fH->GetYaxis()->GetXmin(), fH->GetXaxis()->GetXmax(), fH->GetYaxis()->GetXmax() );
               f2->Paint("surf same");
            } else {
               obj->Paint("cont3 same");
            }
         }
      } else if (obj->InheritsFrom(TF1::Class())) {
         if (obj->TestBit(TF1::kNotDraw) == 0) obj->Paint("lsame");
      } else  {
         //Let's make this 'function' selectable on iOS device (for example, it can be TPaveStat).
         gPad->PushSelectableObject(obj);

         //The next statement is ALWAYS executed on non-iOS platform, on iOS it depends on pad's mode
         //and picked object.
         if (!gPad->PadInHighlightMode() || (gPad->PadInHighlightMode() && obj == gPad->GetSelected()))
            obj->Paint(lnk->GetOption());
      }
      lnk = (TObjOptLink*)lnk->Next();
      padsave->cd();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// [Control routine to draw 1D histograms](#HP01b)

void THistPainter::PaintHist(Option_t *)
{

   //On iOS: do not highlight hist, if part of it was selected.
   //Never executes on non-iOS platform.
   if (gPad->PadInHighlightMode() && gPad->GetSelected() != fH)
      return;

   static char chopth[17];

   Int_t htype, oldhtype;
   Int_t i, j, first, last, nbins, fixbin;
   Double_t c1, yb;
   yb = 0;

   strlcpy(chopth, "                ",17);

   Double_t ymin = Hparam.ymin;
   Double_t ymax = Hparam.ymax;
   Double_t baroffset = fH->GetBarOffset();
   Double_t barwidth  = fH->GetBarWidth();
   Double_t baroffsetsave = gStyle->GetBarOffset();
   Double_t barwidthsave  = gStyle->GetBarWidth();
   gStyle->SetBarOffset(baroffset);
   gStyle->SetBarWidth(barwidth);

   //       Create "LIFE" structure to keep current histogram status

   first = Hparam.xfirst;
   last  = Hparam.xlast;
   nbins = last - first + 1;

   Double_t *keepx = 0;
   Double_t *keepy = 0;
   if (fXaxis->GetXbins()->fN) fixbin = 0;
   else                        fixbin = 1;
   if (fixbin) keepx = new Double_t[2];
   else        keepx = new Double_t[nbins+1];
   keepy = new Double_t[nbins];
   Double_t logymin = 0;
   if (Hoption.Logy) logymin = TMath::Power(10,ymin);

   //      Loop on histogram bins

   for (j=first; j<=last;j++) {
      c1 = Hparam.factor*fH->GetBinContent(j);
      if (TMath::Abs(ymax-ymin) > 0) {
         if (Hoption.Logy) yb = TMath::Log10(TMath::Max(c1,.1*logymin));
         else              yb = c1;
      }
      if (!Hoption.Line) {
         yb = TMath::Max(yb, ymin);
         yb = TMath::Min(yb, ymax);
      }
      keepy[j-first] = yb;
   }

   //              Draw histogram according to value of FillStyle and FillColor

   if (fixbin) { keepx[0] = Hparam.xmin; keepx[1] = Hparam.xmax; }
   else {
      for (i=0; i<nbins; i++) keepx[i] = fXaxis->GetBinLowEdge(i+first);
      keepx[nbins] = fXaxis->GetBinUpEdge(nbins-1+first);
   }

   //         Prepare Fill area (systematic with option "Bar").

   oldhtype = fH->GetFillStyle();
   htype    = oldhtype;
   if (Hoption.Bar) {
      if (htype == 0 || htype == 1000) htype = 1001;
   }

   Width_t lw = (Width_t)fH->GetLineWidth();

   //         Code option for GrapHist

   if (Hoption.Line) chopth[0] = 'L';
   if (Hoption.Star) chopth[1] = '*';
   if (Hoption.Mark) chopth[2] = 'P';
   if (Hoption.Mark == 10) chopth[3] = '0';
   if (Hoption.Line || Hoption.Curve || Hoption.Hist || Hoption.Bar) {
      if (Hoption.Curve)    chopth[3] = 'C';
      if (Hoption.Hist > 0) chopth[4] = 'H';
      else if (Hoption.Bar) chopth[5] = 'B';
      if (fH->GetFillColor() && htype) {
         if (Hoption.Logy) {
            chopth[6] = '1';
         }
         if (Hoption.Hist > 0 || Hoption.Curve || Hoption.Line) {
            chopth[7] = 'F';
         }
      }
   }
   if (!fixbin && strlen(chopth)) {
      chopth[8] = 'N';
   }

   if (Hoption.Fill == 2)    chopth[13] = '2';

   //         Option LOGX

   if (Hoption.Logx) {
      chopth[9]  = 'G';
      chopth[10] = 'X';
      if (fixbin) {
         keepx[0] = TMath::Power(10,keepx[0]);
         keepx[1] = TMath::Power(10,keepx[1]);
      }
   }

   if (Hoption.Off) {
      chopth[11] = ']';
      chopth[12] = '[';
   }

   //         Draw the histogram

   TGraph graph;
   graph.SetLineWidth(lw);
   graph.SetLineStyle(fH->GetLineStyle());
   graph.SetLineColor(fH->GetLineColor());
   graph.SetFillStyle(htype);
   graph.SetFillColor(fH->GetFillColor());
   graph.SetMarkerStyle(fH->GetMarkerStyle());
   graph.SetMarkerSize(fH->GetMarkerSize());
   graph.SetMarkerColor(fH->GetMarkerColor());
   if (!Hoption.Same) graph.ResetBit(TGraph::kClipFrame);

   graph.PaintGrapHist(nbins, keepx, keepy ,chopth);

   delete [] keepx;
   delete [] keepy;
   gStyle->SetBarOffset(baroffsetsave);
   gStyle->SetBarWidth(barwidthsave);

   htype=oldhtype;
}

////////////////////////////////////////////////////////////////////////////////
/// [Control function to draw a 3D histograms.](#HP01d)

void THistPainter::PaintH3(Option_t *option)
{

   char *cmd;
   TString opt = option;
   opt.ToLower();
   Int_t irep;

   if (fCurrentF3) {
      PaintTF3();
      return;
   } else if (Hoption.Box || Hoption.Lego) {
      if (Hoption.Box == 11 || Hoption.Lego == 11) {
         PaintH3Box(1);
      } else if (Hoption.Box == 12 || Hoption.Lego == 12) {
         PaintH3Box(2);
      } else if (Hoption.Box == 13 || Hoption.Lego == 13) {
         PaintH3Box(3);
      } else {
         PaintH3BoxRaster();
      }
      return;
   } else if (strstr(opt,"iso")) {
      PaintH3Iso();
      return;
   } else if (strstr(opt,"tf3")) {
      PaintTF3();
      return;
   } else {
      cmd = Form("TPolyMarker3D::PaintH3((TH1 *)0x%lx,\"%s\");",(Long_t)fH,option);
   }

   if (strstr(opt,"fb")) Hoption.FrontBox = 0;
   if (strstr(opt,"bb")) Hoption.BackBox = 0;

   TView *view = gPad->GetView();
   if (!view) return;
   Double_t thedeg =  90 - gPad->GetTheta();
   Double_t phideg = -90 - gPad->GetPhi();
   Double_t psideg = view->GetPsi();
   view->SetView(phideg, thedeg, psideg, irep);

   // Paint the data
   gROOT->ProcessLine(cmd);

   if (Hoption.Same) return;

   // Draw axis
   view->SetOutlineToCube();
   TSeqCollection *ol = view->GetOutline();
   if (ol && Hoption.BackBox && Hoption.FrontBox) ol->Paint(option);
   Hoption.System = kCARTESIAN;
   TGaxis *axis = new TGaxis();
   if (!Hoption.Axis && !Hoption.Same) PaintLegoAxis(axis, 90);
   delete axis;

   // Draw palette. In case of 4D plot with TTree::Draw() the palette should
   // be painted with the option colz.
   if (fH->GetDrawOption() && strstr(opt,"colz")) {
      Int_t ndiv   = fH->GetContour();
      if (ndiv == 0 ) {
         ndiv = gStyle->GetNumberContours();
         fH->SetContour(ndiv);
      }
      PaintPalette();
   }

   // Draw title
   PaintTitle();

   //Draw stats and fit results
   TF1 *fit  = 0;
   TIter next(fFunctions);
   TObject *obj;
   while ((obj = next())) {
      if (obj->InheritsFrom(TF1::Class())) {
         fit = (TF1*)obj;
         break;
      }
   }
   if ((Hoption.Same%10) != 1) {
      if (!fH->TestBit(TH1::kNoStats)) {  // bit set via TH1::SetStats
         PaintStat3(gStyle->GetOptStat(),fit);
      }
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Compute histogram parameters used by the drawing routines.

Int_t THistPainter::PaintInit()
{

   if (fH->GetDimension() > 1 || Hoption.Lego || Hoption.Surf) return 1;

   Int_t i;
   static const char *where = "PaintInit";
   Double_t yMARGIN = gStyle->GetHistTopMargin();
   Int_t maximum = 0;
   Int_t minimum = 0;
   if (fH->GetMaximumStored() != -1111) maximum = 1;
   if (fH->GetMinimumStored() != -1111) minimum = 1;

   //     Compute X axis parameters

   Int_t last      = fXaxis->GetLast();
   Int_t first     = fXaxis->GetFirst();
   Hparam.xlowedge = fXaxis->GetBinLowEdge(first);
   Hparam.xbinsize = fXaxis->GetBinWidth(first);
   Hparam.xlast    = last;
   Hparam.xfirst   = first;
   Hparam.xmin     = Hparam.xlowedge;
   Hparam.xmax     = fXaxis->GetBinLowEdge(last)+fXaxis->GetBinWidth(last);

   //       if log scale in X, replace xmin,max by the log
   if (Hoption.Logx) {
      if (Hparam.xmax<=0) {
         Error(where, "cannot set X axis to log scale");
         return 0;
      }
      if (Hparam.xlowedge <=0 ) {
         if (Hoption.Same) {
            Hparam.xlowedge = TMath::Power(10, gPad->GetUxmin());
         } else {
            for (i=first; i<=last; i++) {
               Double_t binLow = fXaxis->GetBinLowEdge(i);
               if (binLow>0) {
                  Hparam.xlowedge = binLow;
                  break;
               }
              if (binLow == 0 && fH->GetBinContent(i) !=0) {
                 Hparam.xlowedge = fXaxis->GetBinUpEdge(i)*0.001;
                 break;
              }
            }
            if (Hparam.xlowedge<=0) {
               Error(where, "cannot set X axis to log scale");
               return 0;
            }
         }
         Hparam.xmin  = Hparam.xlowedge;
      }
      Hparam.xfirst= fXaxis->FindFixBin(Hparam.xmin);
      Hparam.xlast = fXaxis->FindFixBin(Hparam.xmax);
      Hparam.xmin  = TMath::Log10(Hparam.xmin);
      Hparam.xmax  = TMath::Log10(Hparam.xmax);
      if (Hparam.xlast  > last)  Hparam.xlast  = last;
      if (Hparam.xfirst < first) Hparam.xfirst = first;
   }

   //     Compute Y axis parameters
   Double_t bigp = TMath::Power(10,32);
   Double_t ymax = -bigp;
   Double_t ymin = bigp;
   Double_t c1, e1;
   Double_t xv[1];
   Double_t fval;
   TObject *f;
   TF1 *f1;
   Double_t allchan = 0;
   Int_t nonNullErrors = 0;
   TIter   next(fFunctions);
   for (i=first; i<=last;i++) {
      c1 = fH->GetBinContent(i);
      ymax = TMath::Max(ymax,c1);
      if (Hoption.Logy) {
         if (c1 > 0) ymin = TMath::Min(ymin,c1);
      } else {
         ymin = TMath::Min(ymin,c1);
      }
      if (Hoption.Error) {
         if (fH->GetBinErrorOption() == TH1::kNormal)
            e1 = fH->GetBinError(i);
         else
            e1 = fH->GetBinErrorUp(i);
         if (e1 > 0) nonNullErrors++;
         ymax = TMath::Max(ymax,c1+e1);
         if (fH->GetBinErrorOption() != TH1::kNormal)
            e1 = fH->GetBinErrorLow(i);

         if (Hoption.Logy) {
            if (c1-e1>0.01*TMath::Abs(c1)) ymin = TMath::Min(ymin,c1-e1);
         } else {
            ymin = TMath::Min(ymin,c1-e1);
         }
      }
      if (Hoption.Func) {
         xv[0] = fXaxis->GetBinCenter(i);
         while ((f = (TObject*) next())) {
            if (f->IsA() == TF1::Class()) {
               f1 = (TF1*)f;
               if (xv[0] < f1->GetXmin() || xv[0] > f1->GetXmax()) continue;
               fval = f1->Eval(xv[0],0,0);
               if (f1->GetMaximumStored() != -1111) fval = TMath::Min(f1->GetMaximumStored(), fval);
               ymax = TMath::Max(ymax,fval);
               if (Hoption.Logy) {
                  if (c1 > 0 && fval > 0.3*c1) ymin = TMath::Min(ymin,fval);
               }
            }
         }
         next.Reset();
      }
      allchan += c1;
   }
   if (!nonNullErrors) {
      if (Hoption.Error) {
         if (!Hoption.Mark && !Hoption.Line && !Hoption.Star && !Hoption.Curve) Hoption.Hist = 2;
         Hoption.Error=0;
      }
   }


   //     Take into account maximum , minimum

   if (Hoption.Logy && ymin <= 0) {
      if (ymax >= 1) ymin = TMath::Max(.005,ymax*1e-10);
      else           ymin = 0.001*ymax;
   }

   Double_t xm = ymin;
   if (maximum) ymax = fH->GetMaximumStored();
   if (minimum) xm   = fH->GetMinimumStored();
   if (Hoption.Logy && xm < 0) {
      Error(where, "log scale requested with a negative argument (%f)", xm);
      return 0;
   } else if (Hoption.Logy && xm>=0 && ymax==0) { // empty histogram in log scale
      ymin = 0.01;
      ymax = 10.;
   } else {
      ymin = xm;
   }

   if (ymin >= ymax) {
      if (Hoption.Logy) {
         if (ymax > 0) ymin = 0.001*ymax;
         else {
            if (!Hoption.Same) Error(where, "log scale is requested but maximum is less or equal 0 (%f)", ymax);
            return 0;
         }
      }
      else {
         if (ymin > 0) {
            ymin = 0;
            ymax *= 2;
         } else if (ymin < 0) {
            ymax = 0;
            ymin *= 2;
         } else {
            ymin = 0;
            ymax = 1;
         }
      }
   }

   // In some cases, mainly because of precision issues, ymin and ymax could almost equal.
   if (TMath::AreEqualRel(ymin,ymax,1E-15)) {
      ymin = ymin*(1-1E-14);
      ymax = ymax*(1+1E-14);
   }

   //     take into account normalization factor
   Hparam.allchan = allchan;
   Double_t factor = allchan;
   if (fH->GetNormFactor() > 0) factor = fH->GetNormFactor();
   if (allchan) factor /= allchan;
   if (factor == 0) factor = 1;
   Hparam.factor = factor;
   ymax = factor*ymax;
   ymin = factor*ymin;
   //just in case the norm factor is negative
   // this may happen with a positive norm factor and a negative integral !
   if (ymax < ymin) {
      Double_t temp = ymax;
      ymax = ymin;
      ymin = temp;
   }

   //         For log scales, histogram coordinates are LOG10(ymin) and
   //         LOG10(ymax). Final adjustment (if not option "Same"
   //         or "+" for ymax) of ymax and ymin for logarithmic scale, if
   //         Maximum and Minimum are not defined.
   if (Hoption.Logy) {
      if (ymin <=0 || ymax <=0) {
         Error(where, "Cannot set Y axis to log scale");
         return 0;
      }
      ymin = TMath::Log10(ymin);
      if (!minimum) ymin += TMath::Log10(0.5);
      ymax = TMath::Log10(ymax);
      if (!maximum) ymax += TMath::Log10(2*(0.9/0.95));
      if (!Hoption.Same) {
         Hparam.ymin = ymin;
         Hparam.ymax = ymax;
      }
      return 1;
   }

   //         final adjustment of ymin for linear scale.
   //         if minimum is not set , then ymin is set to zero if >0
   //         or to ymin - margin if <0.
   if (!minimum) {
      if (Hoption.MinimumZero) {
         if (ymin >= 0) ymin = 0;
         else           ymin -= yMARGIN*(ymax-ymin);
      } else {
         Double_t dymin = yMARGIN*(ymax-ymin);
         if (ymin >= 0 && (ymin-dymin <= 0)) ymin  = 0;
         else                                ymin -= dymin;
      }
   }

   //         final adjustment of YMAXI for linear scale (if not option "Same"):
   //         decrease histogram height to MAX% of allowed height if HMAXIM
   //         has not been called.
   if (!maximum) {
      ymax += yMARGIN*(ymax-ymin);
   }

   Hparam.ymin = ymin;
   Hparam.ymax = ymax;
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute histogram parameters used by the drawing routines for a rotated pad.

Int_t THistPainter::PaintInitH()
{

   static const char *where = "PaintInitH";
   Double_t yMARGIN = gStyle->GetHistTopMargin();
   Int_t maximum = 0;
   Int_t minimum = 0;
   if (fH->GetMaximumStored() != -1111) maximum = 1;
   if (fH->GetMinimumStored() != -1111) minimum = 1;

   //     Compute X axis parameters

   Int_t last      = fXaxis->GetLast();
   Int_t first     = fXaxis->GetFirst();
   Hparam.xlowedge = fXaxis->GetBinLowEdge(first);
   Hparam.xbinsize = fXaxis->GetBinWidth(first);
   Hparam.xlast    = last;
   Hparam.xfirst   = first;
   Hparam.ymin     = Hparam.xlowedge;
   Hparam.ymax     = fXaxis->GetBinLowEdge(last)+fXaxis->GetBinWidth(last);

   //       if log scale in Y, replace ymin,max by the log
   if (Hoption.Logy) {
      if (Hparam.xlowedge <=0 ) {
         Hparam.xlowedge = 0.1*Hparam.xbinsize;
         Hparam.ymin  = Hparam.xlowedge;
      }
      if (Hparam.ymin <=0 || Hparam.ymax <=0) {
         Error(where, "cannot set Y axis to log scale");
         return 0;
      }
      Hparam.xfirst= fXaxis->FindFixBin(Hparam.ymin);
      Hparam.xlast = fXaxis->FindFixBin(Hparam.ymax);
      Hparam.ymin  = TMath::Log10(Hparam.ymin);
      Hparam.ymax  = TMath::Log10(Hparam.ymax);
      if (Hparam.xlast > last) Hparam.xlast = last;
   }

   //     Compute Y axis parameters
   Double_t bigp = TMath::Power(10,32);
   Double_t xmax = -bigp;
   Double_t xmin = bigp;
   Double_t c1, e1;
   Double_t xv[1];
   Double_t fval;
   Int_t i;
   TObject *f;
   TF1 *f1;
   Double_t allchan = 0;
   TIter   next(fFunctions);
   for (i=first; i<=last;i++) {
      c1 = fH->GetBinContent(i);
      xmax = TMath::Max(xmax,c1);
      xmin = TMath::Min(xmin,c1);
      if (Hoption.Error) {
         e1 = fH->GetBinError(i);
         xmax = TMath::Max(xmax,c1+e1);
         xmin = TMath::Min(xmin,c1-e1);
      }
      if (Hoption.Func) {
         xv[0] = fXaxis->GetBinCenter(i);
         while ((f = (TObject*) next())) {
            if (f->IsA() == TF1::Class()) {
               f1 = (TF1*)f;
               if (xv[0] < f1->GetXmin() || xv[0] > f1->GetXmax()) continue;
               fval = f1->Eval(xv[0],0,0);
               xmax = TMath::Max(xmax,fval);
               if (Hoption.Logy) {
                  if (fval > 0.3*c1) xmin = TMath::Min(xmin,fval);
               }
            }
         }
         next.Reset();
      }
      allchan += c1;
   }

   //     Take into account maximum , minimum

   if (Hoption.Logx && xmin <= 0) {
      if (xmax >= 1) xmin = TMath::Max(.5,xmax*1e-10);
      else           xmin = 0.001*xmax;
   }
   Double_t xm = xmin;
   if (maximum) xmax = fH->GetMaximumStored();
   if (minimum) xm   = fH->GetMinimumStored();
   if (Hoption.Logx && xm <= 0) {
      Error(where, "log scale requested with zero or negative argument (%f)", xm);
      return 0;
   }
   else xmin = xm;
   if (xmin >= xmax) {
      if (Hoption.Logx) {
         if (xmax > 0) xmin = 0.001*xmax;
         else {
            if (!Hoption.Same) Error(where, "log scale is requested but maximum is less or equal 0 (%f)", xmax);
            return 0;
         }
      }
      else {
         if (xmin > 0) {
            xmin = 0;
            xmax *= 2;
         } else if (xmin < 0) {
            xmax = 0;
            xmin *= 2;
         } else {
            xmin = -1;
            xmax = 1;
         }
      }
   }

   //     take into account normalization factor
   Hparam.allchan = allchan;
   Double_t factor = allchan;
   if (fH->GetNormFactor() > 0) factor = fH->GetNormFactor();
   if (allchan) factor /= allchan;
   if (factor == 0) factor = 1;
   Hparam.factor = factor;
   xmax = factor*xmax;
   xmin = factor*xmin;

   //         For log scales, histogram coordinates are LOG10(ymin) and
   //         LOG10(ymax). Final adjustment (if not option "Same"
   //         or "+" for ymax) of ymax and ymin for logarithmic scale, if
   //         Maximum and Minimum are not defined.
   if (Hoption.Logx) {
      if (xmin <=0 || xmax <=0) {
         Error(where, "Cannot set Y axis to log scale");
         return 0;
      }
      xmin = TMath::Log10(xmin);
      if (!minimum) xmin += TMath::Log10(0.5);
      xmax = TMath::Log10(xmax);
      if (!maximum) xmax += TMath::Log10(2*(0.9/0.95));
      if (!Hoption.Same) {
         Hparam.xmin = xmin;
         Hparam.xmax = xmax;
      }
      return 1;
   }

   //         final adjustment of ymin for linear scale.
   //         if minimum is not set , then ymin is set to zero if >0
   //         or to ymin - margin if <0.
   if (!minimum) {
      if (xmin >= 0) xmin = 0;
      else           xmin -= yMARGIN*(xmax-xmin);
   }

   //         final adjustment of YMAXI for linear scale (if not option "Same"):
   //         decrease histogram height to MAX% of allowed height if HMAXIM
   //         has not been called.
   if (!maximum) {
      xmax += yMARGIN*(xmax-xmin);
   }
   Hparam.xmin = xmin;
   Hparam.xmax = xmax;
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// [Control function to draw a 3D histogram with boxes.](#HP25)

void THistPainter::PaintH3Box(Int_t iopt)
{
   //       Predefined box structure
   Double_t wxyz[8][3] = { {-1,-1,-1}, {1,-1,-1}, {1,1,-1}, {-1,1,-1},
                           {-1,-1, 1}, {1,-1, 1}, {1,1, 1}, {-1,1, 1} };
   Int_t iface[6][4] = { {0,3,2,1}, {4,5,6,7},
                         {0,1,5,4}, {1,2,6,5}, {2,3,7,6}, {3,0,4,7} };

   //       Define dimensions of world space
   TGaxis *axis = new TGaxis();
   TAxis *xaxis = fH->GetXaxis();
   TAxis *yaxis = fH->GetYaxis();
   TAxis *zaxis = fH->GetZaxis();

   fXbuf[0] = xaxis->GetBinLowEdge(xaxis->GetFirst());
   fYbuf[0] = xaxis->GetBinUpEdge(xaxis->GetLast());
   fXbuf[1] = yaxis->GetBinLowEdge(yaxis->GetFirst());
   fYbuf[1] = yaxis->GetBinUpEdge(yaxis->GetLast());
   fXbuf[2] = zaxis->GetBinLowEdge(zaxis->GetFirst());
   fYbuf[2] = zaxis->GetBinUpEdge(zaxis->GetLast());

   fLego = new TPainter3dAlgorithms(fXbuf, fYbuf);

   //       Set view
   TView *view = gPad->GetView();
   if (!view) {
      Error("PaintH3", "no TView in current pad");
      return;
   }
   Double_t thedeg =  90 - gPad->GetTheta();
   Double_t phideg = -90 - gPad->GetPhi();
   Double_t psideg = view->GetPsi();
   Int_t irep;
   view->SetView(phideg, thedeg, psideg, irep);

   Int_t backcolor = gPad->GetFrameFillColor();
   view->PadRange(backcolor);

   //       Draw back surfaces of frame box
   fLego->InitMoveScreen(-1.1,1.1);
   if (Hoption.BackBox) {
      fLego->DefineGridLevels(fZaxis->GetNdivisions()%100);
      fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMove1);
      fLego->BackBox(90);
   }

   fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMode1);

   //       Define order of drawing
   Double_t *tnorm = view->GetTnorm();
   if (!tnorm) return;
   Int_t incrx = (tnorm[ 8] < 0.) ? -1 : +1;
   Int_t incry = (tnorm[ 9] < 0.) ? -1 : +1;
   Int_t incrz = (tnorm[10] < 0.) ? -1 : +1;
   Int_t ix1 = (incrx == +1) ? xaxis->GetFirst() : xaxis->GetLast();
   Int_t iy1 = (incry == +1) ? yaxis->GetFirst() : yaxis->GetLast();
   Int_t iz1 = (incrz == +1) ? zaxis->GetFirst() : zaxis->GetLast();
   Int_t ix2 = (incrx == +1) ? xaxis->GetLast()  : xaxis->GetFirst();
   Int_t iy2 = (incry == +1) ? yaxis->GetLast()  : yaxis->GetFirst();
   Int_t iz2 = (incrz == +1) ? zaxis->GetLast()  : zaxis->GetFirst();

   //       Set graphic attributes (colour, style, etc.)
   Style_t fillsav   = fH->GetFillStyle();
   Style_t colsav    = fH->GetFillColor();
   Style_t coldark   = TColor::GetColorDark(colsav);
   Style_t colbright = TColor::GetColorBright(colsav);

   fH->SetFillStyle(1001);
   fH->TAttFill::Modify();
   fH->TAttLine::Modify();
   Int_t ncolors  = gStyle->GetNumberOfColors();
   Int_t theColor;

   //       Create bin boxes and draw
   Double_t wmin = TMath::Max(fH->GetMinimum(),0.);
   Double_t wmax = TMath::Max(TMath::Abs(fH->GetMaximum()),
                              TMath::Abs(fH->GetMinimum()));

   Double_t pmin[3], pmax[3], sxyz[8][3];
   for (Int_t ix = ix1; ix !=ix2+incrx; ix += incrx) {
      pmin[0] = xaxis->GetBinLowEdge(ix);
      pmax[0] = xaxis->GetBinUpEdge(ix);
      for (Int_t iy = iy1; iy != iy2+incry; iy += incry) {
         pmin[1] = yaxis->GetBinLowEdge(iy);
         pmax[1] = yaxis->GetBinUpEdge(iy);
         for (Int_t iz = iz1; iz != iz2+incrz; iz += incrz) {
            pmin[2]    = zaxis->GetBinLowEdge(iz);
            pmax[2]    = zaxis->GetBinUpEdge(iz);
            Double_t w = fH->GetBinContent(fH->GetBin(ix,iy,iz));
            Bool_t neg = kFALSE;
            Int_t n = 5;
            if (w<0) {
               w   = -w;
               neg = kTRUE;
            }
            if (w < wmin) continue;
            if (w > wmax) w = wmax;
            Double_t scale = (TMath::Power((w-wmin)/(wmax-wmin),1./3.))/2.;
            if (scale == 0) continue;
            for (Int_t i=0; i<3; ++i) {
               Double_t c = (pmax[i] + pmin[i])*0.5;
               Double_t d = (pmax[i] - pmin[i])*scale;
               for (Int_t k=0; k<8; ++k) { // set bin box vertices
                  sxyz[k][i] = wxyz[k][i]*d + c;
               }
            }
            for (Int_t k=0; k<8; ++k) { // transform to normalized space
               view->WCtoNDC(&sxyz[k][0],&sxyz[k][0]);
            }
            Double_t x[8], y[8]; // draw bin box faces
            for (Int_t k=0; k<6; ++k) {
               for (Int_t i=0; i<4; ++i) {
                  Int_t iv = iface[k][i];
                  x[i]     = sxyz[iv][0];
                  y[i]     = sxyz[iv][1];
               }
               x[4] = x[0] ; y[4] = y[0];
               if (neg) {
                  x[5] = x[2] ; y[5] = y[2];
                  x[6] = x[3] ; y[6] = y[3];
                  x[7] = x[1] ; y[7] = y[1];
                  n = 8;
               } else {
                  n = 5;
               }
               Double_t z = (x[2]-x[0])*(y[3]-y[1]) - (y[2]-y[0])*(x[3]-x[1]);
               if (z <= 0.) continue;
               if (iopt == 2) {
                  theColor = ncolors*((w-wmin)/(wmax-wmin)) -1;
                  fH->SetFillColor(gStyle->GetColorPalette(theColor));
               } else {
                  if (k == 3 || k == 5) {
                     fH->SetFillColor(coldark);
                  } else if (k == 0 || k == 1) {
                     fH->SetFillColor(colbright);
                  } else {
                     fH->SetFillColor(colsav);
                  }
               }
               fH->TAttFill::Modify();
               gPad->PaintFillArea(4, x, y);
               if (iopt != 3)gPad->PaintPolyLine(n, x, y);
            }
         }
      }
   }

   //       Draw front surfaces of frame box
   if (Hoption.FrontBox) fLego->FrontBox(90);

   //       Draw axis and title
   if (!Hoption.Axis && !Hoption.Same) PaintLegoAxis(axis, 90);
   PaintTitle();

   // Draw palette. if needed.
   if (Hoption.Zscale) {
      Int_t ndiv   = fH->GetContour();
      if (ndiv == 0 ) {
         ndiv = gStyle->GetNumberContours();
         fH->SetContour(ndiv);
      }
      PaintPalette();
   }

   delete axis;
   delete fLego; fLego = 0;

   fH->SetFillStyle(fillsav);
   fH->SetFillColor(colsav);
   fH->TAttFill::Modify();
}

////////////////////////////////////////////////////////////////////////////////
/// [Control function to draw a 3D histogram with boxes.](#HP25)

void THistPainter::PaintH3BoxRaster()
{
   //       Predefined box structure
   Double_t wxyz[8][3] = {
      {-1,-1,-1}, {1,-1,-1}, {1,1,-1}, {-1,1,-1}, // bottom vertices
      {-1,-1, 1}, {1,-1, 1}, {1,1, 1}, {-1,1, 1}  // top vertices
   };
   Int_t iface[6][4] = {
      {0,3,2,1}, {4,5,6,7},                       // bottom and top faces
      {0,1,5,4}, {1,2,6,5}, {2,3,7,6}, {3,0,4,7}  // side faces
   };
   Double_t normal[6][3] = {
      {0,0,-1}, {0,0,1},                          // Z-, Z+
      {0,-1,0}, {1,0,0}, {0,1,0}, {-1,0,0}        // Y-, X+, Y+, X-
   };

   //       Define dimensions of world space
   TGaxis *axis = new TGaxis();
   TAxis *xaxis = fH->GetXaxis();
   TAxis *yaxis = fH->GetYaxis();
   TAxis *zaxis = fH->GetZaxis();

   fXbuf[0] = xaxis->GetBinLowEdge(xaxis->GetFirst());
   fYbuf[0] = xaxis->GetBinUpEdge(xaxis->GetLast());
   fXbuf[1] = yaxis->GetBinLowEdge(yaxis->GetFirst());
   fYbuf[1] = yaxis->GetBinUpEdge(yaxis->GetLast());
   fXbuf[2] = zaxis->GetBinLowEdge(zaxis->GetFirst());
   fYbuf[2] = zaxis->GetBinUpEdge(zaxis->GetLast());

   fLego = new TPainter3dAlgorithms(fXbuf, fYbuf);

   //       Set view
   TView *view = gPad->GetView();
   if (!view) {
      Error("PaintH3", "no TView in current pad");
      return;
   }
   Double_t thedeg =  90 - gPad->GetTheta();
   Double_t phideg = -90 - gPad->GetPhi();
   Double_t psideg = view->GetPsi();
   Int_t irep;
   view->SetView(phideg, thedeg, psideg, irep);

   Int_t backcolor = gPad->GetFrameFillColor();
   view->PadRange(backcolor);

   //       Draw front surfaces of frame box
   if (Hoption.FrontBox) {
      fLego->InitMoveScreen(-1.1,1.1);
      fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMove2);
   }

   //       Initialize hidden line removal algorithm "raster screen"
   fLego->InitRaster(-1.1,-1.1,1.1,1.1,1000,800);

   //       Define order of drawing
   Double_t *tnorm = view->GetTnorm();
   if (!tnorm) return;
   Int_t incrx = (tnorm[ 8] < 0.) ? +1 : -1;
   Int_t incry = (tnorm[ 9] < 0.) ? +1 : -1;
   Int_t incrz = (tnorm[10] < 0.) ? +1 : -1;
   Int_t ix1 = (incrx == +1) ? xaxis->GetFirst() : xaxis->GetLast();
   Int_t iy1 = (incry == +1) ? yaxis->GetFirst() : yaxis->GetLast();
   Int_t iz1 = (incrz == +1) ? zaxis->GetFirst() : zaxis->GetLast();
   Int_t ix2 = (incrx == +1) ? xaxis->GetLast()  : xaxis->GetFirst();
   Int_t iy2 = (incry == +1) ? yaxis->GetLast()  : yaxis->GetFirst();
   Int_t iz2 = (incrz == +1) ? zaxis->GetLast()  : zaxis->GetFirst();

   //       Set line attributes (colour, style, etc.)
   fH->TAttLine::Modify();

   //       Create bin boxes and draw
   const Int_t NTMAX = 100;
   Double_t tt[NTMAX][2];
   Double_t wmin = TMath::Max(fH->GetMinimum(),0.);
   Double_t wmax = TMath::Max(TMath::Abs(fH->GetMaximum()),
                              TMath::Abs(fH->GetMinimum()));
   Double_t pmin[3], pmax[3], sxyz[8][3], pp[4][2];
   for (Int_t ix = ix1; ix !=ix2+incrx; ix += incrx) {
      pmin[0] = xaxis->GetBinLowEdge(ix);
      pmax[0] = xaxis->GetBinUpEdge(ix);
      for (Int_t iy = iy1; iy != iy2+incry; iy += incry) {
         pmin[1] = yaxis->GetBinLowEdge(iy);
         pmax[1] = yaxis->GetBinUpEdge(iy);
         for (Int_t iz = iz1; iz != iz2+incrz; iz += incrz) {
            pmin[2] = zaxis->GetBinLowEdge(iz);
            pmax[2] = zaxis->GetBinUpEdge(iz);
            Double_t w = fH->GetBinContent(fH->GetBin(ix,iy,iz));
            Bool_t neg = kFALSE;
            if (w<0) {
               w   = -w;
               neg = kTRUE;
            }
            if (w < wmin) continue;
            if (w > wmax) w = wmax;
            Double_t scale = (TMath::Power((w-wmin)/(wmax-wmin),1./3.))/2.;
            if (scale == 0) continue;
            for (Int_t i=0; i<3; ++i) {
               Double_t c = (pmax[i] + pmin[i])*0.5;
               Double_t d = (pmax[i] - pmin[i])*scale;
               for (Int_t k=0; k<8; ++k) { // set bin box vertices
                  sxyz[k][i] = wxyz[k][i]*d + c;
               }
            }
            for (Int_t k=0; k<8; ++k) { // transform to normalized space
               view->WCtoNDC(&sxyz[k][0],&sxyz[k][0]);
            }
            for (Int_t k=0; k<6; ++k) { // draw box faces
               Double_t zn;
               view->FindNormal(normal[k][0], normal[k][1], normal[k][2], zn);
               if (zn <= 0) continue;
               for (Int_t i=0; i<4; ++i) {
                  Int_t ip = iface[k][i];
                  pp[i][0] = sxyz[ip][0];
                  pp[i][1] = sxyz[ip][1];
               }
               for (Int_t i=0; i<4; ++i) {
                  Int_t i1 = i;
                  Int_t i2 = (i == 3) ? 0 : i + 1;
                  Int_t nt;
                  fLego->FindVisibleLine(&pp[i1][0], &pp[i2][0], NTMAX, nt, &tt[0][0]);
                  Double_t xdel = pp[i2][0] - pp[i1][0];
                  Double_t ydel = pp[i2][1] - pp[i1][1];
                  Double_t x[2], y[2];
                  for (Int_t it = 0; it < nt; ++it) {
                     x[0] = pp[i1][0] + xdel*tt[it][0];
                     y[0] = pp[i1][1] + ydel*tt[it][0];
                     x[1] = pp[i1][0] + xdel*tt[it][1];
                     y[1] = pp[i1][1] + ydel*tt[it][1];
                     gPad->PaintPolyLine(2, x, y);
                  }
               }
               if (neg) {
                  Int_t i1 = 0;
                  Int_t i2 = 2;
                  Int_t nt;
                  fLego->FindVisibleLine(&pp[i1][0], &pp[i2][0], NTMAX, nt, &tt[0][0]);
                  Double_t xdel = pp[i2][0] - pp[i1][0];
                  Double_t ydel = pp[i2][1] - pp[i1][1];
                  Double_t x[2], y[2];
                  for (Int_t it = 0; it < nt; ++it) {
                     x[0] = pp[i1][0] + xdel*tt[it][0];
                     y[0] = pp[i1][1] + ydel*tt[it][0];
                     x[1] = pp[i1][0] + xdel*tt[it][1];
                     y[1] = pp[i1][1] + ydel*tt[it][1];
                     gPad->PaintPolyLine(2, x, y);
                  }
                  i1 = 1;
                  i2 = 3;
                  fLego->FindVisibleLine(&pp[i1][0], &pp[i2][0], NTMAX, nt, &tt[0][0]);
                  xdel = pp[i2][0] - pp[i1][0];
                  ydel = pp[i2][1] - pp[i1][1];
                  for (Int_t it = 0; it < nt; ++it) {
                     x[0] = pp[i1][0] + xdel*tt[it][0];
                     y[0] = pp[i1][1] + ydel*tt[it][0];
                     x[1] = pp[i1][0] + xdel*tt[it][1];
                     y[1] = pp[i1][1] + ydel*tt[it][1];
                     gPad->PaintPolyLine(2, x, y);
                  }
               }
               fLego->FillPolygonBorder(4, &pp[0][0]); // update raster screen
            }
         }
      }
   }

   //       Draw frame box
   if (Hoption.BackBox) {
      fLego->DefineGridLevels(fZaxis->GetNdivisions()%100);
      fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceRaster1);
      fLego->BackBox(90);
   }

   if (Hoption.FrontBox) fLego->FrontBox(90);

   //       Draw axis and title
   if (!Hoption.Axis && !Hoption.Same) PaintLegoAxis(axis, 90);
   PaintTitle();

   delete axis;
   delete fLego; fLego = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// [Control function to draw a 3D histogram with Iso Surfaces.](#HP25)

void THistPainter::PaintH3Iso()
{

   const Double_t ydiff = 1;
   const Double_t yligh1 = 10;
   const Double_t qa = 0.15;
   const Double_t qd = 0.15;
   const Double_t qs = 0.8;
   Double_t fmin, fmax;
   Int_t i, irep;
   Int_t nbcol = 28;
   Int_t icol1 = 201;
   Int_t ic1 = icol1;
   Int_t ic2 = ic1+nbcol;
   Int_t ic3 = ic2+nbcol;

   TGaxis *axis = new TGaxis();
   TAxis *xaxis = fH->GetXaxis();
   TAxis *yaxis = fH->GetYaxis();
   TAxis *zaxis = fH->GetZaxis();

   Int_t nx = fH->GetNbinsX();
   Int_t ny = fH->GetNbinsY();
   Int_t nz = fH->GetNbinsZ();

   Double_t *x = new Double_t[nx];
   Double_t *y = new Double_t[ny];
   Double_t *z = new Double_t[nz];

   for (i=0; i<nx; i++) x[i] = xaxis->GetBinCenter(i+1);
   for (i=0; i<ny; i++) y[i] = yaxis->GetBinCenter(i+1);
   for (i=0; i<nz; i++) z[i] = zaxis->GetBinCenter(i+1);

   fXbuf[0] = xaxis->GetBinLowEdge(xaxis->GetFirst());
   fYbuf[0] = xaxis->GetBinUpEdge(xaxis->GetLast());
   fXbuf[1] = yaxis->GetBinLowEdge(yaxis->GetFirst());
   fYbuf[1] = yaxis->GetBinUpEdge(yaxis->GetLast());
   fXbuf[2] = zaxis->GetBinLowEdge(zaxis->GetFirst());
   fYbuf[2] = zaxis->GetBinUpEdge(zaxis->GetLast());

   Double_t s[3];
   s[0] = fH->GetSumOfWeights()/(fH->GetNbinsX()*fH->GetNbinsY()*fH->GetNbinsZ());
   s[1] = 0.5*s[0];
   s[2] = 1.5*s[0];

   fLego = new TPainter3dAlgorithms(fXbuf, fYbuf);

   TView *view = gPad->GetView();
   if (!view) {
      Error("PaintH3Iso", "no TView in current pad");
      delete [] x;
      delete [] y;
      delete [] z;
      return;
   }
   Double_t thedeg =  90 - gPad->GetTheta();
   Double_t phideg = -90 - gPad->GetPhi();
   Double_t psideg = view->GetPsi();
   view->SetView(phideg, thedeg, psideg, irep);

   Int_t backcolor = gPad->GetFrameFillColor();
   if (Hoption.System != kCARTESIAN) backcolor = 0;
   view->PadRange(backcolor);

   Double_t dcol = 0.5/Double_t(nbcol);
   TColor *colref = gROOT->GetColor(fH->GetFillColor());
   if (!colref) {
      delete [] x;
      delete [] y;
      delete [] z;
      return;
   }
   Float_t r, g, b, hue, light, satur;
   colref->GetRGB(r,g,b);
   TColor::RGBtoHLS(r,g,b,hue,light,satur);
   TColor *acol;
   for (Int_t col=0;col<nbcol;col++) {
      acol = gROOT->GetColor(col+icol1);
      TColor::HLStoRGB(hue, .4+col*dcol, satur, r, g, b);
      if (acol) acol->SetRGB(r, g, b);
   }

   fLego->InitMoveScreen(-1.1,1.1);

   if (Hoption.BackBox) {
      fLego->DefineGridLevels(fZaxis->GetNdivisions()%100);
      fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMove1);
      fLego->BackBox(90);
   }

   fLego->LightSource(0, ydiff, 0, 0, 0, irep);
   fLego->LightSource(1, yligh1, 1, 1, 1, irep);
   fLego->SurfaceProperty(qa, qd, qs, 1, irep);
   fmin = ydiff*qa;
   fmax = ydiff*qa + (yligh1+0.1)*(qd+qs);
   fLego->SetIsoSurfaceParameters(fmin, fmax, nbcol, ic1, ic2, ic3);

   fLego->IsoSurface(1, s, nx, ny, nz, x, y, z, "BF");

   if (Hoption.FrontBox) {
      fLego->InitMoveScreen(-1.1,1.1);
      fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMove2);
      fLego->FrontBox(90);
   }
   if (!Hoption.Axis && !Hoption.Same) PaintLegoAxis(axis, 90);

   PaintTitle();

   delete axis;
   delete fLego; fLego = 0;
   delete [] x;
   delete [] y;
   delete [] z;
}

////////////////////////////////////////////////////////////////////////////////
/// [Control function to draw a 2D histogram as a lego plot.](#HP17)

void THistPainter::PaintLego(Option_t *)
{

   Int_t raster = 1;
   if (Hparam.zmin == 0 && Hparam.zmax == 0) {Hparam.zmin = -1; Hparam.zmax = 1;}
   Int_t   nx      = Hparam.xlast - Hparam.xfirst + 1;
   Int_t   ny      = Hparam.ylast - Hparam.yfirst + 1;
   Double_t zmin   = Hparam.zmin;
   Double_t zmax   = Hparam.zmax;
   Double_t xlab1  = Hparam.xmin;
   Double_t xlab2  = Hparam.xmax;
   Double_t ylab1  = Hparam.ymin;
   Double_t ylab2  = Hparam.ymax;
   Double_t dangle = 10*3.141592/180; //Delta angle for Rapidity option
   Double_t deltaz = TMath::Abs(zmin);
   if (deltaz == 0) deltaz = 1;
   if (zmin >= zmax) {
      zmin -= 0.5*deltaz;
      zmax += 0.5*deltaz;
   }
   Double_t z1c = zmin;
   Double_t z2c = zmin + (zmax-zmin)*(1+gStyle->GetHistTopMargin());

   //     Compute the lego limits and instantiate a lego object
   fXbuf[0] = -1;
   fYbuf[0] =  1;
   fXbuf[1] = -1;
   fYbuf[1] =  1;
   if (Hoption.System == kPOLAR) {
      fXbuf[2] = z1c;
      fYbuf[2] = z2c;
   } else if (Hoption.System == kCYLINDRICAL) {
      if (Hoption.Logy) {
         if (ylab1 > 0) fXbuf[2] = TMath::Log10(ylab1);
         else           fXbuf[2] = 0;
         if (ylab2 > 0) fYbuf[2] = TMath::Log10(ylab2);
         else           fYbuf[2] = 0;
      } else {
         fXbuf[2] = ylab1;
         fYbuf[2] = ylab2;
      }
      z1c = 0; z2c = 1;
   } else if (Hoption.System == kSPHERICAL) {
      fXbuf[2] = -1;
      fYbuf[2] =  1;
      z1c = 0; z2c = 1;
   } else if (Hoption.System == kRAPIDITY) {
      fXbuf[2] = -1/TMath::Tan(dangle);
      fYbuf[2] =  1/TMath::Tan(dangle);
   } else {
      fXbuf[0] = xlab1;
      fYbuf[0] = xlab2;
      fXbuf[1] = ylab1;
      fYbuf[1] = ylab2;
      fXbuf[2] = z1c;
      fYbuf[2] = z2c;
      raster  = 0;
   }

   fLego = new TPainter3dAlgorithms(fXbuf, fYbuf, Hoption.System);

   Int_t nids = -1;
   TH1 * hid = NULL;
   Color_t colormain = -1, colordark = -1;
   Bool_t drawShadowsInLego1 = kTRUE;

   // LEGO3 is like LEGO1 except that the black lines around each lego are not drawn.
   if (Hoption.Lego == 13) {
      Hoption.Lego = 11;
      fLego->SetMesh(0);
   }
   // LEGO4 is like LEGO1 except no shadows are drawn.
   if (Hoption.Lego == 14) {
      Hoption.Lego = 11;
      drawShadowsInLego1 = kFALSE;
   }

   //          Create axis object

   TGaxis *axis = new TGaxis();

   //                  Initialize the levels on the Z axis
   Int_t ndiv   = fH->GetContour();
   if (ndiv == 0 ) {
      ndiv = gStyle->GetNumberContours();
      fH->SetContour(ndiv);
   }
   Int_t ndivz  = TMath::Abs(ndiv);
   if (fH->TestBit(TH1::kUserContour) == 0) fH->SetContour(ndiv);

   //     Initialize colors
   if (!fStack) {
      fLego->SetEdgeAtt(fH->GetLineColor(),fH->GetLineStyle(),fH->GetLineWidth(),0);
   } else {
      for (Int_t id=0;id<=fStack->GetSize();id++) {
         hid = (TH1*)fStack->At((id==0)?id:id-1);
         fLego->SetEdgeAtt(hid->GetLineColor(),hid->GetLineStyle(),hid->GetLineWidth(),id);
      }
   }

   if (Hoption.Lego == 11) {
      nids = 1;
      if (fStack) nids = fStack->GetSize();
      hid = fH;
      for (Int_t id=0;id<=nids;id++) {
         if (id > 0 && fStack) hid = (TH1*)fStack->At(id-1);
         colormain = hid->GetFillColor();
         if (colormain == 1) colormain = 17; //avoid drawing with black
         if (drawShadowsInLego1) colordark = TColor::GetColorDark(colormain);
         else                    colordark = colormain;
         fLego->SetColorMain(colormain,id);
         fLego->SetColorDark(colordark,id);
         if (id <= 1)    fLego->SetColorMain(colormain,-1);  // Set Bottom color
         if (id == nids) fLego->SetColorMain(colormain,99);  // Set Top color
      }
   }

   //     Now ready to draw the lego plot
   Int_t irep = 0;

   TView *view = gPad->GetView();
   if (!view) {
      Error("PaintLego", "no TView in current pad");
      return;
   }

   Double_t thedeg =  90 - gPad->GetTheta();
   Double_t phideg = -90 - gPad->GetPhi();
   Double_t psideg = view->GetPsi();
   view->SetView(phideg, thedeg, psideg, irep);

   fLego->SetLineColor(kBlack);  // zgrid color for lego1 & lego2
   fLego->SetFillStyle(fH->GetFillStyle());

   //     Set color/style for back box
   fLego->SetFillStyle(gPad->GetFrameFillStyle());
   fLego->SetFillColor(gPad->GetFrameFillColor());
   fLego->TAttFill::Modify();

   Int_t backcolor = gPad->GetFrameFillColor();
   if (Hoption.System != kCARTESIAN) backcolor = 0;
   view->PadRange(backcolor);

   fLego->SetFillStyle(fH->GetFillStyle());
   fLego->SetFillColor(fH->GetFillColor());
   fLego->TAttFill::Modify();

   fLego->DefineGridLevels(fZaxis->GetNdivisions()%100);

   if (raster) fLego->InitRaster(-1.1,-1.1,1.1,1.1,1000,800);
   else        fLego->InitMoveScreen(-1.1,1.1);

   if (Hoption.Lego == 19) {
      fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMove1);
      if (Hoption.BackBox)   fLego->BackBox(90);
      if (Hoption.FrontBox)  fLego->FrontBox(90);
      if (!Hoption.Axis)     PaintLegoAxis(axis, 90);
      return;
   }

   if (Hoption.Lego == 11 || Hoption.Lego == 12) {
      if (Hoption.System == kCARTESIAN && Hoption.BackBox) {
         fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMove1);
         fLego->BackBox(90);
      }
   }

   if (Hoption.Lego == 12) DefineColorLevels(ndivz);

   fLego->SetLegoFunction(&TPainter3dAlgorithms::LegoFunction);
   if (Hoption.Lego ==  1) fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceRaster2);
   if (Hoption.Lego == 11) fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMode3);
   if (Hoption.Lego == 12) fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMode2);
   if (Hoption.System == kPOLAR) {
      if (Hoption.Lego ==  1) fLego->LegoPolar(1,nx,ny,"FB");
      if (Hoption.Lego == 11) fLego->LegoPolar(1,nx,ny,"BF");
      if (Hoption.Lego == 12) fLego->LegoPolar(1,nx,ny,"BF");
   } else if (Hoption.System == kCYLINDRICAL) {
      if (Hoption.Lego ==  1) fLego->LegoCylindrical(1,nx,ny,"FB");
      if (Hoption.Lego == 11) fLego->LegoCylindrical(1,nx,ny,"BF");
      if (Hoption.Lego == 12) fLego->LegoCylindrical(1,nx,ny,"BF");
   } else if (Hoption.System == kSPHERICAL) {
      if (Hoption.Lego ==  1) fLego->LegoSpherical(0,1,nx,ny,"FB");
      if (Hoption.Lego == 11) fLego->LegoSpherical(0,1,nx,ny,"BF");
      if (Hoption.Lego == 12) fLego->LegoSpherical(0,1,nx,ny,"BF");
   } else if (Hoption.System == kRAPIDITY) {
      if (Hoption.Lego ==  1) fLego->LegoSpherical(1,1,nx,ny,"FB");
      if (Hoption.Lego == 11) fLego->LegoSpherical(1,1,nx,ny,"BF");
      if (Hoption.Lego == 12) fLego->LegoSpherical(1,1,nx,ny,"BF");
   } else {
      if (Hoption.Lego ==  1) {
                              fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMove2);
                              fLego->LegoCartesian(90,nx,ny,"FB");}
      if (Hoption.Lego == 11) fLego->LegoCartesian(90,nx,ny,"BF");
      if (Hoption.Lego == 12) fLego->LegoCartesian(90,nx,ny,"BF");
   }

   if (Hoption.Lego == 1 || Hoption.Lego == 11) {
      if (Hoption.System == kCARTESIAN && Hoption.BackBox) {
         fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMove1);
         fLego->BackBox(90);
      }
   }
   if (Hoption.System == kCARTESIAN) {
      fLego->InitMoveScreen(-1.1,1.1);
      fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMove2);
      if (Hoption.FrontBox) fLego->FrontBox(90);
   }
   if (!Hoption.Axis && !Hoption.Same) PaintLegoAxis(axis, 90);
   if (Hoption.Zscale) PaintPalette();
   delete axis;
   delete fLego; fLego = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the axis for legos and surface plots.

void THistPainter::PaintLegoAxis(TGaxis *axis, Double_t ang)
{

   static Double_t epsil = 0.001;

   Double_t cosa, sina;
   Double_t bmin, bmax;
   Double_t r[24]        /* was [3][8] */;
   Int_t ndivx, ndivy, ndivz, i;
   Double_t x1[3], x2[3], y1[3], y2[3], z1[3], z2[3], av[24]  /*  was [3][8] */;
   static char chopax[8], chopay[8], chopaz[8];
   Int_t ix1, ix2, iy1, iy2, iz1, iz2;
   Double_t rad;

   TView *view = gPad->GetView();
   if (!view) {
      Error("PaintLegoAxis", "no TView in current pad");
      return;
   }

   // In polar coordinates, draw a short line going from the external circle
   // corresponding to r = 1 up to r = 1.1
   if (Hoption.System == kPOLAR) {
      r[0] = 1;
      r[1] = 0;
      r[2] = 0;
      view->WCtoNDC(r, x1);
      r[0] = 1.1;
      r[1] = 0;
      r[2] = 0;
      view->WCtoNDC(r, x2);
      gPad->PaintLine(x1[0],x1[1],x2[0],x2[1]);
      return;
   }

   if (Hoption.System != kCARTESIAN) return;

   rad = TMath::ATan(1.) * 4. /180.;
   cosa = TMath::Cos(ang*rad);
   sina = TMath::Sin(ang*rad);

   view->AxisVertex(ang, av, ix1, ix2, iy1, iy2, iz1, iz2);
   for (i = 1; i <= 8; ++i) {
      r[i*3 - 3] = av[i*3 - 3] + av[i*3 - 2]*cosa;
      r[i*3 - 2] = av[i*3 - 2]*sina;
      r[i*3 - 1] = av[i*3 - 1];
   }

   view->WCtoNDC(&r[ix1*3 - 3], x1);
   view->WCtoNDC(&r[ix2*3 - 3], x2);
   view->WCtoNDC(&r[iy1*3 - 3], y1);
   view->WCtoNDC(&r[iy2*3 - 3], y2);
   view->WCtoNDC(&r[iz1*3 - 3], z1);
   view->WCtoNDC(&r[iz2*3 - 3], z2);

   view->SetAxisNDC(x1, x2, y1, y2, z1, z2);

   Double_t *rmin = view->GetRmin();
   Double_t *rmax = view->GetRmax();
   if (!rmin || !rmax) return;

   // Initialize the axis options
   if (x1[0] > x2[0]) strlcpy(chopax, "SDH=+",8);
   else               strlcpy(chopax, "SDH=-",8);
   if (y1[0] > y2[0]) strlcpy(chopay, "SDH=+",8);
   else               strlcpy(chopay, "SDH=-",8);
   if (z2[1] > z1[1]) strlcpy(chopaz, "SDH=+",8);
   else               strlcpy(chopaz, "SDH=-",8);

   // Option LOG is required ?
   if (Hoption.Logx) strlcat(chopax,"G",8);
   if (Hoption.Logy) strlcat(chopay,"G",8);
   if (Hoption.Logz) strlcat(chopaz,"G",8);

   // Initialize the number of divisions. If the
   // number of divisions is negative, option 'N' is required.
   ndivx = fXaxis->GetNdivisions();
   ndivy = fYaxis->GetNdivisions();
   ndivz = fZaxis->GetNdivisions();
   if (ndivx < 0) {
      ndivx = TMath::Abs(ndivx);
      strlcat(chopax, "N",8);
   }
   if (ndivy < 0) {
      ndivy = TMath::Abs(ndivy);
      strlcat(chopay, "N",8);
   }
   if (ndivz < 0) {
      ndivz = TMath::Abs(ndivz);
      strlcat(chopaz, "N",8);
   }

   // Set Axis attributes.
   // The variable SCALE  rescales the VSIZ
   // in order to have the same label size for all angles.

   axis->SetLineWidth(1);

   // X axis drawing
   if (TMath::Abs(x1[0] - x2[0]) >= epsil || TMath::Abs(x1[1] - x2[1]) > epsil) {
      axis->ImportAxisAttributes(fXaxis);
      axis->SetLabelOffset(fXaxis->GetLabelOffset()+fXaxis->GetTickLength());
      if (Hoption.Logx && !fH->InheritsFrom(TH3::Class())) {
         bmin = TMath::Power(10, rmin[0]);
         bmax = TMath::Power(10, rmax[0]);
      } else {
         bmin = rmin[0];
         bmax = rmax[0];
      }
      // Option time display is required ?
      if (fXaxis->GetTimeDisplay()) {
         strlcat(chopax,"t",8);
         if (strlen(fXaxis->GetTimeFormatOnly()) == 0) {
            axis->SetTimeFormat(fXaxis->ChooseTimeFormat(bmax-bmin));
         } else {
            axis->SetTimeFormat(fXaxis->GetTimeFormat());
         }
      }
      axis->SetOption(chopax);
      axis->PaintAxis(x1[0], x1[1], x2[0], x2[1], bmin, bmax, ndivx, chopax);
   }

   // Y axis drawing
   if (TMath::Abs(y1[0] - y2[0]) >= epsil || TMath::Abs(y1[1] - y2[1]) > epsil) {
      axis->ImportAxisAttributes(fYaxis);
      axis->SetLabelOffset(fYaxis->GetLabelOffset()+fYaxis->GetTickLength());
      if (fYaxis->GetTitleOffset() == 0) axis->SetTitleOffset(1.5);

      if (fH->GetDimension() < 2) {
         strlcpy(chopay, "V=+UN",8);
         ndivy = 0;
      }
      if (TMath::Abs(y1[0] - y2[0]) < epsil) {
         y2[0] = y1[0];
      }
      if (Hoption.Logy && !fH->InheritsFrom(TH3::Class())) {
         bmin = TMath::Power(10, rmin[1]);
         bmax = TMath::Power(10, rmax[1]);
      } else {
         bmin = rmin[1];
         bmax = rmax[1];
      }
      // Option time display is required ?
      if (fYaxis->GetTimeDisplay()) {
         strlcat(chopay,"t",8);
         if (strlen(fYaxis->GetTimeFormatOnly()) == 0) {
            axis->SetTimeFormat(fYaxis->ChooseTimeFormat(bmax-bmin));
         } else {
            axis->SetTimeFormat(fYaxis->GetTimeFormat());
         }
      }
      axis->SetOption(chopay);
      axis->PaintAxis(y1[0], y1[1], y2[0], y2[1], bmin, bmax, ndivy, chopay);
   }

   // Z axis drawing
   if (TMath::Abs(z1[0] - z2[0]) >= 100*epsil || TMath::Abs(z1[1] - z2[1]) > 100*epsil) {
      axis->ImportAxisAttributes(fZaxis);
      if (Hoption.Logz && !fH->InheritsFrom(TH3::Class())) {
         bmin = TMath::Power(10, rmin[2]);
         bmax = TMath::Power(10, rmax[2]);
      } else {
         bmin = rmin[2];
         bmax = rmax[2];
      }
      // Option time display is required ?
      if (fZaxis->GetTimeDisplay()) {
         strlcat(chopaz,"t",8);
         if (strlen(fZaxis->GetTimeFormatOnly()) == 0) {
            axis->SetTimeFormat(fZaxis->ChooseTimeFormat(bmax-bmin));
         } else {
            axis->SetTimeFormat(fZaxis->GetTimeFormat());
         }
      }
      axis->SetOption(chopaz);
      axis->PaintAxis(z1[0], z1[1], z2[0], z2[1], bmin, bmax, ndivz, chopaz);
   }

   //fH->SetLineStyle(1);  /// otherwise fEdgeStyle[i] gets overwritten!
}

////////////////////////////////////////////////////////////////////////////////
/// [Paint the color palette on the right side of the pad.](#HP22)

void THistPainter::PaintPalette()
{

   TPaletteAxis *palette = (TPaletteAxis*)fFunctions->FindObject("palette");
   TView *view = gPad->GetView();
   if (palette) {
      if (view) {
         if (!palette->TestBit(TPaletteAxis::kHasView)) {
            fFunctions->Remove(palette);
            delete palette; palette = 0;
         }
      } else {
         if (palette->TestBit(TPaletteAxis::kHasView)) {
            fFunctions->Remove(palette);
            delete palette; palette = 0;
         }
      }
      // make sure the histogram member of the palette is setup correctly. It may not be after a Clone()
      if (palette && !palette->GetHistogram()) palette->SetHistogram(fH);
   }

   if (!palette) {
      Double_t xup  = gPad->GetUxmax();
      Double_t x2   = gPad->PadtoX(gPad->GetX2());
      Double_t ymin = gPad->PadtoY(gPad->GetUymin());
      Double_t ymax = gPad->PadtoY(gPad->GetUymax());
      Double_t xr   = 0.05*(gPad->GetX2() - gPad->GetX1());
      Double_t xmin = gPad->PadtoX(xup +0.1*xr);
      Double_t xmax = gPad->PadtoX(xup + xr);
      if (xmax > x2) xmax = gPad->PadtoX(gPad->GetX2()-0.01*xr);
      palette = new TPaletteAxis(xmin,ymin,xmax,ymax,fH);
      fFunctions->AddFirst(palette);
      palette->Paint();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// [Control function to draw a 2D histogram as a scatter plot.](#HP11)

void THistPainter::PaintScatterPlot(Option_t *option)
{

   fH->TAttMarker::Modify();

   Int_t k, marker;
   Double_t dz, z, xk,xstep, yk, ystep;
   Double_t scale = 1;
   Bool_t ltest  = kFALSE;
   Double_t zmax  = fH->GetMaximum();
   Double_t zmin  = fH->GetMinimum();
   if (zmin == 0 && zmax == 0) return;
   if (zmin == zmax) {
      zmax += 0.1*TMath::Abs(zmax);
      zmin -= 0.1*TMath::Abs(zmin);
   }
   Int_t ncells = (Hparam.ylast-Hparam.yfirst)*(Hparam.xlast-Hparam.xfirst);
   if (Hoption.Logz) {
      if (zmin > 0) zmin = TMath::Log10(zmin);
      else          zmin = 0;
      if (zmax > 0) zmax = TMath::Log10(zmax);
      else          zmax = 0;
      if (zmin == 0 && zmax == 0) return;
      dz = zmax - zmin;
      scale = 100/dz;
      if (ncells > 10000) scale /= 5;
      ltest = kTRUE;
   } else {
      dz = zmax - zmin;
      if (dz >= kNMAX || zmax < 1) {
         scale = (kNMAX-1)/dz;
         if (ncells > 10000) scale /= 5;
         ltest = kTRUE;
      }
   }
   if (fH->GetMinimumStored() == -1111) {
      Double_t yMARGIN = gStyle->GetHistTopMargin();
      if (Hoption.MinimumZero) {
         if (zmin >= 0) zmin = 0;
         else           zmin -= yMARGIN*(zmax-zmin);
      } else {
         Double_t dzmin = yMARGIN*(zmax-zmin);
         if (zmin >= 0 && (zmin-dzmin <= 0)) zmin  = 0;
         else                                zmin -= dzmin;
      }
   }

   TString opt = option;
   opt.ToLower();
   if (opt.Contains("scat=")) {
      char optscat[100];
      strlcpy(optscat,opt.Data(),100);
      char *oscat = strstr(optscat,"scat=");
      char *blank = strstr(oscat," "); if (blank) *blank = 0;
      sscanf(oscat+5,"%lg",&scale);
   }
   // use an independent instance of a random generator
   // instead of gRandom to avoid conflicts and
   // to get same random numbers when drawing the same histogram
   TRandom2 random;
   marker=0;
   for (Int_t j=Hparam.yfirst; j<=Hparam.ylast;j++) {
      yk    = fYaxis->GetBinLowEdge(j);
      ystep = fYaxis->GetBinWidth(j);
      for (Int_t i=Hparam.xfirst; i<=Hparam.xlast;i++) {
         Int_t bin  = j*(fXaxis->GetNbins()+2) + i;
         xk    = fXaxis->GetBinLowEdge(i);
         xstep = fXaxis->GetBinWidth(i);
         if (!IsInside(xk+0.5*xstep,yk+0.5*ystep)) continue;
         z     = fH->GetBinContent(bin);
         if (z < zmin) z = zmin;
         if (z > zmax) z = zmax;
         if (Hoption.Logz) {
            if (z > 0) z = TMath::Log10(z) - zmin;
         } else {
            z    -=  zmin;
         }
         if (z <= 0) continue;
         k = Int_t(z*scale);
         if (ltest) k++;
         if (k > 0) {
            for (Int_t loop=0; loop<k; loop++) {
               if (k+marker >= kNMAX) {
                  gPad->PaintPolyMarker(marker, fXbuf, fYbuf);
                  marker=0;
               }
               fXbuf[marker] = (random.Rndm()*xstep) + xk;
               fYbuf[marker] = (random.Rndm()*ystep) + yk;
               if (Hoption.Logx) {
                  if (fXbuf[marker] > 0) fXbuf[marker] = TMath::Log10(fXbuf[marker]);
                  else                   break;
               }
               if (Hoption.Logy) {
                  if (fYbuf[marker] > 0) fYbuf[marker] = TMath::Log10(fYbuf[marker]);
                  else                  break;
               }
               if (fXbuf[marker] < gPad->GetUxmin()) break;
               if (fYbuf[marker] < gPad->GetUymin()) break;
               if (fXbuf[marker] > gPad->GetUxmax()) break;
               if (fYbuf[marker] > gPad->GetUymax()) break;
               marker++;
            }
         }
      }
   }
   if (marker > 0) gPad->PaintPolyMarker(marker, fXbuf, fYbuf);

   if (Hoption.Zscale) PaintPalette();
}

////////////////////////////////////////////////////////////////////////////////
/// Static function to paint special objects like vectors and matrices.
/// This function is called via `gROOT->ProcessLine` to paint these objects
/// without having a direct dependency of the graphics or histogramming
/// system.

void THistPainter::PaintSpecialObjects(const TObject *obj, Option_t *option)
{

   if (!obj) return;
   Bool_t status = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   if (obj->InheritsFrom(TMatrixFBase::Class())) {
      // case TMatrixF
      TH2F *R__TMatrixFBase = new TH2F((TMatrixFBase &)*obj);
      R__TMatrixFBase->SetBit(kCanDelete);
      R__TMatrixFBase->Draw(option);

   } else if (obj->InheritsFrom(TMatrixDBase::Class())) {
      // case TMatrixD
      TH2D *R__TMatrixDBase = new TH2D((TMatrixDBase &)*obj);
      R__TMatrixDBase->SetBit(kCanDelete);
      R__TMatrixDBase->Draw(option);

   } else if (obj->InheritsFrom(TVectorF::Class())) {
      //case TVectorF
      TH1F *R__TVectorF = new TH1F((TVectorF &)*obj);
      R__TVectorF->SetBit(kCanDelete);
      R__TVectorF->Draw(option);

   } else if (obj->InheritsFrom(TVectorD::Class())) {
      //case TVectorD
      TH1D *R__TVectorD = new TH1D((TVectorD &)*obj);
      R__TVectorD->SetBit(kCanDelete);
      R__TVectorD->Draw(option);
   }

   TH1::AddDirectory(status);
}

////////////////////////////////////////////////////////////////////////////////
/// [Draw the statistics box for 1D and profile histograms.](#HP07)

void THistPainter::PaintStat(Int_t dostat, TF1 *fit)
{

   TString tt, tf;
   Int_t dofit;
   TPaveStats *stats  = 0;
   TIter next(fFunctions);
   TObject *obj;
   while ((obj = next())) {
      if (obj->InheritsFrom(TPaveStats::Class())) {
         stats = (TPaveStats*)obj;
         break;
      }
   }

   if (stats && dostat) {
      dofit  = stats->GetOptFit();
      dostat = stats->GetOptStat();
   } else {
      dofit  = gStyle->GetOptFit();
   }
   if (!dofit) fit = 0;
   if (dofit  == 1) dofit  =  111;
   if (dostat == 1) dostat = 1111;
   Int_t print_name    = dostat%10;
   Int_t print_entries = (dostat/10)%10;
   Int_t print_mean    = (dostat/100)%10;
   Int_t print_stddev  = (dostat/1000)%10;
   Int_t print_under   = (dostat/10000)%10;
   Int_t print_over    = (dostat/100000)%10;
   Int_t print_integral= (dostat/1000000)%10;
   Int_t print_skew    = (dostat/10000000)%10;
   Int_t print_kurt    = (dostat/100000000)%10;
   Int_t nlines = print_name + print_entries + print_mean + print_stddev +
                  print_under + print_over + print_integral +
                  print_skew + print_kurt;
   Int_t print_fval    = dofit%10;
   Int_t print_ferrors = (dofit/10)%10;
   Int_t print_fchi2   = (dofit/100)%10;
   Int_t print_fprob   = (dofit/1000)%10;
   Int_t nlinesf = print_fval + print_fchi2 + print_fprob;
   if (fit) {
      if (print_fval < 2) nlinesf += fit->GetNumberFreeParameters();
      else                nlinesf += fit->GetNpar();
   }
   if (fH->InheritsFrom(TProfile::Class())) nlinesf += print_mean + print_stddev;

   // Pavetext with statistics
   Bool_t done = kFALSE;
   if (!dostat && !fit) {
      if (stats) { fFunctions->Remove(stats); delete stats;}
      return;
   }
   Double_t  statw  = gStyle->GetStatW();
   if (fit) statw   = 1.8*gStyle->GetStatW();
   Double_t  stath  = (nlines+nlinesf)*gStyle->GetStatFontSize();
   if (stath <= 0 || 3 == (gStyle->GetStatFont()%10)) {
      stath = 0.25*(nlines+nlinesf)*gStyle->GetStatH();
   }
   if (stats) {
      stats->Clear();
      done = kTRUE;
   } else {
      stats  = new TPaveStats(
               gStyle->GetStatX()-statw,
               gStyle->GetStatY()-stath,
               gStyle->GetStatX(),
               gStyle->GetStatY(),"brNDC");

      stats->SetParent(fH);
      stats->SetOptFit(dofit);
      stats->SetOptStat(dostat);
      stats->SetFillColor(gStyle->GetStatColor());
      stats->SetFillStyle(gStyle->GetStatStyle());
      stats->SetBorderSize(gStyle->GetStatBorderSize());
      stats->SetTextFont(gStyle->GetStatFont());
      if (gStyle->GetStatFont()%10 > 2)
         stats->SetTextSize(gStyle->GetStatFontSize());
      stats->SetFitFormat(gStyle->GetFitFormat());
      stats->SetStatFormat(gStyle->GetStatFormat());
      stats->SetName("stats");

      stats->SetTextColor(gStyle->GetStatTextColor());
      stats->SetTextAlign(12);
      stats->SetBit(kCanDelete);
      stats->SetBit(kMustCleanup);
   }
   if (print_name)  stats->AddText(fH->GetName());
   if (print_entries) {
      if (fH->GetEntries() < 1e7) tt.Form("%s = %-7d",gStringEntries.Data(),Int_t(fH->GetEntries()+0.5));
      else                        tt.Form("%s = %14.7g",gStringEntries.Data(),Float_t(fH->GetEntries()));
      stats->AddText(tt.Data());
   }
   if (print_mean) {
      if (print_mean == 1) {
         tf.Form("%s  = %s%s",gStringMean.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),fH->GetMean(1));
      } else {
         tf.Form("%s  = %s%s #pm %s%s",gStringMean.Data(),"%",stats->GetStatFormat()
                                                  ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),fH->GetMean(1),fH->GetMeanError(1));
      }
      stats->AddText(tt.Data());
      if (fH->InheritsFrom(TProfile::Class())) {
         if (print_mean == 1) {
            tf.Form("%s = %s%s",gStringMeanY.Data(),"%",stats->GetStatFormat());
            tt.Form(tf.Data(),fH->GetMean(2));
         } else {
            tf.Form("%s = %s%s #pm %s%s",gStringMeanY.Data(),"%",stats->GetStatFormat()
                                                      ,"%",stats->GetStatFormat());
            tt.Form(tf.Data(),fH->GetMean(2),fH->GetMeanError(2));
         }
         stats->AddText(tt.Data());
      }
   }
   if (print_stddev) {
      if (print_stddev == 1) {
         tf.Form("%s   = %s%s",gStringStdDev.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),fH->GetStdDev(1));
      } else {
         tf.Form("%s   = %s%s #pm %s%s",gStringStdDev.Data(),"%",stats->GetStatFormat()
                                                  ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),fH->GetStdDev(1),fH->GetStdDevError(1));
      }
      stats->AddText(tt.Data());
      if (fH->InheritsFrom(TProfile::Class())) {
         if (print_stddev == 1) {
            tf.Form("%s = %s%s",gStringStdDevY.Data(),"%",stats->GetStatFormat());
            tt.Form(tf.Data(),fH->GetStdDev(2));
         } else {
            tf.Form("%s = %s%s #pm %s%s",gStringStdDevY.Data(),"%",stats->GetStatFormat()
                                                     ,"%",stats->GetStatFormat());
            tt.Form(tf.Data(),fH->GetStdDev(2),fH->GetStdDevError(2));
         }
         stats->AddText(tt.Data());
      }
   }
   if (print_under) {
      tf.Form("%s = %s%s",gStringUnderflow.Data(),"%",stats->GetStatFormat());
      tt.Form(tf.Data(),fH->GetBinContent(0));
      stats->AddText(tt.Data());
   }
   if (print_over) {
      tf.Form("%s  = %s%s",gStringOverflow.Data(),"%",stats->GetStatFormat());
      tt.Form(tf.Data(),fH->GetBinContent(fXaxis->GetNbins()+1));
      stats->AddText(tt.Data());
   }
   if (print_integral) {
      if (print_integral == 1) {
         tf.Form("%s = %s%s",gStringIntegral.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),fH->Integral());
      } else {
         tf.Form("%s = %s%s",gStringIntegralBinWidth.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),fH->Integral("width"));
      }
      stats->AddText(tt.Data());
   }
   if (print_skew) {
      if (print_skew == 1) {
         tf.Form("%s = %s%s",gStringSkewness.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),fH->GetSkewness(1));
      } else {
         tf.Form("%s = %s%s #pm %s%s",gStringSkewness.Data(),"%",stats->GetStatFormat()
                                                     ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),fH->GetSkewness(1),fH->GetSkewness(11));
      }
      stats->AddText(tt.Data());
   }
   if (print_kurt) {
      if (print_kurt == 1) {
         tf.Form("%s = %s%s",gStringKurtosis.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),fH->GetKurtosis(1));
      } else {
         tf.Form("%s = %s%s #pm %s%s",gStringKurtosis.Data(),"%",stats->GetStatFormat()
                                                     ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),fH->GetKurtosis(1),fH->GetKurtosis(11));
      }
      stats->AddText(tt.Data());
   }

   // Draw Fit parameters
   if (fit) {
      Int_t ndf = fit->GetNDF();
      tf.Form("#chi^{2} / ndf = %s%s / %d","%",stats->GetFitFormat(),ndf);
      tt.Form(tf.Data(),(Float_t)fit->GetChisquare());
      if (print_fchi2) stats->AddText(tt.Data());
      if (print_fprob) {
         tf.Form("Prob  = %s%s","%",stats->GetFitFormat());
         tt.Form(tf.Data(),(Float_t)TMath::Prob(fit->GetChisquare(),ndf));
         stats->AddText(tt.Data());
      }
      if (print_fval || print_ferrors) {
         Double_t parmin,parmax;
         for (Int_t ipar=0;ipar<fit->GetNpar();ipar++) {
            fit->GetParLimits(ipar,parmin,parmax);
            if (print_fval < 2 && parmin*parmax != 0 && parmin >= parmax) continue;
            if (print_ferrors) {
               tf.Form("%-8s = %s%s #pm %s ", fit->GetParName(ipar), "%",stats->GetFitFormat(),
                       GetBestFormat(fit->GetParameter(ipar), fit->GetParError(ipar), stats->GetFitFormat()));
               tt.Form(tf.Data(),(Float_t)fit->GetParameter(ipar)
                               ,(Float_t)fit->GetParError(ipar));
            } else {
               tf.Form("%-8s = %s%s ",fit->GetParName(ipar), "%",stats->GetFitFormat());
               tt.Form(tf.Data(),(Float_t)fit->GetParameter(ipar));
            }
            stats->AddText(tt.Data());
         }
      }
   }

   if (!done) fFunctions->Add(stats);
   stats->Paint();
}

////////////////////////////////////////////////////////////////////////////////
/// [Draw the statistics box for 2D histograms.](#HP07)

void THistPainter::PaintStat2(Int_t dostat, TF1 *fit)
{

   if (fH->GetDimension() != 2) return;
   TH2 *h2 = (TH2*)fH;

   TString tt, tf;
   Int_t dofit;
   TPaveStats *stats  = 0;
   TIter next(fFunctions);
   TObject *obj;
   while ((obj = next())) {
      if (obj->InheritsFrom(TPaveStats::Class())) {
         stats = (TPaveStats*)obj;
         break;
      }
   }
   if (stats && dostat) {
      dofit  = stats->GetOptFit();
      dostat = stats->GetOptStat();
   } else {
      dofit  = gStyle->GetOptFit();
   }
   if (dostat == 1) dostat = 1111;
   Int_t print_name    = dostat%10;
   Int_t print_entries = (dostat/10)%10;
   Int_t print_mean    = (dostat/100)%10;
   Int_t print_stddev  = (dostat/1000)%10;
   Int_t print_under   = (dostat/10000)%10;
   Int_t print_over    = (dostat/100000)%10;
   Int_t print_integral= (dostat/1000000)%10;
   Int_t print_skew    = (dostat/10000000)%10;
   Int_t print_kurt    = (dostat/100000000)%10;
   Int_t nlines = print_name + print_entries + 2*print_mean + 2*print_stddev + print_integral;
   if (print_under || print_over) nlines += 3;

   // Pavetext with statistics
   if (!gStyle->GetOptFit()) fit = 0;
   Bool_t done = kFALSE;
   if (!dostat && !fit) {
      if (stats) { fFunctions->Remove(stats); delete stats;}
      return;
   }
   Double_t  statw  = gStyle->GetStatW();
   if (fit) statw   = 1.8*gStyle->GetStatW();
   Double_t  stath  = nlines*gStyle->GetStatFontSize();
   if (stath <= 0 || 3 == (gStyle->GetStatFont()%10)) {
      stath = 0.25*nlines*gStyle->GetStatH();
   }
   if (fit) stath += gStyle->GetStatH();
   if (stats) {
      stats->Clear();
      done = kTRUE;
   } else {
      stats  = new TPaveStats(
               gStyle->GetStatX()-statw,
               gStyle->GetStatY()-stath,
               gStyle->GetStatX(),
               gStyle->GetStatY(),"brNDC");

      stats->SetParent(fH);
      stats->SetOptFit(dofit);
      stats->SetOptStat(dostat);
      stats->SetFillColor(gStyle->GetStatColor());
      stats->SetFillStyle(gStyle->GetStatStyle());
      stats->SetBorderSize(gStyle->GetStatBorderSize());
      stats->SetName("stats");

      stats->SetTextColor(gStyle->GetStatTextColor());
      stats->SetTextAlign(12);
      stats->SetTextFont(gStyle->GetStatFont());
      if (gStyle->GetStatFont()%10 > 2)
         stats->SetTextSize(gStyle->GetStatFontSize());
      stats->SetFitFormat(gStyle->GetFitFormat());
      stats->SetStatFormat(gStyle->GetStatFormat());
      stats->SetBit(kCanDelete);
      stats->SetBit(kMustCleanup);
   }
   if (print_name)  stats->AddText(h2->GetName());
   if (print_entries) {
      if (h2->GetEntries() < 1e7) tt.Form("%s = %-7d",gStringEntries.Data(),Int_t(h2->GetEntries()+0.5));
      else                        tt.Form("%s = %14.7g",gStringEntries.Data(),Float_t(h2->GetEntries()));
      stats->AddText(tt.Data());
   }
   if (print_mean) {
      if (print_mean == 1) {
         tf.Form("%s = %s%s",gStringMeanX.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h2->GetMean(1));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s",gStringMeanY.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h2->GetMean(2));
         stats->AddText(tt.Data());
      } else {
         tf.Form("%s = %s%s #pm %s%s",gStringMeanX.Data(),"%",stats->GetStatFormat()
                                                   ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h2->GetMean(1),h2->GetMeanError(1));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s #pm %s%s",gStringMeanY.Data(),"%",stats->GetStatFormat()
                                                   ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h2->GetMean(2),h2->GetMeanError(2));
         stats->AddText(tt.Data());
      }
   }
   if (print_stddev) {
      if (print_stddev == 1) {
         tf.Form("%s = %s%s",gStringStdDevX.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h2->GetStdDev(1));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s",gStringStdDevY.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h2->GetStdDev(2));
         stats->AddText(tt.Data());
      } else {
         tf.Form("%s = %s%s #pm %s%s",gStringStdDevX.Data(),"%",stats->GetStatFormat()
                                                  ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h2->GetStdDev(1),h2->GetStdDevError(1));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s #pm %s%s",gStringStdDevY.Data(),"%",stats->GetStatFormat()
                                                  ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h2->GetStdDev(2),h2->GetStdDevError(2));
         stats->AddText(tt.Data());
      }
   }
   if (print_integral) {
      tf.Form("%s = %s%s",gStringIntegral.Data(),"%",stats->GetStatFormat());
      tt.Form(tf.Data(),fH->Integral());
      stats->AddText(tt.Data());
   }
   if (print_skew) {
      if (print_skew == 1) {
         tf.Form("%s = %s%s",gStringSkewnessX.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h2->GetSkewness(1));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s",gStringSkewnessY.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h2->GetSkewness(2));
         stats->AddText(tt.Data());
      } else {
         tf.Form("%s = %s%s #pm %s%s",gStringSkewnessX.Data(),"%",stats->GetStatFormat()
                                                       ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h2->GetSkewness(1),h2->GetSkewness(11));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s #pm %s%s",gStringSkewnessY.Data(),"%",stats->GetStatFormat()
                                                       ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h2->GetSkewness(2),h2->GetSkewness(12));
         stats->AddText(tt.Data());
      }
   }
   if (print_kurt) {
      if (print_kurt == 1) {
         tf.Form("%s = %s%s",gStringKurtosisX.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h2->GetKurtosis(1));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s",gStringKurtosisY.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h2->GetKurtosis(2));
         stats->AddText(tt.Data());
      } else {
         tf.Form("%s = %s%s #pm %s%s",gStringKurtosisX.Data(),"%",stats->GetStatFormat()
                                                       ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h2->GetKurtosis(1),h2->GetKurtosis(11));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s #pm %s%s",gStringKurtosisY.Data(),"%",stats->GetStatFormat()
                                                       ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h2->GetKurtosis(2),h2->GetKurtosis(12));
         stats->AddText(tt.Data());
      }
   }
   if (print_under || print_over) {
      //get 3*3 under/overflows for 2d hist
      Double_t unov[9];

      Int_t cellsX = h2->GetXaxis()->GetNbins() + 1;
      Int_t cellsY = h2->GetYaxis()->GetNbins() + 1;
      Int_t firstX = std::max(1, h2->GetXaxis()->GetFirst());
      Int_t firstY = std::max(1, h2->GetYaxis()->GetFirst());
      Int_t lastX  = std::min(h2->GetXaxis()->GetLast(), h2->GetXaxis()->GetNbins());
      Int_t lastY  = std::min(h2->GetYaxis()->GetLast(), h2->GetYaxis()->GetNbins());

      unov[0] = h2->Integral(      0, firstX-1, lastY+1, cellsY  );
      unov[1] = h2->Integral(firstX , lastX   , lastY+1, cellsY  );
      unov[2] = h2->Integral(lastX+1, cellsX  , lastY+1, cellsY  );
      unov[3] = h2->Integral(      0, firstX-1, firstY , lastY   );
      unov[4] = h2->Integral(firstX , lastX   , firstY , lastY   );
      unov[5] = h2->Integral(lastX+1, cellsX  , firstY , lastY   );
      unov[6] = h2->Integral(      0, firstX-1,       0, firstY-1);
      unov[7] = h2->Integral(firstX, lastX,           0, firstY-1);
      unov[8] = h2->Integral(lastX+1, cellsX  ,       0, firstY-1);

      tt.Form(" %7d|%7d|%7d\n", (Int_t)unov[0], (Int_t)unov[1], (Int_t)unov[2]);
      stats->AddText(tt.Data());
      if (TMath::Abs(unov[4]) < 1.e7)
         tt.Form(" %7d|%7d|%7d\n", (Int_t)unov[3], (Int_t)unov[4], (Int_t)unov[5]);
      else
         tt.Form(" %7d|%14.7g|%7d\n", (Int_t)unov[3], (Float_t)unov[4], (Int_t)unov[5]);
      stats->AddText(tt.Data());
      tt.Form(" %7d|%7d|%7d\n", (Int_t)unov[6], (Int_t)unov[7], (Int_t)unov[8]);
      stats->AddText(tt.Data());
   }

   // Draw Fit parameters
   if (fit) {
      Int_t ndf = fit->GetNDF();
      tt.Form("#chi^{2} / ndf = %6.4g / %d",(Float_t)fit->GetChisquare(),ndf);
      stats->AddText(tt.Data());
      for (Int_t ipar=0;ipar<fit->GetNpar();ipar++) {
         tt.Form("%-8s = %5.4g #pm %5.4g ",fit->GetParName(ipar)
                                   ,(Float_t)fit->GetParameter(ipar)
                                   ,(Float_t)fit->GetParError(ipar));
         stats->AddText(tt.Data());
      }
   }

   if (!done) fFunctions->Add(stats);
   stats->Paint();
}

////////////////////////////////////////////////////////////////////////////////
/// [Draw the statistics box for 3D histograms.](#HP07)

void THistPainter::PaintStat3(Int_t dostat, TF1 *fit)
{

   if (fH->GetDimension() != 3) return;
   TH3 *h3 = (TH3*)fH;

   TString tt, tf;
   Int_t dofit;
   TPaveStats *stats  = 0;
   TIter next(fFunctions);
   TObject *obj;
   while ((obj = next())) {
      if (obj->InheritsFrom(TPaveStats::Class())) {
         stats = (TPaveStats*)obj;
         break;
      }
   }
   if (stats && dostat) {
      dofit  = stats->GetOptFit();
      dostat = stats->GetOptStat();
   } else {
      dofit  = gStyle->GetOptFit();
   }
   if (dostat == 1) dostat = 1111;
   Int_t print_name    = dostat%10;
   Int_t print_entries = (dostat/10)%10;
   Int_t print_mean    = (dostat/100)%10;
   Int_t print_stddev  = (dostat/1000)%10;
   Int_t print_under   = (dostat/10000)%10;
   Int_t print_over    = (dostat/100000)%10;
   Int_t print_integral= (dostat/1000000)%10;
   Int_t print_skew    = (dostat/10000000)%10;
   Int_t print_kurt    = (dostat/100000000)%10;
   Int_t nlines = print_name + print_entries + 3*print_mean + 3*print_stddev + print_integral;
   if (print_under || print_over) nlines += 3;

   // Pavetext with statistics
   if (!gStyle->GetOptFit()) fit = 0;
   Bool_t done = kFALSE;
   if (!dostat && !fit) {
      if (stats) { fFunctions->Remove(stats); delete stats;}
      return;
   }
   Double_t  statw  = gStyle->GetStatW();
   if (fit) statw   = 1.8*gStyle->GetStatW();
   Double_t  stath  = nlines*gStyle->GetStatFontSize();
   if (stath <= 0 || 3 == (gStyle->GetStatFont()%10)) {
      stath = 0.25*nlines*gStyle->GetStatH();
   }
   if (fit) stath += gStyle->GetStatH();
   if (stats) {
      stats->Clear();
      done = kTRUE;
   } else {
      stats  = new TPaveStats(
               gStyle->GetStatX()-statw,
               gStyle->GetStatY()-stath,
               gStyle->GetStatX(),
               gStyle->GetStatY(),"brNDC");

      stats->SetParent(fH);
      stats->SetOptFit(dofit);
      stats->SetOptStat(dostat);
      stats->SetFillColor(gStyle->GetStatColor());
      stats->SetFillStyle(gStyle->GetStatStyle());
      stats->SetBorderSize(gStyle->GetStatBorderSize());
      stats->SetName("stats");

      stats->SetTextColor(gStyle->GetStatTextColor());
      stats->SetTextAlign(12);
      stats->SetTextFont(gStyle->GetStatFont());
      stats->SetFitFormat(gStyle->GetFitFormat());
      stats->SetStatFormat(gStyle->GetStatFormat());
      stats->SetBit(kCanDelete);
      stats->SetBit(kMustCleanup);
   }
   if (print_name)  stats->AddText(h3->GetName());
   if (print_entries) {
      if (h3->GetEntries() < 1e7) tt.Form("%s = %-7d",gStringEntries.Data(),Int_t(h3->GetEntries()+0.5));
      else                        tt.Form("%s = %14.7g",gStringEntries.Data(),Float_t(h3->GetEntries()+0.5));
      stats->AddText(tt.Data());
   }
   if (print_mean) {
      if (print_mean == 1) {
         tf.Form("%s = %s%s",gStringMeanX.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetMean(1));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s",gStringMeanY.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetMean(2));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s",gStringMeanZ.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetMean(3));
         stats->AddText(tt.Data());
      } else {
         tf.Form("%s = %s%s #pm %s%s",gStringMeanX.Data(),"%",stats->GetStatFormat()
                                                   ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetMean(1),h3->GetMeanError(1));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s #pm %s%s",gStringMeanY.Data(),"%",stats->GetStatFormat()
                                                   ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetMean(2),h3->GetMeanError(2));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s #pm %s%s",gStringMeanZ.Data(),"%",stats->GetStatFormat()
                                                   ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetMean(3),h3->GetMeanError(3));
         stats->AddText(tt.Data());
      }
   }
   if (print_stddev) {
      if (print_stddev == 1) {
         tf.Form("%s = %s%s",gStringStdDevX.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetStdDev(1));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s",gStringStdDevY.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetStdDev(2));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s",gStringStdDevZ.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetStdDev(3));
         stats->AddText(tt.Data());
      } else {
         tf.Form("%s = %s%s #pm %s%s",gStringStdDevX.Data(),"%",stats->GetStatFormat()
                                                  ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetStdDev(1),h3->GetStdDevError(1));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s #pm %s%s",gStringStdDevY.Data(),"%",stats->GetStatFormat()
                                                  ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetStdDev(2),h3->GetStdDevError(2));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s #pm %s%s",gStringStdDevZ.Data(),"%",stats->GetStatFormat()
                                                  ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetStdDev(3),h3->GetStdDevError(3));
         stats->AddText(tt.Data());
      }
   }
   if (print_integral) {
      tt.Form("%s  = %6.4g",gStringIntegral.Data(),h3->Integral());
      stats->AddText(tt.Data());
   }
   if (print_skew) {
      if (print_skew == 1) {
         tf.Form("%s = %s%s",gStringSkewnessX.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetSkewness(1));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s",gStringSkewnessY.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetSkewness(2));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s",gStringSkewnessZ.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetSkewness(3));
         stats->AddText(tt.Data());
      } else {
         tf.Form("%s = %s%s #pm %s%s",gStringSkewnessX.Data(),"%",stats->GetStatFormat()
                                                       ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetSkewness(1),h3->GetSkewness(11));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s #pm %s%s",gStringSkewnessY.Data(),"%",stats->GetStatFormat()
                                                       ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetSkewness(2),h3->GetSkewness(12));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s #pm %s%s",gStringSkewnessZ.Data(),"%",stats->GetStatFormat()
                                                       ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetSkewness(3),h3->GetSkewness(13));
         stats->AddText(tt.Data());
      }
   }
   if (print_kurt) {
      if (print_kurt == 1) {
         tf.Form("%s = %s%s",gStringKurtosisX.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetKurtosis(1));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s",gStringKurtosisY.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetKurtosis(2));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s",gStringKurtosisZ.Data(),"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetKurtosis(3));
         stats->AddText(tt.Data());
      } else {
         tf.Form("%s = %s%s #pm %s%s",gStringKurtosisX.Data(),"%",stats->GetStatFormat()
                                                       ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetKurtosis(1),h3->GetKurtosis(11));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s #pm %s%s",gStringKurtosisY.Data(),"%",stats->GetStatFormat()
                                                       ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetKurtosis(2),h3->GetKurtosis(12));
         stats->AddText(tt.Data());
         tf.Form("%s = %s%s #pm %s%s",gStringKurtosisZ.Data(),"%",stats->GetStatFormat()
                                                       ,"%",stats->GetStatFormat());
         tt.Form(tf.Data(),h3->GetKurtosis(3),h3->GetKurtosis(13));
         stats->AddText(tt.Data());
      }
   }
   if (print_under || print_over) {
      // no underflow - overflow printing for a 3D histogram
      // one would need a 3D table
   }

   // Draw Fit parameters
   if (fit) {
      Int_t ndf = fit->GetNDF();
      tt.Form("#chi^{2} / ndf = %6.4g / %d",(Float_t)fit->GetChisquare(),ndf);
      stats->AddText(tt.Data());
      for (Int_t ipar=0;ipar<fit->GetNpar();ipar++) {
         tt.Form("%-8s = %5.4g #pm %5.4g ",fit->GetParName(ipar)
                                   ,(Float_t)fit->GetParameter(ipar)
                                   ,(Float_t)fit->GetParError(ipar));
         stats->AddText(tt.Data());
      }
   }

   if (!done) fFunctions->Add(stats);
   stats->Paint();
}

////////////////////////////////////////////////////////////////////////////////
/// [Control function to draw a 2D histogram as a surface plot.](#HP18)

void THistPainter::PaintSurface(Option_t *)
{

   const Double_t ydiff = 1;
   const Double_t yligh1 = 10;
   const Double_t qa = 0.15;
   const Double_t qd = 0.15;
   const Double_t qs = 0.8;
   Double_t fmin, fmax;
   Int_t raster = 0;
   Int_t irep   = 0;

   if (Hparam.zmin == 0 && Hparam.zmax == 0) {Hparam.zmin = -1; Hparam.zmax = 1;}
   Int_t   nx      = Hparam.xlast - Hparam.xfirst;
   Int_t   ny      = Hparam.ylast - Hparam.yfirst;
   Double_t zmin   = Hparam.zmin;
   Double_t zmax   = Hparam.zmax;
   Double_t xlab1  = Hparam.xmin;
   Double_t xlab2  = Hparam.xmax;
   Double_t ylab1  = Hparam.ymin;
   Double_t ylab2  = Hparam.ymax;
   Double_t dangle = 10*3.141592/180; //Delta angle for Rapidity option
   Double_t deltaz = TMath::Abs(zmin);
   if (deltaz == 0) deltaz = 1;
   if (zmin >= zmax) {
      zmin -= 0.5*deltaz;
      zmax += 0.5*deltaz;
   }
   Double_t z1c = zmin;
   Double_t z2c = zmin + (zmax-zmin)*(1+gStyle->GetHistTopMargin());
   //     Compute the lego limits and instantiate a lego object
   fXbuf[0] = -1;
   fYbuf[0] =  1;
   fXbuf[1] = -1;
   fYbuf[1] =  1;
   if (Hoption.System >= kPOLAR && (Hoption.Surf == 1 || Hoption.Surf == 13)) raster = 1;
   if (Hoption.System == kPOLAR) {
      fXbuf[2] = z1c;
      fYbuf[2] = z2c;
   } else if (Hoption.System == kCYLINDRICAL) {
      if (Hoption.Logy) {
         if (ylab1 > 0) fXbuf[2] = TMath::Log10(ylab1);
         else           fXbuf[2] = 0;
         if (ylab2 > 0) fYbuf[2] = TMath::Log10(ylab2);
         else           fYbuf[2] = 0;
      } else {
         fXbuf[2] = ylab1;
         fYbuf[2] = ylab2;
      }
      z1c = 0; z2c = 1;
   } else if (Hoption.System == kSPHERICAL) {
      fXbuf[2] = -1;
      fYbuf[2] =  1;
      z1c = 0; z2c = 1;
   } else if (Hoption.System == kRAPIDITY) {
      fXbuf[2] = -1/TMath::Tan(dangle);
      fYbuf[2] =  1/TMath::Tan(dangle);
   } else {
      fXbuf[0] = xlab1;
      fYbuf[0] = xlab2;
      fXbuf[1] = ylab1;
      fYbuf[1] = ylab2;
      fXbuf[2] = z1c;
      fYbuf[2] = z2c;
   }

   fLego = new TPainter3dAlgorithms(fXbuf, fYbuf, Hoption.System);
   fLego->SetEdgeAtt(fH->GetLineColor(),fH->GetLineStyle(),fH->GetLineWidth(),0);
   fLego->SetFillColor(fH->GetFillColor());

   //          Create axis object

   TGaxis *axis = new TGaxis();

   //                  Initialize the levels on the Z axis
   Int_t ndiv   = fH->GetContour();
   if (ndiv == 0 ) {
      ndiv = gStyle->GetNumberContours();
      fH->SetContour(ndiv);
   }
   Int_t ndivz  = TMath::Abs(ndiv);
   if (fH->TestBit(TH1::kUserContour) == 0) fH->SetContour(ndiv);

   if (Hoption.Surf == 13 || Hoption.Surf == 15) fLego->SetMesh(3);
   if (Hoption.Surf == 12 || Hoption.Surf == 14 || Hoption.Surf == 17) fLego->SetMesh(0);

   //     Close the surface in case of non cartesian coordinates.

   if (Hoption.System != kCARTESIAN) {nx++; ny++;}

   //     Now ready to draw the surface plot

   TView *view = gPad->GetView();
   if (!view) {
      Error("PaintSurface", "no TView in current pad");
      return;
   }

   Double_t thedeg =  90 - gPad->GetTheta();
   Double_t phideg = -90 - gPad->GetPhi();
   Double_t psideg = view->GetPsi();
   view->SetView(phideg, thedeg, psideg, irep);

   //     Set color/style for back box
   if (Hoption.Same) {
      fLego->SetFillStyle(0);
      fLego->SetFillColor(1);
   } else {
      fLego->SetFillStyle(gPad->GetFrameFillStyle());
      fLego->SetFillColor(gPad->GetFrameFillColor());
   }
   fLego->TAttFill::Modify();

   Int_t backcolor = gPad->GetFrameFillColor();
   if (Hoption.System != kCARTESIAN) backcolor = 0;
   view->PadRange(backcolor);

   fLego->SetFillStyle(fH->GetFillStyle());
   fLego->SetFillColor(fH->GetFillColor());
   fLego->TAttFill::Modify();

   //     Draw the filled contour on top
   Int_t icol1 = fH->GetFillColor();

   Int_t hoption35 = Hoption.Surf;
   if (Hoption.Surf == 13 || Hoption.Surf == 15) {
      DefineColorLevels(ndivz);
      Hoption.Surf = 23;
      fLego->SetSurfaceFunction(&TPainter3dAlgorithms::SurfaceFunction);
      fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMode2);
      if (Hoption.System == kPOLAR)       fLego->SurfacePolar(1,nx,ny,"BF");
      if (Hoption.System == kCYLINDRICAL) fLego->SurfaceCylindrical(1,nx,ny,"BF");
      if (Hoption.System == kSPHERICAL)   fLego->SurfaceSpherical(0,1,nx,ny,"BF");
      if (Hoption.System == kRAPIDITY )   fLego->SurfaceSpherical(1,1,nx,ny,"BF");
      if (Hoption.System == kCARTESIAN)   fLego->SurfaceCartesian(90,nx,ny,"BF");
      Hoption.Surf = hoption35;
      fLego->SetMesh(1);
   }

   if (raster) fLego->InitRaster(-1.1,-1.1,1.1,1.1,1000,800);
   else        fLego->InitMoveScreen(-1.1,1.1);

   if (Hoption.Surf == 11 || Hoption.Surf == 12 || Hoption.Surf == 14 || Hoption.Surf == 17) {
      fLego->DefineGridLevels(fZaxis->GetNdivisions()%100);
      if (Hoption.System == kCARTESIAN && Hoption.BackBox) {
         fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMove1);
         fLego->BackBox(90);
      }
   }

   //     Gouraud Shading surface
   if (Hoption.Surf == 14) {
   //    Set light sources
      fLego->LightSource(0, ydiff, 0,0,0,irep);
      fLego->LightSource(1, yligh1 ,1,1,1,irep);
      fLego->SurfaceProperty(qa, qd, qs, 1, irep);
      fmin = ydiff*qa;
      fmax = fmin + (yligh1+0.1)*(qd+qs);
      Int_t nbcol = 28;
      icol1 = 201;
      Double_t dcol = 0.5/Double_t(nbcol);
      TColor *colref = gROOT->GetColor(fH->GetFillColor());
      if (!colref) return;
      Float_t r,g,b,hue,light,satur;
      colref->GetRGB(r,g,b);
      TColor::RGBtoHLS(r,g,b,hue,light,satur);
      TColor *acol;
      for (Int_t col=0;col<nbcol;col++) {
         acol = gROOT->GetColor(col+icol1);
         TColor::HLStoRGB(hue,.4+col*dcol,satur,r,g,b);
         if (acol) acol->SetRGB(r,g,b);
      }
      fLego->Spectrum(nbcol, fmin, fmax, icol1, 1, irep);
      fLego->SetSurfaceFunction(&TPainter3dAlgorithms::GouraudFunction);
      fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMode2);
      if (Hoption.System == kPOLAR)       fLego->SurfacePolar(1,nx,ny,"BF");
      if (Hoption.System == kCYLINDRICAL) fLego->SurfaceCylindrical(1,nx,ny,"BF");
      if (Hoption.System == kSPHERICAL)   fLego->SurfaceSpherical(0,1,nx,ny,"BF");
      if (Hoption.System == kRAPIDITY )   fLego->SurfaceSpherical(1,1,nx,ny,"BF");
      if (Hoption.System == kCARTESIAN)   fLego->SurfaceCartesian(90,nx,ny,"BF");
   } else if (Hoption.Surf == 15) {
   // The surface is not drawn in this case.
   } else {
   //     Draw the surface
      if (Hoption.Surf == 11 || Hoption.Surf == 12 || Hoption.Surf == 16 || Hoption.Surf == 17) {
         DefineColorLevels(ndivz);
      } else {
         fLego->DefineGridLevels(fZaxis->GetNdivisions()%100);
      }
      fLego->SetSurfaceFunction(&TPainter3dAlgorithms::SurfaceFunction);
      if (Hoption.Surf ==  1 || Hoption.Surf == 13) fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceRaster1);
      if (Hoption.Surf == 11 || Hoption.Surf == 12 || Hoption.Surf == 17) fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMode2);
      if (Hoption.System == kPOLAR) {
         if (Hoption.Surf ==  1 || Hoption.Surf == 13) fLego->SurfacePolar(1,nx,ny,"FB");
         if (Hoption.Surf == 11 || Hoption.Surf == 12 || Hoption.Surf == 17) fLego->SurfacePolar(1,nx,ny,"BF");
      } else if (Hoption.System == kCYLINDRICAL) {
         if (Hoption.Surf ==  1 || Hoption.Surf == 13) fLego->SurfaceCylindrical(1,nx,ny,"FB");
         if (Hoption.Surf == 11 || Hoption.Surf == 12 || Hoption.Surf == 17) fLego->SurfaceCylindrical(1,nx,ny,"BF");
      } else if (Hoption.System == kSPHERICAL) {
         if (Hoption.Surf ==  1 || Hoption.Surf == 13) fLego->SurfaceSpherical(0,1,nx,ny,"FB");
         if (Hoption.Surf == 11 || Hoption.Surf == 12 || Hoption.Surf == 17) fLego->SurfaceSpherical(0,1,nx,ny,"BF");
      } else if (Hoption.System == kRAPIDITY) {
         if (Hoption.Surf ==  1 || Hoption.Surf == 13) fLego->SurfaceSpherical(1,1,nx,ny,"FB");
         if (Hoption.Surf == 11 || Hoption.Surf == 12 || Hoption.Surf == 17) fLego->SurfaceSpherical(1,1,nx,ny,"BF");
      } else {
         if (Hoption.Surf ==  1 || Hoption.Surf == 13) fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMove1);
         if (Hoption.Surf == 16) fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMove3);
         if (Hoption.Surf ==  1 || Hoption.Surf == 13 || Hoption.Surf == 16) fLego->SurfaceCartesian(90,nx,ny,"FB");
         if (Hoption.Surf == 11 || Hoption.Surf == 12 || Hoption.Surf == 17) fLego->SurfaceCartesian(90,nx,ny,"BF");
      }
   }

   // Paint the line contour on top for option SURF7
   if (Hoption.Surf == 17) {
      fLego->InitMoveScreen(-1.1,1.1);
      fLego->DefineGridLevels(fZaxis->GetNdivisions()%100);
      Hoption.Surf = 23;
      fLego->SetSurfaceFunction(&TPainter3dAlgorithms::SurfaceFunction);
      fLego->SetDrawFace(&TPainter3dAlgorithms::DrawLevelLines);
      if (Hoption.System == kPOLAR)       fLego->SurfacePolar(1,nx,ny,"FB");
      if (Hoption.System == kCYLINDRICAL) fLego->SurfaceCylindrical(1,nx,ny,"FB");
      if (Hoption.System == kSPHERICAL)   fLego->SurfaceSpherical(0,1,nx,ny,"FB");
      if (Hoption.System == kRAPIDITY )   fLego->SurfaceSpherical(1,1,nx,ny,"FB");
      if (Hoption.System == kCARTESIAN)   fLego->SurfaceCartesian(90,nx,ny,"FB");
   }

   if ((!Hoption.Same) &&
       (Hoption.Surf == 1 || Hoption.Surf == 13 || Hoption.Surf == 16)) {
      if (Hoption.System == kCARTESIAN && Hoption.BackBox) {
         fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMove1);
         fLego->BackBox(90);
      }
   }
   if (Hoption.System == kCARTESIAN) {
      fLego->InitMoveScreen(-1.1,1.1);
      fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMove2);
      if (Hoption.FrontBox) fLego->FrontBox(90);
   }
   if (!Hoption.Axis && !Hoption.Same) PaintLegoAxis(axis, 90);

   if (Hoption.Zscale) PaintPalette();

   delete axis;
   delete fLego; fLego = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Control function to draw a table using Delaunay triangles.

void THistPainter::PaintTriangles(Option_t *option)
{

   TGraphDelaunay2D *dt = nullptr;
   TGraphDelaunay *dtOld = nullptr;

   // Check if fH contains a TGraphDelaunay2D
   TList *hl = fH->GetListOfFunctions();
   dt = (TGraphDelaunay2D*)hl->FindObject("TGraphDelaunay2D");
   if (!dt) dtOld = (TGraphDelaunay*)hl->FindObject("TGraphDelaunay");
   if (!dt && !dtOld) return;

   // If needed, create a TGraph2DPainter
   if (!fGraph2DPainter) {
      if (dt) fGraph2DPainter = new TGraph2DPainter(dt);
      else fGraph2DPainter = new TGraph2DPainter(dtOld);
   }

   // Define the 3D view
   if (Hparam.zmin == 0 && Hparam.zmax == 0) {Hparam.zmin = -1; Hparam.zmax = 1;}
   if (Hoption.Same) {
      TView *viewsame = gPad->GetView();
      if (!viewsame) {
         Error("PaintTriangles", "no TView in current pad, do not use option SAME");
         return;
      }
      Double_t *rmin = viewsame->GetRmin();
      Double_t *rmax = viewsame->GetRmax();
      if (!rmin || !rmax) return;
      fXbuf[0] = rmin[0];
      fYbuf[0] = rmax[0];
      fXbuf[1] = rmin[1];
      fYbuf[1] = rmax[1];
      fXbuf[2] = rmin[2];
      fYbuf[2] = rmax[2];
   } else {
      fXbuf[0] = Hparam.xmin;
      fYbuf[0] = Hparam.xmax;
      fXbuf[1] = Hparam.ymin;
      fYbuf[1] = Hparam.ymax;
      fXbuf[2] = Hparam.zmin;
      fYbuf[2] = Hparam.zmax;
   }

   fLego = new TPainter3dAlgorithms(fXbuf, fYbuf);
   TView *view = gPad->GetView();
   if (!view) {
      Error("PaintTriangles", "no TView in current pad");
      return;
   }
   Double_t thedeg =  90 - gPad->GetTheta();
   Double_t phideg = -90 - gPad->GetPhi();
   Double_t psideg = view->GetPsi();
   Int_t irep;
   view->SetView(phideg, thedeg, psideg, irep);

   // Set color/style for back box
   fLego->SetFillStyle(gPad->GetFrameFillStyle());
   fLego->SetFillColor(gPad->GetFrameFillColor());
   fLego->TAttFill::Modify();
   Int_t backcolor = gPad->GetFrameFillColor();
   if (Hoption.System != kCARTESIAN) backcolor = 0;
   view->PadRange(backcolor);
   fLego->SetFillStyle(fH->GetFillStyle());
   fLego->SetFillColor(fH->GetFillColor());
   fLego->TAttFill::Modify();

   // Paint the Back Box if needed
   if (Hoption.BackBox && !Hoption.Same) {
      fLego->InitMoveScreen(-1.1,1.1);
      fLego->DefineGridLevels(fZaxis->GetNdivisions()%100);
      fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMove1);
      fLego->BackBox(90);
   }

   // Paint the triangles
   fGraph2DPainter->Paint(option);

   // Paint the Front Box if needed
   if (Hoption.FrontBox) {
      fLego->InitMoveScreen(-1.1,1.1);
      fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMove2);
      fLego->FrontBox(90);
   }

   // Paint the Axis if needed
   if (!Hoption.Axis && !Hoption.Same) {
      TGaxis *axis = new TGaxis();
      PaintLegoAxis(axis, 90);
      delete axis;
   }

   if (Hoption.Zscale) PaintPalette();

   delete fLego; fLego = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Define the color levels used to paint legos, surfaces etc..

void THistPainter::DefineColorLevels(Int_t ndivz)
{

   Int_t i, irep;

   // Initialize the color levels
   if (ndivz >= 100) {
      Warning("PaintSurface", "too many color levels, %d, reset to 8", ndivz);
      ndivz = 8;
   }
   Double_t *funlevel = new Double_t[ndivz+1];
   Int_t *colorlevel = new Int_t[ndivz+1];
   Int_t theColor;
   Int_t ncolors = gStyle->GetNumberOfColors();
   for (i = 0; i < ndivz; ++i) {
      funlevel[i] = fH->GetContourLevelPad(i);
      theColor = Int_t((i+0.99)*Float_t(ncolors)/Float_t(ndivz));
      colorlevel[i] = gStyle->GetColorPalette(theColor);
   }
   colorlevel[ndivz] = gStyle->GetColorPalette(ncolors-1);
   fLego->ColorFunction(ndivz, funlevel, colorlevel, irep);
   delete [] colorlevel;
   delete [] funlevel;
}

////////////////////////////////////////////////////////////////////////////////
/// [Control function to draw 2D/3D histograms (tables).](#HP01c)

void THistPainter::PaintTable(Option_t *option)
{

   // Fill Hparam structure with histo parameters
   if (!TableInit()) return;

   // Draw histogram frame
   PaintFrame();

   // If palette option not specified, delete a possible existing palette
   if (!Hoption.Zscale) {
      TObject *palette = fFunctions->FindObject("palette");
      if (palette) { fFunctions->Remove(palette); delete palette;}
   }

   // Do not draw the histogram. Only the attached functions will be drawn.
   if (Hoption.Func == 2) {
      if (Hoption.Zscale) {
         Int_t ndiv   = fH->GetContour();
         if (ndiv == 0 ) {
            ndiv = gStyle->GetNumberContours();
            fH->SetContour(ndiv);
         }
         PaintPalette();
      }

   // Draw the histogram according to the option
   } else {
      if (fH->InheritsFrom(TH2Poly::Class())) {
         if (Hoption.Fill)         PaintTH2PolyBins("f");
         if (Hoption.Color)        PaintTH2PolyColorLevels(option);
         if (Hoption.Scat)         PaintTH2PolyScatterPlot(option);
         if (Hoption.Text)         PaintTH2PolyText(option);
         if (Hoption.Line)         PaintTH2PolyBins("l");
         if (Hoption.Mark)         PaintTH2PolyBins("P");
      } else if (fH->GetEntries() != 0 && Hoption.Axis<=0) {
         if (Hoption.Scat)         PaintScatterPlot(option);
         if (Hoption.Arrow)        PaintArrows(option);
         if (Hoption.Box)          PaintBoxes(option);
         if (Hoption.Color) {
            if (Hoption.Color == 3) PaintColorLevelsFast(option);
            else                    PaintColorLevels(option);
         }
         if (Hoption.Contour)      PaintContour(option);
         if (Hoption.Text)         PaintText(option);
         if (Hoption.Error >= 100) Paint2DErrors(option);
         if (Hoption.Candle)       PaintCandlePlot(option);
      }
      if (Hoption.Lego)                     PaintLego(option);
      if (Hoption.Surf && !Hoption.Contour) PaintSurface(option);
      if (Hoption.Tri)                      PaintTriangles(option);
   }

   // Draw histogram title
   PaintTitle();

   // Draw the axes
   if (!Hoption.Lego && !Hoption.Surf &&
       !Hoption.Tri  && !(Hoption.Error >= 100)) PaintAxis(kFALSE);

   TF1 *fit  = 0;
   TIter next(fFunctions);
   TObject *obj;
   while ((obj = next())) {
      if (obj->InheritsFrom(TF1::Class())) {
         fit = (TF1*)obj;
         break;
      }
   }
   if ((Hoption.Same%10) != 1) {
      if (!fH->TestBit(TH1::kNoStats)) {  // bit set via TH1::SetStats
         if (!gPad->PadInSelectionMode() && !gPad->PadInHighlightMode()) {
            //ALWAYS executed on non-iOS platform.
            //On iOS, depends on mode.
            PaintStat2(gStyle->GetOptStat(),fit);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Control function to draw a TH2Poly bins' contours.
///
/// - option = "F" draw the bins as filled areas.
/// - option = "L" draw the bins as line.
/// - option = "P" draw the bins as markers.

void THistPainter::PaintTH2PolyBins(Option_t *option)
{

   //Do not highlight the histogram, if its part was picked.
   if (gPad->PadInHighlightMode() && gPad->GetSelected() != fH) return;

   TString opt = option;
   opt.ToLower();
   Bool_t line = kFALSE;
   Bool_t fill = kFALSE;
   Bool_t mark = kFALSE;
   if (opt.Contains("l")) line = kTRUE;
   if (opt.Contains("f")) fill = kTRUE;
   if (opt.Contains("p")) mark = kTRUE;

   TH2PolyBin  *b;
   Double_t z;

   TIter next(((TH2Poly*)fH)->GetBins());
   TObject *obj, *poly;

   while ((obj=next())) {
      b = (TH2PolyBin*)obj;
      z = b->GetContent();
      if (z==0 && Hoption.Zero) continue; // Do not draw empty bins in case of option "COL0 L"
      poly  = b->GetPolygon();

      // Paint the TGraph bins.
      if (poly->IsA() == TGraph::Class()) {
         TGraph *g  = (TGraph*)poly;
         g->TAttLine::Modify();
         g->TAttMarker::Modify();
         g->TAttFill::Modify();
         if (line) {
            Int_t fs = g->GetFillStyle();
            Int_t fc = g->GetFillColor();
            g->SetFillStyle(0);
            g->SetFillColor(g->GetLineColor());
            g->Paint("F");
            g->SetFillStyle(fs);
            g->SetFillColor(fc);
         }
         if (fill) g->Paint("F");
         if (mark) g->Paint("P");
      }

      // Paint the TMultiGraph bins.
      if (poly->IsA() == TMultiGraph::Class()) {
         TMultiGraph *mg = (TMultiGraph*)poly;
         TList *gl = mg->GetListOfGraphs();
         if (!gl) return;
         TGraph *g;
         TIter nextg(gl);
         while ((g = (TGraph*) nextg())) {
            g->TAttLine::Modify();
            g->TAttMarker::Modify();
            g->TAttFill::Modify();
            if (line) {
               Int_t fs = g->GetFillStyle();
               Int_t fc = g->GetFillColor();
               g->SetFillStyle(0);
               g->SetFillColor(g->GetLineColor());
               g->Paint("F");
               g->SetFillStyle(fs);
               g->SetFillColor(fc);
            }
            if (fill) g->Paint("F");
            if (mark) g->Paint("P");
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// [Control function to draw a TH2Poly as a color plot.](#HP20a)

void THistPainter::PaintTH2PolyColorLevels(Option_t *)
{

   //Do not highlight the histogram, if its part was picked.
   if (gPad->PadInHighlightMode() && gPad->GetSelected() != fH)
      return;

   Int_t ncolors, color, theColor;
   Double_t z, zc;
   Double_t zmin = fH->GetMinimum();
   Double_t zmax = fH->GetMaximum();
   if (Hoption.Logz) {
      if (zmax > 0) {
         if (zmin <= 0) zmin = TMath::Min((Double_t)1, (Double_t)0.001*zmax);
         zmin = TMath::Log10(zmin);
         zmax = TMath::Log10(zmax);
      } else {
         return;
      }
   }
   Double_t dz = zmax - zmin;

   // Initialize the levels on the Z axis
   ncolors     = gStyle->GetNumberOfColors();
   Int_t ndiv  = fH->GetContour();
   if (ndiv == 0 ) {
      ndiv = gStyle->GetNumberContours();
      fH->SetContour(ndiv);
   }
   Int_t ndivz  = TMath::Abs(ndiv);
   if (fH->TestBit(TH1::kUserContour) == 0) fH->SetContour(ndiv);
   Double_t scale = ndivz/dz;

   TH2PolyBin  *b;

   TIter next(((TH2Poly*)fH)->GetBins());
   TObject *obj, *poly;

   while ((obj=next())) {
      b     = (TH2PolyBin*)obj;
      poly  = b->GetPolygon();

      z = b->GetContent();
      if (z==0 && Hoption.Zero) continue;
      if (Hoption.Logz) {
         if (z > 0) z = TMath::Log10(z);
         else       z = zmin;
      }
      if (z < zmin) continue;

      // Define the bin color.
      if (fH->TestBit(TH1::kUserContour)) {
         zc = fH->GetContourLevelPad(0);
         if (z < zc) continue;
         color = -1;
         for (Int_t k=0; k<ndiv; k++) {
            zc = fH->GetContourLevelPad(k);
            if (z < zc) {
               continue;
            } else {
               color++;
            }
         }
      } else {
         color = Int_t(0.01+(z-zmin)*scale);
      }
      theColor = Int_t((color+0.99)*Float_t(ncolors)/Float_t(ndivz));
      if (theColor > ncolors-1) theColor = ncolors-1;

      // Paint the TGraph bins.
      if (poly->IsA() == TGraph::Class()) {
         TGraph *g  = (TGraph*)poly;
         g->SetFillColor(gStyle->GetColorPalette(theColor));
         g->TAttFill::Modify();
         g->Paint("F");
      }

      // Paint the TMultiGraph bins.
      if (poly->IsA() == TMultiGraph::Class()) {
         TMultiGraph *mg = (TMultiGraph*)poly;
         TList *gl = mg->GetListOfGraphs();
         if (!gl) return;
         TGraph *g;
         TIter nextg(gl);
         while ((g = (TGraph*) nextg())) {
            g->SetFillColor(gStyle->GetColorPalette(theColor));
            g->TAttFill::Modify();
            g->Paint("F");
         }
      }
   }
   if (Hoption.Zscale) PaintPalette();
}

////////////////////////////////////////////////////////////////////////////////
/// [Control function to draw a TH2Poly as a scatter plot.](#HP20a)

void THistPainter::PaintTH2PolyScatterPlot(Option_t *)
{

   //Do not highlight the histogram, if its part was selected.
   if (gPad->PadInHighlightMode() && gPad->GetSelected() != fH)
      return;

   Int_t k, loop, marker=0;
   Double_t z, xk,xstep, yk, ystep, xp, yp;
   Double_t scale = 1;
   Double_t zmin = fH->GetMinimum();
   Double_t zmax = fH->GetMaximum();
   if (Hoption.Logz) {
      if (zmax > 0) {
         if (zmin <= 0) zmin = TMath::Min((Double_t)1, (Double_t)0.001*zmax);
         zmin = TMath::Log10(zmin);
         zmax = TMath::Log10(zmax);
      } else {
         return;
      }
   }
   Double_t dz = zmax - zmin;
   scale = (kNMAX-1)/dz;


   // use an independent instance of a random generator
   // instead of gRandom to avoid conflicts and
   // to get same random numbers when drawing the same histogram
   TRandom2 random;

   TH2PolyBin  *b;

   TIter next(((TH2Poly*)fH)->GetBins());
   TObject *obj, *poly;

   Double_t maxarea = 0, a;
   while ((obj=next())) {
      b     = (TH2PolyBin*)obj;
      a     = b->GetArea();
      if (a>maxarea) maxarea = a;
   }

   next.Reset();

   while ((obj=next())) {
      b     = (TH2PolyBin*)obj;
      poly  = b->GetPolygon();
      z     = b->GetContent();
      if (z < zmin) z = zmin;
      if (z > zmax) z = zmax;
      if (Hoption.Logz) {
         if (z > 0) z = TMath::Log10(z) - zmin;
      } else {
         z    -=  zmin;
      }
      k     = Int_t((z*scale)*(b->GetArea()/maxarea));
      xk    = b->GetXMin();
      yk    = b->GetYMin();
      xstep = b->GetXMax()-xk;
      ystep = b->GetYMax()-yk;

      // Paint the TGraph bins.
      if (poly->IsA() == TGraph::Class()) {
         TGraph *g  = (TGraph*)poly;
         if (k <= 0 || z <= 0) continue;
         loop = 0;
         while (loop<k) {
            if (k+marker >= kNMAX) {
               gPad->PaintPolyMarker(marker, fXbuf, fYbuf);
               marker=0;
            }
            xp = (random.Rndm()*xstep) + xk;
            yp = (random.Rndm()*ystep) + yk;
            if (g->IsInside(xp,yp)) {
               fXbuf[marker] = xp;
               fYbuf[marker] = yp;
               marker++;
               loop++;
            }
         }
         if (marker > 0) gPad->PaintPolyMarker(marker, fXbuf, fYbuf);
      }

      // Paint the TMultiGraph bins.
      if (poly->IsA() == TMultiGraph::Class()) {
         TMultiGraph *mg = (TMultiGraph*)poly;
         TList *gl = mg->GetListOfGraphs();
         if (!gl) return;
         if (k <= 0 || z <= 0) continue;
         loop = 0;
         while (loop<k) {
            if (k+marker >= kNMAX) {
               gPad->PaintPolyMarker(marker, fXbuf, fYbuf);
               marker=0;
            }
            xp = (random.Rndm()*xstep) + xk;
            yp = (random.Rndm()*ystep) + yk;
            if (mg->IsInside(xp,yp)) {
               fXbuf[marker] = xp;
               fYbuf[marker] = yp;
               marker++;
               loop++;
            }
         }
         if (marker > 0) gPad->PaintPolyMarker(marker, fXbuf, fYbuf);
      }
   }
   PaintTH2PolyBins("l");
}

////////////////////////////////////////////////////////////////////////////////
/// [Control function to draw a TH2Poly as a text plot.](#HP20a)

void THistPainter::PaintTH2PolyText(Option_t *)
{

   TLatex text;
   text.SetTextFont(gStyle->GetTextFont());
   text.SetTextColor(fH->GetMarkerColor());
   text.SetTextSize(0.02*fH->GetMarkerSize());

   Double_t x, y, z, e, angle = 0;
   TString tt, tf;
   tf.Form("%s%s","%",gStyle->GetPaintTextFormat());
   if (Hoption.Text >= 1000) angle = Hoption.Text%1000;
   Int_t opt = (Int_t)Hoption.Text/1000;

   text.SetTextAlign(22);
   if (Hoption.Text ==  1) angle = 0;
   text.SetTextAngle(angle);
   text.TAttText::Modify();

   TH2PolyBin *b;

   TIter next(((TH2Poly*)fH)->GetBins());
   TObject *obj, *p;

   while ((obj=next())) {
      b = (TH2PolyBin*)obj;
      p = b->GetPolygon();
      x = (b->GetXMin()+b->GetXMax())/2;
      if (Hoption.Logx) {
         if (x > 0)  x  = TMath::Log10(x);
         else continue;
      }
      y = (b->GetYMin()+b->GetYMax())/2;
      if (Hoption.Logy) {
         if (y > 0)  y  = TMath::Log10(y);
         else continue;
      }
      z = b->GetContent();
      if (z < fH->GetMinimum() || (z == 0 && !Hoption.MinimumZero)) continue;
      if (opt==2) {
         e = fH->GetBinError(b->GetBinNumber());
         tf.Form("#splitline{%s%s}{#pm %s%s}",
                                    "%",gStyle->GetPaintTextFormat(),
                                    "%",gStyle->GetPaintTextFormat());
         tt.Form(tf.Data(),z,e);
      } else {
         tt.Form(tf.Data(),z);
      }
      if (opt==3) text.PaintLatex(x,y,angle,0.02*fH->GetMarkerSize(),p->GetName());
      else        text.PaintLatex(x,y,angle,0.02*fH->GetMarkerSize(),tt.Data());
   }

   PaintTH2PolyBins("l");
}

////////////////////////////////////////////////////////////////////////////////
/// [Control function to draw a 1D/2D histograms with the bin values.](#HP15)

void THistPainter::PaintText(Option_t *)
{

   TLatex text;
   text.SetTextFont(gStyle->GetTextFont());
   text.SetTextColor(fH->GetMarkerColor());
   text.SetTextSize(0.02*fH->GetMarkerSize());

   Double_t x, y, z, e, angle = 0;
   TString tt, tf;
   tf.Form("%s%s","%",gStyle->GetPaintTextFormat());
   if (Hoption.Text >= 1000) angle = Hoption.Text%1000;

   // 1D histograms
   if (fH->GetDimension() == 1) {
      Bool_t getentries = kFALSE;
      Double_t yt;
      TProfile *hp = (TProfile*)fH;
      if (Hoption.Text>2000 && fH->InheritsFrom(TProfile::Class())) {
         Hoption.Text = Hoption.Text-2000;
         getentries = kTRUE;
      }
      if (Hoption.Text ==  1) angle = 90;
      text.SetTextAlign(11);
      if (angle == 90) text.SetTextAlign(12);
      if (angle ==  0) text.SetTextAlign(21);
      text.TAttText::Modify();
      Double_t dt = 0.02*(gPad->GetY2()-gPad->GetY1());
      for (Int_t i=Hparam.xfirst; i<=Hparam.xlast;i++) {
         if (Hoption.Bar) {
            x  = fH->GetXaxis()->GetBinLowEdge(i)+
                 fH->GetXaxis()->GetBinWidth(i)*
                 (fH->GetBarOffset()+0.5*fH->GetBarWidth());
         } else {
            x  = fH->GetXaxis()->GetBinCenter(i);
         }
         y  = fH->GetBinContent(i);
         yt = y;
         if (Hoption.MinimumZero && y<0) y = 0;
         if (getentries) yt = hp->GetBinEntries(i);
         if (yt == 0.) continue;
         tt.Form(tf.Data(),yt);
         if (Hoption.Logx) {
            if (x > 0)  x  = TMath::Log10(x);
            else continue;
         }
         if (Hoption.Logy) {
            if (y > 0)  y  = TMath::Log10(y);
            else continue;
         }
         if (y >= gPad->GetY2()) continue;
         if (y <= gPad->GetY1()) continue;
         text.PaintLatex(x,y+0.2*dt,angle,0.02*fH->GetMarkerSize(),tt.Data());
      }

   // 2D histograms
   } else {
      text.SetTextAlign(22);
      if (Hoption.Text ==  1) angle = 0;
      text.SetTextAngle(angle);
      text.TAttText::Modify();
      for (Int_t j=Hparam.yfirst; j<=Hparam.ylast;j++) {
         y    = fYaxis->GetBinCenter(j);
         if (Hoption.Logy) {
            if (y > 0)  y  = TMath::Log10(y);
            else continue;
         }
         for (Int_t i=Hparam.xfirst; i<=Hparam.xlast;i++) {
            Int_t bin  = j*(fXaxis->GetNbins()+2) + i;
            x    = fXaxis->GetBinCenter(i);
            if (Hoption.Logx) {
               if (x > 0)  x  = TMath::Log10(x);
               else continue;
            }
            if (!IsInside(x,y)) continue;
            z = fH->GetBinContent(bin);
            if (z < Hparam.zmin || (z == 0 && !Hoption.MinimumZero)) continue;
            if (Hoption.Text>2000) {
               e = fH->GetBinError(bin);
               tf.Form("#splitline{%s%s}{#pm %s%s}",
                                          "%",gStyle->GetPaintTextFormat(),
                                          "%",gStyle->GetPaintTextFormat());
               tt.Form(tf.Data(),z,e);
            } else {
               tt.Form(tf.Data(),z);
            }
            text.PaintLatex(x,y+fH->GetBarOffset()*fYaxis->GetBinWidth(j),
                            angle,0.02*fH->GetMarkerSize(),tt.Data());
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// [Control function to draw a 3D implicit functions.](#HP27)

void THistPainter::PaintTF3()
{

   Int_t irep;

   TGaxis *axis = new TGaxis();
   TAxis *xaxis = fH->GetXaxis();
   TAxis *yaxis = fH->GetYaxis();
   TAxis *zaxis = fH->GetZaxis();

   fXbuf[0] = xaxis->GetBinLowEdge(xaxis->GetFirst());
   fYbuf[0] = xaxis->GetBinUpEdge(xaxis->GetLast());
   fXbuf[1] = yaxis->GetBinLowEdge(yaxis->GetFirst());
   fYbuf[1] = yaxis->GetBinUpEdge(yaxis->GetLast());
   fXbuf[2] = zaxis->GetBinLowEdge(zaxis->GetFirst());
   fYbuf[2] = zaxis->GetBinUpEdge(zaxis->GetLast());

   fLego = new TPainter3dAlgorithms(fXbuf, fYbuf);

   TView *view = gPad->GetView();
   if (!view) {
      Error("PaintTF3", "no TView in current pad");
      return;
   }
   Double_t thedeg =  90 - gPad->GetTheta();
   Double_t phideg = -90 - gPad->GetPhi();
   Double_t psideg = view->GetPsi();
   view->SetView(phideg, thedeg, psideg, irep);

   fLego->InitMoveScreen(-1.1,1.1);

   if (Hoption.BackBox) {
      fLego->DefineGridLevels(fZaxis->GetNdivisions()%100);
      fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMove1);
      fLego->BackBox(90);
   }

   fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMode1);

   fLego->ImplicitFunction(fCurrentF3, fXbuf, fYbuf, fH->GetNbinsX(),
                                       fH->GetNbinsY(),
                                       fH->GetNbinsZ(), "BF");

   if (Hoption.FrontBox) {
      fLego->InitMoveScreen(-1.1,1.1);
      fLego->SetDrawFace(&TPainter3dAlgorithms::DrawFaceMove2);
      fLego->FrontBox(90);
   }
   if (!Hoption.Axis && !Hoption.Same) PaintLegoAxis(axis, 90);

   PaintTitle();

   delete axis;
   delete fLego; fLego = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the histogram title
///
/// The title is drawn according to the title alignment returned by
/// `GetTitleAlign()`. It is a 2 digits integer): hv
///
/// where `h` is the horizontal alignment and `v` is the
/// vertical alignment.
///
/// - `h` can get the values 1 2 3 for left, center, and right
/// - `v` can get the values 1 2 3 for bottom, middle and top
///
/// for instance the default alignment is: 13 (left top)

void THistPainter::PaintTitle()
{
   // probably best place for calls PaintHighlightBin
   // calls after paint histo (1D or 2D) and before paint title and stats
   if (!gPad->GetView()) PaintHighlightBin();

   if (Hoption.Same) return;
   if (fH->TestBit(TH1::kNoTitle)) return;
   Int_t nt = strlen(fH->GetTitle());
   TPaveText *title = 0;
   TObject *obj;
   TIter next(gPad->GetListOfPrimitives());
   while ((obj = next())) {
      if (!obj->InheritsFrom(TPaveText::Class())) continue;
      title = (TPaveText*)obj;
      if (strcmp(title->GetName(),"title")) {title = 0; continue;}
      break;
   }
   if (nt == 0 || gStyle->GetOptTitle() <= 0) {
      if (title) delete title;
      return;
   }
   Double_t ht = gStyle->GetTitleH();
   Double_t wt = gStyle->GetTitleW();

   if (ht <= 0) {
      if (gStyle->GetTitleFont("")%10 == 3) {
         Double_t hw = TMath::Max((Double_t)gPad->XtoPixel(gPad->GetX2()),
                                  (Double_t)gPad->YtoPixel(gPad->GetY1()));
         ht = 1.1*(gStyle->GetTitleSize("")/hw);
      } else {
         ht = 1.1*gStyle->GetTitleFontSize();
      }
   }
   if (ht <= 0) ht = 0.05;
   if (wt <= 0) {
      TLatex l;
      l.SetTextSize(ht);
      l.SetTitle(fH->GetTitle());
      // adjustment in case the title has several lines (#splitline)
      ht = TMath::Max(ht, 1.2*l.GetYsize()/(gPad->GetY2() - gPad->GetY1()));
      Double_t wndc = l.GetXsize()/(gPad->GetX2() - gPad->GetX1());
      wt = TMath::Min(0.7, 0.02+wndc);
   }
   if (title) {
      TText *t0 = (TText*)title->GetLine(0);
      if (t0) {
         if (!strcmp(t0->GetTitle(),fH->GetTitle())) return;
         t0->SetTitle(fH->GetTitle());
         if (wt > 0) title->SetX2NDC(title->GetX1NDC()+wt);
      }
      return;
   }

   Int_t talh = gStyle->GetTitleAlign()/10;
   if (talh < 1) talh = 1; else if (talh > 3) talh = 3;
   Int_t talv = gStyle->GetTitleAlign()%10;
   if (talv < 1) talv = 1; else if (talv > 3) talv = 3;
   Double_t xpos, ypos;
   xpos = gStyle->GetTitleX();
   ypos = gStyle->GetTitleY();
   if (talh == 2) xpos = xpos-wt/2.;
   if (talh == 3) xpos = xpos-wt;
   if (talv == 2) ypos = ypos+ht/2.;
   if (talv == 1) ypos = ypos+ht;

   TPaveText *ptitle = new TPaveText(xpos, ypos-ht, xpos+wt, ypos,"blNDC");

   //     box with the histogram title
   ptitle->SetFillColor(gStyle->GetTitleFillColor());
   ptitle->SetFillStyle(gStyle->GetTitleStyle());
   ptitle->SetName("title");
   ptitle->SetBorderSize(gStyle->GetTitleBorderSize());
   ptitle->SetTextColor(gStyle->GetTitleTextColor());
   ptitle->SetTextFont(gStyle->GetTitleFont(""));
   if (gStyle->GetTitleFont("")%10 > 2)
      ptitle->SetTextSize(gStyle->GetTitleFontSize());
   ptitle->AddText(fH->GetTitle());
   ptitle->SetBit(kCanDelete);
   ptitle->Draw();
   ptitle->Paint();

   if(!gPad->IsEditable()) delete ptitle;
}

////////////////////////////////////////////////////////////////////////////////
/// Process message `mess`.

void THistPainter::ProcessMessage(const char *mess, const TObject *obj)
{
   if (!strcmp(mess,"SetF3")) {
      fCurrentF3 = (TF3 *)obj;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Static function.
///
/// Convert Right Ascension, Declination to X,Y using an AITOFF projection.
/// This procedure can be used to create an all-sky map in Galactic
/// coordinates with an equal-area Aitoff projection.  Output map
/// coordinates are zero longitude centered.
/// Also called Hammer-Aitoff projection (first presented by Ernst von Hammer in 1892)
///
/// source: GMT
///
/// code from  Ernst-Jan Buis

Int_t THistPainter::ProjectAitoff2xy(Double_t l, Double_t b, Double_t &Al, Double_t &Ab)
{

   Double_t x, y;

   Double_t alpha2 = (l/2)*TMath::DegToRad();
   Double_t delta  = b*TMath::DegToRad();
   Double_t r2     = TMath::Sqrt(2.);
   Double_t f      = 2*r2/TMath::Pi();
   Double_t cdec   = TMath::Cos(delta);
   Double_t denom  = TMath::Sqrt(1. + cdec*TMath::Cos(alpha2));
   x      = cdec*TMath::Sin(alpha2)*2.*r2/denom;
   y      = TMath::Sin(delta)*r2/denom;
   x     *= TMath::RadToDeg()/f;
   y     *= TMath::RadToDeg()/f;
   //  x *= -1.; // for a skymap swap left<->right
   Al = x;
   Ab = y;

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function
///
/// Probably the most famous of the various map projections, the Mercator projection
/// takes its name from Mercator who presented it in 1569. It is a cylindrical, conformal projection
/// with no distortion along the equator.
/// The Mercator projection has been used extensively for world maps in which the distortion towards
/// the polar regions grows rather large, thus incorrectly giving the impression that, for example,
/// Greenland is larger than South America. In reality, the latter is about eight times the size of
/// Greenland. Also, the Former Soviet Union looks much bigger than Africa or South America. One may wonder
/// whether this illusion has had any influence on U.S. foreign policy.' (Source: GMT)
/// code from  Ernst-Jan Buis

Int_t THistPainter::ProjectMercator2xy(Double_t l, Double_t b, Double_t &Al, Double_t &Ab)
{

   Al = l;
   Double_t aid = TMath::Tan((TMath::PiOver2() + b*TMath::DegToRad())/2);
   Ab = TMath::Log(aid);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function code from  Ernst-Jan Buis

Int_t THistPainter::ProjectSinusoidal2xy(Double_t l, Double_t b, Double_t &Al, Double_t &Ab)
{

   Al = l*cos(b*TMath::DegToRad());
   Ab = b;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function code from  Ernst-Jan Buis

Int_t THistPainter::ProjectParabolic2xy(Double_t l, Double_t b, Double_t &Al, Double_t &Ab)
{

   Al = l*(2.*TMath::Cos(2*b*TMath::DegToRad()/3) - 1);
   Ab = 180*TMath::Sin(b*TMath::DegToRad()/3);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Recompute the histogram range following graphics operations.

void THistPainter::RecalculateRange()
{

   if (Hoption.Same) return;

   //     Compute x,y range
   Double_t xmin = Hparam.xmin;
   Double_t xmax = Hparam.xmax;
   Double_t ymin = Hparam.ymin;
   Double_t ymax = Hparam.ymax;

   Double_t xmin_aid, ymin_aid, xmax_aid, ymax_aid;
   if (Hoption.Proj ==1) {
      // TODO : check x range not lower than -180 and not higher than 180
      THistPainter::ProjectAitoff2xy(Hparam.xmin, Hparam.ymin, xmin_aid, ymin_aid);
      THistPainter::ProjectAitoff2xy(Hparam.xmin, Hparam.ymax, xmin,     ymax_aid);
      THistPainter::ProjectAitoff2xy(Hparam.xmax, Hparam.ymax, xmax_aid, ymax);
      THistPainter::ProjectAitoff2xy(Hparam.xmax, Hparam.ymin, xmax,     ymin);

      if (xmin > xmin_aid) xmin = xmin_aid;
      if (ymin > ymin_aid) ymin = ymin_aid;
      if (xmax < xmax_aid) xmax = xmax_aid;
      if (ymax < ymax_aid) ymax = ymax_aid;
      if (Hparam.ymin<0 && Hparam.ymax>0) {
         // there is an  'equator', check its range in the plot..
         THistPainter::ProjectAitoff2xy(Hparam.xmin*0.9999, 0, xmin_aid, ymin_aid);
         THistPainter::ProjectAitoff2xy(Hparam.xmax*0.9999, 0, xmax_aid, ymin_aid);
         if (xmin >xmin_aid) xmin = xmin_aid;
         if (xmax <xmax_aid) xmax = xmax_aid;
      }
      if (Hparam.xmin<0 && Hparam.xmax>0) {
         THistPainter::ProjectAitoff2xy(0, Hparam.ymin, xmin_aid, ymin_aid);
         THistPainter::ProjectAitoff2xy(0, Hparam.ymax, xmax_aid, ymax_aid);
         if (ymin >ymin_aid) ymin = ymin_aid;
         if (ymax <ymax_aid) ymax = ymax_aid;
      }
   } else if ( Hoption.Proj ==2) {
      if (Hparam.ymin <= -90 || Hparam.ymax >=90) {
         Warning("Mercator Projection", "Latitude out of range %f or %f", Hparam.ymin, Hparam.ymax);
         Hoption.Proj = 0;
      } else {
         THistPainter::ProjectMercator2xy(Hparam.xmin, Hparam.ymin, xmin, ymin);
         THistPainter::ProjectMercator2xy(Hparam.xmax, Hparam.ymax, xmax, ymax);
      }
   } else if (Hoption.Proj == 3) {
      THistPainter::ProjectSinusoidal2xy(Hparam.xmin, Hparam.ymin, xmin_aid, ymin_aid);
      THistPainter::ProjectSinusoidal2xy(Hparam.xmin, Hparam.ymax, xmin,     ymax_aid);
      THistPainter::ProjectSinusoidal2xy(Hparam.xmax, Hparam.ymax, xmax_aid, ymax);
      THistPainter::ProjectSinusoidal2xy(Hparam.xmax, Hparam.ymin, xmax,     ymin);

      if (xmin > xmin_aid) xmin = xmin_aid;
      if (ymin > ymin_aid) ymin = ymin_aid;
      if (xmax < xmax_aid) xmax = xmax_aid;
      if (ymax < ymax_aid) ymax = ymax_aid;
      if (Hparam.ymin<0 && Hparam.ymax>0) {
         THistPainter::ProjectSinusoidal2xy(Hparam.xmin, 0, xmin_aid, ymin_aid);
         THistPainter::ProjectSinusoidal2xy(Hparam.xmax, 0, xmax_aid, ymin_aid);
         if (xmin >xmin_aid) xmin = xmin_aid;
         if (xmax <xmax_aid) xmax = xmax_aid;
      }
      if (Hparam.xmin<0 && Hparam.xmax>0) {
         THistPainter::ProjectSinusoidal2xy(0,Hparam.ymin, xmin_aid, ymin_aid);
         THistPainter::ProjectSinusoidal2xy(0, Hparam.ymax, xmax_aid, ymin_aid);
         if (ymin >ymin_aid) ymin = ymin_aid;
         if (ymax <ymax_aid) ymax = ymax_aid;
      }
   } else if (Hoption.Proj == 4) {
      THistPainter::ProjectParabolic2xy(Hparam.xmin, Hparam.ymin, xmin_aid, ymin_aid);
      THistPainter::ProjectParabolic2xy(Hparam.xmin, Hparam.ymax, xmin,     ymax_aid);
      THistPainter::ProjectParabolic2xy(Hparam.xmax, Hparam.ymax, xmax_aid, ymax);
      THistPainter::ProjectParabolic2xy(Hparam.xmax, Hparam.ymin, xmax,     ymin);

      if (xmin > xmin_aid) xmin = xmin_aid;
      if (ymin > ymin_aid) ymin = ymin_aid;
      if (xmax < xmax_aid) xmax = xmax_aid;
      if (ymax < ymax_aid) ymax = ymax_aid;
      if (Hparam.ymin<0 && Hparam.ymax>0) {
         THistPainter::ProjectParabolic2xy(Hparam.xmin, 0, xmin_aid, ymin_aid);
         THistPainter::ProjectParabolic2xy(Hparam.xmax, 0, xmax_aid, ymin_aid);
         if (xmin >xmin_aid) xmin = xmin_aid;
         if (xmax <xmax_aid) xmax = xmax_aid;
      }
      if (Hparam.xmin<0 && Hparam.xmax>0) {
         THistPainter::ProjectParabolic2xy(0, Hparam.ymin, xmin_aid, ymin_aid);
         THistPainter::ProjectParabolic2xy(0, Hparam.ymax, xmax_aid, ymin_aid);
         if (ymin >ymin_aid) ymin = ymin_aid;
         if (ymax <ymax_aid) ymax = ymax_aid;
      }
   }
   Hparam.xmin= xmin;
   Hparam.xmax= xmax;
   Hparam.ymin= ymin;
   Hparam.ymax= ymax;

   Double_t dx   = xmax-xmin;
   Double_t dy   = ymax-ymin;
   Double_t dxr  = dx/(1 - gPad->GetLeftMargin()   - gPad->GetRightMargin());
   Double_t dyr  = dy/(1 - gPad->GetBottomMargin() - gPad->GetTopMargin());

   // Range() could change the size of the pad pixmap and therefore should
   // be called before the other paint routines
   gPad->Range(xmin - dxr*gPad->GetLeftMargin(),
                      ymin - dyr*gPad->GetBottomMargin(),
                      xmax + dxr*gPad->GetRightMargin(),
                      ymax + dyr*gPad->GetTopMargin());
   gPad->RangeAxis(xmin, ymin, xmax, ymax);
}

////////////////////////////////////////////////////////////////////////////////
/// Set current histogram to `h`

void THistPainter::SetHistogram(TH1 *h)
{

   if (h == 0)  return;
   fH = h;
   fXaxis = h->GetXaxis();
   fYaxis = h->GetYaxis();
   fZaxis = h->GetZaxis();
   fFunctions = fH->GetListOfFunctions();
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize various options to draw 2D histograms.

Int_t THistPainter::TableInit()
{

   static const char *where = "TableInit";

   Int_t first, last;
   Double_t yMARGIN= gStyle->GetHistTopMargin();
   Double_t zmin, zmax;
   Int_t maximum = 0;
   Int_t minimum = 0;
   if (fH->GetMaximumStored() != -1111) maximum = 1;
   if (fH->GetMinimumStored() != -1111) minimum = 1;

   //    -----------------  Compute X axis parameters
   first           = fXaxis->GetFirst();
   last            = fXaxis->GetLast();
   Hparam.xlast    = last;
   Hparam.xfirst   = first;
   Hparam.xlowedge = fXaxis->GetBinLowEdge(first);
   Hparam.xbinsize = fXaxis->GetBinWidth(first);
   Hparam.xmin     = Hparam.xlowedge;
   Hparam.xmax     = fXaxis->GetBinLowEdge(last)+fXaxis->GetBinWidth(last);

   //       if log scale in X, replace xmin,max by the log
   if (Hoption.Logx) {
   //   find the first edge of a bin that is > 0
      if (Hparam.xlowedge <=0 ) {
         Hparam.xlowedge = fXaxis->GetBinUpEdge(fXaxis->FindFixBin(0.01*Hparam.xbinsize));
         Hparam.xmin  = Hparam.xlowedge;
      }
      if (Hparam.xmin <=0 || Hparam.xmax <=0) {
         Error(where, "cannot set X axis to log scale");
         return 0;
      }
      Hparam.xfirst= fXaxis->FindFixBin(Hparam.xmin);
      if (Hparam.xfirst < first) Hparam.xfirst = first;
      Hparam.xlast = fXaxis->FindFixBin(Hparam.xmax);
      if (Hparam.xlast > last) Hparam.xlast = last;
      Hparam.xmin  = TMath::Log10(Hparam.xmin);
      Hparam.xmax  = TMath::Log10(Hparam.xmax);
   }

   //    -----------------  Compute Y axis parameters
   first           = fYaxis->GetFirst();
   last            = fYaxis->GetLast();
   Hparam.ylast    = last;
   Hparam.yfirst   = first;
   Hparam.ylowedge = fYaxis->GetBinLowEdge(first);
   Hparam.ybinsize = fYaxis->GetBinWidth(first);
   if (!Hparam.ybinsize) Hparam.ybinsize = 1;
   Hparam.ymin     = Hparam.ylowedge;
   Hparam.ymax     = fYaxis->GetBinLowEdge(last)+fYaxis->GetBinWidth(last);

   //       if log scale in Y, replace ymin,max by the log
   if (Hoption.Logy) {
      if (Hparam.ylowedge <=0 ) {
         Hparam.ylowedge = fYaxis->GetBinUpEdge(fYaxis->FindFixBin(0.01*Hparam.ybinsize));
         Hparam.ymin  = Hparam.ylowedge;
      }
      if (Hparam.ymin <=0 || Hparam.ymax <=0) {
         Error(where, "cannot set Y axis to log scale");
         return 0;
      }
      Hparam.yfirst= fYaxis->FindFixBin(Hparam.ymin);
      if (Hparam.yfirst < first) Hparam.yfirst = first;
      Hparam.ylast = fYaxis->FindFixBin(Hparam.ymax);
      if (Hparam.ylast > last) Hparam.ylast = last;
      Hparam.ymin  = TMath::Log10(Hparam.ymin);
      Hparam.ymax  = TMath::Log10(Hparam.ymax);
   }


   //    -----------------  Compute Z axis parameters
   Double_t bigp = TMath::Power(10,32);
   zmax = -bigp;
   zmin = bigp;
   Double_t c1, e1;
   Double_t allchan = 0;
   for (Int_t j=Hparam.yfirst; j<=Hparam.ylast;j++) {
      for (Int_t i=Hparam.xfirst; i<=Hparam.xlast;i++) {
         c1 = fH->GetBinContent(i,j);
         zmax = TMath::Max(zmax,c1);
         if (Hoption.Error) {
            e1 = fH->GetBinError(i,j);
            zmax = TMath::Max(zmax,c1+e1);
         }
         zmin = TMath::Min(zmin,c1);
         allchan += c1;
      }
   }

   //     Take into account maximum , minimum

   if (maximum) zmax = fH->GetMaximumStored();
   if (minimum) zmin = fH->GetMinimumStored();
   if (Hoption.Logz && zmax < 0) {
      if (!Hoption.Same) Error(where, "log scale is requested but maximum is less or equal 0 (%f)", zmax);
      return 0;
   } else if (Hoption.Logz && zmin>=0 && zmax==0) { // empty histogram in log scale
      zmin = 0.01;
      zmax = 10.;
   }
   if (zmin >= zmax) {
      if (Hoption.Logz) {
         if (zmax > 0) zmin = 0.001*zmax;
         else {
            if (!Hoption.Same) Error(where, "log scale is requested but maximum is less or equal 0 (%f)", zmax);
            return 0;
         }
      }
   }

   //     take into account normalization factor
   Hparam.allchan = allchan;
   Double_t factor = allchan;
   if (fH->GetNormFactor() > 0) factor = fH->GetNormFactor();
   if (allchan) factor /= allchan;
   if (factor == 0) factor = 1;
   Hparam.factor = factor;
   zmax = factor*zmax;
   zmin = factor*zmin;
   c1 = zmax;
   if (TMath::Abs(zmin) > TMath::Abs(c1)) c1 = zmin;

   //         For log scales, histogram coordinates are log10(ymin) and
   //         log10(ymax). Final adjustment (if not option "Same")
   //         or "+" for ymax) of ymax and ymin for logarithmic scale, if
   //         Maximum and Minimum are not defined.
   if (Hoption.Logz) {
      if (zmin <= 0) {
         zmin = TMath::Min((Double_t)1, (Double_t)0.001*zmax);
         fH->SetMinimum(zmin);
      }
      zmin = TMath::Log10(zmin);
      if (!minimum) zmin += TMath::Log10(0.5);
      zmax = TMath::Log10(zmax);
      if (!maximum) zmax += TMath::Log10(2*(0.9/0.95));
      goto LZMIN;
   }

   //         final adjustment of YMAXI for linear scale (if not option "Same"):
   //         decrease histogram height to MAX% of allowed height if HMAXIM
   //         has not been called.
   //         MAX% is the value in percent which has been set in HPLSET
   //         (default is 90%).
   if (!maximum) {
      zmax += yMARGIN*(zmax-zmin);
   }

   //         final adjustment of ymin for linear scale.
   //         if minimum is not set , then ymin is set to zero if >0
   //         or to ymin - yMARGIN if <0.
   if (!minimum) {
      if (Hoption.MinimumZero) {
         if (zmin >= 0) zmin = 0;
         else           zmin -= yMARGIN*(zmax-zmin);
      } else {
         Double_t dzmin = yMARGIN*(zmax-zmin);
         if (zmin >= 0 && (zmin-dzmin <= 0)) zmin  = 0;
         else                                zmin -= dzmin;
      }
   }

LZMIN:
   Hparam.zmin = zmin;
   Hparam.zmax = zmax;

   //     Set bar offset and width
   Hparam.baroffset = fH->GetBarOffset();
   Hparam.barwidth  = fH->GetBarWidth();

   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// This function returns the best format to print the error value (e)
/// knowing the parameter value (v) and the format (f) used to print it.

const char * THistPainter::GetBestFormat(Double_t v, Double_t e, const char *f)
{

   static TString ef;
   TString tf, tv;

   // print v with the format f in tv.
   tf.Form("%s%s","%",f);
   tv.Form(tf.Data(),v);

   // Analyse tv.
   int ie = tv.Index("e");
   int iE = tv.Index("E");
   int id = tv.Index(".");

   // v has been printed with the exponent notation.
   // There is 2 cases, the exponent is positive or negative
   if (ie >= 0 || iE >= 0) {
      if (tv.Index("+") >= 0) {
         if (e < 1) {
            ef.Form("%s.1f","%");
         } else {
            if (ie >= 0) {
               ef.Form("%s.%de","%",ie-id-1);
            } else {
               ef.Form("%s.%dE","%",iE-id-1);
            }
         }
      } else {
         if (ie >= 0) {
            ef.Form("%s.%de","%",ie-id-1);
         } else {
            ef.Form("%s.%dE","%",iE-id-1);
         }
      }

   // There is not '.' in tv. e will be printed with one decimal digit.
   } else if (id < 0) {
      ef.Form("%s.1f","%");

   // There is a '.' in tv and no exponent notation. e's decimal part will
   // have the same number of digits as v's one.
   } else {
      ef.Form("%s.%df","%",tv.Length()-id-1);
   }

   return ef.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Set projection.

void THistPainter::SetShowProjection(const char *option,Int_t nbins)
{

   if (fShowProjection) return;
   TString opt = option;
   opt.ToLower();
   Int_t projection = 0;
   if (opt.Contains("x"))  projection = 1;
   if (opt.Contains("y"))  projection = 2;
   if (opt.Contains("z"))  projection = 3;
   if (opt.Contains("xy")) projection = 4;
   if (opt.Contains("yx")) projection = 5;
   if (opt.Contains("xz")) projection = 6;
   if (opt.Contains("zx")) projection = 7;
   if (opt.Contains("yz")) projection = 8;
   if (opt.Contains("zy")) projection = 9;
   if (projection < 4) fShowOption = option+1;
   else                fShowOption = option+2;
   fShowProjection = projection+100*nbins;
   gROOT->MakeDefCanvas();
   gPad->SetName(Form("c_%lx_projection_%d", (ULong_t)fH, fShowProjection));
   gPad->SetGrid();
}

////////////////////////////////////////////////////////////////////////////////
/// Show projection onto X.

void THistPainter::ShowProjectionX(Int_t /*px*/, Int_t py)
{

   Int_t nbins = (Int_t)fShowProjection/100;
   gPad->SetDoubleBuffer(0); // turn off double buffer mode
   gVirtualX->SetDrawMode(TVirtualX::kInvert); // set the drawing mode to XOR mode

   // Erase old position and draw a line at current position
   static int pyold1 = 0;
   static int pyold2 = 0;
   float uxmin = gPad->GetUxmin();
   float uxmax = gPad->GetUxmax();
   int pxmin   = gPad->XtoAbsPixel(uxmin);
   int pxmax   = gPad->XtoAbsPixel(uxmax);
   Float_t upy = gPad->AbsPixeltoY(py);
   Float_t y   = gPad->PadtoY(upy);
   Int_t biny1 = fH->GetYaxis()->FindBin(y);
   Int_t biny2 = TMath::Min(biny1+nbins-1, fH->GetYaxis()->GetNbins());
   Int_t py1   = gPad->YtoAbsPixel(fH->GetYaxis()->GetBinLowEdge(biny1));
   Int_t py2   = gPad->YtoAbsPixel(fH->GetYaxis()->GetBinUpEdge(biny2));

   if (pyold1 || pyold2) gVirtualX->DrawBox(pxmin,pyold1,pxmax,pyold2,TVirtualX::kFilled);
   gVirtualX->DrawBox(pxmin,py1,pxmax,py2,TVirtualX::kFilled);
   pyold1 = py1;
   pyold2 = py2;

   // Create or set the new canvas proj x
   TVirtualPad *padsav = gPad;
   TVirtualPad *c = (TVirtualPad*)gROOT->GetListOfCanvases()->FindObject(Form("c_%lx_projection_%d",
                                                                              (ULong_t)fH, fShowProjection));
   if (c) {
      c->Clear();
   } else {
      fShowProjection = 0;
      pyold1 = 0;
      pyold2 = 0;
      return;
   }
   c->cd();
   c->SetLogy(padsav->GetLogz());
   c->SetLogx(padsav->GetLogx());

   // Draw slice corresponding to mouse position
   TString prjName = TString::Format("slice_px_of_%s",fH->GetName());
   TH1D *hp = ((TH2*)fH)->ProjectionX(prjName, biny1, biny2);
   if (hp) {
      hp->SetFillColor(38);
      // apply a patch from Oliver Freyermuth to set the title in the projection
      // using the range of the projected Y values
      if (biny1 == biny2) {
         Double_t valueFrom   = fH->GetYaxis()->GetBinLowEdge(biny1);
         Double_t valueTo     = fH->GetYaxis()->GetBinUpEdge(biny1);
         // Limit precision to 1 digit more than the difference between upper and lower bound (to also catch 121.5-120.5).
         Int_t valuePrecision = -TMath::Nint(TMath::Log10(valueTo-valueFrom))+1;
         if (fH->GetYaxis()->GetLabels() != NULL) {
            hp->SetTitle(TString::Format("ProjectionX of biny=%d [y=%.*lf..%.*lf] %s", biny1, valuePrecision, valueFrom, valuePrecision, valueTo, fH->GetYaxis()->GetBinLabel(biny1)));
         } else {
            hp->SetTitle(TString::Format("ProjectionX of biny=%d [y=%.*lf..%.*lf]", biny1, valuePrecision, valueFrom, valuePrecision, valueTo));
         }
      } else {
         Double_t valueFrom   = fH->GetYaxis()->GetBinLowEdge(biny1);
         Double_t valueTo     = fH->GetYaxis()->GetBinUpEdge(biny2);
         // Limit precision to 1 digit more than the difference between upper and lower bound (to also catch 121.5-120.5).
         // biny1 is used here to get equal precision no matter how large the binrange is,
         // otherwise precision may change when moving the mouse to the histogram boundaries (limiting effective binrange).
         Int_t valuePrecision = -TMath::Nint(TMath::Log10(fH->GetYaxis()->GetBinUpEdge(biny1)-valueFrom))+1;
         if (fH->GetYaxis()->GetLabels() != NULL) {
            hp->SetTitle(TString::Format("ProjectionX of biny=[%d,%d] [y=%.*lf..%.*lf] [%s..%s]", biny1, biny2, valuePrecision, valueFrom, valuePrecision, valueTo, fH->GetYaxis()->GetBinLabel(biny1), fH->GetYaxis()->GetBinLabel(biny2)));
         } else {
            hp->SetTitle(TString::Format("ProjectionX of biny=[%d,%d] [y=%.*lf..%.*lf]", biny1, biny2, valuePrecision, valueFrom, valuePrecision, valueTo));
         }
      }
      hp->SetXTitle(fH->GetXaxis()->GetTitle());
      hp->SetYTitle("Number of Entries");
      hp->Draw();
      c->Update();
      padsav->cd();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Show projection onto Y.

void THistPainter::ShowProjectionY(Int_t px, Int_t /*py*/)
{

   Int_t nbins = (Int_t)fShowProjection/100;
   gPad->SetDoubleBuffer(0);             // turn off double buffer mode
   gVirtualX->SetDrawMode(TVirtualX::kInvert);  // set the drawing mode to XOR mode

   // Erase old position and draw a line at current position
   static int pxold1 = 0;
   static int pxold2 = 0;
   float uymin = gPad->GetUymin();
   float uymax = gPad->GetUymax();
   int pymin   = gPad->YtoAbsPixel(uymin);
   int pymax   = gPad->YtoAbsPixel(uymax);
   Float_t upx = gPad->AbsPixeltoX(px);
   Float_t x   = gPad->PadtoX(upx);
   Int_t binx1 = fH->GetXaxis()->FindBin(x);
   Int_t binx2 = TMath::Min(binx1+nbins-1, fH->GetXaxis()->GetNbins());
   Int_t px1   = gPad->XtoAbsPixel(fH->GetXaxis()->GetBinLowEdge(binx1));
   Int_t px2   = gPad->XtoAbsPixel(fH->GetXaxis()->GetBinUpEdge(binx2));

   if (pxold1 || pxold2) gVirtualX->DrawBox(pxold1,pymin,pxold2,pymax,TVirtualX::kFilled);
   gVirtualX->DrawBox(px1,pymin,px2,pymax,TVirtualX::kFilled);
   pxold1 = px1;
   pxold2 = px2;

   // Create or set the new canvas proj y
   TVirtualPad *padsav = gPad;
   TVirtualPad *c = (TVirtualPad*)gROOT->GetListOfCanvases()->FindObject(Form("c_%lx_projection_%d",
                                                                              (ULong_t)fH, fShowProjection));
   if (c) {
      c->Clear();
   } else {
      fShowProjection = 0;
      pxold1 = 0;
      pxold2 = 0;
      return;
   }
   c->cd();
   c->SetLogy(padsav->GetLogz());
   c->SetLogx(padsav->GetLogy());

   // Draw slice corresponding to mouse position
   TString prjName = TString::Format("slice_py_of_%s",fH->GetName());
   TH1D *hp = ((TH2*)fH)->ProjectionY(prjName, binx1, binx2);
   if (hp) {
      hp->SetFillColor(38);
      // apply a patch from Oliver Freyermuth to set the title in the projection
      // using the range of the projected X values
      if (binx1 == binx2) {
         Double_t valueFrom   = fH->GetXaxis()->GetBinLowEdge(binx1);
         Double_t valueTo     = fH->GetXaxis()->GetBinUpEdge(binx1);
         // Limit precision to 1 digit more than the difference between upper and lower bound (to also catch 121.5-120.5).
         Int_t valuePrecision = -TMath::Nint(TMath::Log10(valueTo-valueFrom))+1;
         if (fH->GetXaxis()->GetLabels() != NULL) {
            hp->SetTitle(TString::Format("ProjectionY of binx=%d [x=%.*lf..%.*lf] [%s]", binx1, valuePrecision, valueFrom, valuePrecision, valueTo, fH->GetXaxis()->GetBinLabel(binx1)));
         } else {
            hp->SetTitle(TString::Format("ProjectionY of binx=%d [x=%.*lf..%.*lf]", binx1, valuePrecision, valueFrom, valuePrecision, valueTo));
         }
      } else {
         Double_t valueFrom   = fH->GetXaxis()->GetBinLowEdge(binx1);
         Double_t valueTo     = fH->GetXaxis()->GetBinUpEdge(binx2);
         // Limit precision to 1 digit more than the difference between upper and lower bound (to also catch 121.5-120.5).
         // binx1 is used here to get equal precision no matter how large the binrange is,
         // otherwise precision may change when moving the mouse to the histogram boundaries (limiting effective binrange).
         Int_t valuePrecision = -TMath::Nint(TMath::Log10(fH->GetXaxis()->GetBinUpEdge(binx1)-valueFrom))+1;
         if (fH->GetXaxis()->GetLabels() != NULL) {
            hp->SetTitle(TString::Format("ProjectionY of binx=[%d,%d] [x=%.*lf..%.*lf] [%s..%s]", binx1, binx2, valuePrecision, valueFrom, valuePrecision, valueTo, fH->GetXaxis()->GetBinLabel(binx1), fH->GetXaxis()->GetBinLabel(binx2)));
         } else {
            hp->SetTitle(TString::Format("ProjectionY of binx=[%d,%d] [x=%.*lf..%.*lf]", binx1, binx2, valuePrecision, valueFrom, valuePrecision, valueTo));
         }
      }
      hp->SetXTitle(fH->GetYaxis()->GetTitle());
      hp->SetYTitle("Number of Entries");
      hp->Draw();
      c->Update();
      padsav->cd();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Show projection (specified by `fShowProjection`) of a `TH3`.
/// The drawing option for the projection is in `fShowOption`.
///
/// First implementation; R.Brun
///
/// Full implementation: Tim Tran (timtran@jlab.org)  April 2006

void THistPainter::ShowProjection3(Int_t px, Int_t py)
{

   Int_t nbins=(Int_t)fShowProjection/100; //decode nbins
   if (fH->GetDimension() < 3) {
      if (fShowProjection%100 == 1) {ShowProjectionX(px,py); return;}
      if (fShowProjection%100 == 2) {ShowProjectionY(px,py); return;}
   }

   gPad->SetDoubleBuffer(0);             // turn off double buffer mode
   gVirtualX->SetDrawMode(TVirtualX::kInvert);  // set the drawing mode to XOR mode

   // Erase old position and draw a line at current position
   TView *view = gPad->GetView();
   if (!view) return;
   TH3 *h3 = (TH3*)fH;
   TAxis *xaxis = h3->GetXaxis();
   TAxis *yaxis = h3->GetYaxis();
   TAxis *zaxis = h3->GetZaxis();
   Double_t u[3],xx[3];

   static TPoint line1[2];//store end points of a line, initialised 0 by default
   static TPoint line2[2];// second line when slice thickness > 1 bin thickness
   static TPoint line3[2];
   static TPoint line4[2];
   static TPoint endface1[5];
   static TPoint endface2[5];
   static TPoint rect1[5];//store vertices of the polyline (rectangle), initialsed 0 by default
   static TPoint rect2[5];// second rectangle when slice thickness > 1 bin thickness

   Double_t uxmin = gPad->GetUxmin();
   Double_t uxmax = gPad->GetUxmax();
   Double_t uymin = gPad->GetUymin();
   Double_t uymax = gPad->GetUymax();

   int pxmin = gPad->XtoAbsPixel(uxmin);
   int pxmax = gPad->XtoAbsPixel(uxmax);
   if (pxmin==pxmax) return;
   int pymin = gPad->YtoAbsPixel(uymin);
   int pymax = gPad->YtoAbsPixel(uymax);
   if (pymin==pymax) return;
   Double_t cx    = (pxmax-pxmin)/(uxmax-uxmin);
   Double_t cy    = (pymax-pymin)/(uymax-uymin);
   TVirtualPad *padsav = gPad;
   TVirtualPad *c = (TVirtualPad*)gROOT->GetListOfCanvases()->FindObject(Form("c_%lx_projection_%d",
                                                                              (ULong_t)fH, fShowProjection));
   if (!c) {
      fShowProjection = 0;
      return;
   }

   switch ((Int_t)fShowProjection%100) {
      case 1:
         // "x"
         {
            Int_t firstY = yaxis->GetFirst();
            Int_t lastY  = yaxis->GetLast();
            Int_t biny = firstY + Int_t((lastY-firstY)*(px-pxmin)/(pxmax-pxmin));
            Int_t biny2 = TMath::Min(biny+nbins-1,yaxis->GetNbins() );
            yaxis->SetRange(biny,biny2);
            Int_t firstZ = zaxis->GetFirst();
            Int_t lastZ  = zaxis->GetLast();
            Int_t binz = firstZ + Int_t((lastZ-firstZ)*(py-pymin)/(pymax-pymin));
            Int_t binz2 = TMath::Min(binz+nbins-1,zaxis->GetNbins() );
            zaxis->SetRange(binz,binz2);
            if (line1[0].GetX()) gVirtualX->DrawPolyLine(2,line1);
            if (nbins>1 && line1[0].GetX()) {
               gVirtualX->DrawPolyLine(2,line2);
               gVirtualX->DrawPolyLine(2,line3);
               gVirtualX->DrawPolyLine(2,line4);
               gVirtualX->DrawPolyLine(5,endface1);
               gVirtualX->DrawPolyLine(5,endface2);
            }
            xx[0] = xaxis->GetXmin();
            xx[2] = zaxis->GetBinCenter(binz);
            xx[1] = yaxis->GetBinCenter(biny);
            view->WCtoNDC(xx,u);
            line1[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            line1[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
            xx[0] = xaxis->GetXmax();
            view->WCtoNDC(xx,u);
            line1[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            line1[1].SetY(pymin + Int_t((u[1]-uymin)*cy));
            gVirtualX->DrawPolyLine(2,line1);
            if (nbins>1) {
               xx[0] = xaxis->GetXmin();
               xx[2] = zaxis->GetBinCenter(binz+nbins-1);
               xx[1] = yaxis->GetBinCenter(biny);
               view->WCtoNDC(xx,u);
               line2[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               line2[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[0] = xaxis->GetXmax();
               view->WCtoNDC(xx,u);
               line2[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               line2[1].SetY(pymin + Int_t((u[1]-uymin)*cy));

               xx[0] = xaxis->GetXmin();
               xx[2] = zaxis->GetBinCenter(binz+nbins-1);
               xx[1] = yaxis->GetBinCenter(biny+nbins-1);
               view->WCtoNDC(xx,u);
               line3[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               line3[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[0] = xaxis->GetXmax();
               view->WCtoNDC(xx,u);
               line3[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               line3[1].SetY(pymin + Int_t((u[1]-uymin)*cy));

               xx[0] = xaxis->GetXmin();
               xx[2] = zaxis->GetBinCenter(binz);
               xx[1] = yaxis->GetBinCenter(biny+nbins-1);
               view->WCtoNDC(xx,u);
               line4[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               line4[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[0] = xaxis->GetXmax();
               view->WCtoNDC(xx,u);
               line4[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               line4[1].SetY(pymin + Int_t((u[1]-uymin)*cy));

               endface1[0].SetX(line1[0].GetX());
               endface1[0].SetY(line1[0].GetY());
               endface1[1].SetX(line2[0].GetX());
               endface1[1].SetY(line2[0].GetY());
               endface1[2].SetX(line3[0].GetX());
               endface1[2].SetY(line3[0].GetY());
               endface1[3].SetX(line4[0].GetX());
               endface1[3].SetY(line4[0].GetY());
               endface1[4].SetX(line1[0].GetX());
               endface1[4].SetY(line1[0].GetY());

               endface2[0].SetX(line1[1].GetX());
               endface2[0].SetY(line1[1].GetY());
               endface2[1].SetX(line2[1].GetX());
               endface2[1].SetY(line2[1].GetY());
               endface2[2].SetX(line3[1].GetX());
               endface2[2].SetY(line3[1].GetY());
               endface2[3].SetX(line4[1].GetX());
               endface2[3].SetY(line4[1].GetY());
               endface2[4].SetX(line1[1].GetX());
               endface2[4].SetY(line1[1].GetY());

               gVirtualX->DrawPolyLine(2,line2);
               gVirtualX->DrawPolyLine(2,line3);
               gVirtualX->DrawPolyLine(2,line4);
               gVirtualX->DrawPolyLine(5,endface1);
               gVirtualX->DrawPolyLine(5,endface2);
            }
            c->Clear();
            c->cd();
            TH1 *hp = h3->Project3D("x");
            yaxis->SetRange(firstY,lastY);
            zaxis->SetRange(firstZ,lastZ);
            if (hp) {
               hp->SetFillColor(38);
               if (nbins == 1)
                  hp->SetTitle(TString::Format("ProjectionX of biny=%d [y=%.1f..%.1f] binz=%d [z=%.1f..%.1f]", biny, yaxis->GetBinLowEdge(biny), yaxis->GetBinUpEdge(biny),
                                               binz, zaxis->GetBinLowEdge(binz), zaxis->GetBinUpEdge(binz)));
               else {
                  hp->SetTitle(TString::Format("ProjectionX, biny=[%d,%d] [y=%.1f..%.1f], binz=[%d,%d] [z=%.1f..%.1f]", biny, biny2, yaxis->GetBinLowEdge(biny), yaxis->GetBinUpEdge(biny2),
                                               binz, binz2, zaxis->GetBinLowEdge(binz), zaxis->GetBinUpEdge(binz2) ) );
               }
               hp->SetXTitle(fH->GetXaxis()->GetTitle());
               hp->SetYTitle("Number of Entries");
               hp->Draw(fShowOption.Data());
            }
         }
         break;

      case 2:
         // "y"
         {
            Int_t firstX = xaxis->GetFirst();
            Int_t lastX  = xaxis->GetLast();
            Int_t binx = firstX + Int_t((lastX-firstX)*(px-pxmin)/(pxmax-pxmin));
            Int_t binx2 = TMath::Min(binx+nbins-1,xaxis->GetNbins() );
            xaxis->SetRange(binx,binx2);
            Int_t firstZ = zaxis->GetFirst();
            Int_t lastZ  = zaxis->GetLast();
            Int_t binz = firstZ + Int_t((lastZ-firstZ)*(py-pymin)/(pymax-pymin));
            Int_t binz2 = TMath::Min(binz+nbins-1,zaxis->GetNbins() );
            zaxis->SetRange(binz,binz2);
            if (line1[0].GetX()) gVirtualX->DrawPolyLine(2,line1);
            if (nbins>1 && line1[0].GetX()) {
               gVirtualX->DrawPolyLine(2,line2);
               gVirtualX->DrawPolyLine(2,line3);
               gVirtualX->DrawPolyLine(2,line4);
               gVirtualX->DrawPolyLine(5,endface1);
               gVirtualX->DrawPolyLine(5,endface2);
            }
            xx[0]=xaxis->GetBinCenter(binx);
            xx[2] = zaxis->GetBinCenter(binz);
            xx[1] = yaxis->GetXmin();
            view->WCtoNDC(xx,u);
            line1[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            line1[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
            xx[1] = yaxis->GetXmax();
            view->WCtoNDC(xx,u);
            line1[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            line1[1].SetY(pymin + Int_t((u[1]-uymin)*cy));
            gVirtualX->DrawPolyLine(2,line1);
            if (nbins>1) {
               xx[1] = yaxis->GetXmin();
               xx[2] = zaxis->GetBinCenter(binz+nbins-1);
               xx[0] = xaxis->GetBinCenter(binx);
               view->WCtoNDC(xx,u);
               line2[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               line2[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[1] = yaxis->GetXmax();
               view->WCtoNDC(xx,u);
               line2[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               line2[1].SetY(pymin + Int_t((u[1]-uymin)*cy));

               xx[1] = yaxis->GetXmin();
               xx[2] = zaxis->GetBinCenter(binz+nbins-1);
               xx[0] = xaxis->GetBinCenter(binx+nbins-1);
               view->WCtoNDC(xx,u);
               line3[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               line3[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[1] = yaxis->GetXmax();
               view->WCtoNDC(xx,u);
               line3[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               line3[1].SetY(pymin + Int_t((u[1]-uymin)*cy));

               xx[1] = yaxis->GetXmin();
               xx[2] = zaxis->GetBinCenter(binz);
               xx[0] = xaxis->GetBinCenter(binx+nbins-1);
               view->WCtoNDC(xx,u);
               line4[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               line4[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[1] = yaxis->GetXmax();
               view->WCtoNDC(xx,u);
               line4[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               line4[1].SetY(pymin + Int_t((u[1]-uymin)*cy));

               endface1[0].SetX(line1[0].GetX());
               endface1[0].SetY(line1[0].GetY());
               endface1[1].SetX(line2[0].GetX());
               endface1[1].SetY(line2[0].GetY());
               endface1[2].SetX(line3[0].GetX());
               endface1[2].SetY(line3[0].GetY());
               endface1[3].SetX(line4[0].GetX());
               endface1[3].SetY(line4[0].GetY());
               endface1[4].SetX(line1[0].GetX());
               endface1[4].SetY(line1[0].GetY());

               endface2[0].SetX(line1[1].GetX());
               endface2[0].SetY(line1[1].GetY());
               endface2[1].SetX(line2[1].GetX());
               endface2[1].SetY(line2[1].GetY());
               endface2[2].SetX(line3[1].GetX());
               endface2[2].SetY(line3[1].GetY());
               endface2[3].SetX(line4[1].GetX());
               endface2[3].SetY(line4[1].GetY());
               endface2[4].SetX(line1[1].GetX());
               endface2[4].SetY(line1[1].GetY());

               gVirtualX->DrawPolyLine(2,line2);
               gVirtualX->DrawPolyLine(2,line3);
               gVirtualX->DrawPolyLine(2,line4);
               gVirtualX->DrawPolyLine(5,endface1);
               gVirtualX->DrawPolyLine(5,endface2);
            }
            c->Clear();
            c->cd();
            TH1 *hp = h3->Project3D("y");
            xaxis->SetRange(firstX,lastX);
            zaxis->SetRange(firstZ,lastZ);
            if (hp) {
               hp->SetFillColor(38);
               if (nbins == 1)
                  hp->SetTitle(TString::Format("ProjectionY of binx=%d [x=%.1f..%.1f] binz=%d [z=%.1f..%.1f]", binx, xaxis->GetBinLowEdge(binx), xaxis->GetBinUpEdge(binx),
                                               binz, zaxis->GetBinLowEdge(binz), zaxis->GetBinUpEdge(binz)));
               else
                  hp->SetTitle(TString::Format("ProjectionY, binx=[%d,%d] [x=%.1f..%.1f], binz=[%d,%d] [z=%.1f..%.1f]", binx, binx2, xaxis->GetBinLowEdge(binx), xaxis->GetBinUpEdge(binx2),
                                               binz, binz2, zaxis->GetBinLowEdge(binz), zaxis->GetBinUpEdge(binz2) ) );
               hp->SetXTitle(fH->GetYaxis()->GetTitle());
               hp->SetYTitle("Number of Entries");
               hp->Draw(fShowOption.Data());
            }
         }
         break;

      case 3:
         // "z"
         {
            Int_t firstX = xaxis->GetFirst();
            Int_t lastX  = xaxis->GetLast();
            Int_t binx = firstX + Int_t((lastX-firstX)*(px-pxmin)/(pxmax-pxmin));
            Int_t binx2 = TMath::Min(binx+nbins-1,xaxis->GetNbins() );
            xaxis->SetRange(binx,binx2);
            Int_t firstY = yaxis->GetFirst();
            Int_t lastY  = yaxis->GetLast();
            Int_t biny = firstY + Int_t((lastY-firstY)*(py-pymin)/(pymax-pymin));
            Int_t biny2 = TMath::Min(biny+nbins-1,yaxis->GetNbins() );
            yaxis->SetRange(biny,biny2);
            if (line1[0].GetX()) gVirtualX->DrawPolyLine(2,line1);
            if (nbins>1 && line1[0].GetX()) {
               gVirtualX->DrawPolyLine(2,line2);
               gVirtualX->DrawPolyLine(2,line3);
               gVirtualX->DrawPolyLine(2,line4);
               gVirtualX->DrawPolyLine(5,endface1);
               gVirtualX->DrawPolyLine(5,endface2);
            }
            xx[0] = xaxis->GetBinCenter(binx);
            xx[1] = yaxis->GetBinCenter(biny);
            xx[2] = zaxis->GetXmin();
            view->WCtoNDC(xx,u);
            line1[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            line1[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
            xx[2] = zaxis->GetXmax();
            view->WCtoNDC(xx,u);
            line1[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            line1[1].SetY(pymin + Int_t((u[1]-uymin)*cy));
            gVirtualX->DrawPolyLine(2,line1);
            if (nbins>1) {
               xx[2] = zaxis->GetXmin();
               xx[1] = yaxis->GetBinCenter(biny+nbins-1);
               xx[0] = xaxis->GetBinCenter(binx);
               view->WCtoNDC(xx,u);
               line2[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               line2[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[2] = zaxis->GetXmax();
               view->WCtoNDC(xx,u);
               line2[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               line2[1].SetY(pymin + Int_t((u[1]-uymin)*cy));

               xx[2] = zaxis->GetXmin();
               xx[1] = yaxis->GetBinCenter(biny+nbins-1);
               xx[0] = xaxis->GetBinCenter(binx+nbins-1);
               view->WCtoNDC(xx,u);
               line3[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               line3[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[2] = zaxis->GetXmax();
               view->WCtoNDC(xx,u);
               line3[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               line3[1].SetY(pymin + Int_t((u[1]-uymin)*cy));

               xx[2] = zaxis->GetXmin();
               xx[1] = yaxis->GetBinCenter(biny);
               xx[0] = xaxis->GetBinCenter(binx+nbins-1);
               view->WCtoNDC(xx,u);
               line4[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               line4[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[2] = zaxis->GetXmax();
               view->WCtoNDC(xx,u);
               line4[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               line4[1].SetY(pymin + Int_t((u[1]-uymin)*cy));

               endface1[0].SetX(line1[0].GetX());
               endface1[0].SetY(line1[0].GetY());
               endface1[1].SetX(line2[0].GetX());
               endface1[1].SetY(line2[0].GetY());
               endface1[2].SetX(line3[0].GetX());
               endface1[2].SetY(line3[0].GetY());
               endface1[3].SetX(line4[0].GetX());
               endface1[3].SetY(line4[0].GetY());
               endface1[4].SetX(line1[0].GetX());
               endface1[4].SetY(line1[0].GetY());

               endface2[0].SetX(line1[1].GetX());
               endface2[0].SetY(line1[1].GetY());
               endface2[1].SetX(line2[1].GetX());
               endface2[1].SetY(line2[1].GetY());
               endface2[2].SetX(line3[1].GetX());
               endface2[2].SetY(line3[1].GetY());
               endface2[3].SetX(line4[1].GetX());
               endface2[3].SetY(line4[1].GetY());
               endface2[4].SetX(line1[1].GetX());
               endface2[4].SetY(line1[1].GetY());

               gVirtualX->DrawPolyLine(2,line2);
               gVirtualX->DrawPolyLine(2,line3);
               gVirtualX->DrawPolyLine(2,line4);
               gVirtualX->DrawPolyLine(5,endface1);
               gVirtualX->DrawPolyLine(5,endface2);
            }
            c->Clear();
            c->cd();
            TH1 *hp = h3->Project3D("z");
            xaxis->SetRange(firstX,lastX);
            yaxis->SetRange(firstY,lastY);
            if (hp) {
               hp->SetFillColor(38);
               if (nbins == 1)
                  hp->SetTitle(TString::Format("ProjectionZ of binx=%d [x=%.1f..%.1f] biny=%d [y=%.1f..%.1f]", binx, xaxis->GetBinLowEdge(binx), xaxis->GetBinUpEdge(binx),
                                               biny, yaxis->GetBinLowEdge(biny), yaxis->GetBinUpEdge(biny)));
               else
                  hp->SetTitle(TString::Format("ProjectionZ, binx=[%d,%d] [x=%.1f..%.1f], biny=[%d,%d] [y=%.1f..%.1f]", binx, binx2, xaxis->GetBinLowEdge(binx), xaxis->GetBinUpEdge(binx2),
                                               biny, biny2, yaxis->GetBinLowEdge(biny), yaxis->GetBinUpEdge(biny2) ) );
               hp->SetXTitle(fH->GetZaxis()->GetTitle());
               hp->SetYTitle("Number of Entries");
               hp->Draw(fShowOption.Data());
            }
         }
         break;

      case 4:
         // "xy"
         {
            Int_t first = zaxis->GetFirst();
            Int_t last  = zaxis->GetLast();
            Int_t binz  = first + Int_t((last-first)*(py-pymin)/(pymax-pymin));
            Int_t binz2 = TMath::Min(binz+nbins-1,zaxis->GetNbins() );
            zaxis->SetRange(binz,binz2);
            if (rect1[0].GetX())            gVirtualX->DrawPolyLine(5,rect1);
            if (nbins>1 && rect2[0].GetX()) gVirtualX->DrawPolyLine(5,rect2);
            xx[0] = xaxis->GetXmin();
            xx[1] = yaxis->GetXmax();
            xx[2] = zaxis->GetBinCenter(binz);
            view->WCtoNDC(xx,u);
            rect1[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
            rect1[4].SetX(rect1[0].GetX());
            rect1[4].SetY(rect1[0].GetY());
            xx[0] = xaxis->GetXmax();
            view->WCtoNDC(xx,u);
            rect1[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[1].SetY(pymin + Int_t((u[1]-uymin)*cy));
            xx[1] = yaxis->GetXmin();
            view->WCtoNDC(xx,u);
            rect1[2].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[2].SetY(pymin + Int_t((u[1]-uymin)*cy));
            xx[0] = xaxis->GetXmin();
            view->WCtoNDC(xx,u);
            rect1[3].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[3].SetY(pymin + Int_t((u[1]-uymin)*cy));
            gVirtualX->DrawPolyLine(5,rect1);
            if (nbins>1) {
               xx[0] = xaxis->GetXmin();
               xx[1] = yaxis->GetXmax();
               xx[2] = zaxis->GetBinCenter(binz+nbins-1);
               view->WCtoNDC(xx,u);
               rect2[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
               rect2[4].SetX(rect2[0].GetX());
               rect2[4].SetY(rect2[0].GetY());
               xx[0] = xaxis->GetXmax();
               view->WCtoNDC(xx,u);
               rect2[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[1].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[1] = yaxis->GetXmin();
               view->WCtoNDC(xx,u);
               rect2[2].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[2].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[0] = xaxis->GetXmin();
               view->WCtoNDC(xx,u);
               rect2[3].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[3].SetY(pymin + Int_t((u[1]-uymin)*cy));
               gVirtualX->DrawPolyLine(5,rect2);
            }

            c->Clear();
            c->cd();
            TH2 *hp = (TH2*)h3->Project3D("xy");
            zaxis->SetRange(first,last);
            if (hp) {
               hp->SetFillColor(38);
               if (nbins==1)hp->SetTitle(TString::Format("ProjectionXY of binz=%d [z=%.1f..%.f]", binz,zaxis->GetBinLowEdge(binz),zaxis->GetBinUpEdge(binz)));
               else        hp->SetTitle(TString::Format("ProjectionXY, binz=[%d,%d] [z=%.1f..%.1f]", binz,binz2,zaxis->GetBinLowEdge(binz),zaxis->GetBinUpEdge(binz2)));
               hp->SetXTitle(fH->GetYaxis()->GetTitle());
               hp->SetYTitle(fH->GetXaxis()->GetTitle());
               hp->SetZTitle("Number of Entries");
               hp->Draw(fShowOption.Data());
            }
         }
         break;

      case 5:
         // "yx"
         {
            Int_t first = zaxis->GetFirst();
            Int_t last  = zaxis->GetLast();
            Int_t binz = first + Int_t((last-first)*(py-pymin)/(pymax-pymin));
            Int_t binz2 = TMath::Min(binz+nbins-1,zaxis->GetNbins() );
            zaxis->SetRange(binz,binz2);
            if (rect1[0].GetX())            gVirtualX->DrawPolyLine(5,rect1);
            if (nbins>1 && rect2[0].GetX()) gVirtualX->DrawPolyLine(5,rect2);
            xx[0] = xaxis->GetXmin();
            xx[1] = yaxis->GetXmax();
            xx[2] = zaxis->GetBinCenter(binz);
            view->WCtoNDC(xx,u);
            rect1[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
            rect1[4].SetX(rect1[0].GetX());
            rect1[4].SetY(rect1[0].GetY());
            xx[0] = xaxis->GetXmax();
            view->WCtoNDC(xx,u);
            rect1[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[1].SetY(pymin + Int_t((u[1]-uymin)*cy));
            xx[1] = yaxis->GetXmin();
            view->WCtoNDC(xx,u);
            rect1[2].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[2].SetY(pymin + Int_t((u[1]-uymin)*cy));
            xx[0] = xaxis->GetXmin();
            view->WCtoNDC(xx,u);
            rect1[3].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[3].SetY(pymin + Int_t((u[1]-uymin)*cy));
            gVirtualX->DrawPolyLine(5,rect1);
            if (nbins>1) {
               xx[0] = xaxis->GetXmin();
               xx[1] = yaxis->GetXmax();
               xx[2] = zaxis->GetBinCenter(binz+nbins-1);
               view->WCtoNDC(xx,u);
               rect2[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
               rect2[4].SetX(rect2[0].GetX());
               rect2[4].SetY(rect2[0].GetY());
               xx[0] = xaxis->GetXmax();
               view->WCtoNDC(xx,u);
               rect2[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[1].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[1] = yaxis->GetXmin();
               view->WCtoNDC(xx,u);
               rect2[2].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[2].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[0] = xaxis->GetXmin();
               view->WCtoNDC(xx,u);
               rect2[3].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[3].SetY(pymin + Int_t((u[1]-uymin)*cy));
               gVirtualX->DrawPolyLine(5,rect2);
            }
            c->Clear();
            c->cd();
            TH2 *hp = (TH2*)h3->Project3D("yx");
            zaxis->SetRange(first,last);
            if (hp) {
               hp->SetFillColor(38);
               if (nbins==1)hp->SetTitle(TString::Format("ProjectionYX of binz=%d [z=%.1f..%.f]", binz,zaxis->GetBinLowEdge(binz),zaxis->GetBinUpEdge(binz)));
               else        hp->SetTitle(TString::Format("ProjectionYX, binz=[%d,%d] [z=%.1f..%.1f]", binz,binz2,zaxis->GetBinLowEdge(binz),zaxis->GetBinUpEdge(binz2)));
               hp->SetXTitle(fH->GetXaxis()->GetTitle());
               hp->SetYTitle(fH->GetYaxis()->GetTitle());
               hp->SetZTitle("Number of Entries");
               hp->Draw(fShowOption.Data());
            }
         }
         break;

      case 6:
         // "xz"
         {
            Int_t first = yaxis->GetFirst();
            Int_t last  = yaxis->GetLast();
            Int_t biny = first + Int_t((last-first)*(py-pymin)/(pymax-pymin));
            Int_t biny2 = TMath::Min(biny+nbins-1,yaxis->GetNbins() );
            yaxis->SetRange(biny,biny2);
            if (rect1[0].GetX())            gVirtualX->DrawPolyLine(5,rect1);
            if (nbins>1 && rect1[0].GetX()) gVirtualX->DrawPolyLine(5,rect2);
            xx[0] = xaxis->GetXmin();
            xx[2] = zaxis->GetXmax();
            xx[1] = yaxis->GetBinCenter(biny);
            view->WCtoNDC(xx,u);
            rect1[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
            rect1[4].SetX(rect1[0].GetX());
            rect1[4].SetY(rect1[0].GetY());
            xx[0] = xaxis->GetXmax();
            view->WCtoNDC(xx,u);
            rect1[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[1].SetY(pymin + Int_t((u[1]-uymin)*cy));
            xx[2] = zaxis->GetXmin();
            view->WCtoNDC(xx,u);
            rect1[2].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[2].SetY(pymin + Int_t((u[1]-uymin)*cy));
            xx[0] = xaxis->GetXmin();
            view->WCtoNDC(xx,u);
            rect1[3].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[3].SetY(pymin + Int_t((u[1]-uymin)*cy));
            gVirtualX->DrawPolyLine(5,rect1);
            if (nbins>1) {
               xx[0] = xaxis->GetXmin();
               xx[2] = zaxis->GetXmax();
               xx[1] = yaxis->GetBinCenter(biny+nbins-1);
               view->WCtoNDC(xx,u);
               rect2[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
               rect2[4].SetX(rect2[0].GetX());
               rect2[4].SetY(rect2[0].GetY());
               xx[0] = xaxis->GetXmax();
               view->WCtoNDC(xx,u);
               rect2[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[1].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[2] = zaxis->GetXmin();
               view->WCtoNDC(xx,u);
               rect2[2].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[2].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[0] = xaxis->GetXmin();
               view->WCtoNDC(xx,u);
               rect2[3].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[3].SetY(pymin + Int_t((u[1]-uymin)*cy));
               gVirtualX->DrawPolyLine(5,rect2);
            }
            c->Clear();
            c->cd();
            TH2 *hp = (TH2*)h3->Project3D("xz");
            yaxis->SetRange(first,last);
            if (hp) {
               hp->SetFillColor(38);
               if (nbins==1)hp->SetTitle(TString::Format("ProjectionXZ of biny=%d [y=%.1f..%.f]", biny,yaxis->GetBinLowEdge(biny),yaxis->GetBinUpEdge(biny)));
               else        hp->SetTitle(TString::Format("ProjectionXZ, biny=[%d,%d] [y=%.1f..%.1f]", biny,biny2,yaxis->GetBinLowEdge(biny),yaxis->GetBinUpEdge(biny2)));
               hp->SetXTitle(fH->GetZaxis()->GetTitle());
               hp->SetYTitle(fH->GetXaxis()->GetTitle());
               hp->SetZTitle("Number of Entries");
               hp->Draw(fShowOption.Data());
            }
         }
         break;

      case 7:
         // "zx"
         {
            Int_t first = yaxis->GetFirst();
            Int_t last  = yaxis->GetLast();
            Int_t biny = first + Int_t((last-first)*(py-pymin)/(pymax-pymin));
            Int_t biny2 = TMath::Min(biny+nbins-1,yaxis->GetNbins() );
            yaxis->SetRange(biny,biny2);
            if (rect1[0].GetX())            gVirtualX->DrawPolyLine(5,rect1);
            if (nbins>1 && rect1[0].GetX()) gVirtualX->DrawPolyLine(5,rect2);
            xx[0] = xaxis->GetXmin();
            xx[2] = zaxis->GetXmax();
            xx[1] = yaxis->GetBinCenter(biny);
            view->WCtoNDC(xx,u);
            rect1[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
            rect1[4].SetX(rect1[0].GetX());
            rect1[4].SetY(rect1[0].GetY());
            xx[0] = xaxis->GetXmax();
            view->WCtoNDC(xx,u);
            rect1[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[1].SetY(pymin + Int_t((u[1]-uymin)*cy));
            xx[2] = zaxis->GetXmin();
            view->WCtoNDC(xx,u);
            rect1[2].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[2].SetY(pymin + Int_t((u[1]-uymin)*cy));
            xx[0] = xaxis->GetXmin();
            view->WCtoNDC(xx,u);
            rect1[3].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[3].SetY(pymin + Int_t((u[1]-uymin)*cy));
            gVirtualX->DrawPolyLine(5,rect1);
            if (nbins>1) {
               xx[0] = xaxis->GetXmin();
               xx[2] = zaxis->GetXmax();
               xx[1] = yaxis->GetBinCenter(biny+nbins-1);
               view->WCtoNDC(xx,u);
               rect2[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
               rect2[4].SetX(rect2[0].GetX());
               rect2[4].SetY(rect2[0].GetY());
               xx[0] = xaxis->GetXmax();
               view->WCtoNDC(xx,u);
               rect2[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[1].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[2] = zaxis->GetXmin();
               view->WCtoNDC(xx,u);
               rect2[2].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[2].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[0] = xaxis->GetXmin();
               view->WCtoNDC(xx,u);
               rect2[3].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[3].SetY(pymin + Int_t((u[1]-uymin)*cy));
               gVirtualX->DrawPolyLine(5,rect2);
            }
            c->Clear();
            c->cd();
            TH2 *hp = (TH2*)h3->Project3D("zx");
            yaxis->SetRange(first,last);
            if (hp) {
               hp->SetFillColor(38);
               if (nbins==1)hp->SetTitle(TString::Format("ProjectionZX of biny=%d [y=%.1f..%.f]", biny,yaxis->GetBinLowEdge(biny),yaxis->GetBinUpEdge(biny)));
               else        hp->SetTitle(TString::Format("ProjectionZX, biny=[%d,%d] [y=%.1f..%.1f]", biny,biny2,yaxis->GetBinLowEdge(biny),yaxis->GetBinUpEdge(biny2)));
               hp->SetXTitle(fH->GetXaxis()->GetTitle());
               hp->SetYTitle(fH->GetZaxis()->GetTitle());
               hp->SetZTitle("Number of Entries");
               hp->Draw(fShowOption.Data());
            }
         }
         break;

      case 8:
         // "yz"
         {
            Int_t first = xaxis->GetFirst();
            Int_t last  = xaxis->GetLast();
            Int_t binx = first + Int_t((last-first)*(px-pxmin)/(pxmax-pxmin));
            Int_t binx2 = TMath::Min(binx+nbins-1,xaxis->GetNbins() );
            xaxis->SetRange(binx,binx2);
            if (rect1[0].GetX()) gVirtualX->DrawPolyLine(5,rect1);
            if (nbins>1 && rect1[0].GetX()) gVirtualX->DrawPolyLine(5,rect2);
            xx[2] = zaxis->GetXmin();
            xx[1] = yaxis->GetXmax();
            xx[0] = xaxis->GetBinCenter(binx);
            view->WCtoNDC(xx,u);
            rect1[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
            rect1[4].SetX(rect1[0].GetX());
            rect1[4].SetY(rect1[0].GetY());
            xx[2] = zaxis->GetXmax();
            view->WCtoNDC(xx,u);
            rect1[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[1].SetY(pymin + Int_t((u[1]-uymin)*cy));
            xx[1] = yaxis->GetXmin();
            view->WCtoNDC(xx,u);
            rect1[2].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[2].SetY(pymin + Int_t((u[1]-uymin)*cy));
            xx[2] = zaxis->GetXmin();
            view->WCtoNDC(xx,u);
            rect1[3].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[3].SetY(pymin + Int_t((u[1]-uymin)*cy));
            gVirtualX->DrawPolyLine(5,rect1);
            if (nbins>1) {
               xx[2] = zaxis->GetXmin();
               xx[1] = yaxis->GetXmax();
               xx[0] = xaxis->GetBinCenter(binx+nbins-1);
               view->WCtoNDC(xx,u);
               rect2[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
               rect2[4].SetX(rect2[0].GetX());
               rect2[4].SetY(rect2[0].GetY());
               xx[2] = zaxis->GetXmax();
               view->WCtoNDC(xx,u);
               rect2[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[1].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[1] = yaxis->GetXmin();
               view->WCtoNDC(xx,u);
               rect2[2].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[2].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[2] = zaxis->GetXmin();
               view->WCtoNDC(xx,u);
               rect2[3].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[3].SetY(pymin + Int_t((u[1]-uymin)*cy));
               gVirtualX->DrawPolyLine(5,rect2);
            }
            c->Clear();
            c->cd();
            TH2 *hp = (TH2*)h3->Project3D("yz");
            xaxis->SetRange(first,last);
            if (hp) {
               hp->SetFillColor(38);
               if (nbins==1)hp->SetTitle(TString::Format("ProjectionYZ of binx=%d [x=%.1f..%.f]", binx,xaxis->GetBinLowEdge(binx),xaxis->GetBinUpEdge(binx)));
               else         hp->SetTitle(TString::Format("ProjectionYZ, binx=[%d,%d] [x=%.1f..%.1f]", binx,binx2,xaxis->GetBinLowEdge(binx),xaxis->GetBinUpEdge(binx2)));
               hp->SetXTitle(fH->GetZaxis()->GetTitle());
               hp->SetYTitle(fH->GetYaxis()->GetTitle());
               hp->SetZTitle("Number of Entries");
               hp->Draw(fShowOption.Data());
            }
         }
         break;

      case 9:
         // "zy"
         {
            Int_t first = xaxis->GetFirst();
            Int_t last  = xaxis->GetLast();
            Int_t binx = first + Int_t((last-first)*(px-pxmin)/(pxmax-pxmin));
            Int_t binx2 = TMath::Min(binx+nbins-1,xaxis->GetNbins() );
            xaxis->SetRange(binx,binx2);
            if (rect1[0].GetX()) gVirtualX->DrawPolyLine(5,rect1);
            if (nbins>1 && rect1[0].GetX()) gVirtualX->DrawPolyLine(5,rect2);
            xx[2] = zaxis->GetXmin();
            xx[1] = yaxis->GetXmax();
            xx[0] = xaxis->GetBinCenter(binx);
            view->WCtoNDC(xx,u);
            rect1[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
            rect1[4].SetX(rect1[0].GetX());
            rect1[4].SetY(rect1[0].GetY());
            xx[2] = zaxis->GetXmax();
            view->WCtoNDC(xx,u);
            rect1[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[1].SetY(pymin + Int_t((u[1]-uymin)*cy));
            xx[1] = yaxis->GetXmin();
            view->WCtoNDC(xx,u);
            rect1[2].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[2].SetY(pymin + Int_t((u[1]-uymin)*cy));
            xx[2] = zaxis->GetXmin();
            view->WCtoNDC(xx,u);
            rect1[3].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
            rect1[3].SetY(pymin + Int_t((u[1]-uymin)*cy));
            gVirtualX->DrawPolyLine(5,rect1);
            if (nbins>1) {
               xx[2] = zaxis->GetXmin();
               xx[1] = yaxis->GetXmax();
               xx[0] = xaxis->GetBinCenter(binx+nbins-1);
               view->WCtoNDC(xx,u);
               rect2[0].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[0].SetY(pymin + Int_t((u[1]-uymin)*cy));
               rect2[4].SetX(rect2[0].GetX());
               rect2[4].SetY(rect2[0].GetY());
               xx[2] = zaxis->GetXmax();
               view->WCtoNDC(xx,u);
               rect2[1].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[1].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[1] = yaxis->GetXmin();
               view->WCtoNDC(xx,u);
               rect2[2].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[2].SetY(pymin + Int_t((u[1]-uymin)*cy));
               xx[2] = zaxis->GetXmin();
               view->WCtoNDC(xx,u);
               rect2[3].SetX(pxmin + Int_t((u[0]-uxmin)*cx));
               rect2[3].SetY(pymin + Int_t((u[1]-uymin)*cy));
               gVirtualX->DrawPolyLine(5,rect2);
            }
            c->Clear();
            c->cd();
            TH2 *hp = (TH2*)h3->Project3D("zy");
            xaxis->SetRange(first,last);
            if (hp) {
               hp->SetFillColor(38);
               if (nbins==1)hp->SetTitle(TString::Format("ProjectionZY of binx=%d [x=%.1f..%.f]", binx,xaxis->GetBinLowEdge(binx),xaxis->GetBinUpEdge(binx)));
               else         hp->SetTitle(TString::Format("ProjectionZY, binx=[%d,%d] [x=%.1f..%.1f]", binx,binx2,xaxis->GetBinLowEdge(binx),xaxis->GetBinUpEdge(binx2)));
               hp->SetXTitle(fH->GetYaxis()->GetTitle());
               hp->SetYTitle(fH->GetZaxis()->GetTitle());
               hp->SetZTitle("Number of Entries");
               hp->Draw(fShowOption.Data());
            }
         }
         break;
   }
   c->Update();
   padsav->cd();
}
