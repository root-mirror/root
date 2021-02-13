// @(#)root/postscript:$Id$
// Author: Olivier Couet

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifdef WIN32
#pragma optimize("",off)
#endif

#include <cstdlib>
#include <cstring>
#include <cctype>
#include <fstream>

#include "TROOT.h"
#include "TColor.h"
#include "TVirtualPad.h"
#include "TPoints.h"
#include "TTeXDump.h"
#include "TStyle.h"
#include "TMath.h"

ClassImp(TTeXDump);

/** \class TTeXDump
\ingroup PS

\brief Interface to TeX.

This class allow to generate <b>PGF/TikZ</b> vector graphics output
which can be included in TeX and LaTeX documents.

PGF is a TeX macro package for generating graphics. It is platform
and format-independent and works together with the most important TeX
backend drivers, including pdftex and dvips. It comes with a
user-friendly syntax layer called TikZ.

To generate a such file it is enough to do:
~~~ {.cpp}
   gStyle->SetPaperSize(10.,10.);
   hpx->Draw();
   gPad->Print("hpx.tex");
~~~

Then, the generated file (`hpx.tex`) can be included in a
LaTeX document (`simple.tex`) in the following way:
~~~ {.cpp}
\documentclass{article}
\usepackage{tikz}
\usetikzlibrary{patterns}
\usetikzlibrary{plotmarks}
\title{A simple LaTeX example}
\date{July 2013}
\begin{document}
\maketitle
The following image as been generated using the TTeXDump class:
\par
\input{hpx.tex}
\end{document}
~~~

Note the three directives needed at the top of the LaTeX file:
~~~ {.cpp}
\usepackage{tikz}
\usetikzlibrary{patterns}
\usetikzlibrary{plotmarks}
~~~

Then including the picture in the document is done with the
`\input` directive.

 The command `pdflatex simple.tex` will generate the
corresponding pdf file `simple.pdf`.
*/

////////////////////////////////////////////////////////////////////////////////
/// Default TeX constructor

TTeXDump::TTeXDump() : TVirtualPS()
{
   fStream       = 0;
   fType         = 0;
   gVirtualPS    = this;
   fBoundingBox  = kFALSE;
   fRange        = kFALSE;
   fXsize        = 0.;
   fYsize        = 0.;
   fCurrentRed   = -1.;
   fCurrentGreen = -1.;
   fCurrentBlue  = -1.;
   fCurrentAlpha = 1.;
   fLineScale    = 0.;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the TeX interface
///
///  --fname : TeX file name
///  - wtype : TeX workstation type. Not used in the TeX driver. But as TTeXDump
///            inherits from TVirtualPS it should be kept. Anyway it is not
///            necessary to specify this parameter at creation time because it
///            has a default value (which is ignore in the TeX case).

TTeXDump::TTeXDump(const char *fname, Int_t wtype) : TVirtualPS(fname, wtype)
{
   fStream       = 0;
   fType         = 0;
   gVirtualPS    = this;
   fBoundingBox  = kFALSE;
   fRange        = kFALSE;
   fXsize        = 0.;
   fYsize        = 0.;
   fCurrentRed   = -1.;
   fCurrentGreen = -1.;
   fCurrentBlue  = -1.;
   fCurrentAlpha = 1.;
   fLineScale    = 0.;

   Open(fname, wtype);
}

////////////////////////////////////////////////////////////////////////////////
/// Open a TeX file

void TTeXDump::Open(const char *fname, Int_t wtype)
{
   if (fStream) {
      Warning("Open", "TeX file already open");
      return;
   }

   SetLineScale(gStyle->GetLineScalePS());
   fLenBuffer = 0;
   fType      = abs(wtype);

   gStyle->GetPaperSize(fXsize, fYsize);

   Float_t xrange, yrange;
   if (gPad) {
      Double_t ww = gPad->GetWw();
      Double_t wh = gPad->GetWh();
      ww *= gPad->GetWNDC();
      wh *= gPad->GetHNDC();
      Double_t ratio = wh/ww;
      xrange = fXsize;
      yrange = fXsize*ratio;
      if (yrange > fYsize) { yrange = fYsize; xrange = yrange/ratio;}
      fXsize = xrange; fYsize = yrange;
   }

   // Open OS file
   fStream   = new std::ofstream(fname,std::ios::out);
   if (fStream == 0 || !fStream->good()) {
      printf("ERROR in TTeXDump::Open: Cannot open file:%s\n",fname);
      if (fStream == 0) return;
   }

   gVirtualPS = this;

   for (Int_t i=0;i<fSizBuffer;i++) fBuffer[i] = ' ';

   fBoundingBox = kFALSE;
   fRange       = kFALSE;

   // Set a default range
   Range(fXsize, fYsize);

   NewPage();
}

////////////////////////////////////////////////////////////////////////////////
/// Default TeX destructor

TTeXDump::~TTeXDump()
{
   Close();
}

////////////////////////////////////////////////////////////////////////////////
/// Close a TeX file

void TTeXDump::Close(Option_t *)
{
   if (!gVirtualPS) return;
   if (!fStream) return;
   if (gPad) gPad->Update();
   PrintStr("@");
   PrintStr("\\end{tikzpicture}@");

   // Close file stream
   if (fStream) { fStream->close(); delete fStream; fStream = 0;}

   gVirtualPS = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Activate an already open TeX file

void TTeXDump::On()
{
   // fType is used to know if the TeX file is open. Unlike TPostScript, TTeXDump
   // has no "workstation type". In fact there is only one TeX type.

   if (!fType) {
      Error("On", "no TeX file open");
      Off();
      return;
   }
   gVirtualPS = this;
}

////////////////////////////////////////////////////////////////////////////////
/// Deactivate an already open TeX file

void TTeXDump::Off()
{
   gVirtualPS = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a Box

void TTeXDump::DrawBox(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
   Float_t x1c = XtoTeX(x1);
   Float_t y1c = YtoTeX(y1);
   Float_t x2c = XtoTeX(x2);
   Float_t y2c = YtoTeX(y2);

   Int_t fillis = fFillStyle/1000;
   Int_t fillsi = fFillStyle%1000;

   if (fillis==1) {
      SetColor(fFillColor);
      PrintStr("@");
      PrintStr("\\draw [color=c, fill=c");
      if (fCurrentAlpha != 1.) {
         PrintStr(", fill opacity=");
         WriteReal(fCurrentAlpha, kFALSE);
      }
      PrintStr("] (");
      WriteReal(x1c, kFALSE);
      PrintFast(1,",");
      WriteReal(y1c, kFALSE);
      PrintStr(") rectangle (");
      WriteReal(x2c, kFALSE);
      PrintFast(1,",");
      WriteReal(y2c, kFALSE);
      PrintStr(");");
   }
   if (fillis>1 && fillis<4) {
      SetColor(fFillColor);
      PrintStr("@");
      PrintStr("\\draw [pattern=");
      switch (fillsi) {
      case 1 :
         PrintStr("crosshatch dots");
         break;
      case 2 :
      case 3 :
         PrintStr("dots");
         break;
      case 4 :
         PrintStr("north east lines");
         break;
      case 5 :
         PrintStr("north west lines");
         break;
      case 6 :
         PrintStr("vertical lines");
         break;
      case 7 :
         PrintStr("horizontal lines");
         break;
      case 10 :
          PrintStr("bricks");
         break;
      case 13 :
         PrintStr("crosshatch");
         break;
      }
      PrintStr(", draw=none, pattern color=c");
      if (fCurrentAlpha != 1.) {
         PrintStr(", fill opacity=");
         WriteReal(fCurrentAlpha, kFALSE);
      }
      PrintStr("] (");
      WriteReal(x1c, kFALSE);
      PrintFast(1,",");
      WriteReal(y1c, kFALSE);
      PrintStr(") rectangle (");
      WriteReal(x2c, kFALSE);
      PrintFast(1,",");
      WriteReal(y2c, kFALSE);
      PrintStr(");");
   }
   if (fillis == 0) {
      if (fLineWidth<=0) return;
      SetColor(fLineColor);
      PrintStr("@");
      PrintStr("\\draw [c");
      PrintStr(",line width=");
      WriteReal(0.3*fLineScale*fLineWidth, kFALSE);
      if (fCurrentAlpha != 1.) {
         PrintStr(", opacity=");
         WriteReal(fCurrentAlpha, kFALSE);
      }
      PrintStr("] (");
      WriteReal(x1c, kFALSE);
      PrintFast(1,",");
      WriteReal(y1c, kFALSE);
      PrintStr(") -- (");
      WriteReal(x1c, kFALSE);
      PrintFast(1,",");
      WriteReal(y2c, kFALSE);
      PrintStr(") -- (");
      WriteReal(x2c, kFALSE);
      PrintFast(1,",");
      WriteReal(y2c, kFALSE);
      PrintStr(") -- (");
      WriteReal(x2c, kFALSE);
      PrintFast(1,",");
      WriteReal(y1c, kFALSE);
      PrintStr(") -- (");
      WriteReal(x1c, kFALSE);
      PrintFast(1,",");
      WriteReal(y1c, kFALSE);
      PrintStr(");");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a Frame around a box
///
/// mode = -1  the box looks as it is behind the screen
/// mode =  1  the box looks as it is in front of the screen
/// border is the border size in already pre-computed TeX units dark is the
/// color for the dark part of the frame light is the color for the light
/// part of the frame

void TTeXDump::DrawFrame(Double_t, Double_t, Double_t, Double_t,
                         Int_t, Int_t, Int_t, Int_t)
{
   Warning("DrawFrame", "not yet implemented");
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a PolyLine
///
///  Draw a polyline through  the points  xy.
///  - If NN=1 moves only to point x,y.
///  - If NN=0 the x,y are  written in the TeX file
///       according to the current transformation.
///  - If NN>0 the line is clipped as a line.
///  - If NN<0 the line is clipped as a fill area.

void TTeXDump::DrawPolyLine(Int_t, TPoints *)
{
   Warning("DrawPolyLine", "not yet implemented");
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a PolyLine in NDC space
///
///  Draw a polyline through  the points  xy.
///  - If NN=1 moves only to point x,y.
///  - If NN=0 the x,y are  written in the TeX file
///       according to the current transformation.
///  - If NN>0 the line is clipped as a line.
///  - If NN<0 the line is clipped as a fill area.

void TTeXDump::DrawPolyLineNDC(Int_t, TPoints *)
{
   Warning("DrawPolyLineNDC", "not yet implemented");
}

////////////////////////////////////////////////////////////////////////////////
/// Paint PolyMarker

void TTeXDump::DrawPolyMarker(Int_t, Float_t *, Float_t *)
{
   Warning("DrawPolyMarker", "not yet implemented");
}

////////////////////////////////////////////////////////////////////////////////
/// Paint PolyMarker

void TTeXDump::DrawPolyMarker(Int_t n, Double_t *xw, Double_t *yw)
{
   Float_t x, y;

   SetColor(fMarkerColor);

   PrintStr("@");
   PrintStr("\\foreach \\P in {");

   x = XtoTeX(xw[0]);
   y = YtoTeX(yw[0]);

   PrintStr("(");
   WriteReal(x, kFALSE);
   PrintFast(1,",");
   WriteReal(y, kFALSE);
   PrintStr(")");

   for (Int_t i=1;i<n;i++) {
      x = XtoTeX(xw[i]);
      y = YtoTeX(yw[i]);
      PrintFast(3,", (");
      WriteReal(x, kFALSE);
      PrintFast(1,",");
      WriteReal(y, kFALSE);
      PrintFast(1,")");
   }

   PrintStr("}{\\draw[mark options={color=c,fill=c");

   if (fCurrentAlpha != 1.) {
      PrintStr(",opacity=");
      WriteReal(fCurrentAlpha, kFALSE);
   }

   if (TAttMarker::GetMarkerStyleBase(fMarkerStyle) == 23 || TAttMarker::GetMarkerStyleBase(fMarkerStyle) == 32) PrintStr(",rotate=180");

   PrintStr(Form("},mark size=%fpt", 8./3.33*(fMarkerSize - TMath::Floor(TAttMarker::GetMarkerLineWidth(fMarkerStyle)/2.)/4.)));
   PrintStr(Form(", line width=%fpt", 4./3.33*TMath::Floor(TAttMarker::GetMarkerLineWidth(fMarkerStyle)/2.)));
   PrintStr(", mark=");
   switch (TAttMarker::GetMarkerStyleBase(fMarkerStyle)) {
   case 1 :
      PrintStr("*");
      PrintStr(",mark size=1pt");
      break;
   case 2 :
      PrintStr("+");
      break;
   case 3 :
      PrintStr("asterisk");
      break;
   case 4 :
      PrintStr("o");
      break;
   case 5 :
      PrintStr("x");
      break;
   case 20 :
      PrintStr("*");
      break;
   case 21 :
      PrintStr("square*");
      break;
   case 22 :
      PrintStr("triangle*");
      break;
   case 23 :
      PrintStr("triangle*");
      break;
   case 24 :
      PrintStr("o");
      break;
   case 25 :
      PrintStr("square");
      break;
   case 26 :
      PrintStr("triangle");
      break;
   case 27 :
      PrintStr("diamond");
      break;
   case 28 :
      PrintStr("cross");
      break;
   case 29 :
      PrintStr("newstar*");
      break;
   case 30 :
      PrintStr("newstar");
      break;
   case 31 :
      PrintStr("10-pointed star");
      break;
   case 32 :
      PrintStr("triangle");
      break;
   case 33 :
      PrintStr("diamond*");
      break;
   case 34 :
      PrintStr("cross*");
      break;
   }
   PrintStr("] plot coordinates {\\P};}");
}

////////////////////////////////////////////////////////////////////////////////
/// This function defines a path with xw and yw and draw it according the
/// value of nn:
///
///  - If nn>0 a line is drawn.
///  - If nn<0 a closed polygon is drawn.

void TTeXDump::DrawPS(Int_t nn, Double_t *xw, Double_t *yw)
{
   Int_t  n = TMath::Abs(nn);;
   Float_t x, y;

   if( n <= 1) {
      Error("DrawPS", "Two points are needed");
      return;
   }

   x = XtoTeX(xw[0]);
   y = YtoTeX(yw[0]);

   Int_t fillis = fFillStyle/1000;
   Int_t fillsi = fFillStyle%1000;

   if (nn>0) {
      if (fLineWidth<=0) return;
      SetColor(fLineColor);
      PrintStr("@");
      PrintStr("\\draw [c");
      TString spec = gStyle->GetLineStyleString(fLineStyle);
      TString tikzSpec;
      TString stripped = TString{spec.Strip(TString::kBoth)};
      if (stripped.Length()) {
         tikzSpec.Append(",dash pattern=");
         Ssiz_t i{0}, j{0};
         bool on{true}, iterate{true};

         while (iterate){
            j = stripped.Index(" ", 1, i, TString::kExact);
            if (j == kNPOS){
               iterate = false;
               j = stripped.Length();
            }

            if (on) {
               tikzSpec.Append("on ");
               on = false;
            } else {
               tikzSpec.Append("off ");
               on = true;
            }
            int num = TString{stripped(i, j - i)}.Atoi();
            float pt = 0.2*num;
            tikzSpec.Append(TString::Format("%.2fpt ", pt));
            i = j + 1;
         }
         PrintStr(tikzSpec.Data());
      }
      PrintStr(",line width=");
      WriteReal(0.3*fLineScale*fLineWidth, kFALSE);
      if (fCurrentAlpha != 1.) {
         PrintStr(",opacity=");
         WriteReal(fCurrentAlpha, kFALSE);
      }
   } else {
      SetColor(fFillColor);
      if (fillis==1) {
         PrintStr("@");
         PrintStr("\\draw [c, fill=c");
      } else if (fillis==0) {
         PrintStr("@");
         PrintStr("\\draw [c");
      } else {
         PrintStr("\\draw [pattern=");
         switch (fillsi) {
         case 1 :
            PrintStr("crosshatch dots");
            break;
         case 2 :
         case 3 :
            PrintStr("dots");
            break;
         case 4 :
            PrintStr("north east lines");
            break;
         case 5 :
            PrintStr("north west lines");
            break;
         case 6 :
            PrintStr("vertical lines");
            break;
         case 7 :
            PrintStr("horizontal lines");
            break;
         case 10 :
             PrintStr("bricks");
            break;
         case 13 :
            PrintStr("crosshatch");
            break;
         }
         PrintStr(", draw=none, pattern color=c");
      }
      if (fCurrentAlpha != 1.) {
         PrintStr(", fill opacity=");
         WriteReal(fCurrentAlpha, kFALSE);
      }
   }
   PrintStr("] (");
   WriteReal(x, kFALSE);
   PrintFast(1,",");
   WriteReal(y, kFALSE);
   PrintStr(") -- ");

   for (Int_t i=1;i<n;i++) {
      x = XtoTeX(xw[i]);
      y = YtoTeX(yw[i]);
      PrintFast(1,"(");
      WriteReal(x, kFALSE);
      PrintFast(1,",");
      WriteReal(y, kFALSE);
      PrintFast(1,")");
      if (i<n-1) PrintStr(" -- ");
      else PrintStr(";@");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Start the TeX page. This function starts the tikzpicture environment

void TTeXDump::NewPage()
{
   // Compute pad conversion coefficients
   if (gPad) {
      Double_t ww   = gPad->GetWw();
      Double_t wh   = gPad->GetWh();
      fYsize        = fXsize*wh/ww;
   } else {
      fYsize = 27;
   }

   if(!fBoundingBox) {
      PrintStr("\\begin{tikzpicture}@");
      DefineMarkers();
      fBoundingBox = kTRUE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the range for the paper in centimetres

void TTeXDump::Range(Float_t xsize, Float_t ysize)
{
   fXsize = xsize;
   fYsize = ysize;

   fRange = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set color index for fill areas

void TTeXDump::SetFillColor( Color_t cindex )
{
   fFillColor = cindex;
   if (gStyle->GetFillColor() <= 0) cindex = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set color index for lines

void TTeXDump::SetLineColor( Color_t cindex )
{
   fLineColor = cindex;
}

////////////////////////////////////////////////////////////////////////////////
/// Change the line style
///
///  - linestyle = 2 dashed
///  - linestyle = 3 dotted
///  - linestyle = 4 dash-dotted
///  - linestyle = else solid (1 in is used most of the time)

void TTeXDump::SetLineStyle(Style_t linestyle)
{
   fLineStyle = linestyle;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the lines width.

void TTeXDump::SetLineWidth(Width_t linewidth)
{
   fLineWidth = linewidth;
}

////////////////////////////////////////////////////////////////////////////////
/// Set size for markers.

void TTeXDump::SetMarkerSize( Size_t msize)
{
   fMarkerSize = msize;
}

////////////////////////////////////////////////////////////////////////////////
/// Set color index for markers.

void TTeXDump::SetMarkerColor( Color_t cindex)
{
   fMarkerColor = cindex;
}

////////////////////////////////////////////////////////////////////////////////
/// Set color with its color index

void TTeXDump::SetColor(Int_t color)
{
   if (color < 0) color = 0;
   TColor *col = gROOT->GetColor(color);

   if (col) {
      SetColor(col->GetRed(), col->GetGreen(), col->GetBlue());
      fCurrentAlpha = col->GetAlpha();
   } else {
      SetColor(1., 1., 1.);
      fCurrentAlpha = 1.;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set color with its R G B components
///
///  - r: % of red in [0,1]
///  - g: % of green in [0,1]
///  - b: % of blue in [0,1]

void TTeXDump::SetColor(Float_t r, Float_t g, Float_t b)
{
   if (fCurrentRed == r && fCurrentGreen == g && fCurrentBlue == b) return;

   fCurrentRed   = r;
   fCurrentGreen = g;
   fCurrentBlue  = b;
   PrintStr("@");
   PrintStr("\\definecolor{c}{rgb}{");
   WriteReal(r, kFALSE);
   PrintFast(1,",");
   WriteReal(g, kFALSE);
   PrintFast(1,",");
   WriteReal(b, kFALSE);
   PrintFast(2,"};");
}

////////////////////////////////////////////////////////////////////////////////
/// Set color index for text

void TTeXDump::SetTextColor( Color_t cindex )
{
   fTextColor = cindex;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw text
///
///  - xx: x position of the text
///  - yy: y position of the text
///  - chars: text to be drawn

void TTeXDump::Text(Double_t x, Double_t y, const char *chars)
{
   Double_t wh = (Double_t)gPad->XtoPixel(gPad->GetX2());
   Double_t hh = (Double_t)gPad->YtoPixel(gPad->GetY1());
   Float_t tsize, ftsize;
   if (wh < hh) {
      tsize = fTextSize*wh;
      Int_t sizeTTF = (Int_t)(tsize+0.5);
      ftsize = (sizeTTF*fXsize*gPad->GetAbsWNDC())/wh;
   } else {
      tsize = fTextSize*hh;
      Int_t sizeTTF = (Int_t)(tsize+0.5);
      ftsize = (sizeTTF*fYsize*gPad->GetAbsHNDC())/hh;
   }
   ftsize *= 2.22097;
   if (ftsize <= 0) return;

   TString t(chars);
   if (t.Index("\\")>=0 || t.Index("^{")>=0 || t.Index("_{")>=0) {
      t.Prepend("$");
      t.Append("$");
   } else {
      t.ReplaceAll("<","$<$");
      t.ReplaceAll(">","$>$");
      t.ReplaceAll("_","\\_");
   }
   t.ReplaceAll("&","\\&");
   t.ReplaceAll("#","\\#");
   t.ReplaceAll("%","\\%");

   Int_t txalh = fTextAlign/10;
   if (txalh <1) txalh = 1; else if (txalh > 3) txalh = 3;
   Int_t txalv = fTextAlign%10;
   if (txalv <1) txalv = 1; else if (txalv > 3) txalv = 3;
   SetColor(fTextColor);
   PrintStr("@");
   PrintStr("\\draw");
   if (txalh!=2 || txalv!=2) {
      PrintStr(" [anchor=");
      if (txalv==1) PrintStr("base");
      if (txalv==3) PrintStr("north");
      if (txalh==1) PrintStr(" west");
      if (txalh==3) PrintStr(" east");
      PrintFast(1,"]");
   }
   PrintFast(2," (");
   WriteReal(XtoTeX(x), kFALSE);
   PrintFast(1,",");
   WriteReal(YtoTeX(y), kFALSE);
   PrintStr(") node[scale=");
   WriteReal(ftsize, kFALSE);
   PrintStr(", color=c");
   if (fCurrentAlpha != 1.) {
      PrintStr(",opacity=");
      WriteReal(fCurrentAlpha, kFALSE);
   }
   PrintStr(", rotate=");
   WriteReal(fTextAngle, kFALSE);
   PrintFast(2,"]{");
   PrintStr(t.Data());
   PrintFast(2,"};");
}

////////////////////////////////////////////////////////////////////////////////
/// Write a string of characters in NDC

void TTeXDump::TextNDC(Double_t u, Double_t v, const char *chars)
{
   Double_t x = gPad->GetX1() + u*(gPad->GetX2() - gPad->GetX1());
   Double_t y = gPad->GetY1() + v*(gPad->GetY2() - gPad->GetY1());
   Text(x, y, chars);
}

////////////////////////////////////////////////////////////////////////////////
/// Convert U from NDC coordinate to TeX

Float_t TTeXDump::UtoTeX(Double_t u)
{
   Double_t cm = fXsize*(gPad->GetAbsXlowNDC() + u*gPad->GetAbsWNDC());
   return cm;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert V from NDC coordinate to TeX

Float_t TTeXDump::VtoTeX(Double_t v)
{
   Double_t cm = fYsize*(gPad->GetAbsYlowNDC() + v*gPad->GetAbsHNDC());
   return cm;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert X from world coordinate to TeX

Float_t TTeXDump::XtoTeX(Double_t x)
{
   Double_t u = (x - gPad->GetX1())/(gPad->GetX2() - gPad->GetX1());
   return  UtoTeX(u);
}

////////////////////////////////////////////////////////////////////////////////
/// Convert Y from world coordinate to TeX

Float_t TTeXDump::YtoTeX(Double_t y)
{
   Double_t v = (y - gPad->GetY1())/(gPad->GetY2() - gPad->GetY1());
   return  VtoTeX(v);
}

////////////////////////////////////////////////////////////////////////////////
/// Begin the Cell Array painting

void TTeXDump::CellArrayBegin(Int_t, Int_t, Double_t, Double_t, Double_t,
                          Double_t)
{
   Warning("CellArrayBegin", "not yet implemented");
}

////////////////////////////////////////////////////////////////////////////////
/// Paint the Cell Array

void TTeXDump::CellArrayFill(Int_t, Int_t, Int_t)
{
   Warning("CellArrayFill", "not yet implemented");
}

////////////////////////////////////////////////////////////////////////////////
/// End the Cell Array painting

void TTeXDump::CellArrayEnd()
{
   Warning("CellArrayEnd", "not yet implemented");
}

////////////////////////////////////////////////////////////////////////////////
/// Not needed in TeX case

void TTeXDump::DrawPS(Int_t, Float_t *, Float_t *)
{
   Warning("DrawPS", "not yet implemented");
}

////////////////////////////////////////////////////////////////////////////////
/// add additional pgfplotmarks

void TTeXDump::DefineMarkers()
{
  // open cross
  PrintStr("\\pgfdeclareplotmark{cross} {@");
  PrintStr("\\pgfpathmoveto{\\pgfpoint{-0.3\\pgfplotmarksize}{\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{+0.3\\pgfplotmarksize}{\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{+0.3\\pgfplotmarksize}{0.3\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{+1\\pgfplotmarksize}{0.3\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{+1\\pgfplotmarksize}{-0.3\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{+0.3\\pgfplotmarksize}{-0.3\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{+0.3\\pgfplotmarksize}{-1.\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{-0.3\\pgfplotmarksize}{-1.\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{-0.3\\pgfplotmarksize}{-0.3\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{-1.\\pgfplotmarksize}{-0.3\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{-1.\\pgfplotmarksize}{0.3\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{-0.3\\pgfplotmarksize}{0.3\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathclose@");
  PrintStr("\\pgfusepathqstroke@");
  PrintStr("}@");

  // filled cross
  PrintStr("\\pgfdeclareplotmark{cross*} {@");
  PrintStr("\\pgfpathmoveto{\\pgfpoint{-0.3\\pgfplotmarksize}{\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{+0.3\\pgfplotmarksize}{\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{+0.3\\pgfplotmarksize}{0.3\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{+1\\pgfplotmarksize}{0.3\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{+1\\pgfplotmarksize}{-0.3\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{+0.3\\pgfplotmarksize}{-0.3\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{+0.3\\pgfplotmarksize}{-1.\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{-0.3\\pgfplotmarksize}{-1.\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{-0.3\\pgfplotmarksize}{-0.3\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{-1.\\pgfplotmarksize}{-0.3\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{-1.\\pgfplotmarksize}{0.3\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfpoint{-0.3\\pgfplotmarksize}{0.3\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathclose@");
  PrintStr("\\pgfusepathqfillstroke@");
  PrintStr("}@");

  // open star
  PrintStr("\\pgfdeclareplotmark{newstar} {@");
  PrintStr("\\pgfpathmoveto{\\pgfqpoint{0pt}{\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfqpointpolar{44}{0.5\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfqpointpolar{18}{\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfqpointpolar{-20}{0.5\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfqpointpolar{-54}{\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfqpointpolar{-90}{0.5\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfqpointpolar{234}{\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfqpointpolar{198}{0.5\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfqpointpolar{162}{\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfqpointpolar{134}{0.5\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathclose@");
  PrintStr("\\pgfusepathqstroke@");
  PrintStr("}@");

  // filled star
  PrintStr("\\pgfdeclareplotmark{newstar*} {@");
  PrintStr("\\pgfpathmoveto{\\pgfqpoint{0pt}{\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfqpointpolar{44}{0.5\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfqpointpolar{18}{\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfqpointpolar{-20}{0.5\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfqpointpolar{-54}{\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfqpointpolar{-90}{0.5\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfqpointpolar{234}{\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfqpointpolar{198}{0.5\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfqpointpolar{162}{\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathlineto{\\pgfqpointpolar{134}{0.5\\pgfplotmarksize}}@");
  PrintStr("\\pgfpathclose@");
  PrintStr("\\pgfusepathqfillstroke@");
  PrintStr("}@");
}
