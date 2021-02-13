// @(#)root/postscript:$Id: TPDF.cxx,v 1.0
// Author: Olivier Couet

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TPDF
\ingroup PS

\brief Interface to PDF.

Like PostScript, PDF is a vector graphics output format allowing a very high
graphics output quality. The functionalities provided by this class are very
similar to those provided by `TPostScript`.

Compare to PostScript output, the PDF files are usually smaller because some
parts of them can be compressed.

PDF also allows to define table of contents. This facility can be used in ROOT.
The following example shows how to proceed:
~~~ {.cpp}
{
   TCanvas* canvas = new TCanvas("canvas");
   TH1F* histo = new TH1F("histo","test 1",10,0.,10.);
   histo->SetFillColor(2);
   histo->Fill(2.);
   histo->Draw();
   canvas->Print("plots.pdf(","Title:One bin filled");
   histo->Fill(4.);
   histo->Draw();
   canvas->Print("plots.pdf","Title:Two bins filled");
   histo->Fill(6.);
   histo->Draw();
   canvas->Print("plots.pdf","Title:Three bins filled");
   histo->Fill(8.);
   histo->Draw();
   canvas->Print("plots.pdf","Title:Four bins filled");
   histo->Fill(8.);
   histo->Draw();
   canvas->Print("plots.pdf)","Title:The fourth bin content is 2");
}
~~~
Each character string following the keyword "Title:" makes a new entry in
the table of contents.
*/

#ifdef WIN32
#pragma optimize("",off)
#endif

#include <cstdlib>
#include <cstring>
#include <cctype>
#include <fstream>

#include "TROOT.h"
#include "TDatime.h"
#include "TColor.h"
#include "TVirtualPad.h"
#include "TPoints.h"
#include "TPDF.h"
#include "TStyle.h"
#include "TMath.h"
#include "TStorage.h"
#include "TText.h"
#include "zlib.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "snprintf.h"

// To scale fonts to the same size as the old TT version
const Float_t kScale = 0.93376068;

// Objects numbers
const Int_t kObjRoot             =  1; // Root object
const Int_t kObjInfo             =  2; // Info object
const Int_t kObjOutlines         =  3; // Outlines object
const Int_t kObjPages            =  4; // Pages object (pages index)
const Int_t kObjPageResources    =  5; // Pages Resources object
const Int_t kObjContents         =  6; // Table of content
const Int_t kObjFont             =  7; // First Font object (14 in total)
const Int_t kObjColorSpace       = 22; // ColorSpace object
const Int_t kObjPatternResourses = 23; // Pattern Resources object
const Int_t kObjPatternList      = 24; // Pattern list object
const Int_t kObjTransList        = 25; // List of transparencies
const Int_t kObjPattern          = 26; // First pattern object (25 in total)
const Int_t kObjFirstPage        = 51; // First page object

// Number of fonts
const Int_t kNumberOfFonts = 15;

Int_t TPDF::fgLineJoin = 0;
Int_t TPDF::fgLineCap  = 0;

ClassImp(TPDF);

////////////////////////////////////////////////////////////////////////////////
/// Default PDF constructor

TPDF::TPDF() : TVirtualPS()
{
   fStream          = 0;
   fCompress        = kFALSE;
   fPageNotEmpty    = kFALSE;
   gVirtualPS       = this;
   fRed             = 0.;
   fGreen           = 0.;
   fBlue            = 0.;
   fAlpha           = 1.;
   fXsize           = 0.;
   fYsize           = 0.;
   fType            = 0;
   fPageFormat      = 0;
   fPageOrientation = 0;
   fStartStream     = 0;
   fLineScale       = 0.;
   fObjPosSize      = 0;
   fObjPos          = 0;
   fNbObj           = 0;
   fNbPage          = 0;
   fRange           = kFALSE;
   SetTitle("PDF");
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the PDF interface
///
///  - fname : PDF file name
///  - wtype : PDF workstation type. Not used in the PDF driver. But as TPDF
///            inherits from TVirtualPS it should be kept. Anyway it is not
///            necessary to specify this parameter at creation time because it
///            has a default value (which is ignore in the PDF case).

TPDF::TPDF(const char *fname, Int_t wtype) : TVirtualPS(fname, wtype)
{
   fStream          = 0;
   fCompress        = kFALSE;
   fPageNotEmpty    = kFALSE;
   fRed             = 0.;
   fGreen           = 0.;
   fBlue            = 0.;
   fAlpha           = 1.;
   fXsize           = 0.;
   fYsize           = 0.;
   fType            = 0;
   fPageFormat      = 0;
   fPageOrientation = 0;
   fStartStream     = 0;
   fLineScale       = 0.;
   fObjPosSize      = 0;
   fNbObj           = 0;
   fNbPage          = 0;
   fRange           = kFALSE;
   SetTitle("PDF");
   Open(fname, wtype);
}

////////////////////////////////////////////////////////////////////////////////
/// Default PDF destructor

TPDF::~TPDF()
{
   Close();

   if (fObjPos) delete [] fObjPos;
}

////////////////////////////////////////////////////////////////////////////////
/// Begin the Cell Array painting

void TPDF::CellArrayBegin(Int_t, Int_t, Double_t, Double_t, Double_t,
                          Double_t)
{
   Warning("TPDF::CellArrayBegin", "not yet implemented");
}

////////////////////////////////////////////////////////////////////////////////
/// Paint the Cell Array

void TPDF::CellArrayFill(Int_t, Int_t, Int_t)
{
   Warning("TPDF::CellArrayFill", "not yet implemented");
}

////////////////////////////////////////////////////////////////////////////////
/// End the Cell Array painting

void TPDF::CellArrayEnd()
{
   Warning("TPDF::CellArrayEnd", "not yet implemented");
}

////////////////////////////////////////////////////////////////////////////////
/// Close a PDF file

void TPDF::Close(Option_t *)
{
   Int_t i;

   if (!gVirtualPS) return;
   if (!fStream) return;
   if (gPad) gPad->Update();

   // Close the currently opened page
   WriteCompressedBuffer();
   PrintStr("endstream@");
   Int_t streamLength = fNByte-fStartStream-10;
   PrintStr("endobj@");
   NewObject(4*(fNbPage-1)+kObjFirstPage+2);
   WriteInteger(streamLength, 0);
   PrintStr("@");
   PrintStr("endobj@");
   NewObject(4*(fNbPage-1)+kObjFirstPage+3);
   PrintStr("<<@");
   if (!strstr(GetTitle(),"PDF")) {
      PrintStr("/Title (");
      PrintStr(GetTitle());
      PrintStr(")@");
   } else {
      PrintStr("/Title (Page");
      WriteInteger(fNbPage);
      PrintStr(")@");
   }
   PrintStr("/Dest [");
   WriteInteger(4*(fNbPage-1)+kObjFirstPage);
   PrintStr(" 0 R /XYZ null null 0]@");
   PrintStr("/Parent");
   WriteInteger(kObjContents);
   PrintStr(" 0 R");
   PrintStr("@");
   if (fNbPage > 1) {
      PrintStr("/Prev");
      WriteInteger(4*(fNbPage-2)+kObjFirstPage+3);
      PrintStr(" 0 R");
      PrintStr("@");
   }
   PrintStr(">>@");

   NewObject(kObjOutlines);
   PrintStr("<<@");
   PrintStr("/Type /Outlines@");
   PrintStr("/Count");
   WriteInteger(fNbPage+1);
   PrintStr("@");
   PrintStr("/First");
   WriteInteger(kObjContents);
   PrintStr(" 0 R");
   PrintStr("@");
   PrintStr("/Last");
   WriteInteger(kObjContents);
   PrintStr(" 0 R");
   PrintStr("@");
   PrintStr(">>@");
   PrintStr("endobj@");

   NewObject(kObjContents);
   PrintStr("<<@");
   PrintStr("/Title (Contents)@");
   PrintStr("/Dest [");
   WriteInteger(kObjFirstPage);
   PrintStr(" 0 R /XYZ null null 0]@");
   PrintStr("/Count");
   WriteInteger(fNbPage);
   PrintStr("@");
   PrintStr("/Parent");
   WriteInteger(kObjOutlines);
   PrintStr(" 0 R");
   PrintStr("@");
   PrintStr("/First");
   WriteInteger(kObjFirstPage+3);
   PrintStr(" 0 R");
   PrintStr("@");
   PrintStr("/Last");
   WriteInteger(4*(fNbPage-1)+kObjFirstPage+3);
   PrintStr(" 0 R");
   PrintStr("@");
   PrintStr(">>@");

   // List of all the pages
   NewObject(kObjPages);
   PrintStr("<<@");
   PrintStr("/Type /Pages@");
   PrintStr("/Count");
   WriteInteger(fNbPage);
   PrintStr("@");
   PrintStr("/Kids [");
   for (i=1; i<=fNbPage; i++) {
      WriteInteger(4*(i-1)+kObjFirstPage);
      PrintStr(" 0 R");
   }
   PrintStr(" ]");
   PrintStr("@");
   PrintStr(">>@");
   PrintStr("endobj@");

   NewObject(kObjTransList);
   PrintStr("<<@");
   for (i=0; i<(int)fAlphas.size(); i++) {
      PrintStr(
      Form("/ca%3.2f << /Type /ExtGState /ca %3.2f >> /CA%3.2f << /Type /ExtGState /CA %3.2f >>@",
      fAlphas[i],fAlphas[i],fAlphas[i],fAlphas[i]));
   }
   PrintStr(">>@");
   PrintStr("endobj@");
   if (fAlphas.size()) fAlphas.clear();

   // Cross-Reference Table
   Int_t refInd = fNByte;
   PrintStr("xref@");
   PrintStr("0");
   WriteInteger(fNbObj+1);
   PrintStr("@");
   PrintStr("0000000000 65535 f @");
   char str[21];
   for (i=0; i<fNbObj; i++) {
      snprintf(str,21,"%10.10d 00000 n @",fObjPos[i]);
      PrintStr(str);
   }

   // Trailer
   PrintStr("trailer@");
   PrintStr("<<@");
   PrintStr("/Size");
   WriteInteger(fNbObj+1);
   PrintStr("@");
   PrintStr("/Root");
   WriteInteger(kObjRoot);
   PrintStr(" 0 R");
   PrintStr("@");
   PrintStr("/Info");
   WriteInteger(kObjInfo);
   PrintStr(" 0 R@");
   PrintStr(">>@");
   PrintStr("startxref@");
   WriteInteger(refInd, 0);
   PrintStr("@");
   PrintStr("%%EOF@");

   // Close file stream
   if (fStream) { fStream->close(); delete fStream; fStream = 0;}

   gVirtualPS = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a Box

void TPDF::DrawBox(Double_t x1, Double_t y1, Double_t x2, Double_t  y2)
{
   static Double_t x[4], y[4];
   Double_t ix1 = XtoPDF(x1);
   Double_t ix2 = XtoPDF(x2);
   Double_t iy1 = YtoPDF(y1);
   Double_t iy2 = YtoPDF(y2);
   Int_t fillis = fFillStyle/1000;
   Int_t fillsi = fFillStyle%1000;

   if (fillis == 3 || fillis == 2) {
      if (fillsi > 99) {
         x[0] = x1;   y[0] = y1;
         x[1] = x2;   y[1] = y1;
         x[2] = x2;   y[2] = y2;
         x[3] = x1;   y[3] = y2;
         return;
      }
      if (fillsi > 0 && fillsi < 26) {
         x[0] = x1;   y[0] = y1;
         x[1] = x2;   y[1] = y1;
         x[2] = x2;   y[2] = y2;
         x[3] = x1;   y[3] = y2;
         DrawPS(-4, &x[0], &y[0]);
      }
      if (fillsi == -3) {
         SetColor(5);
         if (fAlpha == 1) PrintFast(15," q 0.4 w [] 0 d");
         WriteReal(ix1);
         WriteReal(iy1);
         WriteReal(ix2 - ix1);
         WriteReal(iy2 - iy1);
         if (fAlpha == 1) PrintFast(8," re b* Q");
         else             PrintFast(6," re f*");
      }
   }
   if (fillis == 1) {
      SetColor(fFillColor);
      if (fAlpha == 1) PrintFast(15," q 0.4 w [] 0 d");
      WriteReal(ix1);
      WriteReal(iy1);
      WriteReal(ix2 - ix1);
      WriteReal(iy2 - iy1);
      if (fAlpha == 1) PrintFast(8," re b* Q");
      else             PrintFast(6," re f*");
   }
   if (fillis == 0) {
      if (fLineWidth<=0) return;
      SetColor(fLineColor);
      WriteReal(ix1);
      WriteReal(iy1);
      WriteReal(ix2 - ix1);
      WriteReal(iy2 - iy1);
      PrintFast(5," re S");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a Frame around a box
///
///  - mode = -1  box looks as it is behind the screen
///  - mode =  1  box looks as it is in front of the screen
///  - border is the border size in already precomputed PDF units
///  - dark  is the color for the dark part of the frame
///  - light is the color for the light part of the frame

void TPDF::DrawFrame(Double_t xl, Double_t yl, Double_t xt, Double_t  yt,
                            Int_t mode, Int_t border, Int_t dark, Int_t light)
{
   static Double_t xps[7], yps[7];
   Int_t i;

   // Draw top&left part of the box
   if (mode == -1) SetColor(dark);
   else            SetColor(light);
   xps[0] = XtoPDF(xl);          yps[0] = YtoPDF(yl);
   xps[1] = xps[0] + border;     yps[1] = yps[0] + border;
   xps[2] = xps[1];              yps[2] = YtoPDF(yt) - border;
   xps[3] = XtoPDF(xt) - border; yps[3] = yps[2];
   xps[4] = XtoPDF(xt);          yps[4] = YtoPDF(yt);
   xps[5] = xps[0];              yps[5] = yps[4];
   xps[6] = xps[0];              yps[6] = yps[0];

   MoveTo(xps[0], yps[0]);
   for (i=1;i<7;i++) LineTo(xps[i], yps[i]);
   PrintFast(3," f*");

   // Draw bottom&right part of the box
   if (mode == -1) SetColor(light);
   else            SetColor(dark);
   xps[0] = XtoPDF(xl);          yps[0] = YtoPDF(yl);
   xps[1] = xps[0] + border;     yps[1] = yps[0] + border;
   xps[2] = XtoPDF(xt) - border; yps[2] = yps[1];
   xps[3] = xps[2];              yps[3] = YtoPDF(yt) - border;
   xps[4] = XtoPDF(xt);          yps[4] = YtoPDF(yt);
   xps[5] = xps[4];              yps[5] = yps[0];
   xps[6] = xps[0];              yps[6] = yps[0];

   MoveTo(xps[0], yps[0]);
   for (i=1;i<7;i++) LineTo(xps[i], yps[i]);
   PrintFast(3," f*");
}

////////////////////////////////////////////////////////////////////////////////
/// Draw Fill area with hatch styles

void TPDF::DrawHatch(Float_t, Float_t, Int_t, Float_t *, Float_t *)
{
   Warning("DrawHatch", "hatch fill style not yet implemented");
}

////////////////////////////////////////////////////////////////////////////////
/// Draw Fill area with hatch styles

void TPDF::DrawHatch(Float_t, Float_t, Int_t, Double_t *, Double_t *)
{
   Warning("DrawHatch", "hatch fill style not yet implemented");
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a PolyLine
///
///  Draw a polyline through  the points xy.
///
///  - If NN=1 moves only to point x,y.
///  - If NN=0 the x,y are  written  in the PDF file
///       according to the current transformation.
///  - If NN>0 the line is clipped as a line.
///  - If NN<0 the line is clipped as a fill area.

void TPDF::DrawPolyLine(Int_t nn, TPoints *xy)
{
   Int_t  n;

   Style_t linestylesav = fLineStyle;
   Width_t linewidthsav = fLineWidth;

   if (nn > 0) {
      if (fLineWidth<=0) return;
      n = nn;
      SetLineStyle(fLineStyle);
      SetLineWidth(fLineWidth);
      SetColor(Int_t(fLineColor));
   } else {
      n = -nn;
      SetLineStyle(1);
      SetLineWidth(1);
      SetColor(Int_t(fLineColor));
   }

   WriteReal(XtoPDF(xy[0].GetX()));
   WriteReal(YtoPDF(xy[0].GetY()));
   if (n <= 1) {
      if (n == 0) return;
      PrintFast(2," m");
      return;
   }

   PrintFast(2," m");

   for (Int_t i=1;i<n;i++) LineTo(XtoPDF(xy[i].GetX()), YtoPDF(xy[i].GetY()));

   if (nn > 0) {
      if (xy[0].GetX() == xy[n-1].GetX() && xy[0].GetY() == xy[n-1].GetY()) PrintFast(3," cl");
      PrintFast(2," S");
   } else {
      PrintFast(3," f*");
   }

   SetLineStyle(linestylesav);
   SetLineWidth(linewidthsav);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a PolyLine in NDC space
///
///  Draw a polyline through the points xy.
///
///  - If NN=1 moves only to point x,y.
///  - If NN=0 the x,y are  written in the PDF file
///       according to the current transformation.
///  - If NN>0 the line is clipped as a line.
///  - If NN<0 the line is clipped as a fill area.

void TPDF::DrawPolyLineNDC(Int_t nn, TPoints *xy)
{
   Int_t  n;

   Style_t linestylesav = fLineStyle;
   Width_t linewidthsav = fLineWidth;

   if (nn > 0) {
      if (fLineWidth<=0) return;
      n = nn;
      SetLineStyle(fLineStyle);
      SetLineWidth(fLineWidth);
      SetColor(Int_t(fLineColor));
   } else {
      n = -nn;
      SetLineStyle(1);
      SetLineWidth(1);
      SetColor(Int_t(fLineColor));
   }

   WriteReal(UtoPDF(xy[0].GetX()));
   WriteReal(VtoPDF(xy[0].GetY()));
   if (n <= 1) {
      if (n == 0) return;
      PrintFast(2," m");
      return;
   }

   PrintFast(2," m");

   for (Int_t i=1;i<n;i++) LineTo(UtoPDF(xy[i].GetX()), VtoPDF(xy[i].GetY()));

   if (nn > 0) {
      if (xy[0].GetX() == xy[n-1].GetX() && xy[0].GetY() == xy[n-1].GetY()) PrintFast(3," cl");
      PrintFast(2," S");
   } else {
      PrintFast(3," f*");
   }

   SetLineStyle(linestylesav);
   SetLineWidth(linewidthsav);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw markers at the n WC points xw, yw

void TPDF::DrawPolyMarker(Int_t n, Float_t *xw, Float_t *yw)
{
   fMarkerStyle = TMath::Abs(fMarkerStyle);
   Style_t linestylesav = fLineStyle;
   Width_t linewidthsav = fLineWidth;
   SetLineStyle(1);
   SetLineWidth(TMath::Max(1, Int_t(TAttMarker::GetMarkerLineWidth(fMarkerStyle))));
   SetColor(Int_t(fMarkerColor));
   Int_t ms = TAttMarker::GetMarkerStyleBase(fMarkerStyle);

   if (ms == 4)
      ms = 24;
   else if (ms >= 6 && ms <= 8)
      ms = 20;
   else if (ms >= 9 && ms <= 19)
      ms = 1;

   // Define the marker size
   Float_t msize  = fMarkerSize - TMath::Floor(TAttMarker::GetMarkerLineWidth(fMarkerStyle)/2.)/4.*fLineScale;
   if (fMarkerStyle == 1 || (fMarkerStyle >= 9 && fMarkerStyle <= 19)) {
     msize = 1.;
   } else if (fMarkerStyle == 6) {
     msize = 1.;
   } else if (fMarkerStyle == 7) {
     msize = 1.5;
   } else {
      const Int_t kBASEMARKER = 8;
      Float_t sbase = msize*kBASEMARKER;
      Float_t s2x = sbase / Float_t(gPad->GetWw() * gPad->GetAbsWNDC());
      msize = this->UtoPDF(s2x) - this->UtoPDF(0);
   }

   Double_t m  = msize;
   Double_t m2 = m/2;
   Double_t m3 = m/3;
   Double_t m4 = m2*1.333333333333;
   Double_t m6 = m/6;
   Double_t m0 = m/10.;
   Double_t m8 = m/4;
   Double_t m9 = m/8;

   // Draw the marker according to the type
   Double_t ix,iy;
   for (Int_t i=0;i<n;i++) {
      ix = XtoPDF(xw[i]);
      iy = YtoPDF(yw[i]);
      // Dot (.)
      if (ms == 1) {
         MoveTo(ix-1, iy);
         LineTo(ix  , iy);
      // Plus (+)
      } else if (ms == 2) {
         MoveTo(ix-m2, iy);
         LineTo(ix+m2, iy);
         MoveTo(ix   , iy-m2);
         LineTo(ix   , iy+m2);
      // X shape (X)
      } else if (ms == 5) {
         MoveTo(ix-m2*0.707, iy-m2*0.707);
         LineTo(ix+m2*0.707, iy+m2*0.707);
         MoveTo(ix-m2*0.707, iy+m2*0.707);
         LineTo(ix+m2*0.707, iy-m2*0.707);
      // Asterisk shape (*)
      } else if (ms == 3 || ms == 31) {
         MoveTo(ix-m2, iy);
         LineTo(ix+m2, iy);
         MoveTo(ix   , iy-m2);
         LineTo(ix   , iy+m2);
         MoveTo(ix-m2*0.707, iy-m2*0.707);
         LineTo(ix+m2*0.707, iy+m2*0.707);
         MoveTo(ix-m2*0.707, iy+m2*0.707);
         LineTo(ix+m2*0.707, iy-m2*0.707);
      // Circle
      } else if (ms == 24 || ms == 20) {
         MoveTo(ix-m2, iy);
         WriteReal(ix-m2); WriteReal(iy+m4);
         WriteReal(ix+m2); WriteReal(iy+m4);
         WriteReal(ix+m2); WriteReal(iy)   ; PrintFast(2," c");
         WriteReal(ix+m2); WriteReal(iy-m4);
         WriteReal(ix-m2); WriteReal(iy-m4);
         WriteReal(ix-m2); WriteReal(iy)   ; PrintFast(4," c h");
      // Square
      } else if (ms == 25 || ms == 21) {
         WriteReal(ix-m2); WriteReal(iy-m2);
         WriteReal(m)    ; WriteReal(m)    ; PrintFast(3," re");
      // Down triangle
      } else if (ms == 23 || ms == 32) {
         MoveTo(ix   , iy-m2);
         LineTo(ix+m2, iy+m2);
         LineTo(ix-m2, iy+m2);
         PrintFast(2," h");
      // Up triangle
      } else if (ms == 26 || ms == 22) {
         MoveTo(ix-m2, iy-m2);
         LineTo(ix+m2, iy-m2);
         LineTo(ix   , iy+m2);
         PrintFast(2," h");
      } else if (ms == 27 || ms == 33) {
         MoveTo(ix   , iy-m2);
         LineTo(ix+m3, iy);
         LineTo(ix   , iy+m2);
         LineTo(ix-m3, iy)   ;
         PrintFast(2," h");
      } else if (ms == 28 || ms == 34) {
         MoveTo(ix-m6, iy-m6);
         LineTo(ix-m6, iy-m2);
         LineTo(ix+m6, iy-m2);
         LineTo(ix+m6, iy-m6);
         LineTo(ix+m2, iy-m6);
         LineTo(ix+m2, iy+m6);
         LineTo(ix+m6, iy+m6);
         LineTo(ix+m6, iy+m2);
         LineTo(ix-m6, iy+m2);
         LineTo(ix-m6, iy+m6);
         LineTo(ix-m2, iy+m6);
         LineTo(ix-m2, iy-m6);
         PrintFast(2," h");
      } else if (ms == 29 || ms == 30) {
         MoveTo(ix           , iy+m2);
         LineTo(ix+0.112255*m, iy+0.15451*m);
         LineTo(ix+0.47552*m , iy+0.15451*m);
         LineTo(ix+0.181635*m, iy-0.05902*m);
         LineTo(ix+0.29389*m , iy-0.40451*m);
         LineTo(ix           , iy-0.19098*m);
         LineTo(ix-0.29389*m , iy-0.40451*m);
         LineTo(ix-0.181635*m, iy-0.05902*m);
         LineTo(ix-0.47552*m , iy+0.15451*m);
         LineTo(ix-0.112255*m, iy+0.15451*m);
         PrintFast(2," h");
      } else if (ms == 35 ) {
      // diamond with cross
         MoveTo(ix-m2, iy   );
         LineTo(ix   , iy-m2);
         LineTo(ix+m2, iy   );
         LineTo(ix   , iy+m2);
         LineTo(ix-m2, iy   );
         LineTo(ix+m2, iy   );
         LineTo(ix   , iy+m2);
         LineTo(ix   , iy-m2);
         PrintFast(2," h");
      } else if (ms == 36 ) {
      // square with diagonal cross
         MoveTo(ix-m2, iy-m2);
         LineTo(ix+m2, iy-m2);
         LineTo(ix+m2, iy+m2);
         LineTo(ix-m2, iy+m2);
         LineTo(ix-m2, iy-m2);
         LineTo(ix+m2, iy+m2);
         LineTo(ix-m2, iy+m2);
         LineTo(ix+m2, iy-m2);
         PrintFast(2," h");
      } else if (ms == 37 || ms == 39 ) {
      // square with cross
         MoveTo(ix   , iy   );
         LineTo(ix-m8, iy+m2);
         LineTo(ix-m2, iy   );
         LineTo(ix   , iy   );
         LineTo(ix-m8, iy-m2);
         LineTo(ix+m8, iy-m2);
         LineTo(ix   , iy   );
         LineTo(ix+m2, iy   );
         LineTo(ix+m8, iy+m2);
         LineTo(ix   , iy   );
         PrintFast(2," h");
      } else if (ms == 38 ) {
      // + shaped marker with octagon
         MoveTo(ix-m2, iy   );
         LineTo(ix-m2, iy-m8);
         LineTo(ix-m8, iy-m2);
         LineTo(ix+m8, iy-m2);
         LineTo(ix+m2, iy-m8);
         LineTo(ix+m2, iy+m8);
         LineTo(ix+m8, iy+m2);
         LineTo(ix-m8, iy+m2);
         LineTo(ix-m2, iy+m8);
         LineTo(ix-m2, iy   );
         LineTo(ix+m2, iy   );
         LineTo(ix   , iy   );
         LineTo(ix   , iy-m2);
         LineTo(ix   , iy+m2);
         LineTo(ix   , iy);
         PrintFast(2," h");
      } else if (ms == 40 || ms == 41 ) {
       // four triangles X
         MoveTo(ix   , iy   );
         LineTo(ix+m8, iy+m2);
         LineTo(ix+m2, iy+m8);
         LineTo(ix   , iy   );
         LineTo(ix+m2, iy-m8);
         LineTo(ix+m8, iy-m2);
         LineTo(ix   , iy   );
         LineTo(ix-m8, iy-m2);
         LineTo(ix-m2, iy-m8);
         LineTo(ix   , iy   );
         LineTo(ix-m2, iy+m8);
         LineTo(ix-m8, iy+m2);
         LineTo(ix   , iy   );
         PrintFast(2," h");
      } else if (ms == 42 || ms == 43 ) {
      // double diamonds
         MoveTo(ix   , iy+m2);
         LineTo(ix-m9, iy+m9);
         LineTo(ix-m2, iy   );
         LineTo(ix-m9, iy-m9);
         LineTo(ix   , iy-m2);
         LineTo(ix+m9, iy-m9);
         LineTo(ix+m2, iy   );
         LineTo(ix+m9, iy+m9);
         LineTo(ix   , iy+m2);
         PrintFast(2," h");
      } else if (ms == 44 ) {
      // open four triangles plus
         MoveTo(ix   , iy   );
         LineTo(ix+m8, iy+m2);
         LineTo(ix-m8, iy+m2);
         LineTo(ix+m8, iy-m2);
         LineTo(ix-m8, iy-m2);
         LineTo(ix   , iy   );
         LineTo(ix+m2, iy+m8);
         LineTo(ix+m2, iy-m8);
         LineTo(ix-m2, iy+m8);
         LineTo(ix-m2, iy-m8);
         LineTo(ix   , iy   );
         PrintFast(2," h");
      } else if (ms == 45 ) {
      // filled four triangles plus
         MoveTo(ix+m0, iy+m0);
         LineTo(ix+m8, iy+m2);
         LineTo(ix-m8, iy+m2);
         LineTo(ix-m0, iy+m0);
         LineTo(ix-m2, iy+m8);
         LineTo(ix-m2, iy-m8);
         LineTo(ix-m0, iy-m0);
         LineTo(ix-m8, iy-m2);
         LineTo(ix+m8, iy-m2);
         LineTo(ix+m0, iy-m0);
         LineTo(ix+m2, iy-m8);
         LineTo(ix+m2, iy+m8);
         LineTo(ix+m0, iy+m0);
         PrintFast(2," h");
      } else if (ms == 46 || ms == 47 ) {
      // four triangles X
         MoveTo(ix   , iy+m8);
         LineTo(ix-m8, iy+m2);
         LineTo(ix-m2, iy+m8);
         LineTo(ix-m8, iy   );
         LineTo(ix-m2, iy-m8);
         LineTo(ix-m8, iy-m2);
         LineTo(ix   , iy-m8);
         LineTo(ix+m8, iy-m2);
         LineTo(ix+m2, iy-m8);
         LineTo(ix+m8, iy   );
         LineTo(ix+m2, iy+m8);
         LineTo(ix+m8, iy+m2);
         LineTo(ix   , iy+m8);
         PrintFast(2," h");
      } else if (ms == 48 ) {
      // four filled squares X
         MoveTo(ix   , iy+m8*1.01);
         LineTo(ix-m8, iy+m2);
         LineTo(ix-m2, iy+m8);
         LineTo(ix-m8, iy   );
         LineTo(ix-m2, iy-m8);
         LineTo(ix-m8, iy-m2);
         LineTo(ix   , iy-m8);
         LineTo(ix+m8, iy-m2);
         LineTo(ix+m2, iy-m8);
         LineTo(ix+m8, iy   );
         LineTo(ix+m2, iy+m8);
         LineTo(ix+m8, iy+m2);
         LineTo(ix   , iy+m8*0.99);
         LineTo(ix+m8*0.99, iy   );
         LineTo(ix   , iy-m8*0.99);
         LineTo(ix-m8*0.99, iy   );
         LineTo(ix   , iy+m8*0.99);
         PrintFast(2," h");
      } else if (ms == 49 ) {
      // four filled squares plus
         MoveTo(ix-m6, iy-m6*1.01);
         LineTo(ix-m6, iy-m2);
         LineTo(ix+m6, iy-m2);
         LineTo(ix+m6, iy-m6);
         LineTo(ix+m2, iy-m6);
         LineTo(ix+m2, iy+m6);
         LineTo(ix+m6, iy+m6);
         LineTo(ix+m6, iy+m2);
         LineTo(ix-m6, iy+m2);
         LineTo(ix-m6, iy+m6);
         LineTo(ix-m2, iy+m6);
         LineTo(ix-m2, iy-m6);
         LineTo(ix-m6, iy-m6*0.99);
         LineTo(ix-m6, iy+m6);
         LineTo(ix+m6, iy+m6);
         LineTo(ix+m6, iy-m6);
         PrintFast(2," h");
      } else {
         MoveTo(ix-m6, iy-m6);
         LineTo(ix-m6, iy-m2);
      }

      if ((ms > 19 && ms < 24) || ms == 29 || ms == 33 || ms == 34 ||
          ms == 39 || ms == 41 || ms == 43 || ms == 45 ||
          ms == 47 || ms == 48 || ms == 49) {
         PrintFast(2," f");
      } else {
         PrintFast(2," S");
      }
   }

   SetLineStyle(linestylesav);
   SetLineWidth(linewidthsav);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw markers at the n WC points xw, yw

void TPDF::DrawPolyMarker(Int_t n, Double_t *xw, Double_t *yw)
{
   fMarkerStyle = TMath::Abs(fMarkerStyle);
   Style_t linestylesav = fLineStyle;
   Width_t linewidthsav = fLineWidth;
   SetLineStyle(1);
   SetLineWidth(TMath::Max(1, Int_t(TAttMarker::GetMarkerLineWidth(fMarkerStyle))));
   SetColor(Int_t(fMarkerColor));
   Int_t ms = TAttMarker::GetMarkerStyleBase(fMarkerStyle);

   if (ms == 4)
      ms = 24;
   else if (ms >= 6 && ms <= 8)
      ms = 20;
   else if (ms >= 9 && ms <= 19)
      ms = 1;

   // Define the marker size
   Float_t msize  = fMarkerSize - TMath::Floor(TAttMarker::GetMarkerLineWidth(fMarkerStyle)/2.)/4.*fLineScale;
   if (fMarkerStyle == 1 || (fMarkerStyle >= 9 && fMarkerStyle <= 19)) {
     msize = 1.;
   } else if (fMarkerStyle == 6) {
     msize = 1.5;
   } else if (fMarkerStyle == 7) {
     msize = 3.;
   } else {
      const Int_t kBASEMARKER = 8;
      Float_t sbase = msize*kBASEMARKER;
      Float_t s2x = sbase / Float_t(gPad->GetWw() * gPad->GetAbsWNDC());
      msize = this->UtoPDF(s2x) - this->UtoPDF(0);
   }

   Double_t m  = msize;
   Double_t m2 = m/2;
   Double_t m3 = m/3;
   Double_t m4 = m2*1.333333333333;
   Double_t m6 = m/6;
   Double_t m8 = m/4;
   Double_t m9 = m/8;

   // Draw the marker according to the type
   Double_t ix,iy;
   for (Int_t i=0;i<n;i++) {
      ix = XtoPDF(xw[i]);
      iy = YtoPDF(yw[i]);
      // Dot (.)
      if (ms == 1) {
         MoveTo(ix-1, iy);
         LineTo(ix  , iy);
      // Plus (+)
      } else if (ms == 2) {
         MoveTo(ix-m2, iy);
         LineTo(ix+m2, iy);
         MoveTo(ix   , iy-m2);
         LineTo(ix   , iy+m2);
      // X shape (X)
      } else if (ms == 5) {
         MoveTo(ix-m2*0.707, iy-m2*0.707);
         LineTo(ix+m2*0.707, iy+m2*0.707);
         MoveTo(ix-m2*0.707, iy+m2*0.707);
         LineTo(ix+m2*0.707, iy-m2*0.707);
      // Asterisk shape (*)
      } else if (ms == 3 || ms == 31) {
         MoveTo(ix-m2, iy);
         LineTo(ix+m2, iy);
         MoveTo(ix   , iy-m2);
         LineTo(ix   , iy+m2);
         MoveTo(ix-m2*0.707, iy-m2*0.707);
         LineTo(ix+m2*0.707, iy+m2*0.707);
         MoveTo(ix-m2*0.707, iy+m2*0.707);
         LineTo(ix+m2*0.707, iy-m2*0.707);
      // Circle
      } else if (ms == 24 || ms == 20) {
         MoveTo(ix-m2, iy);
         WriteReal(ix-m2); WriteReal(iy+m4);
         WriteReal(ix+m2); WriteReal(iy+m4);
         WriteReal(ix+m2); WriteReal(iy)   ; PrintFast(2," c");
         WriteReal(ix+m2); WriteReal(iy-m4);
         WriteReal(ix-m2); WriteReal(iy-m4);
         WriteReal(ix-m2); WriteReal(iy)   ; PrintFast(4," c h");
      // Square
      } else if (ms == 25 || ms == 21) {
         WriteReal(ix-m2); WriteReal(iy-m2);
         WriteReal(m)    ; WriteReal(m)    ; PrintFast(3," re");
      // Down triangle
      } else if (ms == 23 || ms == 32) {
         MoveTo(ix   , iy-m2);
         LineTo(ix+m2, iy+m2);
         LineTo(ix-m2, iy+m2);
         PrintFast(2," h");
      // Up triangle
      } else if (ms == 26 || ms == 22) {
         MoveTo(ix-m2, iy-m2);
         LineTo(ix+m2, iy-m2);
         LineTo(ix   , iy+m2);
         PrintFast(2," h");
      } else if (ms == 27 || ms == 33) {
         MoveTo(ix   , iy-m2);
         LineTo(ix+m3, iy);
         LineTo(ix   , iy+m2);
         LineTo(ix-m3, iy)   ;
         PrintFast(2," h");
      } else if (ms == 28 || ms == 34) {
         MoveTo(ix-m6, iy-m6);
         LineTo(ix-m6, iy-m2);
         LineTo(ix+m6, iy-m2);
         LineTo(ix+m6, iy-m6);
         LineTo(ix+m2, iy-m6);
         LineTo(ix+m2, iy+m6);
         LineTo(ix+m6, iy+m6);
         LineTo(ix+m6, iy+m2);
         LineTo(ix-m6, iy+m2);
         LineTo(ix-m6, iy+m6);
         LineTo(ix-m2, iy+m6);
         LineTo(ix-m2, iy-m6);
         PrintFast(2," h");
      } else if (ms == 29 || ms == 30) {
         MoveTo(ix           , iy-m2);
         LineTo(ix-0.112255*m, iy-0.15451*m);
         LineTo(ix-0.47552*m , iy-0.15451*m);
         LineTo(ix-0.181635*m, iy+0.05902*m);
         LineTo(ix-0.29389*m , iy+0.40451*m);
         LineTo(ix           , iy+0.19098*m);
         LineTo(ix+0.29389*m , iy+0.40451*m);
         LineTo(ix+0.181635*m, iy+0.05902*m);
         LineTo(ix+0.47552*m , iy-0.15451*m);
         LineTo(ix+0.112255*m, iy-0.15451*m);
         PrintFast(2," h");
      } else if (ms == 35 ) {
         MoveTo(ix-m2, iy   );
         LineTo(ix   , iy-m2);
         LineTo(ix+m2, iy   );
         LineTo(ix   , iy+m2);
         LineTo(ix-m2, iy   );
         LineTo(ix+m2, iy   );
         LineTo(ix   , iy+m2);
         LineTo(ix   , iy-m2);
         PrintFast(2," h");
      } else if (ms == 36 ) {
         MoveTo(ix-m2, iy-m2);
         LineTo(ix+m2, iy-m2);
         LineTo(ix+m2, iy+m2);
         LineTo(ix-m2, iy+m2);
         LineTo(ix-m2, iy-m2);
         LineTo(ix+m2, iy+m2);
         LineTo(ix-m2, iy+m2);
         LineTo(ix+m2, iy-m2);
         PrintFast(2," h");
      } else if (ms == 37 || ms == 39 ) {
         MoveTo(ix   , iy   );
         LineTo(ix-m8, iy+m2);
         LineTo(ix-m2, iy   );
         LineTo(ix   , iy   );
         LineTo(ix-m8, iy-m2);
         LineTo(ix+m8, iy-m2);
         LineTo(ix   , iy   );
         LineTo(ix+m2, iy   );
         LineTo(ix+m8, iy+m2);
         LineTo(ix   , iy   );
         PrintFast(2," h");
      } else if (ms == 38 ) {
         MoveTo(ix-m2, iy   );
         LineTo(ix-m2, iy-m8);
         LineTo(ix-m8, iy-m2);
         LineTo(ix+m8, iy-m2);
         LineTo(ix+m2, iy-m8);
         LineTo(ix+m2, iy+m8);
         LineTo(ix+m8, iy+m2);
         LineTo(ix-m8, iy+m2);
         LineTo(ix-m2, iy+m8);
         LineTo(ix-m2, iy   );
         LineTo(ix+m2, iy   );
         LineTo(ix   , iy   );
         LineTo(ix   , iy-m2);
         LineTo(ix   , iy+m2);
         LineTo(ix   , iy   );
         PrintFast(2," h");
      } else if (ms == 40 || ms == 41 ) {
         MoveTo(ix   , iy   );
         LineTo(ix+m8, iy+m2);
         LineTo(ix+m2, iy+m8);
         LineTo(ix   , iy   );
         LineTo(ix+m2, iy-m8);
         LineTo(ix+m8, iy-m2);
         LineTo(ix   , iy   );
         LineTo(ix-m8, iy-m2);
         LineTo(ix-m2, iy-m8);
         LineTo(ix   , iy   );
         LineTo(ix-m2, iy+m8);
         LineTo(ix-m8, iy+m2);
         LineTo(ix   , iy   );
         PrintFast(2," h");
      } else if (ms == 42 || ms == 43 ) {
         MoveTo(ix   , iy+m2);
         LineTo(ix-m9, iy+m9);
         LineTo(ix-m2, iy   );
         LineTo(ix-m9, iy-m9);
         LineTo(ix   , iy-m2);
         LineTo(ix+m9, iy-m9);
         LineTo(ix+m2, iy   );
         LineTo(ix+m9, iy+m9);
         LineTo(ix   , iy+m2);
         PrintFast(2," h");
      } else if (ms == 44 ) {
         MoveTo(ix   , iy   );
         LineTo(ix+m8, iy+m2);
         LineTo(ix-m8, iy+m2);
         LineTo(ix+m8, iy-m2);
         LineTo(ix-m8, iy-m2);
         LineTo(ix   , iy   );
         LineTo(ix+m2, iy+m8);
         LineTo(ix+m2, iy-m8);
         LineTo(ix-m2, iy+m8);
         LineTo(ix-m2, iy-m8);
         LineTo(ix   , iy   );
         PrintFast(2," h");
      } else if (ms == 45 ) {
         MoveTo(ix+m6/2., iy+m6/2.);
         LineTo(ix+m8, iy+m2);
         LineTo(ix-m8, iy+m2);
         LineTo(ix-m6/2., iy+m6/2.);
         LineTo(ix-m2, iy+m8);
         LineTo(ix-m2, iy-m8);
         LineTo(ix-m6/2., iy-m6/2.);
         LineTo(ix-m8, iy-m2);
         LineTo(ix+m8, iy-m2);
         LineTo(ix+m6/2., iy-m6/2.);
         LineTo(ix+m2, iy-m8);
         LineTo(ix+m2, iy+m8);
         LineTo(ix+m6/2., iy+m6/2.);
         PrintFast(2," h");
      } else if (ms == 46 || ms == 47 ) {
         MoveTo(ix   , iy+m8);
         LineTo(ix-m8, iy+m2);
         LineTo(ix-m2, iy+m8);
         LineTo(ix-m8, iy   );
         LineTo(ix-m2, iy-m8);
         LineTo(ix-m8, iy-m2);
         LineTo(ix   , iy-m8);
         LineTo(ix+m8, iy-m2);
         LineTo(ix+m2, iy-m8);
         LineTo(ix+m8, iy   );
         LineTo(ix+m2, iy+m8);
         LineTo(ix+m8, iy+m2);
         LineTo(ix   , iy+m8);
         PrintFast(2," h");
      } else if (ms == 48 ) {
         MoveTo(ix   , iy+m8*1.005);
         LineTo(ix-m8, iy+m2);
         LineTo(ix-m2, iy+m8);
         LineTo(ix-m8, iy   );
         LineTo(ix-m2, iy-m8);
         LineTo(ix-m8, iy-m2);
         LineTo(ix   , iy-m8);
         LineTo(ix+m8, iy-m2);
         LineTo(ix+m2, iy-m8);
         LineTo(ix+m8, iy   );
         LineTo(ix+m2, iy+m8);
         LineTo(ix+m8, iy+m2);
         LineTo(ix   , iy+m8*0.995);
         LineTo(ix+m8*0.995, iy   );
         LineTo(ix   , iy-m8*0.995);
         LineTo(ix-m8*0.995, iy   );
         LineTo(ix   , iy+m8*0.995);
         PrintFast(2," h");
      } else if (ms == 49 ) {
         MoveTo(ix-m6, iy-m6*1.01);
         LineTo(ix-m6, iy-m2);
         LineTo(ix+m6, iy-m2);
         LineTo(ix+m6, iy-m6);
         LineTo(ix+m2, iy-m6);
         LineTo(ix+m2, iy+m6);
         LineTo(ix+m6, iy+m6);
         LineTo(ix+m6, iy+m2);
         LineTo(ix-m6, iy+m2);
         LineTo(ix-m6, iy+m6);
         LineTo(ix-m2, iy+m6);
         LineTo(ix-m2, iy-m6);
         LineTo(ix-m6, iy-m6*0.99);
         LineTo(ix-m6, iy+m6);
         LineTo(ix+m6, iy+m6);
         LineTo(ix+m6, iy-m6);
         MoveTo(ix-m6, iy-m6*1.01);
         PrintFast(2," h");
      } else {
         MoveTo(ix-1, iy);
         LineTo(ix  , iy);
      }
      if ((ms > 19 && ms < 24) || ms == 29 || ms == 33 || ms == 34 ||
          ms == 39 || ms == 41 || ms == 43 || ms == 45 ||
          ms == 47 || ms == 48 || ms == 49) {
         PrintFast(2," f");
      } else {
         PrintFast(2," S");
      }
   }

   SetLineStyle(linestylesav);
   SetLineWidth(linewidthsav);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a PolyLine
///
///  Draw a polyline through the points xw,yw.
///
///  - If nn=1 moves only to point xw,yw.
///  - If nn=0 the XW(1) and YW(1) are  written  in the PDF file
///            according to the current NT.
///  - If nn>0 the line is clipped as a line.
///  - If nn<0 the line is clipped as a fill area.

void TPDF::DrawPS(Int_t nn, Float_t *xw, Float_t *yw)
{
   static Float_t dyhatch[24] = {.0075,.0075,.0075,.0075,.0075,.0075,.0075,.0075,
                                 .01  ,.01  ,.01  ,.01  ,.01  ,.01  ,.01  ,.01  ,
                                 .015 ,.015 ,.015 ,.015 ,.015 ,.015 ,.015 ,.015};
   static Float_t anglehatch[24] = {180, 90,135, 45,150, 30,120, 60,
                                    180, 90,135, 45,150, 30,120, 60,
                                    180, 90,135, 45,150, 30,120, 60};
   Int_t  n = 0, fais = 0 , fasi = 0;

   Style_t linestylesav = fLineStyle;
   Width_t linewidthsav = fLineWidth;

   if (nn > 0) {
      if (fLineWidth<=0) return;
      n = nn;
      SetLineStyle(fLineStyle);
      SetLineWidth(fLineWidth);
      SetColor(Int_t(fLineColor));
   }
   if (nn < 0) {
      n = -nn;
      SetLineStyle(1);
      SetLineWidth(1);
      SetColor(Int_t(fFillColor));
      fais = fFillStyle/1000;
      fasi = fFillStyle%1000;
      if (fais == 3 || fais == 2) {
         if (fasi > 100 && fasi <125) {
            DrawHatch(dyhatch[fasi-101],anglehatch[fasi-101], n, xw, yw);
            SetLineStyle(linestylesav);
            SetLineWidth(linewidthsav);
            return;
         }
         if (fasi > 0 && fasi < 26) {
            SetFillPatterns(fasi, Int_t(fFillColor));
         }
      }
   }

   WriteReal(XtoPDF(xw[0]));
   WriteReal(YtoPDF(yw[0]));
   if (n <= 1) {
      if (n == 0) return;
      PrintFast(2," m");
      return;
   }

   PrintFast(2," m");

   for (Int_t i=1;i<n;i++) LineTo(XtoPDF(xw[i]), YtoPDF(yw[i]));

   if (nn > 0) {
      if (xw[0] == xw[n-1] && yw[0] == yw[n-1]) PrintFast(2," h");
      PrintFast(2," S");
   } else {
      if (fais == 0) {PrintFast(2," s"); return;}
      if (fais == 3 || fais == 2) {
         if (fasi > 0 && fasi < 26) {
            PrintFast(3," f*");
            fRed   = -1;
            fGreen = -1;
            fBlue  = -1;
            fAlpha = -1.;
         }
         SetLineStyle(linestylesav);
         SetLineWidth(linewidthsav);
         return;
      }
      PrintFast(3," f*");
   }

   SetLineStyle(linestylesav);
   SetLineWidth(linewidthsav);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a PolyLine
///
/// Draw a polyline through  the points xw,yw.
///
///  - If nn=1 moves only to point xw,yw.
///  - If nn=0 the xw(1) and YW(1) are  written  in the PDF file
///            according to the current NT.
///  - If nn>0 the line is clipped as a line.
///  - If nn<0 the line is clipped as a fill area.

void TPDF::DrawPS(Int_t nn, Double_t *xw, Double_t *yw)
{
   static Float_t dyhatch[24] = {.0075,.0075,.0075,.0075,.0075,.0075,.0075,.0075,
                                 .01  ,.01  ,.01  ,.01  ,.01  ,.01  ,.01  ,.01  ,
                                 .015 ,.015 ,.015 ,.015 ,.015 ,.015 ,.015 ,.015};
   static Float_t anglehatch[24] = {180, 90,135, 45,150, 30,120, 60,
                                    180, 90,135, 45,150, 30,120, 60,
                                    180, 90,135, 45,150, 30,120, 60};
   Int_t  n = 0, fais = 0, fasi = 0;

   Style_t linestylesav = fLineStyle;
   Width_t linewidthsav = fLineWidth;

   if (nn > 0) {
      if (fLineWidth<=0) return;
      n = nn;
      SetLineStyle(fLineStyle);
      SetLineWidth(fLineWidth);
      SetColor(Int_t(fLineColor));
   }
   if (nn < 0) {
      n = -nn;
      SetLineStyle(1);
      SetLineWidth(1);
      SetColor(Int_t(fFillColor));
      fais = fFillStyle/1000;
      fasi = fFillStyle%1000;
      if (fais == 3 || fais == 2) {
         if (fasi > 100 && fasi <125) {
            DrawHatch(dyhatch[fasi-101],anglehatch[fasi-101], n, xw, yw);
            SetLineStyle(linestylesav);
            SetLineWidth(linewidthsav);
            return;
         }
         if (fasi > 0 && fasi < 26) {
            SetFillPatterns(fasi, Int_t(fFillColor));
         }
      }
   }

   WriteReal(XtoPDF(xw[0]));
   WriteReal(YtoPDF(yw[0]));
   if (n <= 1) {
      if (n == 0) return;
      PrintFast(2," m");
      return;
   }

   PrintFast(2," m");

   for (Int_t i=1;i<n;i++) LineTo(XtoPDF(xw[i]), YtoPDF(yw[i]));

   if (nn > 0) {
      if (xw[0] == xw[n-1] && yw[0] == yw[n-1]) PrintFast(2," h");
      PrintFast(2," S");
   } else {
      if (fais == 0) {PrintFast(2," s"); return;}
      if (fais == 3 || fais == 2) {
         if (fasi > 0 && fasi < 26) {
            PrintFast(3," f*");
            fRed   = -1;
            fGreen = -1;
            fBlue  = -1;
            fAlpha = -1.;
         }
         SetLineStyle(linestylesav);
         SetLineWidth(linewidthsav);
         return;
      }
      PrintFast(3," f*");
   }

   SetLineStyle(linestylesav);
   SetLineWidth(linewidthsav);
}

////////////////////////////////////////////////////////////////////////////////
/// Font encoding

void TPDF::FontEncode()
{
   static const char *sdtfonts[] = {
   "/Times-Italic"         , "/Times-Bold"         , "/Times-BoldItalic",
   "/Helvetica"            , "/Helvetica-Oblique"  , "/Helvetica-Bold"  ,
   "/Helvetica-BoldOblique", "/Courier"            , "/Courier-Oblique" ,
   "/Courier-Bold"         , "/Courier-BoldOblique", "/Symbol"          ,
   "/Times-Roman"          , "/ZapfDingbats"       , "/Symbol"};

   for (Int_t i=0; i<kNumberOfFonts; i++) {
      NewObject(kObjFont+i);
      PrintStr("<<@");
      PrintStr("/Type /Font@");
      PrintStr("/Subtype /Type1@");
      PrintStr("/Name /F");
      WriteInteger(i+1,0);
      PrintStr("@");
      PrintStr("/BaseFont ");
      PrintStr(sdtfonts[i]);
      PrintStr("@");
      if (i!=11 && i!=13 && i!=14) {
         PrintStr("/Encoding /WinAnsiEncoding");
         PrintStr("@");
      }
      PrintStr(">>@");
      PrintStr("endobj@");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a line to a new position

void TPDF::LineTo(Double_t x, Double_t y)
{
   WriteReal(x);
   WriteReal(y);
   PrintFast(2," l");
}

////////////////////////////////////////////////////////////////////////////////
/// Move to a new position

void TPDF::MoveTo(Double_t x, Double_t y)
{
   WriteReal(x);
   WriteReal(y);
   PrintFast(2," m");
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new object in the PDF file

void TPDF::NewObject(Int_t n)
{
   if (!fObjPos || n >= fObjPosSize) {
      Int_t newN = TMath::Max(2*fObjPosSize,n+1);
      Int_t *saveo = new Int_t [newN];
      if (fObjPos && fObjPosSize) {
         memcpy(saveo,fObjPos,fObjPosSize*sizeof(Int_t));
         memset(&saveo[fObjPosSize],0,(newN-fObjPosSize)*sizeof(Int_t));
         delete [] fObjPos;
      }
      fObjPos     = saveo;
      fObjPosSize = newN;
   }
   fObjPos[n-1] = fNByte;
   fNbObj       = TMath::Max(fNbObj,n);
   WriteInteger(n, 0);
   PrintStr(" 0 obj");
   PrintStr("@");
}

////////////////////////////////////////////////////////////////////////////////
/// Start a new PDF page.

void TPDF::NewPage()
{
   if (!fPageNotEmpty) return;

   // Compute pad conversion coefficients
   if (gPad) {
      Double_t ww   = gPad->GetWw();
      Double_t wh   = gPad->GetWh();
      fYsize        = fXsize*wh/ww;
   } else {
      fYsize = 27;
   }

   fNbPage++;

   if (fNbPage>1) {
      // Close the currently opened page
      WriteCompressedBuffer();
      PrintStr("endstream@");
      Int_t streamLength = fNByte-fStartStream-10;
      PrintStr("endobj@");
      NewObject(4*(fNbPage-2)+kObjFirstPage+2);
      WriteInteger(streamLength, 0);
      PrintStr("@");
      PrintStr("endobj@");
      NewObject(4*(fNbPage-2)+kObjFirstPage+3);
      PrintStr("<<@");
      if (!strstr(GetTitle(),"PDF")) {
         PrintStr("/Title (");
         PrintStr(GetTitle());
         PrintStr(")@");
      } else {
         PrintStr("/Title (Page");
         WriteInteger(fNbPage-1);
         PrintStr(")@");
      }
      PrintStr("/Dest [");
      WriteInteger(4*(fNbPage-2)+kObjFirstPage);
      PrintStr(" 0 R /XYZ null null 0]@");
      PrintStr("/Parent");
      WriteInteger(kObjContents);
      PrintStr(" 0 R");
      PrintStr("@");
      PrintStr("/Next");
      WriteInteger(4*(fNbPage-1)+kObjFirstPage+3);
      PrintStr(" 0 R");
      PrintStr("@");
      if (fNbPage>2) {
         PrintStr("/Prev");
         WriteInteger(4*(fNbPage-3)+kObjFirstPage+3);
         PrintStr(" 0 R");
         PrintStr("@");
      }
      PrintStr(">>@");
   }

   // Start a new page
   NewObject(4*(fNbPage-1)+kObjFirstPage);
   PrintStr("<<@");
   PrintStr("/Type /Page@");
   PrintStr("@");
   PrintStr("/Parent");
   WriteInteger(kObjPages);
   PrintStr(" 0 R");
   PrintStr("@");

   Double_t xlow=0, ylow=0, xup=1, yup=1;
   if (gPad) {
      xlow = gPad->GetAbsXlowNDC();
      xup  = xlow + gPad->GetAbsWNDC();
      ylow = gPad->GetAbsYlowNDC();
      yup  = ylow + gPad->GetAbsHNDC();
   }

   PrintStr("/MediaBox [");
   Double_t width, height;
   switch (fPageFormat) {
      case 100 :
         width  = 8.5*2.54;
         height = 11.*2.54;
         break;
      case 200 :
         width  = 8.5*2.54;
         height = 14.*2.54;
         break;
      case 300 :
         width  = 11.*2.54;
         height = 17.*2.54;
         break;
      default  :
         width  = 21.0*TMath::Power(TMath::Sqrt(2.), 4-fPageFormat);
         height = 29.7*TMath::Power(TMath::Sqrt(2.), 4-fPageFormat);
   };
   WriteReal(CMtoPDF(fXsize*xlow));
   WriteReal(CMtoPDF(fYsize*ylow));
   WriteReal(CMtoPDF(width));
   WriteReal(CMtoPDF(height));
   PrintStr("]");
   PrintStr("@");

   Double_t xmargin = CMtoPDF(0.7);
   Double_t ymargin = 0;
   if (fPageOrientation == 1) ymargin = CMtoPDF(TMath::Sqrt(2.)*0.7);
   if (fPageOrientation == 2) ymargin = CMtoPDF(height)-CMtoPDF(0.7);

   PrintStr("/CropBox [");
   if (fPageOrientation == 1) {
      WriteReal(xmargin);
      WriteReal(ymargin);
      WriteReal(xmargin+CMtoPDF(fXsize*xup));
      WriteReal(ymargin+CMtoPDF(fYsize*yup));
   }
   if (fPageOrientation == 2) {
      WriteReal(xmargin);
      WriteReal(CMtoPDF(height)-CMtoPDF(fXsize*xup)-xmargin);
      WriteReal(xmargin+CMtoPDF(fYsize*yup));
      WriteReal(CMtoPDF(height)-xmargin);
   }
   PrintStr("]");
   PrintStr("@");

   if (fPageOrientation == 1) PrintStr("/Rotate 0@");
   if (fPageOrientation == 2) PrintStr("/Rotate 90@");

   PrintStr("/Resources");
   WriteInteger(kObjPageResources);
   PrintStr(" 0 R");
   PrintStr("@");

   PrintStr("/Contents");
   WriteInteger(4*(fNbPage-1)+kObjFirstPage+1);
   PrintStr(" 0 R@");
   PrintStr(">>@");
   PrintStr("endobj@");

   NewObject(4*(fNbPage-1)+kObjFirstPage+1);
   PrintStr("<<@");
   PrintStr("/Length");
   WriteInteger(4*(fNbPage-1)+kObjFirstPage+2);
   PrintStr(" 0 R@");
   PrintStr("/Filter [/FlateDecode]@");
   PrintStr(">>@");
   PrintStr("stream@");
   fStartStream = fNByte;
   fCompress = kTRUE;

   // Force the line width definition next time TPDF::SetLineWidth will be called.
   fLineWidth = -1;

   // Force the color definition next time TPDF::SetColor will be called.
   fRed   = -1;
   fGreen = -1;
   fBlue  = -1;
   fAlpha = -1.;

   PrintStr("1 0 0 1");
   if (fPageOrientation == 2) {
      ymargin = CMtoPDF(height)-CMtoPDF(fXsize*xup)-xmargin;
      xmargin = xmargin+CMtoPDF(fYsize*yup);
   }
   WriteReal(xmargin);
   WriteReal(ymargin);
   PrintStr(" cm");
   if (fPageOrientation == 2) PrintStr(" 0 1 -1 0 0 0 cm");
   if (fgLineJoin) {
      WriteInteger(fgLineJoin);
      PrintFast(2," j");
   }
   if (fgLineCap) {
      WriteInteger(fgLineCap);
      PrintFast(2," J");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Deactivate an already open PDF file

void TPDF::Off()
{
   gVirtualPS = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Activate an already open PDF file

void TPDF::On()
{
   // fType is used to know if the PDF file is open. Unlike TPostScript, TPDF
   // has no "workstation type".

   if (!fType) {
      Error("On", "no PDF file open");
      Off();
      return;
   }
   gVirtualPS = this;
}

////////////////////////////////////////////////////////////////////////////////
/// Open a PDF file

void TPDF::Open(const char *fname, Int_t wtype)
{
   Int_t i;

   if (fStream) {
      Warning("Open", "PDF file already open");
      return;
   }

   fLenBuffer = 0;
   fRed       = -1;
   fGreen     = -1;
   fBlue      = -1;
   fAlpha     = -1.;
   fType      = abs(wtype);
   SetLineJoin(gStyle->GetJoinLinePS());
   SetLineCap(gStyle->GetCapLinePS());
   SetLineScale(gStyle->GetLineScalePS()/4.);
   gStyle->GetPaperSize(fXsize, fYsize);
   Float_t xrange, yrange;
   if (gPad) {
      Double_t ww = gPad->GetWw();
      Double_t wh = gPad->GetWh();
      if (fType == 113) {
         ww *= gPad->GetWNDC();
         wh *= gPad->GetHNDC();
      }
      Double_t ratio = wh/ww;
      xrange = fXsize;
      yrange = fXsize*ratio;
      if (yrange > fYsize) { yrange = fYsize; xrange = yrange/ratio;}
      fXsize = xrange; fYsize = yrange;
   }

   // Open OS file
   fStream = new std::ofstream();
#ifdef R__WIN32
      fStream->open(fname, std::ofstream::out | std::ofstream::binary);
#else
      fStream->open(fname, std::ofstream::out);
#endif
   if (fStream == 0 || !fStream->good()) {
      printf("ERROR in TPDF::Open: Cannot open file:%s\n",fname);
      if (fStream == 0) return;
   }

   gVirtualPS = this;

   for (i=0; i<fSizBuffer; i++) fBuffer[i] = ' ';

   // The page orientation is last digit of PDF workstation type
   //  orientation = 1 for portrait
   //  orientation = 2 for landscape
   fPageOrientation = fType%10;
   if (fPageOrientation < 1 || fPageOrientation > 2) {
      Error("Open", "Invalid page orientation %d", fPageOrientation);
      return;
   }

   // format = 0-99 is the European page format (A4,A3 ...)
   // format = 100 is the US format  8.5x11.0 inch
   // format = 200 is the US format  8.5x14.0 inch
   // format = 300 is the US format 11.0x17.0 inch
   fPageFormat = fType/1000;
   if (fPageFormat == 0)  fPageFormat = 4;
   if (fPageFormat == 99) fPageFormat = 0;

   fRange = kFALSE;

   // Set a default range
   Range(fXsize, fYsize);

   fObjPos = 0;
   fObjPosSize = 0;
   fNbObj = 0;
   fNbPage = 0;

   PrintStr("%PDF-1.4@");
   PrintStr("%\342\343\317\323");
   PrintStr("@");

   NewObject(kObjRoot);
   PrintStr("<<@");
   PrintStr("/Type /Catalog@");
   PrintStr("/Pages");
   WriteInteger(kObjPages);
   PrintStr(" 0 R@");
   PrintStr("/Outlines");
   WriteInteger(kObjOutlines);
   PrintStr(" 0 R@");
   PrintStr("/PageMode /UseOutlines@");
   PrintStr(">>@");
   PrintStr("endobj@");

   NewObject(kObjInfo);
   PrintStr("<<@");
   PrintStr("/Creator (ROOT Version ");
   PrintStr(gROOT->GetVersion());
   PrintStr(")");
   PrintStr("@");
   PrintStr("/CreationDate (");
   TDatime t;
   Int_t toff = t.Convert(kFALSE) - t.Convert(kTRUE); // time zone and dst offset
   toff = toff/60;
   char str[24];
   snprintf(str,24,"D:%4.4d%2.2d%2.2d%2.2d%2.2d%2.2d%c%2.2d'%2.2d'",
            t.GetYear()  , t.GetMonth(),
            t.GetDay()   , t.GetHour(),
            t.GetMinute(), t.GetSecond(),
            toff < 0 ? '-' : '+',
            // TMath::Abs(toff/60), TMath::Abs(toff%60)); // format-truncation warning
            TMath::Abs(toff/60) & 0x3F, TMath::Abs(toff%60) & 0x3F); // now 2 digits
   PrintStr(str);
   PrintStr(")");
   PrintStr("@");
   PrintStr("/ModDate (");
   PrintStr(str);
   PrintStr(")");
   PrintStr("@");
   PrintStr("/Title (");
   if (strlen(GetName())<=80) PrintStr(GetName());
   PrintStr(")");
   PrintStr("@");
   PrintStr("/Keywords (ROOT)@");
   PrintStr(">>@");
   PrintStr("endobj@");

   NewObject(kObjPageResources);
   PrintStr("<<@");
   PrintStr("/ProcSet [/PDF /Text]@");

   PrintStr("/Font@");
   PrintStr("<<@");
   for (i=0; i<kNumberOfFonts; i++) {
      PrintStr(" /F");
      WriteInteger(i+1,0);
      WriteInteger(kObjFont+i);
      PrintStr(" 0 R");
   }
   PrintStr("@");
   PrintStr(">>@");

   PrintStr("/ExtGState");
   WriteInteger(kObjTransList);
   PrintStr(" 0 R @");
   if (fAlphas.size()) fAlphas.clear();

   PrintStr("/ColorSpace << /Cs8");
   WriteInteger(kObjColorSpace);
   PrintStr(" 0 R >>");
   PrintStr("@");
   PrintStr("/Pattern");
   WriteInteger(kObjPatternList);
   PrintStr(" 0 R");
   PrintStr("@");
   PrintStr(">>@");
   PrintStr("endobj@");

   FontEncode();
   PatternEncode();

   NewPage();
   fPageNotEmpty = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Output the string str in the output buffer

void TPDF::PrintStr(const char *str)
{
   Int_t len = strlen(str);
   if (len == 0) return;
   fPageNotEmpty = kTRUE;

   if (fCompress) {
      if (fLenBuffer+len >= fSizBuffer) {
         fBuffer  = TStorage::ReAllocChar(fBuffer, 2*fSizBuffer, fSizBuffer);
         fSizBuffer = 2*fSizBuffer;
      }
      strcpy(fBuffer + fLenBuffer, str);
      fLenBuffer += len;
      return;
   }

   TVirtualPS::PrintStr(str);
}

////////////////////////////////////////////////////////////////////////////////
/// Fast version of Print

void TPDF::PrintFast(Int_t len, const char *str)
{
   fPageNotEmpty = kTRUE;
   if (fCompress) {
      if (fLenBuffer+len >= fSizBuffer) {
         fBuffer  = TStorage::ReAllocChar(fBuffer, 2*fSizBuffer, fSizBuffer);
         fSizBuffer = 2*fSizBuffer;
      }
      strcpy(fBuffer + fLenBuffer, str);
      fLenBuffer += len;
      return;
   }

   TVirtualPS::PrintFast(len, str);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the range for the paper in centimetres

void TPDF::Range(Float_t xsize, Float_t ysize)
{
   Float_t xps, yps, xncm, yncm, dxwn, dywn, xwkwn, ywkwn, xymax;

   fXsize = xsize;
   fYsize = ysize;

   xps = xsize;
   yps = ysize;

   if (xsize <= xps && ysize < yps) {
      if ( xps > yps) xymax = xps;
      else            xymax = yps;
      xncm  = xsize/xymax;
      yncm  = ysize/xymax;
      dxwn  = ((xps/xymax)-xncm)/2;
      dywn  = ((yps/xymax)-yncm)/2;
   } else {
      if (xps/yps < 1) xwkwn = xps/yps;
      else             xwkwn = 1;
      if (yps/xps < 1) ywkwn = yps/xps;
      else             ywkwn = 1;

      if (xsize < ysize)  {
         xncm = ywkwn*xsize/ysize;
         yncm = ywkwn;
         dxwn = (xwkwn-xncm)/2;
         dywn = 0;
         if (dxwn < 0) {
            xncm = xwkwn;
            dxwn = 0;
            yncm = xwkwn*ysize/xsize;
            dywn = (ywkwn-yncm)/2;
         }
      } else {
         xncm = xwkwn;
         yncm = xwkwn*ysize/xsize;
         dxwn = 0;
         dywn = (ywkwn-yncm)/2;
         if (dywn < 0) {
            yncm = ywkwn;
            dywn = 0;
            xncm = ywkwn*xsize/ysize;
            dxwn = (xwkwn-xncm)/2;
         }
      }
   }
   fRange = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the alpha channel value.

void TPDF::SetAlpha(Float_t a)
{
   if (a == fAlpha) return;
   fAlpha = a;
   if (fAlpha  <= 0.000001) fAlpha  = 0;

   Bool_t known = kFALSE;
   for (int i=0; i<(int)fAlphas.size(); i++) {
      if (fAlpha == fAlphas[i]) {
         known = kTRUE;
         break;
      }
   }
   if (!known) fAlphas.push_back(fAlpha);
   PrintStr(Form(" /ca%3.2f gs /CA%3.2f gs",fAlpha,fAlpha));
}

////////////////////////////////////////////////////////////////////////////////
/// Set color with its color index.

void TPDF::SetColor(Int_t color)
{
   if (color < 0) color = 0;
   TColor *col = gROOT->GetColor(color);

   if (col) {
      SetColor(col->GetRed(), col->GetGreen(), col->GetBlue());
      SetAlpha(col->GetAlpha());
   } else {
      SetColor(1., 1., 1.);
      SetAlpha(1.);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set color with its R G B components:
///
///  - r: % of red in [0,1]
///  - g: % of green in [0,1]
///  - b: % of blue in [0,1]

void TPDF::SetColor(Float_t r, Float_t g, Float_t b)
{
   if (r == fRed && g == fGreen && b == fBlue) return;

   fRed   = r;
   fGreen = g;
   fBlue  = b;
   if (fRed   <= 0.000001) fRed   = 0;
   if (fGreen <= 0.000001) fGreen = 0;
   if (fBlue  <= 0.000001) fBlue  = 0;

   if (gStyle->GetColorModelPS()) {
      Double_t colCyan, colMagenta, colYellow;
      Double_t colBlack = TMath::Min(TMath::Min(1-fRed,1-fGreen),1-fBlue);
      if (colBlack==1) {
         colCyan    = 0;
         colMagenta = 0;
         colYellow  = 0;
      } else {
         colCyan    = (1-fRed-colBlack)/(1-colBlack);
         colMagenta = (1-fGreen-colBlack)/(1-colBlack);
         colYellow  = (1-fBlue-colBlack)/(1-colBlack);
      }
      if (colCyan    <= 0.000001) colCyan    = 0;
      if (colMagenta <= 0.000001) colMagenta = 0;
      if (colYellow  <= 0.000001) colYellow  = 0;
      if (colBlack   <= 0.000001) colBlack   = 0;
      WriteReal(colCyan);
      WriteReal(colMagenta);
      WriteReal(colYellow);
      WriteReal(colBlack);
      PrintFast(2," K");
      WriteReal(colCyan);
      WriteReal(colMagenta);
      WriteReal(colYellow);
      WriteReal(colBlack);
      PrintFast(2," k");
   } else {
      WriteReal(fRed);
      WriteReal(fGreen);
      WriteReal(fBlue);
      PrintFast(3," RG");
      WriteReal(fRed);
      WriteReal(fGreen);
      WriteReal(fBlue);
      PrintFast(3," rg");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set color index for fill areas

void TPDF::SetFillColor( Color_t cindex )
{
   fFillColor = cindex;
   if (gStyle->GetFillColor() <= 0) cindex = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the fill patterns (1 to 25) for fill areas

void TPDF::SetFillPatterns(Int_t ipat, Int_t color)
{
   char cpat[10];
   TColor *col = gROOT->GetColor(color);
   if (!col) return;
   PrintStr(" /Cs8 cs");
   Double_t colRed   = col->GetRed();
   Double_t colGreen = col->GetGreen();
   Double_t colBlue  = col->GetBlue();
   if (gStyle->GetColorModelPS()) {
      Double_t colBlack = TMath::Min(TMath::Min(1-colRed,1-colGreen),1-colBlue);
      if (colBlack==1) {
         WriteReal(0);
         WriteReal(0);
         WriteReal(0);
         WriteReal(colBlack);
      } else {
         Double_t colCyan    = (1-colRed-colBlack)/(1-colBlack);
         Double_t colMagenta = (1-colGreen-colBlack)/(1-colBlack);
         Double_t colYellow  = (1-colBlue-colBlack)/(1-colBlack);
         WriteReal(colCyan);
         WriteReal(colMagenta);
         WriteReal(colYellow);
         WriteReal(colBlack);
      }
   } else {
      WriteReal(colRed);
      WriteReal(colGreen);
      WriteReal(colBlue);
   }
   snprintf(cpat,10," /P%2.2d scn", ipat);
   PrintStr(cpat);
}

////////////////////////////////////////////////////////////////////////////////
/// Set color index for lines

void TPDF::SetLineColor( Color_t cindex )
{
   fLineColor = cindex;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the value of the global parameter TPDF::fgLineJoin.
/// This parameter determines the appearance of joining lines in a PDF
/// output.
/// It takes one argument which may be:
///   - 0 (miter join)
///   - 1 (round join)
///   - 2 (bevel join)
/// The default value is 0 (miter join).
///
/// \image html postscript_1.png
///
/// To change the line join behaviour just do:
/// ~~~ {.cpp}
/// gStyle->SetJoinLinePS(2); // Set the PDF line join to bevel.
/// ~~~

void TPDF::SetLineJoin( Int_t linejoin )
{
   fgLineJoin = linejoin;
   if (fgLineJoin<0) fgLineJoin=0;
   if (fgLineJoin>2) fgLineJoin=2;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the value of the global parameter TPDF::fgLineCap.
/// This parameter determines the appearance of line caps in a PDF
/// output.
/// It takes one argument which may be:
///   - 0 (butt caps)
///   - 1 (round caps)
///   - 2 (projecting caps)
/// The default value is 0 (butt caps).
///
/// \image html postscript_2.png
///
/// To change the line cap behaviour just do:
/// ~~~ {.cpp}
/// gStyle->SetCapLinePS(2); // Set the PDF line cap to projecting.
/// ~~~

void TPDF::SetLineCap( Int_t linecap )
{
   fgLineCap = linecap;
   if (fgLineCap<0) fgLineCap=0;
   if (fgLineCap>2) fgLineCap=2;
}

////////////////////////////////////////////////////////////////////////////////
/// Change the line style
///
///  - linestyle = 2 dashed
///  - linestyle = 3 dotted
///  - linestyle = 4 dash-dotted
///  - linestyle = else solid (1 in is used most of the time)

void TPDF::SetLineStyle(Style_t linestyle)
{
   if ( linestyle == fLineStyle) return;
   fLineStyle = linestyle;
   TString st = (TString)gStyle->GetLineStyleString(linestyle);
   PrintFast(2," [");
   TObjArray *tokens = st.Tokenize(" ");
   for (Int_t j = 0; j<tokens->GetEntries(); j++) {
      Int_t it;
      sscanf(((TObjString*)tokens->At(j))->GetName(), "%d", &it);
      WriteInteger((Int_t)(it/4));
   }
   delete tokens;
   PrintFast(5,"] 0 d");
}

////////////////////////////////////////////////////////////////////////////////
/// Change the line width

void TPDF::SetLineWidth(Width_t linewidth)
{
   if (linewidth == fLineWidth) return;
   fLineWidth = linewidth;
   if (fLineWidth!=0) {
      WriteReal(fLineScale*fLineWidth);
      PrintFast(2," w");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set color index for markers.

void TPDF::SetMarkerColor( Color_t cindex )
{
   fMarkerColor = cindex;
}

////////////////////////////////////////////////////////////////////////////////
/// Set color index for text

void TPDF::SetTextColor( Color_t cindex )
{
   fTextColor = cindex;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw text
///
///  - xx: x position of the text
///  - yy: y position of the text
///  - chars: text to be drawn

void TPDF::Text(Double_t xx, Double_t yy, const char *chars)
{
   if (fTextSize <= 0) return;

   const Double_t kDEGRAD = TMath::Pi()/180.;
   char str[8];
   Double_t x = xx;
   Double_t y = yy;

   // Font and text size
   Int_t font = abs(fTextFont)/10;
   if (font > kNumberOfFonts || font < 1) font = 1;

   Double_t wh = (Double_t)gPad->XtoPixel(gPad->GetX2());
   Double_t hh = (Double_t)gPad->YtoPixel(gPad->GetY1());
   Float_t tsize, ftsize;
   if (wh < hh) {
      tsize = fTextSize*wh;
      Int_t sizeTTF = (Int_t)(tsize*kScale+0.5); // TTF size
      ftsize = (sizeTTF*fXsize*gPad->GetAbsWNDC())/wh;
   } else {
      tsize = fTextSize*hh;
      Int_t sizeTTF = (Int_t)(tsize*kScale+0.5); // TTF size
      ftsize = (sizeTTF*fYsize*gPad->GetAbsHNDC())/hh;
   }
   Double_t fontsize = 72*(ftsize)/2.54;
   if (fontsize <= 0) return;

   // Text color
   SetColor(Int_t(fTextColor));

   // Clipping
   PrintStr(" q");
   Double_t x1 = XtoPDF(gPad->GetX1());
   Double_t x2 = XtoPDF(gPad->GetX2());
   Double_t y1 = YtoPDF(gPad->GetY1());
   Double_t y2 = YtoPDF(gPad->GetY2());
   WriteReal(x1);
   WriteReal(y1);
   WriteReal(x2 - x1);
   WriteReal(y2 - y1);
   PrintStr(" re W n");

   // Start the text
   if (!fCompress) PrintStr("@");

   // Text alignment
   Float_t tsizex = gPad->AbsPixeltoX(Int_t(tsize))-gPad->AbsPixeltoX(0);
   Float_t tsizey = gPad->AbsPixeltoY(0)-gPad->AbsPixeltoY(Int_t(tsize));
   Int_t txalh = fTextAlign/10;
   if (txalh < 1) txalh = 1; else if (txalh > 3) txalh = 3;
   Int_t txalv = fTextAlign%10;
   if (txalv < 1) txalv = 1; else if (txalv > 3) txalv = 3;
   if (txalv == 3) {
      y -= 0.8*tsizey*TMath::Cos(kDEGRAD*fTextAngle);
      x += 0.8*tsizex*TMath::Sin(kDEGRAD*fTextAngle);
   } else if (txalv == 2) {
      y -= 0.4*tsizey*TMath::Cos(kDEGRAD*fTextAngle);
      x += 0.4*tsizex*TMath::Sin(kDEGRAD*fTextAngle);
   }

   if (txalh > 1) {
      TText t;
      UInt_t w=0, h;
      t.SetTextSize(fTextSize);
      t.SetTextFont(fTextFont);
      t.GetTextExtent(w, h, chars);
      Double_t twx = gPad->AbsPixeltoX(w)-gPad->AbsPixeltoX(0);
      Double_t twy = gPad->AbsPixeltoY(0)-gPad->AbsPixeltoY(w);
      if (txalh == 2) {
         x = x-(twx/2)*TMath::Cos(kDEGRAD*fTextAngle);
         y = y-(twy/2)*TMath::Sin(kDEGRAD*fTextAngle);
      }
      if (txalh == 3) {
         x = x-twx*TMath::Cos(kDEGRAD*fTextAngle);
         y = y-twy*TMath::Sin(kDEGRAD*fTextAngle);
      }
   }

   // Text angle
   if (fTextAngle == 0) {
      PrintStr(" 1 0 0 1");
      WriteReal(XtoPDF(x));
      WriteReal(YtoPDF(y));
   } else if (fTextAngle == 90) {
      PrintStr(" 0 1 -1 0");
      WriteReal(XtoPDF(x));
      WriteReal(YtoPDF(y));
   } else if (fTextAngle == 270) {
      PrintStr(" 0 -1 1 0");
      WriteReal(XtoPDF(x));
      WriteReal(YtoPDF(y));
   } else {
      WriteReal(TMath::Cos(kDEGRAD*fTextAngle));
      WriteReal(TMath::Sin(kDEGRAD*fTextAngle));
      WriteReal(-TMath::Sin(kDEGRAD*fTextAngle));
      WriteReal(TMath::Cos(kDEGRAD*fTextAngle));
      WriteReal(XtoPDF(x));
      WriteReal(YtoPDF(y));
   }
   PrintStr(" cm");

   // Symbol Italic tan(15) = .26794
   if (font == 15) PrintStr(" 1 0 0.26794 1 0 0 cm");

   PrintStr(" BT");

   snprintf(str,8," /F%d",font);
   PrintStr(str);
   WriteReal(fontsize);
   PrintStr(" Tf");

   const Int_t len=strlen(chars);

   // Calculate the individual character placements.
   // Otherwise, if a string is printed in one line the kerning is not
   // performed. In order to measure the precise character positions we need to
   // trick FreeType into rendering high-resolution characters otherwise it will
   // stick to the screen pixel grid which is far worse than we can achieve on
   // print.
   const Float_t scale = 16.0;
   // Save current text attributes.
   TText saveAttText;
   saveAttText.TAttText::operator=(*this);
   TText t;
   t.SetTextSize(fTextSize * scale);
   t.SetTextFont(fTextFont);
   UInt_t wa1=0, wa0=0;
   t.GetTextAdvance(wa0, chars, kFALSE);
   t.GetTextAdvance(wa1, chars);
   t.TAttText::Modify();
   Bool_t kerning;
   if (wa0-wa1 != 0) kerning = kTRUE;
   else              kerning = kFALSE;
   Int_t *charDeltas = 0;
   if (kerning) {
        charDeltas = new Int_t[len];
        for (Int_t i = 0;i < len;i++) {
            UInt_t ww=0;
            t.GetTextAdvance(ww, chars + i);
            charDeltas[i] = wa1 - ww;
        }
        for (Int_t i = len - 1;i > 0;i--) {
            charDeltas[i] -= charDeltas[i-1];
        }
        char tmp[2];
        tmp[1] = 0;
        for (Int_t i = 1;i < len;i++) {
            tmp[0] = chars[i-1];
            UInt_t width=0;
            t.GetTextAdvance(width, &tmp[0], kFALSE);
            Double_t wwl = gPad->AbsPixeltoX(width - charDeltas[i]) - gPad->AbsPixeltoX(0);
            wwl -= 0.5*(gPad->AbsPixeltoX(1) - gPad->AbsPixeltoX(0)); // half a pixel ~ rounding error
            charDeltas[i] = (Int_t)((1000.0/Float_t(fontsize))*(XtoPDF(wwl) - XtoPDF(0))/scale);
        }
   }
   // Restore text attributes.
   saveAttText.TAttText::Modify();

   // Output the text. Escape some characters if needed
   if (kerning) PrintStr(" [");
   else         PrintStr(" (");

   for (Int_t i=0; i<len;i++) {
      if (chars[i]!='\n') {
         if (kerning) PrintStr("(");
         if (chars[i]=='(' || chars[i]==')') {
            snprintf(str,8,"\\%c",chars[i]);
         } else {
            snprintf(str,8,"%c",chars[i]);
         }
         PrintStr(str);
         if (kerning) {
            PrintStr(") ");
            if (i < len-1) {
               WriteInteger(charDeltas[i+1]);
            }
         }
      }
   }

   if (kerning) PrintStr("] TJ ET Q");
   else         PrintStr(") Tj ET Q");
   if (!fCompress) PrintStr("@");
   if (kerning) delete [] charDeltas;
}

////////////////////////////////////////////////////////////////////////////////
/// Write a string of characters
///
/// This method writes the string chars into a PDF file
/// at position xx,yy in world coordinates.

void TPDF::Text(Double_t, Double_t, const wchar_t *)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Write a string of characters in NDC

void TPDF::TextNDC(Double_t u, Double_t v, const char *chars)
{
   Double_t x = gPad->GetX1() + u*(gPad->GetX2() - gPad->GetX1());
   Double_t y = gPad->GetY1() + v*(gPad->GetY2() - gPad->GetY1());
   Text(x, y, chars);
}

////////////////////////////////////////////////////////////////////////////////
/// Write a string of characters in NDC

void TPDF::TextNDC(Double_t u, Double_t v, const wchar_t *chars)
{
   Double_t x = gPad->GetX1() + u*(gPad->GetX2() - gPad->GetX1());
   Double_t y = gPad->GetY1() + v*(gPad->GetY2() - gPad->GetY1());
   Text(x, y, chars);
}

////////////////////////////////////////////////////////////////////////////////
/// Convert U from NDC coordinate to PDF

Double_t TPDF::UtoPDF(Double_t u)
{
   Double_t cm = fXsize*(gPad->GetAbsXlowNDC() + u*gPad->GetAbsWNDC());
   return 72*cm/2.54;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert V from NDC coordinate to PDF

Double_t TPDF::VtoPDF(Double_t v)
{
   Double_t cm = fYsize*(gPad->GetAbsYlowNDC() + v*gPad->GetAbsHNDC());
   return 72*cm/2.54;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert X from world coordinate to PDF

Double_t TPDF::XtoPDF(Double_t x)
{
   Double_t u = (x - gPad->GetX1())/(gPad->GetX2() - gPad->GetX1());
   return  UtoPDF(u);
}

////////////////////////////////////////////////////////////////////////////////
/// Convert Y from world coordinate to PDF

Double_t TPDF::YtoPDF(Double_t y)
{
   Double_t v = (y - gPad->GetY1())/(gPad->GetY2() - gPad->GetY1());
   return  VtoPDF(v);
}

////////////////////////////////////////////////////////////////////////////////
/// Write the buffer in a compressed way

void TPDF::WriteCompressedBuffer()
{
   z_stream stream;
   int err;
   char *out = new char[2*fLenBuffer];

   stream.next_in   = (Bytef*)fBuffer;
   stream.avail_in  = (uInt)fLenBuffer;
   stream.next_out  = (Bytef*)out;
   stream.avail_out = (uInt)2*fLenBuffer;
   stream.zalloc    = (alloc_func)0;
   stream.zfree     = (free_func)0;
   stream.opaque    = (voidpf)0;

   err = deflateInit(&stream, Z_DEFAULT_COMPRESSION);
   if (err != Z_OK) {
      Error("WriteCompressedBuffer", "error in deflateInit (zlib)");
      delete [] out;
      return;
   }

   err = deflate(&stream, Z_FINISH);
   if (err != Z_STREAM_END) {
      deflateEnd(&stream);
      Error("WriteCompressedBuffer", "error in deflate (zlib)");
      delete [] out;
      return;
   }

   err = deflateEnd(&stream);

   fStream->write(out, stream.total_out);

   fNByte += stream.total_out;
   fStream->write("\n",1); fNByte++;
   fLenBuffer = 0;
   delete [] out;
   fCompress = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Write a Real number to the file.
/// This method overwrites TVirtualPS::WriteReal. Some PDF reader like
/// Acrobat do not work when a PDF file contains reals with exponent. This
/// method writes the real number "z" using the format "%f" instead of the
/// format "%g" when writing it with "%g" generates a number with exponent.

void TPDF::WriteReal(Float_t z, Bool_t space)
{
   char str[15];
   if (space) {
      snprintf(str,15," %g", z);
      if (strstr(str,"e") || strstr(str,"E")) snprintf(str,15," %10.8f", z);
   } else {
      snprintf(str,15,"%g", z);
      if (strstr(str,"e") || strstr(str,"E")) snprintf(str,15,"%10.8f", z);
   }
   PrintStr(str);
}

////////////////////////////////////////////////////////////////////////////////
/// Patterns encoding

void TPDF::PatternEncode()
{
   Int_t patternNb = kObjPattern;

   NewObject(kObjColorSpace);
   if (gStyle->GetColorModelPS()) {
      PrintStr("[/Pattern /DeviceCMYK]@");
   } else {
      PrintStr("[/Pattern /DeviceRGB]@");
   }
   PrintStr("endobj@");
   NewObject(kObjPatternResourses);
   PrintStr("<</ProcSet[/PDF]>>@");
   PrintStr("endobj@");

   NewObject(kObjPatternList);
   PrintStr("<<@");
   PrintStr(" /P01");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P02");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P03");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P04");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P05");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P06");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P07");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P08");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P09");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P10");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P11");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P12");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P13");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P14");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P15");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P16");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P17");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P18");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P19");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P20");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P21");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P22");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P23");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P24");
   WriteInteger(patternNb++);
   PrintStr(" 0 R");
   PrintStr(" /P25");
   WriteInteger(patternNb++);
   PrintStr(" 0 R@");
   PrintStr(">>@");
   PrintStr("endobj@");

   patternNb = kObjPattern;

   // P01
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[1 0 0 1 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 98 4]/XStep 98/YStep 4/Length 91/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211*\3442T\310T\3402P0P04\200\340\242T\256p\205<\240\220\027P0K\301P\241\034(\254\340\253\020m\250\020k\240\220\302e\244`\242\220\313ei\t\244r\200\272\215A\034\v \225\003\2241\202\310\030\201e\f!2\206@N0W \027@\200\001\0|c\024\357\n", 93);
   fNByte += 93;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P02
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.75 0 0 0.75 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 96 4]/XStep 96/YStep 4/Length 92/Filter/FlateDecode>>@");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211$\2121\n\2000\024C\367\234\"G\370\277\025\321+\b\016\342\340P\334tP\252\240\213\3277\332!\204\274\227\v\316\2150\032\335J\356\025\023O\241Np\247\363\021f\317\344\214\234\215\v\002+\036h\033U\326/~\243Ve\231PL\370\215\027\343\032#\006\274\002\f\0\242`\025:\n", 94);
   fNByte += 94;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P03
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.5 0 0 0.5 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 96 16]/XStep 96/YStep 16/Length 93/Filter/FlateDecode>>@");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211$\2121\n\2000\024C\367\234\"G\370\261(\366\n\202\20388\210\233\016J\025t\361\372\376\332!\204\274\227\033\342N\030\215\262\222g\303\304\313Q\347\360\240\370:f\317Y\f\\\214+**\360Dls'\177\306\274\032\257\344\256.\252\376\215\212\221\217\021\003>\001\006\0\317\243\025\254\n", 95);
   fNByte += 95;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P04
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.06 0 0 0.06 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 100 100]/XStep 100/YStep 100/Length 63/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211*\3442T\310T\3402P0P04\200\340\242T\256p\205<\240\220\027P0K\301D\241\034(\254\340\253\020\035k\240\220\002V\231\313\005S\233\303\025\314\025\310\005\020`\0\344\270\r\274\n", 65);
   fNByte += 65;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P05
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.06 0 0 0.06 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 100 100]/XStep 100/YStep 100/Length 66/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211*\3442T\310T\3402P0P04\200\340\242T\256p\205<\240\220\027P0K\301D\241\034(\254\340\253\020\035k\240\220\302\005Q\223\313\005\"\r\024r\270\202\271\002\271\0\002\f\0\344\320\r\274\n", 68);
   fNByte += 68;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P06
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.03 0 0 0.03 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 100 100]/XStep 100/YStep 100/Length 66/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211*\3442T\310T\3402P0P04\200\340\242T\256p\205<\240\220\027P0K\301D\241\034(\254\340\253\020\035k\240\220\302e\nR\232\v\242@js\270\202\271\002\271\0\002\f\0\345X\r\305\n", 68);
   fNByte += 68;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P07
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.03 0 0 0.03 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 100 100]/XStep 100/YStep 100/Length 68/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211*\3442T\310T\3402P0P04\200\340\242T\256p\205<\240\220\027P0K\301D\241\034(\254\340\253\020\035k\240\220\002\02465P\310\345\002)\0042r\270\202\271\002\271\0\002\f\0\345=\r\305\n", 70);
   fNByte += 70;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P08
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.06 0 0 0.06 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 101 101]/XStep 100/YStep 100/Length 139/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211D\217\261\016\3020\fDw\177\305\315L6Q\225|\003\022C\305\300Puk+\201\032$\272\360\373\330\265\323\016\271\330\367\234\344\"x\201\030\214\252\232\030+%\353VZ.jd\367\205\003x\241({]\311\324]\323|\342\006\033J\201:\306\325\230Jg\226J\261\275D\257#\337=\220\260\354k\233\351\211\217Z75\337\020\374\324\306\035\303\310\230\342x=\303\371\275\307o\332s\331\223\224\240G\330\a\365\364\027`\0\nX1}\n",141);
   fNByte += 141;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P09
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.06 0 0 0.06 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 100 100]/XStep 100/YStep 100/Length 108/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211*\3442T\310T\3402P0P04\200\340\242T\256p\205<\240\220\027P0K\301D\241\034(\254\340\253\020\035k\240\220\002\02465P\310\005RFFz&\020\002,d\240\220\314en\256g\0065\b,\001b\230\202$\240\232\214@\362\246`\2169H\336\024\2426\231\v&\200,\n\326\030\314\025\310\005\020`\0\f@\036\227\n", 110);
   fNByte += 110;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P10
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.06 0 0 0.06 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 100 100]/XStep 100/YStep 100/Length 93/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211*\3442T\310T\3402P0P04\200\340\242T\256p\205<\240\220\027P0K\301D\241\034(\254\340\253\020\035k\240\220\002\02465P\310\345\002)\0042r\200\332\r\241\\C \017dN.\027L\312\0\302\205\2535\205j6\205X\224\303\025\314\025\310\005\020`\0\2127\031\t\n", 95);
   fNByte += 95;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P11
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.125 0 0 0.125 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 100 100]/XStep 100/YStep 100/Length 164/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211\\\2171\016\3020\fEw\237\342\037\301ip\223^\001\211\001u`@l0\200(\022,\\\037;v\204\332\241\211\336\373\337V\363\246\204;\210\301H\354\337\347F'\274T\355U>\220\360U\215\003\316\027\306\2655\027=\a\306\223\304I\002m\332\330\356&\030\325\333fZ\275F\337\205\235\265O\270\032\004\331\214\336\305\270\004\227`\357i\256\223\342;]\344\255(!\372\356\205j\030\377K\335\220\344\377\210\274\306\022\330\337T{\214,\212;\301\3508\006\346\206\021O=\216|\212|\246#\375\004\030\0\216FF\207\n", 166);
   fNByte += 166;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P12
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.125 0 0 0.125 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 100 100]/XStep 100/YStep 100/Length 226/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211<P;n\3030\f\335y\n\236 \220DK\242\256P\240C\321\241C\221\311\311\220\242\016\220.\275~D\221/\203I\342}\370(?(\363\215)q\342\234\374\373\273\322\027\337'\3646\301\037\316\374?a~\347\357s\342\313\2045\361A9\237\322fc\231\200\236F\263\301\334;\211\017\207\rN\311\252S\\\227{\247\006w\207\244\303\255p+(\205\333\360e/v\356a\315\317\360\272\320b|w\276\203o\340k\b\004\027\v$b\226\235,\242\254t(\024\nu\305Vm\313\021\375\327\272\257\227fuf\226ju\356\222x\030\024\313\261S\215\377\341\274,\203\254\253Z\\\262A\262\205eD\350\210\320\201\225\212\320\036\241\355\025\372JE,\2266\344\366\310U\344\016HFx>\351\203\236\002\f\0d}e\216\n", 228);
   fNByte += 228;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P13
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.06 0 0 0.06 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 100 100]/XStep 100/YStep 100/Length 69/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211*\3442T\310T\3402P0P04\200\340\242T\256p\205<\240\220\027P0K\301D\241\034(\254\340\253\020\035k\240\220\002V\231\313\005S\233\303\005\241!\" ~0W \027@\200\001\0\331\227\020\253\n", 71);
   fNByte += 71;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P14
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.15 0 0 0.15 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 100 100]/XStep 80/YStep 80/Length 114/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\2114\214=\n\2000\f\205\367\234\342\035!-\241\364\f\202\20388\210\233\016J+\350\342\365M\3723\224\327\367}I\036r8A\f\206\343\372\336\203\026\334\212\006\205\027\004\237b\214X7\306\256\33032\331\240~\022y[\315\026\206\222\372\330}\264\036\253\217\335\353\240\030\b%\223\245o=X\227\346\245\355K\341\345@\3613M\364\v0\0\207o\"\261\n", 116);
   fNByte += 116;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P15
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.102 0 0 0.102 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 60 60]/XStep 60/YStep 60/Length 218/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211<\2211\016\3020\fEw\237\302'@\211c\267w@b@\f\f\210\2510\200(\022,\\\037\347\307\256Z\325\221\375\337\377\225\363\241\312\017\246\302\205'\274\337;\235\371\355\215\275\267\236\\\371\307\265\360\201/\327\3027o\233\361J\262\233\247~\362g\336\211zur!A]{\035}\031S\343\006p\241\226dKI\v\326\202\265\3153\331)X)\335fE\205M\235\373\327\r*\374\026\252\022\216u\223\200\361I\211\177\031\022\001#``\342GI\211\004c\221gi\246\231\247\221\247\231\247\233$XM3\315<\215<\315<K\211e\036#\215a4\366\344\035lm\214Z\314b\211Xj\337K\\\201$\332\325\v\365\2659\204\362\242\274'\v\221\r\321\211\216\364\027`\0\212'_\215\n", 220);
   fNByte += 220;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P16
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.1 0 0 0.05 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 100 100]/XStep 100/YStep 100/Length 123/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211*\3442T\310T\3402P0P04\200\340\242T\256p\205<\240\220\027P0K\301D\241\034(\254\340\253\020\035k\240\220\302ej\240\0D\271 \332\314X\317B\301\330\002H\230\233*\030\231\202\310d.CC=#\020\v*\rV\235\214\254\v\210r@\264\261\031P\241\031H5D\253\021H\267\005\3104 \v\344\016\260\002\020\003lB0W \027@\200\001\0hU \305\n", 125);
   fNByte += 125;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P17
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.06 0 0 0.06 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 100 100]/XStep 100/YStep 100/Length 66/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211*\3442T\310T\3402P0P04\200\340\242T\256p\205<\240\220\027P0K\301D\241\034(\254\340\253\020md\242\020k\240\220\002V\234\313\005S\236\303\025\314\025\310\005\020`\0\r\351\016B\n", 68);
   fNByte += 68;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P18
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.06 0 0 0.06 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 100 100]/XStep 100/YStep 100/Length 69/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211*\3442T\310T\3402P0P04\200\340\242T\256p\205<\240\220\027P0K\301D\241\034(\254\340\253\020md\242\020k\240\220\302\005Q\226\313\005\"\r\024r\270\202\271\002\271\0\002\f\0\016\001\016B\n", 71);
   fNByte += 71;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P19
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.117 0 0 0.117 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 100 100]/XStep 100/YStep 100/Length 149/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211L\216;\016\302@\fD{\237bN\020\331+6a\257\200D\201((P\252@\001R\220\240\341\372\370\263\216(\326\266f\336\330\373&\301\003\304`\b\307\373\334\351\202\227J\a\025\237\020|U\306\021\327\231q\243\306\250\214\325\372T\006\336\367\032\262\326\205\3124\264b\243$\"n.\244=\314\250!\2139\033\327\022i=\323\317\2518\332T}\347.\202\346W\373\372j\315\221\344\266\213=\237\241\344\034\361\264!\236w\344\177\271o8\323\211~\002\f\0\366\3026\233\n", 151);
   fNByte += 151;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P20
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.05 0 0 0.1 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 100 100]/XStep 100/YStep 100/Length 122/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211<L;\016\2030\f\335}\212w\002\344$M\2323 1 \006\006\304\224vhU\220`\341\372<\aT\311\366\263\336o\023\207\017D\241pz\355\376\226\021+\251\226\344\027\017\034\244\321a\232\025/\211\n\316r\343ORh\262}\317\210\344\032o\310)\302\2233\245\252[m\274\332\313\277!$\332\371\371\210`N\242\267$\217\263\246\252W\257\245\006\351\345\024`\0o\347 \305\n", 124);
   fNByte += 124;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P21
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.125 0 0 0.125 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 101 101]/XStep 100/YStep 100/Length 117/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211D\2151\n\2000\fE\367\234\342\037!)\224\336Ap\020\a\aq\323A\251\202.^\337$-\025\022^\372\033^n\022\354 \006CX\274\237\215&\\\032u\032\036\020\274\032\243\307\2740V]\027\234\024\242\"\033\2642En\324\312\224bc\262\\\230\377\301\332WM\224\212(U\221\375\265\301\025\016?\350\317P\215\221\033\213o\244\201>\001\006\0\031I'f\n", 119);
   fNByte += 119;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P22
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.125 0 0 0.125 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 101 101]/XStep 100/YStep 100/Length 118/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211<\215=\n\204P\f\204\373\234b\216\220<\b\357\016\302\026ba!vZ(\273\v\332x}\223\274\237\"|\223a\230\271Hp\200\030\fa\211\273w\232\3617k0\363\204\3401\033\037,+c#\3170~\2244\304\327EV\243r\247\272oOcr\337\323]H\t\226\252\334\252r\255\362\257\213(\t\304\250\326\315T\267\032\275q\242\221^\001\006\0\272\367(&\n", 120);
   fNByte += 120;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P23
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.06 0 0 0.06 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 100 100]/XStep 100/YStep 100/Length 169/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211<\220\273\n\0021\020E\373\371\212[[M\326\331\354\344\027\004\v\261\260\020;\025\224D\320\306\337w\036\254p\363\230\223\341$\344M\005\017\020\203Q8\307\347F'\274\f\355\f>Q\3605\214=\316\005\v.\214kt\217\230;)\324\366\245Fa\213e\320v\212r\022X\006\211Fi\3242\250J\224\302\020\367h\212\254I\\\325R\225o\03143\346U\235@a\t[\202Za\tA\202E`\351~O\002\235`\351~S\202\306h.m\253\264)\232K\217t\310\017q\354\a\353\247\364\377C\356\033\372\t0\0\bm:\375\n", 171);
   fNByte += 171;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P24
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.125 0 0 0.125 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 100 100]/XStep 100/YStep 100/Length 280/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\211DQ9N\004A\f\314\373\025\216\211\326\343v\037_@\"@\004\004\210\f\220@\003\022$|\177\335\345j\220v\345\251\303\343*\215\312\273\024\275\\d\375?\361dM\3162\306\337\214\337Y\336n\240m\217\036\301y\343\\<,i\250\0038F\035)\347l\322\026o\377\023\353|[\254\177\343\005;\315\317ky\224\257\240n\203\374\020\225\337\240\345N\236T\272<_\344\245\304^\3238\030\tc\236E\233xO\034\363\204>\251\317\324\233\023{\352\235\376\336S\357Fl\251\017\372\207\247>xoh&_\366Ud\331\253\314D\023\332\241\211\016\205\246\235\326\236*\275\307\204z8!s\031\335\306\\\306C\306\\\225\376\312\\\225\307\252\246\356\364\273Q\347\271:\371\341l\177\311e\210\3571\211\251#\374\302H\037:\342c\241\323\2617\320 \034\250\0\302\323a{\005%\302a\373(Zx\313\026\213@\215p\324}\026=\274e\217E8s\326}\026M\036\312}\271\n0\0\215\263\207\016\n", 282);
   fNByte += 282;
   PrintStr("endstream@");
   PrintStr("endobj@");

   // P25
   NewObject(patternNb++);
   PrintStr("<</Type/Pattern/Matrix[0.125 0 0 0.125 20 28]/PatternType 1/Resources");
   WriteInteger(kObjPatternResourses);
   PrintStr(" 0 R/PaintType 2/TilingType 1/BBox[0 0 101 101]/XStep 100/YStep 100/Length 54/Filter/FlateDecode>>");
   PrintStr("@");
   fStream->write("stream",6); fNByte += 6;
   fStream->write("\r\nH\2112T\310T\3402P0P\310\34526P\0\242\034.s\004m\016\242\r\r\f\024@\030\302\002\321iZP\305`M\346\310\212\201R\0\001\006\0\206\322\017\200\n", 56);
   fNByte += 56;
   PrintStr("endstream@");
   PrintStr("endobj@");
}
