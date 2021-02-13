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
#include "TDatime.h"
#include "TColor.h"
#include "TVirtualPad.h"
#include "TPoints.h"
#include "TSVG.h"
#include "TStyle.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "snprintf.h"

Int_t TSVG::fgLineJoin = 0;
Int_t TSVG::fgLineCap  = 0;

ClassImp(TSVG);

/** \class TSVG
\ingroup PS

\brief Interface to SVG

[SVG](http://www.w3.org/Graphics/SVG/Overview.htm8)
(Scalable Vector Graphics) is a language for describing
two-dimensional graphics in XML. SVG allows high quality vector graphics in
HTML pages.

To print a ROOT canvas "c1" into an SVG file simply do:
~~~ {.cpp}
     c1->Print("c1.svg");
~~~
The result is the ASCII file `c1.svg`.

It can be open directly using a web browser or included in a html document
the following way:
~~~ {.cpp}
<embed width="95%" height="500" src="c1.svg">
~~~
It is best viewed with Internet Explorer and you need the
[Adobe SVG Viewer](http://www.adobe.com/svg/viewer/install/main.html)

To zoom using the Adobe SVG Viewer, position the mouse over
the area you want to zoom and click the right button.

To define the zoom area,
use Control+drag to mark the boundaries of the zoom area.

To pan, use Alt+drag.
By clicking with the right mouse button on the SVG graphics you will get
a pop-up menu giving other ways to interact with the image.

SVG files can be used directly in compressed mode to minimize the time
transfer over the network. Compressed SVG files should be created using
`gzip` on a normal ASCII SVG file and should then be renamed
using the file extension `.svgz`.
*/

////////////////////////////////////////////////////////////////////////////////
/// Default SVG constructor

TSVG::TSVG() : TVirtualPS()
{
   fStream      = 0;
   fType        = 0;
   gVirtualPS   = this;
   fBoundingBox = kFALSE;
   fRange       = kFALSE;
   fXsize       = 0.;
   fYsize       = 0.;
   fYsizeSVG    = 0;
   SetTitle("SVG");
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the SVG interface
///
///  - fname : SVG file name
///  - wtype : SVG workstation type. Not used in the SVG driver. But as TSVG
///            inherits from TVirtualPS it should be kept. Anyway it is not
///            necessary to specify this parameter at creation time because it
///            has a default value (which is ignore in the SVG case).

TSVG::TSVG(const char *fname, Int_t wtype) : TVirtualPS(fname, wtype)
{
   fStream = 0;
   SetTitle("SVG");
   Open(fname, wtype);
}

////////////////////////////////////////////////////////////////////////////////
/// Open a SVG file

void TSVG::Open(const char *fname, Int_t wtype)
{
   if (fStream) {
      Warning("Open", "SVG file already open");
      return;
   }

   fLenBuffer = 0;
   fType      = abs(wtype);
   SetLineJoin(gStyle->GetJoinLinePS());
   SetLineCap(gStyle->GetCapLinePS());
   SetLineScale(gStyle->GetLineScalePS());
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
      printf("ERROR in TSVG::Open: Cannot open file:%s\n",fname);
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
/// Default SVG destructor

TSVG::~TSVG()
{
   Close();
}

////////////////////////////////////////////////////////////////////////////////
/// Close a SVG file

void TSVG::Close(Option_t *)
{
   if (!gVirtualPS) return;
   if (!fStream) return;
   if (gPad) gPad->Update();
   PrintStr("</svg>@");

   // Close file stream
   if (fStream) { fStream->close(); delete fStream; fStream = 0;}

   gVirtualPS = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Activate an already open SVG file

void TSVG::On()
{
   // fType is used to know if the SVG file is open. Unlike TPostScript, TSVG
   // has no "workstation type". In fact there is only one SVG type.

   if (!fType) {
      Error("On", "no SVG file open");
      Off();
      return;
   }
   gVirtualPS = this;
}

////////////////////////////////////////////////////////////////////////////////
/// Deactivate an already open SVG file

void TSVG::Off()
{
   gVirtualPS = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a Box

void TSVG::DrawBox(Double_t x1, Double_t y1, Double_t x2, Double_t  y2)
{
   static Double_t x[4], y[4];
   Double_t ix1 = XtoSVG(TMath::Min(x1,x2));
   Double_t ix2 = XtoSVG(TMath::Max(x1,x2));
   Double_t iy1 = YtoSVG(TMath::Min(y1,y2));
   Double_t iy2 = YtoSVG(TMath::Max(y1,y2));
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
         PrintStr("@");
         PrintFast(9,"<rect x=\"");
         WriteReal(ix1, kFALSE);
         PrintFast(5,"\" y=\"");
         WriteReal(iy2, kFALSE);
         PrintFast(9,"\" width=\"");
         WriteReal(ix2-ix1, kFALSE);
         PrintFast(10,"\" height=\"");
         WriteReal(iy1-iy2, kFALSE);
         PrintFast(7,"\" fill=");
         SetColorAlpha(5);
         PrintFast(2,"/>");
      }
   }
   if (fillis == 1) {
      PrintStr("@");
      PrintFast(9,"<rect x=\"");
      WriteReal(ix1, kFALSE);
      PrintFast(5,"\" y=\"");
      WriteReal(iy2, kFALSE);
      PrintFast(9,"\" width=\"");
      WriteReal(ix2-ix1, kFALSE);
      PrintFast(10,"\" height=\"");
      WriteReal(iy1-iy2, kFALSE);
      PrintFast(7,"\" fill=");
      SetColorAlpha(fFillColor);
      PrintFast(2,"/>");
   }
   if (fillis == 0) {
      if (fLineWidth<=0) return;
      PrintStr("@");
      PrintFast(9,"<rect x=\"");
      WriteReal(ix1, kFALSE);
      PrintFast(5,"\" y=\"");
      WriteReal(iy2, kFALSE);
      PrintFast(9,"\" width=\"");
      WriteReal(ix2-ix1, kFALSE);
      PrintFast(10,"\" height=\"");
      WriteReal(iy1-iy2, kFALSE);
      PrintFast(21,"\" fill=\"none\" stroke=");
      SetColorAlpha(fLineColor);
      PrintFast(2,"/>");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a Frame around a box
///
///  - mode = -1  the box looks as it is behind the screen
///  - mode =  1  the box looks as it is in front of the screen
///  - border is the border size in already pre-computed SVG units dark is the
///    color for the dark part of the frame light is the color for the light
///    part of the frame

void TSVG::DrawFrame(Double_t xl, Double_t yl, Double_t xt, Double_t  yt,
                            Int_t mode, Int_t border, Int_t dark, Int_t light)
{
   static Double_t xps[7], yps[7];
   Int_t i;
   Double_t ixd0, iyd0, ixdi, iydi, ix, iy;
   Int_t idx, idy;

   //- Draw top&left part of the box

   xps[0] = XtoSVG(xl);          yps[0] = YtoSVG(yl);
   xps[1] = xps[0] + border;     yps[1] = yps[0] - border;
   xps[2] = xps[1];              yps[2] = YtoSVG(yt) + border;
   xps[3] = XtoSVG(xt) - border; yps[3] = yps[2];
   xps[4] = XtoSVG(xt);          yps[4] = YtoSVG(yt);
   xps[5] = xps[0];              yps[5] = yps[4];
   xps[6] = xps[0];              yps[6] = yps[0];

   ixd0 = xps[0];
   iyd0 = yps[0];
   PrintStr("@");
   PrintFast(10,"<path d=\"M");
   WriteReal(ixd0, kFALSE);
   PrintFast(1,",");
   WriteReal(iyd0, kFALSE);

   idx = 0;
   idy = 0;
   for (i=1; i<7; i++) {
      ixdi = xps[i];
      iydi = yps[i];
      ix   = ixdi - ixd0;
      iy   = iydi - iyd0;
      ixd0 = ixdi;
      iyd0 = iydi;
      if( ix && iy) {
         if( idx ) { MovePS(idx,0); idx = 0; }
         if( idy ) { MovePS(0,idy); idy = 0; }
         MovePS(ix,iy);
         continue;
      }
      if ( ix ) {
         if( idy )  { MovePS(0,idy); idy = 0; }
         if( !idx ) { idx = ix; continue;}
         if( ix*idx > 0 ) {
            idx += ix;
         } else {
            MovePS(idx,0);
            idx  = ix;
         }
         continue;
      }
      if( iy ) {
         if( idx ) { MovePS(idx,0); idx = 0; }
         if( !idy) { idy = iy; continue;}
         if( iy*idy > 0 ) {
            idy += iy;
         } else {
            MovePS(0,idy);
            idy  = iy;
         }
      }
   }
   if( idx ) MovePS(idx,0);
   if( idy ) MovePS(0,idy);
   PrintFast(8,"z\" fill=");
   if (mode == -1) {
      SetColorAlpha(dark);
   } else {
      SetColorAlpha(light);
   }
   if (fgLineJoin)
      PrintStr(Form(" stroke-linejoin=\"%s\"", fgLineJoin == 1 ? "round" : "bevel"));
   if (fgLineCap)
      PrintStr(Form(" stroke-linecap=\"%s\"", fgLineCap == 1 ? "round" : "square"));
   PrintFast(2,"/>");

   //- Draw bottom&right part of the box
   xps[0] = XtoSVG(xl);          yps[0] = YtoSVG(yl);
   xps[1] = xps[0] + border;     yps[1] = yps[0] - border;
   xps[2] = XtoSVG(xt) - border; yps[2] = yps[1];
   xps[3] = xps[2];              yps[3] = YtoSVG(yt) + border;
   xps[4] = XtoSVG(xt);          yps[4] = YtoSVG(yt);
   xps[5] = xps[4];              yps[5] = yps[0];
   xps[6] = xps[0];              yps[6] = yps[0];

   ixd0 = xps[0];
   iyd0 = yps[0];
   PrintStr("@");
   PrintFast(10,"<path d=\"M");
   WriteReal(ixd0, kFALSE);
   PrintFast(1,",");
   WriteReal(iyd0, kFALSE);

   idx = 0;
   idy = 0;
   for (i=1;i<7;i++) {
      ixdi = xps[i];
      iydi = yps[i];
      ix   = ixdi - ixd0;
      iy   = iydi - iyd0;
      ixd0 = ixdi;
      iyd0 = iydi;
      if( ix && iy) {
         if( idx ) { MovePS(idx,0); idx = 0; }
         if( idy ) { MovePS(0,idy); idy = 0; }
         MovePS(ix,iy);
         continue;
      }
      if ( ix ) {
         if( idy )  { MovePS(0,idy); idy = 0; }
         if( !idx ) { idx = ix; continue;}
         if( ix*idx > 0 ) {
            idx += ix;
         } else {
            MovePS(idx,0);
            idx  = ix;
         }
         continue;
      }
      if( iy ) {
         if( idx ) { MovePS(idx,0); idx = 0; }
         if( !idy) { idy = iy; continue;}
         if( iy*idy > 0 ) {
            idy += iy;
         } else {
            MovePS(0,idy);
            idy  = iy;
         }
      }
   }
   if( idx ) MovePS(idx,0);
   if( idy ) MovePS(0,idy);
   PrintFast(8,"z\" fill=");
   if (mode == -1) {
      SetColorAlpha(light);
   } else {
      SetColorAlpha(dark);
   }
   if (fgLineJoin)
      PrintStr(Form(" stroke-linejoin=\"%s\"", fgLineJoin == 1 ? "round" : "bevel"));
   if (fgLineCap)
      PrintStr(Form(" stroke-linecap=\"%s\"", fgLineCap == 1 ? "round" : "square"));
   PrintFast(2,"/>");
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a PolyLine
///
///  Draw a polyline through  the points  xy.
///  - If NN=1 moves only to point x,y.
///  - If NN=0 the x,y are  written  in the SVG        file
///       according to the current transformation.
///  - If NN>0 the line is clipped as a line.
///  - If NN<0 the line is clipped as a fill area.

void TSVG::DrawPolyLine(Int_t nn, TPoints *xy)
{
   Int_t  n, idx, idy;
   Double_t ixd0, iyd0, ixdi, iydi, ix, iy;

   if (nn > 0) {
      n = nn;
   } else {
      n = -nn;
   }

   ixd0 = XtoSVG(xy[0].GetX());
   iyd0 = YtoSVG(xy[0].GetY());
   if( n <= 1) return;

   PrintFast(2," m");
   idx = 0;
   idy = 0;
   for (Int_t i=1;i<n;i++) {
      ixdi = XtoSVG(xy[i].GetX());
      iydi = YtoSVG(xy[i].GetY());
      ix   = ixdi - ixd0;
      iy   = iydi - iyd0;
      ixd0 = ixdi;
      iyd0 = iydi;
      if( ix && iy) {
         if( idx ) { MovePS(idx,0); idx = 0; }
         if( idy ) { MovePS(0,idy); idy = 0; }
         MovePS(ix,iy);
         continue;
      }
      if ( ix ) {
         if( idy )  { MovePS(0,idy); idy = 0; }
         if( !idx ) { idx = ix; continue;}
         if( ix*idx > 0 ) {
            idx += ix;
         } else {
            MovePS(idx,0);
            idx  = ix;
         }
         continue;
      }
      if( iy ) {
         if( idx ) { MovePS(idx,0); idx = 0; }
         if( !idy) { idy = iy; continue;}
         if( iy*idy > 0 ) {
            idy += iy;
         } else {
            MovePS(0,idy);
            idy  = iy;
         }
      }
   }
   if( idx ) MovePS(idx,0);
   if( idy ) MovePS(0,idy);

   if (nn > 0 ) {
   } else {
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a PolyLine in NDC space
///
///  Draw a polyline through  the points  xy.
///  --If NN=1 moves only to point x,y.
///  --If NN=0 the x,y are  written  in the SVG        file
///       according to the current transformation.
///  --If NN>0 the line is clipped as a line.
///  - If NN<0 the line is clipped as a fill area.

void TSVG::DrawPolyLineNDC(Int_t nn, TPoints *xy)
{
   Int_t  n, idx, idy;
   Double_t ixd0, iyd0, ixdi, iydi, ix, iy;

   if (nn > 0) {
      n = nn;
   } else {
      n = -nn;
   }

   ixd0 = UtoSVG(xy[0].GetX());
   iyd0 = VtoSVG(xy[0].GetY());
   if( n <= 1) return;

   idx = 0;
   idy = 0;
   for (Int_t i=1;i<n;i++) {
      ixdi = UtoSVG(xy[i].GetX());
      iydi = VtoSVG(xy[i].GetY());
      ix   = ixdi - ixd0;
      iy   = iydi - iyd0;
      ixd0 = ixdi;
      iyd0 = iydi;
      if( ix && iy) {
         if( idx ) { MovePS(idx,0); idx = 0; }
         if( idy ) { MovePS(0,idy); idy = 0; }
         MovePS(ix,iy);
         continue;
      }
      if ( ix ) {
         if( idy )  { MovePS(0,idy); idy = 0; }
         if( !idx ) { idx = ix; continue;}
         if( ix*idx > 0 ) {
            idx += ix;
         } else {
            MovePS(idx,0);
            idx  = ix;
         }
         continue;
      }
      if( iy ) {
         if( idx ) { MovePS(idx,0); idx = 0; }
         if( !idy) { idy = iy; continue;}
         if( iy*idy > 0 ) {
            idy += iy;
         } else {
            MovePS(0,idy);
            idy  = iy;
         }
      }
   }
   if( idx ) MovePS(idx,0);
   if( idy ) MovePS(0,idy);

   if (nn > 0 ) {
      if (xy[0].GetX() == xy[n-1].GetX() && xy[0].GetY() == xy[n-1].GetY()) PrintFast(3," cl");
   } else {
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint PolyMarker

void TSVG::DrawPolyMarker(Int_t n, Float_t *xw, Float_t *yw)
{
   fMarkerStyle = TMath::Abs(fMarkerStyle);
   Int_t ms = TAttMarker::GetMarkerStyleBase(fMarkerStyle);

   if (ms == 4)
      ms = 24;
   else if (ms >= 6 && ms <= 8)
      ms = 20;
   else if (ms >= 9 && ms <= 19)
      ms = 1;

   // Define the marker size
   Float_t msize  = fMarkerSize - TMath::Floor(TAttMarker::GetMarkerLineWidth(fMarkerStyle)/2.)/4.;
   if (fMarkerStyle == 1 || (fMarkerStyle >= 9 && fMarkerStyle <= 19)) msize = 0.01;
   if (fMarkerStyle == 6) msize = 0.02;
   if (fMarkerStyle == 7) msize = 0.04;

   const Int_t kBASEMARKER = 8;
   Float_t sbase = msize*kBASEMARKER;
   Float_t s2x = sbase / Float_t(gPad->GetWw() * gPad->GetAbsWNDC());
   msize = this->UtoSVG(s2x) - this->UtoSVG(0);

   Double_t m  = msize;
   Double_t m2 = m/2.;
   Double_t m3 = m/3.;
   Double_t m6 = m/6.;
   Double_t m8 = m/8.;
   Double_t m4 = m/4.;
   Double_t m0 = m/10.;

   // Draw the marker according to the type
   PrintStr("@");
   if ((ms > 19 && ms < 24) || ms == 29 || ms == 33 || ms == 34 ||
       ms == 39 || ms == 41 || ms == 43 || ms == 45 ||
       ms == 47 || ms == 48 || ms == 49) {
      PrintStr("<g fill=");
      SetColorAlpha(Int_t(fMarkerColor));
      PrintStr(">");
   } else {
      PrintStr("<g stroke=");
      SetColorAlpha(Int_t(fMarkerColor));
      PrintStr(" stroke-width=\"");
      WriteReal(TMath::Max(1, Int_t(TAttMarker::GetMarkerLineWidth(fMarkerStyle))), kFALSE);
      PrintStr("\" fill=\"none\"");
      if (fgLineJoin)
         PrintStr(Form(" stroke-linejoin=\"%s\"", fgLineJoin == 1 ? "round" : "bevel"));
      if (fgLineCap)
         PrintStr(Form(" stroke-linecap=\"%s\"", fgLineCap == 1 ? "round" : "square"));
      PrintStr(">");
   }
   Double_t ix,iy;
   for (Int_t i=0;i<n;i++) {
      ix = XtoSVG(xw[i]);
      iy = YtoSVG(yw[i]);
      PrintStr("@");
      // Dot (.)
      if (ms == 1) {
         PrintStr("<line x1=\"");
         WriteReal(ix-1, kFALSE);
         PrintStr("\" y1=\"");
         WriteReal(iy, kFALSE);
         PrintStr("\" x2=\"");
         WriteReal(ix, kFALSE);
         PrintStr("\" y2=\"");
         WriteReal(iy, kFALSE);
         PrintStr("\"/>");
      // Plus (+)
      } else if (ms == 2) {
         PrintStr("<line x1=\"");
         WriteReal(ix-m2, kFALSE);
         PrintStr("\" y1=\"");
         WriteReal(iy, kFALSE);
         PrintStr("\" x2=\"");
         WriteReal(ix+m2, kFALSE);
         PrintStr("\" y2=\"");
         WriteReal(iy, kFALSE);
         PrintStr("\"/>");

         PrintStr("<line x1=\"");
         WriteReal(ix, kFALSE);
         PrintStr("\" y1=\"");
         WriteReal(iy-m2, kFALSE);
         PrintStr("\" x2=\"");
         WriteReal(ix, kFALSE);
         PrintStr("\" y2=\"");
         WriteReal(iy+m2, kFALSE);
         PrintStr("\"/>");
      // X shape (X)
      } else if (ms == 5) {
         PrintStr("<line x1=\"");
         WriteReal(ix-m2*0.707, kFALSE);
         PrintStr("\" y1=\"");
         WriteReal(iy-m2*0.707, kFALSE);
         PrintStr("\" x2=\"");
         WriteReal(ix+m2*0.707, kFALSE);
         PrintStr("\" y2=\"");
         WriteReal(iy+m2*0.707, kFALSE);
         PrintStr("\"/>");

         PrintStr("<line x1=\"");
         WriteReal(ix-m2*0.707, kFALSE);
         PrintStr("\" y1=\"");
         WriteReal(iy+m2*0.707, kFALSE);
         PrintStr("\" x2=\"");
         WriteReal(ix+m2*0.707, kFALSE);
         PrintStr("\" y2=\"");
         WriteReal(iy-m2*0.707, kFALSE);
         PrintStr("\"/>");
      // Asterisk shape (*)
      } else if (ms == 3 || ms == 31) {
         PrintStr("<line x1=\"");
         WriteReal(ix-m2, kFALSE);
         PrintStr("\" y1=\"");
         WriteReal(iy, kFALSE);
         PrintStr("\" x2=\"");
         WriteReal(ix+m2, kFALSE);
         PrintStr("\" y2=\"");
         WriteReal(iy, kFALSE);
         PrintStr("\"/>");

         PrintStr("<line x1=\"");
         WriteReal(ix, kFALSE);
         PrintStr("\" y1=\"");
         WriteReal(iy-m2, kFALSE);
         PrintStr("\" x2=\"");
         WriteReal(ix, kFALSE);
         PrintStr("\" y2=\"");
         WriteReal(iy+m2, kFALSE);
         PrintStr("\"/>");

         PrintStr("<line x1=\"");
         WriteReal(ix-m2*0.707, kFALSE);
         PrintStr("\" y1=\"");
         WriteReal(iy-m2*0.707, kFALSE);
         PrintStr("\" x2=\"");
         WriteReal(ix+m2*0.707, kFALSE);
         PrintStr("\" y2=\"");
         WriteReal(iy+m2*0.707, kFALSE);
         PrintStr("\"/>");

         PrintStr("<line x1=\"");
         WriteReal(ix-m2*0.707, kFALSE);
         PrintStr("\" y1=\"");
         WriteReal(iy+m2*0.707, kFALSE);
         PrintStr("\" x2=\"");
         WriteReal(ix+m2*0.707, kFALSE);
         PrintStr("\" y2=\"");
         WriteReal(iy-m2*0.707, kFALSE);
         PrintStr("\"/>");
      // Circle
      } else if (ms == 24 || ms == 20) {
         PrintStr("<circle cx=\"");
         WriteReal(ix, kFALSE);
         PrintStr("\" cy=\"");
         WriteReal(iy, kFALSE);
         PrintStr("\" r=\"");
         if (m2<=0) m2=1;
         WriteReal(m2, kFALSE);
         PrintStr("\" fill=\"none\"");
         PrintStr("/>");
      // Square
      } else if (ms == 25 || ms == 21) {
         PrintStr("<rect x=\"");
         WriteReal(ix-m2, kFALSE);
         PrintStr("\" y=\"");
         WriteReal(iy-m2, kFALSE);
         PrintStr("\" width=\"");
         WriteReal(m, kFALSE);
         PrintStr("\" height=\"");
         WriteReal(m, kFALSE);
         PrintStr("\" fill=\"none\"");
         PrintStr("/>");
      // Down triangle
      } else if (ms == 26 || ms == 22) {
         PrintStr("<polygon points=\"");
         WriteReal(ix); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m2);
         PrintStr("\"/>");
      // Up triangle
      } else if (ms == 23 || ms == 32) {
         PrintStr("<polygon points=\"");
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix); PrintStr(","); WriteReal(iy+m2);
         PrintStr("\"/>");
      // Diamond
      } else if (ms == 27 || ms == 33) {
         PrintStr("<polygon points=\"");
         WriteReal(ix); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m3); PrintStr(","); WriteReal(iy);
         WriteReal(ix); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m3); PrintStr(","); WriteReal(iy);
         PrintStr("\"/>");
      // Cross
      } else if (ms == 28 || ms == 34) {
         PrintStr("<polygon points=\"");
         WriteReal(ix-m6); PrintStr(","); WriteReal(iy-m6);
         WriteReal(ix-m6); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m6); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m6); PrintStr(","); WriteReal(iy-m6);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m6);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m6);
         WriteReal(ix+m6); PrintStr(","); WriteReal(iy+m6);
         WriteReal(ix+m6); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m6); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m6); PrintStr(","); WriteReal(iy+m6);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m6);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m6);
         PrintStr("\"/>");
      } else if (ms == 29 || ms == 30) {
         PrintStr("<polygon points=\"");
         WriteReal(ix); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix+0.112255*m); PrintStr(","); WriteReal(iy+0.15451*m);
         WriteReal(ix+0.47552*m); PrintStr(","); WriteReal(iy+0.15451*m);
         WriteReal(ix+0.181635*m); PrintStr(","); WriteReal(iy-0.05902*m);
         WriteReal(ix+0.29389*m); PrintStr(","); WriteReal(iy-0.40451*m);
         WriteReal(ix); PrintStr(","); WriteReal(iy-0.19098*m);
         WriteReal(ix-0.29389*m); PrintStr(","); WriteReal(iy-0.40451*m);
         WriteReal(ix-0.181635*m); PrintStr(","); WriteReal(iy-0.05902*m);
         WriteReal(ix-0.47552*m); PrintStr(","); WriteReal(iy+0.15451*m);
         WriteReal(ix-0.112255*m); PrintStr(","); WriteReal(iy+0.15451*m);
         PrintStr("\"/>");
      } else if (ms == 35) {
         PrintStr("<polygon points=\"");
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy   );
         WriteReal(ix   ); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy   );
         WriteReal(ix   ); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy   );
         WriteReal(ix   ); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy-m2);
         PrintStr("\"/>");
      } else if (ms == 36) {
         PrintStr("<polygon points=\"");
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m2);
         PrintStr("\"/>");
      } else if (ms == 37 || ms == 39) {
         PrintStr("<polygon points=\"");
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         PrintStr("\"/>");
      } else if (ms == 38) {
         PrintStr("<polygon points=\"");
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy   );
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy   );
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         WriteReal(ix   ); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy);
         PrintStr("\"/>");
      } else if (ms == 40 || ms == 41) {
         PrintStr("<polygon points=\"");
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         PrintStr("\"/>");
      } else if (ms == 42 || ms == 43) {
         PrintStr("<polygon points=\"");
         WriteReal(ix   ); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m8); PrintStr(","); WriteReal(iy+m8);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy   );
         WriteReal(ix-m8); PrintStr(","); WriteReal(iy-m8);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m8); PrintStr(","); WriteReal(iy-m8);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m8); PrintStr(","); WriteReal(iy+m8);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy+m2);
         PrintStr("\"/>");
      } else if (ms == 44) {
         PrintStr("<polygon points=\"");
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         PrintStr("\"/>");
      } else if (ms == 45) {
         PrintStr("<polygon points=\"");
         WriteReal(ix+m0); PrintStr(","); WriteReal(iy+m0);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m0); PrintStr(","); WriteReal(iy+m0);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix-m0); PrintStr(","); WriteReal(iy-m0);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m0); PrintStr(","); WriteReal(iy-m0);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix+m0); PrintStr(","); WriteReal(iy+m0);
         PrintStr("\"/>");
      } else if (ms == 46 || ms == 47) {
         PrintStr("<polygon points=\"");
         WriteReal(ix   ); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy   );
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy+m4);
         PrintStr("\"/>");
      } else if (ms == 48) {
         PrintStr("<polygon points=\"");
         WriteReal(ix   ); PrintStr(","); WriteReal(iy+m4*1.01);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy   );
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy+m4*0.99);
         WriteReal(ix+m4*0.99); PrintStr(","); WriteReal(iy   );
         WriteReal(ix   ); PrintStr(","); WriteReal(iy-m4*0.99);
         WriteReal(ix-m4*0.99); PrintStr(","); WriteReal(iy   );
         WriteReal(ix   ); PrintStr(","); WriteReal(iy+m4*0.99);
         PrintStr("\"/>");
      } else if (ms == 49) {
         PrintStr("<polygon points=\"");
         WriteReal(ix-m6); PrintStr(","); WriteReal(iy-m6*1.01);
         WriteReal(ix-m6); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m6); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m6); PrintStr(","); WriteReal(iy-m6);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m6);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m6);
         WriteReal(ix+m6); PrintStr(","); WriteReal(iy+m6);
         WriteReal(ix+m6); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m6); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m6); PrintStr(","); WriteReal(iy+m6);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m6);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m6);
         WriteReal(ix-m6); PrintStr(","); WriteReal(iy-m6*0.99);
         WriteReal(ix-m6); PrintStr(","); WriteReal(iy+m6);
         WriteReal(ix+m6); PrintStr(","); WriteReal(iy+m6);
         WriteReal(ix+m6); PrintStr(","); WriteReal(iy-m6);
         PrintStr("\"/>");
      } else {
         PrintStr("<line x1=\"");
         WriteReal(ix-1, kFALSE);
         PrintStr("\" y1=\"");
         WriteReal(iy, kFALSE);
         PrintStr("\" x2=\"");
         WriteReal(ix, kFALSE);
         PrintStr("\" y2=\"");
         WriteReal(iy, kFALSE);
         PrintStr("\"/>");
      }
   }
   PrintStr("@");
   PrintStr("</g>");
}

////////////////////////////////////////////////////////////////////////////////
/// Paint PolyMarker

void TSVG::DrawPolyMarker(Int_t n, Double_t *xw, Double_t *yw)
{
   fMarkerStyle = TMath::Abs(fMarkerStyle);
   Int_t ms = TAttMarker::GetMarkerStyleBase(fMarkerStyle);

   if (ms == 4)
      ms = 24;
   else if (ms >= 6 && ms <= 8)
      ms = 20;
   else if (ms >= 9 && ms <= 19)
      ms = 1;

   // Define the marker size
   Float_t msize  = fMarkerSize - TMath::Floor(TAttMarker::GetMarkerLineWidth(fMarkerStyle)/2.)/4.;
   if (fMarkerStyle == 1 || (fMarkerStyle >= 9 && fMarkerStyle <= 19)) msize = 0.01;
   if (fMarkerStyle == 6) msize = 0.02;
   if (fMarkerStyle == 7) msize = 0.04;

   const Int_t kBASEMARKER = 8;
   Float_t sbase = msize*kBASEMARKER;
   Float_t s2x = sbase / Float_t(gPad->GetWw() * gPad->GetAbsWNDC());
   msize = this->UtoSVG(s2x) - this->UtoSVG(0);

   Double_t m  = msize;
   Double_t m2 = m/2;
   Double_t m3 = m/3;
   Double_t m6 = m/6;
   Double_t m4 = m/4.;
   Double_t m8 = m/8.;
   Double_t m0 = m/10.;

   // Draw the marker according to the type
   PrintStr("@");
   if ((ms > 19 && ms < 24) || ms == 29 || ms == 33 || ms == 34 ||
       ms == 39 || ms == 41 || ms == 43 || ms == 45 ||
       ms == 47 || ms == 48 || ms == 49) {
      PrintStr("<g fill=");
      SetColorAlpha(Int_t(fMarkerColor));
      PrintStr(">");
   } else {
      PrintStr("<g stroke=");
      SetColorAlpha(Int_t(fMarkerColor));
      PrintStr(" stroke-width=\"");
      WriteReal(TMath::Max(1, Int_t(TAttMarker::GetMarkerLineWidth(fMarkerStyle))), kFALSE);
      PrintStr("\" fill=\"none\"");
      if (fgLineJoin)
         PrintStr(Form(" stroke-linejoin=\"%s\"", fgLineJoin == 1 ? "round" : "bevel"));
      if (fgLineCap)
         PrintStr(Form(" stroke-linecap=\"%s\"", fgLineCap == 1 ? "round" : "square"));
      PrintStr(">");
   }
   Double_t ix,iy;
   for (Int_t i=0;i<n;i++) {
      ix = XtoSVG(xw[i]);
      iy = YtoSVG(yw[i]);
      PrintStr("@");
      // Dot (.)
      if (ms == 1) {
         PrintStr("<line x1=\"");
         WriteReal(ix-1, kFALSE);
         PrintStr("\" y1=\"");
         WriteReal(iy, kFALSE);
         PrintStr("\" x2=\"");
         WriteReal(ix, kFALSE);
         PrintStr("\" y2=\"");
         WriteReal(iy, kFALSE);
         PrintStr("\"/>");
      // Plus (+)
      } else if (ms == 2) {
         PrintStr("<line x1=\"");
         WriteReal(ix-m2, kFALSE);
         PrintStr("\" y1=\"");
         WriteReal(iy, kFALSE);
         PrintStr("\" x2=\"");
         WriteReal(ix+m2, kFALSE);
         PrintStr("\" y2=\"");
         WriteReal(iy, kFALSE);
         PrintStr("\"/>");

         PrintStr("<line x1=\"");
         WriteReal(ix, kFALSE);
         PrintStr("\" y1=\"");
         WriteReal(iy-m2, kFALSE);
         PrintStr("\" x2=\"");
         WriteReal(ix, kFALSE);
         PrintStr("\" y2=\"");
         WriteReal(iy+m2, kFALSE);
         PrintStr("\"/>");
      // X shape (X)
      } else if (ms == 5) {
         PrintStr("<line x1=\"");
         WriteReal(ix-m2*0.707, kFALSE);
         PrintStr("\" y1=\"");
         WriteReal(iy-m2*0.707, kFALSE);
         PrintStr("\" x2=\"");
         WriteReal(ix+m2*0.707, kFALSE);
         PrintStr("\" y2=\"");
         WriteReal(iy+m2*0.707, kFALSE);
         PrintStr("\"/>");

         PrintStr("<line x1=\"");
         WriteReal(ix-m2*0.707, kFALSE);
         PrintStr("\" y1=\"");
         WriteReal(iy+m2*0.707, kFALSE);
         PrintStr("\" x2=\"");
         WriteReal(ix+m2*0.707, kFALSE);
         PrintStr("\" y2=\"");
         WriteReal(iy-m2*0.707, kFALSE);
         PrintStr("\"/>");
      // Asterisk shape (*)
      } else if (ms == 3 || ms == 31) {
         PrintStr("<line x1=\"");
         WriteReal(ix-m2, kFALSE);
         PrintStr("\" y1=\"");
         WriteReal(iy, kFALSE);
         PrintStr("\" x2=\"");
         WriteReal(ix+m2, kFALSE);
         PrintStr("\" y2=\"");
         WriteReal(iy, kFALSE);
         PrintStr("\"/>");

         PrintStr("<line x1=\"");
         WriteReal(ix, kFALSE);
         PrintStr("\" y1=\"");
         WriteReal(iy-m2, kFALSE);
         PrintStr("\" x2=\"");
         WriteReal(ix, kFALSE);
         PrintStr("\" y2=\"");
         WriteReal(iy+m2, kFALSE);
         PrintStr("\"/>");

         PrintStr("<line x1=\"");
         WriteReal(ix-m2*0.707, kFALSE);
         PrintStr("\" y1=\"");
         WriteReal(iy-m2*0.707, kFALSE);
         PrintStr("\" x2=\"");
         WriteReal(ix+m2*0.707, kFALSE);
         PrintStr("\" y2=\"");
         WriteReal(iy+m2*0.707, kFALSE);
         PrintStr("\"/>");

         PrintStr("<line x1=\"");
         WriteReal(ix-m2*0.707, kFALSE);
         PrintStr("\" y1=\"");
         WriteReal(iy+m2*0.707, kFALSE);
         PrintStr("\" x2=\"");
         WriteReal(ix+m2*0.707, kFALSE);
         PrintStr("\" y2=\"");
         WriteReal(iy-m2*0.707, kFALSE);
         PrintStr("\"/>");
      // Circle
      } else if (ms == 24 || ms == 20) {
         PrintStr("<circle cx=\"");
         WriteReal(ix, kFALSE);
         PrintStr("\" cy=\"");
         WriteReal(iy, kFALSE);
         PrintStr("\" r=\"");
         if (m2<=0) m2=1;
         WriteReal(m2, kFALSE);
         PrintStr("\"/>");
      // Square
      } else if (ms == 25 || ms == 21) {
         PrintStr("<rect x=\"");
         WriteReal(ix-m2, kFALSE);
         PrintStr("\" y=\"");
         WriteReal(iy-m2, kFALSE);
         PrintStr("\" width=\"");
         WriteReal(m, kFALSE);
         PrintStr("\" height=\"");
         WriteReal(m, kFALSE);
         PrintStr("\"/>");
      // Down triangle
      } else if (ms == 26 || ms == 22) {
         PrintStr("<polygon points=\"");
         WriteReal(ix); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m2);
         PrintStr("\"/>");
      // Up triangle
      } else if (ms == 23 || ms == 32) {
         PrintStr("<polygon points=\"");
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix); PrintStr(","); WriteReal(iy+m2);
         PrintStr("\"/>");
      // Diamond
      } else if (ms == 27 || ms == 33) {
         PrintStr("<polygon points=\"");
         WriteReal(ix); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m3); PrintStr(","); WriteReal(iy);
         WriteReal(ix); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m3); PrintStr(","); WriteReal(iy);
         PrintStr("\"/>");
      // Cross
      } else if (ms == 28 || ms == 34) {
         PrintStr("<polygon points=\"");
         WriteReal(ix-m6); PrintStr(","); WriteReal(iy-m6);
         WriteReal(ix-m6); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m6); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m6); PrintStr(","); WriteReal(iy-m6);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m6);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m6);
         WriteReal(ix+m6); PrintStr(","); WriteReal(iy+m6);
         WriteReal(ix+m6); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m6); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m6); PrintStr(","); WriteReal(iy+m6);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m6);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m6);
         PrintStr("\"/>");
      } else if (ms == 29 || ms == 30) {
         PrintStr("<polygon points=\"");
         WriteReal(ix); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix+0.112255*m); PrintStr(","); WriteReal(iy+0.15451*m);
         WriteReal(ix+0.47552*m); PrintStr(","); WriteReal(iy+0.15451*m);
         WriteReal(ix+0.181635*m); PrintStr(","); WriteReal(iy-0.05902*m);
         WriteReal(ix+0.29389*m); PrintStr(","); WriteReal(iy-0.40451*m);
         WriteReal(ix); PrintStr(","); WriteReal(iy-0.19098*m);
         WriteReal(ix-0.29389*m); PrintStr(","); WriteReal(iy-0.40451*m);
         WriteReal(ix-0.181635*m); PrintStr(","); WriteReal(iy-0.05902*m);
         WriteReal(ix-0.47552*m); PrintStr(","); WriteReal(iy+0.15451*m);
         WriteReal(ix-0.112255*m); PrintStr(","); WriteReal(iy+0.15451*m);
         PrintStr("\"/>");
      } else if (ms == 35) {
         PrintStr("<polygon points=\"");
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy   );
         WriteReal(ix   ); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy   );
         WriteReal(ix   ); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy   );
         WriteReal(ix   ); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy-m2);
         PrintStr("\"/>");
      } else if (ms == 36) {
         PrintStr("<polygon points=\"");
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m2);
         PrintStr("\"/>");
      } else if (ms == 37 || ms == 39) {
         PrintStr("<polygon points=\"");
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         PrintStr("\"/>");
      } else if (ms == 38) {
         PrintStr("<polygon points=\"");
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy   );
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy   );
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         WriteReal(ix   ); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy);
         PrintStr("\"/>");
      } else if (ms == 40 || ms == 41) {
         PrintStr("<polygon points=\"");
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         PrintStr("\"/>");
      } else if (ms == 42 || ms == 43) {
         PrintStr("<polygon points=\"");
         WriteReal(ix   ); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m8); PrintStr(","); WriteReal(iy+m8);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy   );
         WriteReal(ix-m8); PrintStr(","); WriteReal(iy-m8);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m8); PrintStr(","); WriteReal(iy-m8);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m8); PrintStr(","); WriteReal(iy+m8);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy+m2);
         PrintStr("\"/>");
      } else if (ms == 44) {
         PrintStr("<polygon points=\"");
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy   );
         PrintStr("\"/>");
      } else if (ms == 45) {
         PrintStr("<polygon points=\"");
         WriteReal(ix+m0); PrintStr(","); WriteReal(iy+m0);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m0); PrintStr(","); WriteReal(iy+m0);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix-m0); PrintStr(","); WriteReal(iy-m0);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m0); PrintStr(","); WriteReal(iy-m0);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix+m0); PrintStr(","); WriteReal(iy+m0);
         PrintStr("\"/>");
      } else if (ms == 46 || ms == 47) {
         PrintStr("<polygon points=\"");
         WriteReal(ix   ); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy   );
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy+m4);
         PrintStr("\"/>");
      } else if (ms == 48) {
         PrintStr("<polygon points=\"");
         WriteReal(ix   ); PrintStr(","); WriteReal(iy+m4*1.01);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy   );
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix-m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m4);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy   );
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m4);
         WriteReal(ix+m4); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix   ); PrintStr(","); WriteReal(iy+m4*0.99);
         WriteReal(ix+m4*0.99); PrintStr(","); WriteReal(iy   );
         WriteReal(ix   ); PrintStr(","); WriteReal(iy-m4*0.99);
         WriteReal(ix-m4*0.99); PrintStr(","); WriteReal(iy   );
         WriteReal(ix   ); PrintStr(","); WriteReal(iy+m4*0.99);
         PrintStr("\"/>");
      } else if (ms == 49) {
         PrintStr("<polygon points=\"");
         WriteReal(ix-m6); PrintStr(","); WriteReal(iy-m6*1.01);
         WriteReal(ix-m6); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m6); PrintStr(","); WriteReal(iy-m2);
         WriteReal(ix+m6); PrintStr(","); WriteReal(iy-m6);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy-m6);
         WriteReal(ix+m2); PrintStr(","); WriteReal(iy+m6);
         WriteReal(ix+m6); PrintStr(","); WriteReal(iy+m6);
         WriteReal(ix+m6); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m6); PrintStr(","); WriteReal(iy+m2);
         WriteReal(ix-m6); PrintStr(","); WriteReal(iy+m6);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy+m6);
         WriteReal(ix-m2); PrintStr(","); WriteReal(iy-m6);
         WriteReal(ix-m6); PrintStr(","); WriteReal(iy-m6*0.99);
         WriteReal(ix-m6); PrintStr(","); WriteReal(iy+m6);
         WriteReal(ix+m6); PrintStr(","); WriteReal(iy+m6);
         WriteReal(ix+m6); PrintStr(","); WriteReal(iy-m6);
         PrintStr("\"/>");
      } else {
         PrintStr("<line x1=\"");
         WriteReal(ix-1, kFALSE);
         PrintStr("\" y1=\"");
         WriteReal(iy, kFALSE);
         PrintStr("\" x2=\"");
         WriteReal(ix, kFALSE);
         PrintStr("\" y2=\"");
         WriteReal(iy, kFALSE);
         PrintStr("\"/>");
      }
   }
   PrintStr("@");
   PrintStr("</g>");
}

////////////////////////////////////////////////////////////////////////////////
/// This function defines a path with xw and yw and draw it according the
/// value of nn:
///
///  - If nn>0 a line is drawn.
///  - If nn<0 a closed polygon is drawn.

void TSVG::DrawPS(Int_t nn, Double_t *xw, Double_t *yw)
{
   Int_t  n, fais, fasi;
   Double_t ixd0, iyd0, idx, idy, ixdi, iydi, ix, iy;
   fais = fasi = 0;

   if (nn > 0) {
      if (fLineWidth<=0) return;
      n = nn;
   } else {
      n = -nn;
      fais = fFillStyle/1000;
      fasi = fFillStyle%1000;
      if (fais == 3 || fais == 2) {
         if (fasi > 100 && fasi <125) {
            return;
         }
         if (fasi > 0 && fasi < 26) {
         }
      }
   }

   if( n <= 1) {
      Error("DrawPS", "Two points are needed");
      return;
   }

   ixd0 = XtoSVG(xw[0]);
   iyd0 = YtoSVG(yw[0]);

   PrintStr("@");
   PrintFast(10,"<path d=\"M");
   WriteReal(ixd0, kFALSE);
   PrintFast(1,",");
   WriteReal(iyd0, kFALSE);

   idx = idy = 0;
   for (Int_t i=1;i<n;i++) {
      ixdi = XtoSVG(xw[i]);
      iydi = YtoSVG(yw[i]);
      ix   = ixdi - ixd0;
      iy   = iydi - iyd0;
      ixd0 = ixdi;
      iyd0 = iydi;
      if( ix && iy) {
         if( idx ) { MovePS(idx,0); idx = 0; }
         if( idy ) { MovePS(0,idy); idy = 0; }
         MovePS(ix,iy);
      } else if ( ix ) {
         if( idy )  { MovePS(0,idy); idy = 0;}
         if( !idx ) { idx = ix;}
         else if( TMath::Sign(ix,idx) == ix )       idx += ix;
         else { MovePS(idx,0);  idx  = ix;}
      } else if( iy ) {
         if( idx ) { MovePS(idx,0); idx = 0;}
         if( !idy) { idy = iy;}
         else if( TMath::Sign(iy,idy) == iy)         idy += iy;
         else { MovePS(0,idy);    idy  = iy;}
      }
   }
   if (idx) MovePS(idx,0);
   if (idy) MovePS(0,idy);

   if (nn > 0 ) {
      if (xw[0] == xw[n-1] && yw[0] == yw[n-1]) PrintFast(1,"z");
      PrintFast(21,"\" fill=\"none\" stroke=");
      SetColorAlpha(fLineColor);
      if(fLineWidth > 1.) {
         PrintFast(15," stroke-width=\"");
         WriteReal(fLineWidth, kFALSE);
         PrintFast(1,"\"");
      }
      if (fLineStyle > 1) {
         PrintFast(19," stroke-dasharray=\"");
         TString st = (TString)gStyle->GetLineStyleString(fLineStyle);
         TObjArray *tokens = st.Tokenize(" ");
         for (Int_t j = 0; j<tokens->GetEntries(); j++) {
            Int_t it;
            sscanf(((TObjString*)tokens->At(j))->GetName(), "%d", &it);
            if (j>0) PrintFast(1,",");
            WriteReal(it/4);
         }
         delete tokens;
         PrintFast(1,"\"");
      }
   } else {
      PrintFast(8,"z\" fill=");
      if (fais == 0) {
         PrintFast(14,"\"none\" stroke=");
         SetColorAlpha(fFillColor);
      } else {
         SetColorAlpha(fFillColor);
      }
   }
   if (fgLineJoin)
      PrintStr(Form(" stroke-linejoin=\"%s\"", fgLineJoin == 1 ? "round" : "bevel"));
   if (fgLineCap)
      PrintStr(Form(" stroke-linecap=\"%s\"", fgLineCap == 1 ? "round" : "square"));
   PrintFast(2,"/>");
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the SVG file. The main task of the function is to output the
/// SVG header file which consist in <title>, <desc> and <defs>. The
/// HeaderPS provided by the user program is written in the <defs> part.

void TSVG::Initialize()
{
   // Title
   PrintStr("<title>@");
   PrintStr(GetName());
   PrintStr("@");
   PrintStr("</title>@");

   // Description
   PrintStr("<desc>@");
   PrintFast(22,"Creator: ROOT Version ");
   PrintStr(gROOT->GetVersion());
   PrintStr("@");
   PrintFast(14,"CreationDate: ");
   TDatime t;
   PrintStr(t.AsString());
   //Check a special header is defined in the current style
   Int_t nh = strlen(gStyle->GetHeaderPS());
   if (nh) {
      PrintFast(nh,gStyle->GetHeaderPS());
   }
   PrintStr("</desc>@");

   // Definitions
   PrintStr("<defs>@");
   PrintStr("</defs>@");

}

////////////////////////////////////////////////////////////////////////////////
/// Move to a new position (ix, iy). The move is done in relative coordinates
/// which allows to have short numbers which decrease the size of the file.
/// This function use the full power of the SVG's paths by using the
/// horizontal and vertical move whenever it is possible.

void TSVG::MovePS(Double_t ix, Double_t iy)
{
   if (ix != 0 && iy != 0)  {
      PrintFast(1,"l");
      WriteReal(ix);
      PrintFast(1,",");
      WriteReal(iy);
   } else if (ix != 0)  {
      PrintFast(1,"h");
      WriteReal(ix);
   } else if (iy != 0)  {
      PrintFast(1,"v");
      WriteReal(iy);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Start the SVG page. This function initialize the pad conversion
/// coefficients and output the <svg> directive which is close later in the
/// the function Close.

void TSVG::NewPage()
{
   // Compute pad conversion coefficients
   if (gPad) {
      Double_t ww   = gPad->GetWw();
      Double_t wh   = gPad->GetWh();
      fYsize        = fXsize*wh/ww;
   } else {
      fYsize = 27;
   }

   // <svg> directive. It defines the viewBox.
   if(!fBoundingBox) {
      PrintStr("@<?xml version=\"1.0\" standalone=\"no\"?>");
      PrintStr("@<svg width=\"");
      WriteReal(CMtoSVG(fXsize), kFALSE);
      PrintStr("\" height=\"");
      fYsizeSVG = CMtoSVG(fYsize);
      WriteReal(fYsizeSVG, kFALSE);
      PrintStr("\" viewBox=\"0 0");
      WriteReal(CMtoSVG(fXsize));
      WriteReal(fYsizeSVG);
      PrintStr("\" xmlns=\"http://www.w3.org/2000/svg\" shape-rendering=\"crispEdges\">");
      PrintStr("@");
      Initialize();
      fBoundingBox  = kTRUE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the range for the paper in centimetres

void TSVG::Range(Float_t xsize, Float_t ysize)
{
   Float_t xps, yps, xncm, yncm, dxwn, dywn, xwkwn, ywkwn, xymax;

   fXsize = xsize;
   fYsize = ysize;

   xps = xsize;
   yps = ysize;

   if( xsize <= xps && ysize < yps) {
      if ( xps > yps ) xymax = xps;
      else             xymax = yps;
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
         if( dxwn < 0) {
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
         if( dywn < 0) {
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
/// Set color index for fill areas

void TSVG::SetFillColor( Color_t cindex )
{
   fFillColor = cindex;
   if (gStyle->GetFillColor() <= 0) cindex = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set color index for lines

void TSVG::SetLineColor( Color_t cindex )
{
   fLineColor = cindex;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the value of the global parameter TSVG::fgLineJoin.
/// This parameter determines the appearance of joining lines in a SVG
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
/// gStyle->SetJoinLinePS(2); // Set the PS line join to bevel.
/// ~~~

void TSVG::SetLineJoin( Int_t linejoin )
{
   fgLineJoin = linejoin;
   if (fgLineJoin<0) fgLineJoin=0;
   if (fgLineJoin>2) fgLineJoin=2;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the value of the global parameter TSVG::fgLineCap.
/// This parameter determines the appearance of line caps in a SVG
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
/// gStyle->SetCapLinePS(2); // Set the PS line cap to projecting.
/// ~~~

void TSVG::SetLineCap( Int_t linecap )
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

void TSVG::SetLineStyle(Style_t linestyle)
{
   fLineStyle = linestyle;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the lines width.

void TSVG::SetLineWidth(Width_t linewidth)
{
   fLineWidth = linewidth;
}

////////////////////////////////////////////////////////////////////////////////
/// Set color index for markers.

void TSVG::SetMarkerColor( Color_t cindex )
{
   fMarkerColor = cindex;
}

////////////////////////////////////////////////////////////////////////////////
/// Set RGBa color with its color index

void TSVG::SetColorAlpha(Int_t color)
{
   if (color < 0) color = 0;
   TColor *col = gROOT->GetColor(color);
   if (col) {
      SetColor(col->GetRed(), col->GetGreen(), col->GetBlue());
      Float_t a = col->GetAlpha();
      if (a<1.) PrintStr(Form(" fill-opacity=\"%3.2f\" stroke-opacity=\"%3.2f\"",a,a));
   } else {
      SetColor(1., 1., 1.);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set RGB (without alpha channel) color with its color index

void TSVG::SetColor(Int_t color)
{
   if (color < 0) color = 0;
   TColor *col = gROOT->GetColor(color);
   if (col) {
      SetColor(col->GetRed(), col->GetGreen(), col->GetBlue());
   } else {
      SetColor(1., 1., 1.);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set color with its R G B components
///
///  - r: % of red in [0,1]
///  --g: % of green in [0,1]
///  - b: % of blue in [0,1]

void TSVG::SetColor(Float_t r, Float_t g, Float_t b)
{
   if (r <= 0. && g <= 0. && b <= 0. ) {
      PrintFast(7,"\"black\"");
   } else if (r >= 1. && g >= 1. && b >= 1. ) {
      PrintFast(7,"\"white\"");
   } else {
      char str[12];
      snprintf(str,12,"\"#%2.2x%2.2x%2.2x\"",Int_t(255.*r)
                                            ,Int_t(255.*g)
                                            ,Int_t(255.*b));
      PrintStr(str);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set color index for text

void TSVG::SetTextColor( Color_t cindex )
{
   fTextColor = cindex;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw text
///
///  - xx: x position of the text
///  - yy: y position of the text
///  - chars: text to be drawn

void TSVG::Text(Double_t xx, Double_t yy, const char *chars)
{
   static const char *fontFamily[] = {
   "Times"    , "Times"    , "Times",
   "Helvetica", "Helvetica", "Helvetica"   , "Helvetica",
   "Courier"  , "Courier"  , "Courier"     , "Courier",
   "Times"    ,"Times"     , "ZapfDingbats", "Times"};

   static const char *fontWeight[] = {
   "normal", "bold", "bold",
   "normal", "normal", "bold"  , "bold",
   "normal", "normal", "bold"  , "bold",
   "normal", "normal", "normal", "normal"};

   static const char *fontStyle[] = {
   "italic", "normal" , "italic",
   "normal", "oblique", "normal", "oblique",
   "normal", "oblique", "normal", "oblique",
   "normal", "normal" , "normal", "italic"};

   Double_t ix    = XtoSVG(xx);
   Double_t iy    = YtoSVG(yy);
   Double_t txalh = fTextAlign/10;
   if (txalh <1) txalh = 1; else if (txalh > 3) txalh = 3;
   Double_t txalv = fTextAlign%10;
   if (txalv <1) txalv = 1; else if (txalv > 3) txalv = 3;

   Double_t     wh = (Double_t)gPad->XtoPixel(gPad->GetX2());
   Double_t     hh = (Double_t)gPad->YtoPixel(gPad->GetY1());
   Float_t fontrap = 1.09; //scale down compared to X11
   Float_t ftsize;

   Int_t font  = abs(fTextFont)/10;
   if (font > 42 || font < 1) font = 1;
   if (wh < hh) {
      ftsize = fTextSize*fXsize*gPad->GetAbsWNDC();
   } else {
      ftsize = fTextSize*fYsize*gPad->GetAbsHNDC();
   }
   Int_t ifont = font-1;

   Double_t fontsize = CMtoSVG(ftsize/fontrap);
   if( fontsize <= 0) return;

   if (txalv == 3) iy = iy+fontsize;
   if (txalv == 2) iy = iy+(fontsize/2);

   if (fTextAngle != 0.) {
      PrintStr("@");
      PrintFast(21,"<g transform=\"rotate(");
      WriteReal(-fTextAngle, kFALSE);
      PrintFast(1,",");
      WriteReal(ix, kFALSE);
      PrintFast(1,",");
      WriteReal(iy, kFALSE);
      PrintFast(3,")\">");
   }

   PrintStr("@");
   PrintFast(30,"<text xml:space=\"preserve\" x=\"");
   WriteReal(ix, kFALSE);
   PrintFast(5,"\" y=\"");
   WriteReal(iy, kFALSE);
   PrintFast(1,"\"");
   if (txalh == 2) {
      PrintFast(21," text-anchor=\"middle\"");
   } else if (txalh == 3) {
      PrintFast(18," text-anchor=\"end\"");
   }
   PrintFast(6," fill=");
   SetColorAlpha(Int_t(fTextColor));
   PrintFast(12," font-size=\"");
   WriteReal(fontsize, kFALSE);
   PrintFast(15,"\" font-family=\"");
   PrintStr(fontFamily[ifont]);
   if (strcmp(fontWeight[ifont],"normal")) {
      PrintFast(15,"\" font-weight=\"");
      PrintStr(fontWeight[ifont]);
   }
   if (strcmp(fontStyle[ifont],"normal")) {
      PrintFast(14,"\" font-style=\"");
      PrintStr(fontStyle[ifont]);
   }
   PrintFast(2,"\">");

   if (font == 12 || font == 15) {
      Int_t ichar = chars[0]+848;
      Int_t ic    = ichar;

      // Math Symbols (cf: http://www.fileformat.info/info/unicode/category/Sm/list.htm)
      if (ic == 755) ichar =  8804;
      if (ic == 759) ichar =  9827;
      if (ic == 760) ichar =  9830;
      if (ic == 761) ichar =  9829;
      if (ic == 762) ichar =  9824;
      if (ic == 766) ichar =  8594;
      if (ic == 776) ichar =   247;
      if (ic == 757) ichar =  8734;
      if (ic == 758) ichar =   402;
      if (ic == 771) ichar =  8805;
      if (ic == 774) ichar =  8706;
      if (ic == 775) ichar =  8226;
      if (ic == 779) ichar =  8776;
      if (ic == 805) ichar =  8719;
      if (ic == 821) ichar =  8721;
      if (ic == 834) ichar =  8747;
      if (ic == 769) ichar =   177;
      if (ic == 772) ichar =   215;
      if (ic == 768) ichar =   176;
      if (ic == 791) ichar =  8745;
      if (ic == 793) ichar =  8835; // SUPERSET OF
      if (ic == 794) ichar =  8839; // SUPERSET OF OR EQUAL TO
      if (ic == 795) ichar =  8836; // NOT A SUBSET OF
      if (ic == 796) ichar =  8834;
      if (ic == 893) ichar =  8722;
      if (ic == 803) ichar =   169; // COPYRIGHT SIGN
      if (ic == 819) ichar =   169; // COPYRIGHT SIGN
      if (ic == 804) ichar =  8482;
      if (ic == 770) ichar =    34;
      if (ic == 823) ichar = 10072;
      if (ic == 781) ichar = 10072;
      if (ic == 824) ichar =  9117; // LEFT PARENTHESIS LOWER HOOK
      if (ic == 822) ichar =  9115; // LEFT PARENTHESIS UPPER HOOK
      if (ic == 767) ichar =  8595; // DOWNWARDS ARROW
      if (ic == 763) ichar =  8596; // LEFT RIGHT ARROW
      if (ic == 764) ichar =  8592; // LEFTWARDS ARROW
      if (ic == 788) ichar =  8855; // CIRCLED TIMES
      if (ic == 784) ichar =  8501;
      if (ic == 777) ichar =  8800;
      if (ic == 797) ichar =  8838;
      if (ic == 800) ichar =  8736;
      if (ic == 812) ichar =  8656; // LEFTWARDS DOUBLE ARROW
      if (ic == 817) ichar =    60; // LESS-THAN SIGN
      if (ic == 833) ichar =    62; // GREATER-THAN SIGN
      if (ic == 778) ichar =  8803; // STRICTLY EQUIVALENT TO
      if (ic == 809) ichar =  8743; // LOGICAL AND
      if (ic == 802) ichar =  9415; // CIRCLED LATIN CAPITAL LETTER R
      if (ic == 780) ichar =  8230; // HORIZONTAL ELLIPSIS
      if (ic == 801) ichar =  8711; // NABLA
      if (ic == 783) ichar =  8629; // DOWNWARDS ARROW WITH CORNER LEFTWARDS
      if (ic == 782) ichar =  8213;
      if (ic == 799) ichar =  8713;
      if (ic == 792) ichar =  8746;
      if (ic == 828) ichar =  9127;
      if (ic == 765) ichar =  8593; // UPWARDS ARROW
      if (ic == 789) ichar =  8853; // CIRCLED PLUS
      if (ic == 813) ichar =  8657; // UPWARDS DOUBLE ARROW
      if (ic == 773) ichar =  8733; // PROPORTIONAL TO
      if (ic == 790) ichar =  8709; // EMPTY SET
      if (ic == 810) ichar =  8744;
      if (ic == 756) ichar =  8260;
      if (ic == 807) ichar =  8231;
      if (ic == 808) ichar =  8989; // TOP RIGHT CORNER
      if (ic == 814) ichar =  8658; // RIGHTWARDS DOUBLE ARROW
      if (ic == 806) ichar =  8730; // SQUARE ROOT
      if (ic == 827) ichar =  9123;
      if (ic == 829) ichar =  9128;
      if (ic == 786) ichar =  8476;
      if (ic == 785) ichar =  8465;
      if (ic == 787) ichar =  8472;

      // Greek characters
      if (ic == 918) ichar = 934;
      if (ic == 919) ichar = 915;
      if (ic == 920) ichar = 919;
      if (ic == 923) ichar = 922;
      if (ic == 924) ichar = 923;
      if (ic == 925) ichar = 924;
      if (ic == 926) ichar = 925;
      if (ic == 929) ichar = 920;
      if (ic == 930) ichar = 929;
      if (ic == 936) ichar = 926;
      if (ic == 915) ichar = 935;
      if (ic == 937) ichar = 936;
      if (ic == 935) ichar = 937;
      if (ic == 938) ichar = 918;
      if (ic == 951) ichar = 947;
      if (ic == 798) ichar = 949;
      if (ic == 970) ichar = 950;
      if (ic == 952) ichar = 951;
      if (ic == 961) ichar = 952;
      if (ic == 955) ichar = 954;
      if (ic == 956) ichar = 955;
      if (ic == 957) ichar = 956;
      if (ic == 958) ichar = 957;
      if (ic == 968) ichar = 958;
      if (ic == 934) ichar = 962;
      if (ic == 962) ichar = 961;
      if (ic == 966) ichar = 969;
      if (ic == 950) ichar = 966;
      if (ic == 947) ichar = 967;
      if (ic == 969) ichar = 968;
      if (ic == 967) ichar = 969;
      if (ic == 954) ichar = 966;
      if (ic == 922) ichar = 952;
      if (ic == 753) ichar = 965;
      PrintStr(Form("&#%4.4d;",ichar));
   } else {
      Int_t len=strlen(chars);
      for (Int_t i=0; i<len;i++) {
         if (chars[i]!='\n') {
            if (chars[i]=='<') {
               PrintFast(4,"&lt;");
            } else if (chars[i]=='>') {
               PrintFast(4,"&gt;");
            } else if (chars[i]=='\305') {
               PrintFast(7,"&#8491;"); // ANGSTROM SIGN
            } else if (chars[i]=='\345') {
               PrintFast(6,"&#229;");
            } else if (chars[i]=='&') {
               PrintFast(5,"&amp;");
            } else {
               PrintFast(1,&chars[i]);
            }
         }
      }
   }

   PrintFast(7,"</text>");

   if (fTextAngle != 0.) {
      PrintStr("@");
      PrintFast(4,"</g>");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write a string of characters in NDC

void TSVG::TextNDC(Double_t u, Double_t v, const char *chars)
{
   Double_t x = gPad->GetX1() + u*(gPad->GetX2() - gPad->GetX1());
   Double_t y = gPad->GetY1() + v*(gPad->GetY2() - gPad->GetY1());
   Text(x, y, chars);
}

////////////////////////////////////////////////////////////////////////////////
/// Convert U from NDC coordinate to SVG

Double_t TSVG::UtoSVG(Double_t u)
{
   Double_t cm = fXsize*(gPad->GetAbsXlowNDC() + u*gPad->GetAbsWNDC());
   return 0.5 + 72*cm/2.54;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert V from NDC coordinate to SVG

Double_t TSVG::VtoSVG(Double_t v)
{
   Double_t cm = fYsize*(gPad->GetAbsYlowNDC() + v*gPad->GetAbsHNDC());
   return 0.5 + 72*cm/2.54;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert X from world coordinate to SVG

Double_t TSVG::XtoSVG(Double_t x)
{
   Double_t u = (x - gPad->GetX1())/(gPad->GetX2() - gPad->GetX1());
   return  UtoSVG(u);
}

////////////////////////////////////////////////////////////////////////////////
/// Convert Y from world coordinate to SVG

Double_t TSVG::YtoSVG(Double_t y)
{
   Double_t v = (y - gPad->GetY1())/(gPad->GetY2() - gPad->GetY1());
   return  fYsizeSVG-VtoSVG(v);
}

////////////////////////////////////////////////////////////////////////////////
/// Begin the Cell Array painting

void TSVG::CellArrayBegin(Int_t, Int_t, Double_t, Double_t, Double_t,
                          Double_t)
{
   Warning("TSVG::CellArrayBegin", "not yet implemented");
}

////////////////////////////////////////////////////////////////////////////////
/// Paint the Cell Array

void TSVG::CellArrayFill(Int_t, Int_t, Int_t)
{
   Warning("TSVG::CellArrayFill", "not yet implemented");
}

////////////////////////////////////////////////////////////////////////////////
/// End the Cell Array painting

void TSVG::CellArrayEnd()
{
   Warning("TSVG::CellArrayEnd", "not yet implemented");
}

////////////////////////////////////////////////////////////////////////////////
/// Not needed in SVG case

void TSVG::DrawPS(Int_t, Float_t *, Float_t *)
{
   Warning("TSVG::DrawPS", "not yet implemented");
}
