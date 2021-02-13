// @(#)root/postscript:$Id$
// Author: Olivier Couet

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TSVG
#define ROOT_TSVG


#include "TVirtualPS.h"

class TPoints;

class TSVG : public TVirtualPS {

protected:
   Float_t      fXsize;           ///< Page size along X
   Float_t      fYsize;           ///< Page size along Y
   Int_t        fType;            ///< Workstation type used to know if the SVG is open
   Bool_t       fBoundingBox;     ///< True when the SVG header is printed
   Bool_t       fRange;           ///< True when a range has been defined
   Double_t     fYsizeSVG;        ///< Page's Y size in SVG units

   static Int_t fgLineJoin;       ///< Appearance of joining lines
   static Int_t fgLineCap;        ///< Appearance of line caps

public:
   TSVG();
   TSVG(const char *filename, Int_t type=-113);
   virtual ~TSVG();

   void  CellArrayBegin(Int_t W, Int_t H, Double_t x1, Double_t x2, Double_t y1, Double_t y2);
   void  CellArrayFill(Int_t r, Int_t g, Int_t b);
   void  CellArrayEnd();
   void  Close(Option_t *opt="");
   Double_t CMtoSVG(Double_t u) {return 0.5 + 72*u/2.54;}
   void  DrawBox(Double_t x1, Double_t y1,Double_t x2, Double_t  y2);
   void  DrawFrame(Double_t xl, Double_t yl, Double_t xt, Double_t  yt,
                   Int_t mode, Int_t border, Int_t dark, Int_t light);
   void  DrawPolyLine(Int_t n, TPoints *xy);
   void  DrawPolyLineNDC(Int_t n, TPoints *uv);
   void  DrawPolyMarker(Int_t n, Float_t *x, Float_t *y);
   void  DrawPolyMarker(Int_t n, Double_t *x, Double_t *y);
   void  DrawPS(Int_t n, Float_t *xw, Float_t *yw);
   void  DrawPS(Int_t n, Double_t *xw, Double_t *yw);
   void  Initialize();
   void  MovePS(Double_t x, Double_t y);
   void  NewPage();
   void  Off();
   void  On();
   void  Open(const char *filename, Int_t type=-111);
   void  Range(Float_t xrange, Float_t yrange);
   void  SetColorAlpha(Int_t color = 1);
   void  SetColor(Int_t color = 1);
   void  SetColor(Float_t r, Float_t g, Float_t b);
   void  SetFillColor( Color_t cindex=1);
   void  SetLineColor( Color_t cindex=1);
   void  SetLineJoin(Int_t linejoin=0);
   void  SetLineCap(Int_t linecap=0);
   void  SetLineStyle(Style_t linestyle = 1);
   void  SetLineWidth(Width_t linewidth = 1);
   void  SetLineScale(Float_t =3) { }
   void  SetMarkerColor( Color_t cindex=1);
   void  SetTextColor( Color_t cindex=1);
   void  Text(Double_t x, Double_t y, const char *string);
   void  Text(Double_t, Double_t, const wchar_t *){}
   void  TextNDC(Double_t u, Double_t v, const char *string);
   void  TextNDC(Double_t, Double_t, const wchar_t *){}
   Double_t UtoSVG(Double_t u);
   Double_t VtoSVG(Double_t v);
   Double_t XtoSVG(Double_t x);
   Double_t YtoSVG(Double_t y);

   ClassDef(TSVG, 0);  //SVG driver
};

#endif
