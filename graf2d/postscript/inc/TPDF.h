// @(#)root/postscript:$Id: TPDF.h,v 1.0
// Author: Olivier Couet

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPDF
#define ROOT_TPDF


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TPDF                                                                 //
//                                                                      //
// PDF driver.                                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TVirtualPS.h"
#include <vector>


class TPoints;

class TPDF : public TVirtualPS {

protected:
   Float_t            fRed;             ///< Per cent of red
   Float_t            fGreen;           ///< Per cent of green
   Float_t            fBlue;            ///< Per cent of blue
   Float_t            fAlpha;           ///< Per cent of transparency
   std::vector<float> fAlphas;          ///< List of alpha values used
   Float_t            fXsize;           ///< Page size along X
   Float_t            fYsize;           ///< Page size along Y
   Int_t              fType;            ///< Workstation type used to know if the PDF is open
   Int_t              fPageFormat;      ///< Page format (A4, Letter etc ...)
   Int_t              fPageOrientation; ///< Page orientation (Portrait, Landscape)
   Int_t              fStartStream;     ///<
   Float_t            fLineScale;       ///< Line width scale factor
   Int_t             *fObjPos;          ///< Objets position
   Int_t              fObjPosSize;      ///< Real size of fObjPos
   Int_t              fNbObj;           ///< Number of objects
   Int_t              fNbPage;          ///< Number of pages
   Bool_t             fPageNotEmpty;    ///< True if the current page is not empty
   Bool_t             fCompress;        ///< True when fBuffer must be compressed
   Bool_t             fRange;           ///< True when a range has been defined

   static Int_t       fgLineJoin;       ///< Appearance of joining lines
   static Int_t       fgLineCap;        ///< Appearance of line caps

public:
   TPDF();
   TPDF(const char *filename, Int_t type=-111);
   virtual ~TPDF();

   void     CellArrayBegin(Int_t W, Int_t H, Double_t x1, Double_t x2, Double_t y1, Double_t y2);
   void     CellArrayFill(Int_t r, Int_t g, Int_t b);
   void     CellArrayEnd();
   void     Close(Option_t *opt="");
   Double_t CMtoPDF(Double_t u) {return Int_t(0.5 + 72*u/2.54);}
   void     DrawBox(Double_t x1, Double_t y1,Double_t x2, Double_t  y2);
   void     DrawFrame(Double_t xl, Double_t yl, Double_t xt, Double_t  yt,
                      Int_t mode, Int_t border, Int_t dark, Int_t light);
   void     DrawHatch(Float_t dy, Float_t angle, Int_t n, Float_t *x, Float_t *y);
   void     DrawHatch(Float_t dy, Float_t angle, Int_t n, Double_t *x, Double_t *y);
   void     DrawPolyLine(Int_t n, TPoints *xy);
   void     DrawPolyLineNDC(Int_t n, TPoints *uv);
   void     DrawPolyMarker(Int_t n, Float_t *x, Float_t *y);
   void     DrawPolyMarker(Int_t n, Double_t *x, Double_t *y);
   void     DrawPS(Int_t n, Float_t *xw, Float_t *yw);
   void     DrawPS(Int_t n, Double_t *xw, Double_t *yw);
   void     LineTo(Double_t x, Double_t y);
   void     MoveTo(Double_t x, Double_t y);
   void     FontEncode();
   void     NewObject(Int_t n);
   void     NewPage();
   void     Off();
   void     On();
   void     Open(const char *filename, Int_t type=-111);
   void     PatternEncode();
   void     PrintFast(Int_t nch, const char *string="");
   void     PrintStr(const char *string="");
   void     Range(Float_t xrange, Float_t yrange);
   void     SetAlpha(Float_t alpha = 1.);
   void     SetColor(Int_t color = 1);
   void     SetColor(Float_t r, Float_t g, Float_t b);
   void     SetFillColor( Color_t cindex=1);
   void     SetFillPatterns(Int_t ipat, Int_t color);
   void     SetLineColor( Color_t cindex=1);
   void     SetLineJoin(Int_t linejoin=0);
   void     SetLineCap(Int_t linecap=0);
   void     SetLineScale(Float_t scale=1) {fLineScale = scale;}
   void     SetLineStyle(Style_t linestyle = 1);
   void     SetLineWidth(Width_t linewidth = 1);
   void     SetMarkerColor( Color_t cindex=1);
   void     SetTextColor( Color_t cindex=1);
   void     Text(Double_t x, Double_t y, const char *string);
   void     Text(Double_t, Double_t, const wchar_t *);
   void     TextNDC(Double_t u, Double_t v, const char *string);
   void     TextNDC(Double_t, Double_t, const wchar_t *);
   void     WriteCompressedBuffer();
   virtual  void WriteReal(Float_t r, Bool_t space=kTRUE);
   Double_t UtoPDF(Double_t u);
   Double_t VtoPDF(Double_t v);
   Double_t XtoPDF(Double_t x);
   Double_t YtoPDF(Double_t y);

   ClassDef(TPDF, 0);  //PDF driver
};

#endif
