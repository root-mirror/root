// @(#)root/postscript:$Id$
// Author: O.Couet   16/07/99

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPostScript
#define ROOT_TPostScript


#include "TVirtualPS.h"

class TPoints;

class TPostScript : public TVirtualPS {

protected:
   Float_t fX1v;             ///< X bottom left corner of paper
   Float_t fY1v;             ///< Y bottom left corner of paper
   Float_t fX2v;             ///< X top right corner of paper
   Float_t fY2v;             ///< Y top right corner of paper
   Float_t fX1w;             ///<
   Float_t fY1w;             ///<
   Float_t fX2w;             ///<
   Float_t fY2w;             ///<
   Float_t fDXC;             ///<
   Float_t fDYC;             ///<
   Float_t fXC;              ///<
   Float_t fYC;              ///<
   Float_t fFX;              ///<
   Float_t fFY;              ///<
   Float_t fXVP1;            ///<
   Float_t fXVP2;            ///<
   Float_t fYVP1;            ///<
   Float_t fYVP2;            ///<
   Float_t fXVS1;            ///<
   Float_t fXVS2;            ///<
   Float_t fYVS1;            ///<
   Float_t fYVS2;            ///<
   Float_t fXsize;           ///< Page size along X
   Float_t fYsize;           ///< Page size along Y
   Float_t fMaxsize;         ///< Largest dimension of X and Y
   Float_t fRed;             ///< Per cent of red
   Float_t fGreen;           ///< Per cent of green
   Float_t fBlue;            ///< Per cent of blue
   Float_t fLineScale;       ///< Line width scale factor
   Int_t   fSave;            ///< Number of gsave for restore
   Int_t   fNXzone;          ///< Number of zones along X
   Int_t   fNYzone;          ///< Number of zones along Y
   Int_t   fIXzone;          ///< Current zone along X
   Int_t   fIYzone;          ///< Current zone along Y
   Float_t fMarkerSizeCur;   ///< current transformed value of marker size
   Int_t   fCurrentColor;    ///< current Postscript color index
   Int_t   fNpages;          ///< number of pages
   Int_t   fType;            ///< PostScript workstation type
   Int_t   fMode;            ///< PostScript mode
   Int_t   fClip;            ///< Clipping mode
   Bool_t  fBoundingBox;     ///< True for Encapsulated PostScript
   Bool_t  fClear;           ///< True when page must be cleared
   Bool_t  fClipStatus;      ///< Clipping Indicator
   Bool_t  fRange;           ///< True when a range has been defined
   Bool_t  fZone;            ///< Zone indicator
   char    fPatterns[32];    ///< Indicate if pattern n is defined
   Int_t   fNbinCT;          ///< Number of entries in the current Cell Array
   Int_t   fNbCellW;         ///< Number of boxes per line
   Int_t   fNbCellLine;      ///< Number of boxes in the current line
   Int_t   fMaxLines;        ///< Maximum number of lines in a PS array
   Int_t   fLastCellRed;     ///< Last red value
   Int_t   fLastCellGreen;   ///< Last green value
   Int_t   fLastCellBlue;    ///< Last blue value
   Int_t   fNBSameColorCell; ///< Number of boxes with the same color
   TString fFileName;        ///< PS file name
   Bool_t  fFontEmbed;       ///< True is FontEmbed has been called

   static Int_t fgLineJoin;  ///< Appearance of joining lines
   static Int_t fgLineCap;   ///< Appearance of line caps

public:
   TPostScript();
   TPostScript(const char *filename, Int_t type=-111);
   virtual ~TPostScript();

   void  CellArrayBegin(Int_t W, Int_t H, Double_t x1, Double_t x2,
                                          Double_t y1, Double_t y2);
   void  CellArrayFill(Int_t r, Int_t g, Int_t b);
   void  CellArrayEnd();
   void  Close(Option_t *opt="");
   Int_t CMtoPS(Double_t u) {return Int_t(0.5 + 72*u/2.54);}
   void  DefineMarkers();
   void  DrawBox(Double_t x1, Double_t y1,Double_t x2, Double_t  y2);
   void  DrawFrame(Double_t xl, Double_t yl, Double_t xt, Double_t  yt,
                   Int_t mode, Int_t border, Int_t dark, Int_t light);
   void  DrawHatch(Float_t dy, Float_t angle, Int_t n, Float_t *x,
                   Float_t *y);
   void  DrawHatch(Float_t dy, Float_t angle, Int_t n, Double_t *x,
                   Double_t *y);
   void  DrawPolyLine(Int_t n, TPoints *xy);
   void  DrawPolyLineNDC(Int_t n, TPoints *uv);
   void  DrawPolyMarker(Int_t n, Float_t *x, Float_t *y);
   void  DrawPolyMarker(Int_t n, Double_t *x, Double_t *y);
   void  DrawPS(Int_t n, Float_t *xw, Float_t *yw);
   void  DrawPS(Int_t n, Double_t *xw, Double_t *yw);
   bool  FontEmbedType1(const char *filename);
   bool  FontEmbedType2(const char *filename);
   bool  FontEmbedType42(const char *filename);
   void  FontEmbed();
   void  FontEncode();
   void  Initialize();
   void  NewPage();
   void  Off();
   void  On();
   void  Open(const char *filename, Int_t type=-111);
   void  SaveRestore(Int_t flag);
   void  SetFillColor( Color_t cindex=1);
   void  SetFillPatterns(Int_t ipat, Int_t color);
   void  SetLineColor( Color_t cindex=1);
   void  SetLineJoin(Int_t linejoin=0);
   void  SetLineCap(Int_t linecap=0);
   void  SetLineStyle(Style_t linestyle = 1);
   void  SetLineWidth(Width_t linewidth = 1);
   void  SetLineScale(Float_t scale=3) {fLineScale = scale;}
   void  SetMarkerColor( Color_t cindex=1);
   void  SetTextColor( Color_t cindex=1);
   void  MovePS(Int_t x, Int_t y);
   void  Range(Float_t xrange, Float_t yrange);
   void  SetColor(Int_t color = 1);
   void  SetColor(Float_t r, Float_t g, Float_t b);
   void  Text(Double_t x, Double_t y, const char *string);
   void  Text(Double_t x, Double_t y, const wchar_t *string);
   void  TextNDC(Double_t u, Double_t v, const char *string);
   void  TextNDC(Double_t u, Double_t v, const wchar_t *string);
   Int_t UtoPS(Double_t u);
   Int_t VtoPS(Double_t v);
   Int_t XtoPS(Double_t x);
   Int_t YtoPS(Double_t y);
   void  Zone();

   ClassDef(TPostScript,0)  //PostScript driver
};

#endif
