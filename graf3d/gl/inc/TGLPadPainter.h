// @(#)root/gl:$Id$
// Author:  Timur Pocheptsov  06/05/2009

/*************************************************************************
 * Copyright (C) 1995-2009, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGLPadPainter
#define ROOT_TGLPadPainter

#include "TVirtualPadPainter.h"
#include "TGLFontManager.h"
#include "TGLPadUtils.h"
#include "TPoint.h"

#include <vector>

class TLinearGradient;
class TRadialGradient;
/*
The _main_ purpose of TGLPadPainter is to enable 2d gl raphics
inside standard TPad/TCanvas.
*/
class TGLPadPainter : public TVirtualPadPainter {
private:
   Rgl::Pad::PolygonStippleSet fSSet;
   Rgl::Pad::Tesselator        fTess;
   Rgl::Pad::MarkerPainter     fMarker;
   Rgl::Pad::GLLimits          fLimits;

   std::vector<Double_t>       fVs;//Vertex buffer for tesselator.

   TGLFontManager              fFM;
   TGLFont                     fF;

   Int_t                       fVp[4];

   std::vector<TPoint>         fPoly;
   Bool_t                      fIsHollowArea;

   Bool_t                      fLocked;

   template<class Char_t>
   void DrawTextHelper(Double_t x, Double_t y, const Char_t *text, ETextMode mode);
public:
   TGLPadPainter();

   //Final overriders for TVirtualPadPainter pure virtual functions.
   //1. Part, which simply delegates to TVirtualX.
   //Line attributes.
   Color_t  GetLineColor() const;
   Style_t  GetLineStyle() const;
   Width_t  GetLineWidth() const;

   void     SetLineColor(Color_t lcolor);
   void     SetLineStyle(Style_t lstyle);
   void     SetLineWidth(Width_t lwidth);
   //Fill attributes.
   Color_t  GetFillColor() const;
   Style_t  GetFillStyle() const;
   Bool_t   IsTransparent() const;

   void     SetFillColor(Color_t fcolor);
   void     SetFillStyle(Style_t fstyle);
   void     SetOpacity(Int_t percent);
   //Text attributes.
   Short_t  GetTextAlign() const;
   Float_t  GetTextAngle() const;
   Color_t  GetTextColor() const;
   Font_t   GetTextFont()  const;
   Float_t  GetTextSize()  const;
   Float_t  GetTextMagnitude() const;

   void     SetTextAlign(Short_t align);
   void     SetTextAngle(Float_t tangle);
   void     SetTextColor(Color_t tcolor);
   void     SetTextFont(Font_t tfont);
   void     SetTextSize(Float_t tsize);
   void     SetTextSizePixels(Int_t npixels);

   //2. "Off-screen management" part.
   Int_t    CreateDrawable(UInt_t w, UInt_t h);
   void     ClearDrawable();
   void     CopyDrawable(Int_t device, Int_t px, Int_t py);
   void     DestroyDrawable(Int_t device);
   void     SelectDrawable(Int_t device);

   void     InitPainter();
   void     InvalidateCS();
   void     LockPainter();

   void     DrawLine(Double_t x1, Double_t y1, Double_t x2, Double_t y2);
   void     DrawLineNDC(Double_t u1, Double_t v1, Double_t u2, Double_t v2);

   void     DrawBox(Double_t x1, Double_t y1, Double_t x2, Double_t y2, EBoxMode mode);
   //TPad needs double and float versions.
   void     DrawFillArea(Int_t n, const Double_t *x, const Double_t *y);
   void     DrawFillArea(Int_t n, const Float_t *x, const Float_t *y);

   //TPad needs both double and float versions of DrawPolyLine.
   void     DrawPolyLine(Int_t n, const Double_t *x, const Double_t *y);
   void     DrawPolyLine(Int_t n, const Float_t *x, const Float_t *y);
   void     DrawPolyLineNDC(Int_t n, const Double_t *u, const Double_t *v);

   //TPad needs both versions.
   void     DrawPolyMarker(Int_t n, const Double_t *x, const Double_t *y);
   void     DrawPolyMarker(Int_t n, const Float_t *x, const Float_t *y);

   void     DrawText(Double_t x, Double_t y, const char *text, ETextMode mode);
   void     DrawText(Double_t, Double_t, const wchar_t *, ETextMode);
   void     DrawTextNDC(Double_t x, Double_t y, const char *text, ETextMode mode);
   void     DrawTextNDC(Double_t, Double_t, const wchar_t *, ETextMode);

   //jpg, png, gif and bmp output.
   void     SaveImage(TVirtualPad *pad, const char *fileName, Int_t type) const;

   //TASImage support.
   void     DrawPixels(const unsigned char *pixelData, UInt_t width, UInt_t height,
                       Int_t dstX, Int_t dstY, Bool_t enableBlending);


private:

   //Attention! GL_PROJECTION will become
   //the current matrix after these calls.
   void     SaveProjectionMatrix()const;
   void     RestoreProjectionMatrix()const;

   //Attention! GL_MODELVIEW will become the
   //current matrix after these calls.
   void     SaveModelviewMatrix()const;
   void     RestoreModelviewMatrix()const;

   void     SaveViewport();
   void     RestoreViewport();

   void     DrawPolyMarker();

   //Aux. functions for a gradient and solid fill:
   void DrawPolygonWithGradient(Int_t n, const Double_t *x, const Double_t *y);
   //
   void DrawGradient(const TLinearGradient *gradient, Int_t n, const Double_t *x, const Double_t *y);
   void DrawGradient(const TRadialGradient *gradient, Int_t n, const Double_t *x, const Double_t *y);
   //
   void DrawTesselation(Int_t n, const Double_t *x, const Double_t *y);

   TGLPadPainter(const TGLPadPainter &rhs);
   TGLPadPainter & operator = (const TGLPadPainter &rhs);

   ClassDef(TGLPadPainter, 0)
};

#endif

