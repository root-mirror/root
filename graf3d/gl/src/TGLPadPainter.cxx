// @(#)root/gl:$Id$
// Author:  Timur Pocheptsov  06/05/2009

/*************************************************************************
 * Copyright (C) 1995-2009, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <stdexcept>
#include <cassert>
#include <limits>
#include <memory>
#include <vector>

#include "TAttMarker.h"
#include "TVirtualX.h"
#include "TError.h"
#include "TImage.h"
#include "TROOT.h"
#include "TPad.h"

#include "TColorGradient.h"
#include "TGLPadPainter.h"
#include "TGLIncludes.h"
#include "TGLUtil.h"
#include "TMath.h"

namespace {

////////////////////////////////////////////////////////////////////////////////
///Not a bad idea to assert on gVirtualX != nullptr

bool IsGradientFill(Color_t fillColorIndex)
{
   return dynamic_cast<TColorGradient *>(gROOT->GetColor(fillColorIndex));
}

}

/** \class TGLPadPainter
\ingroup opengl
"Delegating" part of TGLPadPainter. Line/fill/etc. attributes can be
set inside TPad, but not only there:
many of them are set by base sub-objects of 2d primitives
(2d primitives usually inherit TAttLine or TAttFill etc.).  And these sub-objects
call gVirtualX->SetLineWidth ... etc. So, if I save some attributes in my painter,
it will be mess - at any moment I do not know, where to take line attribute - from
gVirtualX or from my own member. So! All attributed, _ALL_ go to/from gVirtualX.
*/

ClassImp(TGLPadPainter);

////////////////////////////////////////////////////////////////////////////////

TGLPadPainter::TGLPadPainter()
                  : fIsHollowArea(kFALSE),
                    fLocked(kTRUE)
{
   fVp[0] = fVp[1] = fVp[2] = fVp[3] = 0;
}


////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

Color_t TGLPadPainter::GetLineColor() const
{
   return gVirtualX->GetLineColor();
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

Style_t TGLPadPainter::GetLineStyle() const
{
   return gVirtualX->GetLineStyle();
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

Width_t TGLPadPainter::GetLineWidth() const
{
   return gVirtualX->GetLineWidth();
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

void TGLPadPainter::SetLineColor(Color_t lcolor)
{
   gVirtualX->SetLineColor(lcolor);
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

void TGLPadPainter::SetLineStyle(Style_t lstyle)
{
   gVirtualX->SetLineStyle(lstyle);
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

void TGLPadPainter::SetLineWidth(Width_t lwidth)
{
   gVirtualX->SetLineWidth(lwidth);
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

Color_t TGLPadPainter::GetFillColor() const
{
   return gVirtualX->GetFillColor();
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

Style_t TGLPadPainter::GetFillStyle() const
{
   return gVirtualX->GetFillStyle();
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.
///IsTransparent is implemented as inline function in TAttFill.

Bool_t TGLPadPainter::IsTransparent() const
{
   return gVirtualX->IsTransparent();
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

void TGLPadPainter::SetFillColor(Color_t fcolor)
{
   gVirtualX->SetFillColor(fcolor);
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

void TGLPadPainter::SetFillStyle(Style_t fstyle)
{
   gVirtualX->SetFillStyle(fstyle);
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

void TGLPadPainter::SetOpacity(Int_t percent)
{
   gVirtualX->SetOpacity(percent);
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

Short_t TGLPadPainter::GetTextAlign() const
{
   return gVirtualX->GetTextAlign();
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

Float_t TGLPadPainter::GetTextAngle() const
{
   return gVirtualX->GetTextAngle();
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

Color_t TGLPadPainter::GetTextColor() const
{
   return gVirtualX->GetTextColor();
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

Font_t TGLPadPainter::GetTextFont() const
{
   return gVirtualX->GetTextFont();
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

Float_t TGLPadPainter::GetTextSize() const
{
   return gVirtualX->GetTextSize();
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

Float_t TGLPadPainter::GetTextMagnitude() const
{
   return gVirtualX->GetTextMagnitude();
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

void TGLPadPainter::SetTextAlign(Short_t align)
{
   gVirtualX->SetTextAlign(align);
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

void TGLPadPainter::SetTextAngle(Float_t tangle)
{
   gVirtualX->SetTextAngle(tangle);
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

void TGLPadPainter::SetTextColor(Color_t tcolor)
{
   gVirtualX->SetTextColor(tcolor);
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

void TGLPadPainter::SetTextFont(Font_t tfont)
{
   gVirtualX->SetTextFont(tfont);
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

void TGLPadPainter::SetTextSize(Float_t tsize)
{
   gVirtualX->SetTextSize(tsize);
}

////////////////////////////////////////////////////////////////////////////////
///Delegate to gVirtualX.

void TGLPadPainter::SetTextSizePixels(Int_t npixels)
{
   gVirtualX->SetTextSizePixels(npixels);
}

/*
"Pixmap" part of TGLPadPainter.
*/

////////////////////////////////////////////////////////////////////////////////
///Not required at the moment.

Int_t TGLPadPainter::CreateDrawable(UInt_t/*w*/, UInt_t/*h*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
///Not required at the moment.

void TGLPadPainter::ClearDrawable()
{
}

////////////////////////////////////////////////////////////////////////////////
///Not required at the moment.

void TGLPadPainter::CopyDrawable(Int_t /*device*/, Int_t /*px*/, Int_t /*py*/)
{
}

////////////////////////////////////////////////////////////////////////////////
///Not required at the moment.

void TGLPadPainter::DestroyDrawable(Int_t /*device*/)
{
}

////////////////////////////////////////////////////////////////////////////////
///For gVirtualX this means select pixmap (or window)
///and all subsequent drawings will go into
///this pixmap. For OpenGL this means the change of
///coordinate system and viewport.

void TGLPadPainter::SelectDrawable(Int_t /*device*/)
{
   if (fLocked)
      return;

   if (TPad *pad = dynamic_cast<TPad *>(gPad)) {
      Int_t px = 0, py = 0;

      pad->XYtoAbsPixel(pad->GetX1(), pad->GetY1(), px, py);

      py = gPad->GetWh() - py;
      //
      TGLUtil::InitializeIfNeeded();
      const Float_t scale = TGLUtil::GetScreenScalingFactor();

      glViewport(GLint(px * scale), GLint(py * scale),
                 GLsizei(gPad->GetWw() * pad->GetAbsWNDC() * scale),
                 GLsizei(gPad->GetWh() * pad->GetAbsHNDC() * scale));

      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      glOrtho(pad->GetX1(), pad->GetX2(), pad->GetY1(), pad->GetY2(), -10., 10.);

      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
      glTranslated(0., 0., -1.);
   } else {
      ::Error("TGLPadPainter::SelectDrawable",
               "function was called not from TPad or TCanvas code\n");
      throw std::runtime_error("");
   }
}

////////////////////////////////////////////////////////////////////////////////
///Init gl-pad painter:
///1. 2D painter does not use depth test, should not modify
///   depth-buffer content (except initial cleanup).
///2. Disable cull face.
///3. Disable lighting.
///4. Set viewport (to the whole canvas area).
///5. Set camera.
///6. Unlock painter.

void TGLPadPainter::InitPainter()
{
   glDisable(GL_DEPTH_TEST);
   glDisable(GL_CULL_FACE);
   glDisable(GL_LIGHTING);

   //Clear the buffer
   glViewport(0, 0, GLsizei(gPad->GetWw()), GLsizei(gPad->GetWh()));

   glDepthMask(GL_TRUE);
   glClearColor(1.,1.,1.,1.);
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   glDepthMask(GL_FALSE);

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();

   glOrtho(gPad->GetX1(), gPad->GetX2(), gPad->GetY1(), gPad->GetY2(), -10., 10.);

   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   glTranslated(0., 0., -1.);

   fLocked = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
///When TPad::Range for gPad is called, projection
///must be changed in OpenGL.

void TGLPadPainter::InvalidateCS()
{
   if (fLocked) return;

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();

   glOrtho(gPad->GetX1(), gPad->GetX2(), gPad->GetY1(), gPad->GetY2(), -10., 10.);

   glMatrixMode(GL_MODELVIEW);
}

////////////////////////////////////////////////////////////////////////////////
///Locked state of painter means, that
///GL context can be invalid, so no GL calls
///can be executed.

void TGLPadPainter::LockPainter()
{
   if (fLocked) return;

   glFinish();
   fLocked = kTRUE;
}

/*
2D primitives.
*/

const Double_t lineWidthTS = 3.;

////////////////////////////////////////////////////////////////////////////////
///Draw line segment.

void TGLPadPainter::DrawLine(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
   if (fLocked) {
      //GL pad painter can be called in non-standard situation:
      //not from TPad::Paint, but
      //from TView3D::ExecuteRotateView. This means in fact,
      //that TView3D wants to draw itself in a XOR mode, via
      //gVirtualX.
      if (gVirtualX->GetDrawMode() == TVirtualX::kInvert) {
         gVirtualX->DrawLine(gPad->XtoAbsPixel(x1), gPad->YtoAbsPixel(y1),
                             gPad->XtoAbsPixel(x2), gPad->YtoAbsPixel(y2));
      }

      return;
   }

   const Rgl::Pad::LineAttribSet lineAttribs(kTRUE, gVirtualX->GetLineStyle(), fLimits.GetMaxLineWidth(), kTRUE);

   glBegin(GL_LINES);
   glVertex2d(x1, y1);
   glVertex2d(x2, y2);
   glEnd();

   if (gVirtualX->GetLineWidth() > lineWidthTS) {
      Double_t pointSize = gVirtualX->GetLineWidth();
      if (pointSize > fLimits.GetMaxPointSize())
         pointSize = fLimits.GetMaxPointSize();
      glPointSize((GLfloat)pointSize);
      const TGLEnableGuard pointSmooth(GL_POINT_SMOOTH);
      glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
      glBegin(GL_POINTS);

      glVertex2d(x1, y1);
      glVertex2d(x2, y2);

      glEnd();
      glPointSize(1.f);
   }

}

////////////////////////////////////////////////////////////////////////////////
///Draw line segment in NDC coordinates.

void TGLPadPainter::DrawLineNDC(Double_t u1, Double_t v1, Double_t u2, Double_t v2)
{
   if (fLocked) return;

   const Rgl::Pad::LineAttribSet lineAttribs(kTRUE, gVirtualX->GetLineStyle(), fLimits.GetMaxLineWidth(), kTRUE);
   const Double_t xRange = gPad->GetX2() - gPad->GetX1();
   const Double_t yRange = gPad->GetY2() - gPad->GetY1();

   glBegin(GL_LINES);
   glVertex2d(gPad->GetX1() + u1 * xRange, gPad->GetY1() + v1 * yRange);
   glVertex2d(gPad->GetX1() + u2 * xRange, gPad->GetY1() + v2 * yRange);
   glEnd();
}

////////////////////////////////////////////////////////////////////////////////
///Draw filled or hollow box.

void TGLPadPainter::DrawBox(Double_t x1, Double_t y1, Double_t x2, Double_t y2, EBoxMode mode)
{
   if (fLocked) return;

   if (IsGradientFill(gVirtualX->GetFillColor())) {
      Double_t xs[] = {x1, x2, x2, x1};
      Double_t ys[] = {y1, y1, y2, y2};
      DrawPolygonWithGradient(4, xs, ys);
      return;
   }

   if (mode == kHollow) {
      const Rgl::Pad::LineAttribSet lineAttribs(kTRUE, 0, fLimits.GetMaxLineWidth(), kTRUE);
      //
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      glRectd(x1, y1, x2, y2);
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      glLineWidth(1.f);
   } else {
      const Rgl::Pad::FillAttribSet fillAttribs(fSSet, kFALSE);//Set filling parameters.
      glRectd(x1, y1, x2, y2);
   }
}

////////////////////////////////////////////////////////////////////////////////
///Draw tesselated polygon (probably, outline only).

void TGLPadPainter::DrawFillArea(Int_t n, const Double_t *x, const Double_t *y)
{
   assert(x != 0 && "DrawFillArea, parameter 'x' is null");
   assert(y != 0 && "DrawFillArea, parameter 'y' is null");

   if (fLocked)
      return;

   if (n < 3) {
      ::Error("TGLPadPainter::DrawFillArea",
              "invalid number of points in a polygon");
      return;
   }

   if (IsGradientFill(gVirtualX->GetFillColor()))
      return DrawPolygonWithGradient(n, x, y);

   if (!gVirtualX->GetFillStyle()) {
      fIsHollowArea = kTRUE;
      return DrawPolyLine(n, x, y);
   }

   const Rgl::Pad::FillAttribSet fillAttribs(fSSet, kFALSE);
   DrawTesselation(n, x, y);
}

////////////////////////////////////////////////////////////////////////////////
///Draw tesselated polygon (never called, probably, since TPad::PaintFillArea for floats
///is deprecated).

void TGLPadPainter::DrawFillArea(Int_t n, const Float_t *x, const Float_t *y)
{
   if (fLocked) return;

   if (!gVirtualX->GetFillStyle()) {
      fIsHollowArea = kTRUE;
      return DrawPolyLine(n, x, y);
   }

   fVs.resize(n * 3);

   for (Int_t i = 0; i < n; ++i) {
      fVs[i * 3]     = x[i];
      fVs[i * 3 + 1] = y[i];
   }

   const Rgl::Pad::FillAttribSet fillAttribs(fSSet, kFALSE);

   GLUtesselator *t = (GLUtesselator *)fTess.GetTess();
   gluBeginPolygon(t);
   gluNextContour(t, (GLenum)GLU_UNKNOWN);

   for (Int_t i = 0; i < n; ++i)
      gluTessVertex(t, &fVs[i * 3], &fVs[i * 3]);


   gluEndPolygon(t);
}

////////////////////////////////////////////////////////////////////////////////
///Draw poly-line in user coordinates.

void TGLPadPainter::DrawPolyLine(Int_t n, const Double_t *x, const Double_t *y)
{
   if (fLocked) return;

   const Rgl::Pad::LineAttribSet lineAttribs(kTRUE, gVirtualX->GetLineStyle(), fLimits.GetMaxLineWidth(), kTRUE);

   glBegin(GL_LINE_STRIP);

   for (Int_t i = 0; i < n; ++i)
      glVertex2d(x[i], y[i]);

   if (fIsHollowArea) {
      glVertex2d(x[0], y[0]);
      fIsHollowArea = kFALSE;
   }
   glEnd();

   if (gVirtualX->GetLineWidth() > lineWidthTS) {
      Double_t pointSize = gVirtualX->GetLineWidth();
      if (pointSize > fLimits.GetMaxPointSize())
         pointSize = fLimits.GetMaxPointSize();
      glPointSize((GLfloat)pointSize);
      const TGLEnableGuard pointSmooth(GL_POINT_SMOOTH);
      glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
      glBegin(GL_POINTS);

      for (Int_t i = 0; i < n; ++i)
         glVertex2d(x[i], y[i]);

      glEnd();
      glPointSize(1.f);
   }
}

////////////////////////////////////////////////////////////////////////////////
///Never called?

void TGLPadPainter::DrawPolyLine(Int_t n, const Float_t *x, const Float_t *y)
{
   if (fLocked) return;

   const Rgl::Pad::LineAttribSet lineAttribs(kTRUE, gVirtualX->GetLineStyle(), fLimits.GetMaxLineWidth(), kTRUE);

   glBegin(GL_LINE_STRIP);

   for (Int_t i = 0; i < n; ++i)
      glVertex2f(x[i], y[i]);

   if (fIsHollowArea) {
      glVertex2f(x[0], y[0]);
      fIsHollowArea = kFALSE;
   }

   glEnd();
}

////////////////////////////////////////////////////////////////////////////////
///Poly line in NDC.

void TGLPadPainter::DrawPolyLineNDC(Int_t n, const Double_t *u, const Double_t *v)
{
   if (fLocked) return;

   const Rgl::Pad::LineAttribSet lineAttribs(kTRUE, gVirtualX->GetLineStyle(), fLimits.GetMaxLineWidth(), kTRUE);
   const Double_t xRange = gPad->GetX2() - gPad->GetX1();
   const Double_t yRange = gPad->GetY2() - gPad->GetY1();
   const Double_t x1 = gPad->GetX1(), y1 = gPad->GetY1();

   glBegin(GL_LINE_STRIP);

   for (Int_t i = 0; i < n; ++i)
      glVertex2d(x1 + u[i] * xRange, y1 + v[i] * yRange);

   glEnd();
}

namespace {

//Aux. function.
template<class ValueType>
void ConvertMarkerPoints(Int_t n, const ValueType *x, const ValueType *y, std::vector<TPoint> & dst);

}

////////////////////////////////////////////////////////////////////////////////
///Poly-marker.

void TGLPadPainter::DrawPolyMarker(Int_t n, const Double_t *x, const Double_t *y)
{
   if (fLocked) return;

   ConvertMarkerPoints(n, x, y, fPoly);
   DrawPolyMarker();
}

////////////////////////////////////////////////////////////////////////////////
///Poly-marker.

void TGLPadPainter::DrawPolyMarker(Int_t n, const Float_t *x, const Float_t *y)
{
   if (fLocked) return;

   ConvertMarkerPoints(n, x, y, fPoly);
   DrawPolyMarker();
}

////////////////////////////////////////////////////////////////////////////////
///Poly-marker.

void TGLPadPainter::DrawPolyMarker()
{
   if (fLocked) return;

   SaveProjectionMatrix();
   glLoadIdentity();
   //
   glOrtho(0, gPad->GetAbsWNDC() * gPad->GetWw(), 0, gPad->GetAbsHNDC() * gPad->GetWh(), -10., 10.);
   //
   glMatrixMode(GL_MODELVIEW);
   //
   const TGLEnableGuard blendGuard(GL_BLEND);

   Float_t rgba[4] = {};
   Rgl::Pad::ExtractRGBA(gVirtualX->GetMarkerColor(), rgba);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   glColor4fv(rgba);

   const Width_t w = TMath::Max(1, Int_t(TAttMarker::GetMarkerLineWidth(gVirtualX->GetMarkerStyle())));
   glLineWidth(w > fLimits.GetMaxLineWidth() ? fLimits.GetMaxLineWidth() : !w ? 1.f : w);

   const TPoint *xy = &fPoly[0];
   const Style_t markerStyle = TAttMarker::GetMarkerStyleBase(gVirtualX->GetMarkerStyle());
   const UInt_t n = UInt_t(fPoly.size());
   switch (markerStyle) {
   case kDot:
      fMarker.DrawDot(n, xy);
      break;
   case kPlus:
      fMarker.DrawPlus(n, xy);
      break;
   case kStar:
   case 31:
      fMarker.DrawStar(n, xy);
      break;
   case kCircle:
   case kOpenCircle:
      fMarker.DrawCircle(n, xy);
      break;
   case kMultiply:
      fMarker.DrawX(n, xy);
      break;
   case kFullDotSmall://"Full dot small"
      fMarker.DrawFullDotSmall(n, xy);
      break;
   case kFullDotMedium:
      fMarker.DrawFullDotMedium(n, xy);
      break;
   case kFullDotLarge:
   case kFullCircle:
      fMarker.DrawFullDotLarge(n, xy);
      break;
   case kFullSquare:
      fMarker.DrawFullSquare(n, xy);
      break;
   case kFullTriangleUp:
      fMarker.DrawFullTrianlgeUp(n, xy);
      break;
   case kFullTriangleDown:
      fMarker.DrawFullTrianlgeDown(n, xy);
      break;
   case kOpenSquare:
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      fMarker.DrawFullSquare(n, xy);
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      break;
   case kOpenTriangleUp:
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      fMarker.DrawFullTrianlgeUp(n, xy);
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      break;
   case kOpenDiamond:
      fMarker.DrawDiamond(n, xy);
      break;
   case kOpenCross:
      fMarker.DrawOpenCross(n, xy);
      break;
   case kFullStar:
      fMarker.DrawFullStar(n, xy);
      break;
   case kOpenStar:
      fMarker.DrawOpenStar(n, xy);
      break;
   case kOpenTriangleDown:
      fMarker.DrawOpenTrianlgeDown(n, xy);
      break;
   case kFullDiamond:
      fMarker.DrawFullDiamond(n, xy);
      break;
   case kFullCross:
      fMarker.DrawFullCross(n, xy);
      break;
   case kOpenDiamondCross:
      fMarker.DrawOpenDiamondCross(n, xy);
      break;
   case kOpenSquareDiagonal:
      fMarker.DrawOpenSquareDiagonal(n, xy);
      break;
   case kOpenThreeTriangles:
      fMarker.DrawOpenThreeTriangles(n, xy);
      break;
   case kOctagonCross:
      fMarker.DrawOctagonCross(n, xy);
      break;
   case kFullThreeTriangles:
      fMarker.DrawFullThreeTriangles(n, xy);
      break;
   case kOpenFourTrianglesX:
      fMarker.DrawOpenFourTrianglesX(n, xy);
      break;
   case kFullFourTrianglesX:
      fMarker.DrawFullFourTrianglesX(n, xy);
      break;
   case kOpenDoubleDiamond:
      fMarker.DrawOpenDoubleDiamond(n, xy);
      break;
   case kFullDoubleDiamond:
      fMarker.DrawFullDoubleDiamond(n, xy);
      break;
   case kOpenFourTrianglesPlus:
      fMarker.DrawOpenFourTrianglesPlus(n, xy);
      break;
   case kFullFourTrianglesPlus:
      fMarker.DrawFullFourTrianglesPlus(n, xy);
      break;
   case kOpenCrossX:
      fMarker.DrawOpenCrossX(n, xy);
      break;
   case kFullCrossX:
      fMarker.DrawFullCrossX(n, xy);
      break;
   case kFourSquaresX:
      fMarker.DrawFourSquaresX(n, xy);
      break;
   case kFourSquaresPlus:
      fMarker.DrawFourSquaresPlus(n, xy);
      break;
   }

   RestoreProjectionMatrix();
   glMatrixMode(GL_MODELVIEW);
   glLineWidth(1.f);
}

////////////////////////////////////////////////////////////////////////////////

template<class Char>
void TGLPadPainter::DrawTextHelper(Double_t x, Double_t y, const Char *text, ETextMode /*mode*/)
{
   SaveProjectionMatrix();

   glLoadIdentity();
   //
   glOrtho(0, gPad->GetAbsWNDC() * gPad->GetWw(), 0, gPad->GetAbsHNDC() * gPad->GetWh(), -10., 10.);
   //
   glMatrixMode(GL_MODELVIEW);

   Float_t rgba[4] = {};
   Rgl::Pad::ExtractRGBA(gVirtualX->GetTextColor(), rgba);
   glColor4fv(rgba);

   //10 is the first valid font index.
   //20 is FreeSerifBold, as in TTF.cxx and in TGLFontManager.cxx.
   //shift - is the shift to access "extended" fonts.
   const Int_t shift = TGLFontManager::GetExtendedFontStartIndex();

   Int_t fontIndex = TMath::Max(Short_t(10), gVirtualX->GetTextFont());
   if (fontIndex / 10 + shift > TGLFontManager::GetFontFileArray()->GetEntries())
      fontIndex = 20 + shift * 10;
   else
      fontIndex += shift * 10;

   fFM.RegisterFont(TMath::Max(Int_t(gVirtualX->GetTextSize()) - 1, 10),//kTexture does not work if size < 10.
                               TGLFontManager::GetFontNameFromId(fontIndex),
                               TGLFont::kTexture, fF);
   fF.PreRender();

   const UInt_t padH = UInt_t(gPad->GetAbsHNDC() * gPad->GetWh());
   fF.Render(text, gPad->XtoPixel(x), padH - gPad->YtoPixel(y), GetTextAngle(), GetTextMagnitude());

   fF.PostRender();
   RestoreProjectionMatrix();

   glMatrixMode(GL_MODELVIEW);
}

////////////////////////////////////////////////////////////////////////////////
///Draw text. This operation is especially
///dangerous if in locked state -
///ftgl will assert on zero texture size
///(which is result of bad GL context).

void TGLPadPainter::DrawText(Double_t x, Double_t y, const char *text, ETextMode mode)
{
   if (fLocked) return;

   if (!gVirtualX->GetTextSize())
      return;

   DrawTextHelper(x, y, text, mode);
}

////////////////////////////////////////////////////////////////////////////////
///Draw text. This operation is especially
///dangerous if in locked state -
///ftgl will assert on zero texture size
///(which is result of bad GL context).

void TGLPadPainter::DrawText(Double_t x, Double_t y, const wchar_t *text, ETextMode mode)
{
   if (fLocked) return;

   if (!gVirtualX->GetTextSize())
      return;

   DrawTextHelper(x, y, text, mode);
}

////////////////////////////////////////////////////////////////////////////////
///Draw text in NDC. This operation is especially
///dangerous if in locked state -
///ftgl will assert on zero texture size
///(which is result of bad GL context).

void TGLPadPainter::DrawTextNDC(Double_t u, Double_t v, const char *text, ETextMode mode)
{
   if (fLocked) return;

   const Double_t xRange = gPad->GetX2() - gPad->GetX1();
   const Double_t yRange = gPad->GetY2() - gPad->GetY1();
   DrawText(gPad->GetX1() + u * xRange, gPad->GetY1() + v * yRange, text, mode);
}

////////////////////////////////////////////////////////////////////////////////
///Draw text in NDC. This operation is especially
///dangerous if in locked state -
///ftgl will assert on zero texture size
///(which is result of bad GL context).

void TGLPadPainter::DrawTextNDC(Double_t u, Double_t v, const wchar_t *text, ETextMode mode)
{
   if (fLocked) return;

   const Double_t xRange = gPad->GetX2() - gPad->GetX1();
   const Double_t yRange = gPad->GetY2() - gPad->GetY1();
   DrawText(gPad->GetX1() + u * xRange, gPad->GetY1() + v * yRange, text, mode);
}

////////////////////////////////////////////////////////////////////////////////
///Save the projection matrix.
///Attention! GL_PROJECTION will become the current matrix
///after this call!

void TGLPadPainter::SaveProjectionMatrix()const
{
   glMatrixMode(GL_PROJECTION);
   glPushMatrix();
}

////////////////////////////////////////////////////////////////////////////////
///Restore the projection matrix.
///Attention! GL_PROJECTION will become the current matrix
///after this call!

void TGLPadPainter::RestoreProjectionMatrix()const
{
   glMatrixMode(GL_PROJECTION);
   glPopMatrix();
}

////////////////////////////////////////////////////////////////////////////////
///Save the modelview matrix.
///Attention! GL_MODELVIEW will become the current matrix
///after this call!

void TGLPadPainter::SaveModelviewMatrix()const
{
   glMatrixMode(GL_MODELVIEW);
   glPushMatrix();
}

////////////////////////////////////////////////////////////////////////////////
///Restore the modelview matrix.
///Attention! GL_MODELVIEW will become the current matrix
///after this call!

void TGLPadPainter::RestoreModelviewMatrix()const
{
   glMatrixMode(GL_MODELVIEW);
   glPopMatrix();
}

////////////////////////////////////////////////////////////////////////////////
///Extract and save the current viewport.

void TGLPadPainter::SaveViewport()
{
   glGetIntegerv(GL_VIEWPORT, fVp);
}

////////////////////////////////////////////////////////////////////////////////
///Restore the saved viewport.

void TGLPadPainter::RestoreViewport()
{
   glViewport(fVp[0], fVp[1], fVp[2], fVp[3]);
}

////////////////////////////////////////////////////////////////////////////////
/// Using TImage save frame-buffer contents as a picture.

void TGLPadPainter::SaveImage(TVirtualPad *pad, const char *fileName, Int_t type) const
{
   TVirtualPad *canvas = (TVirtualPad *)pad->GetCanvas();
   if (!canvas)
      return;

   gROOT->ProcessLine(Form("((TCanvas *)0x%lx)->Flush();", (ULong_t)canvas));

   std::vector<unsigned> buff(canvas->GetWw() * canvas->GetWh());
   glPixelStorei(GL_PACK_ALIGNMENT, 1);
   glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
   //In case GL_BGRA is not in gl.h (old windows' gl) - comment/uncomment lines.
   //glReadPixels(0, 0, canvas->GetWw(), canvas->GetWh(), GL_BGRA, GL_UNSIGNED_BYTE, (char *)&buff[0]);
   glReadPixels(0, 0, canvas->GetWw(), canvas->GetWh(), GL_RGBA, GL_UNSIGNED_BYTE, (char *)&buff[0]);

   std::unique_ptr<TImage> image(TImage::Create());
   if (!image.get()) {
      ::Error("TGLPadPainter::SaveImage", "TImage creation failed");
      return;
   }

   image->DrawRectangle(0, 0, canvas->GetWw(), canvas->GetWh());
   UInt_t *argb = image->GetArgbArray();

   if (!argb) {
      ::Error("TGLPadPainter::SaveImage", "null argb array in TImage object");
      return;
   }

   const Int_t nLines  = canvas->GetWh();
   const Int_t nPixels = canvas->GetWw();

   for (Int_t i = 0; i < nLines; ++i) {
     Int_t base = (nLines - 1 - i) * nPixels;
     for (Int_t j = 0; j < nPixels; ++j, ++base) {
        //Uncomment/comment if you don't have GL_BGRA.

        const UInt_t pix  = buff[base];
        const UInt_t bgra = ((pix & 0xff) << 16) | (pix & 0xff00) |
                            ((pix & 0xff0000) >> 16) | (pix & 0xff000000);

        //argb[i * nPixels + j] = buff[base];
        argb[i * nPixels + j] = bgra;
     }
   }

   image->WriteImage(fileName, (TImage::EImageFileTypes)type);
}

////////////////////////////////////////////////////////////////////////////////

void TGLPadPainter::DrawPixels(const unsigned char *pixelData, UInt_t width, UInt_t height,
                               Int_t dstX, Int_t dstY, Bool_t enableBlending)
{
   if (fLocked)
      return;

   if (!pixelData) {
      //I'd prefer an assert.
      ::Error("TGLPadPainter::DrawPixels", "pixel data is null");
      return;
   }

   if (std::numeric_limits<UInt_t>::digits >= 32) {
      //TASImage uses bit 31 as ...
      //alpha channel flag! FUUUUUUUUUUUUU .....   !!!
      CLRBIT(width, 31);
      CLRBIT(height, 31);
   }

   if (!width) {
      //Assert is better.
      ::Error("TGLPadPainter::DrawPixels", "invalid width");
      return;
   }

   if (!height) {
      //Assert is better.
      ::Error("TGLPadPainter::DrawPixels", "invalid height");
      return;
   }

   if (TPad *pad = dynamic_cast<TPad *>(gPad)) {
      //TASImage passes pixel coordinates in pad's pixmap coordinate space.
      //While glRasterPosX said to work with 'window' coordinates,
      //that's a lie :) it does not :)

      const Double_t rasterX = Double_t(dstX) / (pad->GetAbsWNDC() * pad->GetWw()) *
                                (pad->GetX2() - pad->GetX1()) + pad->GetX1();

      const Double_t yRange = pad->GetY2() - pad->GetY1();
      const Double_t rasterY = yRange - Double_t(dstY + height) / (pad->GetAbsHNDC() * pad->GetWh()) * yRange +
                               pad->GetY1();

      GLdouble oldPos[4] = {};
      //Save the previous raster pos.
      glGetDoublev(GL_CURRENT_RASTER_POSITION, oldPos);

      glRasterPos2d(rasterX, rasterY);
      //Stupid asimage provides us upside-down image.
      std::vector<unsigned char> upsideDownImage(4 * width * height);
      const unsigned char *srcLine = pixelData + 4 * width * (height - 1);
      unsigned char *dstLine = &upsideDownImage[0];
      for (UInt_t i = 0; i < height; ++i, srcLine -= 4 * width, dstLine += 4 * width)
         std::copy(srcLine, srcLine + 4 * width, dstLine);

      if (enableBlending) {
         glEnable(GL_BLEND);
         glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      }

      glDrawPixels(width, height, GL_BGRA, GL_UNSIGNED_BYTE, &upsideDownImage[0]);

      if (enableBlending)
         glDisable(GL_BLEND);

      //Restore raster pos.
      glRasterPos2d(oldPos[0], oldPos[1]);
   } else
      ::Error("TGLPadPainter::DrawPixels", "no pad found to draw");
}

//Aux. functions - gradient and solid fill of arbitrary area.

////////////////////////////////////////////////////////////////////////////////
///At the moment I assume both linear and radial gradients will work the same way  -
///using a stencil buffer and some big rectangle(s) to fill with a gradient.
///Thus I have a 'common' part - the part responsible for a stencil test.

void TGLPadPainter::DrawPolygonWithGradient(Int_t n, const Double_t *x, const Double_t *y)
{
   assert(n > 2 && "DrawPolygonWithGradient, invalid number of points");
   assert(x != 0 && "DrawPolygonWithGradient, parameter 'x' is null");
   assert(y != 0 && "DrawPolygonWithGradient, parameter 'y' is null");

   assert(dynamic_cast<TColorGradient *>(gROOT->GetColor(gVirtualX->GetFillColor())) != 0 &&
          "DrawPolygonWithGradient, the current fill color is not a gradient fill");
   const TColorGradient * const grad =
         dynamic_cast<TColorGradient *>(gROOT->GetColor(gVirtualX->GetFillColor()));

   if (fLocked)
      return;

   //Now, some magic!
   const TGLEnableGuard stencilGuard(GL_STENCIL_TEST);

   //TODO: check that the state is restored back correctly after
   //      we done with a gradient.
   //TODO: make sure that we have glDepthMask set to false in general!
   glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

   glStencilFunc(GL_NEVER, 1, 0xFF);
   glStencilOp(GL_REPLACE, GL_KEEP, GL_KEEP);// draw 1s on test fail (always)
   //Draw stencil pattern
   glStencilMask(0xFF);
   glClear(GL_STENCIL_BUFFER_BIT);

   //Draw our polygon into the stencil buffer:
   DrawTesselation(n, x, y);

   glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
   glStencilMask(0x00);
   //Draw where stencil's value is 0
   glStencilFunc(GL_EQUAL, 0, 0xFF);
   //Draw only where stencil's value is 1
   glStencilFunc(GL_EQUAL, 1, 0xFF);

   //At the moment radial gradient is derived from linear - it was convenient
   //at some point, but in fact it was a bad idea. And now I have to
   //first check radial gradient.
   //TODO: TRadialGradient must inherit TColorGradient directly.
   const TRadialGradient * const rGrad = dynamic_cast<const TRadialGradient *>(grad);
   if (rGrad)
      DrawGradient(rGrad, n, x, y);
   else {
      const TLinearGradient * const lGrad = dynamic_cast<const TLinearGradient *>(grad);
      assert(lGrad != 0 && "DrawPolygonWithGradient, unknown gradient type");
      DrawGradient(lGrad, n, x, y);
   }
}

////////////////////////////////////////////////////////////////////////////////

void TGLPadPainter::DrawGradient(const TRadialGradient *grad, Int_t nPoints,
                                 const Double_t *xs, const Double_t *ys)
{
   assert(grad != 0 && "DrawGradient, parameter 'grad' is null");
   assert(nPoints > 2 && "DrawGradient, invalid number of points");
   assert(xs != 0 && "DrawGradient, parameter 'xs' is null");
   assert(ys != 0 && "DrawGradient, parameter 'ys' is null");

   if (grad->GetGradientType() != TRadialGradient::kSimple) {
      ::Warning("TGLPadPainter::DrawGradient",
                "extended radial gradient is not supported");//yet?
      return;
   }

   //TODO: check the polygon's bbox!
   const auto &bbox = Rgl::Pad::FindBoundingRect(nPoints, xs, ys);
   //
   auto center = grad->GetCenter();
   auto radius = grad->GetRadius();
   //Adjust the center and radius depending on coordinate mode.
   if (grad->GetCoordinateMode() == TColorGradient::kObjectBoundingMode) {
      radius *= TMath::Max(bbox.fWidth, bbox.fHeight);
      center.fX = bbox.fWidth * center.fX + bbox.fXMin;
      center.fY = bbox.fHeight * center.fY + bbox.fYMin;
   } else {
      const auto w = gPad->GetX2() - gPad->GetX1();
      const auto h = gPad->GetY2() - gPad->GetY1();

      radius *= TMath::Max(w, h);
      center.fX *= w;
      center.fY *= h;
   }
   //Now for the gradient fill we switch into pixel coordinates:
   const auto pixelW = gPad->GetAbsWNDC() * gPad->GetWw();
   const auto pixelH = gPad->GetAbsHNDC() * gPad->GetWh();
   //
   SaveProjectionMatrix();
   SaveModelviewMatrix();
   //A new ortho projection:
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   //
   glOrtho(0., pixelW, 0., pixelH, -10., 10.);
   //
   radius *= TMath::Max(pixelH, pixelW);
   center.fX = gPad->XtoPixel(center.fX);
   center.fY = pixelH - gPad->YtoPixel(center.fY);

   Double_t maxR = 0.;
   {
   const Double_t xMin = gPad->XtoPixel(bbox.fXMin);
   const Double_t xMax = gPad->XtoPixel(bbox.fXMax);
   const Double_t yMin = pixelH - gPad->YtoPixel(bbox.fYMin);
   const Double_t yMax = pixelH - gPad->YtoPixel(bbox.fYMax);
   //Get the longest distance from the center to the bounding box vertices
   //(this will be the maximum possible radius):
   const Double_t maxDistX = TMath::Max(TMath::Abs(center.fX - xMin),
                                        TMath::Abs(center.fX - xMax));
   const Double_t maxDistY = TMath::Max(TMath::Abs(center.fY - yMin),
                                        TMath::Abs(center.fY - yMax));
   maxR = TMath::Sqrt(maxDistX * maxDistX + maxDistY * maxDistY);
   }

   //If gradient 'stops inside the polygon', we use
   //the solid fill for the area outside of radial gradient:
   const Bool_t solidFillAfter = maxR > radius;
   //We emulate a radial gradient using triangles and linear gradient:
   //TODO: Can be something smarter? (btw even 100 seems to be enough)
   const UInt_t nSlices = 500;

   const auto nColors = grad->GetNumberOfSteps();
   //+1 - the strip from the last color's position to radius,
   //and (probably) + 1 for solidFillAfter.
   const auto nCircles = nColors + 1 + solidFillAfter;

   //TODO: can locations be outside of [0., 1.] ???
   //at the moment I assume the answer is NO, NEVER.
   const auto locations = grad->GetColorPositions();
   // * 2 below == x,y
   std::vector<Double_t> circles(nSlices * nCircles * 2);
   const Double_t angle = TMath::TwoPi() / nSlices;

   //"Main" circles (for colors at locations[i]).
   for (UInt_t i = 0; i < nColors; ++i) {
      const auto circle = &circles[i * nSlices * 2];
      //TODO: either check locations here or somewhere else.
      const auto r = radius * locations[i];
      for (UInt_t j = 0, e = nSlices * 2 - 2; j < e; j += 2) {
         circle[j] = center.fX + r * TMath::Cos(angle * j);
         circle[j + 1] = center.fY + r * TMath::Sin(angle * j);
      }
      //The "closing" vertices:
      circle[(nSlices - 1) * 2] = circle[0];
      circle[(nSlices - 1) * 2 + 1] = circle[1];
   }

   {
   //The strip between lastPos and radius:
   const auto circle = &circles[nColors * nSlices * 2];
   for (UInt_t j = 0, e = nSlices * 2 - 2; j < e; j += 2) {
      circle[j] = center.fX + radius * TMath::Cos(angle * j);
      circle[j + 1] = center.fY + radius * TMath::Sin(angle * j);
   }

   circle[(nSlices - 1) * 2] = circle[0];
   circle[(nSlices - 1) * 2 + 1] = circle[1];
   }

   if (solidFillAfter) {
      //The strip after the radius:
      const auto circle = &circles[(nCircles - 1) * nSlices * 2];
      for (UInt_t j = 0, e = nSlices * 2 - 2; j < e; j += 2) {
         circle[j] = center.fX + maxR * TMath::Cos(angle * j);
         circle[j + 1] = center.fY + maxR * TMath::Sin(angle * j);
      }

      circle[(nSlices - 1) * 2] = circle[0];
      circle[(nSlices - 1) * 2 + 1] = circle[1];
   }

   //Now we draw:
   //1) triangle fan in the center (from center to the locations[1],
   //   with a solid fill).
   //2) quad strips for colors.
   //3) additional quad strip from the lastLocation to the radius
   //4) additional quad strip (if any) from the radius to maxR.

   //RGBA values:
   const auto rgba = grad->GetColors();

   const TGLEnableGuard alphaGuard(GL_BLEND);
   //TODO?
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   //Probably a degenerated case. Maybe not.
   glBegin(GL_TRIANGLE_FAN);
   glColor4dv(rgba);
   glVertex2d(center.fX, center.fY);

   for (UInt_t i = 0, e = nSlices * 2; i < e; i += 2)
      glVertex2dv(&circles[i]);

   glEnd();

   //No auto for circles, explicit types to have const Double_t * const, not Duble_t * const.
   for (UInt_t i = 0; i < nColors - 1; ++i) {
      const Double_t * const inner = &circles[i * nSlices * 2];
      const auto innerRGBA = rgba + i * 4;
      const auto outerRGBA = rgba + (i + 1) * 4;
      const Double_t * const outer = &circles[(i + 1) * nSlices * 2];

      Rgl::DrawQuadStripWithRadialGradientFill(nSlices, inner, innerRGBA, outer, outerRGBA);
   }

   //Probably degenerated strip.
   {
   glBegin(GL_QUAD_STRIP);
   const Double_t * const inner = &circles[nSlices * (nColors - 1) * 2];
   const auto solidRGBA = rgba + (nColors - 1) * 4;
   const Double_t * const outer = &circles[nSlices * nColors * 2];

   Rgl::DrawQuadStripWithRadialGradientFill(nSlices, inner, solidRGBA, outer, solidRGBA);
   }

   if (solidFillAfter) {
      glBegin(GL_QUAD_STRIP);
      const Double_t * const inner = &circles[nSlices * nColors * 2];
      const auto solidRGBA = rgba + (nColors - 1) * 4;
      const Double_t * const outer = &circles[nSlices * (nColors + 1) * 2];

      Rgl::DrawQuadStripWithRadialGradientFill(nSlices, inner, solidRGBA, outer, solidRGBA);
   }

   RestoreProjectionMatrix();
   RestoreModelviewMatrix();
}

////////////////////////////////////////////////////////////////////////////////

void TGLPadPainter::DrawGradient(const TLinearGradient *grad, Int_t n,
                                 const Double_t *x, const Double_t *y)
{
   assert(grad != 0 && "DrawGradient, parameter 'grad' is null");
   assert(n > 2 && "DrawGradient, invalid number of points");
   assert(x != 0 && "DrawGradient, parameter 'x' is null");
   assert(y != 0 && "DrawGradient, parameter 'y' is null");

   //Now we fill the whole scene with one big rectangle
   //(group of rectangles) with a gradient fill using
   //stencil test.

   //Find a bounding rect.
   const auto &bbox = Rgl::Pad::FindBoundingRect(n, x, y);
   //TODO: check the bbox??

   //For the gradient fill we switch into the
   //pixel coordinates.
   SaveProjectionMatrix();
   SaveModelviewMatrix();

   //A new ortho projection:
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();

   const Double_t pixelW = gPad->GetAbsWNDC() * gPad->GetWw();
   const Double_t pixelH = gPad->GetAbsHNDC() * gPad->GetWh();
   glOrtho(0., pixelW, 0., pixelH, -10., 10.);

   //A new modelview:
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   //
   TColorGradient::Point start = grad->GetStart();
   TColorGradient::Point end   = grad->GetEnd();

   //Change gradient coordinates from 'NDC' to pad coords:
   if (grad->GetCoordinateMode() == TColorGradient::kPadMode)
   {
      const Double_t w = gPad->GetX2() - gPad->GetX1();
      const Double_t h = gPad->GetY2() - gPad->GetY1();

      start.fX = start.fX * w;
      start.fY = start.fY * h;
      end.fX   = end.fX * w;
      end.fY   = end.fY * h;
   } else {
      start.fX = start.fX * bbox.fWidth + bbox.fXMin;
      start.fY = start.fY * bbox.fHeight + bbox.fYMin;
      end.fX   = end.fX * bbox.fWidth + bbox.fXMin;
      end.fY   = end.fY * bbox.fHeight + bbox.fYMin;
   }

   //TODO: with a radial fill we'll have to extract the code
   //      below into the separate function/and have additional function
   //      for a radial gradient.
   //Now from pad to pixels:
   start.fX = gPad->XtoPixel(start.fX);
   start.fY = pixelH - gPad->YtoPixel(start.fY);
   end.fX = gPad->XtoPixel(end.fX);
   end.fY = pixelH - gPad->YtoPixel(end.fY);
   const Double_t xMin = gPad->XtoPixel(bbox.fXMin);
   const Double_t xMax = gPad->XtoPixel(bbox.fXMax);
   const Double_t yMin = pixelH - gPad->YtoPixel(bbox.fYMin);
   const Double_t yMax = pixelH - gPad->YtoPixel(bbox.fYMax);
   //

   //TODO: check all calculations!

   //Get the longest distance from the start point to the bounding box vertices:
   const Double_t maxDistX = TMath::Max(TMath::Abs(start.fX - xMin), TMath::Abs(start.fX - xMax));
   const Double_t maxDistY = TMath::Max(TMath::Abs(start.fY - yMin), TMath::Abs(start.fY - yMax));

   const Double_t startEndLength = TMath::Sqrt((end.fX - start.fX) * (end.fX - start.fX) +
                                               (end.fY - start.fY) * (end.fY - start.fY));
   const Double_t h = TMath::Max(TMath::Sqrt(maxDistX * maxDistX + maxDistY * maxDistY),
                                 startEndLength);

   //Boxes with a gradients to emulate gradient fill with many colors:
   const Double_t * const colorPositions = grad->GetColorPositions();
   std::vector<Double_t> gradBoxes(grad->GetNumberOfSteps() + 2);
   gradBoxes[0] = start.fY - h;
   for (unsigned i = 1; i <= grad->GetNumberOfSteps(); ++i)
      gradBoxes[i] = startEndLength * colorPositions[i - 1] + start.fY;

   gradBoxes[grad->GetNumberOfSteps() + 1] = start.fY + h;

   //Rotation angle - gradient's axis:
   Double_t angle = TMath::ACos((startEndLength * (end.fY - start.fY)) /
                                (startEndLength * startEndLength)) * TMath::RadToDeg();
   if (end.fX > start.fX)
      angle *= -1;

   glTranslated(start.fX, start.fY, 0.);
   glRotated(angle, 0., 0., 1.);
   glTranslated(-start.fX, -start.fY, 0.);
   //
   const Double_t * const rgba = grad->GetColors();

   const unsigned nEdges = gradBoxes.size();
   const unsigned nColors = grad->GetNumberOfSteps();
   const Double_t xLeft = start.fX - h, xRight = start.fX + h;

   const TGLEnableGuard alphaGuard(GL_BLEND);
   //TODO?
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   Rgl::DrawBoxWithGradientFill(gradBoxes[0], gradBoxes[1], xLeft, xRight, rgba, rgba);
   Rgl::DrawBoxWithGradientFill(gradBoxes[nEdges - 2], gradBoxes[nEdges - 1], xLeft, xRight,
                           rgba + (nColors - 1) * 4, rgba + (nColors - 1) * 4);

   for (unsigned i = 1; i < nEdges - 2; ++i)
      Rgl::DrawBoxWithGradientFill(gradBoxes[i], gradBoxes[i + 1], xLeft,
                                   xRight, rgba + (i - 1) * 4, rgba + i * 4);

   RestoreProjectionMatrix();
   RestoreModelviewMatrix();
}

////////////////////////////////////////////////////////////////////////////////

void TGLPadPainter::DrawTesselation(Int_t n, const Double_t *x, const Double_t *y)
{
   assert(n > 2 && "DrawTesselation, invalid number of points");
   assert(x != 0 && "DrawTesselation, parameter 'x' is null");
   assert(y != 0 && "DrawTesselation, parameter 'y' is null");

   //Data for a tesselator:
   fVs.resize(n * 3);

   for (Int_t i = 0; i < n; ++i) {
      fVs[i * 3]     = x[i];
      fVs[i * 3 + 1] = y[i];
      fVs[i * 3 + 2] = 0.;
   }

   //TODO: A very primitive way to tesselate - check what
   //kind of polygons we can really have from TPad/TCanvas.
   GLUtesselator *t = (GLUtesselator *)fTess.GetTess();
   gluBeginPolygon(t);
   gluNextContour(t, (GLenum)GLU_UNKNOWN);

   for (Int_t i = 0; i < n; ++i)
      gluTessVertex(t, &fVs[i * 3], &fVs[i * 3]);

   gluEndPolygon(t);
}


//Aux. functions.
namespace {

template<class ValueType>
void ConvertMarkerPoints(Int_t n, const ValueType *x, const ValueType *y, std::vector<TPoint> & dst)
{
   const UInt_t padH = UInt_t(gPad->GetAbsHNDC() * gPad->GetWh());

   dst.resize(n);
   for (Int_t i = 0; i < n; ++i) {
      dst[i].fX = gPad->XtoPixel(x[i]);
      dst[i].fY = padH - gPad->YtoPixel(y[i]);
   }
}

}

