// @(#)root/graf2d:$Id$
// Author: Olivier Couet, Timur Pocheptsov 23/01/2012


/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//#define NDEBUG

#include <stdexcept>
#include <iostream>
#include <cstring>
#include <cassert>
#include <limits>

#include <Cocoa/Cocoa.h>

#  include <ft2build.h>
#  include FT_FREETYPE_H
#  include FT_GLYPH_H

#include "QuartzFillArea.h"
#include "TColorGradient.h"
#include "QuartzMarker.h"
#include "CocoaPrivate.h"
#include "QuartzWindow.h"
#include "QuartzPixmap.h"
#include "QuartzUtils.h"
#include "X11Drawable.h"
#include "QuartzText.h"
#include "QuartzLine.h"
#include "CocoaUtils.h"
#include "TGQuartz.h"
#include "TString.h"
#include "TPoint.h"
#include "TColor.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TEnv.h"
#include "TMath.h"

// To scale fonts to the same size as the TTF version
const Float_t kScale = 0.93376068;

ClassImp(TGQuartz)

namespace X11 = ROOT::MacOSX::X11;
namespace Quartz = ROOT::Quartz;
namespace Util = ROOT::MacOSX::Util;

namespace {

//______________________________________________________________________________
void ConvertPointsROOTToCocoa(Int_t nPoints, const TPoint *xy, std::vector<TPoint> &dst,
                              NSObject<X11Drawable> *drawable)
{
   assert(nPoints != 0 && "ConvertPointsROOTToCocoa, nPoints parameter is 0");
   assert(xy != 0 && "ConvertPointsROOTToCocoa, xy parameter is null");
   assert(drawable != 0 && "ConvertPointsROOTToCocoa, drawable parameter is null");

   const auto scaleFactor = drawable.fScaleFactor;

   dst.resize(nPoints);
   for (Int_t i = 0; i < nPoints; ++i) {
      dst[i].fX = SCoord_t(xy[i].fX * scaleFactor);
      dst[i].fY = SCoord_t(X11::LocalYROOTToCocoa(drawable, xy[i].fY) * scaleFactor);
   }
}

}

//______________________________________________________________________________
TGQuartz::TGQuartz()
            : fUseAA(true), fUseFAAA(false)
{
   //Default ctor.


   if (!TTF::IsInitialized())
      TTF::Init();

   //I do not know why TTF::Init returns void and I have to check IsInitialized() again.
   if (!TTF::IsInitialized())
      Error("TGQuartz", "TTF::Init() failed");

   fAlign.x = 0;
   fAlign.y = 0;

   SetAA();
}


//______________________________________________________________________________
TGQuartz::TGQuartz(const char *name, const char *title)
            : TGCocoa(name, title),
              fUseAA(true), fUseFAAA(false)
{
   //Constructor.
   if (!TTF::IsInitialized())
      TTF::Init();

   //I do not know why TTF::Init returns void and I have to check IsInitialized() again.
   if (!TTF::IsInitialized())
      Error("TGQuartz", "TTF::Init() failed");

   fAlign.x = 0;
   fAlign.y = 0;

   SetAA();
}


//______________________________________________________________________________
void TGQuartz::DrawBox(Int_t x1, Int_t y1, Int_t x2, Int_t y2, EBoxMode mode)
{
   //Check some conditions first.
   if (fDirectDraw) {
      if (!fPimpl->GetDrawable(fSelectedDrawable).fIsPixmap)
         fPimpl->fX11CommandBuffer.AddDrawBoxXor(fSelectedDrawable, x1, y1, x2, y2);
      return;
   }

   NSObject<X11Drawable> * const drawable = (NSObject<X11Drawable> *)GetSelectedDrawableChecked("DrawBox");
   if (!drawable)
      return;

   CGContextRef ctx = drawable.fContext;
   const Quartz::CGStateGuard ctxGuard(ctx);
   //AA flag is not a part of a state.
   const Quartz::CGAAStateGuard aaCtxGuard(ctx, fUseAA);

   //Go to low-left-corner system.
   y1 = Int_t(X11::LocalYROOTToCocoa(drawable, y1));
   y2 = Int_t(X11::LocalYROOTToCocoa(drawable, y2));

   if (const TColorGradient * const gradient = dynamic_cast<TColorGradient *>(gROOT->GetColor(GetFillColor()))) {
      //Draw a box with a gradient fill and a shadow.
      //Ignore all fill styles and EBoxMode, use a gradient fill.
      TPoint polygon[4];
      polygon[0].fX = x1, polygon[0].fY = y1;
      polygon[1].fX = x2, polygon[1].fY = y1;
      polygon[2].fX = x2, polygon[2].fY = y2;
      polygon[3].fX = x1, polygon[3].fY = y2;

      Quartz::DrawPolygonWithGradientFill(ctx, gradient, CGSizeMake(drawable.fWidth, drawable.fHeight),
                                          4, polygon, kFALSE); //kFALSE == don't draw a shadow.
   } else {
      const bool isHollow = mode == kHollow || GetFillStyle() / 1000 == 2;

      //Note! Pattern index (and its address) MUST live
      //long enough to be valid at the point of Quartz::DrawBox call!
      unsigned patternIndex = 0;
      if (isHollow) {
         if (!Quartz::SetLineColor(ctx, GetLineColor())) {
            Error("DrawBox", "Can not find color for index %d", int(GetLineColor()));
            return;
         }
      } else {
         if (!Quartz::SetFillAreaParameters(ctx, &patternIndex)) {
            Error("DrawBox", "SetFillAreaParameters failed");
            return;
         }
      }
      Quartz::SetLineStyle(ctx, GetLineStyle());
      Quartz::SetLineWidth(ctx, GetLineWidth());
      Quartz::DrawBox(ctx, x1, y1, x2, y2, isHollow);
   }
}


//______________________________________________________________________________
void TGQuartz::DrawFillArea(Int_t n, TPoint *xy)
{
   //Comment from TVirtualX:

   // Draw a filled area through all points.
   // n         : number of points
   // xy        : array of points

   //End of comment.
   if (n < 3)
      return;

   //Do some checks first.
   if (fDirectDraw)//To avoid warnings from Quartz - no context at the moment!
      return;

   NSObject<X11Drawable> * const drawable =
               (NSObject<X11Drawable> *)GetSelectedDrawableChecked("DrawFillArea");

   if (!drawable)
      return;

   CGContextRef ctx = drawable.fContext;

   //Convert points to bottom-left system:
   ConvertPointsROOTToCocoa(n, xy, fConvertedPoints, drawable);

   const Quartz::CGStateGuard ctxGuard(ctx);
   //AA flag is not a part of a state.
   const Quartz::CGAAStateGuard aaCtxGuard(ctx, fUseFAAA);

   if (drawable.fScaleFactor > 1.) {
      // The CTM will be restored by 'ctxGuard'.
      CGContextScaleCTM(ctx, 1. / drawable.fScaleFactor, 1. / drawable.fScaleFactor);
   }

   const TColor * const fillColor = gROOT->GetColor(GetFillColor());
   if (!fillColor) {
      Error("DrawFillArea", "Could not find TColor for index %d", GetFillColor());
      return;
   }

   if (const TColorGradient * const gradient = dynamic_cast<const TColorGradient *>(fillColor)) {
      Quartz::DrawPolygonWithGradientFill(ctx, gradient, CGSizeMake(drawable.fWidth, drawable.fHeight),
                                          n, &fConvertedPoints[0], kFALSE);//kFALSE == don't draw a shadow.
   } else {
      unsigned patternIndex = 0;
      if (!Quartz::SetFillAreaParameters(ctx, &patternIndex)) {
         Error("DrawFillArea", "SetFillAreaParameters failed");
         return;
      }

      Quartz::DrawFillArea(ctx, n, &fConvertedPoints[0], kFALSE);//The last argument - do not draw shadows.
   }
}


//______________________________________________________________________________
void TGQuartz::DrawCellArray(Int_t /*x1*/, Int_t /*y1*/, Int_t /*x2*/, Int_t /*y2*/,
                             Int_t /*nx*/, Int_t /*ny*/, Int_t */*ic*/)
{
   //Noop.
}


//______________________________________________________________________________
void TGQuartz::DrawLine(Int_t x1, Int_t y1, Int_t x2, Int_t y2)
{
   // Draw a line.
   // x1,y1        : begin of line
   // x2,y2        : end of line

   if (fDirectDraw) {
      if (!fPimpl->GetDrawable(fSelectedDrawable).fIsPixmap) {
         QuartzView * const view = (QuartzView *)fPimpl->GetWindow(fSelectedDrawable).fContentView;
         if (!view) {
             ::Warning("DrawLine", "Invalid view/window for XOR-mode");
             return;
         }

         [view.fQuartzWindow addCrosshairWindow];
         fPimpl->fX11CommandBuffer.AddDrawLineXor(fSelectedDrawable, x1, y1, x2, y2);
      }

      return;
   }

   //Do some checks first:
   assert(fSelectedDrawable > fPimpl->GetRootWindowID() && "DrawLine, bad drawable is selected");
   NSObject<X11Drawable> * const drawable =
                     (NSObject<X11Drawable> *)GetSelectedDrawableChecked("DrawLine");
   if (!drawable)
      return;

   CGContextRef ctx = drawable.fContext;
   const Quartz::CGStateGuard ctxGuard(ctx);
   //AA flag is not a part of a state.
   const Quartz::CGAAStateGuard aaCtxGuard(ctx, fUseAA);

   if (!Quartz::SetLineColor(ctx, GetLineColor())) {
      Error("DrawLine", "Could not set line color for index %d", int(GetLineColor()));
      return;
   }

   Quartz::SetLineStyle(ctx, GetLineStyle());
   Quartz::SetLineWidth(ctx, GetLineWidth());

   Quartz::DrawLine(ctx, x1, X11::LocalYROOTToCocoa(drawable, y1), x2,
                    X11::LocalYROOTToCocoa(drawable, y2));
}


//______________________________________________________________________________
void TGQuartz::DrawPolyLine(Int_t n, TPoint *xy)
{
   //Comment from TVirtualX:
   // Draw a line through all points.
   // n         : number of points
   // xy        : list of points
   //End of comment.

   //Some checks first.
   if (fDirectDraw)//To avoid warnings from Quartz - no context at the moment!
      return;

   NSObject<X11Drawable> * const drawable =
                     (NSObject<X11Drawable> *)GetSelectedDrawableChecked("DrawPolyLine");
   if (!drawable)
      return;

   CGContextRef ctx = drawable.fContext;
   const Quartz::CGStateGuard ctxGuard(ctx);
   //AA flag is not a part of a state.
   const Quartz::CGAAStateGuard aaCtxGuard(ctx, fUseAA);

   if (!Quartz::SetLineColor(ctx, GetLineColor())) {
      Error("DrawPolyLine", "Could not find TColor for index %d", GetLineColor());
      return;
   }

   Quartz::SetLineStyle(ctx, GetLineStyle());
   Quartz::SetLineWidth(ctx, GetLineWidth());

   //Convert to bottom-left-corner system.
   ConvertPointsROOTToCocoa(n, xy, fConvertedPoints, drawable);

   if (drawable.fScaleFactor > 1.)
      CGContextScaleCTM(ctx, 1. / drawable.fScaleFactor, 1. / drawable.fScaleFactor);

   Quartz::DrawPolyLine(ctx, n, &fConvertedPoints[0]);

   // CTM (current transformation matrix) is restored by 'ctxGuard's dtor.
}


//______________________________________________________________________________
void TGQuartz::DrawPolyMarker(Int_t n, TPoint *xy)
{
   //Comment from TVirtualX:
   // Draw PolyMarker
   // n         : number of points
   // xy        : list of points
   //End of comment.

   //Do some checks first.
   if (fDirectDraw)//To avoid warnings from Quartz - no context at the moment!
      return;

   NSObject<X11Drawable> * const drawable =
                        (NSObject<X11Drawable> *)GetSelectedDrawableChecked("DrawPolyMarker");
   if (!drawable)
      return;

   CGContextRef ctx = drawable.fContext;
   const Quartz::CGStateGuard ctxGuard(ctx);
   //AA flag is not a part of a state.
   const Quartz::CGAAStateGuard aaCtxGuard(ctx, fUseAA);

   if (!Quartz::SetFillColor(ctx, GetMarkerColor())) {
      Error("DrawPolyMarker", "Could not find TColor for index %d", GetMarkerColor());
      return;
   }

   Quartz::SetLineColor(ctx, GetMarkerColor());//Can not fail (for coverity).
   Quartz::SetLineStyle(ctx, 1);
   Quartz::SetLineWidth(ctx, TMath::Max(1, Int_t(TAttMarker::GetMarkerLineWidth(GetMarkerStyle()))));

   ConvertPointsROOTToCocoa(n, xy, fConvertedPoints, drawable);

   if (drawable.fScaleFactor > 1.)
      CGContextScaleCTM(ctx, 1. / drawable.fScaleFactor, 1. / drawable.fScaleFactor);

   Style_t markerstyle = TAttMarker::GetMarkerStyleBase(GetMarkerStyle());

   // The fast pixel markers need to be treated separately
   if (markerstyle == 1 || markerstyle == 6 || markerstyle == 7) {
       CGContextSetLineJoin(ctx, kCGLineJoinMiter);
       CGContextSetLineCap(ctx, kCGLineCapButt);
   } else {
       CGContextSetLineJoin(ctx, kCGLineJoinRound);
       CGContextSetLineCap(ctx, kCGLineCapRound);
   }

   Float_t MarkerSizeReduced = GetMarkerSize() - TMath::Floor(TAttMarker::GetMarkerLineWidth(GetMarkerStyle())/2.)/4.;
   Quartz::DrawPolyMarker(ctx, n, &fConvertedPoints[0], MarkerSizeReduced * drawable.fScaleFactor, markerstyle);

   CGContextSetLineJoin(ctx, kCGLineJoinMiter);
   CGContextSetLineCap(ctx, kCGLineCapButt);
}


//______________________________________________________________________________
void TGQuartz::DrawText(Int_t x, Int_t y, Float_t /*angle*/, Float_t /*mgn*/,
                        const char *text, ETextMode /*mode*/)
{
   if (fDirectDraw)//To avoid warnings from Quartz - no context at the moment!
      return;

   if (!text || !text[0])//Can this ever happen? TPad::PaintText does not check this.
      return;

   if (GetTextSize()<1.5)//Do not draw anything, or CoreText will create some small (but not of size 0 font).
      return;

   NSObject<X11Drawable> * const drawable =
                     (NSObject<X11Drawable> *)GetSelectedDrawableChecked("DrawText");
   if (!drawable)
      return;

   CGContextRef ctx = drawable.fContext;
   const Quartz::CGStateGuard ctxGuard(ctx);

   //Before any core text drawing operations, reset text matrix.
   CGContextSetTextMatrix(ctx, CGAffineTransformIdentity);

   try {
      if (CTFontRef currentFont = fPimpl->fFontManager.SelectFont(GetTextFont(), kScale*GetTextSize())) {
         const unsigned fontIndex = GetTextFont() / 10;
         if (fontIndex == 12 || fontIndex == 15) {//Greek and math symbols.
            //This is a hack. Correct way is to extract glyphs from symbol.ttf,
            //find correct mapping, place this glyphs. This requires manual layout though (?),
            //and as usually, I have to many things to do, may be, one day I'll fix text rendering also.
            //This hack work only on MacOSX 10.7.3, does not work on iOS and I'm not sure about future/previous
            //versions of MacOSX.
            typedef std::vector<UniChar>::size_type size_type;

            std::vector<UniChar> unichars(std::strlen(text));
            for (size_type i = 0, len = unichars.size(); i < len; ++i)
               unichars[i] = 0xF000 + (unsigned char)text[i];

            Quartz::TextLine ctLine(unichars, currentFont, GetTextColor());
            ctLine.DrawLine(ctx, x, X11::LocalYROOTToCocoa(drawable, y));
         } else {
            const Quartz::TextLine ctLine(text, currentFont, GetTextColor());
            ctLine.DrawLine(ctx, x, X11::LocalYROOTToCocoa(drawable, y));
         }
      }
   } catch (const std::exception &e) {
      Error("DrawText", "Exception from Quartz::TextLine: %s", e.what());
   }
}

//______________________________________________________________________________
void TGQuartz::DrawText(Int_t x, Int_t y, Float_t angle, Float_t /*mgn*/, const wchar_t *text, ETextMode mode)
{
   if (!text || !text[0])
      return;

   if (!TTF::IsInitialized()) {
      Error("DrawText", "wchar_t string to draw, but TTF initialization failed");
      return;
   }

   if (!GetTextSize())//Do not draw anything, or CoreText will create some small (but not of size 0 font).
      return;

   (void)x;
   (void)y;
   (void)angle;
   (void)mode;

   TTF::SetSmoothing(kTRUE);
   TTF::SetRotationMatrix(angle);
   TTF::PrepareString(text);
   TTF::LayoutGlyphs();

   AlignTTFString();
   RenderTTFString(x, y, mode);
}

//______________________________________________________________________________
void TGQuartz::GetTextExtent(UInt_t &w, UInt_t &h, char *text)
{
   // Returns the size of the specified character string "mess".
   //
   // w    - the text width
   // h    - the text height
   // text - the string

   if (!text || !text[0]) {
      w = 0;
      h = 0;
      return;
   }

   if (fPimpl->fFontManager.SelectFont(GetTextFont(), kScale*GetTextSize())) {
      const unsigned fontIndex = GetTextFont() / 10;
      if (fontIndex == 12 || fontIndex == 15) {//Greek and math symbols.
         typedef std::vector<UniChar>::size_type size_type;

         std::vector<UniChar> unichars(std::strlen(text));
         for (size_type i = 0, len = unichars.size(); i < len; ++i)
            unichars[i] = 0xF000 + (unsigned char)text[i];

         fPimpl->fFontManager.GetTextBounds(w, h, unichars);
      } else {
         fPimpl->fFontManager.GetTextBounds(w, h, text);
      }
   }
}

//______________________________________________________________________________
Int_t TGQuartz::GetFontAscent() const
{
   // Returns the ascent of the current font (in pixels).
   // The ascent of a font is the distance from the baseline
   // to the highest position characters extend to.
   if (fPimpl->fFontManager.SelectFont(GetTextFont(), kScale*GetTextSize()))
      return Int_t(fPimpl->fFontManager.GetAscent());

   return 0;
}

//______________________________________________________________________________
Int_t TGQuartz::GetFontAscent(const char *text) const
{
   // Returns the ascent of the current font (in pixels).
   // The ascent of a font is the distance from the baseline
   // to the highest position characters extend to.

   //In case of any problem we can always resort to the old version:
   if (!text || !text[0])//How it's usually tested in ROOT
      return GetFontAscent();

   if (fPimpl->fFontManager.SelectFont(GetTextFont(), kScale*GetTextSize())) {
      const unsigned fontIndex = GetTextFont() / 10;
      if (fontIndex == 12 || fontIndex == 15) {//Greek and math symbols.
         //That's an ugly hack :)
         typedef std::vector<UniChar>::size_type size_type;

         std::vector<UniChar> unichars(std::strlen(text));
         for (size_type i = 0, len = unichars.size(); i < len; ++i)
            unichars[i] = 0xF000 + (unsigned char)text[i];

         return Int_t(fPimpl->fFontManager.GetAscent(unichars));
      } else
         return Int_t(fPimpl->fFontManager.GetAscent(text));
   }

   return 0;
}

//______________________________________________________________________________
Int_t TGQuartz::GetFontDescent() const
{
   // Returns the descent of the current font (in pixels.
   // The descent is the distance from the base line
   // to the lowest point characters extend to.
   if (fPimpl->fFontManager.SelectFont(GetTextFont(), kScale*GetTextSize()))
      return Int_t(fPimpl->fFontManager.GetDescent());

   return 0;
}

//______________________________________________________________________________
Int_t TGQuartz::GetFontDescent(const char *text) const
{
   // Returns the descent of the current font (in pixels.
   // The descent is the distance from the base line
   // to the lowest point characters extend to.

   //That's how it's tested in ROOT:
   if (!text || !text[0])
      return GetFontDescent();

   if (fPimpl->fFontManager.SelectFont(GetTextFont(), kScale*GetTextSize())) {
      const unsigned fontIndex = GetTextFont() / 10;
      if (fontIndex == 12 || fontIndex == 15) {//Greek and math symbols.
         //That's an ugly hack :)
         typedef std::vector<UniChar>::size_type size_type;

         std::vector<UniChar> unichars(std::strlen(text));
         for (size_type i = 0, len = unichars.size(); i < len; ++i)
            unichars[i] = 0xF000 + (unsigned char)text[i];

         return Int_t(fPimpl->fFontManager.GetDescent(unichars));
      } else
         return Int_t(fPimpl->fFontManager.GetDescent(text));
   }

   return 0;
}


//______________________________________________________________________________
Float_t TGQuartz::GetTextMagnitude()
{
   // Returns the current font magnification factor
   return 0;
}

//______________________________________________________________________________
void TGQuartz::SetLineColor(Color_t cindex)
{
   // Set color index "cindex" for drawing lines.
   TAttLine::SetLineColor(cindex);
}


//______________________________________________________________________________
void TGQuartz::SetLineStyle(Style_t lstyle)
{
   // Set line style.
   TAttLine::SetLineStyle(lstyle);
}


//______________________________________________________________________________
void TGQuartz::SetLineWidth(Width_t width)
{
   // Set the line width.

   TAttLine::SetLineWidth(width);
}


//______________________________________________________________________________
void TGQuartz::SetFillColor(Color_t cindex)
{
   // Set color index "cindex" for fill areas.

   TAttFill::SetFillColor(cindex);
}


//______________________________________________________________________________
void TGQuartz::SetFillStyle(Style_t style)
{
   // Set fill area style.
   TAttFill::SetFillStyle(style);
}


//______________________________________________________________________________
void TGQuartz::SetMarkerColor(Color_t cindex)
{
   // Set color index "cindex" for markers.
   TAttMarker::SetMarkerColor(cindex);
}


//______________________________________________________________________________
void TGQuartz::SetMarkerSize(Float_t markersize)
{
   // Set marker size index.
   //
   // markersize - the marker scale factor
   TAttMarker::SetMarkerSize(markersize);
}


//______________________________________________________________________________
void TGQuartz::SetMarkerStyle(Style_t markerstyle)
{
   // Set marker style.

   TAttMarker::SetMarkerStyle(markerstyle);
}


//______________________________________________________________________________
void TGQuartz::SetTextAlign(Short_t talign)
{
   // Set the text alignment.
   //
   // talign = txalh horizontal text alignment
   // talign = txalv vertical text alignment

   TAttText::SetTextAlign(talign);
}

//______________________________________________________________________________
void TGQuartz::SetTextColor(Color_t cindex)
{
   // Set the color index "cindex" for text.

   TAttText::SetTextColor(cindex);
}


//______________________________________________________________________________
void TGQuartz::SetTextFont(Font_t fontNumber)
{
   // Set the current text font number.

   TAttText::SetTextFont(fontNumber);

   if (!TTF::IsInitialized()) {
      Error("SetTextFont", "TTF is not initialized");
      return;
   }

   TTF::SetTextFont(fontNumber);
}

//______________________________________________________________________________
Int_t TGQuartz::SetTextFont(char *fontName, ETextSetMode /*mode*/)
{
   //This function is never used in gPad (in normal text rendering,
   //so I'm not setting anything for CoreText).
   if (!TTF::IsInitialized()) {
      Error("SetTextFont", "TTF is not initialized");
      return 0;
   }

   return TTF::SetTextFont(fontName);
}

//______________________________________________________________________________
void TGQuartz::SetTextSize(Float_t textsize)
{
   // Set the current text size to "textsize"

   TAttText::SetTextSize(textsize);

   if (!TTF::IsInitialized()) {
      Error("SetTextSize", "TTF is not initialized");
      return;
   }

   TTF::SetTextSize(textsize);
}


//______________________________________________________________________________
void TGQuartz::SetOpacity(Int_t /*percent*/)
{
   // Set opacity of the current window. This image manipulation routine
   // works by adding to a percent amount of neutral to each pixels RGB.
   // Since it requires quite some additional color map entries is it
   // only supported on displays with more than > 8 color planes (> 256
   // colors).
}

//TTF related part.

//______________________________________________________________________________
void TGQuartz::AlignTTFString()
{
   //Comment from TGX11TTF:
   // Compute alignment variables. The alignment is done on the horizontal string
   // then the rotation is applied on the alignment variables.
   // SetRotation and LayoutGlyphs should have been called before.
   //End of comment.

   //This code is from TGX11TTF (with my fixes).
   //It looks like align can not be both X and Y align?

   const EAlign align = EAlign(fTextAlign);

   // vertical alignment
   if (align == kTLeft || align == kTCenter || align == kTRight) {
      fAlign.y = TTF::GetAscent();
   } else if (align == kMLeft || align == kMCenter || align == kMRight) {
      fAlign.y = TTF::GetAscent() / 2;
   } else {
      fAlign.y = 0;
   }

   // horizontal alignment
   if (align == kTRight || align == kMRight || align == kBRight) {
      fAlign.x = TTF::GetWidth();
   } else if (align == kTCenter || align == kMCenter || align == kBCenter) {
      fAlign.x = TTF::GetWidth() / 2;
   } else {
      fAlign.x = 0;
   }

   FT_Vector_Transform(&fAlign, TTF::GetRotMatrix());
   //This shift is from the original code.
   fAlign.x = fAlign.x >> 6;
   fAlign.y = fAlign.y >> 6;
}

//______________________________________________________________________________
Bool_t TGQuartz::IsTTFStringVisible(Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   //Comment from TGX11TTF:
   // Test if there is really something to render.
   //End of comment.

   //This code is from TGX11TTF (with modifications).

   //Comment from TGX11TTF:
   // If w or h is 0, very likely the string is only blank characters
   if (!w || !h)
      return kFALSE;

   UInt_t width = 0;
   UInt_t height = 0;
   Int_t xy = 0;

   GetWindowSize(GetCurrentWindow(), xy, xy, width, height);

   // If string falls outside window, there is probably no need to draw it.
   if (x + int(w) <= 0 || x >= int(width))
      return kFALSE;

   if (y + int(h) <= 0 || y >= int(height))
      return kFALSE;

   return kTRUE;
}

//______________________________________________________________________________
void TGQuartz::RenderTTFString(Int_t x, Int_t y, ETextMode mode)
{
   //Comment from TGX11TTF:
   // Perform the string rendering in the pad.
   // LayoutGlyphs should have been called before.
   //End of comment.

   //This code is a modified (for Quartz) version of TG11TTF::RenderString.

   NSObject<X11Drawable> * const drawable = (NSObject<X11Drawable> *)GetSelectedDrawableChecked("DrawText");
   if (!drawable)
      return;

   QuartzPixmap *dstPixmap = nil;
   if ([drawable isKindOfClass : [QuartzPixmap class]])
      dstPixmap = (QuartzPixmap *)drawable;
   else if ([drawable isKindOfClass : [QuartzView class]] || [drawable isKindOfClass : [QuartzWindow class]])
      dstPixmap = ((NSObject<X11Window> *)drawable).fBackBuffer;

   if (!dstPixmap) {
      //I can not read pixels from a window (I can, but this is too slow and unreliable).
      Error("DrawText", "fSelectedDrawable is neither QuartzPixmap nor a double buffered window");
      return;
   }

   //Comment from TGX11TTF:
   // compute the size and position of the XImage that will contain the text
   const Int_t xOff = TTF::GetBox().xMin < 0 ? -TTF::GetBox().xMin : 0;
   const Int_t yOff = TTF::GetBox().yMin < 0 ? -TTF::GetBox().yMin : 0;

   const Int_t w = TTF::GetBox().xMax + xOff;
   const Int_t h = TTF::GetBox().yMax + yOff;

   const Int_t x1 = x - xOff - fAlign.x;
   const Int_t y1 = y + yOff + fAlign.y - h;

   if (!IsTTFStringVisible(x1, y1, w, h))
      return;

   //By default, all pixels are set to 0 (all components, that's what code in TGX11TTF also does here).
   Util::NSScopeGuard<QuartzPixmap> pixmap([[QuartzPixmap alloc] initWithW : w H : h scaleFactor : 1.f]);
   if (!pixmap.Get()) {
      Error("DrawText", "pixmap creation failed");
      return;
   }

   const unsigned char defaultBackgroundPixel[] = {255, 255, 255, 255};
   Util::ScopedArray<unsigned char> arrayGuard;
   if (mode == kClear) {
      //For this mode, TGX11TTF does some work to: a) preserve pixels under symbols
      //b) calculate (interpolate) pixel for glyphs.

      X11::Rectangle bbox(x1, y1, w, h);
      //We already check IsVisible, so, in principle, bbox at least has intersection with
      //the current selected drawable.
      if (X11::AdjustCropArea(dstPixmap, bbox))
         arrayGuard.Reset([dstPixmap readColorBits : bbox]);

      if (!arrayGuard.Get()) {
         Error("DrawText", "problem with reading background pixels");
         return;
      }

      // This is copy & paste from TGX11TTF:
      const Int_t xo = x1 < 0 ? -x1 : 0;
      const Int_t yo = y1 < 0 ? -y1 : 0;

      for (int yp = 0; yp < int(bbox.fHeight) && yo + yp < h; ++yp) {
         const unsigned char *srcBase = arrayGuard.Get() + bbox.fWidth * yp * 4;
         for (int xp = 0; xp < int(bbox.fWidth) && xo + xp < w; ++xp) {
            const unsigned char * const pixel = srcBase + xp * 4;
            [pixmap.Get() putPixel : pixel X : xo + xp Y : yo + yp];
         }
      }
   } else {
      //Find background color and set for all pixels.
      [pixmap.Get() addPixel : defaultBackgroundPixel];
   }

   CGContextRef ctx = drawable.fContext;
   const Quartz::CGStateGuard ctxGuard(ctx);

   CGContextSetRGBStrokeColor(ctx, 0., 0., 1., 1.);
   // paint the glyphs in the pixmap.
   TTF::TTGlyph *glyph = TTF::GetGlyphs();
   for (int n = 0; n < TTF::GetNumGlyphs(); ++n, ++glyph) {
      if (FT_Glyph_To_Bitmap(&glyph->fImage, TTF::GetSmoothing() ? ft_render_mode_normal : ft_render_mode_mono, 0, 1 ))
         continue;

      FT_BitmapGlyph bitmap = (FT_BitmapGlyph)glyph->fImage;
      FT_Bitmap *source = &bitmap->bitmap;
      const Int_t bx = bitmap->left + xOff;
      const Int_t by = h - bitmap->top - yOff;

      DrawFTGlyphIntoPixmap(pixmap.Get(), source, TGCocoa::GetPixel(GetTextColor()),
                            mode == kClear ? ULong_t(-1) : 0xffffff, bx, by);
   }

   const X11::Rectangle copyArea(0, 0, w, h);
   const X11::Point dstPoint(x1, y1);
   [dstPixmap copy : pixmap.Get() area : copyArea withMask : nil clipOrigin : X11::Point() toPoint : dstPoint];
}

//______________________________________________________________________________
void TGQuartz::DrawFTGlyphIntoPixmap(void *pHack, FT_Bitmap *source, ULong_t fore, ULong_t back, Int_t bx, Int_t by)
{
   //This function is a "remake" of TGX11FFT::DrawImage.

   //I'm using this code to reproduce the same text as generated by TGX11TTF.
   //It's quite sloppy, as in original version. I tried to make it not so ugly and
   //more or less readable.

   QuartzPixmap *pixmap = (QuartzPixmap *)pHack;
   assert(pixmap != nil && "DrawFTGlyphIntoPixmap, pixmap parameter is nil");
   assert(source != 0 && "DrawFTGlyphIntoPixmap, source parameter is null");

   if (TTF::GetSmoothing()) {
      static ColorStruct_t col[5];
      // background kClear, i.e. transparent, we take as background color
      // the average of the rgb values of all pixels covered by this character
      if (back == ULong_t(-1) && source->width) {
         const int maxDots = 50000;
         int dots = Int_t(source->width * source->rows);
         if (dots > maxDots)
            dots = maxDots;

         //In original code, they first have to extract
         //pixels and call XQueryColors.
         //I have only one loop here.
         ULong_t r = 0, g = 0, b = 0;
         for (int y = 0, dotCnt = 0; y < int(source->rows); y++) {
            for (int x = 0; x < int(source->width); x++) {
               if (x + bx < int(pixmap.fWidth) && y + by < int(pixmap.fHeight)) {
                  const unsigned char * const pixels = pixmap.fData + (y + by) * pixmap.fWidth * 4 + (x + bx) * 4;
                  r += UShort_t(pixels[0] / 255. * 0xffff);
                  g += UShort_t(pixels[1] / 255. * 0xffff);
                  b += UShort_t(pixels[2] / 255. * 0xffff);
               }

               if (++dotCnt >= maxDots)
                  break;
            }
         }

         if (dots) {
            r /= dots;
            g /= dots;
            b /= dots;
         }

         if (col[0].fRed == r && col[0].fGreen == g && col[0].fBlue == b) {
            col[0].fPixel = back;
         } else {
            col[0].fPixel = ~back;//???
            col[0].fRed = (UShort_t) r;
            col[0].fGreen = (UShort_t) g;
            col[0].fBlue = (UShort_t) b;
         }
      }

      // if fore or background have changed from previous character
      // recalculate the 3 smoothing colors (interpolation between fore-
      // and background colors)
      if (fore != col[4].fPixel || back != col[0].fPixel) {
         col[4].fPixel = fore;
         TGCocoa::QueryColor(kNone, col[4]);//calculate fRed/fGreen/fBlue triple from fPixel.
         if (back != (ULong_t)-1) {
            col[0].fPixel = back;
            TGCocoa::QueryColor(kNone, col[0]);
         }

         // interpolate between fore and background colors
         for (int x = 3; x > 0; --x) {
            col[x].fRed   = (col[4].fRed   * x + col[0].fRed   * (4 - x)) / 4;
            col[x].fGreen = (col[4].fGreen * x + col[0].fGreen * (4 - x)) / 4;
            col[x].fBlue  = (col[4].fBlue  * x + col[0].fBlue  * (4 - x)) / 4;
            TGCocoa::AllocColor(kNone, col[x]);//Calculate fPixel from fRed/fGreen/fBlue triplet.
         }
      }

      // put smoothed character, character pixmap values are an index
      // into the 5 colors used for aliasing (4 = foreground, 0 = background)
      const unsigned char *s = source->buffer;
      for (int y = 0; y < (int) source->rows; ++y) {
         for (int x = 0; x < (int) source->width; ++x) {
            unsigned char d = *s++ & 0xff;//???
            d = ((d + 10) * 5) / 256;//???
            if (d > 4)
               d = 4;
            if (d && x < (int) source->width) {
               const UChar_t pixel[] = {UChar_t(double(col[d].fRed) / 0xffff * 255),
                                        UChar_t(double(col[d].fGreen) / 0xffff * 255),
                                        UChar_t(double(col[d].fBlue) / 0xffff * 255), 255};
               [pixmap putPixel : pixel X : bx + x Y : by + y];
            }
         }
      }
   } else {
      // no smoothing, just put character using foreground color
      unsigned char rgba[4] = {};
      rgba[3] = 255;
      X11::PixelToRGB(fore, rgba);
      unsigned char d = 0;

      const unsigned char *row = source->buffer;
      for (int y = 0; y < int(source->rows); ++y) {
         int n = 0;
         const unsigned char *s = row;
         for (int x = 0; x < int(source->width); ++x) {
            if (!n)
               d = *s++;

            if (TESTBIT(d,7 - n))
               [pixmap putPixel : rgba X : bx + x Y : by + y];

            if (++n == int(kBitsPerByte))
               n = 0;
         }

         row += source->pitch;
      }
   }
}

//Aux. functions.

//______________________________________________________________________________
void TGQuartz::SetAA()
{
   if (gEnv) {
      const TString value(TString(gEnv->GetValue("Cocoa.EnableAntiAliasing", "auto")).Strip());
      if (value == "auto") {
         [[NSScreen mainScreen] backingScaleFactor] > 1. ? fUseAA = true : fUseAA = false;
      } else if (value == "no")
         fUseAA = false;
      else {
         assert(value == "yes" && "SetAA, value must be 'yes', 'no' or 'auto'");
         fUseAA = true;
      }
      const TString valuefa(TString(gEnv->GetValue("Cocoa.EnableFillAreaAntiAliasing", "auto")).Strip());
      if (valuefa == "auto") {
         [[NSScreen mainScreen] backingScaleFactor] > 1. ? fUseFAAA = true : fUseFAAA = false;
      } else if (valuefa == "no")
         fUseFAAA = false;
      else {
         assert(valuefa == "yes" && "SetAA, value must be 'yes', 'no' or 'auto'");
         fUseFAAA = true;
      }
   }
}

//______________________________________________________________________________
void *TGQuartz::GetSelectedDrawableChecked(const char *calledFrom) const
{
   assert(calledFrom != 0 && "GetSelectedDrawableChecked, calledFrom parameter is null");
   assert(fSelectedDrawable > fPimpl->GetRootWindowID() && "GetSelectedDrawableChecked, bad drawable is selected");

   NSObject<X11Drawable> *drawable = fPimpl->GetDrawable(fSelectedDrawable);
   if (!drawable.fIsPixmap) {
      //TPad/TCanvas ALWAYS draw only into a pixmap.
      if ([drawable isKindOfClass : [QuartzView class]]) {
         QuartzView *view = (QuartzView *)drawable;
         if (!view.fBackBuffer) {
            Error(calledFrom, "Selected window is not double buffered");
            return 0;
         }

         drawable = view.fBackBuffer;
      } else {
         Error(calledFrom, "Selected drawable is neither a pixmap, nor a double buffered window");
         return 0;
      }
   }

   if (!drawable.fContext) {
      Error(calledFrom, "Context is null");
      return 0;
   }

   return drawable;
}
