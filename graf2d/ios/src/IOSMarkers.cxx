// @(#)root/graf2d:$Id$
// Author: Timur Pocheptsov, 14/8/2011

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TAttMarker.h"

#include "IOSGraphicUtils.h"
#include "IOSMarkers.h"

namespace ROOT {
namespace iOS {
namespace GraphicUtils {

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerDot(CGContextRef ctx, unsigned n, const TPoint *xy)
{
   for (unsigned i = 0; i < n; ++i)
      CGContextFillRect(ctx, CGRectMake(xy[i].fX, xy[i].fY, 1.f, 1.f));
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerPlus(CGContextRef ctx, unsigned n, const TPoint *xy, Size_t markerSize)
{
   const Double_t im = 4 * markerSize + 0.5;

   for (UInt_t i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, -im + x, y);
      CGContextAddLineToPoint(ctx, im + x, y);
      CGContextStrokePath(ctx);

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x, -im + y);
      CGContextAddLineToPoint(ctx, x, im + y);
      CGContextStrokePath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerStar(CGContextRef ctx, unsigned n, const TPoint *xy, Size_t markerSize)
{
   Double_t im = 4 * markerSize + 0.5;

   TPoint star[8];
   star[0].fX = -im;  star[0].fY = 0;
   star[1].fX =  im;  star[1].fY = 0;
   star[2].fX = 0  ;  star[2].fY = -im;
   star[3].fX = 0  ;  star[3].fY = im;

   im = 0.707 * im + 0.5;
   star[4].fX = -im;  star[4].fY = -im;
   star[5].fX =  im;  star[5].fY = im;
   star[6].fX = -im;  star[6].fY = im;
   star[7].fX =  im;  star[7].fY = -im;

   for (UInt_t i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, star[0].fX + x, star[0].fY + y);
      CGContextAddLineToPoint(ctx, star[1].fX + x, star[1].fY + y);
      CGContextStrokePath(ctx);

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, star[2].fX + x, star[2].fY + y);
      CGContextAddLineToPoint(ctx, star[3].fX + x, star[3].fY + y);
      CGContextStrokePath(ctx);

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, star[4].fX + x, star[4].fY + y);
      CGContextAddLineToPoint(ctx, star[5].fX + x, star[5].fY + y);
      CGContextStrokePath(ctx);

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, star[6].fX + x, star[6].fY + y);
      CGContextAddLineToPoint(ctx, star[7].fX + x, star[7].fY + y);
      CGContextStrokePath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerOpenCircle(CGContextRef ctx, unsigned n, const TPoint *xy, Size_t markerSize)
{
   Double_t r = 4 * markerSize + 0.5;
   if (r > 100.)
      r = 100.;//as in TGX11.

   const Double_t d = 2 * r;

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      const CGRect rect = CGRectMake(x - r, y - r, d, d);
      CGContextStrokeEllipseInRect(ctx, rect);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerX(CGContextRef ctx, unsigned n, const TPoint *xy, Size_t markerSize)
{
   const Double_t im = 0.707 * (4 * markerSize + 0.5) + 0.5;
   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, -im + x, -im + y);
      CGContextAddLineToPoint(ctx, im + x, im + y);
      CGContextStrokePath(ctx);

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, -im + x, im + y);
      CGContextAddLineToPoint(ctx, im + x, -im + y);
      CGContextStrokePath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerFullDotSmall(CGContextRef ctx, unsigned n, const TPoint *xy)
{
   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, -1. + x, y);
      CGContextAddLineToPoint(ctx, x + 1., y);
      CGContextStrokePath(ctx);

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x, -1. + y);
      CGContextAddLineToPoint(ctx, x, 1. + y);
      CGContextStrokePath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerFullDotMedium(CGContextRef ctx, unsigned n, const TPoint *xy)
{
   for (unsigned i = 0; i < n; ++i)
      CGContextFillRect(ctx, CGRectMake(xy[i].fX - 1, xy[i].fY - 1, 3.f, 3.f));
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerFullDotLarge(CGContextRef ctx, unsigned n, const TPoint *xy, Size_t markerSize)
{
   Double_t radius = 4 * markerSize + 0.5;
   if (radius > 100.)
      radius = 100;//as in TGX11.

   const Double_t d = 2 * radius;

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      const CGRect rect = CGRectMake(x - radius, y - radius, d, d);
      CGContextFillEllipseInRect(ctx, rect);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerFullSquare(CGContextRef ctx, unsigned n, const TPoint *xy, Size_t markerSize)
{
   const Double_t im = 4 * markerSize + 0.5;
   for (unsigned i = 0; i < n; ++i) {
      const CGRect rect = CGRectMake(xy[i].fX - im, xy[i].fY - im, im * 2, im * 2);
      CGContextFillRect(ctx, rect);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerOpenSquare(CGContextRef ctx, unsigned n, const TPoint *xy, Size_t markerSize)
{
   const Double_t im = 4 * markerSize + 0.5;
   for (unsigned i = 0; i < n; ++i) {
      const CGRect rect = CGRectMake(xy[i].fX - im, xy[i].fY - im, im * 2, im * 2);
      CGContextStrokeRect(ctx, rect);
   }
}


////////////////////////////////////////////////////////////////////////////////

void DrawMarkerFullTriangleUp(CGContextRef ctx, unsigned n, const TPoint *xy, Size_t markerSize)
{
   const Double_t im = 4 * markerSize + 0.5;
   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;
      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x - im, y - im);
      CGContextAddLineToPoint(ctx, x + im, y - im);
      CGContextAddLineToPoint(ctx, x, im + y);
      CGContextFillPath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerOpenTriangleUp(CGContextRef ctx, unsigned n, const TPoint *xy, Size_t markerSize)
{
   const Double_t im = 4 * markerSize + 0.5;
   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;
      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x - im, y - im);
      CGContextAddLineToPoint(ctx, x + im, y - im);
      CGContextAddLineToPoint(ctx, x, im + y);
      CGContextAddLineToPoint(ctx, x - im, y - im);
      CGContextStrokePath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerOpenTriangleDown(CGContextRef ctx, unsigned n, const TPoint *xy, Size_t markerSize)
{
   const Int_t im = Int_t(4 * markerSize + 0.5);
   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x - im, y + im);
      CGContextAddLineToPoint(ctx, x, y - im);
      CGContextAddLineToPoint(ctx, im + x, y + im);
      CGContextAddLineToPoint(ctx, x - im, y + im);
      CGContextStrokePath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerFullTriangleDown(CGContextRef ctx, unsigned n, const TPoint *xy, Size_t markerSize)
{
   const Int_t im = Int_t(4 * markerSize + 0.5);
   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x - im, y + im);
      CGContextAddLineToPoint(ctx, x, y - im);
      CGContextAddLineToPoint(ctx, im + x, y + im);
      CGContextFillPath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerFullDiamond(CGContextRef ctx, unsigned n, const TPoint *xy, Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);
   const Int_t imx = Int_t(2.66 * markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x - imx,  y);
      CGContextAddLineToPoint(ctx, x, y - im);
      CGContextAddLineToPoint(ctx, x + imx, y);
      CGContextAddLineToPoint(ctx, x, y + im);
      CGContextDrawPath(ctx, kCGPathFillStroke);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerOpenDiamond(CGContextRef ctx, unsigned n, const TPoint *xy, Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);
   const Int_t imx = Int_t(2.66 * markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x - imx,  y);
      CGContextAddLineToPoint(ctx, x, y - im);
      CGContextAddLineToPoint(ctx, x + imx, y);
      CGContextAddLineToPoint(ctx, x, y + im);
      CGContextAddLineToPoint(ctx, x - imx,  y);
      CGContextStrokePath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerFullCross(CGContextRef ctx, unsigned n, const TPoint *xy, Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);
   const Int_t imx = Int_t(1.33 * markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x - im, y - imx);
      CGContextAddLineToPoint(ctx, x - imx, y - imx);
      CGContextAddLineToPoint(ctx, x - imx, y - im);
      CGContextAddLineToPoint(ctx, x + imx, y - im);
      CGContextAddLineToPoint(ctx, x + imx, y - imx);
      CGContextAddLineToPoint(ctx, x + im, y - imx);
      CGContextAddLineToPoint(ctx, x + im, y + imx);
      CGContextAddLineToPoint(ctx, x + imx, y + imx);
      CGContextAddLineToPoint(ctx, x + imx, y + im);
      CGContextAddLineToPoint(ctx, x - imx, y + im);
      CGContextAddLineToPoint(ctx, x - imx, y + imx);
      CGContextAddLineToPoint(ctx, x - im, y + imx);
      CGContextAddLineToPoint(ctx, x - im, y - imx);
      CGContextFillPath(ctx);
   }
}


////////////////////////////////////////////////////////////////////////////////

void DrawMarkerOpenCross(CGContextRef ctx, unsigned n, const TPoint *xy, Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);
   const Int_t imx = Int_t(1.33 * markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x - im, y - imx);
      CGContextAddLineToPoint(ctx, x - imx, y - imx);
      CGContextAddLineToPoint(ctx, x - imx, y - im);
      CGContextAddLineToPoint(ctx, x + imx, y - im);
      CGContextAddLineToPoint(ctx, x + imx, y - imx);
      CGContextAddLineToPoint(ctx, x + im, y - imx);
      CGContextAddLineToPoint(ctx, x + im, y + imx);
      CGContextAddLineToPoint(ctx, x + imx, y + imx);
      CGContextAddLineToPoint(ctx, x + imx, y + im);
      CGContextAddLineToPoint(ctx, x - imx, y + im);
      CGContextAddLineToPoint(ctx, x - imx, y + imx);
      CGContextAddLineToPoint(ctx, x - im, y + imx);
      CGContextAddLineToPoint(ctx, x - im, y - imx);
      CGContextStrokePath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// HIGZ full star pentagone

void DrawMarkerFullStar(CGContextRef ctx, unsigned n, const TPoint *xy, Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);
   const Int_t im1 = Int_t(0.66 * markerSize + 0.5);
   const Int_t im2 = Int_t(2.00 * markerSize + 0.5);
   const Int_t im3 = Int_t(2.66 * markerSize + 0.5);
   const Int_t im4 = Int_t(1.33 * markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x - im, y - im4);
      CGContextAddLineToPoint(ctx, x - im2, y + im1);
      CGContextAddLineToPoint(ctx, x - im4, y - im4);
      CGContextFillPath(ctx);

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x - im2, y + im1);//1
      CGContextAddLineToPoint(ctx, x - im3, y + im);//2
      CGContextAddLineToPoint(ctx, x, y + im2);//3
      CGContextFillPath(ctx);

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x, y + im2);//3
      CGContextAddLineToPoint(ctx, x + im3, y + im);//4
      CGContextAddLineToPoint(ctx, x + im2, y + im1);//5
      CGContextFillPath(ctx);

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x + im2, y + im1);//5
      CGContextAddLineToPoint(ctx, x + im, y - im4);//6
      CGContextAddLineToPoint(ctx,x + im4, y - im4);//7
      CGContextFillPath(ctx);

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x + im4, y - im4);//7
      CGContextAddLineToPoint(ctx, x, y - im);//8
      CGContextAddLineToPoint(ctx, x - im4, y - im4);//9
      CGContextFillPath(ctx);

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x - im4, y - im4);//9
      CGContextAddLineToPoint(ctx, x - im2, y + im1);//1
      CGContextAddLineToPoint(ctx, x, y + im2);//3
      CGContextFillPath(ctx);

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x - im4, y - im4);//9
      CGContextAddLineToPoint(ctx, x, y + im2);//3
      CGContextAddLineToPoint(ctx, x + im2, y + im1);//5
      CGContextFillPath(ctx);

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x - im4, y - im4);//9
      CGContextAddLineToPoint(ctx, x + im2, y + im1);//5
      CGContextAddLineToPoint(ctx, x + im4, y - im4);//7
      CGContextFillPath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerOpenStar(CGContextRef ctx, unsigned n, const TPoint *xy, Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);
   const Int_t im1 = Int_t(0.66 * markerSize + 0.5);
   const Int_t im2 = Int_t(2.00 * markerSize + 0.5);
   const Int_t im3 = Int_t(2.66 * markerSize + 0.5);
   const Int_t im4 = Int_t(1.33 * markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x - im, y - im4);
      CGContextAddLineToPoint(ctx, x - im2, y + im1);
      CGContextAddLineToPoint(ctx, x - im3, y + im);
      CGContextAddLineToPoint(ctx, x, y + im2);
      CGContextAddLineToPoint(ctx, x + im3, y + im);
      CGContextAddLineToPoint(ctx, x + im2, y + im1);
      CGContextAddLineToPoint(ctx, x + im, y - im4);
      CGContextAddLineToPoint(ctx, x + im4, y - im4);
      CGContextAddLineToPoint(ctx, x, y - im);
      CGContextAddLineToPoint(ctx, x - im4, y - im4);
      CGContextAddLineToPoint(ctx, x - im, y - im4);
      CGContextStrokePath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerOpenSquareDiagonal(CGContextRef ctx, unsigned n, const TPoint *xy,
                          Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x - im, y - im);
      CGContextAddLineToPoint(ctx, x + im, y - im);
      CGContextAddLineToPoint(ctx, x + im, y + im);
      CGContextAddLineToPoint(ctx, x - im, y + im);
      CGContextAddLineToPoint(ctx, x - im, y - im);
      CGContextAddLineToPoint(ctx, x + im, y + im);
      CGContextStrokePath(ctx);

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x - im, y + im);
      CGContextAddLineToPoint(ctx, x + im, y - im);
      CGContextStrokePath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerOpenDiamondCross(CGContextRef ctx, unsigned n, const TPoint *xy,
                          Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x - im, y    );
      CGContextAddLineToPoint(ctx, x     , y - im);
      CGContextAddLineToPoint(ctx, x + im, y     );
      CGContextAddLineToPoint(ctx, x     , y + im);
      CGContextAddLineToPoint(ctx, x - im, y     );
      CGContextAddLineToPoint(ctx, x + im, y     );
      CGContextStrokePath(ctx);

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x     , y + im);
      CGContextAddLineToPoint(ctx, x     , y - im);
      CGContextStrokePath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerOpenThreeTriangles(CGContextRef ctx, unsigned n, const TPoint *xy,
                          Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);
   const Int_t im2 = Int_t(2.00 * markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x     , y     );
      CGContextAddLineToPoint(ctx, x -im2, y + im);
      CGContextAddLineToPoint(ctx, x - im, y     );
      CGContextAddLineToPoint(ctx, x     , y     );
      CGContextAddLineToPoint(ctx, x -im2, y - im);
      CGContextAddLineToPoint(ctx, x +im2, y - im);
      CGContextAddLineToPoint(ctx, x     , y     );
      CGContextAddLineToPoint(ctx, x + im, y     );
      CGContextAddLineToPoint(ctx, x +im2, y + im);
      CGContextAddLineToPoint(ctx, x     , y     );
      CGContextStrokePath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerOctagonCross(CGContextRef ctx, unsigned n, const TPoint *xy,
                          Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);
   const Int_t im2 = Int_t(2.00 * markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x - im, y     );
      CGContextAddLineToPoint(ctx, x - im, y -im2);
      CGContextAddLineToPoint(ctx, x -im2, y - im);
      CGContextAddLineToPoint(ctx, x +im2, y - im);
      CGContextAddLineToPoint(ctx, x + im, y -im2);
      CGContextAddLineToPoint(ctx, x + im, y +im2);
      CGContextAddLineToPoint(ctx, x +im2, y + im);
      CGContextAddLineToPoint(ctx, x -im2, y + im);
      CGContextAddLineToPoint(ctx, x - im, y +im2);
      CGContextAddLineToPoint(ctx, x - im, y     );
      CGContextAddLineToPoint(ctx, x + im, y     );
      CGContextStrokePath(ctx);

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x     , y - im);
      CGContextAddLineToPoint(ctx, x    , y + im);
      CGContextStrokePath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerFilledThreeTriangles(CGContextRef ctx, unsigned n, const TPoint *xy,
                          Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);
   const Int_t im2 = Int_t(2.00 * markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x     , y     );
      CGContextAddLineToPoint(ctx, x -im2, y + im);
      CGContextAddLineToPoint(ctx, x - im, y     );
      CGContextAddLineToPoint(ctx, x     , y     );
      CGContextAddLineToPoint(ctx, x -im2, y - im);
      CGContextAddLineToPoint(ctx, x +im2, y - im);
      CGContextAddLineToPoint(ctx, x     , y     );
      CGContextAddLineToPoint(ctx, x + im, y     );
      CGContextAddLineToPoint(ctx, x +im2, y + im);
      CGContextAddLineToPoint(ctx, x     , y     );
      CGContextFillPath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerOpenFourTrianglesX(CGContextRef ctx, unsigned n, const TPoint *xy,
                          Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);
   const Int_t im2 = Int_t(2.00 * markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x     , y     );
      CGContextAddLineToPoint(ctx, x +im2, y + im);
      CGContextAddLineToPoint(ctx, x + im, y +im2);
      CGContextAddLineToPoint(ctx, x     , y     );
      CGContextAddLineToPoint(ctx, x + im, y -im2);
      CGContextAddLineToPoint(ctx, x +im2, y - im);
      CGContextAddLineToPoint(ctx, x     , y     );
      CGContextAddLineToPoint(ctx, x -im2, y - im);
      CGContextAddLineToPoint(ctx, x - im, y -im2);
      CGContextAddLineToPoint(ctx, x     , y     );
      CGContextAddLineToPoint(ctx, x - im, y +im2);
      CGContextAddLineToPoint(ctx, x -im2, y + im);
      CGContextAddLineToPoint(ctx, x     , y     );
      CGContextStrokePath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerFilledFourTrianglesX(CGContextRef ctx, unsigned n, const TPoint *xy,
                          Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);
   const Int_t im2 = Int_t(2.00 * markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x     , y     );
      CGContextAddLineToPoint(ctx, x +im2, y + im);
      CGContextAddLineToPoint(ctx, x + im, y +im2);
      CGContextAddLineToPoint(ctx, x     , y     );
      CGContextAddLineToPoint(ctx, x + im, y -im2);
      CGContextAddLineToPoint(ctx, x +im2, y - im);
      CGContextAddLineToPoint(ctx, x     , y     );
      CGContextAddLineToPoint(ctx, x -im2, y - im);
      CGContextAddLineToPoint(ctx, x - im, y -im2);
      CGContextAddLineToPoint(ctx, x     , y     );
      CGContextAddLineToPoint(ctx, x - im, y +im2);
      CGContextAddLineToPoint(ctx, x -im2, y + im);
      CGContextAddLineToPoint(ctx, x     , y     );
      CGContextFillPath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerOpenDoubleDiamond(CGContextRef ctx, unsigned n, const TPoint *xy,
                          Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);
   const Int_t im4 = Int_t(markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x     , y + im);
      CGContextAddLineToPoint(ctx, x -im4, y +im4);
      CGContextAddLineToPoint(ctx, x - im, y     );
      CGContextAddLineToPoint(ctx, x -im4, y -im4);
      CGContextAddLineToPoint(ctx, x     , y - im);
      CGContextAddLineToPoint(ctx, x +im4, y -im4);
      CGContextAddLineToPoint(ctx, x + im, y     );
      CGContextAddLineToPoint(ctx, x +im4, y +im4);
      CGContextAddLineToPoint(ctx, x     , y + im);
      CGContextStrokePath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerFilledDoubleDiamond(CGContextRef ctx, unsigned n, const TPoint *xy,
                          Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);
   const Int_t im4 = Int_t(markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x     , y + im);
      CGContextAddLineToPoint(ctx, x -im4, y +im4);
      CGContextAddLineToPoint(ctx, x - im, y     );
      CGContextAddLineToPoint(ctx, x -im4, y -im4);
      CGContextAddLineToPoint(ctx, x     , y - im);
      CGContextAddLineToPoint(ctx, x +im4, y -im4);
      CGContextAddLineToPoint(ctx, x + im, y     );
      CGContextAddLineToPoint(ctx, x +im4, y +im4);
      CGContextAddLineToPoint(ctx, x     , y + im);
      CGContextFillPath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerOpenFourTrianglesPlus(CGContextRef ctx, unsigned n, const TPoint *xy,
                          Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);
   const Int_t im2 = Int_t(2.00 * markerSize + 0.5);
   const Int_t im4 = Int_t(1.33 * markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x     , y     );
      CGContextAddLineToPoint(ctx, x +im2, y + im);
      CGContextAddLineToPoint(ctx, x -im2, y + im);
      CGContextAddLineToPoint(ctx, x +im2, y - im);
      CGContextAddLineToPoint(ctx, x -im2, y - im);
      CGContextAddLineToPoint(ctx, x     , y     );
      CGContextAddLineToPoint(ctx, x + im, y +im2);
      CGContextAddLineToPoint(ctx, x + im, y -im2);
      CGContextAddLineToPoint(ctx, x - im, y +im2);
      CGContextAddLineToPoint(ctx, x - im, y -im2);
      CGContextAddLineToPoint(ctx, x     , y     );
      CGContextStrokePath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerFilledFourTrianglesPlus(CGContextRef ctx, unsigned n, const TPoint *xy,
                          Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);
   const Int_t im2 = Int_t(2.00 * markerSize + 0.5);
   const Int_t im4 = Int_t(0.2 * markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x +im4, y +im4);
      CGContextAddLineToPoint(ctx, x +im2, y + im);
      CGContextAddLineToPoint(ctx, x -im2, y + im);
      CGContextAddLineToPoint(ctx, x -im4, y +im4);
      CGContextAddLineToPoint(ctx, x - im, y +im2);
      CGContextAddLineToPoint(ctx, x - im, y -im2);
      CGContextAddLineToPoint(ctx, x -im4, y -im4);
      CGContextAddLineToPoint(ctx, x -im2, y - im);
      CGContextAddLineToPoint(ctx, x +im2, y - im);
      CGContextAddLineToPoint(ctx, x +im4, y -im4);
      CGContextAddLineToPoint(ctx, x + im, y -im2);
      CGContextAddLineToPoint(ctx, x + im, y +im2);
      CGContextAddLineToPoint(ctx, x +im4, y +im4);
      CGContextFillPath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerOpenCrossX(CGContextRef ctx, unsigned n, const TPoint *xy,
                          Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);
   const Int_t im2 = Int_t(2.00 * markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x     , y + im2);
      CGContextAddLineToPoint(ctx, x -im2, y + im);
      CGContextAddLineToPoint(ctx, x - im, y +im2);
      CGContextAddLineToPoint(ctx, x -im2, y     );
      CGContextAddLineToPoint(ctx, x - im, y -im2);
      CGContextAddLineToPoint(ctx, x -im2, y - im);
      CGContextAddLineToPoint(ctx, x     , y -im2);
      CGContextAddLineToPoint(ctx, x +im2, y - im);
      CGContextAddLineToPoint(ctx, x + im, y -im2);
      CGContextAddLineToPoint(ctx, x +im2, y     );
      CGContextAddLineToPoint(ctx, x + im, y +im2);
      CGContextAddLineToPoint(ctx, x +im2, y + im);
      CGContextAddLineToPoint(ctx, x     , y +im2);
      CGContextStrokePath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerFilledCrossX(CGContextRef ctx, unsigned n, const TPoint *xy,
                          Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);
   const Int_t im2 = Int_t(2.0 * markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x -im2, y - im2*1.005);
      CGContextAddLineToPoint(ctx, x -im2, y - im);
      CGContextAddLineToPoint(ctx, x +im2, y - im);
      CGContextAddLineToPoint(ctx, x +im2, y -im2);
      CGContextAddLineToPoint(ctx, x + im, y -im2);
      CGContextAddLineToPoint(ctx, x + im, y +im2);
      CGContextAddLineToPoint(ctx, x +im2, y +im2);
      CGContextAddLineToPoint(ctx, x +im2, y + im);
      CGContextAddLineToPoint(ctx, x -im2, y + im);
      CGContextAddLineToPoint(ctx, x -im2, y +im2);
      CGContextAddLineToPoint(ctx, x - im, y +im2);
      CGContextAddLineToPoint(ctx, x - im, y -im2);
      CGContextAddLineToPoint(ctx, x -im2, y -im2*0.995);
      CGContextAddLineToPoint(ctx, x -im2, y +im2);
      CGContextAddLineToPoint(ctx, x +im2, y +im2);
      CGContextAddLineToPoint(ctx, x +im2, y -im2);
      CGContextAddLineToPoint(ctx, x -im2, y -im2*1.005);
      CGContextFillPath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerFourSuaresX(CGContextRef ctx, unsigned n, const TPoint *xy,
                          Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);
   const Int_t im2 = Int_t(2.00 * markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x     , y + im2*1.01);
      CGContextAddLineToPoint(ctx, x -im2, y + im);
      CGContextAddLineToPoint(ctx, x - im, y +im2);
      CGContextAddLineToPoint(ctx, x -im2, y     );
      CGContextAddLineToPoint(ctx, x - im, y -im2);
      CGContextAddLineToPoint(ctx, x -im2, y - im);
      CGContextAddLineToPoint(ctx, x     , y -im2);
      CGContextAddLineToPoint(ctx, x +im2, y - im);
      CGContextAddLineToPoint(ctx, x + im, y -im2);
      CGContextAddLineToPoint(ctx, x +im2, y     );
      CGContextAddLineToPoint(ctx, x + im, y +im2);
      CGContextAddLineToPoint(ctx, x +im2, y + im);
      CGContextAddLineToPoint(ctx, x     , y +im2*0.99);
      CGContextAddLineToPoint(ctx, x +im2*0.99, y     );
      CGContextAddLineToPoint(ctx, x     , y -im2*0.99);
      CGContextAddLineToPoint(ctx, x -im2*0.99, y     );
      CGContextAddLineToPoint(ctx, x     , y +im2*0.99);
      CGContextFillPath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawMarkerFourSuaresPlus(CGContextRef ctx, unsigned n, const TPoint *xy,
                          Size_t markerSize)
{
   const Int_t im  = Int_t(4 * markerSize + 0.5);
   const Int_t im2 = Int_t(1.33 * markerSize + 0.5);

   for (unsigned i = 0; i < n; ++i) {
      const Double_t x = xy[i].fX;
      const Double_t y = xy[i].fY;

      CGContextBeginPath(ctx);
      CGContextMoveToPoint(ctx, x -im2, y - im2*1.01);
      CGContextAddLineToPoint(ctx, x -im2, y + im);
      CGContextAddLineToPoint(ctx, x - im, y +im2);
      CGContextAddLineToPoint(ctx, x -im2, y     );
      CGContextAddLineToPoint(ctx, x - im, y -im2);
      CGContextAddLineToPoint(ctx, x -im2, y - im);
      CGContextAddLineToPoint(ctx, x     , y -im2);
      CGContextAddLineToPoint(ctx, x +im2, y - im);
      CGContextAddLineToPoint(ctx, x + im, y -im2);
      CGContextAddLineToPoint(ctx, x +im2, y     );
      CGContextAddLineToPoint(ctx, x + im, y +im2);
      CGContextAddLineToPoint(ctx, x +im2, y + im);
      CGContextAddLineToPoint(ctx, x     , y +im2*0.99);
      CGContextAddLineToPoint(ctx, x +im2*0.99, y     );
      CGContextAddLineToPoint(ctx, x     , y -im2*0.99);
      CGContextAddLineToPoint(ctx, x -im2*0.99, y     );
      CGContextAddLineToPoint(ctx, x     , y +im2*0.99);
      CGContextFillPath(ctx);
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawPolyMarker(CGContextRef ctx, unsigned nPoints, const TPoint *xy, Size_t markerSize, Style_t markerStyle)
{
   switch (markerStyle) {
   case kDot:
      DrawMarkerDot(ctx, nPoints, xy);
      break;
   case kPlus:
      DrawMarkerPlus(ctx, nPoints, xy, markerSize);
      break;
   case kStar:
      DrawMarkerStar(ctx, nPoints, xy, markerSize);
      break;
   case kCircle:
   case kOpenCircle:
      DrawMarkerOpenCircle(ctx, nPoints, xy, markerSize);
      break;
   case kMultiply:
      DrawMarkerX(ctx, nPoints, xy, markerSize);
      break;
   case kFullDotSmall:
      DrawMarkerFullDotSmall(ctx, nPoints, xy);
      break;
   case kFullDotMedium:
      DrawMarkerFullDotMedium(ctx, nPoints, xy);
      break;
   case kFullDotLarge:
   case kFullCircle:
      DrawMarkerFullDotLarge(ctx, nPoints, xy, markerSize);
      break;
   case kFullSquare:
      DrawMarkerFullSquare(ctx, nPoints, xy, markerSize);
      break;
   case kFullTriangleUp:
      DrawMarkerFullTriangleUp(ctx, nPoints, xy, markerSize);
      break;
   case kFullTriangleDown:
      DrawMarkerFullTriangleDown(ctx, nPoints, xy, markerSize);
      break;
   case kOpenSquare:
      DrawMarkerOpenSquare(ctx, nPoints, xy, markerSize);
      break;
   case kOpenTriangleUp:
      DrawMarkerOpenTriangleUp(ctx, nPoints, xy, markerSize);
      break;
   case kOpenTriangleDown:
      DrawMarkerOpenTriangleDown(ctx, nPoints, xy, markerSize);
      break;
   case kOpenDiamond:
      DrawMarkerOpenDiamond(ctx, nPoints, xy, markerSize);
      break;
   case kFullDiamond:
      DrawMarkerFullDiamond(ctx, nPoints, xy, markerSize);
      break;
   case kOpenCross:
      DrawMarkerOpenCross(ctx, nPoints, xy, markerSize);
      break;
   case kFullCross:
      DrawMarkerFullCross(ctx, nPoints, xy, markerSize);
      break;
   case kFullStar:
      DrawMarkerFullStar(ctx, nPoints, xy, markerSize);
      break;
   case kOpenStar:
      DrawMarkerOpenStar(ctx, nPoints, xy, markerSize);
      break;
   case kOpenDiamondCross:     // YS
      DrawMarkerOpenDiamondCross(ctx, nPoints, xy, markerSize);
      break;
   case kOpenSquareDiagonal:   // YS
      DrawMarkerOpenSquareDiagonal(ctx, nPoints, xy, markerSize);
      break;
   case kOpenThreeTriangles:   // YS
      DrawMarkerOpenThreeTriangles(ctx, nPoints, xy, markerSize);
      break;
   case kOctagonCross:   // YS
      DrawMarkerOctagonCross(ctx, nPoints, xy, markerSize);
      break;
   case kFullThreeTriangles:   // YS
      DrawMarkerFilledThreeTriangles(ctx, nPoints, xy, markerSize);
      break;
   case kOpenFourTrianglesX:   // YS
      DrawMarkerOpenFourTrianglesX(ctx, nPoints, xy, markerSize);
      break;
   case kFullFourTrianglesX:   // YS
      DrawMarkerFilledFourTrianglesX(ctx, nPoints, xy, markerSize);
      break;
   case kOpenDoubleDiamond:   // YS
      DrawMarkerOpenDoubleDiamond(ctx, nPoints, xy, markerSize);
      break;
   case kFullDoubleDiamond:   // YS
      DrawMarkerFilledDoubleDiamond(ctx, nPoints, xy, markerSize);
      break;
   case kOpenFourTrianglesPlus:   // YS
      DrawMarkerOpenFourTrianglesPlus(ctx, nPoints, xy, markerSize);
      break;
   case kFullFourTrianglesPlus:   // YS
      DrawMarkerFilledFourTrianglesPlus(ctx, nPoints, xy, markerSize);
      break;
   case kOpenCrossX:   // YS
      DrawMarkerOpenCrossX(ctx, nPoints, xy, markerSize);
      break;
   case kFullCrossX:   // YS
      DrawMarkerFilledCrossX(ctx, nPoints, xy, markerSize);
      break;
   case kFourSquaresX:   // YS
      DrawMarkerFourSquaresX(ctx, nPoints, xy, markerSize);
      break;
   case kFourSquaresPlus:   // YS
      DrawMarkerFourSquaresPlus(ctx, nPoints, xy, markerSize);
      break;
   }
}

////////////////////////////////////////////////////////////////////////////////

void DrawPolyMarker(CGContextRef ctx, const std::vector<TPoint> &xy, Size_t markerSize, Style_t markerStyle)
{
   DrawPolyMarker(ctx, xy.size(), &xy[0], markerSize, markerStyle);
}

}//namespace GraphicUtils
}//namespace iOS
}//namespace ROOT
