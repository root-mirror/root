// @(#)root/graf2d:$Id$
// Author: Olivier Couet, 23/01/2012

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_QuartzFillArea
#define ROOT_QuartzFillArea

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// QuartzFillArea                                                       //
//                                                                      //
// Aux. functions to draw fill area.                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <Cocoa/Cocoa.h>

#include "TAttFill.h"
#include "Rtypes.h"

#include "TPoint.h"

class TColorGradient;

namespace ROOT {
namespace Quartz {

Bool_t SetFillColor(CGContextRef ctx, Color_t colorIndex);
Bool_t SetFillAreaParameters(CGContextRef ctx, unsigned *patternIndex);

void DrawBox(CGContextRef ctx, Int_t x1, Int_t y1, Int_t x2, Int_t y2, bool hollow);
void DrawFillArea(CGContextRef ctx, Int_t n, TPoint *xy, Bool_t drawShadow);

void DrawPolygonWithGradientFill(CGContextRef ctx, const TColorGradient *extendedColor, const CGSize &sizeOfDrawable,
                                 Int_t nPoints, const TPoint *xy, Bool_t drawShadow);

}
}

#endif
