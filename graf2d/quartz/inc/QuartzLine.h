// @(#)root/graf2d:$Id$
// Author: Olivier Couet, 23/01/2012

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_QuartzLine
#define ROOT_QuartzLine

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// QuartzLine                                                           //
//                                                                      //
// Aux. functions to draw line.                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <Cocoa/Cocoa.h>

#include "Rtypes.h"

#include "TPoint.h"

namespace ROOT {
namespace Quartz {

Bool_t SetLineColor(CGContextRef ctx, Color_t colorIndex);
void SetLineStyle(CGContextRef ctx, Int_t lstyle);
void SetLineWidth(CGContextRef ctx, Int_t width);

void DrawLine(CGContextRef ctx, Int_t x1, Int_t y1, Int_t x2, Int_t y2);
void DrawPolyLine(CGContextRef ctx, Int_t n, TPoint * xy);

}
}

#endif
