// @(#)root/graf2d:$Id$
// Author: Timur Pocheptsov, 14/8/2011

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_IOSMarkers
#define ROOT_IOSMarkers

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// IOSMarkers                                                           //
//                                                                      //
// Aux. functions to draw poly-markers.                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>

#include <CoreGraphics/CGContext.h>

#include "Rtypes.h"

#include "TPoint.h"

namespace ROOT {
namespace iOS {
namespace GraphicUtils {

void DrawPolyMarker(CGContextRef ctx, const std::vector<TPoint> &marker, Size_t markerSize, Style_t markerStyle);
void DrawPolyMarker(CGContextRef ctx, unsigned nPoints, const TPoint *marker, Size_t markerSize, Style_t markerStyle);

}
}
}

#endif
