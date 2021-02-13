// @(#)root/gui:$Id$
// Author: Fons Rademakers   02/01/98

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
/**************************************************************************

    This source is based on Xclass95, a Win95-looking GUI toolkit.
    Copyright (C) 1996, 1997 David Barth, Ricky Ralston, Hector Peraza.

    Xclass95 is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

**************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGDimension, TGPosition, TGLongPosition, TGInsets and TGRectangle    //
//                                                                      //
// Several small geometry classes that implement dimensions             //
// (width and height), positions (x and y), insets and rectangles.      //
// They are trivial and their members are public.                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGDimension.h"
#include "TMathBase.h"

void TGRectangle::Merge(const TGRectangle &r)
{
   // Merge parameters
   Int_t max_x = TMath::Max(fX + (Int_t) fW, r.fX + (Int_t) r.fW);
   fX = TMath::Min(fX, r.fX);
   Int_t max_y = TMath::Max(fY + (Int_t) fH, r.fY + (Int_t) r.fH);
   fY = TMath::Min(fY, r.fY);
   fW = max_x - fX;
   fH = max_y - fY;
}
