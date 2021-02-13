/* @(#)root/histpainter:$Id$ */

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_Hparam
#define ROOT_Hparam


////////////////////////////////////////////////////////////////////////////////
/*! \struct Hparam_t
\brief Histogram parameters structure.

Structure to store current histogram parameters.
Used internally by `THistpainter` to manage histogram parameters.

*/

#include "RtypesCore.h"

typedef struct Hparam_t {
   Double_t  xbinsize;      ///< bin size in case of equidistant bins
   Double_t  xlowedge;      ///< low edge of axis
   Double_t  xmin;          ///< minimum value along X
   Double_t  xmax;          ///< maximum value along X
   Double_t  ybinsize;      ///< bin size in case of equidistant bins
   Double_t  ylowedge;      ///< low edge of axis
   Double_t  ymin;          ///< minimum value along y
   Double_t  ymax;          ///< maximum value along y
   Double_t  zbinsize;      ///< bin size in case of equidistant bins
   Double_t  zlowedge;      ///< low edge of axis
   Double_t  zmin;          ///< minimum value along Z
   Double_t  zmax;          ///< maximum value along Z
   Double_t  factor;        ///< multiplication factor (normalization)
   Double_t  allchan;       ///< integrated sum of contents
   Double_t  baroffset;     ///< offset of bin for bars or legos [0,1]
   Double_t  barwidth;      ///< width of bin for bars and legos [0,1]
   Int_t     xfirst;        ///< first bin number along X
   Int_t     xlast;         ///< last bin number along X
   Int_t     yfirst;        ///< first bin number along Y
   Int_t     ylast;         ///< last bin number along Y
   Int_t     zfirst;        ///< first bin number along Z
   Int_t     zlast;         ///< last bin number along Z
} Hparam_t;

#endif
