// @(#)root/g3d:$Id$
// Author: Nenad Buncic 17/09/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TBRIK.h"
#include "TNode.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TGeometry.h"

ClassImp(TBRIK);

/** \class TBRIK
\ingroup g3d
A box with faces perpendicular to the axes.

\image html g3d_brik.png
It has 6 parameters:

  - name:       name of the shape
  - title:      shape's title
  - material:  (see TMaterial)
  - dx:         half-length of the box along the x-axis
  - dy:         half-length of the box along the y-axis
  - dz:         half-length of the box along the z-axis
*/

////////////////////////////////////////////////////////////////////////////////
/// BRIK shape default constructor

TBRIK::TBRIK()
{
   fDx = 0.;
   fDy = 0.;
   fDz = 0.;
}

////////////////////////////////////////////////////////////////////////////////
/// BRIK shape normal constructor

TBRIK::TBRIK(const char *name, const char *title, const char *material, Float_t dx, Float_t dy, Float_t dz)
      : TShape (name, title,material)
{
   fDx = dx;
   fDy = dy;
   fDz = dz;
}

////////////////////////////////////////////////////////////////////////////////
/// BRIK shape default destructor

TBRIK::~TBRIK()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from point px,py to a BRIK
///
/// Compute the closest distance of approach from point px,py to each corner
/// point of the BRIK.

Int_t TBRIK::DistancetoPrimitive(Int_t px, Int_t py)
{
   const Int_t numPoints = 8;
   return ShapeDistancetoPrimitive(numPoints,px,py);
}

////////////////////////////////////////////////////////////////////////////////
/// Create BRIK points

void TBRIK::SetPoints(Double_t *points) const
{
   if (points) {
      points[ 0] = -fDx ; points[ 1] = -fDy ; points[ 2] = -fDz;
      points[ 3] = -fDx ; points[ 4] =  fDy ; points[ 5] = -fDz;
      points[ 6] =  fDx ; points[ 7] =  fDy ; points[ 8] = -fDz;
      points[ 9] =  fDx ; points[10] = -fDy ; points[11] = -fDz;
      points[12] = -fDx ; points[13] = -fDy ; points[14] =  fDz;
      points[15] = -fDx ; points[16] =  fDy ; points[17] =  fDz;
      points[18] =  fDx ; points[19] =  fDy ; points[20] =  fDz;
      points[21] =  fDx ; points[22] = -fDy ; points[23] =  fDz;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return total X3D needed by TNode::ls (when called with option "x")

void TBRIK::Sizeof3D() const
{
   gSize3D.numPoints += 8;
   gSize3D.numSegs   += 12;
   gSize3D.numPolys  += 6;
}

////////////////////////////////////////////////////////////////////////////////
/// Get buffer 3D

const TBuffer3D & TBRIK::GetBuffer3D(Int_t reqSections) const
{
   static TBuffer3D buffer(TBuffer3DTypes::kGeneric);

   TShape::FillBuffer3D(buffer, reqSections);

   // No kShapeSpecific or kBoundingBox

   if (reqSections & TBuffer3D::kRawSizes) {
      Int_t nbPnts = 8;
      Int_t nbSegs = 12;
      Int_t nbPols = 6;
      if (buffer.SetRawSizes(nbPnts, nbPnts*3, nbSegs, nbSegs*3, nbPols, nbPols*6)) {
         buffer.SetSectionsValid(TBuffer3D::kRawSizes);
      }
   }
   if ((reqSections & TBuffer3D::kRaw) && buffer.SectionsValid(TBuffer3D::kRawSizes)) {
      // Points
      SetPoints(buffer.fPnts);
      if (!buffer.fLocalFrame) {
         TransformPoints(buffer.fPnts, buffer.NbPnts());
      }

      Int_t c = GetBasicColor();

      // Segments
      buffer.fSegs[ 0] = c   ; buffer.fSegs[ 1] = 0   ; buffer.fSegs[ 2] = 1   ;
      buffer.fSegs[ 3] = c+1 ; buffer.fSegs[ 4] = 1   ; buffer.fSegs[ 5] = 2   ;
      buffer.fSegs[ 6] = c+1 ; buffer.fSegs[ 7] = 2   ; buffer.fSegs[ 8] = 3   ;
      buffer.fSegs[ 9] = c   ; buffer.fSegs[10] = 3   ; buffer.fSegs[11] = 0   ;
      buffer.fSegs[12] = c+2 ; buffer.fSegs[13] = 4   ; buffer.fSegs[14] = 5   ;
      buffer.fSegs[15] = c+2 ; buffer.fSegs[16] = 5   ; buffer.fSegs[17] = 6   ;
      buffer.fSegs[18] = c+3 ; buffer.fSegs[19] = 6   ; buffer.fSegs[20] = 7   ;
      buffer.fSegs[21] = c+3 ; buffer.fSegs[22] = 7   ; buffer.fSegs[23] = 4   ;
      buffer.fSegs[24] = c   ; buffer.fSegs[25] = 0   ; buffer.fSegs[26] = 4   ;
      buffer.fSegs[27] = c+2 ; buffer.fSegs[28] = 1   ; buffer.fSegs[29] = 5   ;
      buffer.fSegs[30] = c+1 ; buffer.fSegs[31] = 2   ; buffer.fSegs[32] = 6   ;
      buffer.fSegs[33] = c+3 ; buffer.fSegs[34] = 3   ; buffer.fSegs[35] = 7   ;

      // Polygons
      buffer.fPols[ 0] = c   ; buffer.fPols[ 1] = 4   ;  buffer.fPols[ 2] = 0  ;
      buffer.fPols[ 3] = 9   ; buffer.fPols[ 4] = 4   ;  buffer.fPols[ 5] = 8  ;
      buffer.fPols[ 6] = c+1 ; buffer.fPols[ 7] = 4   ;  buffer.fPols[ 8] = 1  ;
      buffer.fPols[ 9] = 10  ; buffer.fPols[10] = 5   ;  buffer.fPols[11] = 9  ;
      buffer.fPols[12] = c   ; buffer.fPols[13] = 4   ;  buffer.fPols[14] = 2  ;
      buffer.fPols[15] = 11  ; buffer.fPols[16] = 6   ;  buffer.fPols[17] = 10 ;
      buffer.fPols[18] = c+1 ; buffer.fPols[19] = 4   ;  buffer.fPols[20] = 3  ;
      buffer.fPols[21] = 8   ; buffer.fPols[22] = 7   ;  buffer.fPols[23] = 11 ;
      buffer.fPols[24] = c+2 ; buffer.fPols[25] = 4   ;  buffer.fPols[26] = 0  ;
      buffer.fPols[27] = 3   ; buffer.fPols[28] = 2   ;  buffer.fPols[29] = 1  ;
      buffer.fPols[30] = c+3 ; buffer.fPols[31] = 4   ;  buffer.fPols[32] = 4  ;
      buffer.fPols[33] = 5   ; buffer.fPols[34] = 6   ;  buffer.fPols[35] = 7  ;

      buffer.SetSectionsValid(TBuffer3D::kRaw);
   }
   return buffer;
}
