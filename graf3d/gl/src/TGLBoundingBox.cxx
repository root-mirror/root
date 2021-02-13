// @(#)root/gl:$Id$
// Author:  Richard Maunder  25/05/2005

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <iostream>
#include "TGLBoundingBox.h"
#include "TGLIncludes.h"
#include "TMathBase.h"

using namespace std;

/** \class TGLBoundingBox
\ingroup opengl
Concrete class describing an orientated (free) or axis aligned box
of 8 vertices. Supports methods for setting aligned or orientated
boxes, find volume, axes, extents, centers, face planes etc.
Also tests for overlap testing of planes and other bounding boxes,
with fast sphere approximation.
*/

ClassImp(TGLBoundingBox);

////////////////////////////////////////////////////////////////////////////////
/// Construct an empty bounding box

TGLBoundingBox::TGLBoundingBox()
{
   SetEmpty();
}

////////////////////////////////////////////////////////////////////////////////
/// Construct a bounding box from provided 8 vertices

TGLBoundingBox::TGLBoundingBox(const TGLVertex3 vertex[8])
{
   Set(vertex);
}

////////////////////////////////////////////////////////////////////////////////
/// Construct a bounding box from provided 8 vertices

TGLBoundingBox::TGLBoundingBox(const Double_t vertex[8][3])
{
   Set(vertex);
}

////////////////////////////////////////////////////////////////////////////////
/// Construct an global axis ALIGNED bounding box from provided low/high vertex pair

TGLBoundingBox::TGLBoundingBox(const TGLVertex3 & lowVertex, const TGLVertex3 & highVertex)
{
   SetAligned(lowVertex, highVertex);
}

////////////////////////////////////////////////////////////////////////////////
/// Construct a bounding box as copy of existing one

TGLBoundingBox::TGLBoundingBox(const TGLBoundingBox & other)
{
   Set(other);
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy bounding box

TGLBoundingBox::~TGLBoundingBox()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Update the internally cached volume and axes vectors - these are retained
/// for efficiency - many more reads than modifications

void TGLBoundingBox::UpdateCache()
{
   //    y
   //    |
   //    |
   //    |________x
   //   /  3-------2
   //  /  /|      /|
   // z  7-------6 |
   //    | 0-----|-1
   //    |/      |/
   //    4-------5
   //

   // Do axes first so Extents() is correct
   fAxes[0].Set(fVertex[1] - fVertex[0]);
   fAxes[1].Set(fVertex[3] - fVertex[0]);
   fAxes[2].Set(fVertex[4] - fVertex[0]);

   // Sometimes have zero volume BB due to single zero magnitude
   // axis record and try to fix below
   Bool_t fixZeroMagAxis = kFALSE;
   Int_t zeroMagAxisInd = -1;
   for (UInt_t i = 0; i<3; i++) {
      fAxesNorm[i] = fAxes[i];
      Double_t mag = fAxesNorm[i].Mag();
      if (mag > 0.0) {
         fAxesNorm[i] /= mag;
      } else {
         if (!fixZeroMagAxis && zeroMagAxisInd == -1) {
            zeroMagAxisInd = i;
            fixZeroMagAxis = kTRUE;
         } else if (fixZeroMagAxis) {
            fixZeroMagAxis = kFALSE;
         }
      }
   }

   // Try to cope with a zero volume bounding box where one
   // axis is zero by using cross product of other two
   if (fixZeroMagAxis) {
      fAxesNorm[zeroMagAxisInd] = Cross(fAxesNorm[(zeroMagAxisInd+1)%3],
                                        fAxesNorm[(zeroMagAxisInd+2)%3]);
   }

   TGLVector3 extents = Extents();
   fVolume   = TMath::Abs(extents.X() * extents.Y() * extents.Z());
   fDiagonal = extents.Mag();
}

////////////////////////////////////////////////////////////////////////////////
/// Set a bounding box from provided 8 vertices

void TGLBoundingBox::Set(const TGLVertex3 vertex[8])
{
   for (UInt_t v = 0; v < 8; v++) {
      fVertex[v] = vertex[v];
   }
   // Could change cached volume/axes
   UpdateCache();
}

////////////////////////////////////////////////////////////////////////////////
/// Set a bounding box from provided 8 vertices

void TGLBoundingBox::Set(const Double_t vertex[8][3])
{
   for (UInt_t v = 0; v < 8; v++) {
      for (UInt_t a = 0; a < 3; a++) {
         fVertex[v][a] = vertex[v][a];
      }
   }
   // Could change cached volume/axes
   UpdateCache();
}

////////////////////////////////////////////////////////////////////////////////
/// Set a bounding box from vertices of other

void TGLBoundingBox::Set(const TGLBoundingBox & other)
{
   for (UInt_t v = 0; v < 8; v++) {
      fVertex[v].Set(other.fVertex[v]);
   }
   // Could change cached volume/axes
   UpdateCache();
}

////////////////////////////////////////////////////////////////////////////////
/// Set bounding box empty - all vertices at (0,0,0)

void TGLBoundingBox::SetEmpty()
{
   for (UInt_t v = 0; v < 8; v++) {
      fVertex[v].Fill(0.0);
   }
   // Could change cached volume/axes
   UpdateCache();
}

////////////////////////////////////////////////////////////////////////////////
/// Set ALIGNED box from two low/high vertices. Box axes are aligned with
/// global frame axes that vertices are specified in.

void TGLBoundingBox::SetAligned(const TGLVertex3 & lowVertex, const TGLVertex3 & highVertex)
{
   // lowVertex = vertex[0]
   // highVertex = vertex[6]
   //
   //    y
   //    |
   //    |
   //    |________x
   //   /  3-------2
   //  /  /|      /|
   // z  7-------6 |
   //    | 0-----|-1
   //    |/      |/
   //    4-------5
   //

   TGLVector3 diff = highVertex - lowVertex;
   if (diff.X() < 0.0 || diff.Y() < 0.0 || diff.Z() < 0.0) {
      Error("TGLBoundingBox::SetAligned", "low/high vertex range error");
   }
   fVertex[0] = lowVertex;
   fVertex[1] = lowVertex;  fVertex[1].X() += diff.X();
   fVertex[2] = lowVertex;  fVertex[2].X() += diff.X(); fVertex[2].Y() += diff.Y();
   fVertex[3] = lowVertex;  fVertex[3].Y() += diff.Y();
   fVertex[4] = highVertex; fVertex[4].X() -= diff.X(); fVertex[4].Y() -= diff.Y();
   fVertex[5] = highVertex; fVertex[5].Y() -= diff.Y();
   fVertex[6] = highVertex;
   fVertex[7] = highVertex; fVertex[7].X() -= diff.X();
   // Could change cached volume/axes
   UpdateCache();
}

////////////////////////////////////////////////////////////////////////////////
/// Set ALIGNED box from one or more points. Box axes are aligned with
/// global frame axes that points are specified in.

void TGLBoundingBox::SetAligned(UInt_t nbPnts, const Double_t * pnts)
{
   if (nbPnts < 1 || !pnts) {
      assert(false);
      return;
   }

   // Single point gives a zero volume BB
   TGLVertex3 low(pnts[0], pnts[1], pnts[2]);
   TGLVertex3 high(pnts[0], pnts[1], pnts[2]);

   for (UInt_t p = 1; p < nbPnts; p++) {
      for (UInt_t i = 0; i < 3; i++) {
         if (pnts[3*p + i] < low[i]) {
            low[i] = pnts[3*p + i] ;
         }
         if (pnts[3*p + i] > high[i]) {
            high[i] = pnts[3*p + i] ;
         }
      }
   }

   SetAligned(low, high);
}

////////////////////////////////////////////////////////////////////////////////
/// Expand current bbox so that it includes other's bbox.
/// This make the bbox axis-aligned.

void TGLBoundingBox::MergeAligned(const TGLBoundingBox & other)
{
   if (other.IsEmpty()) return;
   if (IsEmpty())
   {
      Set(other);
   }
   else
   {
      TGLVertex3 low (other.MinAAVertex());
      TGLVertex3 high(other.MaxAAVertex());

      low .Minimum(MinAAVertex());
      high.Maximum(MaxAAVertex());
      SetAligned(low, high);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Expand current bbox so that it includes the point.
/// This make the bbox axis-aligned.

void TGLBoundingBox::ExpandAligned(const TGLVertex3 & point)
{
   TGLVertex3 low (MinAAVertex());
   TGLVertex3 high(MaxAAVertex());

   low .Minimum(point);
   high.Maximum(point);

   SetAligned(low, high);
}

////////////////////////////////////////////////////////////////////////////////
/// Isotropically scale bounding box along it's LOCAL axes, preserving center

void TGLBoundingBox::Scale(Double_t factor)
{
   Scale(factor, factor, factor);
   // Could change cached volume/axes
   UpdateCache();
}

////////////////////////////////////////////////////////////////////////////////
/// Asymmetrically scale box along it's LOCAL x,y,z axes, preserving center

void TGLBoundingBox::Scale(Double_t xFactor, Double_t yFactor, Double_t zFactor)
{
   // Get x,y,z edges (non-normalised axis) and scale
   // them by factors
   const TGLVector3 xOffset = Axis(0, kFALSE)*(xFactor - 1.0) / 2.0;
   const TGLVector3 yOffset = Axis(1, kFALSE)*(yFactor - 1.0) / 2.0;
   const TGLVector3 zOffset = Axis(2, kFALSE)*(zFactor - 1.0) / 2.0;

   //    y
   //    |
   //    |
   //    |________x
   //   /  3-------2
   //  /  /|      /|
   // z  7-------6 |
   //    | 0-----|-1
   //    |/      |/
   //    4-------5
   //
   fVertex[0] += -xOffset - yOffset - zOffset;
   fVertex[1] +=  xOffset - yOffset - zOffset;
   fVertex[2] +=  xOffset + yOffset - zOffset;
   fVertex[3] += -xOffset + yOffset - zOffset;

   fVertex[4] += -xOffset - yOffset + zOffset;
   fVertex[5] +=  xOffset - yOffset + zOffset;
   fVertex[6] +=  xOffset + yOffset + zOffset;
   fVertex[7] += -xOffset + yOffset + zOffset;

   // Could change cached volume/axes
   UpdateCache();
}

////////////////////////////////////////////////////////////////////////////////
/// Translate all vertices by offset

void TGLBoundingBox::Translate(const TGLVector3 & offset)
{
   for (UInt_t v = 0; v < 8; v++) {
      fVertex[v] = fVertex[v] + offset;
   }

   // No cache change - volume and axes vectors remain same
}

////////////////////////////////////////////////////////////////////////////////
/// Transform all vertices with matrix.

void TGLBoundingBox::Transform(const TGLMatrix & matrix)
{
   for (UInt_t v = 0; v < 8; v++) {
      matrix.TransformVertex(fVertex[v]);
   }

   // Could change cached volume/axes
   UpdateCache();
}

////////////////////////////////////////////////////////////////////////////////
///return a vector of face vertices
///    y
///    |
///    |
///    |________x
///   /  3-------2
///  /  /|      /|
/// z  7-------6 |
///    | 0-----|-1
///    |/      |/
///    4-------5
///

const std::vector<UInt_t> & TGLBoundingBox::FaceVertices(EFace face) const
{
   static Bool_t init = kFALSE;
   static std::vector<UInt_t> faceIndexes[kFaceCount];
   if (!init) {
      // Low X - 7403
      faceIndexes[kFaceLowX].push_back(7);
      faceIndexes[kFaceLowX].push_back(4);
      faceIndexes[kFaceLowX].push_back(0);
      faceIndexes[kFaceLowX].push_back(3);
      // High X - 2156
      faceIndexes[kFaceHighX].push_back(2);
      faceIndexes[kFaceHighX].push_back(1);
      faceIndexes[kFaceHighX].push_back(5);
      faceIndexes[kFaceHighX].push_back(6);
      // Low Y - 5104
      faceIndexes[kFaceLowY].push_back(5);
      faceIndexes[kFaceLowY].push_back(1);
      faceIndexes[kFaceLowY].push_back(0);
      faceIndexes[kFaceLowY].push_back(4);
      // High Y - 2673
      faceIndexes[kFaceHighY].push_back(2);
      faceIndexes[kFaceHighY].push_back(6);
      faceIndexes[kFaceHighY].push_back(7);
      faceIndexes[kFaceHighY].push_back(3);
      // Low Z - 3012
      faceIndexes[kFaceLowZ].push_back(3);
      faceIndexes[kFaceLowZ].push_back(0);
      faceIndexes[kFaceLowZ].push_back(1);
      faceIndexes[kFaceLowZ].push_back(2);
      // High Z - 6547
      faceIndexes[kFaceHighZ].push_back(6);
      faceIndexes[kFaceHighZ].push_back(5);
      faceIndexes[kFaceHighZ].push_back(4);
      faceIndexes[kFaceHighZ].push_back(7);
      init= kTRUE;
   }
   return faceIndexes[face];
}

////////////////////////////////////////////////////////////////////////////////
/// Fill out supplied plane set vector with TGLPlane objects
/// representing six faces of box

void TGLBoundingBox::PlaneSet(TGLPlaneSet_t & planeSet) const
{
   assert(planeSet.empty());

   //    y
   //    |
   //    |
   //    |________x
   //   /  3-------2
   //  /  /|      /|
   // z  7-------6 |
   //    | 0-----|-1
   //    |/      |/
   //    4-------5
   //
   // Construct plane set using axis + vertices
   planeSet.push_back(TGLPlane( fAxesNorm[2], fVertex[4])); // Near
   planeSet.push_back(TGLPlane(-fAxesNorm[2], fVertex[0])); // Far
   planeSet.push_back(TGLPlane(-fAxesNorm[0], fVertex[0])); // Left
   planeSet.push_back(TGLPlane( fAxesNorm[0], fVertex[1])); // Right
   planeSet.push_back(TGLPlane(-fAxesNorm[1], fVertex[0])); // Bottom
   planeSet.push_back(TGLPlane( fAxesNorm[1], fVertex[3])); // Top
}

////////////////////////////////////////////////////////////////////////////////
/// Return the near-plane.

TGLPlane TGLBoundingBox::GetNearPlane() const
{
   return TGLPlane(fAxesNorm[2], fVertex[4]);
}

////////////////////////////////////////////////////////////////////////////////
/// Find overlap (Inside, Outside, Partial) of plane c.f. bounding box.

Rgl::EOverlap TGLBoundingBox::Overlap(const TGLPlane & plane) const
{
   using namespace Rgl;

   // First : cheap square approximation test. If distance of our
   // center to plane > our half extent length we are outside plane
   if (plane.DistanceTo(Center()) + (Extents().Mag()/2.0) < 0.0) {
      return kOutside;
   }

   // Second : test all 8 box vertices against plane
   Int_t verticesInsidePlane = 8;
   for (UInt_t v = 0; v < 8; v++) {
      if (plane.DistanceTo(fVertex[v]) < 0.0) {
         verticesInsidePlane--;
      }
   }

   if ( verticesInsidePlane == 0 ) {
      return kOutside;
   } else if ( verticesInsidePlane == 8 ) {
      return kInside;
   } else {
      return kPartial;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Find overlap (Inside, Outside, Partial) of other bounding box c.f. us.

Rgl::EOverlap TGLBoundingBox::Overlap(const TGLBoundingBox & other) const
{
   using namespace Rgl;

   // Simplify code with refs
   const TGLBoundingBox & a = *this;
   const TGLBoundingBox & b = other;

   TGLVector3 aHL = a.Extents() / 2.0; // Half length extents
   TGLVector3 bHL = b.Extents() / 2.0; // Half length extents

   // Following tests are greatly simplified
   // if we convert into our local frame

   // Find translation in parent frame
   TGLVector3 parentT = b.Center() - a.Center();

   // First: Do a simple & cheap sphere approximation containment test.
   // In many uses b will be completely contained by a and very much smaller
   // these cases short circuited here

   // We need the inner sphere for the container (box a) - radius = shortest box half length
   Double_t aSphereRadius = aHL[0] < aHL[1] ? aHL[0] : aHL[1];
   if (aHL[2] < aSphereRadius) {
      aSphereRadius = aHL[2];
   }
   // and the outer sphere for container (box b) - radius = box diagonal
   Double_t bSphereRadius = bHL.Mag();

   // If b sphere radius + translation mag is smaller than b sphere radius
   // b is complete contained by a
   if (bSphereRadius + parentT.Mag() < aSphereRadius) {
      return kInside;
   }

   // Second: Perform more expensive 15 separating axes test

   // Find translation in A's frame
   TGLVector3 aT(Dot(parentT, a.Axis(0)), Dot(parentT, a.Axis(1)), Dot(parentT, a.Axis(2)));

   // Find B's basis with respect to A's local frame
   // Get rotation matrix
   Double_t   roaT[3][3];
   UInt_t     i, k;
   for (i=0 ; i<3 ; i++) {
      for (k=0; k<3; k++) {
         roaT[i][k] = Dot(a.Axis(i), b.Axis(k));
         // Force very small components to zero to avoid rounding errors
         if (fabs(roaT[i][k]) < 1e-14) {
            roaT[i][k] = 0.0;
         }
      }
      // Normalise columns to avoid rounding errors
      Double_t norm = sqrt(roaT[i][0]*roaT[i][0] + roaT[i][1]*roaT[i][1] + roaT[i][2]*roaT[i][2]);
      roaT[i][0] /= norm; roaT[i][1] /= norm; roaT[i][2] /= norm;
   }

   // Perform separating axis test for all 15 potential
   // axes. If no separating axes found, the two boxes overlap.
   Double_t ra, rb, t;

   // A's 3 basis vectors
   for (i=0; i<3; i++) {
      ra = aHL[i];
      rb = bHL[0]*fabs(roaT[i][0]) + bHL[1]*fabs(roaT[i][1]) + bHL[2]*fabs(roaT[i][2]);
      t = fabs(aT[i]);
      if (t > ra + rb)
         return kOutside;
      else if (ra < t + rb)
         return kPartial;
   }

   // B's 3 basis vectors
   for (k=0; k<3; k++) {
      ra = aHL[0]*fabs(roaT[0][k]) + aHL[1]*fabs(roaT[1][k]) + aHL[2]*fabs(roaT[2][k]);
      rb = bHL[k];
      t = fabs(aT[0]*roaT[0][k] + aT[1]*roaT[1][k] + aT[2]*roaT[2][k]);
      if (t > ra + rb)
         return kOutside;
      else if (ra < t + rb)
         return kPartial;
   }

   // Now the 9 cross products

   // A0 x B0
   ra = aHL[1]*fabs(roaT[2][0]) + aHL[2]*fabs(roaT[1][0]);
   rb = bHL[1]*fabs(roaT[0][2]) + bHL[2]*fabs(roaT[0][1]);
   t = fabs(aT[2]*roaT[1][0] - aT[1]*roaT[2][0]);
   if (t > ra + rb)
      return kOutside;
   else if (ra < t + rb)
      return kPartial;

   // A0 x B1
   ra = aHL[1]*fabs(roaT[2][1]) + aHL[2]*fabs(roaT[1][1]);
   rb = bHL[0]*fabs(roaT[0][2]) + bHL[2]*fabs(roaT[0][0]);
   t = fabs(aT[2]*roaT[1][1] - aT[1]*roaT[2][1]);
   if (t > ra + rb)
      return kOutside;
   else if (ra < t + rb)
      return kPartial;

   // A0 x B2
   ra = aHL[1]*fabs(roaT[2][2]) + aHL[2]*fabs(roaT[1][2]);
   rb = bHL[0]*fabs(roaT[0][1]) + bHL[1]*fabs(roaT[0][0]);
   t = fabs(aT[2]*roaT[1][2] - aT[1]*roaT[2][2]);
   if (t > ra + rb)
      return kOutside;
   else if (ra < t + rb)
      return kPartial;

   // A1 x B0
   ra = aHL[0]*fabs(roaT[2][0]) + aHL[2]*fabs(roaT[0][0]);
   rb = bHL[1]*fabs(roaT[1][2]) + bHL[2]*fabs(roaT[1][1]);
   t = fabs(aT[0]*roaT[2][0] - aT[2]*roaT[0][0]);
   if (t > ra + rb)
      return kOutside;
   else if (ra < t + rb)
      return kPartial;

   // A1 x B1
   ra = aHL[0]*fabs(roaT[2][1]) + aHL[2]*fabs(roaT[0][1]);
   rb = bHL[0]*fabs(roaT[1][2]) + bHL[2]*fabs(roaT[1][0]);
   t = fabs(aT[0]*roaT[2][1] - aT[2]*roaT[0][1]);
   if (t > ra + rb)
      return kOutside;
   else if (ra < t + rb)
      return kPartial;

   // A1 x B2
   ra = aHL[0]*fabs(roaT[2][2]) + aHL[2]*fabs(roaT[0][2]);
   rb = bHL[0]*fabs(roaT[1][1]) + bHL[1]*fabs(roaT[1][0]);
   t = fabs(aT[0]*roaT[2][2] - aT[2]*roaT[0][2]);
   if (t > ra + rb)
      return kOutside;
   else if (ra < t + rb)
      return kPartial;

   // A2 x B0
   ra = aHL[0]*fabs(roaT[1][0]) + aHL[1]*fabs(roaT[0][0]);
   rb = bHL[1]*fabs(roaT[2][2]) + bHL[2]*fabs(roaT[2][1]);
   t = fabs(aT[1]*roaT[0][0] - aT[0]*roaT[1][0]);
   if (t > ra + rb)
      return kOutside;
   else if (ra < t + rb)
      return kPartial;

   // A2 x B1
   ra = aHL[0]*fabs(roaT[1][1]) + aHL[1]*fabs(roaT[0][1]);
   rb = bHL[0]*fabs(roaT[2][2]) + bHL[2]*fabs(roaT[2][0]);
   t = fabs(aT[1]*roaT[0][1] - aT[0]*roaT[1][1]);
   if (t > ra + rb)
      return kOutside;
   else if (ra < t + rb)
      return kPartial;

   // A2 x B2
   ra = aHL[0]*fabs(roaT[1][2]) + aHL[1]*fabs(roaT[0][2]);
   rb = bHL[0]*fabs(roaT[2][1]) + bHL[1]*fabs(roaT[2][0]);
   t = fabs(aT[1]*roaT[0][2] - aT[0]*roaT[1][2]);
   if (t > ra + rb)
      return kOutside;
   else if (ra < t + rb)
      return kPartial;

   // No separating axis - b is inside a
   return kInside;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the bounding box as either wireframe (default) of solid
/// using current GL color.

void TGLBoundingBox::Draw(Bool_t solid) const
{
   if (!solid) {
      glBegin(GL_LINE_LOOP);
      glVertex3dv(fVertex[0].CArr());
      glVertex3dv(fVertex[1].CArr());
      glVertex3dv(fVertex[2].CArr());
      glVertex3dv(fVertex[3].CArr());
      glVertex3dv(fVertex[7].CArr());
      glVertex3dv(fVertex[6].CArr());
      glVertex3dv(fVertex[5].CArr());
      glVertex3dv(fVertex[4].CArr());
      glEnd();
      glBegin(GL_LINES);
      glVertex3dv(fVertex[1].CArr());
      glVertex3dv(fVertex[5].CArr());
      glVertex3dv(fVertex[2].CArr());
      glVertex3dv(fVertex[6].CArr());
      glVertex3dv(fVertex[0].CArr());
      glVertex3dv(fVertex[3].CArr());
      glVertex3dv(fVertex[4].CArr());
      glVertex3dv(fVertex[7].CArr());
      glEnd();
   } else {
   //    y
   //    |
   //    |
   //    |________x
   //   /  3-------2
   //  /  /|      /|
   // z  7-------6 |
   //    | 0-----|-1
   //    |/      |/
   //    4-------5
      // Clockwise winding
      glBegin(GL_QUADS);
      // Near
      glNormal3d ( fAxesNorm[2].X(),  fAxesNorm[2].Y(),  fAxesNorm[2].Z());
      glVertex3dv(fVertex[4].CArr());
      glVertex3dv(fVertex[7].CArr());
      glVertex3dv(fVertex[6].CArr());
      glVertex3dv(fVertex[5].CArr());
      // Far
      glNormal3d (-fAxesNorm[2].X(), -fAxesNorm[2].Y(), -fAxesNorm[2].Z());
      glVertex3dv(fVertex[0].CArr());
      glVertex3dv(fVertex[1].CArr());
      glVertex3dv(fVertex[2].CArr());
      glVertex3dv(fVertex[3].CArr());
      // Left
      glNormal3d (-fAxesNorm[0].X(), -fAxesNorm[0].Y(), -fAxesNorm[0].Z());
      glVertex3dv(fVertex[0].CArr());
      glVertex3dv(fVertex[3].CArr());
      glVertex3dv(fVertex[7].CArr());
      glVertex3dv(fVertex[4].CArr());
      // Right
      glNormal3d ( fAxesNorm[0].X(),  fAxesNorm[0].Y(),  fAxesNorm[0].Z());
      glVertex3dv(fVertex[6].CArr());
      glVertex3dv(fVertex[2].CArr());
      glVertex3dv(fVertex[1].CArr());
      glVertex3dv(fVertex[5].CArr());
      // Top
      glNormal3d ( fAxesNorm[1].X(),  fAxesNorm[1].Y(),  fAxesNorm[1].Z());
      glVertex3dv(fVertex[3].CArr());
      glVertex3dv(fVertex[2].CArr());
      glVertex3dv(fVertex[6].CArr());
      glVertex3dv(fVertex[7].CArr());
      // Bottom
      glNormal3d (-fAxesNorm[1].X(), -fAxesNorm[1].Y(), -fAxesNorm[1].Z());
      glVertex3dv(fVertex[4].CArr());
      glVertex3dv(fVertex[5].CArr());
      glVertex3dv(fVertex[1].CArr());
      glVertex3dv(fVertex[0].CArr());

      glEnd();
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Find minimum vertex value for axis of index X(0), Y(1), Z(2)

Double_t TGLBoundingBox::Min(UInt_t index) const
{
   Double_t min = fVertex[0][index];
   for (UInt_t v = 1; v < 8; v++) {
      if (fVertex[v][index] < min) {
         min = fVertex[v][index];
      }
   }
   return min;
}

////////////////////////////////////////////////////////////////////////////////
/// Find maximum vertex value for axis of index X(0), Y(1), Z(2)

Double_t TGLBoundingBox::Max(UInt_t index) const
{
   Double_t max = fVertex[0][index];
   for (UInt_t v = 1; v < 8; v++) {
      if (fVertex[v][index] > max) {
         max = fVertex[v][index];
      }
   }
   return max;
}

////////////////////////////////////////////////////////////////////////////////
/// Find minimum vertex values.

TGLVertex3 TGLBoundingBox::MinAAVertex() const
{
   return TGLVertex3(Min(0), Min(1), Min(2));
}

////////////////////////////////////////////////////////////////////////////////
/// Find maximum vertex values.

TGLVertex3 TGLBoundingBox::MaxAAVertex() const
{
   return TGLVertex3(Max(0), Max(1), Max(2));
}

////////////////////////////////////////////////////////////////////////////////
/// Output to std::cout the vertices, center and volume of box

void TGLBoundingBox::Dump() const
{
   for (UInt_t i = 0; i<8; i++) {
      std::cout << "[" << i << "] (" << fVertex[i].X() << "," << fVertex[i].Y() << "," << fVertex[i].Z() << ")" << std::endl;
   }
   std::cout << "Center:  ";   Center().Dump();
   std::cout << "Extents: ";   Extents().Dump();
   std::cout << "Volume:  " << Volume() << std::endl;
}

