// @(#)root/geom:$Id$
// Author: Andrei Gheata   26/09/05

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#include <iostream>

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoScaledShape.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"

/** \class TGeoScaledShape
\ingroup Geometry_classes

A shape scaled by a TGeoScale transformation
\image html geom_scaledshape.png
*/

ClassImp(TGeoScaledShape);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoScaledShape::TGeoScaledShape()
{
   fShape = 0;
   fScale = 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor

TGeoScaledShape::TGeoScaledShape(const char *name, TGeoShape *shape, TGeoScale *scale)
                :TGeoBBox(name,0,0,0)
{
   fShape = shape;
   fScale = scale;
   if (!fScale->IsRegistered()) fScale->RegisterYourself();
   TGeoScaledShape::ComputeBBox();

}

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TGeoScaledShape::TGeoScaledShape(TGeoShape *shape, TGeoScale *scale)
{
   fShape = shape;
   fScale = scale;
   if (!fScale->IsRegistered()) fScale->RegisterYourself();
   TGeoScaledShape::ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoScaledShape::~TGeoScaledShape()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Computes capacity of this shape [length^3]

Double_t TGeoScaledShape::Capacity() const
{
   Double_t capacity = fShape->Capacity();
   const Double_t *scale = fScale->GetScale();
   capacity *= scale[0]*scale[1]*scale[2];
   return capacity;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute bounding box of the scaled shape

void TGeoScaledShape::ComputeBBox()
{
   if (!fShape) {
      Error("ComputeBBox", "Scaled shape %s without shape", GetName());
      return;
   }
   if (fShape->IsAssembly()) fShape->ComputeBBox();
   TGeoBBox *box = (TGeoBBox*)fShape;
   const Double_t *orig = box->GetOrigin();
   Double_t point[3], master[3];
   point[0] = box->GetDX();
   point[1] = box->GetDY();
   point[2] = box->GetDZ();

   fScale->LocalToMaster(orig, fOrigin);
   fScale->LocalToMaster(point, master);
   fDX = TMath::Abs(master[0]);
   fDY = TMath::Abs(master[1]);
   fDZ = TMath::Abs(master[2]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute normal to closest surface from POINT.

void TGeoScaledShape::ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm)
{
   Double_t local[3], ldir[3], lnorm[3];
   fScale->MasterToLocal(point,local);
   fScale->MasterToLocalVect(dir,ldir);
   TGeoMatrix::Normalize(ldir);
   fShape->ComputeNormal(local,ldir,lnorm);
//   fScale->LocalToMasterVect(lnorm, norm);
   fScale->MasterToLocalVect(lnorm, norm);
   TGeoMatrix::Normalize(norm);
}

////////////////////////////////////////////////////////////////////////////////
/// Test if point is inside the scaled shape

Bool_t TGeoScaledShape::Contains(const Double_t *point) const
{
   Double_t local[3];
   fScale->MasterToLocal(point,local);
   return fShape->Contains(local);
}

////////////////////////////////////////////////////////////////////////////////
/// compute closest distance from point px,py to each vertex. Should not be called.

Int_t TGeoScaledShape::DistancetoPrimitive(Int_t px, Int_t py)
{
   Int_t n = fShape->GetNmeshVertices();
   return ShapeDistancetoPrimitive(n, px, py);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from inside point to surface of the scaled shape.

Double_t TGeoScaledShape::DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
{
   Double_t local[3], ldir[3];
   Double_t lstep;
   fScale->MasterToLocal(point,local);
   lstep = fScale->MasterToLocal(step, dir);
   fScale->MasterToLocalVect(dir,ldir);
   TGeoMatrix::Normalize(ldir);
   Double_t dist = fShape->DistFromInside(local,ldir, iact, lstep, safe);
   if (iact<3 && safe) *safe = fScale->LocalToMaster(*safe);
   dist = fScale->LocalToMaster(dist, ldir);
   return dist;
}


////////////////////////////////////////////////////////////////////////////////
/// Compute distance from outside point to surface of the scaled shape.

Double_t TGeoScaledShape::DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
{
   Double_t local[3], ldir[3];
   Double_t lstep;
//   printf("DistFromOutside(%f,%f,%f,  %f,%f,%f)\n", point[0], point[1], point[2], dir[0], dir[1],dir[2]);
   fScale->MasterToLocal(point,local);
//   printf("local: %f,%f,%f\n", local[0],local[1], local[2]);
   lstep = fScale->MasterToLocal(step, dir);
   fScale->MasterToLocalVect(dir,ldir);
   TGeoMatrix::Normalize(ldir);
//   printf("localdir: %f,%f,%f\n",ldir[0],ldir[1],ldir[2]);
   Double_t dist = fShape->DistFromOutside(local,ldir, iact, lstep, safe);
//   printf("local distance: %f\n", dist);
   if (safe) *safe = fScale->LocalToMaster(*safe);
   dist = fScale->LocalToMaster(dist, ldir);
//   printf("converted distance: %f\n",dist);
   return dist;
}

////////////////////////////////////////////////////////////////////////////////
/// Cannot divide assemblies.

TGeoVolume *TGeoScaledShape::Divide(TGeoVolume * /*voldiv*/, const char *divname, Int_t /*iaxis*/, Int_t /*ndiv*/,
                             Double_t /*start*/, Double_t /*step*/)
{
   Error("Divide", "Scaled shapes cannot be divided. Division volume %s not created", divname);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills a static 3D buffer and returns a reference.

const TBuffer3D & TGeoScaledShape::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
{
   TBuffer3D &buffer = (TBuffer3D &)fShape->GetBuffer3D(reqSections, localFrame);

//   TGeoBBox::FillBuffer3D(buffer, reqSections, localFrame);
   Double_t halfLengths[3] = { fDX, fDY, fDZ };
   buffer.SetAABoundingBox(fOrigin, halfLengths);
   if (!buffer.fLocalFrame) {
      TransformPoints(buffer.fBBVertex[0], 8);
   }

   if ((reqSections & TBuffer3D::kRaw) && buffer.SectionsValid(TBuffer3D::kRawSizes)) {
      SetPoints(buffer.fPnts);
      if (!buffer.fLocalFrame) {
         TransformPoints(buffer.fPnts, buffer.NbPnts());
      }
   }

   return buffer;
}
////////////////////////////////////////////////////////////////////////////////
/// in case shape has some negative parameters, these has to be computed
/// in order to fit the mother

TGeoShape *TGeoScaledShape::GetMakeRuntimeShape(TGeoShape * /*mother*/, TGeoMatrix * /*mat*/) const
{
   Error("GetMakeRuntimeShape", "Scaled shapes cannot be parametrized.");
   return NULL;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoScaledShape::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   fShape->GetMeshNumbers(nvert, nsegs, npols);
}

////////////////////////////////////////////////////////////////////////////////
/// print shape parameters

void TGeoScaledShape::InspectShape() const
{
   printf("*** Shape %s: TGeoScaledShape ***\n", GetName());
   fScale->Print();
   fShape->InspectShape();
   TGeoBBox::InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true if the scaled shape is an assembly.

Bool_t TGeoScaledShape::IsAssembly() const
{
   return fShape->IsAssembly();
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the scale transformation is a reflection.

Bool_t TGeoScaledShape::IsReflected() const
{
   return fScale->IsReflection();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a TBuffer3D describing *this* shape.
/// Coordinates are in local reference frame.

TBuffer3D *TGeoScaledShape::MakeBuffer3D() const
{
   TBuffer3D *buff = fShape->MakeBuffer3D();
   if (buff) SetPoints(buff->fPnts);
   return buff;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a scaled shape starting from a non-scaled one.

TGeoShape *TGeoScaledShape::MakeScaledShape(const char *name, TGeoShape *shape, TGeoScale *scale)
{
   TGeoShape *old_shape = shape;
   TGeoShape *new_shape = nullptr;
   if (shape->IsA() == TGeoScaledShape::Class()) {
      TGeoScaledShape *sshape = (TGeoScaledShape*)shape;
      TGeoScale *old_scale = sshape->GetScale();
      old_shape = sshape->GetShape();
      scale->SetScale(scale->GetScale()[0]*old_scale->GetScale()[0],
                      scale->GetScale()[1]*old_scale->GetScale()[1],
                      scale->GetScale()[2]*old_scale->GetScale()[2]);
   }
   if (old_shape->IsAssembly()) {
      // The shape is owned by the assembly, so make sure it does not register to TGeoManager
      new_shape = new TGeoScaledShape(old_shape, scale);
      new_shape->SetName(name);
   } else {
      new_shape = new TGeoScaledShape(name, old_shape, scale);
   }
   return new_shape;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill TBuffer3D structure for segments and polygons.

void TGeoScaledShape::SetSegsAndPols(TBuffer3D &buff) const
{
   fShape->SetSegsAndPols(buff);
}

////////////////////////////////////////////////////////////////////////////////
/// computes the closest distance from given point to this shape, according
/// to option. The matching point on the shape is stored in spoint.

Double_t TGeoScaledShape::Safety(const Double_t *point, Bool_t in) const
{
   Double_t local[3];
   fScale->MasterToLocal(point,local);
   Double_t safe = fShape->Safety(local,in);
   safe = fScale->LocalToMaster(safe);
   return safe;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoScaledShape::SavePrimitive(std::ostream &out, Option_t *option)
{
   if (TObject::TestBit(kGeoSavePrimitive)) return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   if (!fShape || !fScale) {
      out << "##### Invalid shape or scale !. Aborting. #####" << std::endl;
      return;
   }
   fShape->SavePrimitive(out, option);
   TString sname = fShape->GetPointerName();
   const Double_t *sc = fScale->GetScale();
   out << "   // Scale factor:" << std::endl;
   out << "   TGeoScale *pScale = new TGeoScale(\"" << fScale->GetName()
       << "\"," << sc[0] << "," << sc[1] << "," << sc[2] << ");" << std::endl;
   out << "   TGeoScaledShape *" << GetPointerName() << " = new TGeoScaledShape(\""
       << GetName() << "\"," << sname << ", pScale);" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Mesh points for scaled shapes.

void TGeoScaledShape::SetPoints(Double_t *points) const
{
   Int_t npts = fShape->GetNmeshVertices();
   fShape->SetPoints(points);
   Double_t master[3];
   for (Int_t i=0; i<npts; i++) {
      fScale->LocalToMaster(&points[3*i], master);
      memcpy(&points[3*i], master, 3*sizeof(Double_t));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Mesh points for scaled shapes.

void TGeoScaledShape::SetPoints(Float_t *points) const
{
   Int_t npts = fShape->GetNmeshVertices();
   fShape->SetPoints(points);
   Double_t master[3];
   Double_t local[3];
   Int_t index;
   for (Int_t i=0; i<npts; i++) {
      index = 3*i;
      local[0] = points[index];
      local[1] = points[index+1];
      local[2] = points[index+2];
      fScale->LocalToMaster(local, master);
      points[index] = master[0];
      points[index+1] = master[1];
      points[index+2] = master[2];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check the inside status for each of the points in the array.
/// Input: Array of point coordinates + vector size
/// Output: Array of Booleans for the inside of each point

void TGeoScaledShape::Contains_v(const Double_t *points, Bool_t *inside, Int_t vecsize) const
{
   for (Int_t i=0; i<vecsize; i++) inside[i] = Contains(&points[3*i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the normal for an array o points so that norm.dot.dir is positive
/// Input: Arrays of point coordinates and directions + vector size
/// Output: Array of normal directions

void TGeoScaledShape::ComputeNormal_v(const Double_t *points, const Double_t *dirs, Double_t *norms, Int_t vecsize)
{
   for (Int_t i=0; i<vecsize; i++) ComputeNormal(&points[3*i], &dirs[3*i], &norms[3*i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoScaledShape::DistFromInside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromInside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoScaledShape::DistFromOutside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromOutside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute safe distance from each of the points in the input array.
/// Input: Array of point coordinates, array of statuses for these points, size of the arrays
/// Output: Safety values

void TGeoScaledShape::Safety_v(const Double_t *points, const Bool_t *inside, Double_t *safe, Int_t vecsize) const
{
   for (Int_t i=0; i<vecsize; i++) safe[i] = Safety(&points[3*i], inside[i]);
}
