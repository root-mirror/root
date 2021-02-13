// @(#):$Id$
// Author: Andrei Gheata   30/05/02
// TGeoBoolNode::Contains and parser implemented by Mihaela Gheata


/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#include "TGeoBoolNode.h"

#include <iostream>

#include "TVirtualPad.h"
#include "TVirtualViewer3D.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"
#include "TGeoCompositeShape.h"
#include "TGeoMatrix.h"
#include "TGeoManager.h"

/** \class TGeoBoolNode
\ingroup Geometry_classes

Base class for Boolean operations between two shapes.

A Boolean node describes a Boolean operation between 'left' and 'right'
shapes positioned with respect to an ARBITRARY reference frame. The boolean
node is referenced by a mother composite shape and its shape components may
be primitive but also composite shapes. The later situation leads to a binary
tree hierarchy. When the parent composite shape is used to create a volume,
the reference frame of the volume is chosen to match the frame in which
node shape components were defined.

The positioned shape components may or may not be disjoint. The specific
implementations for Boolean nodes are:

  - TGeoUnion - representing the Boolean  union of two positioned shapes
  - TGeoSubtraction - representing the Boolean subtraction of two positioned shapes
  - TGeoIntersection - representing the Boolean intersection of two positioned shapes
*/

ClassImp(TGeoBoolNode);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGeoBoolNode::ThreadData_t::ThreadData_t() :
   fSelected(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGeoBoolNode::ThreadData_t::~ThreadData_t()
{
}

////////////////////////////////////////////////////////////////////////////////

TGeoBoolNode::ThreadData_t& TGeoBoolNode::GetThreadData() const
{
   Int_t tid = TGeoManager::ThreadId();
/*
   std::lock_guard<std::mutex> guard(fMutex);
   if (tid >= fThreadSize) {
      Error("GetThreadData", "Thread id=%d bigger than maximum declared thread number %d. \nUse TGeoManager::SetMaxThreads properly !!!",
             tid, fThreadSize);
   }
   if (tid >= fThreadSize)
   {
      fThreadData.resize(tid + 1);
      fThreadSize = tid + 1;
   }
   if (fThreadData[tid] == 0)
   {
   if (fThreadData[tid] == 0)
      fThreadData[tid] = new ThreadData_t;
   }
*/
   return *fThreadData[tid];
}

////////////////////////////////////////////////////////////////////////////////

void TGeoBoolNode::ClearThreadData() const
{
   std::lock_guard<std::mutex> guard(fMutex);
   std::vector<ThreadData_t*>::iterator i = fThreadData.begin();
   while (i != fThreadData.end())
   {
      delete *i;
      ++i;
   }
   fThreadData.clear();
   fThreadSize = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create thread data for n threads max.

void TGeoBoolNode::CreateThreadData(Int_t nthreads)
{
   std::lock_guard<std::mutex> guard(fMutex);
   fThreadData.resize(nthreads);
   fThreadSize = nthreads;
   for (Int_t tid=0; tid<nthreads; tid++) {
      if (fThreadData[tid] == 0) {
         fThreadData[tid] = new ThreadData_t;
      }
   }
   // Propagate to components
   if (fLeft)  fLeft->CreateThreadData(nthreads);
   if (fRight) fRight->CreateThreadData(nthreads);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the selected branch.

void TGeoBoolNode::SetSelected(Int_t sel)
{
   GetThreadData().fSelected = sel;
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoBoolNode::TGeoBoolNode()
{
   fLeft     = nullptr;
   fRight    = nullptr;
   fLeftMat  = nullptr;
   fRightMat = nullptr;
   fNpoints  = 0;
   fPoints   = nullptr;
   fThreadSize = 0;
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor called by TGeoCompositeShape providing 2 subexpressions for the 2 branches.

TGeoBoolNode::TGeoBoolNode(const char *expr1, const char *expr2)
{
   fLeft     = nullptr;
   fRight    = nullptr;
   fLeftMat  = nullptr;
   fRightMat = nullptr;
   fNpoints  = 0;
   fPoints   = nullptr;
   fThreadSize = 0;
   CreateThreadData(1);
   if (!MakeBranch(expr1, kTRUE)) {
      return;
   }
   if (!MakeBranch(expr2, kFALSE)) {
      return;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor providing left and right shapes and matrices (in the Boolean operation).

TGeoBoolNode::TGeoBoolNode(TGeoShape *left, TGeoShape *right, TGeoMatrix *lmat, TGeoMatrix *rmat)
{
   fLeft = left;
   fRight = right;
   fLeftMat = lmat;
   fNpoints  = 0;
   fPoints   = nullptr;
   fThreadSize = 0;
   CreateThreadData(1);
   if (!fLeftMat) fLeftMat = gGeoIdentity;
   else fLeftMat->RegisterYourself();
   fRightMat = rmat;
   if (!fRightMat) fRightMat = gGeoIdentity;
   else fRightMat->RegisterYourself();
   if (!fLeft) {
      Error("ctor", "left shape is NULL");
      return;
   }
   if (!fRight) {
      Error("ctor", "right shape is NULL");
      return;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.
/// --- deletion of components handled by TGeoManager class.

TGeoBoolNode::~TGeoBoolNode()
{
   if (fPoints) delete [] fPoints;
   ClearThreadData();
}

////////////////////////////////////////////////////////////////////////////////
/// Set fPoints array

void TGeoBoolNode::AssignPoints(Int_t npoints, Double_t *points)
{
   if (fPoints) {
      delete [] fPoints;
      fPoints = nullptr;
      fNpoints = 0;
   }
   if (points) {
      fNpoints = npoints;
      fPoints = new Double_t[3*fNpoints];
      memcpy(fPoints, points, 3*fNpoints*sizeof(Double_t));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Expands the boolean expression either on left or right branch, creating
/// component elements (composite shapes and boolean nodes). Returns true on success.

Bool_t TGeoBoolNode::MakeBranch(const char *expr, Bool_t left)
{
   TString sleft, sright, stransf;
   Int_t boolop = TGeoManager::Parse(expr, sleft, sright, stransf);
   if (boolop<0) {
      Error("MakeBranch", "invalid expression");
      return kFALSE;
   }
   TGeoShape *shape = nullptr;
   TGeoMatrix *mat;
   TString newshape;

   if (stransf.Length() == 0) {
      mat = gGeoIdentity;
   } else {
      mat = (TGeoMatrix*)gGeoManager->GetListOfMatrices()->FindObject(stransf.Data());
   }
   if (!mat) {
      Error("MakeBranch", "transformation %s not found", stransf.Data());
      return kFALSE;
   }
   switch (boolop) {
      case 0:
         // elementary shape
         shape = (TGeoShape*)gGeoManager->GetListOfShapes()->FindObject(sleft.Data());
         if (!shape) {
            Error("MakeBranch", "shape %s not found", sleft.Data());
            return kFALSE;
         }
         break;
      case 1:
         // composite shape - union
         newshape = sleft;
         newshape += "+";
         newshape += sright;
         shape = new TGeoCompositeShape(newshape.Data());
         break;
      case 2:
         // composite shape - difference
         newshape = sleft;
         newshape += "-";
         newshape += sright;
         shape = new TGeoCompositeShape(newshape.Data());
         break;
      case 3:
         // composite shape - intersection
         newshape = sleft;
         newshape += "*";
         newshape += sright;
         shape = new TGeoCompositeShape(newshape.Data());
         break;
   }
   if (boolop && (!shape || !shape->IsValid())) {
      Error("MakeBranch", "Shape %s not valid", newshape.Data());
      if (shape) delete shape;
      return kFALSE;
   }
   if (left) {
      fLeft = shape;
      fLeftMat = mat;
   } else {
      fRight = shape;
      fRightMat = mat;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Special schema for feeding the 3D buffers to the painter client.

void TGeoBoolNode::Paint(Option_t * option)
{
   TVirtualViewer3D * viewer = gPad->GetViewer3D();
   if (!viewer) return;

   // Components of composite shape hierarchies for local frame viewers are painted
   // in coordinate frame of the top level composite shape. So we force
   // conversion to this.  See TGeoPainter::PaintNode for loading of GLMatrix.
   Bool_t localFrame = kFALSE; //viewer->PreferLocalFrame();

   TGeoHMatrix *glmat = (TGeoHMatrix*)TGeoShape::GetTransform();
   TGeoHMatrix mat;
   mat = glmat; // keep a copy

   // Now perform fetch and add of the two components buffers.
   // Note we assume that composite shapes are always completely added
   // so don't bother to get addDaughters flag from viewer->AddObject()

   // Setup matrix and fetch/add the left component buffer
   glmat->Multiply(fLeftMat);
   //fLeft->Paint(option);
   if (TGeoCompositeShape *left = dynamic_cast<TGeoCompositeShape *>(fLeft)) left->PaintComposite(option);
   else {
      const TBuffer3D & leftBuffer = fLeft->GetBuffer3D(TBuffer3D::kAll, localFrame);
      viewer->AddObject(leftBuffer);
   }

   // Setup matrix and fetch/add the right component buffer
   *glmat = &mat;
   glmat->Multiply(fRightMat);
   //fRight->Paint(option);
   if (TGeoCompositeShape *right = dynamic_cast<TGeoCompositeShape *>(fRight)) right->PaintComposite(option);
   else {
      const TBuffer3D & rightBuffer = fRight->GetBuffer3D(TBuffer3D::kAll, localFrame);
      viewer->AddObject(rightBuffer);
   }

   *glmat = &mat;
}

////////////////////////////////////////////////////////////////////////////////
/// Register all matrices of the boolean node and descendents.

void TGeoBoolNode::RegisterMatrices()
{
   if (!fLeftMat->IsIdentity()) fLeftMat->RegisterYourself();
   if (!fRightMat->IsIdentity()) fRightMat->RegisterYourself();
   if (fLeft->IsComposite()) ((TGeoCompositeShape*)fLeft)->GetBoolNode()->RegisterMatrices();
   if (fRight->IsComposite()) ((TGeoCompositeShape*)fRight)->GetBoolNode()->RegisterMatrices();
}

////////////////////////////////////////////////////////////////////////////////
/// Replace one of the matrices. Does not work with TGeoIdentity. Returns true
/// if replacement was successful.

Bool_t TGeoBoolNode::ReplaceMatrix(TGeoMatrix *mat, TGeoMatrix *newmat)
{
   if (mat==gGeoIdentity || newmat==gGeoIdentity) {
      Error("ReplaceMatrix", "Matrices should not be gGeoIdentity. Use default matrix constructor to represent identities.");
      return kFALSE;
   }
   if (!mat || !newmat) {
      Error("ReplaceMatrix", "Matrices should not be null pointers.");
      return kFALSE;
   }
   Bool_t replaced = kFALSE;
   if (fLeftMat == mat) {
      fLeftMat  = newmat;
      replaced = kTRUE;
   }
   if (fRightMat == mat) {
      fRightMat  = newmat;
      replaced = kTRUE;
   }
   return replaced;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoBoolNode::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   fLeft->SavePrimitive(out,option);
   fRight->SavePrimitive(out,option);
   if (!fLeftMat->IsIdentity()) {
      fLeftMat->RegisterYourself();
      fLeftMat->SavePrimitive(out,option);
   }
   if (!fRightMat->IsIdentity()) {
      fRightMat->RegisterYourself();
      fRightMat->SavePrimitive(out,option);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill buffer with shape vertices.

void TGeoBoolNode::SetPoints(Double_t *points) const
{
   TGeoBoolNode *bn = (TGeoBoolNode*)this;
   Int_t npoints = bn->GetNpoints();
   memcpy(points, fPoints, 3*npoints*sizeof(Double_t));
}

////////////////////////////////////////////////////////////////////////////////
/// Fill buffer with shape vertices.

void TGeoBoolNode::SetPoints(Float_t *points) const
{
   TGeoBoolNode *bn = (TGeoBoolNode*)this;
   Int_t npoints = bn->GetNpoints();
   for (Int_t i=0; i<3*npoints; i++) points[i] = fPoints[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Register size of this 3D object

void TGeoBoolNode::Sizeof3D() const
{
   fLeft->Sizeof3D();
   fRight->Sizeof3D();
}

ClassImp(TGeoUnion);

////////////////////////////////////////////////////////////////////////////////
/// Make a clone of this. Pointers are preserved.

TGeoBoolNode *TGeoUnion::MakeClone() const
{
   return new TGeoUnion(fLeft, fRight, fLeftMat, fRightMat);
}

////////////////////////////////////////////////////////////////////////////////
/// Paint method.

void TGeoUnion::Paint(Option_t *option)
{
   TVirtualViewer3D *viewer = gPad->GetViewer3D();

   if (!viewer) {
      Error("Paint", "gPad->GetViewer3D() returned 0, cannot work with composite!\n");
      return;
   }

   viewer->AddCompositeOp(TBuffer3D::kCSUnion);

   TGeoBoolNode::Paint(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoUnion::TGeoUnion()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TGeoUnion::TGeoUnion(const char *expr1, const char *expr2)
          :TGeoBoolNode(expr1, expr2)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor providing pointers to components

TGeoUnion::TGeoUnion(TGeoShape *left, TGeoShape *right, TGeoMatrix *lmat, TGeoMatrix *rmat)
          :TGeoBoolNode(left,right,lmat,rmat)
{
   if (left->TestShapeBit(TGeoShape::kGeoHalfSpace) || right->TestShapeBit(TGeoShape::kGeoHalfSpace)) {
      Fatal("TGeoUnion", "Unions with a half-space (%s + %s) not allowed", left->GetName(), right->GetName());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor
/// --- deletion of components handled by TGeoManager class.

TGeoUnion::~TGeoUnion()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Compute bounding box corresponding to a union of two shapes.

void TGeoUnion::ComputeBBox(Double_t &dx, Double_t &dy, Double_t &dz, Double_t *origin)
{
   if (((TGeoBBox*)fLeft)->IsNullBox()) fLeft->ComputeBBox();
   if (((TGeoBBox*)fRight)->IsNullBox()) fRight->ComputeBBox();
   Double_t vert[48];
   Double_t pt[3];
   Int_t i;
   Double_t xmin, xmax, ymin, ymax, zmin, zmax;
   xmin = ymin = zmin = TGeoShape::Big();
   xmax = ymax = zmax = -TGeoShape::Big();
   ((TGeoBBox*)fLeft)->SetBoxPoints(&vert[0]);
   ((TGeoBBox*)fRight)->SetBoxPoints(&vert[24]);
   for (i=0; i<8; i++) {
      fLeftMat->LocalToMaster(&vert[3*i], &pt[0]);
      if (pt[0]<xmin) xmin=pt[0];
      if (pt[0]>xmax) xmax=pt[0];
      if (pt[1]<ymin) ymin=pt[1];
      if (pt[1]>ymax) ymax=pt[1];
      if (pt[2]<zmin) zmin=pt[2];
      if (pt[2]>zmax) zmax=pt[2];
   }
   for (i=8; i<16; i++) {
      fRightMat->LocalToMaster(&vert[3*i], &pt[0]);
      if (pt[0]<xmin) xmin=pt[0];
      if (pt[0]>xmax) xmax=pt[0];
      if (pt[1]<ymin) ymin=pt[1];
      if (pt[1]>ymax) ymax=pt[1];
      if (pt[2]<zmin) zmin=pt[2];
      if (pt[2]>zmax) zmax=pt[2];
   }
   dx = 0.5*(xmax-xmin);
   origin[0] = 0.5*(xmin+xmax);
   dy = 0.5*(ymax-ymin);
   origin[1] = 0.5*(ymin+ymax);
   dz = 0.5*(zmax-zmin);
   origin[2] = 0.5*(zmin+zmax);
}

////////////////////////////////////////////////////////////////////////////////
/// Find if a union of two shapes contains a given point

Bool_t TGeoUnion::Contains(const Double_t *point) const
{
   Double_t local[3];
   fLeftMat->MasterToLocal(point, &local[0]);
   Bool_t inside = fLeft->Contains(&local[0]);
   if (inside) return kTRUE;
   fRightMat->MasterToLocal(point, &local[0]);
   inside = fRight->Contains(&local[0]);
   return inside;
}

////////////////////////////////////////////////////////////////////////////////
/// Normal computation in POINT. The orientation is chosen so that DIR.dot.NORM>0.

void TGeoUnion::ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm)
{
   ThreadData_t& td = GetThreadData();
   norm[0] = norm[1] = 0.;
   norm[2] = 1.;
   Double_t local[3];
   Double_t ldir[3], lnorm[3];
   if (td.fSelected == 1) {
      fLeftMat->MasterToLocal(point, local);
      fLeftMat->MasterToLocalVect(dir, ldir);
      fLeft->ComputeNormal(local,ldir,lnorm);
      fLeftMat->LocalToMasterVect(lnorm, norm);
      return;
   }
   if (td.fSelected == 2) {
      fRightMat->MasterToLocal(point, local);
      fRightMat->MasterToLocalVect(dir, ldir);
      fRight->ComputeNormal(local,ldir,lnorm);
      fRightMat->LocalToMasterVect(lnorm, norm);
      return;
   }
   fLeftMat->MasterToLocal(point, local);
   if (fLeft->Contains(local)) {
      fLeftMat->MasterToLocalVect(dir, ldir);
      fLeft->ComputeNormal(local,ldir,lnorm);
      fLeftMat->LocalToMasterVect(lnorm, norm);
      return;
   }
   fRightMat->MasterToLocal(point, local);
   if (fRight->Contains(local)) {
      fRightMat->MasterToLocalVect(dir, ldir);
      fRight->ComputeNormal(local,ldir,lnorm);
      fRightMat->LocalToMasterVect(lnorm, norm);
      return;
   }
   // Propagate forward/backward to see which of the components is intersected first
   local[0] = point[0] + 1E-5*dir[0];
   local[1] = point[1] + 1E-5*dir[1];
   local[2] = point[2] + 1E-5*dir[2];

   if (!Contains(local)) {
      local[0] = point[0] - 1E-5*dir[0];
      local[1] = point[1] - 1E-5*dir[1];
      local[2] = point[2] - 1E-5*dir[2];
      if (!Contains(local)) return;
   }
   ComputeNormal(local,dir,norm);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute minimum distance to shape vertices.

Int_t TGeoUnion::DistanceToPrimitive(Int_t /*px*/, Int_t /*py*/)
{
   return 9999;
}

////////////////////////////////////////////////////////////////////////////////
/// Computes distance from a given point inside the shape to its boundary.

Double_t TGeoUnion::DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact,
                                   Double_t step, Double_t *safe) const
{
   if (iact<3 && safe) {
      // compute safe distance
      *safe = Safety(point,kTRUE);
      if (iact==0) return TGeoShape::Big();
      if (iact==1 && step<*safe) return TGeoShape::Big();
   }

   Double_t local[3], local1[3], master[3], ldir[3], rdir[3], pushed[3];
   memcpy(master, point, 3*sizeof(Double_t));
   Int_t i;
   TGeoBoolNode *node = (TGeoBoolNode*)this;
   Double_t d1=0., d2=0., snxt=0., eps=0.;
   fLeftMat->MasterToLocalVect(dir, ldir);
   fRightMat->MasterToLocalVect(dir, rdir);
   fLeftMat->MasterToLocal(point, local);
   Bool_t inside1 = fLeft->Contains(local);
   if (inside1) d1 = fLeft->DistFromInside(local, ldir, 3);
   else memcpy(local1, local, 3*sizeof(Double_t));
   fRightMat->MasterToLocal(point, local);
   Bool_t inside2 = fRight->Contains(local);
   if (inside2) d2 = fRight->DistFromInside(local, rdir, 3);
   if (!(inside1 | inside2)) {
   // This is a pathological case when the point is on the boundary
      d1 = fLeft->DistFromOutside(local1, ldir, 3);
      if (d1<1.E-3) {
         eps = d1+TGeoShape::Tolerance();
         for (i=0; i<3; i++) local1[i] += eps*ldir[i];
         inside1 = kTRUE;
         d1 = fLeft->DistFromInside(local1, ldir, 3);
         d1 += eps;
      } else {
         d2 = fRight->DistFromOutside(local, rdir, 3);
         if (d2<1.E-3) {
           eps = d2+TGeoShape::Tolerance();
           for (i=0; i<3; i++) local[i] += eps*rdir[i];
           inside2 = kTRUE;
           d2 = fRight->DistFromInside(local, rdir, 3);
           d2 += eps;
        }
     }
   }
   while (inside1 || inside2) {
      if (inside1 && inside2) {
         if (d1<d2) {
            snxt += d1;
            node->SetSelected(1);
            // propagate to exit of left shape
            inside1 = kFALSE;
            for (i=0; i<3; i++) master[i] += d1*dir[i];
            // check if propagated point is in right shape
            fRightMat->MasterToLocal(master, local);
            inside2 = fRight->Contains(local);
            if (!inside2) return snxt;
            d2 = fRight->DistFromInside(local, rdir, 3);
            if (d2 < TGeoShape::Tolerance()) return snxt;
         } else {
            snxt += d2;
            node->SetSelected(2);
            // propagate to exit of right shape
            inside2 = kFALSE;
            for (i=0; i<3; i++) master[i] += d2*dir[i];
            // check if propagated point is in left shape
            fLeftMat->MasterToLocal(master, local);
            inside1 = fLeft->Contains(local);
            if (!inside1) return snxt;
            d1 = fLeft->DistFromInside(local, ldir, 3);
            if (d1 < TGeoShape::Tolerance()) return snxt;
         }
      }
      if (inside1) {
         snxt += d1;
         node->SetSelected(1);
         // propagate to exit of left shape
         inside1 = kFALSE;
         for (i=0; i<3; i++) {
            master[i] += d1*dir[i];
            pushed[i] = master[i]+(1.+d1)*TGeoShape::Tolerance()*dir[i];
         }
         // check if propagated point is in right shape
         fRightMat->MasterToLocal(pushed, local);
         inside2 = fRight->Contains(local);
         if (!inside2) return snxt;
         d2 = fRight->DistFromInside(local, rdir, 3);
         if (d2 < TGeoShape::Tolerance()) return snxt;
         d2 += (1.+d1)*TGeoShape::Tolerance();
      }
      if (inside2) {
         snxt += d2;
         node->SetSelected(2);
         // propagate to exit of right shape
         inside2 = kFALSE;
         for (i=0; i<3; i++) {
            master[i] += d2*dir[i];
            pushed[i] = master[i]+(1.+d2)*TGeoShape::Tolerance()*dir[i];
         }
         // check if propagated point is in left shape
         fLeftMat->MasterToLocal(pushed, local);
         inside1 = fLeft->Contains(local);
         if (!inside1) return snxt;
         d1 = fLeft->DistFromInside(local, ldir, 3);
         if (d1 < TGeoShape::Tolerance()) return snxt;
         d1 += (1.+d2)*TGeoShape::Tolerance();
      }
   }
   return snxt;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from a given outside point to the shape.

Double_t TGeoUnion::DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact,
                                    Double_t step, Double_t *safe) const
{
   if (iact<3 && safe) {
      // compute safe distance
      *safe = Safety(point,kFALSE);
      if (iact==0) return TGeoShape::Big();
      if (iact==1 && step<*safe) return TGeoShape::Big();
   }
   TGeoBoolNode *node = (TGeoBoolNode*)this;
   Double_t local[3], ldir[3], rdir[3];
   Double_t d1, d2, snxt;
   fLeftMat->MasterToLocal(point, &local[0]);
   fLeftMat->MasterToLocalVect(dir, &ldir[0]);
   fRightMat->MasterToLocalVect(dir, &rdir[0]);
   d1 = fLeft->DistFromOutside(&local[0], &ldir[0], iact, step, safe);
   fRightMat->MasterToLocal(point, &local[0]);
   d2 = fRight->DistFromOutside(&local[0], &rdir[0], iact, step, safe);
   if (d1<d2) {
      snxt = d1;
      node->SetSelected(1);
   } else {
      snxt = d2;
      node->SetSelected(2);
   }
   return snxt;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns number of vertices for the composite shape described by this union.

Int_t TGeoUnion::GetNpoints()
{
   Int_t itot=0;
   Double_t point[3];
   Double_t tolerance = TGeoShape::Tolerance();
   if (fNpoints) return fNpoints;
   // Local points for the left shape
   Int_t nleft = fLeft->GetNmeshVertices();
   Double_t *points1 = new Double_t[3*nleft];
   fLeft->SetPoints(points1);
   // Local points for the right shape
   Int_t nright = fRight->GetNmeshVertices();
   Double_t *points2 = new Double_t[3*nright];
   fRight->SetPoints(points2);
   Double_t *points = new Double_t[3*(nleft+nright)];
   for (Int_t i=0; i<nleft; i++) {
      if (TMath::Abs(points1[3*i])<tolerance && TMath::Abs(points1[3*i+1])<tolerance) continue;
      fLeftMat->LocalToMaster(&points1[3*i], &points[3*itot]);
      fRightMat->MasterToLocal(&points[3*itot], point);
      if (!fRight->Contains(point)) itot++;
   }
   for (Int_t i=0; i<nright; i++) {
      if (TMath::Abs(points2[3*i])<tolerance && TMath::Abs(points2[3*i+1])<tolerance) continue;
      fRightMat->LocalToMaster(&points2[3*i], &points[3*itot]);
      fLeftMat->MasterToLocal(&points[3*itot], point);
      if (!fLeft->Contains(point)) itot++;
   }

   AssignPoints(itot, points);

   delete [] points1;
   delete [] points2;
   delete [] points;
   return fNpoints;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute safety distance for a union node;

Double_t TGeoUnion::Safety(const Double_t *point, Bool_t in) const
{
   Double_t local1[3], local2[3];
   fLeftMat->MasterToLocal(point,local1);
   Bool_t in1 = fLeft->Contains(local1);
   fRightMat->MasterToLocal(point,local2);
   Bool_t in2 = fRight->Contains(local2);
   Bool_t intrue = in1 | in2;
   if (intrue^in) return 0.0;
   Double_t saf1 = fLeft->Safety(local1, in1);
   Double_t saf2 = fRight->Safety(local2, in2);
   if (in1 && in2) return TMath::Min(saf1, saf2);
   if (in1)        return saf1;
   if (in2)        return saf2;
   return TMath::Min(saf1,saf2);
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoUnion::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   TGeoBoolNode::SavePrimitive(out,option);
   out << "   pBoolNode = new TGeoUnion(";
   out << fLeft->GetPointerName() << ",";
   out << fRight->GetPointerName() << ",";
   if (!fLeftMat->IsIdentity()) out << fLeftMat->GetPointerName() << ",";
   else                         out << "0,";
   if (!fRightMat->IsIdentity()) out << fRightMat->GetPointerName() << ");" << std::endl;
   else                         out << "0);" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Register 3D size of this shape.

void TGeoUnion::Sizeof3D() const
{
   TGeoBoolNode::Sizeof3D();
}


ClassImp(TGeoSubtraction);

////////////////////////////////////////////////////////////////////////////////
/// Make a clone of this. Pointers are preserved.

TGeoBoolNode *TGeoSubtraction::MakeClone() const
{
   return new TGeoSubtraction(fLeft, fRight, fLeftMat, fRightMat);
}

////////////////////////////////////////////////////////////////////////////////
/// Paint method.

void TGeoSubtraction::Paint(Option_t *option)
{
   TVirtualViewer3D *viewer = gPad->GetViewer3D();

   if (!viewer) {
      Error("Paint", "gPad->GetViewer3D() returned 0, cannot work with composite!\n");
      return;
   }

   viewer->AddCompositeOp(TBuffer3D::kCSDifference);

   TGeoBoolNode::Paint(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoSubtraction::TGeoSubtraction()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TGeoSubtraction::TGeoSubtraction(const char *expr1, const char *expr2)
          :TGeoBoolNode(expr1, expr2)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor providing pointers to components

TGeoSubtraction::TGeoSubtraction(TGeoShape *left, TGeoShape *right, TGeoMatrix *lmat, TGeoMatrix *rmat)
                :TGeoBoolNode(left,right,lmat,rmat)
{
   if (left->TestShapeBit(TGeoShape::kGeoHalfSpace)) {
      Fatal("TGeoSubstraction", "Subtractions from a half-space (%s) not allowed", left->GetName());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor
/// --- deletion of components handled by TGeoManager class.

TGeoSubtraction::~TGeoSubtraction()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Compute bounding box corresponding to a subtraction of two shapes.

void TGeoSubtraction::ComputeBBox(Double_t &dx, Double_t &dy, Double_t &dz, Double_t *origin)
{
   TGeoBBox *box = (TGeoBBox*)fLeft;
   if (box->IsNullBox()) fLeft->ComputeBBox();
   Double_t vert[24];
   Double_t pt[3];
   Int_t i;
   Double_t xmin, xmax, ymin, ymax, zmin, zmax;
   xmin = ymin = zmin = TGeoShape::Big();
   xmax = ymax = zmax = -TGeoShape::Big();
   box->SetBoxPoints(&vert[0]);
   for (i=0; i<8; i++) {
      fLeftMat->LocalToMaster(&vert[3*i], &pt[0]);
      if (pt[0]<xmin) xmin=pt[0];
      if (pt[0]>xmax) xmax=pt[0];
      if (pt[1]<ymin) ymin=pt[1];
      if (pt[1]>ymax) ymax=pt[1];
      if (pt[2]<zmin) zmin=pt[2];
      if (pt[2]>zmax) zmax=pt[2];
   }
   dx = 0.5*(xmax-xmin);
   origin[0] = 0.5*(xmin+xmax);
   dy = 0.5*(ymax-ymin);
   origin[1] = 0.5*(ymin+ymax);
   dz = 0.5*(zmax-zmin);
   origin[2] = 0.5*(zmin+zmax);
}

////////////////////////////////////////////////////////////////////////////////
/// Normal computation in POINT. The orientation is chosen so that DIR.dot.NORM>0.

void TGeoSubtraction::ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm)
{
   ThreadData_t& td = GetThreadData();
   norm[0] = norm[1] = 0.;
   norm[2] = 1.;
   Double_t local[3], ldir[3], lnorm[3];
   if (td.fSelected == 1) {
      fLeftMat->MasterToLocal(point, local);
      fLeftMat->MasterToLocalVect(dir, ldir);
      fLeft->ComputeNormal(local,ldir,lnorm);
      fLeftMat->LocalToMasterVect(lnorm, norm);
      return;
   }
   if (td.fSelected == 2) {
      fRightMat->MasterToLocal(point, local);
      fRightMat->MasterToLocalVect(dir, ldir);
      fRight->ComputeNormal(local,ldir,lnorm);
      fRightMat->LocalToMasterVect(lnorm, norm);
      return;
   }
   fRightMat->MasterToLocal(point,local);
   if (fRight->Contains(local)) {
      fRightMat->MasterToLocalVect(dir,ldir);
      fRight->ComputeNormal(local,ldir, lnorm);
      fRightMat->LocalToMasterVect(lnorm,norm);
      return;
   }
   fLeftMat->MasterToLocal(point,local);
   if (!fLeft->Contains(local)) {
      fLeftMat->MasterToLocalVect(dir,ldir);
      fLeft->ComputeNormal(local,ldir, lnorm);
      fLeftMat->LocalToMasterVect(lnorm,norm);
      return;
   }
   // point is inside left shape, but not inside the right
   local[0] = point[0]+1E-5*dir[0];
   local[1] = point[1]+1E-5*dir[1];
   local[2] = point[2]+1E-5*dir[2];
   if (Contains(local)) {
      local[0] = point[0]-1E-5*dir[0];
      local[1] = point[1]-1E-5*dir[1];
      local[2] = point[2]-1E-5*dir[2];
      if (Contains(local)) return;
   }
   ComputeNormal(local,dir,norm);
}

////////////////////////////////////////////////////////////////////////////////
/// Find if a subtraction of two shapes contains a given point

Bool_t TGeoSubtraction::Contains(const Double_t *point) const
{
   Double_t local[3];
   fLeftMat->MasterToLocal(point, &local[0]);
   Bool_t inside = fLeft->Contains(&local[0]);
   if (!inside) return kFALSE;
   fRightMat->MasterToLocal(point, &local[0]);
   inside = !fRight->Contains(&local[0]);
   return inside;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute minimum distance to shape vertices

Int_t TGeoSubtraction::DistanceToPrimitive(Int_t /*px*/, Int_t /*py*/)
{
   return 9999;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from a given point inside to the shape boundary.

Double_t TGeoSubtraction::DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact,
                              Double_t step, Double_t *safe) const
{
   if (iact<3 && safe) {
      // compute safe distance
      *safe = Safety(point,kTRUE);
      if (iact==0) return TGeoShape::Big();
      if (iact==1 && step<*safe) return TGeoShape::Big();
   }
   TGeoBoolNode *node = (TGeoBoolNode*)this;
   Double_t local[3], ldir[3], rdir[3];
   Double_t d1, d2, snxt=0.;
   fLeftMat->MasterToLocal(point, &local[0]);
   fLeftMat->MasterToLocalVect(dir, &ldir[0]);
   fRightMat->MasterToLocalVect(dir, &rdir[0]);
   d1 = fLeft->DistFromInside(&local[0], &ldir[0], iact, step, safe);
   fRightMat->MasterToLocal(point, &local[0]);
   d2 = fRight->DistFromOutside(&local[0], &rdir[0], iact, step, safe);
   if (d1<d2) {
      snxt = d1;
      node->SetSelected(1);
   } else {
      snxt = d2;
      node->SetSelected(2);
   }
   return snxt;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from a given point outside to the shape.

Double_t TGeoSubtraction::DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact,
                              Double_t step, Double_t *safe) const
{
   if (iact<3 && safe) {
      // compute safe distance
      *safe = Safety(point,kFALSE);
      if (iact==0) return TGeoShape::Big();
      if (iact==1 && step<*safe) return TGeoShape::Big();
   }
   TGeoBoolNode *node = (TGeoBoolNode*)this;
   Double_t local[3], master[3], ldir[3], rdir[3];
   memcpy(&master[0], point, 3*sizeof(Double_t));
   Int_t i;
   Double_t d1, d2, snxt=0.;
   fRightMat->MasterToLocal(point, &local[0]);
   fLeftMat->MasterToLocalVect(dir, &ldir[0]);
   fRightMat->MasterToLocalVect(dir, &rdir[0]);
   // check if inside '-'
   Bool_t inside = fRight->Contains(&local[0]);
   Double_t epsil = 0.;
   while (1) {
      if (inside) {
         // propagate to outside of '-'
         node->SetSelected(2);
         d1 = fRight->DistFromInside(&local[0], &rdir[0], iact, step, safe);
         snxt += d1+epsil;
         for (i=0; i<3; i++) master[i] += (d1+1E-8)*dir[i];
         epsil = 1.E-8;
         // now master outside '-'; check if inside '+'
         fLeftMat->MasterToLocal(&master[0], &local[0]);
         if (fLeft->Contains(&local[0])) return snxt;
      }
      // master outside '-' and outside '+' ;  find distances to both
      node->SetSelected(1);
      fLeftMat->MasterToLocal(&master[0], &local[0]);
      d2 = fLeft->DistFromOutside(&local[0], &ldir[0], iact, step, safe);
      if (d2>1E20) return TGeoShape::Big();

      fRightMat->MasterToLocal(&master[0], &local[0]);
      d1 = fRight->DistFromOutside(&local[0], &rdir[0], iact, step, safe);
      if (d2<d1-TGeoShape::Tolerance()) {
         snxt += d2+epsil;
         return snxt;
      }
      // propagate to '-'
      snxt += d1+epsil;
      for (i=0; i<3; i++) master[i] += (d1+1E-8)*dir[i];
      epsil = 1.E-8;
      // now inside '-' and not inside '+'
      fRightMat->MasterToLocal(&master[0], &local[0]);
      inside = kTRUE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns number of vertices for the composite shape described by this subtraction.

Int_t TGeoSubtraction::GetNpoints()
{
   Int_t itot=0;
   Double_t point[3];
   Double_t tolerance = TGeoShape::Tolerance();
   if (fNpoints) return fNpoints;
   Int_t nleft = fLeft->GetNmeshVertices();
   Int_t nright = fRight->GetNmeshVertices();
   Double_t *points = new Double_t[3*(nleft+nright)];
   Double_t *points1 = new Double_t[3*nleft];
   fLeft->SetPoints(points1);
   for (Int_t i=0; i<nleft; i++) {
      if (TMath::Abs(points1[3*i])<tolerance && TMath::Abs(points1[3*i+1])<tolerance) continue;
      fLeftMat->LocalToMaster(&points1[3*i], &points[3*itot]);
      fRightMat->MasterToLocal(&points[3*itot], point);
      if (!fRight->Contains(point)) itot++;
   }
   Double_t *points2 = new Double_t[3*nright];
   fRight->SetPoints(points2);
   for (Int_t i=0; i<nright; i++) {
      if (TMath::Abs(points2[3*i])<tolerance && TMath::Abs(points2[3*i+1])<tolerance) continue;
      fRightMat->LocalToMaster(&points2[3*i], &points[3*itot]);
      fLeftMat->MasterToLocal(&points[3*itot], point);
      if (fLeft->Contains(point)) itot++;
   }

   AssignPoints(itot, points);

   delete [] points1;
   delete [] points2;
   delete [] points;
   return fNpoints;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute safety distance for a union node;

Double_t TGeoSubtraction::Safety(const Double_t *point, Bool_t in) const
{
   Double_t local1[3], local2[3];
   fLeftMat->MasterToLocal(point,local1);
   Bool_t in1 = fLeft->Contains(local1);
   fRightMat->MasterToLocal(point,local2);
   Bool_t in2 = fRight->Contains(local2);
   Bool_t intrue = in1 && (!in2);
   if (in^intrue) return 0.0;
   Double_t saf1 = fLeft->Safety(local1, in1);
   Double_t saf2 = fRight->Safety(local2, in2);
   if (in1 && in2) return saf2;
   if (in1)        return TMath::Min(saf1,saf2);
   if (in2)        return TMath::Max(saf1,saf2);
   return saf1;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoSubtraction::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   TGeoBoolNode::SavePrimitive(out,option);
   out << "   pBoolNode = new TGeoSubtraction(";
   out << fLeft->GetPointerName() << ",";
   out << fRight->GetPointerName() << ",";
   if (!fLeftMat->IsIdentity()) out << fLeftMat->GetPointerName() << ",";
   else                         out << "0,";
   if (!fRightMat->IsIdentity()) out << fRightMat->GetPointerName() << ");" << std::endl;
   else                         out << "0);" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Register 3D size of this shape.

void TGeoSubtraction::Sizeof3D() const
{
   TGeoBoolNode::Sizeof3D();
}

ClassImp(TGeoIntersection);

////////////////////////////////////////////////////////////////////////////////
/// Make a clone of this. Pointers are preserved.

TGeoBoolNode *TGeoIntersection::MakeClone() const
{
   return new TGeoIntersection(fLeft, fRight, fLeftMat, fRightMat);
}

////////////////////////////////////////////////////////////////////////////////
/// Paint method.

void TGeoIntersection::Paint(Option_t *option)
{
   TVirtualViewer3D *viewer = gPad->GetViewer3D();

   if (!viewer) {
      Error("Paint", "gPad->GetViewer3D() returned 0, cannot work with composite!\n");
      return;
   }

   viewer->AddCompositeOp(TBuffer3D::kCSIntersection);

   TGeoBoolNode::Paint(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoIntersection::TGeoIntersection()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TGeoIntersection::TGeoIntersection(const char *expr1, const char *expr2)
          :TGeoBoolNode(expr1, expr2)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor providing pointers to components

TGeoIntersection::TGeoIntersection(TGeoShape *left, TGeoShape *right, TGeoMatrix *lmat, TGeoMatrix *rmat)
                 :TGeoBoolNode(left,right,lmat,rmat)
{
   Bool_t hs1 = (fLeft->TestShapeBit(TGeoShape::kGeoHalfSpace))?kTRUE:kFALSE;
   Bool_t hs2 = (fRight->TestShapeBit(TGeoShape::kGeoHalfSpace))?kTRUE:kFALSE;
   if (hs1 && hs2) Fatal("ctor", "cannot intersect two half-spaces: %s * %s", left->GetName(), right->GetName());
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor
/// --- deletion of components handled by TGeoManager class.

TGeoIntersection::~TGeoIntersection()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Compute bounding box corresponding to a intersection of two shapes.

void TGeoIntersection::ComputeBBox(Double_t &dx, Double_t &dy, Double_t &dz, Double_t *origin)
{
   Bool_t hs1 = (fLeft->TestShapeBit(TGeoShape::kGeoHalfSpace))?kTRUE:kFALSE;
   Bool_t hs2 = (fRight->TestShapeBit(TGeoShape::kGeoHalfSpace))?kTRUE:kFALSE;
   Double_t vert[48];
   Double_t pt[3];
   Int_t i;
   Double_t xmin1, xmax1, ymin1, ymax1, zmin1, zmax1;
   Double_t xmin2, xmax2, ymin2, ymax2, zmin2, zmax2;
   xmin1 = ymin1 = zmin1 = xmin2 = ymin2 = zmin2 = TGeoShape::Big();
   xmax1 = ymax1 = zmax1 = xmax2 = ymax2 = zmax2 =  -TGeoShape::Big();
   if (!hs1) {
      if (((TGeoBBox*)fLeft)->IsNullBox()) fLeft->ComputeBBox();
      ((TGeoBBox*)fLeft)->SetBoxPoints(&vert[0]);
      for (i=0; i<8; i++) {
         fLeftMat->LocalToMaster(&vert[3*i], &pt[0]);
         if (pt[0]<xmin1) xmin1=pt[0];
         if (pt[0]>xmax1) xmax1=pt[0];
         if (pt[1]<ymin1) ymin1=pt[1];
         if (pt[1]>ymax1) ymax1=pt[1];
         if (pt[2]<zmin1) zmin1=pt[2];
         if (pt[2]>zmax1) zmax1=pt[2];
      }
   }
   if (!hs2) {
      if (((TGeoBBox*)fRight)->IsNullBox()) fRight->ComputeBBox();
      ((TGeoBBox*)fRight)->SetBoxPoints(&vert[24]);
      for (i=8; i<16; i++) {
         fRightMat->LocalToMaster(&vert[3*i], &pt[0]);
         if (pt[0]<xmin2) xmin2=pt[0];
         if (pt[0]>xmax2) xmax2=pt[0];
         if (pt[1]<ymin2) ymin2=pt[1];
         if (pt[1]>ymax2) ymax2=pt[1];
         if (pt[2]<zmin2) zmin2=pt[2];
         if (pt[2]>zmax2) zmax2=pt[2];
      }
   }
   if (hs1) {
      dx = 0.5*(xmax2-xmin2);
      origin[0] = 0.5*(xmax2+xmin2);
      dy = 0.5*(ymax2-ymin2);
      origin[1] = 0.5*(ymax2+ymin2);
      dz = 0.5*(zmax2-zmin2);
      origin[2] = 0.5*(zmax2+zmin2);
      return;
   }
   if (hs2) {
      dx = 0.5*(xmax1-xmin1);
      origin[0] = 0.5*(xmax1+xmin1);
      dy = 0.5*(ymax1-ymin1);
      origin[1] = 0.5*(ymax1+ymin1);
      dz = 0.5*(zmax1-zmin1);
      origin[2] = 0.5*(zmax1+zmin1);
      return;
   }
   Double_t sort[4];
   Int_t isort[4];
   sort[0] = xmin1;
   sort[1] = xmax1;
   sort[2] = xmin2;
   sort[3] = xmax2;
   TMath::Sort(4, &sort[0], &isort[0], kFALSE);
   if (isort[1]%2) {
      Warning("ComputeBBox", "shapes %s and %s do not intersect", fLeft->GetName(), fRight->GetName());
      dx = dy = dz = 0;
      memset(origin, 0, 3*sizeof(Double_t));
      return;
   }
   dx = 0.5*(sort[isort[2]]-sort[isort[1]]);
   origin[0] = 0.5*(sort[isort[1]]+sort[isort[2]]);
   sort[0] = ymin1;
   sort[1] = ymax1;
   sort[2] = ymin2;
   sort[3] = ymax2;
   TMath::Sort(4, &sort[0], &isort[0], kFALSE);
   if (isort[1]%2) {
      Warning("ComputeBBox", "shapes %s and %s do not intersect", fLeft->GetName(), fRight->GetName());
      dx = dy = dz = 0;
      memset(origin, 0, 3*sizeof(Double_t));
      return;
   }
   dy = 0.5*(sort[isort[2]]-sort[isort[1]]);
   origin[1] = 0.5*(sort[isort[1]]+sort[isort[2]]);
   sort[0] = zmin1;
   sort[1] = zmax1;
   sort[2] = zmin2;
   sort[3] = zmax2;
   TMath::Sort(4, &sort[0], &isort[0], kFALSE);
   if (isort[1]%2) {
      Warning("ComputeBBox", "shapes %s and %s do not intersect", fLeft->GetName(), fRight->GetName());
      dx = dy = dz = 0;
      memset(origin, 0, 3*sizeof(Double_t));
      return;
   }
   dz = 0.5*(sort[isort[2]]-sort[isort[1]]);
   origin[2] = 0.5*(sort[isort[1]]+sort[isort[2]]);
}

////////////////////////////////////////////////////////////////////////////////
/// Normal computation in POINT. The orientation is chosen so that DIR.dot.NORM>0.

void TGeoIntersection::ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm)
{
   ThreadData_t& td = GetThreadData();
   Double_t local[3], ldir[3], lnorm[3];
   norm[0] = norm[1] = 0.;
   norm[2] = 1.;
   if (td.fSelected == 1) {
      fLeftMat->MasterToLocal(point, local);
      fLeftMat->MasterToLocalVect(dir, ldir);
      fLeft->ComputeNormal(local,ldir,lnorm);
      fLeftMat->LocalToMasterVect(lnorm, norm);
      return;
   }
   if (td.fSelected == 2) {
      fRightMat->MasterToLocal(point, local);
      fRightMat->MasterToLocalVect(dir, ldir);
      fRight->ComputeNormal(local,ldir,lnorm);
      fRightMat->LocalToMasterVect(lnorm, norm);
      return;
   }
   fLeftMat->MasterToLocal(point,local);
   if (!fLeft->Contains(local)) {
      fLeftMat->MasterToLocalVect(dir,ldir);
      fLeft->ComputeNormal(local,ldir,lnorm);
      fLeftMat->LocalToMasterVect(lnorm,norm);
      return;
   }
   fRightMat->MasterToLocal(point,local);
   if (!fRight->Contains(local)) {
      fRightMat->MasterToLocalVect(dir,ldir);
      fRight->ComputeNormal(local,ldir,lnorm);
      fRightMat->LocalToMasterVect(lnorm,norm);
      return;
   }
   // point is inside intersection.
   local[0] = point[0] + 1E-5*dir[0];
   local[1] = point[1] + 1E-5*dir[1];
   local[2] = point[2] + 1E-5*dir[2];
   if (Contains(local)) {
      local[0] = point[0] - 1E-5*dir[0];
      local[1] = point[1] - 1E-5*dir[1];
      local[2] = point[2] - 1E-5*dir[2];
      if (Contains(local)) return;
   }
   ComputeNormal(local,dir,norm);
}

////////////////////////////////////////////////////////////////////////////////
/// Find if a intersection of two shapes contains a given point

Bool_t TGeoIntersection::Contains(const Double_t *point) const
{
   Double_t local[3];
   fLeftMat->MasterToLocal(point, &local[0]);
   Bool_t inside = fLeft->Contains(&local[0]);
   if (!inside) return kFALSE;
   fRightMat->MasterToLocal(point, &local[0]);
   inside = fRight->Contains(&local[0]);
   return inside;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute minimum distance to shape vertices

Int_t TGeoIntersection::DistanceToPrimitive(Int_t /*px*/, Int_t /*py*/)
{
   return 9999;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from a given point inside to the shape boundary.

Double_t TGeoIntersection::DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact,
                              Double_t step, Double_t *safe) const
{
   if (iact<3 && safe) {
      // compute safe distance
      *safe = Safety(point,kTRUE);
      if (iact==0) return TGeoShape::Big();
      if (iact==1 && step<*safe) return TGeoShape::Big();
   }
   TGeoBoolNode *node = (TGeoBoolNode*)this;
   Double_t local[3], ldir[3], rdir[3];
   Double_t d1, d2, snxt=0.;
   fLeftMat->MasterToLocal(point, &local[0]);
   fLeftMat->MasterToLocalVect(dir, &ldir[0]);
   fRightMat->MasterToLocalVect(dir, &rdir[0]);
   d1 = fLeft->DistFromInside(&local[0], &ldir[0], iact, step, safe);
   fRightMat->MasterToLocal(point, &local[0]);
   d2 = fRight->DistFromInside(&local[0], &rdir[0], iact, step, safe);
   if (d1<d2) {
      snxt = d1;
      node->SetSelected(1);
   } else {
      snxt = d2;
      node->SetSelected(2);
   }
   return snxt;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from a given point outside to the shape.

Double_t TGeoIntersection::DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact,
                              Double_t step, Double_t *safe) const
{
   Double_t tol = TGeoShape::Tolerance();
   if (iact<3 && safe) {
      // compute safe distance
      *safe = Safety(point,kFALSE);
      if (iact==0) return TGeoShape::Big();
      if (iact==1 && step<*safe) return TGeoShape::Big();
   }
   TGeoBoolNode *node = (TGeoBoolNode*)this;
   Double_t lpt[3], rpt[3], master[3], ldir[3], rdir[3];
   memcpy(master, point, 3*sizeof(Double_t));
   Int_t i;
   Double_t d1 = 0.;
   Double_t d2 = 0.;
   fLeftMat->MasterToLocal(point, lpt);
   fRightMat->MasterToLocal(point, rpt);
   fLeftMat->MasterToLocalVect(dir, ldir);
   fRightMat->MasterToLocalVect(dir, rdir);
   Bool_t inleft = fLeft->Contains(lpt);
   Bool_t inright = fRight->Contains(rpt);
   node->SetSelected(0);
   Double_t snext = 0.0;
   if (inleft && inright) {
      // It is vey likely to have a numerical issue and the point should
      // be logically outside one of the shapes
      d1 = fLeft->DistFromInside(lpt,ldir,3);
      d2 = fRight->DistFromInside(rpt,rdir,3);
      if (d1<1.E-3) inleft = kFALSE;
      if (d2<1.E-3) inright = kFALSE;
      if (inleft && inright) return snext;
   }

   while (1) {
      d1 = d2 = 0;
      if (!inleft)  {
         d1 = fLeft->DistFromOutside(lpt,ldir,3);
         d1 = TMath::Max(d1,tol);
         if (d1>1E20) return TGeoShape::Big();
      }
      if (!inright) {
         d2 = fRight->DistFromOutside(rpt,rdir,3);
         d2 = TMath::Max(d2,tol);
         if (d2>1E20) return TGeoShape::Big();
      }

      if (d1>d2) {
         // propagate to left shape
         snext += d1;
         node->SetSelected(1);
         inleft = kTRUE;
         for (i=0; i<3; i++) master[i] += d1*dir[i];
         fRightMat->MasterToLocal(master,rpt);
         // Push rpt to avoid a bad boundary condition
         for (i=0; i<3; i++) rpt[i] += tol*rdir[i];
         // check if propagated point is inside right shape
         inright = fRight->Contains(rpt);
         if (inright) return snext;
         // here inleft=true, inright=false
      } else {
         // propagate to right shape
         snext += d2;
         node->SetSelected(2);
         inright = kTRUE;
         for (i=0; i<3; i++) master[i] += d2*dir[i];
         fLeftMat->MasterToLocal(master,lpt);
         // Push lpt to avoid a bad boundary condition
         for (i=0; i<3; i++) lpt[i] += tol*ldir[i];
         // check if propagated point is inside left shape
         inleft = fLeft->Contains(lpt);
         if (inleft) return snext;
         // here inleft=false, inright=true
      }
   }
   return snext;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns number of vertices for the composite shape described by this intersection.

Int_t TGeoIntersection::GetNpoints()
{
   Int_t itot=0;
   Double_t point[3];
   Double_t tolerance = TGeoShape::Tolerance();
   if (fNpoints) return fNpoints;
   Int_t nleft = fLeft->GetNmeshVertices();
   Int_t nright = fRight->GetNmeshVertices();
   Double_t *points = new Double_t[3*(nleft+nright)];
   Double_t *points1 = new Double_t[3*nleft];
   fLeft->SetPoints(points1);
   for (Int_t i=0; i<nleft; i++) {
      if (TMath::Abs(points1[3*i])<tolerance && TMath::Abs(points1[3*i+1])<tolerance) continue;
      fLeftMat->LocalToMaster(&points1[3*i], &points[3*itot]);
      fRightMat->MasterToLocal(&points[3*itot], point);
      if (fRight->Contains(point)) itot++;
   }
   Double_t *points2 = new Double_t[3*nright];
   fRight->SetPoints(points2);
   for (Int_t i=0; i<nright; i++) {
      if (TMath::Abs(points2[3*i])<tolerance && TMath::Abs(points2[3*i+1])<tolerance) continue;
      fRightMat->LocalToMaster(&points2[3*i], &points[3*itot]);
      fLeftMat->MasterToLocal(&points[3*itot], point);
      if (fLeft->Contains(point)) itot++;
   }

   AssignPoints(itot, points);

   delete [] points1;
   delete [] points2;
   delete [] points;
   return fNpoints;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute safety distance for a union node;

Double_t TGeoIntersection::Safety(const Double_t *point, Bool_t in) const
{
   Double_t local1[3], local2[3];
   fLeftMat->MasterToLocal(point,local1);
   Bool_t in1 = fLeft->Contains(local1);
   fRightMat->MasterToLocal(point,local2);
   Bool_t in2 = fRight->Contains(local2);
   Bool_t intrue = in1 & in2;
   if (in^intrue) return 0.0;
   Double_t saf1 = fLeft->Safety(local1, in1);
   Double_t saf2 = fRight->Safety(local2, in2);
   if (in1 && in2) return TMath::Min(saf1, saf2);
   if (in1)        return saf2;
   if (in2)        return saf1;
   return TMath::Max(saf1,saf2);
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoIntersection::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   TGeoBoolNode::SavePrimitive(out,option);
   out << "   pBoolNode = new TGeoIntersection(";
   out << fLeft->GetPointerName() << ",";
   out << fRight->GetPointerName() << ",";
   if (!fLeftMat->IsIdentity()) out << fLeftMat->GetPointerName() << ",";
   else                         out << "0,";
   if (!fRightMat->IsIdentity()) out << fRightMat->GetPointerName() << ");" << std::endl;
   else                         out << "0);" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Register 3D size of this shape.

void TGeoIntersection::Sizeof3D() const
{
   TGeoBoolNode::Sizeof3D();
}
