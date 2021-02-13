// @(#)root/eve7:$Id$
// Author: Matevz Tadel 2007, 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "Rtypes.h"
#include <cassert>


#include <ROOT/REveGeoPolyShape.hxx>
#include <ROOT/REveGeoShape.hxx>
#include <ROOT/REveUtil.hxx>
#include <ROOT/REveGluTess.hxx>
#include <ROOT/REveRenderData.hxx>

#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "CsgOps.h"

#include "TGeoBoolNode.h"
#include "TGeoCompositeShape.h"
#include "TGeoMatrix.h"

using namespace ROOT::Experimental;

/** \class REveGeoPolyShape
\ingroup REve
Description of REveGeoPolyShape
*/

Bool_t REveGeoPolyShape::fgAutoEnforceTriangles = kTRUE;
Bool_t REveGeoPolyShape::fgAutoCalculateNormals = kFALSE;

void   REveGeoPolyShape::SetAutoEnforceTriangles(Bool_t f) { fgAutoEnforceTriangles = f; }
Bool_t REveGeoPolyShape::GetAutoEnforceTriangles()         { return fgAutoEnforceTriangles; }
void   REveGeoPolyShape::SetAutoCalculateNormals(Bool_t f) { fgAutoCalculateNormals = f; }
Bool_t REveGeoPolyShape::GetAutoCalculateNormals()         { return fgAutoCalculateNormals; }

////////////////////////////////////////////////////////////////////////
/// Function produces mesh for provided shape, applying matrix to the result

std::unique_ptr<RootCsg::TBaseMesh> MakeGeoMesh(TGeoMatrix *matr, TGeoShape *shape)
{
   TGeoCompositeShape *comp = dynamic_cast<TGeoCompositeShape *> (shape);

   std::unique_ptr<RootCsg::TBaseMesh> res;

   if (!comp) {
      std::unique_ptr<TBuffer3D> b3d(shape->MakeBuffer3D());

      if (matr) {
         Double_t *v = b3d->fPnts;
         Double_t buf[3];
         for (UInt_t i = 0; i < b3d->NbPnts(); ++i) {
            buf[0] = v[i*3];
            buf[1] = v[i*3+1];
            buf[2] = v[i*3+2];
            matr->LocalToMaster(buf, &v[i*3]);
         }
      }

      res.reset(RootCsg::ConvertToMesh(*b3d.get()));
   } else {
      auto node = comp->GetBoolNode();

      TGeoHMatrix mleft, mright;
      if (matr) { mleft = *matr; mright = *matr; }

      mleft.Multiply(node->GetLeftMatrix());
      auto left = MakeGeoMesh(&mleft, node->GetLeftShape());

      mright.Multiply(node->GetRightMatrix());
      auto right = MakeGeoMesh(&mright, node->GetRightShape());

      if (node->IsA() == TGeoUnion::Class()) res.reset(RootCsg::BuildUnion(left.get(), right.get()));
      if (node->IsA() == TGeoIntersection::Class()) res.reset(RootCsg::BuildIntersection(left.get(), right.get()));
      if (node->IsA() == TGeoSubtraction::Class()) res.reset(RootCsg::BuildDifference(left.get(), right.get()));
   }

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Produce all polygons from composite shape

void REveGeoPolyShape::BuildFromComposite(TGeoCompositeShape *cshape, Int_t n_seg)
{
   fOrigin[0] = cshape->GetOrigin()[0];
   fOrigin[1] = cshape->GetOrigin()[1];
   fOrigin[2] = cshape->GetOrigin()[2];
   fDX = cshape->GetDX();
   fDY = cshape->GetDY();
   fDZ = cshape->GetDZ();

   REveGeoManagerHolder gmgr(REveGeoShape::GetGeoManager(), n_seg);

   auto mesh = MakeGeoMesh(nullptr, cshape);

   Int_t nv = mesh->NumberOfVertices();
   fVertices.reserve(3 * nv);

   for (Int_t i = 0; i < nv; ++i) {
      auto v = mesh->GetVertex(i);
      fVertices.insert(fVertices.end(), v, v + 3);
   }

   fNbPols = mesh->NumberOfPolys();

   Int_t descSize = 0;

   for (Int_t i = 0; i < fNbPols; ++i) descSize += mesh->SizeOfPoly(i) + 1;

   fPolyDesc.reserve(descSize);

   for (Int_t polyIndex = 0; polyIndex < fNbPols; ++polyIndex) {
      Int_t polySize = mesh->SizeOfPoly(polyIndex);

      fPolyDesc.push_back(polySize);

      for (Int_t i = 0; i < polySize; ++i)
         fPolyDesc.push_back(mesh->GetVertexIndex(polyIndex, i));
   }

   if (fgAutoEnforceTriangles) EnforceTriangles();
   if (fgAutoCalculateNormals) CalculateNormals();
}

////////////////////////////////////////////////////////////////////////////////
/// Produce all polygons from normal shape

void REveGeoPolyShape::BuildFromShape(TGeoShape *shape, Int_t n_seg)
{
   TGeoBBox *box = dynamic_cast<TGeoBBox *> (shape);

   if (box) {
      fOrigin[0] = box->GetOrigin()[0];
      fOrigin[1] = box->GetOrigin()[1];
      fOrigin[2] = box->GetOrigin()[2];
      fDX = box->GetDX();
      fDY = box->GetDY();
      fDZ = box->GetDZ();
   }

   REveGeoManagerHolder gmgr(REveGeoShape::GetGeoManager(), n_seg);

   std::unique_ptr<TBuffer3D> b3d(shape->MakeBuffer3D());

   SetFromBuff3D(*b3d.get());
}

////////////////////////////////////////////////////////////////////////////////

void REveGeoPolyShape::FillRenderData(REveRenderData &rd)
{
   // We know all elements are triangles. Or at least they should be.

   rd.Reserve(fVertices.size(), fNormals.size(), 2 + fNbPols * 3);

   for (auto &v: fVertices)
      rd.PushV(v);

   for (auto &n: fNormals)
      rd.PushN(n);

   rd.PushI(REveRenderData::GL_TRIANGLES);
   rd.PushI(fNbPols);

   // count number of index entries etc
   for (Int_t i = 0, j = 0; i < fNbPols; ++i) {
      assert(fPolyDesc[j] == 3);

      rd.PushI(fPolyDesc[j + 1], fPolyDesc[j + 2], fPolyDesc[j + 3]);
      j += 1 + fPolyDesc[j];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set data-members from a Csg mesh.


void REveGeoPolyShape::SetFromBuff3D(const TBuffer3D& buffer)
{
   fNbPols = (Int_t) buffer.NbPols();

   if (fNbPols == 0) return;

   fVertices.insert(fVertices.end(), buffer.fPnts, buffer.fPnts + 3 * buffer.NbPnts());

   Int_t *segs = buffer.fSegs;
   Int_t *pols = buffer.fPols;

   Int_t descSize = 0;

   for (Int_t i = 0, j = 1; i < fNbPols; ++i, ++j)
   {
      descSize += pols[j] + 1;
      j += pols[j] + 1;
   }

   fPolyDesc.resize(descSize);

   for (Int_t numPol = 0, currInd = 0, j = 1; numPol < fNbPols; ++numPol)
   {
      Int_t segmentInd = pols[j] + j;
      Int_t segmentCol = pols[j];
      Int_t s1 = pols[segmentInd];
      segmentInd--;
      Int_t s2 = pols[segmentInd];
      segmentInd--;
      Int_t segEnds[] = {segs[s1 * 3 + 1], segs[s1 * 3 + 2],
                         segs[s2 * 3 + 1], segs[s2 * 3 + 2]};
      Int_t numPnts[3];

      if (segEnds[0] == segEnds[2]) {
         numPnts[0] = segEnds[1]; numPnts[1] = segEnds[0]; numPnts[2] = segEnds[3];
      } else if (segEnds[0] == segEnds[3]) {
         numPnts[0] = segEnds[1]; numPnts[1] = segEnds[0]; numPnts[2] = segEnds[2];
      } else if (segEnds[1] == segEnds[2]) {
         numPnts[0] = segEnds[0]; numPnts[1] = segEnds[1]; numPnts[2] = segEnds[3];
      } else {
         numPnts[0] = segEnds[0]; numPnts[1] = segEnds[1]; numPnts[2] = segEnds[2];
      }

      fPolyDesc[currInd] = 3;
      Int_t sizeInd = currInd++;
      fPolyDesc[currInd++] = numPnts[0];
      fPolyDesc[currInd++] = numPnts[1];
      fPolyDesc[currInd++] = numPnts[2];
      Int_t lastAdded = numPnts[2];

      Int_t end = j + 1;
      for (; segmentInd != end; segmentInd--) {
         segEnds[0] = segs[pols[segmentInd] * 3 + 1];
         segEnds[1] = segs[pols[segmentInd] * 3 + 2];
         if (segEnds[0] == lastAdded) {
            fPolyDesc[currInd++] = segEnds[1];
            lastAdded = segEnds[1];
         } else {
            fPolyDesc[currInd++] = segEnds[0];
            lastAdded = segEnds[0];
         }
         ++fPolyDesc[sizeInd];
      }
      j += segmentCol + 2;
   }

   if (fgAutoEnforceTriangles) EnforceTriangles();
   if (fgAutoCalculateNormals) CalculateNormals();
}

////////////////////////////////////////////////////////////////////////////////
/// Use GLU tesselator to replace all polygons with N > 3 with triangles.
/// After this call polygon descriptions are changed.
/// New vertices are not expected -- exception is thrown if this is
/// requested by the triangulator. Support for adding of new vertices can be
/// provided.

void REveGeoPolyShape::EnforceTriangles()
{
   EveGlu::TriangleCollector tc;

   tc.ProcessData(fVertices, fPolyDesc, fNbPols);

   fPolyDesc.swap(tc.RefPolyDesc());
   fNbPols = tc.GetNTrianlges();
}

////////////////////////////////////////////////////////////////////////////////
/// CalculateNormals per polygon (flat shading)

void REveGeoPolyShape::CalculateNormals()
{
   fNormals.resize(3 * fNbPols);
   if (fNbPols == 0) return;
   Double_t *pnts = &fVertices[0];
   for (Int_t i = 0, j = 0; i < fNbPols; ++i)
   {
      Int_t polEnd = fPolyDesc[j] + j + 1;
      UInt_t norm[] = {fPolyDesc[j + 1], fPolyDesc[j + 2], fPolyDesc[j + 3]};
      j += 4;
      Int_t check = CheckPoints(norm, norm);
      Int_t ngood = check;
      if (check == 3) {
         TMath::Normal2Plane(pnts + norm[0] * 3, pnts + norm[1] * 3,
                             pnts + norm[2] * 3, &fNormals[i * 3]);
         j = polEnd;
         continue;
      }
      while (j < polEnd)
      {
         norm[ngood++] = fPolyDesc[j++];
         if (ngood == 3) {
            ngood = CheckPoints(norm, norm);
            if (ngood == 3) {
               TMath::Normal2Plane(pnts + norm[0] * 3, pnts + norm[1] * 3,
                                   pnts + norm[2] * 3, &fNormals[i * 3]);
               j = polEnd;
               break;
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// CheckPoints

Int_t REveGeoPolyShape::CheckPoints(const UInt_t *source, UInt_t *dest) const
{
   const Double_t * p1 = &fVertices[source[0] * 3];
   const Double_t * p2 = &fVertices[source[1] * 3];
   const Double_t * p3 = &fVertices[source[2] * 3];
   Int_t retVal = 1;

   if (Eq(p1, p2)) {
      dest[0] = source[0];
      if (!Eq(p1, p3) ) {
         dest[1] = source[2];
         retVal = 2;
      }
   } else if (Eq(p1, p3)) {
      dest[0] = source[0];
      dest[1] = source[1];
      retVal = 2;
   } else {
      dest[0] = source[0];
      dest[1] = source[1];
      retVal = 2;
      if (!Eq(p2, p3)) {
         dest[2] = source[2];
         retVal = 3;
      }
   }

   return retVal;
}

////////////////////////////////////////////////////////////////////////////////
/// Test equality of points with epsilon 1e-10.

Bool_t REveGeoPolyShape::Eq(const Double_t *p1, const Double_t *p2)
{
   Double_t dx = TMath::Abs(p1[0] - p2[0]);
   Double_t dy = TMath::Abs(p1[1] - p2[1]);
   Double_t dz = TMath::Abs(p1[2] - p2[2]);
   return (dx < 1e-10) && (dy < 1e-10) && (dz < 1e-10);
}

////////////////////////////////////////////////////////////////////////////////
/// Fill the passed buffer 3D.

void REveGeoPolyShape::FillBuffer3D(TBuffer3D& b, Int_t reqSections, Bool_t) const
{
   if (reqSections & TBuffer3D::kCore)
   {
      // If writing core section all others will be invalid
      b.ClearSectionsValid();

      b.fID = const_cast<REveGeoPolyShape*>(this);
      b.fColor = kMagenta;
      b.fTransparency = 0;
      b.fLocalFrame = kFALSE;
      b.fReflection = kTRUE;

      b.SetSectionsValid(TBuffer3D::kCore);
   }

   if ((reqSections & TBuffer3D::kRawSizes) || (reqSections & TBuffer3D::kRaw))
   {
      Int_t nvrt = fVertices.size() / 3;
      Int_t nseg = 0;

      std::map<Edge_t, Int_t> edges;

      const UInt_t *pd = &fPolyDesc[0];
      for (Int_t i = 0; i < fNbPols; ++i) {
         Int_t nv = pd[0];
         ++pd;
         for (Int_t j = 0; j < nv; ++j) {
            Edge_t e(pd[j], (j != nv - 1) ? pd[j + 1] : pd[0]);
            if (edges.find(e) == edges.end()) {
               edges.insert(std::make_pair(e, 0));
               ++nseg;
            }
         }
         pd += nv;
      }

      b.SetRawSizes(nvrt, 3*nvrt, nseg, 3*nseg, fNbPols, fNbPols+fPolyDesc.size());

      memcpy(b.fPnts, &fVertices[0], sizeof(Double_t)*fVertices.size());

      Int_t si = 0, scnt = 0;
      for (auto &edge : edges) {
         b.fSegs[si++] = 0;
         b.fSegs[si++] = edge.first.fI;
         b.fSegs[si++] = edge.first.fJ;
         edge.second = scnt++;
      }

      Int_t pi = 0;
      pd = &fPolyDesc[0];
      for (Int_t i = 0; i < fNbPols; ++i) {
         Int_t nv = pd[0];
         ++pd;
         b.fPols[pi++] = 0;
         b.fPols[pi++] = nv;
         for (Int_t j = 0; j < nv; ++j) {
            b.fPols[pi++] = edges[Edge_t(pd[j], (j != nv - 1) ? pd[j + 1] : pd[0])];
         }
         pd += nv;
      }

      b.SetSectionsValid(TBuffer3D::kRawSizes | TBuffer3D::kRaw);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill static buffer 3D.

const TBuffer3D& REveGeoPolyShape::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
{
   static TBuffer3D buf(TBuffer3DTypes::kGeneric);

   FillBuffer3D(buf, reqSections, localFrame);

   return buf;
}

////////////////////////////////////////////////////////////////////////////////
/// Create buffer 3D and fill it with point/segment/poly data.

TBuffer3D *REveGeoPolyShape::MakeBuffer3D() const
{
   auto *buf = new TBuffer3D(TBuffer3DTypes::kGeneric);

   FillBuffer3D(*buf, TBuffer3D::kCore | TBuffer3D::kRawSizes | TBuffer3D::kRaw, kFALSE);

   return buf;
}
