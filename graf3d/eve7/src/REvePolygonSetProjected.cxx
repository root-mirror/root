// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/REvePolygonSetProjected.hxx>
#include <ROOT/REveGeoShape.hxx>
#include <ROOT/REveProjectionManager.hxx>
#include <ROOT/REveGluTess.hxx>
#include <ROOT/REveRenderData.hxx>

#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"

#include <cassert>

using namespace ROOT::Experimental;

namespace
{
   struct Seg_t
   {
      // Helper class for building 2D polygons from TBuffer3D.
      Int_t fV1;
      Int_t fV2;

      Seg_t(Int_t i1=-1, Int_t i2=-1) : fV1(i1), fV2(i2) {}
   };

   typedef std::list<Seg_t>           LSeg_t;
}

/** \class REvePolygonSetProjected
\ingroup REve
A set of projected polygons.
Used for storage of projected geometrical shapes.

Internal struct Polygon_t holds only indices into the master vertex
array in REvePolygonSetProjected.
*/

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REvePolygonSetProjected::REvePolygonSetProjected(const std::string &n, const std::string &t) :
   REveShape(n, t),
   fBuff(),
   fPnts()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

REvePolygonSetProjected::~REvePolygonSetProjected()
{
   fPols.clear();
}

////////////////////////////////////////////////////////////////////////////////
/// Fill core part of JSON representation.

Int_t REvePolygonSetProjected::WriteCoreJson(nlohmann::json& j, Int_t rnr_offset)
{
   Int_t ret = REveElement::WriteCoreJson(j, rnr_offset);

   j["fNPnts"] = fPnts.size();

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Crates representation for rendering.
/// This is complicated as we need to:
/// - generate outlines;
/// - convert polygons to triangles.
/// ??? Should we check if polygons are front facing? It was not done in old EVE,
/// just GL_FRONT_AND_BACK was used on top of gluTess.

void REvePolygonSetProjected::BuildRenderData()
{
   fRenderData = std::make_unique<REveRenderData>("makePolygonSetProjected", 3 * fPnts.size());

   Int_t n_pols = fPols.size();
   Int_t n_poly_info = 0;
   for (auto &p : fPols) n_poly_info += 1 + p.NPoints();

   std::vector<Double_t> verts;
   verts.reserve(3 * fPnts.size());
   std::vector<UInt_t>    polys;
   polys.reserve(n_poly_info);

   for (auto &p : fPols)
   {
      polys.emplace_back(p.NPoints());
      polys.insert(polys.end(), p.fPnts.begin(), p.fPnts.end());
   }

   for (unsigned i = 0; i < fPnts.size(); ++i)
   {
      verts.push_back(fPnts[i].fX);
      verts.push_back(fPnts[i].fY);
      verts.push_back(fPnts[i].fZ);
      fRenderData->PushV(fPnts[i]);
   }

   Int_t n_trings = 0;
   {
      EveGlu::TriangleCollector tc;

      tc.ProcessData(verts, polys, n_pols);

      polys.swap(tc.RefPolyDesc());
      n_trings = tc.GetNTrianlges();
   }

   // Calculate size of index buffer.
   Int_t n_idxbuff = 2 + 3 * n_trings + n_pols + n_poly_info;
   fRenderData->Reserve(0,0,n_idxbuff);

   assert(n_trings * 4 == (int)polys.size());

   // Export triangles.
   fRenderData->PushI(REveRenderData::GL_TRIANGLES);
   fRenderData->PushI(n_trings);
   for (int i = 0; i < n_trings; ++i)
   {
      fRenderData->PushI(&polys[i*4 + 1], 3);
   }

   assert (fRenderData->SizeI() == 2 + 3 * n_trings);

   // Export outlines.
   for (auto &p : fPols)
   {
      fRenderData->PushI(REveRenderData::GL_LINE_LOOP);
      fRenderData->PushI(p.NPoints());
      fRenderData->PushI(p.fPnts);
   }

   assert (fRenderData->SizeI() == n_idxbuff);
}

////////////////////////////////////////////////////////////////////////////////
/// Override of virtual method from TAttBBox.

void REvePolygonSetProjected::ComputeBBox()
{
   if (fPnts.size() > 0) {
      BBoxInit();
      for (unsigned pi = 0; pi < fPnts.size(); ++pi)
         BBoxCheckPoint(fPnts[pi].fX, fPnts[pi].fY, fPnts[pi].fZ);
   } else {
      BBoxZero();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// This is virtual method from base-class REveProjected.

void REvePolygonSetProjected::SetProjection(REveProjectionManager* mng,
                                            REveProjectable* model)
{
   REveProjected::SetProjection(mng, model);

   REveGeoShape* gre = dynamic_cast<REveGeoShape*>(model);
   fBuff = gre->MakeBuffer3D();
   CopyVizParams(gre);
}

////////////////////////////////////////////////////////////////////////////////
/// Set depth (z-coordinate) of the projected points.

void REvePolygonSetProjected::SetDepthLocal(Float_t d)
{
   SetDepthCommon(d, this, fBBox);

   for (unsigned i = 0; i < fPnts.size(); ++i)
      fPnts[i].fZ = fDepth;
}

////////////////////////////////////////////////////////////////////////////////
/// This is virtual method from base-class REveProjected.

void REvePolygonSetProjected::UpdateProjection()
{
   if (!fBuff) return;

   // drop polygons and projected/reduced points
   fPols.clear();
   ProjectBuffer3D();
}

////////////////////////////////////////////////////////////////////////////////
/// Compare the two segments and check if the first index of first segment is starting.

Bool_t REvePolygonSetProjected::IsFirstIdxHead(Int_t s0, Int_t s1)
{
   Int_t v0 = fBuff->fSegs[3*s0 + 1];
   Int_t v2 = fBuff->fSegs[3*s1 + 1];
   Int_t v3 = fBuff->fSegs[3*s1 + 2];
   return v0 != v2 && v0 != v3;
}

////////////////////////////////////////////////////////////////////////////////
/// Project and reduce buffer points.

std::vector<UInt_t> REvePolygonSetProjected::ProjectAndReducePoints()
{
   REveProjection* projection = fManager->GetProjection();

   Int_t buffN = fBuff->NbPnts();
   std::vector<REveVector> pnts; pnts.resize(buffN);
   for (Int_t i = 0; i < buffN; ++i)
   {
      pnts[i].Set(fBuff->fPnts[3*i],fBuff->fPnts[3*i+1], fBuff->fPnts[3*i+2]);
      projection->ProjectPoint(pnts[i].fX, pnts[i].fY, pnts[i].fZ, 0,
                               REveProjection::kPP_Plane);
   }

   int npoints = 0;
   std::vector<UInt_t> idxMap;
   idxMap.resize(buffN);

   std::vector<int> ra;
   ra.resize(buffN);  // list of reduced vertices
   for (UInt_t v = 0; v < (UInt_t)buffN; ++v)
   {
      bool found = false;
      for (Int_t k = 0; k < npoints; ++k)
      {
         if (pnts[v].SquareDistance(pnts[ra[k]]) < REveProjection::fgEpsSqr)
         {
            idxMap[v] = k;
            found = true;
            break;
         }
      }
      // have not found a point inside epsilon, add new point in scaled array
      if (!found)
      {
         idxMap[v] = npoints;
         ra[npoints] = v;
         ++npoints;
      }
   }

   // write the array of scaled points
   fPnts.resize(npoints);
   for (Int_t idx = 0; idx < npoints; ++idx)
   {
      Int_t i = ra[idx];
      projection->ProjectPoint(pnts[i].fX, pnts[i].fY, pnts[i].fZ, fDepth,
                               REveProjection::kPP_Distort);
      fPnts[idx].Set(pnts[i]);
   }
   // printf("reduced %d points of %d\n", fNPnts, N);

   return idxMap;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if polygon has dimensions above REveProjection::fgEps and add it
/// to a list if it is not a duplicate.

Float_t REvePolygonSetProjected::AddPolygon(std::list<UInt_t> &pp, vpPolygon_t &pols)
{
   if (pp.size() <= 2) return 0;

   Float_t bbox[4] = { 1e6, -1e6, 1e6, -1e6 };
   for (auto &&idx: pp)
   {
      if (fPnts[idx].fX < bbox[0]) bbox[0] = fPnts[idx].fX;
      if (fPnts[idx].fX > bbox[1]) bbox[1] = fPnts[idx].fX;

      if (fPnts[idx].fY < bbox[2]) bbox[2] = fPnts[idx].fY;
      if (fPnts[idx].fY > bbox[3]) bbox[3] = fPnts[idx].fY;
   }
   Float_t eps = 2*REveProjection::fgEps;
   if ((bbox[1]-bbox[0]) < eps || (bbox[3]-bbox[2]) < eps) return 0;

   // Duplication
   for (auto &&refP : pols)
   {
      if ((Int_t) pp.size() != refP.NPoints())
         continue;

      int start_idx = refP.FindPoint(pp.front());
      if (start_idx < 0)
            continue;
      if (++start_idx >= refP.NPoints()) start_idx = 0;

      // Same orientation duplicate
      {
         auto u = ++pp.begin();
         Int_t pidx = start_idx;
         while (u != pp.end())
         {
            if ((*u) != refP.fPnts[pidx])
               break;
            ++u;
            if (++pidx >= refP.NPoints()) pidx = 0;
         }
         if (u == pp.end()) return 0;
      }
      // Inverse orientation duplicate
      {
         auto u = --pp.end();
         Int_t pidx = start_idx;
         while (u != pp.begin())
         {
            if ((*u) != refP.fPnts[pidx])
               break;
            --u;
            if (++pidx >= refP.NPoints()) pidx = 0;
         }
         if (u == pp.begin()) return 0;
      }
   }

   std::vector<UInt_t> pv(pp.size(), 0);
   int  count = 0;
   for (auto &&u : pp) {
      pv[count++] = u;
   }

   pols.emplace_back(std::move(pv));

   return (bbox[1]-bbox[0]) * (bbox[3]-bbox[2]);
}

////////////////////////////////////////////////////////////////////////////////
/// Build polygons from list of buffer polygons.

Float_t REvePolygonSetProjected::MakePolygonsFromBP(std::vector<UInt_t> &idxMap)
{
   REveProjection* projection = fManager->GetProjection();
   Int_t   *bpols = fBuff->fPols;
   Float_t  surf  = 0; // surface of projected polygons
   for (UInt_t pi = 0; pi < fBuff->NbPols(); ++pi)
   {
      std::list<UInt_t> pp; // points in current polygon
      UInt_t  segN =  bpols[1];
      Int_t  *seg  = &bpols[2];
      // start idx in the fist segment depends of second segment
      UInt_t   tail, head;
      if (IsFirstIdxHead(seg[0], seg[1]))
      {
         head = idxMap[fBuff->fSegs[3*seg[0] + 1]];
         tail = idxMap[fBuff->fSegs[3*seg[0] + 2]];
      }
      else
      {
         head = idxMap[fBuff->fSegs[3*seg[0] + 2]];
         tail = idxMap[fBuff->fSegs[3*seg[0] + 1]];
      }
      pp.emplace_back(head);
      // printf("start idx head %d, tail %d\n", head, tail);
      LSeg_t segs;
      for (UInt_t s = 1; s < segN; ++s)
         segs.emplace_back(fBuff->fSegs[3*seg[s] + 1],fBuff->fSegs[3*seg[s] + 2]);

      for (auto &it: segs)
      {
         UInt_t mv1 = idxMap[it.fV1];
         UInt_t mv2 = idxMap[it.fV2];

         if ( ! projection->AcceptSegment(fPnts[mv1], fPnts[mv2], REveProjection::fgEps))
         {
            pp.clear();
            break;
         }
         if (tail != pp.back()) pp.push_back(tail);
         tail = (mv1 == tail) ? mv2 : mv1;
      }

      if ( ! pp.empty())
      {
         // DirectDraw() implementation: last and first vertices should not be equal
         if (pp.front() == pp.back()) pp.pop_front();
         surf += AddPolygon(pp, fPolsBP);
      }
      bpols += (segN+2);
   }
   return surf;
}

////////////////////////////////////////////////////////////////////////////////
/// Build polygons from the set of buffer segments.
/// First creates a segment pool according to reduced and projected points
/// and then build polygons from the pool.

Float_t REvePolygonSetProjected::MakePolygonsFromBS(std::vector<UInt_t> &idxMap)
{
   LSeg_t   segs;
   Float_t  surf = 0; // surface of projected polygons
   REveProjection *projection = fManager->GetProjection();
   for (UInt_t s = 0; s < fBuff->NbSegs(); ++s)
   {
      Bool_t duplicate = kFALSE;
      Int_t vo1,  vo2;  // idx from fBuff segment
      Int_t vor1, vor2; // mapped idx
      vo1 =  fBuff->fSegs[3*s + 1];
      vo2 =  fBuff->fSegs[3*s + 2]; //... skip color info
      vor1 = idxMap[vo1];
      vor2 = idxMap[vo2];
      if (vor1 == vor2) continue;
      // check duplicate
      for (auto &seg: segs)
      {
         Int_t vv1 = seg.fV1;
         Int_t vv2 = seg.fV2;
         if((vv1 == vor1 && vv2 == vor2) || (vv1 == vor2 && vv2 == vor1))
         {
            duplicate = kTRUE;
            continue;
         }
      }
      if (duplicate == kFALSE && projection->AcceptSegment(fPnts[vor1], fPnts[vor2], REveProjection::fgEps))
         segs.emplace_back(vor1, vor2);
   }

   while (!segs.empty())
   {
      std::list<UInt_t> pp; // points in current polygon
      pp.push_back(segs.front().fV1);
      UInt_t tail = segs.front().fV2;
      segs.pop_front();
      Bool_t match = kTRUE;
      while (match && ! segs.empty())
      {
         for (auto k = segs.begin(); k != segs.end(); ++k)
         {
            UInt_t cv1 = (*k).fV1;
            UInt_t cv2 = (*k).fV2;
            if (cv1 == tail || cv2 == tail)
            {
               pp.emplace_back(tail);
               tail = (cv1 == tail) ? cv2 : cv1;
               segs.erase(k);
               match = kTRUE;
               break;
            }
            else
            {
               match = kFALSE;
            }
         } // end for loop in the segment pool
         if (tail == pp.front())
            break;
      }
      surf += AddPolygon(pp, fPolsBS);
   }
   return surf;
}

////////////////////////////////////////////////////////////////////////////////
/// Project current buffer.

void  REvePolygonSetProjected::ProjectBuffer3D()
{
   // create map from original to projected and reduced point needed only for geometry
   auto idxMap = ProjectAndReducePoints();

   REveProjection::EGeoMode_e mode = fManager->GetProjection()->GetGeoMode();
   switch (mode) {
      case REveProjection::kGM_Polygons: {
         MakePolygonsFromBP(idxMap);
         fPolsBP.swap(fPols);
         break;
      }
      case REveProjection::kGM_Segments: {
         MakePolygonsFromBS(idxMap);
         fPolsBS.swap(fPols);
         break;
      }
      case REveProjection::kGM_Unknown: {
         // take projection with largest surface
         Float_t surfBP = MakePolygonsFromBP(idxMap);
         Float_t surfBS = MakePolygonsFromBS(idxMap);
         if (surfBS < surfBP) {
            fPolsBP.swap(fPols);
            fPolsBS.clear();
         } else {
            fPolsBS.swap(fPols);
            fPolsBP.clear();
         }
         break;
      }
      default: break;
   }

   ResetBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate XY surface of a polygon.

Float_t REvePolygonSetProjected::PolygonSurfaceXY(const REvePolygonSetProjected::Polygon_t &p) const
{
   Float_t surf = 0;
   Int_t nPnts = p.NPoints();
   for (Int_t i = 0; i < nPnts - 1; ++i)
   {
      Int_t a = p.fPnts[i];
      Int_t b = p.fPnts[i+1];
      surf += fPnts[a].fX * fPnts[b].fY - fPnts[a].fY * fPnts[b].fX;
   }
   return 0.5f * TMath::Abs(surf);
}

////////////////////////////////////////////////////////////////////////////////
/// Dump information about built polygons.

void REvePolygonSetProjected::DumpPolys() const
{
   printf("REvePolygonSetProjected %d polygons\n", (Int_t)fPols.size());
   Int_t cnt = 0;
   for ( auto &pol : fPols)
   {
      Int_t nPnts = pol.NPoints();
      printf("Points of polygon %d [Np = %d]:\n", ++cnt, nPnts);
      for (Int_t vi = 0; vi<nPnts; ++vi) {
         Int_t pi = pol.fPnts[vi];
         printf("  (%f, %f, %f)", fPnts[pi].fX, fPnts[pi].fY, fPnts[pi].fZ);
      }
      printf(", surf=%f\n", PolygonSurfaceXY(pol));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Dump information about currently projected buffer.

void REvePolygonSetProjected::DumpBuffer3D()
{
   Int_t* bpols = fBuff->fPols;

   for (UInt_t pi = 0; pi< fBuff->NbPols(); ++pi)
   {
      UInt_t segN = bpols[1];
      printf("%d polygon of %d has %d segments \n", pi, fBuff->NbPols(), segN);

      Int_t* seg =  &bpols[2];
      for (UInt_t a=0; a<segN; ++a)
      {
         Int_t a1 = fBuff->fSegs[3*seg[a] + 1];
         Int_t a2 = fBuff->fSegs[3*seg[a] + 2];
         printf("(%d, %d) \n", a1, a2);
         printf("ORIG points :(%f, %f, %f)  (%f, %f, %f)\n",
                fBuff->fPnts[3*a1],fBuff->fPnts[3*a1+1], fBuff->fPnts[3*a1+2],
                fBuff->fPnts[3*a2],fBuff->fPnts[3*a2+1], fBuff->fPnts[3*a2+2]);
      }
      printf("\n");
      bpols += (segN+2);
   }
}
