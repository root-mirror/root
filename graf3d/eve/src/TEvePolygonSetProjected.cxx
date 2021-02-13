// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEvePolygonSetProjected.h"
#include "TEveGeoShape.h"
#include "TEveProjectionManager.h"

#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TVirtualViewer3D.h"

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
   typedef std::list<Seg_t>::iterator LSegIt_t;
}

/** \class TEvePolygonSetProjected
\ingroup TEve
A set of projected polygons.
Used for storage of projected geometrical shapes.

Internal struct Polygon_t holds only indices into the master vertex
array in TEvePolygonSetProjected.
*/

ClassImp(TEvePolygonSetProjected);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEvePolygonSetProjected::TEvePolygonSetProjected(const char* n, const char* t) :
   TEveShape(n, t),
   fBuff(0),
   fNPnts(0),
   fPnts(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TEvePolygonSetProjected::~TEvePolygonSetProjected()
{
   fPols.clear();
   if (fPnts) delete [] fPnts;
   if (fBuff) delete fBuff;
}

////////////////////////////////////////////////////////////////////////////////
/// Override of virtual method from TAttBBox.

void TEvePolygonSetProjected::ComputeBBox()
{
   if (fNPnts > 0) {
      BBoxInit();
      for (Int_t pi = 0; pi < fNPnts; ++pi)
         BBoxCheckPoint(fPnts[pi].fX, fPnts[pi].fY, fPnts[pi].fZ);
   } else {
      BBoxZero();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// This is virtual method from base-class TEveProjected.

void TEvePolygonSetProjected::SetProjection(TEveProjectionManager* mng,
                                            TEveProjectable* model)
{
   TEveProjected::SetProjection(mng, model);

   TEveGeoShape* gre = dynamic_cast<TEveGeoShape*>(model);
   fBuff = gre->MakeBuffer3D();
   CopyVizParams(gre);
}

////////////////////////////////////////////////////////////////////////////////
/// Set depth (z-coordinate) of the projected points.

void TEvePolygonSetProjected::SetDepthLocal(Float_t d)
{
   SetDepthCommon(d, this, fBBox);

   for (Int_t i = 0; i < fNPnts; ++i)
      fPnts[i].fZ = fDepth;
}

////////////////////////////////////////////////////////////////////////////////
/// This is virtual method from base-class TEveProjected.

void TEvePolygonSetProjected::UpdateProjection()
{
   if (fBuff == 0) return;

   // drop polygons and projected/reduced points
   fPols.clear();
   ProjectBuffer3D();
}

////////////////////////////////////////////////////////////////////////////////
/// Compare the two segments and check if the first index of first segment is starting.

Bool_t TEvePolygonSetProjected::IsFirstIdxHead(Int_t s0, Int_t s1)
{
   Int_t v0 = fBuff->fSegs[3*s0 + 1];
   Int_t v2 = fBuff->fSegs[3*s1 + 1];
   Int_t v3 = fBuff->fSegs[3*s1 + 2];
   return v0 != v2 && v0 != v3;
}

////////////////////////////////////////////////////////////////////////////////
/// Project and reduce buffer points.

Int_t* TEvePolygonSetProjected::ProjectAndReducePoints()
{
   TEveProjection* projection = fManager->GetProjection();

   Int_t buffN = fBuff->NbPnts();
   TEveVector*  pnts  = new TEveVector[buffN];
   for (Int_t i = 0; i < buffN; ++i)
   {
      pnts[i].Set(fBuff->fPnts[3*i],fBuff->fPnts[3*i+1], fBuff->fPnts[3*i+2]);
      projection->ProjectPoint(pnts[i].fX, pnts[i].fY, pnts[i].fZ, 0,
                               TEveProjection::kPP_Plane);
   }

   if (fPnts) delete [] fPnts;
   fNPnts=0;
   Int_t *idxMap = new Int_t[buffN];
   Int_t *ra     = new Int_t[buffN];  // list of reduced vertices
   for (UInt_t v = 0; v < (UInt_t)buffN; ++v)
   {
      idxMap[v] = -1;
      for (Int_t k = 0; k < fNPnts; ++k)
      {
         if (pnts[v].SquareDistance(pnts[ra[k]]) < TEveProjection::fgEpsSqr)
         {
            idxMap[v] = k;
            break;
         }
      }
      // have not found a point inside epsilon, add new point in scaled array
      if (idxMap[v] == -1)
      {
         idxMap[v] = fNPnts;
         ra[fNPnts] = v;
         ++fNPnts;
      }
   }

   // write the array of scaled points
   fPnts = new TEveVector[fNPnts];
   for (Int_t idx = 0; idx < fNPnts; ++idx)
   {
      Int_t i = ra[idx];
      projection->ProjectPoint(pnts[i].fX, pnts[i].fY, pnts[i].fZ, fDepth,
                               TEveProjection::kPP_Distort);
      fPnts[idx].Set(pnts[i]);
   }
   delete [] ra;
   delete [] pnts;
   // printf("reduced %d points of %d\n", fNPnts, N);

   return idxMap;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if polygon has dimensions above TEveProjection::fgEps and add it
/// to a list if it is not a duplicate.

Float_t TEvePolygonSetProjected::AddPolygon(std::list<Int_t>& pp, vpPolygon_t& pols)
{
   if (pp.size() <= 2) return 0;

   Float_t bbox[4] = { 1e6, -1e6, 1e6, -1e6 };
   for (std::list<Int_t>::iterator u = pp.begin(); u != pp.end(); ++u)
   {
      Int_t idx = *u;
      if (fPnts[idx].fX < bbox[0]) bbox[0] = fPnts[idx].fX;
      if (fPnts[idx].fX > bbox[1]) bbox[1] = fPnts[idx].fX;

      if (fPnts[idx].fY < bbox[2]) bbox[2] = fPnts[idx].fY;
      if (fPnts[idx].fY > bbox[3]) bbox[3] = fPnts[idx].fY;
   }
   Float_t eps = 2*TEveProjection::fgEps;
   if ((bbox[1]-bbox[0]) < eps || (bbox[3]-bbox[2]) < eps) return 0;

   // Duplication
   for (vpPolygon_i poi = pols.begin(); poi != pols.end(); ++poi)
   {
      Polygon_t& refP = *poi;

      if ((Int_t) pp.size() != refP.fNPnts)
         continue;

      Int_t start_idx = refP.FindPoint(pp.front());
      if (start_idx < 0)
            continue;
      if (++start_idx >= refP.fNPnts) start_idx = 0;

      // Same orientation duplicate
      {
         std::list<Int_t>::iterator u = ++pp.begin();
         Int_t pidx = start_idx;
         while (u != pp.end())
         {
            if ((*u) != refP.fPnts[pidx])
               break;
            ++u;
            if (++pidx >= refP.fNPnts) pidx = 0;
         }
         if (u == pp.end()) return 0;
      }
      // Inverse orientation duplicate
      {
         std::list<Int_t>::iterator u = --pp.end();
         Int_t pidx = start_idx;
         while (u != pp.begin())
         {
            if ((*u) != refP.fPnts[pidx])
               break;
            --u;
            if (++pidx >= refP.fNPnts) pidx = 0;
         }
         if (u == pp.begin()) return 0;
      }
   }

   Int_t *pv    = new Int_t[pp.size()];
   Int_t  count = 0;
   for (std::list<Int_t>::iterator u = pp.begin(); u != pp.end(); ++u)
   {
      pv[count] = *u;
      ++count;
   }

   pols.push_back(Polygon_t());
   pols.back().fNPnts =  pp.size();
   pols.back().fPnts  = &pv[0];

   return (bbox[1]-bbox[0]) * (bbox[3]-bbox[2]);
}

////////////////////////////////////////////////////////////////////////////////
/// Build polygons from list of buffer polygons.

Float_t TEvePolygonSetProjected::MakePolygonsFromBP(Int_t* idxMap)
{
   TEveProjection* projection = fManager->GetProjection();
   Int_t   *bpols = fBuff->fPols;
   Float_t  surf  = 0; // surface of projected polygons
   for (UInt_t pi = 0; pi < fBuff->NbPols(); ++pi)
   {
      std::list<Int_t> pp; // points in current polygon
      UInt_t  segN =  bpols[1];
      Int_t  *seg  = &bpols[2];
      // start idx in the fist segment depends of second segment
      Int_t   tail, head;
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
      pp.push_back(head);
      // printf("start idx head %d, tail %d\n", head, tail);
      LSeg_t segs;
      for (UInt_t s = 1; s < segN; ++s)
         segs.push_back(Seg_t(fBuff->fSegs[3*seg[s] + 1],fBuff->fSegs[3*seg[s] + 2]));

      for (LSegIt_t it = segs.begin(); it != segs.end(); ++it)
      {
         Int_t mv1 = idxMap[(*it).fV1];
         Int_t mv2 = idxMap[(*it).fV2];

         if ( ! projection->AcceptSegment(fPnts[mv1], fPnts[mv2], TEveProjection::fgEps))
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

Float_t TEvePolygonSetProjected::MakePolygonsFromBS(Int_t* idxMap)
{
   LSeg_t   segs;
   LSegIt_t it;
   Float_t  surf = 0; // surface of projected polygons
   TEveProjection *projection = fManager->GetProjection();
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
      for (it = segs.begin(); it != segs.end(); ++it)
      {
         Int_t vv1 = (*it).fV1;
         Int_t vv2 = (*it).fV2;
         if((vv1 == vor1 && vv2 == vor2) || (vv1 == vor2 && vv2 == vor1))
         {
            duplicate = kTRUE;
            continue;
         }
      }
      if (duplicate == kFALSE && projection->AcceptSegment(fPnts[vor1], fPnts[vor2], TEveProjection::fgEps))
         segs.push_back(Seg_t(vor1, vor2));
   }

   while ( ! segs.empty())
   {
      std::list<Int_t> pp; // points in current polygon
      pp.push_back(segs.front().fV1);
      Int_t tail = segs.front().fV2;
      segs.pop_front();
      Bool_t match = kTRUE;
      while (match && ! segs.empty())
      {
         for (LSegIt_t k = segs.begin(); k != segs.end(); ++k)
         {
            Int_t cv1 = (*k).fV1;
            Int_t cv2 = (*k).fV2;
            if (cv1 == tail || cv2 == tail)
            {
               pp.push_back(tail);
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

void  TEvePolygonSetProjected::ProjectBuffer3D()
{
   // create map from original to projected and reduced point needed only for geometry
   Int_t* idxMap = ProjectAndReducePoints();

   TEveProjection::EGeoMode_e mode = fManager->GetProjection()->GetGeoMode();
   switch (mode)
   {
      case TEveProjection::kGM_Polygons :
      {
         MakePolygonsFromBP(idxMap);
         fPolsBP.swap(fPols);
         break;
      }
      case TEveProjection::kGM_Segments :
      {
         MakePolygonsFromBS(idxMap);
         fPolsBS.swap(fPols);
         break;
      }
      case TEveProjection::kGM_Unknown:
      {
         // take projection with largest surface
        Float_t surfBP = MakePolygonsFromBP(idxMap);
        Float_t surfBS = MakePolygonsFromBS(idxMap);
         if (surfBS < surfBP)
         {
            fPolsBP.swap(fPols);
            fPolsBS.clear();
         }
         else
         {
            fPolsBS.swap(fPols);
            fPolsBP.clear();
         }
      }
      default:
         break;
   }

   delete [] idxMap;
   ResetBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate XY surface of a polygon.

Float_t TEvePolygonSetProjected::PolygonSurfaceXY(const TEvePolygonSetProjected::Polygon_t& p) const
{
   Float_t surf = 0;
   Int_t nPnts = p.fNPnts;
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

void TEvePolygonSetProjected::DumpPolys() const
{
   printf("TEvePolygonSetProjected %d polygons\n", (Int_t)fPols.size());
   Int_t cnt = 0;
   for (vpPolygon_ci i = fPols.begin(); i!= fPols.end(); i++)
   {
      Int_t nPnts = (*i).fNPnts;
      printf("Points of polygon %d [Np = %d]:\n", ++cnt, nPnts);
      for (Int_t vi = 0; vi<nPnts; ++vi) {
         Int_t pi = (*i).fPnts[vi];
         printf("  (%f, %f, %f)", fPnts[pi].fX, fPnts[pi].fY, fPnts[pi].fZ);
      }
      printf(", surf=%f\n", PolygonSurfaceXY(*i));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Dump information about currently projected buffer.

void TEvePolygonSetProjected::DumpBuffer3D()
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
