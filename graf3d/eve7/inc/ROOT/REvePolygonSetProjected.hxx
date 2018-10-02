// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_REvePolygonSetProjected_hxx
#define ROOT_REvePolygonSetProjected_hxx

#include "ROOT/REveVector.hxx"
#include "ROOT/REveShape.hxx"

class TBuffer3D;

namespace ROOT {
namespace Experimental {

class REvePolygonSetProjected : public REveShape, public REveProjected {
private:
   REvePolygonSetProjected(const REvePolygonSetProjected &);            // Not implemented
   REvePolygonSetProjected &operator=(const REvePolygonSetProjected &); // Not implemented

protected:
   struct Polygon_t {
      Int_t fNPnts; // number of points
      Int_t *fPnts; // point indices

      Polygon_t() : fNPnts(0), fPnts(0) {}
      ~Polygon_t() { delete[] fPnts; }

      Polygon_t &operator=(const Polygon_t &x)
      {
         fNPnts = x.fNPnts;
         fPnts = x.fPnts;
         return *this;
      }

      Int_t FindPoint(Int_t pi)
      {
         for (Int_t i = 0; i < fNPnts; ++i) {
            if (fPnts[i] == pi)
               return i;
         }
         return -1;
      }
   };

   typedef std::list<Polygon_t> vpPolygon_t;
   typedef vpPolygon_t::iterator vpPolygon_i;
   typedef vpPolygon_t::const_iterator vpPolygon_ci;

private:
   std::unique_ptr<TBuffer3D> fBuff; // buffer of projectable object

   Bool_t IsFirstIdxHead(Int_t s0, Int_t s1);
   Float_t AddPolygon(std::list<Int_t, std::allocator<Int_t>> &pp, std::list<Polygon_t, std::allocator<Polygon_t>> &p);

   Int_t *ProjectAndReducePoints();
   Float_t MakePolygonsFromBP(Int_t *idxMap);
   Float_t MakePolygonsFromBS(Int_t *idxMap);

protected:
   vpPolygon_t fPols;   // polygons
   vpPolygon_t fPolsBS; // polygons build from TBuffer3D segments
   vpPolygon_t fPolsBP; // polygons build from TBuffer3D polygons

   Int_t fNPnts;      // number of reduced and projected points
   REveVector *fPnts; // reduced and projected points

   virtual void SetDepthLocal(Float_t d);

   Float_t PolygonSurfaceXY(const Polygon_t &poly) const;

public:
   REvePolygonSetProjected(const char *n = "REvePolygonSetProjected", const char *t = "");
   virtual ~REvePolygonSetProjected();

   void ComputeBBox(); // override;
   // TClass* ProjectedClass() same as for REveShape

   virtual void SetProjection(REveProjectionManager *mng, REveProjectable *model);
   virtual void UpdateProjection();
   virtual REveElement *GetProjectedAsElement() { return this; }

   void ProjectBuffer3D();

   virtual void DumpPolys() const;
   void DumpBuffer3D();

   Int_t WriteCoreJson(nlohmann::json &j, Int_t rnr_offset); // override;
   void BuildRenderData();                                   // override;

   ClassDef(REvePolygonSetProjected, 0); // Set of projected polygons with outline; typically produced from a TBuffer3D.
};

} // namespace Experimental
} // namespace ROOT

#endif
