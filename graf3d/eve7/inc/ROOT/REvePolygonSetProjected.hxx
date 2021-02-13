// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REvePolygonSetProjected
#define ROOT7_REvePolygonSetProjected

#include <ROOT/REveVector.hxx>
#include <ROOT/REveShape.hxx>

#include <vector>

class TBuffer3D;

namespace ROOT {
namespace Experimental {

// =========================================================================
// REvePolygonSetProjected
// Set of projected polygons with outline; typically produced from a TBuffer3D.
// =========================================================================

class REvePolygonSetProjected : public REveShape,
                                public REveProjected
{
private:
   REvePolygonSetProjected(const REvePolygonSetProjected &) = delete;
   REvePolygonSetProjected &operator=(const REvePolygonSetProjected &) = delete;

protected:
   struct Polygon_t {
      std::vector<UInt_t> fPnts; // point indices

      Polygon_t() = default;
      Polygon_t(std::vector<UInt_t> &&pnts) : fPnts(pnts){};
      ~Polygon_t() = default;

      int NPoints() const { return (UInt_t)fPnts.size(); }

      int FindPoint(int pi)
      {
         auto dist = std::distance(fPnts.begin(), std::find(fPnts.begin(), fPnts.end(), pi));
         return (dist >= (int) fPnts.size()) ? -1 : (int) dist;
      }
   };

   typedef std::list<Polygon_t> vpPolygon_t;

private:
   std::unique_ptr<TBuffer3D> fBuff; // buffer of projectable object

   Bool_t IsFirstIdxHead(Int_t s0, Int_t s1);
   Float_t AddPolygon(std::list<UInt_t> &pp, std::list<Polygon_t> &p);

   std::vector<UInt_t> ProjectAndReducePoints();
   Float_t MakePolygonsFromBP(std::vector<UInt_t> &idxMap);
   Float_t MakePolygonsFromBS(std::vector<UInt_t> &idxMap);

protected:
   vpPolygon_t fPols;   ///<! polygons
   vpPolygon_t fPolsBS; ///<! polygons build from TBuffer3D segments
   vpPolygon_t fPolsBP; ///<! polygons build from TBuffer3D polygons

   std::vector<REveVector> fPnts; ///<! reduced and projected points

   void SetDepthLocal(Float_t d) override;

   Float_t PolygonSurfaceXY(const Polygon_t &poly) const;

public:
   REvePolygonSetProjected(const std::string &n = "REvePolygonSetProjected", const std::string &t = "");
   virtual ~REvePolygonSetProjected();

   void ComputeBBox() override;
   // TClass* ProjectedClass() same as for REveShape

   void SetProjection(REveProjectionManager *mng, REveProjectable *model) override;
   void UpdateProjection() override;
   REveElement *GetProjectedAsElement() override { return this; }

   void ProjectBuffer3D();

   virtual void DumpPolys() const;
   void DumpBuffer3D();

   Int_t WriteCoreJson(nlohmann::json &j, Int_t rnr_offset) override;
   void BuildRenderData() override;
};

} // namespace Experimental
} // namespace ROOT

#endif
