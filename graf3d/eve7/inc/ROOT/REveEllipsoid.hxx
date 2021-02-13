// @(#)root/eve7:$Id$
// Author: Matevz Tadel, Alja Tadel

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveEllipsoid
#define ROOT7_REveEllipsoid

#include <ROOT/REveShape.hxx>
#include <ROOT/REveVector.hxx>
#include <ROOT/REveProjectionBases.hxx>
#include <ROOT/REveStraightLineSet.hxx>

namespace ROOT {
namespace Experimental {

//------------------------------------------------------------------------------
// REveEllipsoid
//------------------------------------------------------------------------------

class REveEllipsoid : public REveStraightLineSet
{
   friend class REveEllipsoidProjected;

private:
   REveEllipsoid(const REveEllipsoid &) = delete;
   REveEllipsoid &operator=(const REveEllipsoid &) = delete;

protected:
   REveVector fV0;
   REveVector fV1;
   REveVector fV2;

   float      fPhiStep;
   void DrawArch(float phiStart, float phiEnd, float phiStep, REveVector& v0,  REveVector& v1, REveVector& v2);

public:
   REveEllipsoid(const std::string &n = "REveJetConeProjected", const std::string& t = "");
   virtual ~REveEllipsoid() {};

   virtual void Outline();
   void SetBaseVectors(REveVector& v0, REveVector& v1, REveVector& v3);
   void SetPhiStep(float ps);

   Int_t WriteCoreJson(nlohmann::json &j, Int_t rnr_offset) override;

   TClass *ProjectedClass(const REveProjection *p) const override;
};

//------------------------------------------------------------------------------
// REveEllipsoidProjected
//------------------------------------------------------------------------------

class REveEllipsoidProjected :  public REveStraightLineSetProjected
{
private:
   REveEllipsoidProjected(const REveEllipsoidProjected &) = delete;
   REveEllipsoidProjected &operator=(const REveEllipsoidProjected &) = delete;

   void DrawArchProjected(float phiStart, float phiEnd, float phiStep, REveVector& v0,  REveVector& v1, REveVector& v2);
   void GetSurfaceSize(REveVector& p1, REveVector& p2);

   REveVector fMV0;
   REveVector fMV1;

   std::vector <REveVector> fArchPnts;
   float GetEllipseSurface (const REveVector& v1, const REveVector& v2);

public:
   REveEllipsoidProjected(const std::string &n = "REveEllipsoidProjected", const std::string& t = "");
   virtual ~REveEllipsoidProjected();

   void BuildRenderData() override;

   Int_t WriteCoreJson(nlohmann::json &j, Int_t rnr_offset) override;

   virtual void OutlineProjected();
   virtual void SetProjection(REveProjectionManager *mng, REveProjectable *model) override;
   void UpdateProjection() override;
   REveElement *GetProjectedAsElement() override { return this; }
};

} // namespace Experimental
} // namespace ROOT

#endif
