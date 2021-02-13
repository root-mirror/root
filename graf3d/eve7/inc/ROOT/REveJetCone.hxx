// @(#)root/eve7:$Id$
// Author: Matevz Tadel, Jochen Thaeder 2009, 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveJetCone
#define ROOT7_REveJetCone

#include <ROOT/REveShape.hxx>
#include <ROOT/REveVector.hxx>

namespace ROOT {
namespace Experimental {

//------------------------------------------------------------------------------
// REveJetCone
//------------------------------------------------------------------------------

class REveJetCone : public REveShape,
                    public REveProjectable
{
   friend class REveJetConeProjected;

private:
   REveJetCone(const REveJetCone &) = delete;
   REveJetCone &operator=(const REveJetCone &) = delete;

protected:
   REveVector fApex;   // Apex of the cone.
   REveVector fAxis;   // Axis of the cone.
   REveVector fLimits; // Border of Barrel/Cylinder to cut the cone.
   Float_t fThetaC;    // Transition theta
   Float_t fEta, fPhi;
   Float_t fDEta, fDPhi;
   Int_t fNDiv;

   REveVector CalcEtaPhiVec(Float_t eta, Float_t phi) const;
   REveVector CalcBaseVec(Float_t eta, Float_t phi) const;
   REveVector CalcBaseVec(Float_t alpha) const;
   Bool_t IsInTransitionRegion() const;

public:
   REveJetCone(const Text_t *n = "REveJetCone", const Text_t *t = "");
   virtual ~REveJetCone() {}

   Int_t WriteCoreJson(nlohmann::json &j, Int_t rnr_offset) override;
   void BuildRenderData() override;

   void ComputeBBox() override;
   TClass *ProjectedClass(const REveProjection *p) const override;

   void SetApex(const REveVector &a) { fApex = a; }
   void SetCylinder(Float_t r, Float_t z)
   {
      fLimits.Set(0, r, z);
      fThetaC = fLimits.Theta();
   }
   void SetRadius(Float_t r)
   {
      fLimits.Set(r, 0, 0);
      fThetaC = 10;
   }

   Int_t GetNDiv() const { return fNDiv; }
   void SetNDiv(Int_t n);

   Int_t AddCone(Float_t eta, Float_t phi, Float_t cone_r, Float_t length = 0);
   Int_t AddEllipticCone(Float_t eta, Float_t phi, Float_t reta, Float_t rphi, Float_t length = 0);
};


//------------------------------------------------------------------------------
// REveJetConeProjected
//------------------------------------------------------------------------------

class REveJetConeProjected : public REveShape,
                             public REveProjected
{
private:
   REveJetConeProjected(const REveJetConeProjected &) = delete;
   REveJetConeProjected &operator=(const REveJetConeProjected &) = delete;

protected:
   void SetDepthLocal(Float_t d) override;

public:
   REveJetConeProjected(const std::string &n = "REveJetConeProjected", const std::string& t = "");
   virtual ~REveJetConeProjected();

   void BuildRenderData() override;

   // For TAttBBox:
   void ComputeBBox() override;

   // Projected:
   void SetProjection(REveProjectionManager *mng, REveProjectable *model) override;
   void UpdateProjection() override;

   REveElement *GetProjectedAsElement() override { return this; }
};

} // namespace Experimental
} // namespace ROOT

#endif
