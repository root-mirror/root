// @(#)root/eve7:$Id$
// Author: Matevz Tadel, Jochen Thaeder 2009, 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/REveJetCone.hxx>
#include <ROOT/REveTrans.hxx>
#include <ROOT/REveProjectionManager.hxx>
#include <ROOT/REveRenderData.hxx>

#include "TMath.h"
#include "TClass.h"

#include <cassert>

using namespace ROOT::Experimental;

/** \class REveJetCone
\ingroup REve
Draws a jet cone with leading particle is specified in (eta,phi) and
cone radius is given.

If Apex is not set, default is (0.,0.,0.)
In case of cylinder was set, cone is cut at the cylinder edges.

Example :
~~~ {.cpp}
  Float_t coneEta    = r.Uniform(-0.9, 0.9);
  Float_t conePhi    = r.Uniform(0.0, TwoPi() );
  Float_t coneRadius = 0.4;

  REveJetCone* jetCone = new REveJetCone("JetCone");
  jetCone->SetCylinder(250, 250);
  if (jetCone->AddCone(coneEta, conePhi, coneRadius) != -1)
    gEve->AddElement(jetCone);
~~~

#### Implementation notes

REveVector fLimits encodes the following information:
  - fY, fZ:  barrel radius and endcap z-position;
             if both are 0, fX encodes the spherical radius
  - fX    :  scaling for length of the cone
*/

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REveJetCone::REveJetCone(const Text_t* n, const Text_t* t) :
   REveShape(n, t),
   fApex(),
   fLimits(), fThetaC(10),
   fEta(0), fPhi(0), fDEta(0), fDPhi(0), fNDiv(36)
{
   fFillColor = kGreen;
}

////////////////////////////////////////////////////////////////////////////////
/// Set Number of Divisions
////////////////////////////////////////////////////////////////////////////////
void  REveJetCone::SetNDiv(Int_t n)
{
   fNDiv = TMath::Max(4, n);
   if (fNDiv % 4 > 0) { fNDiv += 4 - fNDiv % 4; };
   StampObjProps();
}

////////////////////////////////////////////////////////////////////////////////
/// Fill core part of JSON representation.

Int_t REveJetCone::WriteCoreJson(nlohmann::json &j, Int_t rnr_offset)
{
   Int_t ret = REveElement::WriteCoreJson(j, rnr_offset);

   j["fMainColor"] = GetFillColor();
   j["fLineColor"] = GetLineColor();
   j["fNDiv"] = GetNDiv();

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Crates 3D point array for rendering.

void REveJetCone::BuildRenderData()
{
   assert(fNDiv > 2);

   const Int_t  NP = 1 + fNDiv;

   fRenderData = std::make_unique<REveRenderData>("makeJet", 3 * NP);

   fRenderData->PushV(fApex);

   Float_t angle_step = TMath::TwoPi() / fNDiv;
   Float_t angle      = 0;
   for (Int_t i = 0; i < fNDiv; ++i, angle += angle_step)
   {
      fRenderData->PushV( CalcBaseVec(angle) );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute bounding-box of the data.

void REveJetCone::ComputeBBox()
{
   BBoxInit();
   BBoxCheckPoint(fApex);
   BBoxCheckPoint(CalcBaseVec(0));
   BBoxCheckPoint(CalcBaseVec(TMath::PiOver2()));
   BBoxCheckPoint(CalcBaseVec(TMath::Pi()));
   BBoxCheckPoint(CalcBaseVec(TMath::Pi() + TMath::PiOver2()));
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual from REveProjectable, returns REveJetConeProjected class.

TClass* REveJetCone::ProjectedClass(const REveProjection*) const
{
   return TClass::GetClass<REveJetConeProjected>();
}

////////////////////////////////////////////////////////////////////////////////
/// Add jet cone.
/// parameters are :
/// - (eta,phi)    : of the center/leading particle
/// - cone_r       : cone radius in eta-phi space
/// - length       : length of the cone
///   - if cylinder is set and length is adapted to cylinder.
///      - if length is given, it will be used as scalar factor
///   - if cylinder is not set, length is used as length of the cone
/// Return 0 on success.

Int_t REveJetCone::AddCone(Float_t eta, Float_t phi, Float_t cone_r, Float_t length)
{
   return AddEllipticCone(eta, phi, cone_r, cone_r, length);
}

////////////////////////////////////////////////////////////////////////////////
/// Add jet cone.
/// parameters are :
/// - (eta,phi)    : of the center/leading particle
/// - (reta, rphi) : radius of cone in eta-phi space
/// - length       : length of the cone
///   - if cylinder is set and length is adapted to cylinder.
///      - if length is given, it will be used as scalar factor
///   - if cylinder is not set, length is used as length of the cone
/// Returns 0 on success.

Int_t REveJetCone::AddEllipticCone(Float_t eta, Float_t phi, Float_t reta, Float_t rphi, Float_t length)
{
   if (length != 0) fLimits.fX = length;

   if (fLimits.IsZero())
      return -1;

   fEta = eta; fPhi = phi; fDEta = reta; fDPhi = rphi;

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill REveVector with eta and phi, magnitude 1.

REveVector REveJetCone::CalcEtaPhiVec(Float_t eta, Float_t phi) const
{
   using namespace TMath;

   return REveVector(Cos(phi) / CosH(eta), Sin(phi) / CosH(eta), TanH(eta));
}

////////////////////////////////////////////////////////////////////////////////
/// Returns point on the base of the cone with given eta and phi.

REveVector REveJetCone::CalcBaseVec(Float_t eta, Float_t phi) const
{
   using namespace TMath;

   REveVector vec = CalcEtaPhiVec(eta, phi);

   // -- Set length of the contourPoint
   if (fLimits.fY != 0 && fLimits.fZ != 0)
   {
      Float_t theta = vec.Theta();
      if (theta < fThetaC)
         vec *= fLimits.fZ / Cos(theta);
      else if (theta > Pi() - fThetaC)
         vec *= fLimits.fZ / Cos(theta - Pi());
      else
         vec *= fLimits.fY / Sin(theta);

      if (fLimits.fX != 0) vec *= fLimits.fX;
   }
   else
   {
      vec *= fLimits.fX;
   }

   return vec;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns point on the base of the cone with internal angle alpha:
/// alpha = 0 -> max eta,  alpha = pi/2 -> max phi, ...

REveVector REveJetCone::CalcBaseVec(Float_t alpha) const
{
   using namespace TMath;

   return CalcBaseVec(fEta + fDEta * Cos(alpha), fPhi + fDPhi * Sin(alpha));
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true if the cone is in barrel / endcap transition region.

Bool_t REveJetCone::IsInTransitionRegion() const
{
   using namespace TMath;

   Float_t tm = CalcBaseVec(0).Theta();
   Float_t tM = CalcBaseVec(Pi()).Theta();

   return (tM > fThetaC        && tm < fThetaC) ||
          (tM > Pi() - fThetaC && tm < Pi() - fThetaC);
}


/** \class REveJetConeProjected
\ingroup REve
Projection of REveJetCone.
*/

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REveJetConeProjected::REveJetConeProjected(const std::string& n, const std::string& t) :
   REveShape(n, t)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

REveJetConeProjected::~REveJetConeProjected()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Crates 3D point array for rendering.

void REveJetConeProjected::BuildRenderData()
{
   static const REveException kEH("REveJetConeProjected::BuildRenderData ");

   REveProjection *P = GetManager()->GetProjection();
   REveJetCone    *C = dynamic_cast<REveJetCone*>(GetProjectable());

   fRenderData = std::make_unique<REveRenderData>("makeJetProjected", 4);

   std::vector<REveVector> V;
   V.reserve(4);
   V.resize(3);

   switch (P->GetType())
   {
      case REveProjection::kPT_RPhi:
      {
         V[0] = C->fApex;
         V[1] = C->CalcBaseVec(TMath::Pi() + TMath::PiOver2());
         V[2] = C->CalcBaseVec(TMath::PiOver2());

         for (Int_t i = 0; i < 3; ++i)
            P->ProjectVector(V[i], fDepth);

         break;
      }

      case REveProjection::kPT_RhoZ:
      {
         V[0] = C->fApex;
         V[1] = C->CalcBaseVec(0);
         V[2] = C->CalcBaseVec(TMath::Pi());

         Float_t tm = V[1].Theta();
         Float_t tM = V[2].Theta();

         if (tM > C->fThetaC && tm < C->fThetaC)
         {
            REveVector v(0, C->fLimits.fY, C->fLimits.fZ);

            V.push_back(C->CalcBaseVec(v.Eta(), C->fPhi));
         }

         if (tM > TMath::Pi() - C->fThetaC && tm < TMath::Pi() - C->fThetaC)
         {
            REveVector v(0, C->fLimits.fY, -C->fLimits.fZ);

            V.push_back(C->CalcBaseVec(v.Eta(), C->fPhi));
         }

         for (auto &v : V) P->ProjectVector(v, fDepth);

         std::sort(V.begin() + 1, V.end(),
                   [](const auto& a, const auto &b) -> bool
                      { return a.Phi() < b.Phi(); });

         break;
      }

      default:
         throw kEH + "Unsupported projection type.";
   }

   for (auto &v : V) fRenderData->PushV(v);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute bounding-box, virtual from TAttBBox.

void REveJetConeProjected::ComputeBBox()
{
   BBoxInit();

   REveJetCone    *cone = dynamic_cast<REveJetCone*>(fProjectable);
   REveProjection *proj = GetManager()->GetProjection();
   REveVector v;
   v = cone->fApex;                                       proj->ProjectVector(v, fDepth); BBoxCheckPoint(v);
   v = cone->CalcBaseVec(0);                              proj->ProjectVector(v, fDepth); BBoxCheckPoint(v);
   v = cone->CalcBaseVec(TMath::PiOver2());               proj->ProjectVector(v, fDepth); BBoxCheckPoint(v);
   v = cone->CalcBaseVec(TMath::Pi());                    proj->ProjectVector(v, fDepth); BBoxCheckPoint(v);
   v = cone->CalcBaseVec(TMath::Pi() + TMath::PiOver2()); proj->ProjectVector(v, fDepth); BBoxCheckPoint(v);
}

////////////////////////////////////////////////////////////////////////////////
/// This is virtual method from base-class REveProjected.

void REveJetConeProjected::SetDepthLocal(Float_t d)
{
   SetDepthCommon(d, this, fBBox);
}

////////////////////////////////////////////////////////////////////////////////
/// This is virtual method from base-class REveProjected.

void REveJetConeProjected::SetProjection(REveProjectionManager* mng, REveProjectable* model)
{
   REveProjected::SetProjection(mng, model);
   CopyVizParams(dynamic_cast<REveElement*>(model));
}

////////////////////////////////////////////////////////////////////////////////
/// Re-project the jet-cone.

void REveJetConeProjected::UpdateProjection()
{
}
