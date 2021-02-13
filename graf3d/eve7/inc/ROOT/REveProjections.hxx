// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveProjections
#define ROOT7_REveProjections

#include <ROOT/REveVector.hxx>

#include <vector>
#include <string>

namespace ROOT {
namespace Experimental {

class REveTrans;

///////////////////////////////////////////////////////////////////////////////
/// REveProjection
/// Base for specific classes that implement non-linear projections.
///////////////////////////////////////////////////////////////////////////////

class REveProjection {
public:
   enum EPType_e { kPT_Unknown, kPT_RPhi, kPT_RhoZ, kPT_3D, kPT_End }; // projection type
   enum EPProc_e { kPP_Plane, kPP_Distort, kPP_Full };                 // projection procedure
   enum EGeoMode_e { kGM_Unknown, kGM_Polygons, kGM_Segments };        // strategy for geometry projections

   struct PreScaleEntry_t {
      Float_t fMin{0}, fMax{0};
      Float_t fOffset{0};
      Float_t fScale{1};

      PreScaleEntry_t() = default;

      PreScaleEntry_t(Float_t min, Float_t max, Float_t off, Float_t scale)
         : fMin(min), fMax(max), fOffset(off), fScale(scale)
      {
      }
   };

   typedef std::vector<PreScaleEntry_t> vPreScale_t;

protected:
   EPType_e fType;       // type
   EGeoMode_e fGeoMode;  // strategy of polygon projection (what to try first)
   std::string fName;    // name

   REveVector fCenter; // center of distortion

   bool fDisplaceOrigin; // displace point before projection

   Bool_t fUsePreScale;       // use pre-scaling
   vPreScale_t fPreScales[3]; // scaling before the distortion

   Float_t fDistortion;    // distortion
   Float_t fFixR;          // radius from which scaling remains constant
   Float_t fFixZ;          // z-coordinate from which scaling remains constant
   Float_t fPastFixRFac;   // relative scaling factor beyond fFixR as 10^x
   Float_t fPastFixZFac;   // relative scaling factor beyond fFixZ as 10^x
   Float_t fScaleR;        // scale factor to keep projected radius at fFixR fixed
   Float_t fScaleZ;        // scale factor to keep projected z-coordinate at fFixZ fixed
   Float_t fPastFixRScale; // relative scaling beyond fFixR
   Float_t fPastFixZScale; // relative scaling beyond fFixZ
   Float_t fMaxTrackStep;  // maximum distance between two points on a track

   void PreScaleVariable(Int_t dim, Float_t &v);

public:
   REveProjection();
   virtual ~REveProjection() {}

   virtual Bool_t Is2D() const = 0;
   virtual Bool_t Is3D() const = 0;

   virtual void ProjectPoint(Float_t &x, Float_t &y, Float_t &z, Float_t d, EPProc_e p = kPP_Full) = 0;

   void ProjectPointfv(Float_t *v, Float_t d);
   void ProjectPointdv(Double_t *v, Float_t d);
   void ProjectVector(REveVector &v, Float_t d);

   void ProjectPointfv(const REveTrans *t, const Float_t *p, Float_t *v, Float_t d);
   void ProjectPointdv(const REveTrans *t, const Double_t *p, Double_t *v, Float_t d);
   void ProjectVector(const REveTrans *t, REveVector &v, Float_t d);

   const char *GetName() const { return fName.c_str(); }
   void SetName(const char *txt) { fName = txt; }

   const REveVector &RefCenter() const { return fCenter; }
   virtual void SetCenter(REveVector &v) { fCenter = v; }
   virtual Float_t *GetProjectedCenter();

   void SetDisplaceOrigin(bool);
   Bool_t GetDisplaceOrigin() const { return fDisplaceOrigin; }

   void SetType(EPType_e t) { fType = t; }
   EPType_e GetType() const { return fType; }

   void SetGeoMode(EGeoMode_e m) { fGeoMode = m; }
   EGeoMode_e GetGeoMode() const { return fGeoMode; }

   Bool_t GetUsePreScale() const { return fUsePreScale; }
   void SetUsePreScale(Bool_t x) { fUsePreScale = x; }

   void PreScalePoint(Float_t &x, Float_t &y);
   void PreScalePoint(Float_t &x, Float_t &y, Float_t &z);
   void AddPreScaleEntry(Int_t coord, Float_t max_val, Float_t scale);
   void ChangePreScaleEntry(Int_t coord, Int_t entry, Float_t new_scale);
   void ClearPreScales();

   void SetDistortion(Float_t d);
   Float_t GetDistortion() const { return fDistortion; }
   Float_t GetFixR() const { return fFixR; }
   Float_t GetFixZ() const { return fFixZ; }
   void SetFixR(Float_t x);
   void SetFixZ(Float_t x);
   Float_t GetPastFixRFac() const { return fPastFixRFac; }
   Float_t GetPastFixZFac() const { return fPastFixZFac; }
   void SetPastFixRFac(Float_t x);
   void SetPastFixZFac(Float_t x);
   Float_t GetMaxTrackStep() const { return fMaxTrackStep; }
   void SetMaxTrackStep(Float_t x) { fMaxTrackStep = TMath::Max(x, 1.0f); }

   virtual Bool_t HasSeveralSubSpaces() const { return kFALSE; }
   virtual Bool_t AcceptSegment(REveVector &, REveVector &, Float_t /*tolerance*/) const { return kTRUE; }
   virtual Int_t SubSpaceId(const REveVector &) const { return 0; }
   virtual Bool_t IsOnSubSpaceBoundrary(const REveVector &) const { return kFALSE; }
   virtual void BisectBreakPoint(REveVector &vL, REveVector &vR, Float_t eps_sqr);
   virtual void BisectBreakPoint(REveVector &vL, REveVector &vR, Bool_t project_result = kFALSE, Float_t depth = 0);
   virtual void SetDirectionalVector(Int_t screenAxis, REveVector &vec);

   // utils to draw axis
   REveVector GetOrthogonalCenter(int idx, REveVector &out);
   virtual Float_t GetValForScreenPos(Int_t ax, Float_t value);
   virtual Float_t GetScreenVal(Int_t ax, Float_t value);
   Float_t GetScreenVal(Int_t i, Float_t x, REveVector &dirVec, REveVector &oCenter);
   Float_t GetLimit(Int_t i, Bool_t pos);

   static Float_t fgEps;    // resolution of projected points
   static Float_t fgEpsSqr; // square of resolution of projected points
};

//==============================================================================
// REveRhoZProjection
// Rho/Z non-linear projection.
//==============================================================================

class REveRhoZProjection : public REveProjection {
private:
   REveVector fProjectedCenter; // projected center of distortion.

public:
   REveRhoZProjection();
   virtual ~REveRhoZProjection() {}

   Bool_t Is2D() const override { return kTRUE; }
   Bool_t Is3D() const override { return kFALSE; }

   void ProjectPoint(Float_t &x, Float_t &y, Float_t &z, Float_t d, EPProc_e proc = kPP_Full) override;

   void SetCenter(REveVector &v) override;
   Float_t *GetProjectedCenter() override { return fProjectedCenter.Arr(); }

   Bool_t HasSeveralSubSpaces() const override { return kTRUE; }
   Bool_t AcceptSegment(REveVector &v1, REveVector &v2, Float_t tolerance) const override;
   Int_t SubSpaceId(const REveVector &v) const override;
   Bool_t IsOnSubSpaceBoundrary(const REveVector &v) const override;
   void SetDirectionalVector(Int_t screenAxis, REveVector &vec) override;
};

//==============================================================================
// REveRPhiProjection
// XY non-linear projection.
//==============================================================================

class REveRPhiProjection : public REveProjection {
public:
   REveRPhiProjection();
   virtual ~REveRPhiProjection() {}

   Bool_t Is2D() const override { return kTRUE; }
   Bool_t Is3D() const override { return kFALSE; }

   void ProjectPoint(Float_t &x, Float_t &y, Float_t &z, Float_t d, EPProc_e proc = kPP_Full) override;
};

//==============================================================================
// REve3DProjection
// 3D scaling "projection"
//==============================================================================

class REve3DProjection : public REveProjection {
public:
   REve3DProjection();
   virtual ~REve3DProjection() {}

   Bool_t Is2D() const override { return kFALSE; }
   Bool_t Is3D() const override { return kTRUE; }

   void ProjectPoint(Float_t &x, Float_t &y, Float_t &z, Float_t d, EPProc_e proc = kPP_Full) override;
};

} // namespace Experimental
} // namespace ROOT

#endif
