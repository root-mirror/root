// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEveTrackPropagator
#define ROOT_TEveTrackPropagator

#include "TEveVector.h"
#include "TEvePathMark.h"
#include "TEveUtil.h"
#include "TEveElement.h"
#include "TMarker.h"

#include <vector>

class TEvePointSet;


//==============================================================================
// TEveMagField
//==============================================================================

class TEveMagField
{
protected:
   Bool_t  fFieldConstant;

public:
   TEveMagField() : fFieldConstant(kFALSE) {}
   virtual ~TEveMagField() {}

   virtual Bool_t IsConst() const { return fFieldConstant; }

   virtual void  PrintField(Double_t x, Double_t y, Double_t z) const
   {
      TEveVector b = GetField(x, y, z);
      printf("v(%f, %f, %f) B(%f, %f, %f) \n", x, y, z, b.fX, b.fY, b.fZ);
   }

   // Track propagator uses only GetFieldD() and GetMaxFieldMagD().
   //
   // Here, in this base-class, we have to keep float versions GetField() and GetMaxFieldMag() as
   // primary sources of field data (for backward compatibility in ALICE and CMS).
   //
   // One only needs to redefine the double versions in subclasses.

   virtual Double_t    GetMaxFieldMagD() const { return GetMaxFieldMag(); }
   virtual TEveVectorD GetFieldD(Double_t x, Double_t y, Double_t z) const { return GetField(x, y, z); }

   TEveVectorD GetFieldD(const TEveVectorD &v) const { return GetFieldD(v.fX, v.fY, v.fZ); }

   virtual Float_t    GetMaxFieldMag() const { return 4; }
   virtual TEveVector GetField(Float_t, Float_t, Float_t) const { return TEveVector(); }

   ClassDef(TEveMagField, 0); // Abstract interface to magnetic field
};


//==============================================================================
// TEveMagFieldConst
//==============================================================================

class TEveMagFieldConst : public TEveMagField
{
protected:
   TEveVectorD fB;

public:
   TEveMagFieldConst(Double_t x, Double_t y, Double_t z) :
      TEveMagField(), fB(x, y, z)
   { fFieldConstant = kTRUE; }
   virtual ~TEveMagFieldConst() {}

   virtual Double_t    GetMaxFieldMagD() const { return fB.Mag(); };
   virtual TEveVectorD GetFieldD(Double_t /*x*/, Double_t /*y*/, Double_t /*z*/) const { return fB; }

   ClassDef(TEveMagFieldConst, 0); // Interface to constant magnetic field.
};


//==============================================================================
// TEveMagFieldDuo
//==============================================================================

class TEveMagFieldDuo : public TEveMagField
{
protected:
   TEveVectorD fBIn;
   TEveVectorD fBOut;
   Double_t    fR2;

public:
   TEveMagFieldDuo(Double_t r, Double_t bIn, Double_t bOut) :
      TEveMagField(),
      fBIn(0,0,bIn), fBOut(0,0,bOut), fR2(r*r)
   {
      fFieldConstant = kFALSE;
   }
   virtual ~TEveMagFieldDuo() {}

   virtual Double_t GetMaxFieldMagD() const { return std::max(fBIn.Mag(), fBOut.Mag()); }

   virtual TEveVectorD GetFieldD(Double_t x, Double_t y, Double_t /*z*/) const
   { return  ((x*x+y*y)<fR2) ? fBIn : fBOut; }

   ClassDef(TEveMagFieldDuo, 0); // Interface to magnetic field with two different values depending on radius.
};


//==============================================================================
// TEveTrackPropagator
//==============================================================================

class TEveTrackPropagator : public TEveElementList,
                            public TEveRefBackPtr
{
   friend class TEveTrackPropagatorSubEditor;

public:
   enum EStepper_e    { kHelix, kRungeKutta };

   enum EProjTrackBreaking_e { kPTB_Break, kPTB_UseFirstPointPos, kPTB_UseLastPointPos };

protected:
   struct Helix_t
   {
      Int_t    fCharge;   // Charge of tracked particle.
      Double_t fMaxAng;   // Maximum step angle.
      Double_t fMaxStep;  // Maximum allowed step size.
      Double_t fDelta;    // Maximum error in the middle of the step.

      Double_t fPhi;      // Accumulated angle to check fMaxOrbs by propagator.
      Bool_t   fValid;    // Corner case pT~0 or B~0, possible in variable mag field.

      // ----------------------------------------------------------------

      // helix parameters
      Double_t fLam;         // Momentum ratio pT/pZ.
      Double_t fR;           // Helix radius in cm.
      Double_t fPhiStep;     // Caluclated from fMinAng and fDelta.
      Double_t fSin, fCos;   // Current sin/cos(phistep).

      // Runge-Kutta parameters
      Double_t fRKStep;      // Step for Runge-Kutta.

      // cached
      TEveVectorD fB;        // Current magnetic field, cached.
      TEveVectorD fE1, fE2, fE3; // Base vectors: E1 -> B dir, E2->pT dir, E3 = E1xE2.
      TEveVectorD fPt, fPl;  // Transverse and longitudinal momentum.
      Double_t fPtMag;       // Magnitude of pT.
      Double_t fPlMag;       // Momentum parallel to mag field.
      Double_t fLStep;       // Transverse step arc-length in cm.

      // ----------------------------------------------------------------

      Helix_t();

      void UpdateCommon(const TEveVectorD & p, const TEveVectorD& b);
      void UpdateHelix (const TEveVectorD & p, const TEveVectorD& b, Bool_t full_update, Bool_t enforce_max_step);
      void UpdateRK    (const TEveVectorD & p, const TEveVectorD& b);

      void Step(const TEveVector4D& v, const TEveVectorD& p, TEveVector4D& vOut, TEveVectorD& pOut);

      Double_t GetStep()  { return fLStep * TMath::Sqrt(1 + fLam*fLam); }
      Double_t GetStep2() { return fLStep * fLStep * (1 + fLam*fLam);   }
   };


private:
   TEveTrackPropagator(const TEveTrackPropagator&);            // Not implemented
   TEveTrackPropagator& operator=(const TEveTrackPropagator&); // Not implemented

   void DistributeOffset(const TEveVectorD& off, Int_t first_point, Int_t np, TEveVectorD& p);

protected:
   EStepper_e               fStepper;

   TEveMagField*            fMagFieldObj;
   Bool_t                   fOwnMagFiledObj;

   // Track extrapolation limits
   Double_t                  fMaxR;          // Max radius for track extrapolation
   Double_t                  fMaxZ;          // Max z-coordinate for track extrapolation.
   Int_t                     fNMax;          // Max steps
   // Helix limits
   Double_t                  fMaxOrbs;       // Maximal angular path of tracks' orbits (1 ~ 2Pi).

   // Path-mark / first-vertex control
   Bool_t                   fEditPathMarks; // Show widgets for path-mark control in GUI editor.
   Bool_t                   fFitDaughters;  // Pass through daughter creation points when extrapolating a track.
   Bool_t                   fFitReferences; // Pass through given track-references when extrapolating a track.
   Bool_t                   fFitDecay;      // Pass through decay point when extrapolating a track.
   Bool_t                   fFitCluster2Ds; // Pass through 2D-clusters when extrapolating a track.
   Bool_t                   fFitLineSegments; // Pass through line when extrapolating a track.
   Bool_t                   fRnrDaughters;  // Render daughter path-marks.
   Bool_t                   fRnrReferences; // Render track-reference path-marks.
   Bool_t                   fRnrDecay;      // Render decay path-marks.
   Bool_t                   fRnrCluster2Ds; // Render 2D-clusters.
   Bool_t                   fRnrFV;         // Render first vertex.
   TMarker                  fPMAtt;         // Marker attributes for rendering of path-marks.
   TMarker                  fFVAtt;         // Marker attributes for fits vertex.

   // Handling of discontinuities in projections
   UChar_t                  fProjTrackBreaking; // Handling of projected-track breaking.
   Bool_t                   fRnrPTBMarkers;     // Render break-points on tracks.
   TMarker                  fPTBAtt;            // Marker attributes for track break-points.

   // ----------------------------------------------------------------

   // Propagation, state of current track
   std::vector<TEveVector4D> fPoints;        // Calculated point.
   std::vector<TEveVector4D> fLastPoints;    // Copy of the latest calculated points.
   TEveVectorD               fV;             // Start vertex.
   Helix_t                   fH;             // Helix.

   void    RebuildTracks();
   void    Update(const TEveVector4D& v, const TEveVectorD& p, Bool_t full_update=kFALSE, Bool_t enforce_max_step=kFALSE);
   void    Step(const TEveVector4D &v, const TEveVectorD &p, TEveVector4D &vOut, TEveVectorD &pOut);

   Bool_t  LoopToVertex(TEveVectorD& v, TEveVectorD& p);
   Bool_t  LoopToLineSegment(const TEveVectorD& s, const TEveVectorD& r, TEveVectorD& p);
   void    LoopToBounds(TEveVectorD& p);

   Bool_t  LineToVertex (TEveVectorD& v);
   void    LineToBounds (TEveVectorD& p);

   void    StepRungeKutta(Double_t step, Double_t* vect, Double_t* vout);

   Bool_t  HelixIntersectPlane(const TEveVectorD& p, const TEveVectorD& point, const TEveVectorD& normal,
                               TEveVectorD&itsect);
   Bool_t  LineIntersectPlane(const TEveVectorD& p, const TEveVectorD& point, const TEveVectorD& normal,
                              TEveVectorD& itsect);
   Bool_t  PointOverVertex(const TEveVector4D& v0, const TEveVector4D& v, Double_t* p=0);

   void    ClosestPointFromVertexToLineSegment(const TEveVectorD& v, const TEveVectorD& s, const TEveVectorD& r, Double_t rMagInv, TEveVectorD& c);
   Bool_t  ClosestPointBetweenLines(const TEveVectorD&, const TEveVectorD&, const TEveVectorD&, const TEveVectorD&, TEveVectorD& out);

public:
   TEveTrackPropagator(const char* n="TEveTrackPropagator", const char* t="",
                       TEveMagField* field=0, Bool_t own_field=kTRUE);
   virtual ~TEveTrackPropagator();

   virtual void OnZeroRefCount();

   virtual void CheckReferenceCount(const TEveException& eh="TEveElement::CheckReferenceCount ");

   virtual void ElementChanged(Bool_t update_scenes=kTRUE, Bool_t redraw=kFALSE);

   // propagation
   void   InitTrack(const TEveVectorD& v, Int_t charge);
   void   ResetTrack();

   Int_t    GetCurrentPoint() const;
   Double_t GetTrackLength(Int_t start_point=0, Int_t end_point=-1) const;

   virtual void   GoToBounds(TEveVectorD& p);
   virtual Bool_t GoToVertex(TEveVectorD& v, TEveVectorD& p);
   virtual Bool_t GoToLineSegment(const TEveVectorD& s, const TEveVectorD& r, TEveVectorD& p);

   // TEveVectorF wrappers
   void   InitTrack(const TEveVectorF& v, Int_t charge);
   void   GoToBounds(TEveVectorF& p);
   Bool_t GoToVertex(TEveVectorF& v, TEveVectorF&p);
   Bool_t GoToLineSegment(const TEveVectorF& s, const TEveVectorF& r, TEveVectorF& p);

   Bool_t IntersectPlane(const TEveVectorD& p, const TEveVectorD& point, const TEveVectorD& normal,
                         TEveVectorD& itsect);

   void   FillPointSet(TEvePointSet* ps) const;

   void   SetStepper(EStepper_e s) { fStepper = s; }

   void   SetMagField(Double_t bX, Double_t bY, Double_t bZ);
   void   SetMagField(Double_t b) { SetMagField(0, 0, b); }
   void   SetMagFieldObj(TEveMagField* field, Bool_t own_field=kTRUE);

   void   SetMaxR(Double_t x);
   void   SetMaxZ(Double_t x);
   void   SetMaxOrbs(Double_t x);
   void   SetMinAng(Double_t x);
   void   SetMaxAng(Double_t x);
   void   SetMaxStep(Double_t x);
   void   SetDelta(Double_t x);

   void   SetEditPathMarks(Bool_t x) { fEditPathMarks = x; }
   void   SetRnrDaughters(Bool_t x);
   void   SetRnrReferences(Bool_t x);
   void   SetRnrDecay(Bool_t x);
   void   SetRnrCluster2Ds(Bool_t x);
   void   SetFitDaughters(Bool_t x);
   void   SetFitReferences(Bool_t x);
   void   SetFitDecay(Bool_t x);
   void   SetFitCluster2Ds(Bool_t x);
   void   SetFitLineSegments(Bool_t x);
   void   SetRnrFV(Bool_t x);
   void   SetProjTrackBreaking(UChar_t x);
   void   SetRnrPTBMarkers(Bool_t x);

   TEveVectorD GetMagField(Double_t x, Double_t y, Double_t z) { return fMagFieldObj->GetField(x, y, z); }
   void PrintMagField(Double_t x, Double_t y, Double_t z) const;

   EStepper_e   GetStepper()  const { return fStepper;}

   Double_t GetMaxR()     const { return fMaxR;     }
   Double_t GetMaxZ()     const { return fMaxZ;     }
   Double_t GetMaxOrbs()  const { return fMaxOrbs;  }
   Double_t GetMinAng()   const;
   Double_t GetMaxAng()   const { return fH.fMaxAng;   }
   Double_t GetMaxStep()  const { return fH.fMaxStep;  }
   Double_t GetDelta()    const { return fH.fDelta;    }

   Bool_t  GetEditPathMarks() const { return fEditPathMarks; }
   Bool_t  GetRnrDaughters()  const { return fRnrDaughters;  }
   Bool_t  GetRnrReferences() const { return fRnrReferences; }
   Bool_t  GetRnrDecay()      const { return fRnrDecay;      }
   Bool_t  GetRnrCluster2Ds() const { return fRnrCluster2Ds; }
   Bool_t  GetFitDaughters()  const { return fFitDaughters;  }
   Bool_t  GetFitReferences() const { return fFitReferences; }
   Bool_t  GetFitDecay()      const { return fFitDecay;      }
   Bool_t  GetFitCluster2Ds() const { return fFitCluster2Ds; }
   Bool_t  GetFitLineSegments() const { return fFitLineSegments; }
   Bool_t  GetRnrFV()         const { return fRnrFV;         }
   UChar_t GetProjTrackBreaking() const { return fProjTrackBreaking; }
   Bool_t  GetRnrPTBMarkers()     const { return fRnrPTBMarkers; }

   TMarker& RefPMAtt()  { return fPMAtt; }
   TMarker& RefFVAtt()  { return fFVAtt; }
   TMarker& RefPTBAtt() { return fPTBAtt; }

   const std::vector<TEveVector4D>& GetLastPoints() const { return fLastPoints; }

   static Bool_t IsOutsideBounds(const TEveVectorD& point, Double_t maxRsqr, Double_t maxZ);

   static Double_t             fgDefMagField; // Default value for constant solenoid magnetic field.
   static const Double_t       fgkB2C;        // Constant for conversion of momentum to curvature.
   static TEveTrackPropagator  fgDefault;     // Default track propagator.

   static Double_t             fgEditorMaxR;  // Max R that can be set in GUI editor.
   static Double_t             fgEditorMaxZ;  // Max Z that can be set in GUI editor.

   ClassDef(TEveTrackPropagator, 0); // Calculates path of a particle taking into account special path-marks and imposed boundaries.
};

//______________________________________________________________________________
inline Bool_t TEveTrackPropagator::IsOutsideBounds(const TEveVectorD& point,
                                                   Double_t           maxRsqr,
                                                   Double_t           maxZ)
{
   // Return true if point% is outside of cylindrical bounds detrmined by
   // square radius and z.

   return TMath::Abs(point.fZ) > maxZ ||
          point.fX*point.fX + point.fY*point.fY > maxRsqr;
}

//______________________________________________________________________________
inline Bool_t TEveTrackPropagator::PointOverVertex(const TEveVector4D &v0,
                                                   const TEveVector4D &v,
                                                   Double_t           *p)
{
   static const Double_t kMinPl = 1e-5;

   TEveVectorD dv; dv.Sub(v0, v);

   Double_t dotV;

   if (TMath::Abs(fH.fPlMag) > kMinPl)
   {
      // Use longitudinal momentum to determine crossing point.
      // Works ok for spiraling helices, also for loopers.

      dotV = fH.fE1.Dot(dv);
      if (fH.fPlMag < 0)
         dotV = -dotV;
   }
   else
   {
      // Use full momentum, which is pT, under this conditions.

      dotV = fH.fE2.Dot(dv);
   }

   if (p)
      *p = dotV;

   return dotV < 0;
}

#endif
