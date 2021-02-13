// @(#)root/geom:$Id$
// Author: Andrei Gheata   24/10/01

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoPcon
#define ROOT_TGeoPcon

#include "TGeoBBox.h"

class TGeoPcon : public TGeoBBox
{
protected:
   // data members
   Int_t                 fNz = 0;            // number of z planes (at least two)
   Double_t              fPhi1 = 0;          // lower phi limit (converted to [0,2*pi)
   Double_t              fDphi = 0;          // phi range
   Double_t             *fRmin = nullptr;    //[fNz] pointer to array of inner radii
   Double_t             *fRmax = nullptr;    //[fNz] pointer to array of outer radii
   Double_t             *fZ = nullptr;       //[fNz] pointer to array of Z planes positions
   Bool_t                fFullPhi = false;   //! Full phi range flag
   Double_t              fC1 = 0;            //! Cosine of phi1
   Double_t              fS1 = 0;            //! Sine of phi1
   Double_t              fC2 = 0;            //! Cosine of phi1+dphi
   Double_t              fS2 = 0;            //! Sine of phi1+dphi
   Double_t              fCm = 0;            //! Cosine of (phi1+phi2)/2
   Double_t              fSm = 0;            //! Sine of (phi1+phi2)/2
   Double_t              fCdphi = 0;         //! Cosine of dphi

   // methods
   TGeoPcon(const TGeoPcon&) = delete;
   TGeoPcon& operator=(const TGeoPcon&) = delete;

   Bool_t                HasInsideSurface() const;
   void                  SetSegsAndPolsNoInside(TBuffer3D &buff) const;

public:
   // constructors
   TGeoPcon();
   TGeoPcon(const char *name) : TGeoBBox(name, 0, 0, 0) {}
   TGeoPcon(Double_t phi, Double_t dphi, Int_t nz);
   TGeoPcon(const char *name, Double_t phi, Double_t dphi, Int_t nz);
   TGeoPcon(Double_t *params);
   // destructor
   virtual ~TGeoPcon();
   // methods
   virtual Double_t      Capacity() const;
   virtual void          ComputeBBox();
   virtual void          ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm);
   virtual void          ComputeNormal_v(const Double_t *points, const Double_t *dirs, Double_t *norms, Int_t vecsize);
   virtual Bool_t        Contains(const Double_t *point) const;
   virtual void          Contains_v(const Double_t *points, Bool_t *inside, Int_t vecsize) const;
   virtual void          DefineSection(Int_t snum, Double_t z, Double_t rmin, Double_t rmax);
   virtual Double_t      DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact=1,
                                   Double_t step=TGeoShape::Big(), Double_t *safe=0) const;
   virtual void          DistFromInside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t *step) const;
   virtual Double_t      DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact=1,
                                   Double_t step=TGeoShape::Big(), Double_t *safe=0) const;
   virtual void          DistFromOutside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t *step) const;
   Double_t              DistToSegZ(const Double_t *point, const Double_t *dir, Int_t &iz) const;
   virtual Int_t         DistancetoPrimitive(Int_t px, Int_t py);
   virtual TGeoVolume   *Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv,
                                Double_t start, Double_t step);
   virtual const char   *GetAxisName(Int_t iaxis) const;
   virtual Double_t      GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const;
   virtual void          GetBoundingCylinder(Double_t *param) const;
   virtual const TBuffer3D &GetBuffer3D(Int_t reqSections, Bool_t localFrame) const;
   virtual Int_t         GetByteCount() const {return 60+12*fNz;}
   Double_t              GetPhi1() const {return fPhi1;}
   Double_t              GetDphi() const {return fDphi;}
   Int_t                 GetNz() const   {return fNz;}
   virtual Int_t         GetNsegments() const;
   Double_t             *GetRmin() const {return fRmin;}
   Double_t              GetRmin(Int_t ipl) const;
   Double_t             *GetRmax() const {return fRmax;}
   Double_t              GetRmax(Int_t ipl) const;
   Double_t             *GetZ() const    {return fZ;}
   Double_t              GetZ(Int_t ipl) const;
   virtual TGeoShape    *GetMakeRuntimeShape(TGeoShape * /*mother*/, TGeoMatrix * /*mat*/) const {return 0;}
   virtual Int_t         GetNmeshVertices() const;
   virtual Bool_t        GetPointsOnSegments(Int_t /*npoints*/, Double_t * /*array*/) const {return kFALSE;}
   virtual void          GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const;
   virtual void          InspectShape() const;
   virtual Bool_t        IsCylType() const {return kTRUE;}
   virtual TBuffer3D    *MakeBuffer3D() const;
   Double_t             &Phi1()          {return fPhi1;}
   Double_t             &Dphi()          {return fDphi;}
   Double_t             &Rmin(Int_t ipl) {return fRmin[ipl];}
   Double_t             &Rmax(Int_t ipl) {return fRmax[ipl];}
   Double_t             &Z(Int_t ipl) {return fZ[ipl];}
   virtual Double_t      Safety(const Double_t *point, Bool_t in=kTRUE) const;
   virtual void          Safety_v(const Double_t *points, const Bool_t *inside, Double_t *safe, Int_t vecsize) const;
   Double_t              SafetyToSegment(const Double_t *point, Int_t ipl, Bool_t in=kTRUE, Double_t safmin=TGeoShape::Big()) const;
   virtual void          SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void          SetDimensions(Double_t *param);
   virtual void          SetPoints(Double_t *points) const;
   virtual void          SetPoints(Float_t *points) const;
   virtual void          SetSegsAndPols(TBuffer3D &buff) const;
   virtual void          Sizeof3D() const;

   ClassDef(TGeoPcon, 1)         // polycone class
};

#endif
