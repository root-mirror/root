// @(#)root/geom:$Id$
// Author: Andrei Gheata   31/01/02

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoPara
#define ROOT_TGeoPara

#include "TGeoBBox.h"

class TGeoPara : public TGeoBBox
{
protected :
// data members
   Double_t              fX;        // X half-length
   Double_t              fY;        // Y half-length
   Double_t              fZ;        // Z half-length
   Double_t              fAlpha;    // angle w.r.t Y from the center of low Y to the high Y
   Double_t              fTheta;    // polar angle of segment between low and hi Z surfaces
   Double_t              fPhi;      // azimuthal angle of segment between low and hi Z surfaces
   Double_t              fTxy;      // tangent of XY section angle
   Double_t              fTxz;      // tangent of XZ section angle
   Double_t              fTyz;      // tangent of XZ section angle

// methods
public:
   // constructors
   TGeoPara();
   TGeoPara(Double_t dx, Double_t dy, Double_t dz, Double_t alpha, Double_t theta, Double_t phi);
   TGeoPara(const char *name, Double_t dx, Double_t dy, Double_t dz, Double_t alpha, Double_t theta, Double_t phi);
   TGeoPara(Double_t *param);
   // destructor
   ~TGeoPara() override;
   // methods
   Double_t      Capacity() const override;
   void          ComputeBBox() override;
   void          ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm) override;
   void          ComputeNormal_v(const Double_t *points, const Double_t *dirs, Double_t *norms, Int_t vecsize) override;
   Bool_t        Contains(const Double_t *point) const override;
   void          Contains_v(const Double_t *points, Bool_t *inside, Int_t vecsize) const override;
   Double_t      DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact=1,
                                   Double_t step=TGeoShape::Big(), Double_t *safe=0) const override;
   void          DistFromInside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t *step) const override;
   Double_t      DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact=1,
                                   Double_t step=TGeoShape::Big(), Double_t *safe=0) const override;
   void          DistFromOutside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t *step) const override;
   TGeoVolume   *Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv,
                                Double_t start, Double_t step) override;
   Double_t      GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const override;
   void          GetBoundingCylinder(Double_t *param) const override;
   Int_t         GetByteCount() const override {return 48;}
   Int_t         GetFittingBox(const TGeoBBox *parambox, TGeoMatrix *mat, Double_t &dx, Double_t &dy, Double_t &dz) const override;
   TGeoShape    *GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix *mat) const override;
   Int_t         GetNmeshVertices() const override {return 8;}
   Double_t              GetX() const  {return fX;}
   Double_t              GetY() const  {return fY;}
   Double_t              GetZ() const  {return fZ;}
   Double_t              GetAlpha() const {return fAlpha;}
   Double_t              GetTheta() const {return fTheta;}
   Double_t              GetPhi() const   {return fPhi;}
   Double_t              GetTxy() const {return fTxy;}
   Double_t              GetTxz() const {return fTxz;}
   Double_t              GetTyz() const {return fTyz;}
   void          InspectShape() const override;
   Bool_t        IsCylType() const override {return kFALSE;}
   Double_t      Safety(const Double_t *point, Bool_t in=kTRUE) const override;
   void          Safety_v(const Double_t *points, const Bool_t *inside, Double_t *safe, Int_t vecsize) const override;
   void          SavePrimitive(std::ostream &out, Option_t *option = "") override;
   void          SetDimensions(Double_t *param) override;
   void          SetPoints(Double_t *points) const override;
   void          SetPoints(Float_t *points) const override;
   void          Sizeof3D() const override;

   ClassDef(TGeoPara, 1)         // box primitive
};

#endif
