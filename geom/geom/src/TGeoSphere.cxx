// @(#)root/geom:$Id$
// Author: Andrei Gheata   31/01/02
// TGeoSphere::Contains() DistFromOutside/Out() implemented by Mihaela Gheata

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoSphere
\ingroup Geometry_classes

Spherical shell class. It takes 6 parameters :

   - spherical shell class. It takes 6 parameters :
   - inner and outer radius Rmin, Rmax
   - the theta limits Tmin, Tmax
   - the phi limits Pmin, Pmax (the sector in phi is considered
     starting from Pmin to Pmax counter-clockwise

Begin_Macro(source)
{
   TCanvas *c = new TCanvas("c", "c",0,0,600,600);
   new TGeoManager("sphere", "poza7");
   TGeoMaterial *mat = new TGeoMaterial("Al", 26.98,13,2.7);
   TGeoMedium *med = new TGeoMedium("MED",1,mat);
   TGeoVolume *top = gGeoManager->MakeBox("TOP",med,100,100,100);
   gGeoManager->SetTopVolume(top);
   TGeoVolume *vol = gGeoManager->MakeSphere("SPHERE",med, 30,40,60,120,30,240);
   vol->SetLineWidth(2);
   top->AddNode(vol,1);
   gGeoManager->CloseGeometry();
   gGeoManager->SetNsegments(30);
   top->Draw();
   TView *view = gPad->GetView();
   view->ShowAxis();
}
End_Macro
*/

#include <iostream>

#include "TGeoCone.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TVirtualGeoPainter.h"
#include "TGeoSphere.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"

ClassImp(TGeoSphere);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoSphere::TGeoSphere()
{
   SetShapeBit(TGeoShape::kGeoSph);
   fNz = 0;
   fNseg = 0;
   fRmin = 0.0;
   fRmax = 0.0;
   fTheta1 = 0.0;
   fTheta2 = 180.0;
   fPhi1 = 0.0;
   fPhi2 = 360.0;
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius

TGeoSphere::TGeoSphere(Double_t rmin, Double_t rmax, Double_t theta1,
                       Double_t theta2, Double_t phi1, Double_t phi2)
           :TGeoBBox(0, 0, 0)
{
   SetShapeBit(TGeoShape::kGeoSph);
   SetSphDimensions(rmin, rmax, theta1, theta2, phi1, phi2);
   ComputeBBox();
   SetNumberOfDivisions(20);
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius

TGeoSphere::TGeoSphere(const char *name, Double_t rmin, Double_t rmax, Double_t theta1,
                       Double_t theta2, Double_t phi1, Double_t phi2)
           :TGeoBBox(name, 0, 0, 0)
{
   SetShapeBit(TGeoShape::kGeoSph);
   SetSphDimensions(rmin, rmax, theta1, theta2, phi1, phi2);
   ComputeBBox();
   SetNumberOfDivisions(20);
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius
/// param[0] = Rmin
/// param[1] = Rmax
/// param[2] = theta1
/// param[3] = theta2
/// param[4] = phi1
/// param[5] = phi2

TGeoSphere::TGeoSphere(Double_t *param, Int_t nparam)
           :TGeoBBox(0, 0, 0)
{
   SetShapeBit(TGeoShape::kGeoSph);
   SetDimensions(param, nparam);
   ComputeBBox();
   SetNumberOfDivisions(20);
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoSphere::~TGeoSphere()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Computes capacity of the shape in [length^3]

Double_t TGeoSphere::Capacity() const
{
   Double_t th1 = fTheta1*TMath::DegToRad();
   Double_t th2 = fTheta2*TMath::DegToRad();
   Double_t ph1 = fPhi1*TMath::DegToRad();
   Double_t ph2 = fPhi2*TMath::DegToRad();
   Double_t capacity = (1./3.)*(fRmax*fRmax*fRmax-fRmin*fRmin*fRmin)*
                       TMath::Abs(TMath::Cos(th1)-TMath::Cos(th2))*
                       TMath::Abs(ph2-ph1);
   return capacity;
}

////////////////////////////////////////////////////////////////////////////////
/// compute bounding box of the sphere

void TGeoSphere::ComputeBBox()
{
   if (TGeoShape::IsSameWithinTolerance(TMath::Abs(fTheta2-fTheta1),180)) {
      if (TGeoShape::IsSameWithinTolerance(TMath::Abs(fPhi2-fPhi1),360)) {
         TGeoBBox::SetBoxDimensions(fRmax, fRmax, fRmax);
         memset(fOrigin, 0, 3*sizeof(Double_t));
         return;
      }
   }
   Double_t st1 = TMath::Sin(fTheta1*TMath::DegToRad());
   Double_t st2 = TMath::Sin(fTheta2*TMath::DegToRad());
   Double_t r1min, r1max, r2min, r2max, rmin, rmax;
   r1min = TMath::Min(fRmax*st1, fRmax*st2);
   r1max = TMath::Max(fRmax*st1, fRmax*st2);
   r2min = TMath::Min(fRmin*st1, fRmin*st2);
   r2max = TMath::Max(fRmin*st1, fRmin*st2);
   if (((fTheta1<=90) && (fTheta2>=90)) || ((fTheta2<=90) && (fTheta1>=90))) {
      r1max = fRmax;
      r2max = fRmin;
   }
   rmin = TMath::Min(r1min, r2min);
   rmax = TMath::Max(r1max, r2max);

   Double_t xc[4];
   Double_t yc[4];
   xc[0] = rmax*TMath::Cos(fPhi1*TMath::DegToRad());
   yc[0] = rmax*TMath::Sin(fPhi1*TMath::DegToRad());
   xc[1] = rmax*TMath::Cos(fPhi2*TMath::DegToRad());
   yc[1] = rmax*TMath::Sin(fPhi2*TMath::DegToRad());
   xc[2] = rmin*TMath::Cos(fPhi1*TMath::DegToRad());
   yc[2] = rmin*TMath::Sin(fPhi1*TMath::DegToRad());
   xc[3] = rmin*TMath::Cos(fPhi2*TMath::DegToRad());
   yc[3] = rmin*TMath::Sin(fPhi2*TMath::DegToRad());

   Double_t xmin = xc[TMath::LocMin(4, &xc[0])];
   Double_t xmax = xc[TMath::LocMax(4, &xc[0])];
   Double_t ymin = yc[TMath::LocMin(4, &yc[0])];
   Double_t ymax = yc[TMath::LocMax(4, &yc[0])];
   Double_t dp = fPhi2-fPhi1;
   if (dp<0) dp+=360;
   Double_t ddp = -fPhi1;
   if (ddp<0) ddp+= 360;
   if (ddp>360) ddp-=360;
   if (ddp<=dp) xmax = rmax;
   ddp = 90-fPhi1;
   if (ddp<0) ddp+= 360;
   if (ddp>360) ddp-=360;
   if (ddp<=dp) ymax = rmax;
   ddp = 180-fPhi1;
   if (ddp<0) ddp+= 360;
   if (ddp>360) ddp-=360;
   if (ddp<=dp) xmin = -rmax;
   ddp = 270-fPhi1;
   if (ddp<0) ddp+= 360;
   if (ddp>360) ddp-=360;
   if (ddp<=dp) ymin = -rmax;
   xc[0] = fRmax*TMath::Cos(fTheta1*TMath::DegToRad());
   xc[1] = fRmax*TMath::Cos(fTheta2*TMath::DegToRad());
   xc[2] = fRmin*TMath::Cos(fTheta1*TMath::DegToRad());
   xc[3] = fRmin*TMath::Cos(fTheta2*TMath::DegToRad());
   Double_t zmin = xc[TMath::LocMin(4, &xc[0])];
   Double_t zmax = xc[TMath::LocMax(4, &xc[0])];


   fOrigin[0] = (xmax+xmin)/2;
   fOrigin[1] = (ymax+ymin)/2;
   fOrigin[2] = (zmax+zmin)/2;;
   fDX = (xmax-xmin)/2;
   fDY = (ymax-ymin)/2;
   fDZ = (zmax-zmin)/2;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute normal to closest surface from POINT.

void TGeoSphere::ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm)
{
   Double_t rxy2 = point[0]*point[0]+point[1]*point[1];
   Double_t r2 = rxy2+point[2]*point[2];
   Double_t r=TMath::Sqrt(r2);
   Bool_t rzero=kFALSE;
   if (r<=1E-20) rzero=kTRUE;
   //localize theta
   Double_t phi=0;
   Double_t th=0.;
   if (!rzero) th = TMath::ACos(point[2]/r);

   //localize phi
   phi=TMath::ATan2(point[1], point[0]);

   Double_t saf[4];
   saf[0]=(TGeoShape::IsSameWithinTolerance(fRmin,0) && !TestShapeBit(kGeoThetaSeg) && !TestShapeBit(kGeoPhiSeg))?TGeoShape::Big():TMath::Abs(r-fRmin);
   saf[1]=TMath::Abs(fRmax-r);
   saf[2]=saf[3]= TGeoShape::Big();
   if (TestShapeBit(kGeoThetaSeg)) {
      if (fTheta1>0) {
         saf[2] = r*TMath::Abs(TMath::Sin(th-fTheta1*TMath::DegToRad()));
      }
      if (fTheta2<180) {
         saf[3] = r*TMath::Abs(TMath::Sin(fTheta2*TMath::DegToRad()-th));
      }
   }
   Int_t i = TMath::LocMin(4,saf);
   if (TestShapeBit(kGeoPhiSeg)) {
      Double_t c1 = TMath::Cos(fPhi1*TMath::DegToRad());
      Double_t s1 = TMath::Sin(fPhi1*TMath::DegToRad());
      Double_t c2 = TMath::Cos(fPhi2*TMath::DegToRad());
      Double_t s2 = TMath::Sin(fPhi2*TMath::DegToRad());
      if (TGeoShape::IsCloseToPhi(saf[i], point,c1,s1,c2,s2)) {
         TGeoShape::NormalPhi(point,dir,norm,c1,s1,c2,s2);
         return;
      }
   }
   if (i>1) {
      if (i==2) th=(fTheta1<90)?(fTheta1+90):(fTheta1-90);
      else      th=(fTheta2<90)?(fTheta2+90):(fTheta2-90);
      th *= TMath::DegToRad();
   }

   norm[0] = TMath::Sin(th)*TMath::Cos(phi);
   norm[1] = TMath::Sin(th)*TMath::Sin(phi);
   norm[2] = TMath::Cos(th);
   if (norm[0]*dir[0]+norm[1]*dir[1]+norm[2]*dir[2]<0) {
      norm[0] = -norm[0];
      norm[1] = -norm[1];
      norm[2] = -norm[2];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check if a point in local sphere coordinates is close to a boundary within
/// shape tolerance. Return values:
///  - 0 - not close to boundary
///  - 1 - close to Rmin boundary
///  - 2 - close to Rmax boundary
///  - 3,4 - close to phi1/phi2 boundary
///  - 5,6 - close to theta1/theta2 boundary

Int_t TGeoSphere::IsOnBoundary(const Double_t *point) const
{
   Int_t icode = 0;
   Double_t tol = TGeoShape::Tolerance();
   Double_t r2 = point[0]*point[0]+point[1]*point[1]+point[2]*point[2];
   Double_t drsqout = r2-fRmax*fRmax;
   // Test if point is on fRmax boundary
   if (TMath::Abs(drsqout)<2.*fRmax*tol) return 2;
   Double_t drsqin = r2;
   // Test if point is on fRmin boundary
   if (TestShapeBit(kGeoRSeg)) {
      drsqin -= fRmin*fRmin;
      if (TMath::Abs(drsqin)<2.*fRmin*tol) return 1;
   }
   if (TestShapeBit(kGeoPhiSeg)) {
      Double_t phi = TMath::ATan2(point[1], point[0]);
      if (phi<0) phi+=2*TMath::Pi();
      Double_t phi1 = fPhi1*TMath::DegToRad();
      Double_t phi2 = fPhi2*TMath::DegToRad();
      Double_t ddp = phi-phi1;
      if (r2*ddp*ddp < tol*tol) return 3;
      ddp = phi - phi2;
      if (r2*ddp*ddp < tol*tol) return 4;
   }
   if (TestShapeBit(kGeoThetaSeg)) {
      Double_t r = TMath::Sqrt(r2);
      Double_t theta = TMath::ACos(point[2]/r2);
      Double_t theta1 = fTheta1*TMath::DegToRad();
      Double_t theta2 = fTheta2*TMath::DegToRad();
      Double_t ddt;
      if (fTheta1>0) {
         ddt = TMath::Abs(theta-theta1);
         if (r*ddt < tol) return 5;
      }
      if (fTheta2<180) {
         ddt = TMath::Abs(theta-theta2);
         if (r*ddt < tol) return 6;
      }
   }
   return icode;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if a point is inside radius/theta/phi ranges for the spherical sector.

Bool_t TGeoSphere::IsPointInside(const Double_t *point, Bool_t checkR, Bool_t checkTh, Bool_t checkPh) const
{
   Double_t r2 = point[0]*point[0]+point[1]*point[1]+point[2]*point[2];
   if (checkR) {
      if (TestShapeBit(kGeoRSeg) && (r2<fRmin*fRmin)) return kFALSE;
      if (r2>fRmax*fRmax) return kFALSE;
   }
   if (r2<1E-20) return kTRUE;
   if (checkPh && TestShapeBit(kGeoPhiSeg)) {
      Double_t phi = TMath::ATan2(point[1], point[0]) * TMath::RadToDeg();
      while (phi < fPhi1) phi+=360.;
      Double_t dphi = fPhi2 -fPhi1;
      Double_t ddp = phi - fPhi1;
      if (ddp > dphi) return kFALSE;
   }
   if (checkTh && TestShapeBit(kGeoThetaSeg)) {
      r2=TMath::Sqrt(r2);
      // check theta range
      Double_t theta = TMath::ACos(point[2]/r2)*TMath::RadToDeg();
      if ((theta<fTheta1) || (theta>fTheta2)) return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// test if point is inside this sphere
/// check Rmin<=R<=Rmax

Bool_t TGeoSphere::Contains(const Double_t *point) const
{
   Double_t r2=point[0]*point[0]+point[1]*point[1]+point[2]*point[2];
   if (TestShapeBit(kGeoRSeg) && (r2<fRmin*fRmin)) return kFALSE;
   if (r2>fRmax*fRmax) return kFALSE;
   if (r2<1E-20) return kTRUE;
   // check phi range
   if (TestShapeBit(kGeoPhiSeg)) {
      Double_t phi = TMath::ATan2(point[1], point[0]) * TMath::RadToDeg();
      if (phi < 0 ) phi+=360.;
      Double_t dphi = fPhi2 -fPhi1;
      if (dphi < 0) dphi+=360.;
      Double_t ddp = phi - fPhi1;
      if (ddp < 0) ddp += 360.;
      if (ddp > dphi) return kFALSE;
   }
   if (TestShapeBit(kGeoThetaSeg)) {
      r2=TMath::Sqrt(r2);
      // check theta range
      Double_t theta = TMath::ACos(point[2]/r2)*TMath::RadToDeg();
      if ((theta<fTheta1) || (theta>fTheta2)) return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// compute closest distance from point px,py to each corner

Int_t TGeoSphere::DistancetoPrimitive(Int_t px, Int_t py)
{
   Int_t n = fNseg+1;
   Int_t nz = fNz+1;
   const Int_t numPoints = 2*n*nz;
   return ShapeDistancetoPrimitive(numPoints, px, py);
}

////////////////////////////////////////////////////////////////////////////////
/// compute distance from outside point to surface of the sphere
/// Check if the bounding box is crossed within the requested distance

Double_t TGeoSphere::DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
{
   Double_t sdist = TGeoBBox::DistFromOutside(point,dir, fDX, fDY, fDZ, fOrigin, step);
   if (sdist>=step) return TGeoShape::Big();
   Double_t saf[6];
   Double_t r1,r2,z1,z2,dz,si,ci;
   Double_t rxy2 = point[0]*point[0]+point[1]*point[1];
   Double_t rxy = TMath::Sqrt(rxy2);
   r2 = rxy2+point[2]*point[2];
   Double_t r=TMath::Sqrt(r2);
   Bool_t rzero=kFALSE;
   Double_t phi=0;
   if (r<1E-20) rzero=kTRUE;
   //localize theta
   Double_t th=0.;
   if (TestShapeBit(kGeoThetaSeg) && (!rzero)) {
      th = TMath::ACos(point[2]/r)*TMath::RadToDeg();
   }
   //localize phi
   if (TestShapeBit(kGeoPhiSeg)) {
      phi=TMath::ATan2(point[1], point[0])*TMath::RadToDeg();
      if (phi<0) phi+=360.;
   }
   if (iact<3 && safe) {
      saf[0]=(r<fRmin)?fRmin-r:TGeoShape::Big();
      saf[1]=(r>fRmax)?(r-fRmax):TGeoShape::Big();
      saf[2]=saf[3]=saf[4]=saf[5]= TGeoShape::Big();
      if (TestShapeBit(kGeoThetaSeg)) {
         if (th < fTheta1) {
            saf[2] = r*TMath::Sin((fTheta1-th)*TMath::DegToRad());
         }
         if (th > fTheta2) {
            saf[3] = r*TMath::Sin((th-fTheta2)*TMath::DegToRad());
         }
      }
      if (TestShapeBit(kGeoPhiSeg)) {
         Double_t dph1=phi-fPhi1;
         if (dph1<0) dph1+=360.;
         if (dph1<=90.) saf[4]=rxy*TMath::Sin(dph1*TMath::DegToRad());
         Double_t dph2=fPhi2-phi;
         if (dph2<0) dph2+=360.;
         if (dph2>90.) saf[5]=rxy*TMath::Sin(dph2*TMath::DegToRad());
      }
      *safe = saf[TMath::LocMin(6, &saf[0])];
      if (iact==0) return TGeoShape::Big();
      if (iact==1 && step<*safe) return TGeoShape::Big();
   }
   // compute distance to shape
   // first check if any crossing at all
   Double_t snxt = TGeoShape::Big();
   Double_t rdotn = point[0]*dir[0]+point[1]*dir[1]+point[2]*dir[2];
   Bool_t fullsph = (!TestShapeBit(kGeoThetaSeg) && !TestShapeBit(kGeoPhiSeg))?kTRUE:kFALSE;
   if (r>fRmax) {
      Double_t b = rdotn;
      Double_t c = r2-fRmax*fRmax;
      Double_t d=b*b-c;
      if (d<0) return TGeoShape::Big();
   }
   if (fullsph) {
      Bool_t inrmax = kFALSE;
      Bool_t inrmin = kFALSE;
      if (r<=fRmax+TGeoShape::Tolerance()) inrmax = kTRUE;
      if (r>=fRmin-TGeoShape::Tolerance()) inrmin = kTRUE;
      if (inrmax && inrmin) {
         if ((fRmax-r) < (r-fRmin)) {
         // close to Rmax
            if (rdotn>=0) return TGeoShape::Big();
            return 0.0; // already in
         }
         // close to Rmin
         if (TGeoShape::IsSameWithinTolerance(fRmin,0) || rdotn>=0) return 0.0;
         // check second crossing of Rmin
         return DistToSphere(point, dir, fRmin, kFALSE, kFALSE);
      }
   }

   // do rmin, rmax,  checking phi and theta ranges
   if (r<fRmin) {
      // check first cross of rmin
      snxt = DistToSphere(point, dir, fRmin, kTRUE);
      if (snxt<1E20) return snxt;
   } else {
      if (r>fRmax) {
         // point outside rmax, check first cross of rmax
         snxt = DistToSphere(point, dir, fRmax, kTRUE);
         if (snxt<1E20) return snxt;
         // now check second crossing of rmin
         if (fRmin>0) snxt = DistToSphere(point, dir, fRmin, kTRUE, kFALSE);
      } else {
         // point between rmin and rmax, check second cross of rmin
         if (fRmin>0) snxt = DistToSphere(point, dir, fRmin, kTRUE, kFALSE);
      }
   }
   // check theta conical surfaces
   Double_t ptnew[3];
   Double_t b,delta,xnew,ynew,znew, phi0, ddp;
   Double_t snext = snxt;
   Double_t st1=TGeoShape::Big(), st2=TGeoShape::Big();

   if (TestShapeBit(kGeoThetaSeg)) {
      if (fTheta1>0) {
         if (TGeoShape::IsSameWithinTolerance(fTheta1,90)) {
         // surface is a plane
            if (point[2]*dir[2]<0) {
               snxt = -point[2]/dir[2];
               ptnew[0] = point[0]+snxt*dir[0];
               ptnew[1] = point[1]+snxt*dir[1];
               ptnew[2] = 0;
               // check range
               if (IsPointInside(&ptnew[0], kTRUE, kFALSE, kTRUE)) return TMath::Min(snxt,snext);
            }
         } else {
            si = TMath::Sin(fTheta1*TMath::DegToRad());
            ci = TMath::Cos(fTheta1*TMath::DegToRad());
            if (ci>0) {
               r1 = fRmin*si;
               z1 = fRmin*ci;
               r2 = fRmax*si;
               z2 = fRmax*ci;
            } else {
               r1 = fRmax*si;
               z1 = fRmax*ci;
               r2 = fRmin*si;
               z2 = fRmin*ci;
            }
            dz = 0.5*(z2-z1);
            ptnew[0] = point[0];
            ptnew[1] = point[1];
            ptnew[2] = point[2]-0.5*(z1+z2);
            Double_t zinv = 1./dz;
            Double_t rin = 0.5*(r1+r2+(r2-r1)*ptnew[2]*zinv);
            // Protection in case point is outside
            Double_t sigz = TMath::Sign(1.,point[2]);
            TGeoCone::DistToCone(ptnew, dir, dz, r1, r2, b, delta);
            Bool_t skip = kFALSE;
            if (delta<0) skip = kTRUE;
            if (!skip) {
               snxt = -b-delta;
               if (sigz*ci>0 && sigz*rxy2 > sigz*rin*(rin-sigz*TGeoShape::Tolerance())) {
                  Double_t ddotn = ptnew[0]*dir[0]+ptnew[1]*dir[1]+0.5*(r1-r2)*dir[2]*zinv*TMath::Sqrt(rxy2);
                  if (sigz*ddotn>=0 || -b+delta<1.E-9) skip = kTRUE;
                  snxt = TGeoShape::Big();
               }
               if (snxt<1E10) {
                  znew = ptnew[2]+snxt*dir[2];
                  if (snxt>0 && TMath::Abs(znew)<dz) {
                     if (!TestShapeBit(kGeoPhiSeg)) st1=snxt;
                     else {
                        xnew = ptnew[0]+snxt*dir[0];
                        ynew = ptnew[1]+snxt*dir[1];
                        phi0 = TMath::ATan2(ynew, xnew) * TMath::RadToDeg();
                        ddp = phi0-fPhi1;
                        while (ddp<0) ddp+=360.;
                        if (ddp<=fPhi2-fPhi1) st1 = snxt;
                     }
                  }
               }
               if (!skip && st1>1E10) {
                  snxt = -b+delta;
                  znew = ptnew[2]+snxt*dir[2];
                  if (snxt>0 && TMath::Abs(znew)<dz) {
                     if (!TestShapeBit(kGeoPhiSeg)) st1=snxt;
                     else {
                        xnew = ptnew[0]+snxt*dir[0];
                        ynew = ptnew[1]+snxt*dir[1];
                        phi0 = TMath::ATan2(ynew, xnew) * TMath::RadToDeg();
                        ddp = phi0-fPhi1;
                        while (ddp<0) ddp+=360.;
                        if (ddp<=fPhi2-fPhi1) st1 = snxt;
                     }
                  }
               }
            }
         }
      }

      if (fTheta2<180) {
         if (TGeoShape::IsSameWithinTolerance(fTheta2,90)) {
            // surface is a plane
            if (point[2]*dir[2]<0) {
               snxt = -point[2]/dir[2];
               ptnew[0] = point[0]+snxt*dir[0];
               ptnew[1] = point[1]+snxt*dir[1];
               ptnew[2] = 0;
               // check range
               if (IsPointInside(&ptnew[0], kTRUE, kFALSE, kTRUE)) return TMath::Min(snxt,snext);
            }
         } else {
            si = TMath::Sin(fTheta2*TMath::DegToRad());
            ci = TMath::Cos(fTheta2*TMath::DegToRad());
            if (ci>0) {
               r1 = fRmin*si;
               z1 = fRmin*ci;
               r2 = fRmax*si;
               z2 = fRmax*ci;
            } else {
               r1 = fRmax*si;
               z1 = fRmax*ci;
               r2 = fRmin*si;
               z2 = fRmin*ci;
            }
            dz = 0.5*(z2-z1);
            ptnew[0] = point[0];
            ptnew[1] = point[1];
            ptnew[2] = point[2]-0.5*(z1+z2);
            Double_t zinv = 1./dz;
            Double_t rin = 0.5*(r1+r2+(r2-r1)*ptnew[2]*zinv);
            // Protection in case point is outside
            Double_t sigz = TMath::Sign(1.,point[2]);
            TGeoCone::DistToCone(ptnew, dir, dz, r1, r2, b, delta);
            Bool_t skip = kFALSE;
            if (delta<0) skip = kTRUE;
            if (!skip) {
               snxt = -b-delta;
               if (sigz*ci>0 && sigz*rxy2 < sigz*rin*(rin+sigz*TGeoShape::Tolerance())) {
                  Double_t ddotn = ptnew[0]*dir[0]+ptnew[1]*dir[1]+0.5*(r1-r2)*dir[2]*zinv*TMath::Sqrt(rxy2);
                  if (sigz*ddotn<=0 || -b+delta<1.E-9) skip = kTRUE;
                  snxt = TGeoShape::Big();
               }
               if (snxt<1E10) {
                  znew = ptnew[2]+snxt*dir[2];
                  if (snxt>0 && TMath::Abs(znew)<dz) {
                     if (!TestShapeBit(kGeoPhiSeg)) st2=snxt;
                     else {
                        xnew = ptnew[0]+snxt*dir[0];
                        ynew = ptnew[1]+snxt*dir[1];
                        phi0 = TMath::ATan2(ynew, xnew) * TMath::RadToDeg();
                        ddp = phi0-fPhi1;
                        while (ddp<0) ddp+=360.;
                        if (ddp<=fPhi2-fPhi1) st2 = snxt;
                     }
                  }
               }
               if (!skip && st2>1E10) {
                  snxt = -b+delta;
                  znew = ptnew[2]+snxt*dir[2];
                  if (snxt>0 && TMath::Abs(znew)<dz) {
                     if (!TestShapeBit(kGeoPhiSeg)) st2=snxt;
                     else {
                        xnew = ptnew[0]+snxt*dir[0];
                        ynew = ptnew[1]+snxt*dir[1];
                        phi0 = TMath::ATan2(ynew, xnew) * TMath::RadToDeg();
                        ddp = phi0-fPhi1;
                        while (ddp<0) ddp+=360.;
                        if (ddp<=fPhi2-fPhi1) st2 = snxt;
                     }
                  }
               }
            }
         }
      }
   }
   snxt = TMath::Min(st1, st2);
   snxt = TMath::Min(snxt,snext);
//   if (snxt<1E20) return snxt;
   if (TestShapeBit(kGeoPhiSeg)) {
      Double_t s1 = TMath::Sin(fPhi1*TMath::DegToRad());
      Double_t c1 = TMath::Cos(fPhi1*TMath::DegToRad());
      Double_t s2 = TMath::Sin(fPhi2*TMath::DegToRad());
      Double_t c2 = TMath::Cos(fPhi2*TMath::DegToRad());
      Double_t phim = 0.5*(fPhi1+fPhi2);
      Double_t sm = TMath::Sin(phim*TMath::DegToRad());
      Double_t cm = TMath::Cos(phim*TMath::DegToRad());
      Double_t sfi1=TGeoShape::Big();
      Double_t sfi2=TGeoShape::Big();
      Double_t s=0;
      Double_t safety, un;
      safety = point[0]*s1-point[1]*c1;
      if (safety>0) {
         un = dir[0]*s1-dir[1]*c1;
         if (un<0) {
            s=-safety/un;
            ptnew[0] = point[0]+s*dir[0];
            ptnew[1] = point[1]+s*dir[1];
            ptnew[2] = point[2]+s*dir[2];
            if ((ptnew[1]*cm-ptnew[0]*sm)<=0) {
               sfi1=s;
               if (IsPointInside(&ptnew[0], kTRUE, kTRUE, kFALSE) && sfi1<snxt) return sfi1;
            }
         }
      }
      safety = -point[0]*s2+point[1]*c2;
      if (safety>0) {
         un = -dir[0]*s2+dir[1]*c2;
         if (un<0) {
            s=-safety/un;
            ptnew[0] = point[0]+s*dir[0];
            ptnew[1] = point[1]+s*dir[1];
            ptnew[2] = point[2]+s*dir[2];
            if ((ptnew[1]*cm-ptnew[0]*sm)>=0) {
               sfi2=s;
               if (IsPointInside(&ptnew[0], kTRUE, kTRUE, kFALSE) && sfi2<snxt) return sfi2;
            }
         }
      }
   }
   return snxt;
}

////////////////////////////////////////////////////////////////////////////////
/// compute distance from inside point to surface of the sphere

Double_t TGeoSphere::DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
{
   Double_t saf[6];
   Double_t rxy2 = point[0]*point[0]+point[1]*point[1];
   Double_t rxy = TMath::Sqrt(rxy2);
   Double_t rad2 = rxy2+point[2]*point[2];
   Double_t r=TMath::Sqrt(rad2);
   Bool_t rzero=kFALSE;
   if (r<=1E-20) rzero=kTRUE;
   //localize theta
   Double_t phi=0;;
   Double_t th=0.;
   if (TestShapeBit(kGeoThetaSeg) && (!rzero)) {
      th = TMath::ACos(point[2]/r)*TMath::RadToDeg();
   }
   //localize phi
   if (TestShapeBit(kGeoPhiSeg)) {
      phi=TMath::ATan2(point[1], point[0])*TMath::RadToDeg();
      if (phi<0) phi+=360.;
   }
   if (iact<3 && safe) {
      saf[0]=(TGeoShape::IsSameWithinTolerance(fRmin,0))?TGeoShape::Big():r-fRmin;
      saf[1]=fRmax-r;
      saf[2]=saf[3]=saf[4]=saf[5]= TGeoShape::Big();
      if (TestShapeBit(kGeoThetaSeg)) {
         if (fTheta1>0) {
            saf[2] = r*TMath::Sin((th-fTheta1)*TMath::DegToRad());
         }
         if (fTheta2<180) {
            saf[3] = r*TMath::Sin((fTheta2-th)*TMath::DegToRad());
         }
      }
      if (TestShapeBit(kGeoPhiSeg)) {
         Double_t dph1=phi-fPhi1;
         if (dph1<0) dph1+=360.;
         if (dph1<=90.) saf[4]=rxy*TMath::Sin(dph1*TMath::DegToRad());
         Double_t dph2=fPhi2-phi;
         if (dph2<0) dph2+=360.;
         if (dph2<=90.) saf[5]=rxy*TMath::Sin(dph2*TMath::DegToRad());
      }
      *safe = saf[TMath::LocMin(6, &saf[0])];
      if (iact==0) return TGeoShape::Big();
      if (iact==1 && step<*safe) return TGeoShape::Big();
   }
   // compute distance to shape
   Double_t snxt = TGeoShape::Big();
   if (rzero) {
//      gGeoManager->SetNormalChecked(1.);
      return fRmax;
   }
   // first do rmin, rmax
   Double_t b,delta, xnew,ynew,znew, phi0, ddp;
   Double_t rdotn = point[0]*dir[0]+point[1]*dir[1]+point[2]*dir[2];
   Double_t sn1 = TGeoShape::Big();
   // Inner sphere
   if (fRmin>0) {
      // Protection in case point is actually outside the sphere
      if (r <= fRmin+TGeoShape::Tolerance()) {
         if (rdotn<0) return 0.0;
      } else {
         if (rdotn<0) sn1 = DistToSphere(point, dir, fRmin, kFALSE);
      }
   }
   Double_t sn2 = TGeoShape::Big();
   // Outer sphere
   if (r >= fRmax-TGeoShape::Tolerance()) {
      if (rdotn>=0) return 0.0;
   }
   sn2 = DistToSphere(point, dir, fRmax, kFALSE, kFALSE);
   Double_t sr = TMath::Min(sn1, sn2);
   // check theta conical surfaces
   sn1 = TGeoShape::Big();
   sn2 = TGeoShape::Big();
   if (TestShapeBit(kGeoThetaSeg)) {
      if (TGeoShape::IsSameWithinTolerance(fTheta1,90)) {
      // surface is a plane
         if (point[2]*dir[2]<0)  sn1 = -point[2]/dir[2];
      } else {
         if (fTheta1>0) {
            Double_t r1,r2,z1,z2,dz,ptnew[3];
            Double_t si = TMath::Sin(fTheta1*TMath::DegToRad());
            Double_t ci = TMath::Cos(fTheta1*TMath::DegToRad());
            if (ci>0) {
               r1 = fRmin*si;
               z1 = fRmin*ci;
               r2 = fRmax*si;
               z2 = fRmax*ci;
            } else {
               r1 = fRmax*si;
               z1 = fRmax*ci;
               r2 = fRmin*si;
               z2 = fRmin*ci;
            }
            dz = 0.5*(z2-z1);
            ptnew[0] = point[0];
            ptnew[1] = point[1];
            ptnew[2] = point[2]-0.5*(z1+z2);
            Double_t zinv = 1./dz;
            Double_t rin = 0.5*(r1+r2+(r2-r1)*ptnew[2]*zinv);
            // Protection in case point is outside
            Double_t sigz = TMath::Sign(1.,point[2]);
            if (sigz*ci>0 && sigz*rxy2 < sigz*rin*(rin+sigz*TGeoShape::Tolerance())) {
               Double_t ddotn = ptnew[0]*dir[0]+ptnew[1]*dir[1]+0.5*(r1-r2)*dir[2]*zinv*TMath::Sqrt(rxy2);
               if (sigz*ddotn<=0) return 0.0;
            } else {
               TGeoCone::DistToCone(ptnew, dir, dz, r1, r2, b, delta);
               if (delta>0) {
                  snxt = -b-delta;
                  znew = ptnew[2]+snxt*dir[2];
                  if (snxt>0 && TMath::Abs(znew)<dz) {
                     if (!TestShapeBit(kGeoPhiSeg)) sn1=snxt;
                     else {
                        xnew = ptnew[0]+snxt*dir[0];
                        ynew = ptnew[1]+snxt*dir[1];
                        phi0 = TMath::ATan2(ynew, xnew) * TMath::RadToDeg();
                        ddp = phi0-fPhi1;
                        while (ddp<0) ddp+=360.;
                        if (ddp<=fPhi2-fPhi1) sn1 = snxt;
                     }
                  }
                  if (sn1>1E10) {
                     snxt = -b+delta;
                     znew = ptnew[2]+snxt*dir[2];
                     if (snxt>0 && TMath::Abs(znew)<dz) {
                        if (!TestShapeBit(kGeoPhiSeg)) sn1=snxt;
                        else {
                           xnew = ptnew[0]+snxt*dir[0];
                           ynew = ptnew[1]+snxt*dir[1];
                           phi0 = TMath::ATan2(ynew, xnew) * TMath::RadToDeg();
                           ddp = phi0-fPhi1;
                           while (ddp<0) ddp+=360.;
                           if (ddp<=fPhi2-fPhi1) sn1 = snxt;
                        }
                     }
                  }
               }
            }
         }
      }
      if (TGeoShape::IsSameWithinTolerance(fTheta2,90)) {
         // surface is a plane
         if (point[2]*dir[2]<0)  sn1 = -point[2]/dir[2];
      } else {
         if (fTheta2<180) {
            Double_t r1,r2,z1,z2,dz,ptnew[3];
            Double_t si = TMath::Sin(fTheta2*TMath::DegToRad());
            Double_t ci = TMath::Cos(fTheta2*TMath::DegToRad());
            if (ci>0) {
               r1 = fRmin*si;
               z1 = fRmin*ci;
               r2 = fRmax*si;
               z2 = fRmax*ci;
            } else {
               r1 = fRmax*si;
               z1 = fRmax*ci;
               r2 = fRmin*si;
               z2 = fRmin*ci;
            }
            dz = 0.5*(z2-z1);
            ptnew[0] = point[0];
            ptnew[1] = point[1];
            ptnew[2] = point[2]-0.5*(z1+z2);
            Double_t zinv = 1./dz;
            Double_t rin = 0.5*(r1+r2+(r2-r1)*ptnew[2]*zinv);
            // Protection in case point is outside
            Double_t sigz = TMath::Sign(1.,point[2]);
            if (sigz*ci>0 && sigz*rxy2 > sigz*rin*(rin-sigz*TGeoShape::Tolerance())) {
               Double_t ddotn = ptnew[0]*dir[0]+ptnew[1]*dir[1]+0.5*(r1-r2)*dir[2]*zinv*TMath::Sqrt(rxy2);
               if (sigz*ddotn>=0) return 0.0;
            } else {
               TGeoCone::DistToCone(ptnew, dir, dz, r1, r2, b, delta);
               if (delta>0) {
                  snxt = -b-delta;
                  znew = ptnew[2]+snxt*dir[2];
                  if (snxt>0 && TMath::Abs(znew)<dz) {
                     if (!TestShapeBit(kGeoPhiSeg)) sn2=snxt;
                     else {
                        xnew = ptnew[0]+snxt*dir[0];
                        ynew = ptnew[1]+snxt*dir[1];
                        phi0 = TMath::ATan2(ynew, xnew) * TMath::RadToDeg();
                        ddp = phi0-fPhi1;
                        while (ddp<0) ddp+=360.;
                        if (ddp<=fPhi2-fPhi1) sn2 = snxt;
                     }
                  }
                  if (sn2>1E10) {
                     snxt = -b+delta;
                     znew = ptnew[2]+snxt*dir[2];
                     if (snxt>0 && TMath::Abs(znew)<dz) {
                        if (!TestShapeBit(kGeoPhiSeg)) sn2=snxt;
                        else {
                           xnew = ptnew[0]+snxt*dir[0];
                           ynew = ptnew[1]+snxt*dir[1];
                           phi0 = TMath::ATan2(ynew, xnew) * TMath::RadToDeg();
                           ddp = phi0-fPhi1;
                           while (ddp<0) ddp+=360.;
                           if (ddp<=fPhi2-fPhi1) sn2 = snxt;
                        }
                     }
                  }
               }
            }
         }
      }
   }
   Double_t st = TMath::Min(sn1,sn2);
   Double_t sp = TGeoShape::Big();
   if (TestShapeBit(kGeoPhiSeg)) {
      Double_t s1 = TMath::Sin(fPhi1*TMath::DegToRad());
      Double_t c1 = TMath::Cos(fPhi1*TMath::DegToRad());
      Double_t s2 = TMath::Sin(fPhi2*TMath::DegToRad());
      Double_t c2 = TMath::Cos(fPhi2*TMath::DegToRad());
      Double_t phim = 0.5*(fPhi1+fPhi2);
      Double_t sm = TMath::Sin(phim*TMath::DegToRad());
      Double_t cm = TMath::Cos(phim*TMath::DegToRad());
      sp = TGeoShape::DistToPhiMin(point, dir, s1, c1, s2, c2, sm, cm);
   }
   snxt = TMath::Min(sr, st);
   snxt = TMath::Min(snxt, sp);
   return snxt;
}

////////////////////////////////////////////////////////////////////////////////
/// compute distance to sphere of radius rsph. Direction has to be a unit vector

Double_t TGeoSphere::DistToSphere(const Double_t *point, const Double_t *dir, Double_t rsph, Bool_t check, Bool_t firstcross) const
{
   if (rsph<=0) return TGeoShape::Big();
   Double_t s=TGeoShape::Big();
   Double_t r2 = point[0]*point[0]+point[1]*point[1]+point[2]*point[2];
   Double_t b = point[0]*dir[0]+point[1]*dir[1]+point[2]*dir[2];
   Double_t c = r2-rsph*rsph;
   Bool_t in = (c<=0)?kTRUE:kFALSE;
   Double_t d;

   d=b*b-c;
   if (d<0) return TGeoShape::Big();
   Double_t pt[3];
   Int_t i;
   d = TMath::Sqrt(d);
   if (in) {
      s=-b+d;
   } else {
      s = (firstcross)?(-b-d):(-b+d);
   }
   if (s<0) return TGeoShape::Big();
   if (!check) return s;
   for (i=0; i<3; i++) pt[i]=point[i]+s*dir[i];
   // check theta and phi ranges
   if (IsPointInside(&pt[0], kFALSE)) return s;
   return TGeoShape::Big();
}

////////////////////////////////////////////////////////////////////////////////

TGeoVolume *TGeoSphere::Divide(TGeoVolume * voldiv, const char * divname, Int_t iaxis, Int_t ndiv,
                               Double_t start, Double_t step)
{
   TGeoShape *shape;           //--- shape to be created
   TGeoVolume *vol;            //--- division volume to be created
   TGeoVolumeMulti *vmulti;    //--- generic divided volume
   TGeoPatternFinder *finder;  //--- finder to be attached
   TString opt = "";           //--- option to be attached
   Int_t id;
   Double_t end = start+ndiv*step;
   switch (iaxis) {
      case 1:  //---                R division
         finder = new TGeoPatternSphR(voldiv, ndiv, start, end);
         vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
         voldiv->SetFinder(finder);
         finder->SetDivIndex(voldiv->GetNdaughters());
         for (id=0; id<ndiv; id++) {
            shape = new TGeoSphere(start+id*step, start+(id+1)*step, fTheta1,fTheta2,fPhi1,fPhi2);
            vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
            vmulti->AddVolume(vol);
            opt = "R";
            voldiv->AddNodeOffset(vol, id, 0, opt.Data());
            ((TGeoNodeOffset*)voldiv->GetNodes()->At(voldiv->GetNdaughters()-1))->SetFinder(finder);
         }
         return vmulti;
      case 2:  //---                Phi division
         finder = new TGeoPatternSphPhi(voldiv, ndiv, start, end);
         voldiv->SetFinder(finder);
         finder->SetDivIndex(voldiv->GetNdaughters());
         shape = new TGeoSphere(fRmin, fRmax, fTheta1, fTheta2, -step/2, step/2);
         vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
         vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
         vmulti->AddVolume(vol);
         opt = "Phi";
         for (id=0; id<ndiv; id++) {
            voldiv->AddNodeOffset(vol, id, start+id*step+step/2, opt.Data());
            ((TGeoNodeOffset*)voldiv->GetNodes()->At(voldiv->GetNdaughters()-1))->SetFinder(finder);
         }
         return vmulti;
      case 3:  //---                Theta division
         finder = new TGeoPatternSphTheta(voldiv, ndiv, start, end);
         vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
         voldiv->SetFinder(finder);
         finder->SetDivIndex(voldiv->GetNdaughters());
         for (id=0; id<ndiv; id++) {
            shape = new TGeoSphere(fRmin,fRmax,start+id*step, start+(id+1)*step,fPhi1,fPhi2);
            vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
            vmulti->AddVolume(vol);
            opt = "Theta";
            voldiv->AddNodeOffset(vol, id, 0, opt.Data());
            ((TGeoNodeOffset*)voldiv->GetNodes()->At(voldiv->GetNdaughters()-1))->SetFinder(finder);
         }
         return vmulti;
      default:
         Error("Divide", "In shape %s wrong axis type for division", GetName());
         return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns name of axis IAXIS.

const char *TGeoSphere::GetAxisName(Int_t iaxis) const
{
   switch (iaxis) {
      case 1:
         return "R";
      case 2:
         return "PHI";
      case 3:
         return "THETA";
      default:
         return "UNDEFINED";
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get range of shape for a given axis.

Double_t TGeoSphere::GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const
{
   xlo = 0;
   xhi = 0;
   Double_t dx = 0;
   switch (iaxis) {
      case 1:
         xlo = fRmin;
         xhi = fRmax;
         dx = xhi-xlo;
         return dx;
      case 2:
         xlo = fPhi1;
         xhi = fPhi2;
         dx = xhi-xlo;
         return dx;
      case 3:
         xlo = fTheta1;
         xhi = fTheta2;
         dx = xhi-xlo;
         return dx;
   }
   return dx;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill vector param[4] with the bounding cylinder parameters. The order
/// is the following : Rmin, Rmax, Phi1, Phi2

void TGeoSphere::GetBoundingCylinder(Double_t *param) const
{
   Double_t smin = TMath::Sin(fTheta1*TMath::DegToRad());
   Double_t smax = TMath::Sin(fTheta2*TMath::DegToRad());
   if (smin>smax) {
      Double_t a = smin;
      smin = smax;
      smax = a;
   }
   param[0] = fRmin*smin; // Rmin
   param[0] *= param[0];
   if (((90.-fTheta1)*(fTheta2-90.))>=0) smax = 1.;
   param[1] = fRmax*smax; // Rmax
   param[1] *= param[1];
   param[2] = (fPhi1<0)?(fPhi1+360.):fPhi1; // Phi1
   param[3] = fPhi2;
   if (TGeoShape::IsSameWithinTolerance(param[3]-param[2],360)) {         // Phi2
      param[2] = 0.;
      param[3] = 360.;
   }
   while (param[3]<param[2]) param[3]+=360.;
}

////////////////////////////////////////////////////////////////////////////////
/// print shape parameters

void TGeoSphere::InspectShape() const
{
   printf("*** Shape %s: TGeoSphere ***\n", GetName());
   printf("    Rmin = %11.5f\n", fRmin);
   printf("    Rmax = %11.5f\n", fRmax);
   printf("    Th1  = %11.5f\n", fTheta1);
   printf("    Th2  = %11.5f\n", fTheta2);
   printf("    Ph1  = %11.5f\n", fPhi1);
   printf("    Ph2  = %11.5f\n", fPhi2);
   printf(" Bounding box:\n");
   TGeoBBox::InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a TBuffer3D describing *this* shape.
/// Coordinates are in local reference frame.

TBuffer3D *TGeoSphere::MakeBuffer3D() const
{
   Bool_t full = kTRUE;
   if (TestShapeBit(kGeoThetaSeg) || TestShapeBit(kGeoPhiSeg)) full = kFALSE;
   Int_t ncenter = 1;
   if (full || TestShapeBit(kGeoRSeg)) ncenter = 0;
   Int_t nup = (fTheta1>0)?0:1;
   Int_t ndown = (fTheta2<180)?0:1;
   // number of different latitudes, excluding 0 and 180 degrees
   Int_t nlat = fNz+1-(nup+ndown);
   // number of different longitudes
   Int_t nlong = fNseg;
   if (TestShapeBit(kGeoPhiSeg)) nlong++;

   Int_t nbPnts = nlat*nlong+nup+ndown+ncenter;
   if (TestShapeBit(kGeoRSeg)) nbPnts *= 2;

   Int_t nbSegs = nlat*fNseg + (nlat-1+nup+ndown)*nlong; // outer sphere
   if (TestShapeBit(kGeoRSeg)) nbSegs *= 2; // inner sphere
   if (TestShapeBit(kGeoPhiSeg)) nbSegs += 2*nlat+nup+ndown; // 2 phi planes
   nbSegs += nlong * (2-nup - ndown);  // connecting cones

   Int_t nbPols = fNz*fNseg; // outer
   if (TestShapeBit(kGeoRSeg)) nbPols *=2;  // inner
   if (TestShapeBit(kGeoPhiSeg)) nbPols += 2*fNz; // 2 phi planes
   nbPols += (2-nup-ndown)*fNseg; // connecting

   TBuffer3D* buff = new TBuffer3D(TBuffer3DTypes::kGeneric,
                                   nbPnts, 3*nbPnts, nbSegs, 3*nbSegs, nbPols, 6*nbPols);

   if (buff)
   {
      SetPoints(buff->fPnts);
      SetSegsAndPols(*buff);
   }

   return buff;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill TBuffer3D structure for segments and polygons.

void TGeoSphere::SetSegsAndPols(TBuffer3D & buff) const
{
   Bool_t full = kTRUE;
   if (TestShapeBit(kGeoThetaSeg) || TestShapeBit(kGeoPhiSeg)) full = kFALSE;
   Int_t ncenter = 1;
   if (full || TestShapeBit(kGeoRSeg)) ncenter = 0;
   Int_t nup = (fTheta1>0)?0:1;
   Int_t ndown = (fTheta2<180)?0:1;
   // number of different latitudes, excluding 0 and 180 degrees
   Int_t nlat = fNz+1-(nup+ndown);
   // number of different longitudes
   Int_t nlong = fNseg;
   if (TestShapeBit(kGeoPhiSeg)) nlong++;

   Int_t nbPnts = nlat*nlong+nup+ndown+ncenter;
   if (TestShapeBit(kGeoRSeg)) nbPnts *= 2;

   Int_t nbSegs = nlat*fNseg + (nlat-1+nup+ndown)*nlong; // outer sphere
   if (TestShapeBit(kGeoRSeg)) nbSegs *= 2; // inner sphere
   if (TestShapeBit(kGeoPhiSeg)) nbSegs += 2*nlat+nup+ndown; // 2 phi planes
   nbSegs += nlong * (2-nup - ndown);  // connecting cones

   Int_t nbPols = fNz*fNseg; // outer
   if (TestShapeBit(kGeoRSeg)) nbPols *=2;  // inner
   if (TestShapeBit(kGeoPhiSeg)) nbPols += 2*fNz; // 2 phi planes
   nbPols += (2-nup-ndown)*fNseg; // connecting

   Int_t c = GetBasicColor();
   Int_t i, j;
   Int_t indx;
   indx = 0;
   // outside sphere
   // loop all segments on latitudes (except 0 and 180 degrees)
   // [0, nlat*fNseg)
   Int_t indpar = 0;
   for (i=0; i<nlat; i++) {
      for (j=0; j<fNseg; j++) {
         buff.fSegs[indx++]   = c;
         buff.fSegs[indx++] = i*nlong+j;
         buff.fSegs[indx++] = i*nlong+(j+1)%nlong;
      }
   }
   // loop all segments on longitudes
   // nlat*fNseg + [0, (nlat-1)*nlong)
   Int_t indlong = indpar + nlat*fNseg;
   for (i=0; i<nlat-1; i++) {
      for (j=0; j<nlong; j++) {
         buff.fSegs[indx++]   = c;
         buff.fSegs[indx++] = i*nlong+j;
         buff.fSegs[indx++] = (i+1)*nlong+j;
      }
   }
   Int_t indup = indlong + (nlat-1)*nlong;
   // extra longitudes on top
   // nlat*fNseg+(nlat-1)*nlong + [0, nlong)
   if (nup) {
      Int_t indpup = nlat*nlong;
      for (j=0; j<nlong; j++) {
         buff.fSegs[indx++]   = c;
         buff.fSegs[indx++] = j;
         buff.fSegs[indx++] = indpup;
      }
   }
   Int_t inddown = indup + nup*nlong;
   // extra longitudes on bottom
   // nlat*fNseg+(nlat+nup-1)*nlong + [0, nlong)
   if (ndown) {
      Int_t indpdown = nlat*nlong+nup;
      for (j=0; j<nlong; j++) {
         buff.fSegs[indx++]   = c;
         buff.fSegs[indx++] = (nlat-1)*nlong+j;
         buff.fSegs[indx++] = indpdown;
      }
   }
   Int_t indparin = inddown + ndown*nlong;
   Int_t indlongin = indparin;
   Int_t indupin = indparin;
   Int_t inddownin = indparin;
   Int_t indphi = indparin;
   // inner sphere
   Int_t indptin = nlat*nlong + nup + ndown;
   Int_t iptcenter = indptin;
   // nlat*fNseg+(nlat+nup+ndown-1)*nlong
   if (TestShapeBit(kGeoRSeg)) {
      indlongin = indparin + nlat*fNseg;
      indupin   = indlongin + (nlat-1)*nlong;
      inddownin = indupin + nup*nlong;
      // loop all segments on latitudes (except 0 and 180 degrees)
      // indsegin + [0, nlat*fNseg)
      for (i=0; i<nlat; i++) {
         for (j=0; j<fNseg; j++) {
            buff.fSegs[indx++]   = c+1;
            buff.fSegs[indx++] = indptin + i*nlong+j;
            buff.fSegs[indx++] = indptin + i*nlong+(j+1)%nlong;
         }
      }
      // loop all segments on longitudes
      // indsegin + nlat*fNseg + [0, (nlat-1)*nlong)
      for (i=0; i<nlat-1; i++) {
         for (j=0; j<nlong; j++) {
            buff.fSegs[indx++]   = c+1;
            buff.fSegs[indx++] = indptin + i*nlong+j;
            buff.fSegs[indx++] = indptin + (i+1)*nlong+j;
         }
      }
      // extra longitudes on top
      // indsegin + nlat*fNseg+(nlat-1)*nlong + [0, nlong)
      if (nup) {
         Int_t indupltop = indptin + nlat*nlong;
         for (j=0; j<nlong; j++) {
            buff.fSegs[indx++]   = c+1;
            buff.fSegs[indx++] = indptin + j;
            buff.fSegs[indx++] = indupltop;
         }
      }
      // extra longitudes on bottom
      // indsegin + nlat*fNseg+(nlat+nup-1)*nlong + [0, nlong)
      if (ndown) {
         Int_t indpdown = indptin + nlat*nlong+nup;
         for (j=0; j<nlong; j++) {
            buff.fSegs[indx++]   = c+1;
            buff.fSegs[indx++] = indptin + (nlat-1)*nlong+j;
            buff.fSegs[indx++] = indpdown;
         }
      }
      indphi = inddownin + ndown*nlong;
   }
   Int_t indtheta = indphi;
   // Segments on phi planes
   if (TestShapeBit(kGeoPhiSeg)) {
      indtheta += 2*nlat + nup + ndown;
      for (j=0; j<nlat; j++) {
         buff.fSegs[indx++]   = c+2;
         buff.fSegs[indx++] = j*nlong;
         if (TestShapeBit(kGeoRSeg)) buff.fSegs[indx++] = indptin + j*nlong;
         else buff.fSegs[indx++] = iptcenter;
      }
      for (j=0; j<nlat; j++) {
         buff.fSegs[indx++]   = c+2;
         buff.fSegs[indx++] = (j+1)*nlong-1;
         if (TestShapeBit(kGeoRSeg)) buff.fSegs[indx++] = indptin + (j+1)*nlong-1;
         else buff.fSegs[indx++] = iptcenter;
      }
      if (nup) {
         buff.fSegs[indx++]   = c+2;
         buff.fSegs[indx++] = nlat*nlong;
         if (TestShapeBit(kGeoRSeg)) buff.fSegs[indx++] = indptin + nlat*nlong;
         else buff.fSegs[indx++] = iptcenter;
      }
      if (ndown) {
         buff.fSegs[indx++]   = c+2;
         buff.fSegs[indx++] = nlat*nlong+nup;
         if (TestShapeBit(kGeoRSeg)) buff.fSegs[indx++] = indptin + nlat*nlong+nup;
         else buff.fSegs[indx++] = iptcenter;
      }
   }
   // Segments on cones
   if (!nup) {
      for (j=0; j<nlong; j++) {
         buff.fSegs[indx++]   = c+2;
         buff.fSegs[indx++] = j;
         if (TestShapeBit(kGeoRSeg)) buff.fSegs[indx++] = indptin + j;
         else buff.fSegs[indx++] = iptcenter;
      }
   }
   if (!ndown) {
      for (j=0; j<nlong; j++) {
         buff.fSegs[indx++]   = c+2;
         buff.fSegs[indx++] = (nlat-1)*nlong + j;
         if (TestShapeBit(kGeoRSeg)) buff.fSegs[indx++] = indptin + (nlat-1)*nlong +j;
         else buff.fSegs[indx++] = iptcenter;
      }
   }

   indx = 0;
   // Fill polygons for outside sphere (except 0/180)
   for (i=0; i<nlat-1; i++) {
      for (j=0; j<fNseg; j++) {
         buff.fPols[indx++] = c;
         buff.fPols[indx++] = 4;
         buff.fPols[indx++] = indpar+i*fNseg+j;
         buff.fPols[indx++] = indlong+i*nlong+(j+1)%nlong;
         buff.fPols[indx++] = indpar+(i+1)*fNseg+j;
         buff.fPols[indx++] = indlong+i*nlong+j;
      }
   }
   // upper
   if (nup) {
      for (j=0; j<fNseg; j++) {
         buff.fPols[indx++] = c;
         buff.fPols[indx++] = 3;
         buff.fPols[indx++] = indup + j;
         buff.fPols[indx++] = indup + (j+1)%nlong;
         buff.fPols[indx++] = indpar + j;
      }
   }
   // lower
   if (ndown) {
      for (j=0; j<fNseg; j++) {
         buff.fPols[indx++] = c;
         buff.fPols[indx++] = 3;
         buff.fPols[indx++] = inddown + j;
         buff.fPols[indx++] = indpar + (nlat-1)*fNseg + j;
         buff.fPols[indx++] = inddown + (j+1)%nlong;
      }
   }
   // Fill polygons for inside sphere (except 0/180)

   if (TestShapeBit(kGeoRSeg)) {
      for (i=0; i<nlat-1; i++) {
         for (j=0; j<fNseg; j++) {
            buff.fPols[indx++] = c+1;
            buff.fPols[indx++] = 4;
            buff.fPols[indx++] = indparin+i*fNseg+j;
            buff.fPols[indx++] = indlongin+i*nlong+j;
            buff.fPols[indx++] = indparin+(i+1)*fNseg+j;
            buff.fPols[indx++] = indlongin+i*nlong+(j+1)%nlong;
         }
      }
      // upper
      if (nup) {
         for (j=0; j<fNseg; j++) {
            buff.fPols[indx++] = c+1;
            buff.fPols[indx++] = 3;
            buff.fPols[indx++] = indupin + j;
            buff.fPols[indx++] = indparin + j;
            buff.fPols[indx++] = indupin + (j+1)%nlong;
         }
      }
      // lower
      if (ndown) {
         for (j=0; j<fNseg; j++) {
            buff.fPols[indx++] = c+1;
            buff.fPols[indx++] = 3;
            buff.fPols[indx++] = inddownin + j;
            buff.fPols[indx++] = inddownin + (j+1)%nlong;
            buff.fPols[indx++] = indparin + (nlat-1)*fNseg + j;
         }
      }
   }
   // Polygons on phi planes
   if (TestShapeBit(kGeoPhiSeg)) {
      for (i=0; i<nlat-1; i++) {
         buff.fPols[indx++]   = c+2;
         if (TestShapeBit(kGeoRSeg)) {
            buff.fPols[indx++] = 4;
            buff.fPols[indx++] = indlong + i*nlong;
            buff.fPols[indx++] = indphi + i + 1;
            buff.fPols[indx++] = indlongin + i*nlong;
            buff.fPols[indx++] = indphi + i;
         } else {
            buff.fPols[indx++] = 3;
            buff.fPols[indx++] = indlong + i*nlong;
            buff.fPols[indx++] = indphi + i + 1;
            buff.fPols[indx++] = indphi + i;
         }
      }
      for (i=0; i<nlat-1; i++) {
         buff.fPols[indx++]   = c+2;
         if (TestShapeBit(kGeoRSeg)) {
            buff.fPols[indx++] = 4;
            buff.fPols[indx++] = indlong + (i+1)*nlong-1;
            buff.fPols[indx++] = indphi + nlat + i;
            buff.fPols[indx++] = indlongin + (i+1)*nlong-1;
            buff.fPols[indx++] = indphi + nlat + i + 1;
         } else {
            buff.fPols[indx++] = 3;
            buff.fPols[indx++] = indlong + (i+1)*nlong-1;
            buff.fPols[indx++] = indphi + nlat + i;
            buff.fPols[indx++] = indphi + nlat + i + 1;
         }
      }
      if (nup) {
         buff.fPols[indx++]   = c+2;
         if (TestShapeBit(kGeoRSeg)) {
            buff.fPols[indx++] = 4;
            buff.fPols[indx++] = indup;
            buff.fPols[indx++] = indphi;
            buff.fPols[indx++] = indupin;
            buff.fPols[indx++] = indphi + 2*nlat;
         } else {
            buff.fPols[indx++] = 3;
            buff.fPols[indx++] = indup;
            buff.fPols[indx++] = indphi;
            buff.fPols[indx++] = indphi + 2*nlat;
         }
         buff.fPols[indx++]   = c+2;
         if (TestShapeBit(kGeoRSeg)) {
            buff.fPols[indx++] = 4;
            buff.fPols[indx++] = indup+nlong-1;
            buff.fPols[indx++] = indphi + 2*nlat;
            buff.fPols[indx++] = indupin+nlong-1;
            buff.fPols[indx++] = indphi + nlat;
         } else {
            buff.fPols[indx++] = 3;
            buff.fPols[indx++] = indup+nlong-1;
            buff.fPols[indx++] = indphi + 2*nlat;
            buff.fPols[indx++] = indphi + nlat;
         }
      }
      if (ndown) {
         buff.fPols[indx++]   = c+2;
         if (TestShapeBit(kGeoRSeg)) {
            buff.fPols[indx++] = 4;
            buff.fPols[indx++] = inddown;
            buff.fPols[indx++] = indphi + 2*nlat + nup;
            buff.fPols[indx++] = inddownin;
            buff.fPols[indx++] = indphi + nlat-1;
         } else {
            buff.fPols[indx++] = 3;
            buff.fPols[indx++] = inddown;
            buff.fPols[indx++] = indphi + 2*nlat + nup;
            buff.fPols[indx++] = indphi + nlat-1;
         }
         buff.fPols[indx++]   = c+2;
         if (TestShapeBit(kGeoRSeg)) {
            buff.fPols[indx++] = 4;
            buff.fPols[indx++] = inddown+nlong-1;
            buff.fPols[indx++] = indphi + 2*nlat-1;
            buff.fPols[indx++] = inddownin+nlong-1;
            buff.fPols[indx++] = indphi + 2*nlat+nup;
         } else {
            buff.fPols[indx++] = 3;
            buff.fPols[indx++] = inddown+nlong-1;
            buff.fPols[indx++] = indphi + 2*nlat-1;
            buff.fPols[indx++] = indphi + 2*nlat+nup;
         }
      }
   }
   // Polygons on cones
   if (!nup) {
      for (j=0; j<fNseg; j++) {
         buff.fPols[indx++] = c+2;
         if (TestShapeBit(kGeoRSeg)) {
            buff.fPols[indx++] = 4;
            buff.fPols[indx++] = indpar+j;
            buff.fPols[indx++] = indtheta + j;
            buff.fPols[indx++] = indparin + j;
            buff.fPols[indx++] = indtheta + (j+1)%nlong;
         } else {
            buff.fPols[indx++] = 3;
            buff.fPols[indx++] = indpar+j;
            buff.fPols[indx++] = indtheta + j;
            buff.fPols[indx++] = indtheta + (j+1)%nlong;
         }
      }
   }
   if (!ndown) {
      for (j=0; j<fNseg; j++) {
         buff.fPols[indx++] = c+2;
         if (TestShapeBit(kGeoRSeg)) {
            buff.fPols[indx++] = 4;
            buff.fPols[indx++] = indpar+(nlat-1)*fNseg+j;
            buff.fPols[indx++] = indtheta + (1-nup)*nlong +(j+1)%nlong;
            buff.fPols[indx++] = indparin + (nlat-1)*fNseg + j;
            buff.fPols[indx++] = indtheta + (1-nup)*nlong + j;
         } else {
            buff.fPols[indx++] = 3;
            buff.fPols[indx++] = indpar+(nlat-1)*fNseg+j;
            buff.fPols[indx++] = indtheta + (1-nup)*nlong +(j+1)%nlong;
            buff.fPols[indx++] = indtheta + (1-nup)*nlong + j;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// computes the closest distance from given point to this shape, according
/// to option. The matching point on the shape is stored in spoint.

Double_t TGeoSphere::Safety(const Double_t *point, Bool_t in) const
{
   Double_t r2 = point[0]*point[0]+point[1]*point[1]+point[2]*point[2];
   Double_t r=TMath::Sqrt(r2);
   Bool_t rzero=kFALSE;
   if (r<=1E-20) rzero=kTRUE;
   //localize theta
   Double_t th=0.;
   if (TestShapeBit(kGeoThetaSeg) && (!rzero)) {
      th = TMath::ACos(point[2]/r)*TMath::RadToDeg();
   }
   Double_t saf[4];
   saf[0]=(TGeoShape::IsSameWithinTolerance(fRmin,0) && !TestShapeBit(kGeoThetaSeg) && !TestShapeBit(kGeoPhiSeg))?TGeoShape::Big():r-fRmin;
   saf[1]=fRmax-r;
   saf[2]=saf[3]= TGeoShape::Big();
   if (TestShapeBit(kGeoThetaSeg)) {
      if (fTheta1>0)    saf[2] = r*TMath::Sin((th-fTheta1)*TMath::DegToRad());
      if (fTheta2<180)  saf[3] = r*TMath::Sin((fTheta2-th)*TMath::DegToRad());
   }
   Double_t safphi = TGeoShape::Big();
   Double_t safe = TGeoShape::Big();
   if (TestShapeBit(kGeoPhiSeg)) safphi = TGeoShape::SafetyPhi(point,in,fPhi1,fPhi2);
   if (in) {
      safe = saf[TMath::LocMin(4,saf)];
      return TMath::Min(safe,safphi);
   }
   for (Int_t i=0; i<4; i++) saf[i]=-saf[i];
   safe = saf[TMath::LocMax(4, saf)];
   if (TestShapeBit(kGeoPhiSeg)) return TMath::Max(safe, safphi);
   return safe;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoSphere::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive)) return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   rmin   = " << fRmin << ";" << std::endl;
   out << "   rmax   = " << fRmax << ";" << std::endl;
   out << "   theta1 = " << fTheta1<< ";" << std::endl;
   out << "   theta2 = " << fTheta2 << ";" << std::endl;
   out << "   phi1   = " << fPhi1 << ";" << std::endl;
   out << "   phi2   = " << fPhi2 << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoSphere(\"" << GetName() << "\",rmin,rmax,theta1, theta2,phi1,phi2);" << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Set spherical segment dimensions.

void TGeoSphere::SetSphDimensions(Double_t rmin, Double_t rmax, Double_t theta1,
                               Double_t theta2, Double_t phi1, Double_t phi2)
{
   if (rmin >= rmax) {
      Error("SetDimensions", "invalid parameters rmin/rmax");
      return;
   }
   fRmin = rmin;
   fRmax = rmax;
   if (rmin>0) SetShapeBit(kGeoRSeg);
   if (theta1 >= theta2 || theta1<0 || theta1>180 || theta2>180) {
      Error("SetDimensions", "invalid parameters theta1/theta2");
      return;
   }
   fTheta1 = theta1;
   fTheta2 = theta2;
   if ((theta2-theta1)<180.) SetShapeBit(kGeoThetaSeg);
   fPhi1 = phi1;
   if (phi1<0) fPhi1+=360.;
   fPhi2 = phi2;
   while (fPhi2<=fPhi1) fPhi2+=360.;
   if (!TGeoShape::IsSameWithinTolerance(TMath::Abs(phi2-phi1),360)) SetShapeBit(kGeoPhiSeg);
}

////////////////////////////////////////////////////////////////////////////////
/// Set dimensions of the spherical segment starting from a list of parameters.

void TGeoSphere::SetDimensions(Double_t *param, Int_t nparam)
{
   Double_t rmin = param[0];
   Double_t rmax = param[1];
   Double_t theta1 = 0;
   Double_t theta2 = 180.;
   Double_t phi1 = 0;
   Double_t phi2 = 360.;
   if (nparam > 2) theta1 = param[2];
   if (nparam > 3) theta2 = param[3];
   if (nparam > 4) phi1   = param[4];
   if (nparam > 5) phi2   = param[5];
   SetSphDimensions(rmin, rmax, theta1, theta2, phi1, phi2);
}

////////////////////////////////////////////////////////////////////////////////
/// Set dimensions of the spherical segment starting from a list of parameters.
/// Only takes rmin and rmax

void TGeoSphere::SetDimensions(Double_t *param)
{
   SetDimensions(param,2);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the number of divisions of mesh circles keeping aspect ratio.

void TGeoSphere::SetNumberOfDivisions(Int_t p)
{
   fNseg = p;
   Double_t dphi = fPhi2 - fPhi1;
   if (dphi<0) dphi+=360;
   Double_t dtheta = TMath::Abs(fTheta2-fTheta1);
   fNz = Int_t(fNseg*dtheta/dphi) +1;
   if (fNz<2) fNz=2;
}

////////////////////////////////////////////////////////////////////////////////
/// create sphere mesh points

void TGeoSphere::SetPoints(Double_t *points) const
{
   if (!points) {
      Error("SetPoints", "Input array is NULL");
      return;
   }
   Bool_t full = kTRUE;
   if (TestShapeBit(kGeoThetaSeg) || TestShapeBit(kGeoPhiSeg)) full = kFALSE;
   Int_t ncenter = 1;
   if (full || TestShapeBit(kGeoRSeg)) ncenter = 0;
   Int_t nup = (fTheta1>0)?0:1;
   Int_t ndown = (fTheta2<180)?0:1;
   // number of different latitudes, excluding 0 and 180 degrees
   Int_t nlat = fNz+1-(nup+ndown);
   // number of different longitudes
   Int_t nlong = fNseg;
   if (TestShapeBit(kGeoPhiSeg)) nlong++;
   // total number of points on mesh is:
   //    nlat*nlong + nup + ndown + ncenter;    // in case rmin=0
   //   2*(nlat*nlong + nup + ndown);           // in case rmin>0
   Int_t i,j ;
   Double_t phi1 = fPhi1*TMath::DegToRad();
   Double_t phi2 = fPhi2*TMath::DegToRad();
   Double_t dphi = (phi2-phi1)/fNseg;
   Double_t theta1 = fTheta1*TMath::DegToRad();
   Double_t theta2 = fTheta2*TMath::DegToRad();
   Double_t dtheta = (theta2-theta1)/fNz;
   Double_t z,zi,theta,phi,cphi,sphi;
   Int_t indx=0;
   // FILL ALL POINTS ON OUTER SPHERE
   // (nlat * nlong) points
   // loop all latitudes except 0/180 degrees (nlat times)
   // ilat = [0,nlat]   jlong = [0,nlong]
   // Index(ilat, jlong) = 3*(ilat*nlat + jlong)
   for (i = 0; i < nlat; i++) {
      theta = theta1+(nup+i)*dtheta;
      z =  fRmax * TMath::Cos(theta);
      zi = fRmax * TMath::Sin(theta);
      // loop all different longitudes (nlong times)
      for (j = 0; j < nlong; j++) {
         phi = phi1+j*dphi;
         cphi = TMath::Cos(phi);
         sphi = TMath::Sin(phi);
         points[indx++] = zi * cphi;
         points[indx++] = zi * sphi;
         points[indx++] = z;
      }
   }
   // upper/lower points (if they exist) for outer sphere
   if (nup) {
      // ind_up = 3*nlat*nlong
      points[indx++] = 0.;
      points[indx++] = 0.;
      points[indx++] = fRmax;
   }
   if (ndown) {
      // ind_down = 3*(nlat*nlong+nup)
      points[indx++] = 0.;
      points[indx++] = 0.;
      points[indx++] = -fRmax;
   }
   // do the same for inner sphere if it exist
   // Start_index = 3*(nlat*nlong + nup + ndown)
   if (TestShapeBit(kGeoRSeg)) {
   // Index(ilat, jlong) = start_index + 3*(ilat*nlat + jlong)
      for (i = 0; i < nlat; i++) {
         theta = theta1+(nup+i)*dtheta;
         z =  fRmin * TMath::Cos(theta);
         zi = fRmin * TMath::Sin(theta);
         // loop all different longitudes (nlong times)
         for (j = 0; j < nlong; j++) {
            phi = phi1+j*dphi;
            cphi = TMath::Cos(phi);
            sphi = TMath::Sin(phi);
            points[indx++] = zi * cphi;
            points[indx++] = zi * sphi;
            points[indx++] = z;
         }
      }
      // upper/lower points (if they exist) for inner sphere
      if (nup) {
      // ind_up = start_index + 3*nlat*nlong
         points[indx++] = 0.;
         points[indx++] = 0.;
         points[indx++] = fRmin;
      }
      if (ndown) {
      // ind_down = start_index + 3*(nlat*nlong+nup)
         points[indx++] = 0.;
         points[indx++] = 0.;
         points[indx++] = -fRmin;
      }
   }
   // Add center of sphere if needed
   if (ncenter) {
      // ind_center = 6*(nlat*nlong + nup + ndown)
      points[indx++] = 0.;
      points[indx++] = 0.;
      points[indx++] = 0.;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// create sphere mesh points

void TGeoSphere::SetPoints(Float_t *points) const
{
   if (!points) {
      Error("SetPoints", "Input array is NULL");
      return;
   }
   Bool_t full = kTRUE;
   if (TestShapeBit(kGeoThetaSeg) || TestShapeBit(kGeoPhiSeg)) full = kFALSE;
   Int_t ncenter = 1;
   if (full || TestShapeBit(kGeoRSeg)) ncenter = 0;
   Int_t nup = (fTheta1>0)?0:1;
   Int_t ndown = (fTheta2<180)?0:1;
   // number of different latitudes, excluding 0 and 180 degrees
   Int_t nlat = fNz+1-(nup+ndown);
   // number of different longitudes
   Int_t nlong = fNseg;
   if (TestShapeBit(kGeoPhiSeg)) nlong++;
   // total number of points on mesh is:
   //    nlat*nlong + nup + ndown + ncenter;    // in case rmin=0
   //   2*(nlat*nlong + nup + ndown);           // in case rmin>0
   Int_t i,j ;
   Double_t phi1 = fPhi1*TMath::DegToRad();
   Double_t phi2 = fPhi2*TMath::DegToRad();
   Double_t dphi = (phi2-phi1)/fNseg;
   Double_t theta1 = fTheta1*TMath::DegToRad();
   Double_t theta2 = fTheta2*TMath::DegToRad();
   Double_t dtheta = (theta2-theta1)/fNz;
   Double_t z,zi,theta,phi,cphi,sphi;
   Int_t indx=0;
   // FILL ALL POINTS ON OUTER SPHERE
   // (nlat * nlong) points
   // loop all latitudes except 0/180 degrees (nlat times)
   // ilat = [0,nlat]   jlong = [0,nlong]
   // Index(ilat, jlong) = 3*(ilat*nlat + jlong)
   for (i = 0; i < nlat; i++) {
      theta = theta1+(nup+i)*dtheta;
      z =  fRmax * TMath::Cos(theta);
      zi = fRmax * TMath::Sin(theta);
      // loop all different longitudes (nlong times)
      for (j = 0; j < nlong; j++) {
         phi = phi1+j*dphi;
         cphi = TMath::Cos(phi);
         sphi = TMath::Sin(phi);
         points[indx++] = zi * cphi;
         points[indx++] = zi * sphi;
         points[indx++] = z;
      }
   }
   // upper/lower points (if they exist) for outer sphere
   if (nup) {
      // ind_up = 3*nlat*nlong
      points[indx++] = 0.;
      points[indx++] = 0.;
      points[indx++] = fRmax;
   }
   if (ndown) {
      // ind_down = 3*(nlat*nlong+nup)
      points[indx++] = 0.;
      points[indx++] = 0.;
      points[indx++] = -fRmax;
   }
   // do the same for inner sphere if it exist
   // Start_index = 3*(nlat*nlong + nup + ndown)
   if (TestShapeBit(kGeoRSeg)) {
   // Index(ilat, jlong) = start_index + 3*(ilat*nlat + jlong)
      for (i = 0; i < nlat; i++) {
         theta = theta1+(nup+i)*dtheta;
         z =  fRmin * TMath::Cos(theta);
         zi = fRmin * TMath::Sin(theta);
         // loop all different longitudes (nlong times)
         for (j = 0; j < nlong; j++) {
            phi = phi1+j*dphi;
            cphi = TMath::Cos(phi);
            sphi = TMath::Sin(phi);
            points[indx++] = zi * cphi;
            points[indx++] = zi * sphi;
            points[indx++] = z;
         }
      }
      // upper/lower points (if they exist) for inner sphere
      if (nup) {
      // ind_up = start_index + 3*nlat*nlong
         points[indx++] = 0.;
         points[indx++] = 0.;
         points[indx++] = fRmin;
      }
      if (ndown) {
      // ind_down = start_index + 3*(nlat*nlong+nup)
         points[indx++] = 0.;
         points[indx++] = 0.;
         points[indx++] = -fRmin;
      }
   }
   // Add center of sphere if needed
   if (ncenter) {
      // ind_center = 6*(nlat*nlong + nup + ndown)
      points[indx++] = 0.;
      points[indx++] = 0.;
      points[indx++] = 0.;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoSphere::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   TGeoSphere * localThis = const_cast<TGeoSphere *>(this);
   localThis->SetNumberOfDivisions(gGeoManager->GetNsegments());
   Bool_t full = kTRUE;
   if (TestShapeBit(kGeoThetaSeg) || TestShapeBit(kGeoPhiSeg)) full = kFALSE;
   Int_t ncenter = 1;
   if (full || TestShapeBit(kGeoRSeg)) ncenter = 0;
   Int_t nup = (fTheta1>0)?0:1;
   Int_t ndown = (fTheta2<180)?0:1;
   // number of different latitudes, excluding 0 and 180 degrees
   Int_t nlat = fNz+1-(nup+ndown);
   // number of different longitudes
   Int_t nlong = fNseg;
   if (TestShapeBit(kGeoPhiSeg)) nlong++;

   nvert = nlat*nlong+nup+ndown+ncenter;
   if (TestShapeBit(kGeoRSeg)) nvert *= 2;

   nsegs = nlat*fNseg + (nlat-1+nup+ndown)*nlong; // outer sphere
   if (TestShapeBit(kGeoRSeg)) nsegs *= 2; // inner sphere
   if (TestShapeBit(kGeoPhiSeg)) nsegs += 2*nlat+nup+ndown; // 2 phi planes
   nsegs += nlong * (2-nup - ndown);  // connecting cones

   npols = fNz*fNseg; // outer
   if (TestShapeBit(kGeoRSeg)) npols *=2;  // inner
   if (TestShapeBit(kGeoPhiSeg)) npols += 2*fNz; // 2 phi planes
   npols += (2-nup-ndown)*fNseg; // connecting
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of vertices of the mesh representation

Int_t TGeoSphere::GetNmeshVertices() const
{
   Bool_t full = kTRUE;
   if (TestShapeBit(kGeoThetaSeg) || TestShapeBit(kGeoPhiSeg)) full = kFALSE;
   Int_t ncenter = 1;
   if (full || TestShapeBit(kGeoRSeg)) ncenter = 0;
   Int_t nup = (fTheta1>0)?0:1;
   Int_t ndown = (fTheta2<180)?0:1;
   // number of different latitudes, excluding 0 and 180 degrees
   Int_t nlat = fNz+1-(nup+ndown);
   // number of different longitudes
   Int_t nlong = fNseg;
   if (TestShapeBit(kGeoPhiSeg)) nlong++;
   // total number of points on mesh is:
   //    nlat*nlong + nup + ndown + ncenter;    // in case rmin=0
   //   2*(nlat*nlong + nup + ndown);           // in case rmin>0
   Int_t numPoints = 0;
   if (TestShapeBit(kGeoRSeg)) numPoints = 2*(nlat*nlong+nup+ndown);
   else numPoints = nlat*nlong+nup+ndown+ncenter;
   return numPoints;
}

////////////////////////////////////////////////////////////////////////////////
////// obsolete - to be removed

void TGeoSphere::Sizeof3D() const
{
}

////////////////////////////////////////////////////////////////////////////////
/// Fills a static 3D buffer and returns a reference.

const TBuffer3D & TGeoSphere::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
{
   static TBuffer3DSphere buffer;

   TGeoBBox::FillBuffer3D(buffer, reqSections, localFrame);

   if (reqSections & TBuffer3D::kShapeSpecific) {
      buffer.fRadiusInner  = fRmin;
      buffer.fRadiusOuter  = fRmax;
      buffer.fThetaMin     = fTheta1;
      buffer.fThetaMax     = fTheta2;
      buffer.fPhiMin       = fPhi1;
      buffer.fPhiMax       = fPhi2;
      buffer.SetSectionsValid(TBuffer3D::kShapeSpecific);
   }
   if (reqSections & TBuffer3D::kRawSizes) {
      // We want FillBuffer to be const
      TGeoSphere * localThis = const_cast<TGeoSphere *>(this);
      localThis->SetNumberOfDivisions(gGeoManager->GetNsegments());

      Bool_t full = kTRUE;
      if (TestShapeBit(kGeoThetaSeg) || TestShapeBit(kGeoPhiSeg)) full = kFALSE;
      Int_t ncenter = 1;
      if (full || TestShapeBit(kGeoRSeg)) ncenter = 0;
      Int_t nup = (fTheta1>0)?0:1;
      Int_t ndown = (fTheta2<180)?0:1;
      // number of different latitudes, excluding 0 and 180 degrees
      Int_t nlat = fNz+1-(nup+ndown);
      // number of different longitudes
      Int_t nlong = fNseg;
      if (TestShapeBit(kGeoPhiSeg)) nlong++;

      Int_t nbPnts = nlat*nlong+nup+ndown+ncenter;
      if (TestShapeBit(kGeoRSeg)) nbPnts *= 2;

      Int_t nbSegs = nlat*fNseg + (nlat-1+nup+ndown)*nlong; // outer sphere
      if (TestShapeBit(kGeoRSeg)) nbSegs *= 2; // inner sphere
      if (TestShapeBit(kGeoPhiSeg)) nbSegs += 2*nlat+nup+ndown; // 2 phi planes
      nbSegs += nlong * (2-nup - ndown);  // connecting cones

      Int_t nbPols = fNz*fNseg; // outer
      if (TestShapeBit(kGeoRSeg)) nbPols *=2;  // inner
      if (TestShapeBit(kGeoPhiSeg)) nbPols += 2*fNz; // 2 phi planes
      nbPols += (2-nup-ndown)*fNseg; // connecting

      if (buffer.SetRawSizes(nbPnts, 3*nbPnts, nbSegs, 3*nbSegs, nbPols, 6*nbPols)) {
         buffer.SetSectionsValid(TBuffer3D::kRawSizes);
      }
   }
   if ((reqSections & TBuffer3D::kRaw) && buffer.SectionsValid(TBuffer3D::kRawSizes)) {
      SetPoints(buffer.fPnts);
      if (!buffer.fLocalFrame) {
         TransformPoints(buffer.fPnts, buffer.NbPnts());
      }
      SetSegsAndPols(buffer);
      buffer.SetSectionsValid(TBuffer3D::kRaw);
   }

   return buffer;
}

////////////////////////////////////////////////////////////////////////////////
/// Check the inside status for each of the points in the array.
/// Input: Array of point coordinates + vector size
/// Output: Array of Booleans for the inside of each point

void TGeoSphere::Contains_v(const Double_t *points, Bool_t *inside, Int_t vecsize) const
{
   for (Int_t i=0; i<vecsize; i++) inside[i] = Contains(&points[3*i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the normal for an array o points so that norm.dot.dir is positive
/// Input: Arrays of point coordinates and directions + vector size
/// Output: Array of normal directions

void TGeoSphere::ComputeNormal_v(const Double_t *points, const Double_t *dirs, Double_t *norms, Int_t vecsize)
{
   for (Int_t i=0; i<vecsize; i++) ComputeNormal(&points[3*i], &dirs[3*i], &norms[3*i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoSphere::DistFromInside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromInside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoSphere::DistFromOutside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromOutside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute safe distance from each of the points in the input array.
/// Input: Array of point coordinates, array of statuses for these points, size of the arrays
/// Output: Safety values

void TGeoSphere::Safety_v(const Double_t *points, const Bool_t *inside, Double_t *safe, Int_t vecsize) const
{
   for (Int_t i=0; i<vecsize; i++) safe[i] = Safety(&points[3*i], inside[i]);
}
