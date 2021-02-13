// @(#)root/geom:$Id$
// Author: Andrei Gheata   24/10/01
// TGeoPcon::Contains() implemented by Mihaela Gheata

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoPcon
\ingroup Geometry_classes
A polycone.

It has at least 9 parameters:
           - the lower phi limit;
           - the range in phi;
           - the number of z planes (at least two) where the inner/outer
             radii are changing;
           - z coordinate, inner and outer radius for each z plane

Begin_Macro(source)
{
   TCanvas *c = new TCanvas("c", "c",0,0,600,600);
   new TGeoManager("pcon", "poza10");
   TGeoMaterial *mat = new TGeoMaterial("Al", 26.98,13,2.7);
   TGeoMedium *med = new TGeoMedium("MED",1,mat);
   TGeoVolume *top = gGeoManager->MakeBox("TOP",med,100,100,100);
   gGeoManager->SetTopVolume(top);
   TGeoVolume *vol = gGeoManager->MakePcon("PCON",med, -30.0,300,4);
   TGeoPcon *pcon = (TGeoPcon*)(vol->GetShape());
   pcon->DefineSection(0,0,15,20);
   pcon->DefineSection(1,20,15,20);
   pcon->DefineSection(2,20,15,25);
   pcon->DefineSection(3,50,15,20);
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

#include "TGeoPcon.h"

#include <iostream>

#include "TBuffer.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TVirtualGeoPainter.h"
#include "TGeoTube.h"
#include "TGeoCone.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"

ClassImp(TGeoPcon);

////////////////////////////////////////////////////////////////////////////////
/// dummy ctor

TGeoPcon::TGeoPcon()
         :TGeoBBox(),
          fNz(0),
          fPhi1(0.),
          fDphi(0.),
          fRmin(nullptr),
          fRmax(nullptr),
          fZ(nullptr),
          fFullPhi(kFALSE),
          fC1(0.),
          fS1(0.),
          fC2(0.),
          fS2(0.),
          fCm(0.),
          fSm(0.),
          fCdphi(0.)
{
   SetShapeBit(TGeoShape::kGeoPcon);
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPcon::TGeoPcon(Double_t phi, Double_t dphi, Int_t nz)
         :TGeoBBox(0, 0, 0),
          fNz(nz),
          fPhi1(phi),
          fDphi(dphi),
          fRmin(nullptr),
          fRmax(nullptr),
          fZ(nullptr),
          fFullPhi(kFALSE),
          fC1(0.),
          fS1(0.),
          fC2(0.),
          fS2(0.),
          fCm(0.),
          fSm(0.),
          fCdphi(0.)
{
   SetShapeBit(TGeoShape::kGeoPcon);
   while (fPhi1<0) fPhi1+=360.;
   fRmin = new Double_t [nz];
   fRmax = new Double_t [nz];
   fZ    = new Double_t [nz];
   memset(fRmin, 0, nz*sizeof(Double_t));
   memset(fRmax, 0, nz*sizeof(Double_t));
   memset(fZ, 0, nz*sizeof(Double_t));
   if (TGeoShape::IsSameWithinTolerance(fDphi,360)) fFullPhi = kTRUE;
   Double_t phi1 = fPhi1;
   Double_t phi2 = phi1+fDphi;
   Double_t phim = 0.5*(phi1+phi2);
   fC1 = TMath::Cos(phi1*TMath::DegToRad());
   fS1 = TMath::Sin(phi1*TMath::DegToRad());
   fC2 = TMath::Cos(phi2*TMath::DegToRad());
   fS2 = TMath::Sin(phi2*TMath::DegToRad());
   fCm = TMath::Cos(phim*TMath::DegToRad());
   fSm = TMath::Sin(phim*TMath::DegToRad());
   fCdphi = TMath::Cos(0.5*fDphi*TMath::DegToRad());
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoPcon::TGeoPcon(const char *name, Double_t phi, Double_t dphi, Int_t nz)
         :TGeoBBox(name, 0, 0, 0),
          fNz(nz),
          fPhi1(phi),
          fDphi(dphi),
          fRmin(nullptr),
          fRmax(nullptr),
          fZ(nullptr),
          fFullPhi(kFALSE),
          fC1(0.),
          fS1(0.),
          fC2(0.),
          fS2(0.),
          fCm(0.),
          fSm(0.),
          fCdphi(0.)
{
   SetShapeBit(TGeoShape::kGeoPcon);
   while (fPhi1<0) fPhi1+=360.;
   fRmin = new Double_t [nz];
   fRmax = new Double_t [nz];
   fZ    = new Double_t [nz];
   memset(fRmin, 0, nz*sizeof(Double_t));
   memset(fRmax, 0, nz*sizeof(Double_t));
   memset(fZ, 0, nz*sizeof(Double_t));
   if (TGeoShape::IsSameWithinTolerance(fDphi,360)) fFullPhi = kTRUE;
   Double_t phi1 = fPhi1;
   Double_t phi2 = phi1+fDphi;
   Double_t phim = 0.5*(phi1+phi2);
   fC1 = TMath::Cos(phi1*TMath::DegToRad());
   fS1 = TMath::Sin(phi1*TMath::DegToRad());
   fC2 = TMath::Cos(phi2*TMath::DegToRad());
   fS2 = TMath::Sin(phi2*TMath::DegToRad());
   fCm = TMath::Cos(phim*TMath::DegToRad());
   fSm = TMath::Sin(phim*TMath::DegToRad());
   fCdphi = TMath::Cos(0.5*fDphi*TMath::DegToRad());
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor in GEANT3 style
///  - param[0] = phi1
///  - param[1] = dphi
///  - param[2] = nz
///  - param[3] = z1
///  - param[4] = Rmin1
///  - param[5] = Rmax1
/// ...

TGeoPcon::TGeoPcon(Double_t *param)
         :TGeoBBox(0, 0, 0),
          fNz(0),
          fPhi1(0.),
          fDphi(0.),
          fRmin(nullptr),
          fRmax(nullptr),
          fZ(nullptr),
          fFullPhi(kFALSE),
          fC1(0.),
          fS1(0.),
          fC2(0.),
          fS2(0.),
          fCm(0.),
          fSm(0.),
          fCdphi(0.)
{
   SetShapeBit(TGeoShape::kGeoPcon);
   SetDimensions(param);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoPcon::~TGeoPcon()
{
   if (fRmin) { delete[] fRmin; fRmin = nullptr; }
   if (fRmax) { delete[] fRmax; fRmax = nullptr; }
   if (fZ)    { delete[] fZ; fZ = nullptr; }
}

////////////////////////////////////////////////////////////////////////////////
/// Computes capacity of the shape in [length^3]

Double_t TGeoPcon::Capacity() const
{
   Int_t ipl;
   Double_t rmin1, rmax1, rmin2, rmax2, phi1, phi2, dz;
   Double_t capacity = 0.;
   phi1 = fPhi1;
   phi2 = fPhi1 + fDphi;
   for (ipl=0; ipl<fNz-1; ipl++) {
      dz    = 0.5*(fZ[ipl+1]-fZ[ipl]);
      if (dz < TGeoShape::Tolerance()) continue;
      rmin1 = fRmin[ipl];
      rmax1 = fRmax[ipl];
      rmin2 = fRmin[ipl+1];
      rmax2 = fRmax[ipl+1];
      capacity += TGeoConeSeg::Capacity(dz,rmin1,rmax1,rmin2,rmax2,phi1,phi2);
   }
   return capacity;
}

////////////////////////////////////////////////////////////////////////////////
/// compute bounding box of the pcon
/// Check if the sections are in increasing Z order

void TGeoPcon::ComputeBBox()
{
   for (Int_t isec=0; isec<fNz-1; isec++) {
      if (TMath::Abs(fZ[isec]-fZ[isec+1]) < TGeoShape::Tolerance()) {
         fZ[isec+1]=fZ[isec];
         if (IsSameWithinTolerance(fRmin[isec], fRmin[isec+1]) &&
             IsSameWithinTolerance(fRmax[isec], fRmax[isec+1])) {
            InspectShape();
            Error("ComputeBBox", "Duplicated section %d/%d for shape %s", isec, isec+1, GetName());
         }
      }
      if (fZ[isec]>fZ[isec+1]) {
         InspectShape();
         Fatal("ComputeBBox", "Wrong section order");
      }
   }
   // Check if the last sections are valid
   if (TMath::Abs(fZ[1]-fZ[0]) < TGeoShape::Tolerance() ||
       TMath::Abs(fZ[fNz-1]-fZ[fNz-2]) < TGeoShape::Tolerance()) {
      InspectShape();
      Fatal("ComputeBBox","Shape %s at index %d: Not allowed first two or last two sections at same Z",
             GetName(), gGeoManager->GetListOfShapes()->IndexOf(this));
   }
   Double_t zmin = TMath::Min(fZ[0], fZ[fNz-1]);
   Double_t zmax = TMath::Max(fZ[0], fZ[fNz-1]);
   // find largest rmax an smallest rmin
   Double_t rmin, rmax;
   rmin = fRmin[TMath::LocMin(fNz, fRmin)];
   rmax = fRmax[TMath::LocMax(fNz, fRmax)];

   Double_t xc[4];
   Double_t yc[4];
   xc[0] = rmax*fC1;
   yc[0] = rmax*fS1;
   xc[1] = rmax*fC2;
   yc[1] = rmax*fS2;
   xc[2] = rmin*fC1;
   yc[2] = rmin*fS1;
   xc[3] = rmin*fC2;
   yc[3] = rmin*fS2;

   Double_t xmin = xc[TMath::LocMin(4, &xc[0])];
   Double_t xmax = xc[TMath::LocMax(4, &xc[0])];
   Double_t ymin = yc[TMath::LocMin(4, &yc[0])];
   Double_t ymax = yc[TMath::LocMax(4, &yc[0])];

   Double_t ddp = -fPhi1;
   if (ddp<0) ddp+= 360;
   if (ddp<=fDphi) xmax = rmax;
   ddp = 90-fPhi1;
   if (ddp<0) ddp+= 360;
   if (ddp<=fDphi) ymax = rmax;
   ddp = 180-fPhi1;
   if (ddp<0) ddp+= 360;
   if (ddp<=fDphi) xmin = -rmax;
   ddp = 270-fPhi1;
   if (ddp<0) ddp+= 360;
   if (ddp<=fDphi) ymin = -rmax;
   fOrigin[0] = (xmax+xmin)/2;
   fOrigin[1] = (ymax+ymin)/2;
   fOrigin[2] = (zmax+zmin)/2;
   fDX = (xmax-xmin)/2;
   fDY = (ymax-ymin)/2;
   fDZ = (zmax-zmin)/2;
   SetShapeBit(kGeoClosedShape);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute normal to closest surface from POINT.

void TGeoPcon::ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm)
{
   memset(norm,0,3*sizeof(Double_t));
   Double_t r;
   Double_t ptnew[3];
   Double_t dz, rmin1, rmax1, rmin2, rmax2;
   Bool_t is_tube;
   Int_t ipl = TMath::BinarySearch(fNz, fZ, point[2]);
   if (ipl==(fNz-1) || ipl<0) {
      // point outside Z range
      norm[2] = TMath::Sign(1., dir[2]);
      return;
   }
   Int_t iplclose = ipl;
   if ((fZ[ipl+1]-point[2])<(point[2]-fZ[ipl])) iplclose++;
   dz = TMath::Abs(fZ[iplclose]-point[2]);
   if (dz<1E-5) {
      if (iplclose==0 || iplclose==(fNz-1)) {
         norm[2] = TMath::Sign(1., dir[2]);
         return;
      }
      if (iplclose==ipl && TGeoShape::IsSameWithinTolerance(fZ[ipl],fZ[ipl-1])) {
         r = TMath::Sqrt(point[0]*point[0]+point[1]*point[1]);
         if (r<TMath::Max(fRmin[ipl],fRmin[ipl-1]) || r>TMath::Min(fRmax[ipl],fRmax[ipl-1])) {
            norm[2] = TMath::Sign(1., dir[2]);
            return;
         }
      } else {
         if (TGeoShape::IsSameWithinTolerance(fZ[iplclose],fZ[iplclose+1])) {
            r = TMath::Sqrt(point[0]*point[0]+point[1]*point[1]);
            if (r<TMath::Max(fRmin[iplclose],fRmin[iplclose+1]) || r>TMath::Min(fRmax[iplclose],fRmax[iplclose+1])) {
               norm[2] = TMath::Sign(1., dir[2]);
               return;
            }
         }
      }
   } //-> Z done
   memcpy(ptnew, point, 3*sizeof(Double_t));
   dz = 0.5*(fZ[ipl+1]-fZ[ipl]);
   if (TGeoShape::IsSameWithinTolerance(dz,0.)) {
      norm[2] = TMath::Sign(1., dir[2]);
      return;
   }
   ptnew[2] -= 0.5*(fZ[ipl]+fZ[ipl+1]);
   rmin1 = fRmin[ipl];
   rmax1 = fRmax[ipl];
   rmin2 = fRmin[ipl+1];
   rmax2 = fRmax[ipl+1];
   is_tube = (TGeoShape::IsSameWithinTolerance(rmin1,rmin2) && TGeoShape::IsSameWithinTolerance(rmax1,rmax2))?kTRUE:kFALSE;
   if (!fFullPhi) {
      if (is_tube) TGeoTubeSeg::ComputeNormalS(ptnew,dir,norm,rmin1,rmax1,dz,fC1,fS1,fC2,fS2);
      else         TGeoConeSeg::ComputeNormalS(ptnew,dir,norm,dz,rmin1,rmax1,rmin2,rmax2,fC1,fS1,fC2,fS2);
   } else {
      if (is_tube) TGeoTube::ComputeNormalS(ptnew,dir,norm,rmin1,rmax1,dz);
      else         TGeoCone::ComputeNormalS(ptnew,dir,norm,dz,rmin1,rmax1,rmin2,rmax2);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// test if point is inside this shape
/// check total z range

Bool_t TGeoPcon::Contains(const Double_t *point) const
{
   if ((point[2]<fZ[0]) || (point[2]>fZ[fNz-1])) return kFALSE;
   // check R squared
   Double_t r2 = point[0]*point[0]+point[1]*point[1];

   Int_t izl = 0;
   Int_t izh = fNz-1;
   Int_t izt = (fNz-1)/2;
   while ((izh-izl)>1) {
      if (point[2] > fZ[izt]) izl = izt;
      else izh = izt;
      izt = (izl+izh)>>1;
   }
   // the point is in the section bounded by izl and izh Z planes

   // compute Rmin and Rmax and test the value of R squared
   Double_t rmin, rmax;
   if (TGeoShape::IsSameWithinTolerance(fZ[izl],fZ[izh]) && TGeoShape::IsSameWithinTolerance(point[2],fZ[izl])) {
      rmin = TMath::Min(fRmin[izl], fRmin[izh]);
      rmax = TMath::Max(fRmax[izl], fRmax[izh]);
   } else {
      Double_t dz = fZ[izh] - fZ[izl];
      Double_t dz1 = point[2] - fZ[izl];
      rmin = (fRmin[izl]*(dz-dz1)+fRmin[izh]*dz1)/dz;
      rmax = (fRmax[izl]*(dz-dz1)+fRmax[izh]*dz1)/dz;
   }
   if ((r2<rmin*rmin) || (r2>rmax*rmax)) return kFALSE;
   // now check phi
   if (TGeoShape::IsSameWithinTolerance(fDphi,360)) return kTRUE;
   if (r2<1E-10) return kTRUE;
   Double_t phi = TMath::ATan2(point[1], point[0]) * TMath::RadToDeg();
   if (phi < 0) phi+=360.0;
   Double_t ddp = phi-fPhi1;
   if (ddp<0) ddp+=360.;
   if (ddp<=fDphi) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// compute closest distance from point px,py to each corner

Int_t TGeoPcon::DistancetoPrimitive(Int_t px, Int_t py)
{
   Int_t n = gGeoManager->GetNsegments()+1;
   const Int_t numPoints = 2*n*fNz;
   return ShapeDistancetoPrimitive(numPoints, px, py);
}

////////////////////////////////////////////////////////////////////////////////
/// compute distance from inside point to surface of the polycone

Double_t TGeoPcon::DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
{
   if (iact<3 && safe) {
      *safe = Safety(point, kTRUE);
      if (iact==0) return TGeoShape::Big();
      if ((iact==1) && (*safe>step)) return TGeoShape::Big();
   }
   Double_t snxt = TGeoShape::Big();
   Double_t sstep = 1E-6;
   Double_t point_new[3];
   // determine which z segment contains the point
   Int_t ipl = TMath::BinarySearch(fNz, fZ, point[2]+TMath::Sign(1.E-10,dir[2]));
   if (ipl<0) ipl=0;
   if (ipl==(fNz-1)) ipl--;
   Double_t dz = 0.5*(fZ[ipl+1]-fZ[ipl]);
   Bool_t special_case = kFALSE;
   if (dz<1e-9) {
      // radius changing segment, make sure track is not in the XY plane
      if (TGeoShape::IsSameWithinTolerance(dir[2], 0)) {
         special_case = kTRUE;
      } else {
         //check if a close point is still contained
         point_new[0] = point[0]+sstep*dir[0];
         point_new[1] = point[1]+sstep*dir[1];
         point_new[2] = point[2]+sstep*dir[2];
         if (!Contains(point_new)) return 0.;
         return (DistFromInside(point_new,dir,iact,step,safe)+sstep);
      }
   }
   // determine if the current segment is a tube or a cone
   Bool_t intub = kTRUE;
   if (!TGeoShape::IsSameWithinTolerance(fRmin[ipl],fRmin[ipl+1])) intub=kFALSE;
   else if (!TGeoShape::IsSameWithinTolerance(fRmax[ipl],fRmax[ipl+1])) intub=kFALSE;
   // determine phi segmentation
   memcpy(point_new, point, 2*sizeof(Double_t));
   // new point in reference system of the current segment
   point_new[2] = point[2]-0.5*(fZ[ipl]+fZ[ipl+1]);

   if (special_case) {
      if (!fFullPhi) snxt = TGeoTubeSeg::DistFromInsideS(point_new, dir,
               TMath::Min(fRmin[ipl],fRmin[ipl+1]), TMath::Max(fRmax[ipl],fRmax[ipl+1]),
               dz, fC1,fS1,fC2,fS2,fCm,fSm,fCdphi);
      else       snxt = TGeoTube::DistFromInsideS(point_new, dir,
               TMath::Min(fRmin[ipl],fRmin[ipl+1]), TMath::Max(fRmax[ipl],fRmax[ipl+1]),dz);
      return snxt;
   }
   if (intub) {
      if (!fFullPhi) snxt=TGeoTubeSeg::DistFromInsideS(point_new, dir, fRmin[ipl], fRmax[ipl],dz, fC1,fS1,fC2,fS2,fCm,fSm,fCdphi);
      else snxt=TGeoTube::DistFromInsideS(point_new, dir, fRmin[ipl], fRmax[ipl],dz);
   } else {
      if (!fFullPhi) snxt=TGeoConeSeg::DistFromInsideS(point_new,dir,dz,fRmin[ipl],fRmax[ipl],fRmin[ipl+1],fRmax[ipl+1],fC1,fS1,fC2,fS2,fCm,fSm,fCdphi);
      else snxt=TGeoCone::DistFromInsideS(point_new,dir,dz,fRmin[ipl],fRmax[ipl],fRmin[ipl+1], fRmax[ipl+1]);
   }

   for (Int_t i=0; i<3; i++) point_new[i]=point[i]+(snxt+1E-6)*dir[i];
   if (!Contains(&point_new[0])) return snxt;

   snxt += DistFromInside(&point_new[0], dir, 3) + 1E-6;
   return snxt;
}

////////////////////////////////////////////////////////////////////////////////
/// compute distance to a pcon Z slice. Segment iz must be valid

Double_t TGeoPcon::DistToSegZ(const Double_t *point, const Double_t *dir, Int_t &iz) const
{
   Double_t zmin=fZ[iz];
   Double_t zmax=fZ[iz+1];
   if (TGeoShape::IsSameWithinTolerance(zmin,zmax)) {
      if (TGeoShape::IsSameWithinTolerance(dir[2],0)) return TGeoShape::Big();
      Int_t istep=(dir[2]>0)?1:-1;
      iz+=istep;
      if (iz<0 || iz>(fNz-2)) return TGeoShape::Big();
      return DistToSegZ(point,dir,iz);
   }
   Double_t dz=0.5*(zmax-zmin);
   Double_t local[3];
   memcpy(&local[0], point, 3*sizeof(Double_t));
   local[2]=point[2]-0.5*(zmin+zmax);
   Double_t snxt;
   Double_t rmin1=fRmin[iz];
   Double_t rmax1=fRmax[iz];
   Double_t rmin2=fRmin[iz+1];
   Double_t rmax2=fRmax[iz+1];

   if (TGeoShape::IsSameWithinTolerance(rmin1,rmin2) && TGeoShape::IsSameWithinTolerance(rmax1,rmax2)) {
      if (fFullPhi) snxt=TGeoTube::DistFromOutsideS(local, dir, rmin1, rmax1, dz);
      else snxt=TGeoTubeSeg::DistFromOutsideS(local,dir,rmin1,rmax1,dz,fC1,fS1,fC2,fS2,fCm,fSm,fCdphi);
   } else {
      if (fFullPhi) snxt=TGeoCone::DistFromOutsideS(local,dir,dz,rmin1, rmax1,rmin2,rmax2);
      else snxt=TGeoConeSeg::DistFromOutsideS(local,dir,dz,rmin1,rmax1,rmin2,rmax2,fC1,fS1,fC2,fS2,fCm,fSm,fCdphi);
   }
   if (snxt<1E20) return snxt;
   // check next segment
   if (TGeoShape::IsSameWithinTolerance(dir[2],0)) return TGeoShape::Big();
   Int_t istep=(dir[2]>0)?1:-1;
   iz+=istep;
   if (iz<0 || iz>(fNz-2)) return TGeoShape::Big();
   return DistToSegZ(point,dir,iz);
}

////////////////////////////////////////////////////////////////////////////////
/// compute distance from outside point to surface of the tube

Double_t TGeoPcon::DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
{
   if ((iact<3) && safe) {
      *safe = Safety(point, kFALSE);
      if ((iact==1) && (*safe>step)) return TGeoShape::Big();
      if (iact==0) return TGeoShape::Big();
   }
   // check if ray intersect outscribed cylinder
   if ((point[2]<fZ[0]) && (dir[2]<=0)) return TGeoShape::Big();
   if ((point[2]>fZ[fNz-1]) && (dir[2]>=0)) return TGeoShape::Big();
// Check if the bounding box is crossed within the requested distance
   Double_t sdist = TGeoBBox::DistFromOutside(point,dir, fDX, fDY, fDZ, fOrigin, step);
   if (sdist>=step) return TGeoShape::Big();

   Double_t r2 = point[0]*point[0]+point[1]*point[1];
   Double_t radmax=0;
   radmax=fRmax[TMath::LocMax(fNz, fRmax)];
   if (r2>(radmax*radmax)) {
      Double_t rpr=-point[0]*dir[0]-point[1]*dir[1];
      Double_t nxy=dir[0]*dir[0]+dir[1]*dir[1];
      if (rpr<TMath::Sqrt((r2-radmax*radmax)*nxy)) return TGeoShape::Big();
   }

   // find in which Z segment we are
   Int_t ipl = TMath::BinarySearch(fNz, fZ, point[2]);
   Int_t ifirst = ipl;
   if (ifirst<0) {
      ifirst=0;
   } else if (ifirst>=(fNz-1)) ifirst=fNz-2;
   // find if point is in the phi gap
   Double_t phi=0;
   if (!fFullPhi) {
      phi=TMath::ATan2(point[1], point[0]);
      if (phi<0) phi+=2.*TMath::Pi();
   }

   // compute distance to boundary
   return DistToSegZ(point,dir,ifirst);
}

////////////////////////////////////////////////////////////////////////////////
/// Defines z position of a section plane, rmin and rmax at this z. Sections
/// should be defined in increasing or decreasing Z order and the last section
/// HAS to be snum = fNz-1

void TGeoPcon::DefineSection(Int_t snum, Double_t z, Double_t rmin, Double_t rmax)
{
   if ((snum<0) || (snum>=fNz)) return;
   fZ[snum]    = z;
   fRmin[snum] = rmin;
   fRmax[snum] = rmax;
   if (rmin>rmax)
      Warning("DefineSection", "Shape %s: invalid rmin=%g rmax=%g", GetName(), rmin, rmax);
   if (snum==(fNz-1)) {
      // Reorder sections in increasing Z order
      if (fZ[0] > fZ[snum]) {
         Int_t iz = 0;
         Int_t izi = fNz-1;
         Double_t temp;
         while (iz<izi) {
            temp = fZ[iz];
            fZ[iz] = fZ[izi];
            fZ[izi] = temp;
            temp = fRmin[iz];
            fRmin[iz] = fRmin[izi];
            fRmin[izi] = temp;
            temp = fRmax[iz];
            fRmax[iz] = fRmax[izi];
            fRmax[izi] = temp;
            iz++;
            izi--;
         }
      }
      ComputeBBox();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns number of segments on each mesh circle segment.

Int_t TGeoPcon::GetNsegments() const
{
   return gGeoManager->GetNsegments();
}

////////////////////////////////////////////////////////////////////////////////
/// Divide this polycone shape belonging to volume "voldiv" into ndiv volumes
/// called divname, from start position with the given step. Returns pointer
/// to created division cell volume in case of Z divisions. Z divisions can be
/// performed if the divided range is in between two consecutive Z planes.
/// In case a wrong division axis is supplied, returns pointer to
/// volume that was divided.

TGeoVolume *TGeoPcon::Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv,
                             Double_t start, Double_t step)
{
   TGeoShape *shape;           //--- shape to be created
   TGeoVolume *vol;            //--- division volume to be created
   TGeoVolumeMulti *vmulti;    //--- generic divided volume
   TGeoPatternFinder *finder;  //--- finder to be attached
   TString opt = "";           //--- option to be attached
   Double_t zmin = start;
   Double_t zmax = start+ndiv*step;
   Int_t isect = -1;
   Int_t is, id, ipl;
   switch (iaxis) {
      case 1:  //---               R division
         Error("Divide", "Shape %s: cannot divide a pcon on radius", GetName());
         return 0;
      case 2:  //---               Phi division
         finder = new TGeoPatternCylPhi(voldiv, ndiv, start, start+ndiv*step);
         vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
         voldiv->SetFinder(finder);
         finder->SetDivIndex(voldiv->GetNdaughters());
         shape = new TGeoPcon(-step/2, step, fNz);
         for (is=0; is<fNz; is++)
            ((TGeoPcon*)shape)->DefineSection(is, fZ[is], fRmin[is], fRmax[is]);
         vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
         vmulti->AddVolume(vol);
         opt = "Phi";
         for (id=0; id<ndiv; id++) {
            voldiv->AddNodeOffset(vol, id, start+id*step+step/2, opt.Data());
            ((TGeoNodeOffset*)voldiv->GetNodes()->At(voldiv->GetNdaughters()-1))->SetFinder(finder);
         }
         return vmulti;
      case 3: //---                Z division
         // find start plane
         for (ipl=0; ipl<fNz-1; ipl++) {
            if (start<fZ[ipl]) continue;
            else {
               if ((start+ndiv*step)>fZ[ipl+1]) continue;
            }
            isect = ipl;
            zmin = fZ[isect];
            zmax= fZ[isect+1];
            break;
         }
         if (isect<0) {
            Error("Divide", "Shape %s: cannot divide pcon on Z if divided region is not between 2 planes", GetName());
            return 0;
         }
         finder = new TGeoPatternZ(voldiv, ndiv, start, start+ndiv*step);
         vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
         voldiv->SetFinder(finder);
         finder->SetDivIndex(voldiv->GetNdaughters());
         opt = "Z";
         for (id=0; id<ndiv; id++) {
            Double_t z1 = start+id*step;
            Double_t z2 = start+(id+1)*step;
            Double_t rmin1 = (fRmin[isect]*(zmax-z1)-fRmin[isect+1]*(zmin-z1))/(zmax-zmin);
            Double_t rmax1 = (fRmax[isect]*(zmax-z1)-fRmax[isect+1]*(zmin-z1))/(zmax-zmin);
            Double_t rmin2 = (fRmin[isect]*(zmax-z2)-fRmin[isect+1]*(zmin-z2))/(zmax-zmin);
            Double_t rmax2 = (fRmax[isect]*(zmax-z2)-fRmax[isect+1]*(zmin-z2))/(zmax-zmin);
            Bool_t is_tube = (TGeoShape::IsSameWithinTolerance(fRmin[isect],fRmin[isect+1]) && TGeoShape::IsSameWithinTolerance(fRmax[isect],fRmax[isect+1]))?kTRUE:kFALSE;
            Bool_t is_seg = (fDphi<360)?kTRUE:kFALSE;
            if (is_seg) {
               if (is_tube) shape=new TGeoTubeSeg(fRmin[isect],fRmax[isect],step/2, fPhi1, fPhi1+fDphi);
               else shape=new TGeoConeSeg(step/2, rmin1, rmax1, rmin2, rmax2, fPhi1, fPhi1+fDphi);
            } else {
               if (is_tube) shape=new TGeoTube(fRmin[isect],fRmax[isect],step/2);
               else shape = new TGeoCone(step/2,rmin1,rmax1,rmin2,rmax2);
            }
            vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
            vmulti->AddVolume(vol);
            voldiv->AddNodeOffset(vol, id, start+id*step+step/2, opt.Data());
            ((TGeoNodeOffset*)voldiv->GetNodes()->At(voldiv->GetNdaughters()-1))->SetFinder(finder);
         }
         return vmulti;
      default:
         Error("Divide", "Shape %s: Wrong axis %d for division",GetName(), iaxis);
         return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns name of axis IAXIS.

const char *TGeoPcon::GetAxisName(Int_t iaxis) const
{
   switch (iaxis) {
      case 1:
         return "R";
      case 2:
         return "PHI";
      case 3:
         return "Z";
      default:
         return "UNDEFINED";
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get range of shape for a given axis.

Double_t TGeoPcon::GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const
{
   xlo = 0;
   xhi = 0;
   Double_t dx = 0;
   switch (iaxis) {
      case 2:
         xlo = fPhi1;
         xhi = fPhi1 + fDphi;
         dx = fDphi;
         return dx;
      case 3:
         xlo = fZ[0];
         xhi = fZ[fNz-1];
         dx = xhi-xlo;
         return dx;
   }
   return dx;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill vector param[4] with the bounding cylinder parameters. The order
/// is the following : Rmin, Rmax, Phi1, Phi2

void TGeoPcon::GetBoundingCylinder(Double_t *param) const
{
   param[0] = fRmin[0];           // Rmin
   param[1] = fRmax[0];           // Rmax
   for (Int_t i=1; i<fNz; i++) {
      if (fRmin[i] < param[0]) param[0] = fRmin[i];
      if (fRmax[i] > param[1]) param[1] = fRmax[i];
   }
   param[0] *= param[0];
   param[1] *= param[1];
   if (TGeoShape::IsSameWithinTolerance(fDphi,360.)) {
      param[2] = 0.;
      param[3] = 360.;
      return;
   }
   param[2] = (fPhi1<0)?(fPhi1+360.):fPhi1;     // Phi1
   param[3] = param[2]+fDphi;                   // Phi2
}

////////////////////////////////////////////////////////////////////////////////
/// Returns Rmin for Z segment IPL.

Double_t TGeoPcon::GetRmin(Int_t ipl) const
{
   if (ipl<0 || ipl>(fNz-1)) {
      Error("GetRmin","ipl=%i out of range (0,%i) in shape %s",ipl,fNz-1,GetName());
      return 0.;
   }
   return fRmin[ipl];
}

////////////////////////////////////////////////////////////////////////////////
/// Returns Rmax for Z segment IPL.

Double_t TGeoPcon::GetRmax(Int_t ipl) const
{
   if (ipl<0 || ipl>(fNz-1)) {
      Error("GetRmax","ipl=%i out of range (0,%i) in shape %s",ipl,fNz-1,GetName());
      return 0.;
   }
   return fRmax[ipl];
}

////////////////////////////////////////////////////////////////////////////////
/// Returns Z for segment IPL.

Double_t TGeoPcon::GetZ(Int_t ipl) const
{
   if (ipl<0 || ipl>(fNz-1)) {
      Error("GetZ","ipl=%i out of range (0,%i) in shape %s",ipl,fNz-1,GetName());
      return 0.;
   }
   return fZ[ipl];
}

////////////////////////////////////////////////////////////////////////////////
/// print shape parameters

void TGeoPcon::InspectShape() const
{
   printf("*** Shape %s: TGeoPcon ***\n", GetName());
   printf("    Nz    = %i\n", fNz);
   printf("    phi1  = %11.5f\n", fPhi1);
   printf("    dphi  = %11.5f\n", fDphi);
   for (Int_t ipl=0; ipl<fNz; ipl++)
      printf("     plane %i: z=%11.5f Rmin=%11.5f Rmax=%11.5f\n", ipl, fZ[ipl], fRmin[ipl], fRmax[ipl]);
   printf(" Bounding box:\n");
   TGeoBBox::InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a TBuffer3D describing *this* shape.
/// Coordinates are in local reference frame.

TBuffer3D *TGeoPcon::MakeBuffer3D() const
{
   Int_t nbPnts, nbSegs, nbPols;
   GetMeshNumbers(nbPnts, nbSegs, nbPols);
   if (nbPnts <= 0) return nullptr;

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

void TGeoPcon::SetSegsAndPols(TBuffer3D &buff) const
{
   if (!HasInsideSurface()) {
      SetSegsAndPolsNoInside(buff);
      return;
   }

   Int_t i, j;
   const Int_t n = gGeoManager->GetNsegments()+1;
   Int_t nz = GetNz();
   if (nz < 2) return;
   Int_t nbPnts = nz*2*n;
   if (nbPnts <= 0) return;
   Double_t dphi = GetDphi();

   Bool_t specialCase = TGeoShape::IsSameWithinTolerance(dphi,360);
   Int_t c = GetBasicColor();

   Int_t indx, indx2, k;
   indx = indx2 = 0;

   //inside & outside circles, number of segments: 2*nz*(n-1)
   //             special case number of segments: 2*nz*n
   for (i = 0; i < nz*2; i++) {
      indx2 = i*n;
      for (j = 1; j < n; j++) {
         buff.fSegs[indx++] = c;
         buff.fSegs[indx++] = indx2+j-1;
         buff.fSegs[indx++] = indx2+j;
      }
      if (specialCase) {
         buff.fSegs[indx++] = c;
         buff.fSegs[indx++] = indx2+j-1;
         buff.fSegs[indx++] = indx2;
      }
   }

   //bottom & top lines, number of segments: 2*n
   for (i = 0; i < 2; i++) {
      indx2 = i*(nz-1)*2*n;
      for (j = 0; j < n; j++) {
         buff.fSegs[indx++] = c;
         buff.fSegs[indx++] = indx2+j;
         buff.fSegs[indx++] = indx2+n+j;
      }
   }

   //inside & outside cilindres, number of segments: 2*(nz-1)*n
   for (i = 0; i < (nz-1); i++) {
      //inside cilinder
      indx2 = i*n*2;
      for (j = 0; j < n; j++) {
         buff.fSegs[indx++] = c+2;
         buff.fSegs[indx++] = indx2+j;
         buff.fSegs[indx++] = indx2+n*2+j;
      }
      //outside cilinder
      indx2 = i*n*2+n;
      for (j = 0; j < n; j++) {
         buff.fSegs[indx++] = c+3;
         buff.fSegs[indx++] = indx2+j;
         buff.fSegs[indx++] = indx2+n*2+j;
      }
   }

   //left & right sections, number of segments: 2*(nz-2)
   //          special case number of segments: 0
   if (!specialCase) {
      for (i = 1; i < (nz-1); i++) {
         for (j = 0; j < 2; j++) {
            buff.fSegs[indx++] = c;
            buff.fSegs[indx++] =  2*i    * n + j*(n-1);
            buff.fSegs[indx++] = (2*i+1) * n + j*(n-1);
         }
      }
   }

   Int_t m = n - 1 + (specialCase ?  1 : 0);
   indx = 0;

   //bottom & top, number of polygons: 2*(n-1)
   // special case number of polygons: 2*n
   for (j = 0; j < n-1; j++) {
      buff.fPols[indx++] = c+3;
      buff.fPols[indx++] = 4;
      buff.fPols[indx++] = 2*nz*m+j;
      buff.fPols[indx++] = m+j;
      buff.fPols[indx++] = 2*nz*m+j+1;
      buff.fPols[indx++] = j;
   }
   for (j = 0; j < n-1; j++) {
      buff.fPols[indx++] = c+3;
      buff.fPols[indx++] = 4;
      buff.fPols[indx++] = 2*nz*m+n+j;
      buff.fPols[indx++] = (nz*2-2)*m+j;
      buff.fPols[indx++] = 2*nz*m+n+j+1;
      buff.fPols[indx++] = (nz*2-2)*m+m+j;
   }
   if (specialCase) {
      buff.fPols[indx++] = c+3;
      buff.fPols[indx++] = 4;
      buff.fPols[indx++] = 2*nz*m+j;
      buff.fPols[indx++] = m+j;
      buff.fPols[indx++] = 2*nz*m;
      buff.fPols[indx++] = j;

      buff.fPols[indx++] = c+3;
      buff.fPols[indx++] = 4;
      buff.fPols[indx++] = 2*nz*m+n+j;
      buff.fPols[indx++] = (nz*2-2)*m+m+j;
      buff.fPols[indx++] = 2*nz*m+n;
      buff.fPols[indx++] = (nz*2-2)*m+j;
   }

   //inside & outside, number of polygons: (nz-1)*2*(n-1)
   for (k = 0; k < (nz-1); k++) {
      for (j = 0; j < n-1; j++) {
         buff.fPols[indx++] = c;
         buff.fPols[indx++] = 4;
         buff.fPols[indx++] = 2*k*m+j;
         buff.fPols[indx++] = nz*2*m+(2*k+2)*n+j+1;
         buff.fPols[indx++] = (2*k+2)*m+j;
         buff.fPols[indx++] = nz*2*m+(2*k+2)*n+j;
      }
      for (j = 0; j < n-1; j++) {
         buff.fPols[indx++] = c+1;
         buff.fPols[indx++] = 4;
         buff.fPols[indx++] = (2*k+1)*m+j;
         buff.fPols[indx++] = nz*2*m+(2*k+3)*n+j;
         buff.fPols[indx++] = (2*k+3)*m+j;
         buff.fPols[indx++] = nz*2*m+(2*k+3)*n+j+1;
      }
      if (specialCase) {
         buff.fPols[indx++] = c;
         buff.fPols[indx++] = 4;
         buff.fPols[indx++] = 2*k*m+j;
         buff.fPols[indx++] = nz*2*m+(2*k+2)*n;
         buff.fPols[indx++] = (2*k+2)*m+j;
         buff.fPols[indx++] = nz*2*m+(2*k+2)*n+j;

         buff.fPols[indx++] = c+1;
         buff.fPols[indx++] = 4;
         buff.fPols[indx++] = (2*k+1)*m+j;
         buff.fPols[indx++] = nz*2*m+(2*k+3)*n+j;
         buff.fPols[indx++] = (2*k+3)*m+j;
         buff.fPols[indx++] = nz*2*m+(2*k+3)*n;
      }
   }

   //left & right sections, number of polygons: 2*(nz-1)
   //          special case number of polygons: 0
   if (!specialCase) {
      indx2 = nz*2*(n-1);
      for (k = 0; k < (nz-1); k++) {
         buff.fPols[indx++] = c+2;
         buff.fPols[indx++] = 4;
         buff.fPols[indx++] = k==0 ? indx2 : indx2+2*nz*n+2*(k-1);
         buff.fPols[indx++] = indx2+2*(k+1)*n;
         buff.fPols[indx++] = indx2+2*nz*n+2*k;
         buff.fPols[indx++] = indx2+(2*k+3)*n;

         buff.fPols[indx++] = c+2;
         buff.fPols[indx++] = 4;
         buff.fPols[indx++] = k==0 ? indx2+n-1 : indx2+2*nz*n+2*(k-1)+1;
         buff.fPols[indx++] = indx2+(2*k+3)*n+n-1;
         buff.fPols[indx++] = indx2+2*nz*n+2*k+1;
         buff.fPols[indx++] = indx2+2*(k+1)*n+n-1;
      }
      buff.fPols[indx-8] = indx2+n;
      buff.fPols[indx-2] = indx2+2*n-1;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Fill TBuffer3D structure for segments and polygons, when no inner surface exists

void TGeoPcon::SetSegsAndPolsNoInside(TBuffer3D &buff) const
{
   const Int_t n = gGeoManager->GetNsegments() + 1;
   const Int_t nz = GetNz();
   const Int_t nbPnts = nz * n + 2;

   if ((nz < 2) || (nbPnts <= 0) || (n < 2)) return;

   Int_t c = GetBasicColor();

   Int_t indx = 0, indx1 = 0, indx2 = 0, i, j;

   //  outside circles, number of segments: nz*n
   for (i = 0; i < nz; i++) {
      indx2 = i * n;
      for (j = 1; j < n; j++) {
         buff.fSegs[indx++] = c;
         buff.fSegs[indx++] = indx2 + j - 1;
         buff.fSegs[indx++] = indx2 + j % (n-1);
      }
   }

   indx2 = 0;
   // bottom lines
   for (j = 0; j < n; j++) {
      buff.fSegs[indx++] = c;
      buff.fSegs[indx++] = indx2 + j % (n-1);
      buff.fSegs[indx++] = nbPnts - 2;
   }

   indx2 = (nz-1) * n;
   // top lines
   for (j = 0; j < n; j++) {
      buff.fSegs[indx++] = c;
      buff.fSegs[indx++] = indx2 + j % (n-1);
      buff.fSegs[indx++] = nbPnts - 1;
   }

   // outside cylinders, number of segments: (nz-1)*n
   for (i = 0; i < (nz - 1); i++) {
      // outside cylinder
      indx2 = i * n;
      for (j = 0; j < n; j++) {
         buff.fSegs[indx++] = c;
         buff.fSegs[indx++] = indx2 + j % (n-1);
         buff.fSegs[indx++] = indx2 + n + j % (n-1);
      }
   }

   indx = 0;

   // bottom cap
   indx1 = 0; // start of first z layer
   indx2 = nz*(n-1);
   for (j = 0; j < n - 1; j++) {
      buff.fPols[indx++] = c;
      buff.fPols[indx++] = 3;
      buff.fPols[indx++] = indx1 + j;
      buff.fPols[indx++] = indx2 + j + 1;
      buff.fPols[indx++] = indx2 + j;
   }

   // top cap
   indx1 = (nz-1)*(n-1); // start last z layer
   indx2 = nz*(n-1) + n;
   for (j = 0; j < n - 1; j++) {
      buff.fPols[indx++] = c;
      buff.fPols[indx++] = 3;
      buff.fPols[indx++] = indx1 + j; // last z layer
      buff.fPols[indx++] = indx2 + j;
      buff.fPols[indx++] = indx2 + j + 1;
   }

   // outside, number of polygons: (nz-1)*(n-1)
   for (Int_t k = 0; k < (nz - 1); k++) {
      indx1 = k*(n-1);
      indx2 = nz*(n-1) + n*2 + k*n;
      for (j = 0; j < n-1; j++) {
         buff.fPols[indx++] = c;
         buff.fPols[indx++] = 4;
         buff.fPols[indx++] = indx1 + j;
         buff.fPols[indx++] = indx2 + j;
         buff.fPols[indx++] = indx1 + j + (n-1);
         buff.fPols[indx++] = indx2 + j + 1;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute safety from POINT to segment between planes ipl, ipl+1 within safmin.

Double_t TGeoPcon::SafetyToSegment(const Double_t *point, Int_t ipl, Bool_t in, Double_t safmin) const
{
   Double_t safe = TGeoShape::Big();
   if (ipl<0 || ipl>fNz-2) return (safmin+1.); // error in input plane
// Get info about segment.
   Double_t dz = 0.5*(fZ[ipl+1]-fZ[ipl]);
   if (dz<1E-9) return 1E9; // radius-changing segment
   Double_t ptnew[3];
   memcpy(ptnew, point, 3*sizeof(Double_t));
   ptnew[2] -= 0.5*(fZ[ipl]+fZ[ipl+1]);
   safe = TMath::Abs(ptnew[2])-dz;
   if (safe>safmin) return TGeoShape::Big(); // means: stop checking further segments
   Double_t rmin1 = fRmin[ipl];
   Double_t rmax1 = fRmax[ipl];
   Double_t rmin2 = fRmin[ipl+1];
   Double_t rmax2 = fRmax[ipl+1];
   Bool_t   is_tube = (TGeoShape::IsSameWithinTolerance(rmin1,rmin2) && TGeoShape::IsSameWithinTolerance(rmax1,rmax2))?kTRUE:kFALSE;
   if (!fFullPhi) {
      if (is_tube) safe = TGeoTubeSeg::SafetyS(ptnew,in,rmin1,rmax1, dz,fPhi1,fPhi1+fDphi,0);
      else         safe = TGeoConeSeg::SafetyS(ptnew,in,dz,rmin1,rmax1,rmin2,rmax2,fPhi1,fPhi1+fDphi,0);
   } else {
      if (is_tube) safe = TGeoTube::SafetyS(ptnew,in,rmin1,rmax1,dz,0);
      else         safe = TGeoCone::SafetyS(ptnew,in,dz,rmin1,rmax1,rmin2,rmax2,0);
   }
   if (safe<0) safe=0;
   return safe;
}

////////////////////////////////////////////////////////////////////////////////
/// computes the closest distance from given point to this shape, according
/// to option. The matching point on the shape is stored in spoint.
/// localize the Z segment

Double_t TGeoPcon::Safety(const Double_t *point, Bool_t in) const
{
   Double_t safmin, saftmp;
   Double_t dz;
   Int_t ipl, iplane;

   if (in) {
   //---> point is inside pcon
      ipl = TMath::BinarySearch(fNz, fZ, point[2]);
      if (ipl==(fNz-1)) return 0;   // point on last Z boundary
      if (ipl<0) return 0;          // point on first Z boundary
      if (ipl>0 && TGeoShape::IsSameWithinTolerance(fZ[ipl-1],fZ[ipl]) && TGeoShape::IsSameWithinTolerance(point[2],fZ[ipl-1])) ipl--;
      dz = 0.5*(fZ[ipl+1]-fZ[ipl]);
      if (dz<1E-8) {
         // Point on a segment-changing plane
         safmin = TMath::Min(point[2]-fZ[ipl-1],fZ[ipl+2]-point[2]);
         saftmp = TGeoShape::Big();
         if (fDphi<360) saftmp = TGeoShape::SafetyPhi(point,in,fPhi1,fPhi1+fDphi);
         if (saftmp<safmin) safmin = saftmp;
         Double_t radius = TMath::Sqrt(point[0]*point[0]+point[1]*point[1]);
         if (fRmin[ipl]>0) safmin = TMath::Min(safmin, radius-fRmin[ipl]);
         if (fRmin[ipl+1]>0) safmin = TMath::Min(safmin, radius-fRmin[ipl+1]);
         safmin = TMath::Min(safmin, fRmax[ipl]-radius);
         safmin = TMath::Min(safmin, fRmax[ipl+1]-radius);
         if (safmin<0) safmin = 0;
         return safmin;
      }
      // Check safety for current segment
      safmin = SafetyToSegment(point, ipl);
      if (safmin>1E10) {
         //  something went wrong - point is not inside current segment
         return 0.;
      }
      if (safmin<1E-6) return TMath::Abs(safmin); // point on radius-changing plane
      // check increasing iplanes
      iplane = ipl+1;
      saftmp = 0.;
/*
      while ((iplane<fNz-1) && saftmp<1E10) {
         saftmp = TMath::Abs(SafetyToSegment(point,iplane,kFALSE,safmin));
         if (saftmp<safmin) safmin=saftmp;
         iplane++;
      }
      // now decreasing nplanes
      iplane = ipl-1;
      saftmp = 0.;
      while ((iplane>=0) && saftmp<1E10) {
         saftmp = TMath::Abs(SafetyToSegment(point,iplane,kFALSE,safmin));
         if (saftmp<safmin) safmin=saftmp;
         iplane--;
      }
*/
      return safmin;
   }
   //---> point is outside pcon
   ipl = TMath::BinarySearch(fNz, fZ, point[2]);
   if (ipl<0) ipl=0;
   else if (ipl==fNz-1) ipl=fNz-2;
   dz = 0.5*(fZ[ipl+1]-fZ[ipl]);
   if (dz<1E-8 && (ipl+2<fNz)) {
      ipl++;
      dz = 0.5*(fZ[ipl+1]-fZ[ipl]);
   }
   // Check safety for current segment
   safmin = SafetyToSegment(point, ipl, kFALSE);
   if (safmin<1E-6) return TMath::Abs(safmin); // point on radius-changing plane
   saftmp = 0.;
   // check increasing iplanes
   iplane = ipl+1;
   saftmp = 0.;
   while ((iplane<fNz-1) && saftmp<1E10) {
      saftmp = TMath::Abs(SafetyToSegment(point,iplane,kFALSE,safmin));
      if (saftmp<safmin) safmin=saftmp;
      iplane++;
   }
   // now decreasing nplanes
   iplane = ipl-1;
   saftmp = 0.;
   while ((iplane>=0) && saftmp<1E10) {
      saftmp = TMath::Abs(SafetyToSegment(point,iplane,kFALSE,safmin));
      if (saftmp<safmin) safmin=saftmp;
      iplane--;
   }
   return safmin;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoPcon::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive)) return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   phi1  = " << fPhi1 << ";" << std::endl;
   out << "   dphi  = " << fDphi << ";" << std::endl;
   out << "   nz    = " << fNz << ";" << std::endl;
   out << "   TGeoPcon *pcon = new TGeoPcon(\"" << GetName() << "\",phi1,dphi,nz);" << std::endl;
   for (Int_t i=0; i<fNz; i++) {
      out << "      z     = " << fZ[i] << ";" << std::endl;
      out << "      rmin  = " << fRmin[i] << ";" << std::endl;
      out << "      rmax  = " << fRmax[i] << ";" << std::endl;
      out << "   pcon->DefineSection(" << i << ", z,rmin,rmax);" << std::endl;
   }
   out << "   TGeoShape *" << GetPointerName() << " = pcon;" << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Set polycone dimensions starting from an array.

void TGeoPcon::SetDimensions(Double_t *param)
{
   fPhi1    = param[0];
   while (fPhi1<0) fPhi1 += 360.;
   fDphi    = param[1];
   fNz      = (Int_t)param[2];
   if (fNz<2) {
      Error("SetDimensions","Pcon %s: Number of Z sections must be > 2", GetName());
      return;
   }
   if (fRmin) delete [] fRmin;
   if (fRmax) delete [] fRmax;
   if (fZ) delete [] fZ;
   fRmin = new Double_t [fNz];
   fRmax = new Double_t [fNz];
   fZ    = new Double_t [fNz];
   memset(fRmin, 0, fNz*sizeof(Double_t));
   memset(fRmax, 0, fNz*sizeof(Double_t));
   memset(fZ, 0, fNz*sizeof(Double_t));
   if (TGeoShape::IsSameWithinTolerance(fDphi,360)) fFullPhi = kTRUE;
   Double_t phi1 = fPhi1;
   Double_t phi2 = phi1+fDphi;
   Double_t phim = 0.5*(phi1+phi2);
   fC1 = TMath::Cos(phi1*TMath::DegToRad());
   fS1 = TMath::Sin(phi1*TMath::DegToRad());
   fC2 = TMath::Cos(phi2*TMath::DegToRad());
   fS2 = TMath::Sin(phi2*TMath::DegToRad());
   fCm = TMath::Cos(phim*TMath::DegToRad());
   fSm = TMath::Sin(phim*TMath::DegToRad());
   fCdphi = TMath::Cos(0.5*fDphi*TMath::DegToRad());

   for (Int_t i=0; i<fNz; i++)
      DefineSection(i, param[3+3*i], param[4+3*i], param[5+3*i]);
}

////////////////////////////////////////////////////////////////////////////////
/// create polycone mesh points

void TGeoPcon::SetPoints(Double_t *points) const
{
   Double_t phi, dphi;
   Int_t n = gGeoManager->GetNsegments() + 1;
   dphi = fDphi/(n-1);
   Int_t i, j;
   Int_t indx = 0;

   Bool_t hasInside = HasInsideSurface();

   if (points) {
      for (i = 0; i < fNz; i++) {
         if (hasInside)
            for (j = 0; j < n; j++) {
               phi = (fPhi1+j*dphi)*TMath::DegToRad();
               points[indx++] = fRmin[i] * TMath::Cos(phi);
               points[indx++] = fRmin[i] * TMath::Sin(phi);
               points[indx++] = fZ[i];
            }
         for (j = 0; j < n; j++) {
            phi = (fPhi1+j*dphi)*TMath::DegToRad();
            points[indx++] = fRmax[i] * TMath::Cos(phi);
            points[indx++] = fRmax[i] * TMath::Sin(phi);
            points[indx++] = fZ[i];
         }
      }
      if (!hasInside) {
         points[indx++] = 0;
         points[indx++] = 0;
         points[indx++] = fZ[0];

         points[indx++] = 0;
         points[indx++] = 0;
         points[indx++] = fZ[GetNz()-1];
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// create polycone mesh points

void TGeoPcon::SetPoints(Float_t *points) const
{
   Double_t phi, dphi;
   Int_t n = gGeoManager->GetNsegments() + 1;
   dphi = fDphi/(n-1);
   Int_t i, j;
   Int_t indx = 0;

   Bool_t hasInside = HasInsideSurface();

   if (points) {
      for (i = 0; i < fNz; i++) {
         if (hasInside)
            for (j = 0; j < n; j++) {
               phi = (fPhi1+j*dphi)*TMath::DegToRad();
               points[indx++] = fRmin[i] * TMath::Cos(phi);
               points[indx++] = fRmin[i] * TMath::Sin(phi);
               points[indx++] = fZ[i];
            }
         for (j = 0; j < n; j++) {
            phi = (fPhi1+j*dphi)*TMath::DegToRad();
            points[indx++] = fRmax[i] * TMath::Cos(phi);
            points[indx++] = fRmax[i] * TMath::Sin(phi);
            points[indx++] = fZ[i];
         }
      }
      if (!hasInside) {
         points[indx++] = 0;
         points[indx++] = 0;
         points[indx++] = fZ[0];

         points[indx++] = 0;
         points[indx++] = 0;
         points[indx++] = fZ[GetNz()-1];
      }
   }
}
////////////////////////////////////////////////////////////////////////////////
/// Return number of vertices of the mesh representation

Int_t TGeoPcon::GetNmeshVertices() const
{
   Int_t nvert, nsegs, npols;
   GetMeshNumbers(nvert, nsegs, npols);
   return nvert;
}

////////////////////////////////////////////////////////////////////////////////
/// fill size of this 3-D object

void TGeoPcon::Sizeof3D() const
{
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true when pgon has internal surface
/// It will be only disabled when all Rmin values are 0

Bool_t TGeoPcon::HasInsideSurface() const
{
   // only when full 360 is used, internal part can be excluded
   Bool_t specialCase = TGeoShape::IsSameWithinTolerance(GetDphi(), 360);
   if (!specialCase) return kTRUE;

   for (Int_t i = 0; i < GetNz(); i++)
      if (fRmin[i] > 0.)
         return kTRUE;

   return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoPcon::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   nvert = nsegs = npols = 0;

   Int_t n = gGeoManager->GetNsegments()+1;
   Int_t nz = GetNz();
   if (nz < 2) return;

   if (HasInsideSurface()) {
      Bool_t specialCase = TGeoShape::IsSameWithinTolerance(GetDphi(),360);
      nvert = nz*2*n;
      nsegs = 4*(nz*n-1+(specialCase ?  1 : 0));
      npols = 2*(nz*n-1+(specialCase ?  1 : 0));
   } else {
      nvert = nz * n + 2;
      nsegs = nz * (n - 1) + n * 2 + (nz - 1) * n;
      npols = 2 * (n - 1) + (nz - 1) * (n - 1);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fills a static 3D buffer and returns a reference.

const TBuffer3D & TGeoPcon::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
{
   static TBuffer3D buffer(TBuffer3DTypes::kGeneric);

   TGeoBBox::FillBuffer3D(buffer, reqSections, localFrame);

   if (reqSections & TBuffer3D::kRawSizes) {
      Int_t nbPnts, nbSegs, nbPols;
      GetMeshNumbers(nbPnts, nbSegs, nbPols);
      if (nbPnts > 0) {
         if (buffer.SetRawSizes(nbPnts, 3*nbPnts, nbSegs, 3*nbSegs, nbPols, 6*nbPols)) {
            buffer.SetSectionsValid(TBuffer3D::kRawSizes);
         }
      }
   }
   // TODO: Push down to TGeoShape?? Would have to do raw sizes set first..
   // can rest of TGeoShape be deferred until after this?
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
/// Stream an object of class TGeoPcon.

void TGeoPcon::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TGeoPcon::Class(),this);
      if (TGeoShape::IsSameWithinTolerance(fDphi,360)) fFullPhi = kTRUE;
      Double_t phi1 = fPhi1;
      Double_t phi2 = phi1+fDphi;
      Double_t phim = 0.5*(phi1+phi2);
      fC1 = TMath::Cos(phi1*TMath::DegToRad());
      fS1 = TMath::Sin(phi1*TMath::DegToRad());
      fC2 = TMath::Cos(phi2*TMath::DegToRad());
      fS2 = TMath::Sin(phi2*TMath::DegToRad());
      fCm = TMath::Cos(phim*TMath::DegToRad());
      fSm = TMath::Sin(phim*TMath::DegToRad());
      fCdphi = TMath::Cos(0.5*fDphi*TMath::DegToRad());
   } else {
      R__b.WriteClassBuffer(TGeoPcon::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check the inside status for each of the points in the array.
/// Input: Array of point coordinates + vector size
/// Output: Array of Booleans for the inside of each point

void TGeoPcon::Contains_v(const Double_t *points, Bool_t *inside, Int_t vecsize) const
{
   for (Int_t i=0; i<vecsize; i++) inside[i] = Contains(&points[3*i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the normal for an array o points so that norm.dot.dir is positive
/// Input: Arrays of point coordinates and directions + vector size
/// Output: Array of normal directions

void TGeoPcon::ComputeNormal_v(const Double_t *points, const Double_t *dirs, Double_t *norms, Int_t vecsize)
{
   for (Int_t i=0; i<vecsize; i++) ComputeNormal(&points[3*i], &dirs[3*i], &norms[3*i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoPcon::DistFromInside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromInside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoPcon::DistFromOutside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromOutside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute safe distance from each of the points in the input array.
/// Input: Array of point coordinates, array of statuses for these points, size of the arrays
/// Output: Safety values

void TGeoPcon::Safety_v(const Double_t *points, const Bool_t *inside, Double_t *safe, Int_t vecsize) const
{
   for (Int_t i=0; i<vecsize; i++) safe[i] = Safety(&points[3*i], inside[i]);
}
