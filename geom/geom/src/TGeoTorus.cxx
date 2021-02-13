// @(#)root/geom:$Id$
// Author: Andrei Gheata   28/07/03

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoTorus
\ingroup Geometry_classes

Torus segment class. A torus has 5 parameters :
  - R    - axial radius
  - Rmin - inner radius
  - Rmax - outer radius
  - Phi1 - starting phi
  - Dphi - phi extent


Begin_Macro(source)
{
   TCanvas *c = new TCanvas("c", "c",0,0,600,600);
   new TGeoManager("torus", "poza2");
   TGeoMaterial *mat = new TGeoMaterial("Al", 26.98,13,2.7);
   TGeoMedium *med = new TGeoMedium("MED",1,mat);
   TGeoVolume *top = gGeoManager->MakeBox("TOP",med,100,100,100);
   gGeoManager->SetTopVolume(top);
   TGeoVolume *vol = gGeoManager->MakeTorus("TORUS",med, 40,20,25,0,270);
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

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoTube.h"
#include "TVirtualGeoPainter.h"
#include "TGeoTorus.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"

ClassImp(TGeoTorus);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoTorus::TGeoTorus()
{
   SetShapeBit(TGeoShape::kGeoTorus);
   fR    = 0.0;
   fRmin = 0.0;
   fRmax = 0.0;
   fPhi1 = 0.0;
   fDphi = 0.0;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor without name.

TGeoTorus::TGeoTorus(Double_t r, Double_t rmin, Double_t rmax, Double_t phi1, Double_t dphi)
          :TGeoBBox(0, 0, 0)
{
   SetShapeBit(TGeoShape::kGeoTorus);
   SetTorusDimensions(r, rmin, rmax, phi1, dphi);
   if ((fRmin<0) || (fRmax<0))
      SetShapeBit(kGeoRunTimeShape);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor with name.

TGeoTorus::TGeoTorus(const char *name, Double_t r, Double_t rmin, Double_t rmax, Double_t phi1, Double_t dphi)
          :TGeoBBox(name, 0, 0, 0)
{
   SetShapeBit(TGeoShape::kGeoTorus);
   SetTorusDimensions(r, rmin, rmax, phi1, dphi);
   if ((fRmin<0) || (fRmax<0))
      SetShapeBit(kGeoRunTimeShape);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor based on an array of parameters.
///  - param[0] = R
///  - param[1] = Rmin
///  - param[2] = Rmax
///  - param[3] = Phi1
///  - param[4] = Dphi

TGeoTorus::TGeoTorus(Double_t *param)
          :TGeoBBox(0, 0, 0)
{
   SetShapeBit(TGeoShape::kGeoTorus);
   SetDimensions(param);
   if (fRmin<0 || fRmax<0) SetShapeBit(kGeoRunTimeShape);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Computes capacity of the shape in [length^3]

Double_t TGeoTorus::Capacity() const
{
   Double_t capacity = (fDphi/180.)*TMath::Pi()*TMath::Pi()*fR*(fRmax*fRmax-fRmin*fRmin);
   return capacity;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute bounding box of the torus.

void TGeoTorus::ComputeBBox()
{
   fDZ = fRmax;
   if (TGeoShape::IsSameWithinTolerance(fDphi,360)) {
      fDX = fDY = fR+fRmax;
      return;
   }
   Double_t xc[4];
   Double_t yc[4];
   xc[0] = (fR+fRmax)*TMath::Cos(fPhi1*TMath::DegToRad());
   yc[0] = (fR+fRmax)*TMath::Sin(fPhi1*TMath::DegToRad());
   xc[1] = (fR+fRmax)*TMath::Cos((fPhi1+fDphi)*TMath::DegToRad());
   yc[1] = (fR+fRmax)*TMath::Sin((fPhi1+fDphi)*TMath::DegToRad());
   xc[2] = (fR-fRmax)*TMath::Cos(fPhi1*TMath::DegToRad());
   yc[2] = (fR-fRmax)*TMath::Sin(fPhi1*TMath::DegToRad());
   xc[3] = (fR-fRmax)*TMath::Cos((fPhi1+fDphi)*TMath::DegToRad());
   yc[3] = (fR-fRmax)*TMath::Sin((fPhi1+fDphi)*TMath::DegToRad());

   Double_t xmin = xc[TMath::LocMin(4, &xc[0])];
   Double_t xmax = xc[TMath::LocMax(4, &xc[0])];
   Double_t ymin = yc[TMath::LocMin(4, &yc[0])];
   Double_t ymax = yc[TMath::LocMax(4, &yc[0])];
   Double_t ddp = -fPhi1;
   if (ddp<0) ddp+= 360;
   if (ddp<=fDphi) xmax = fR+fRmax;
   ddp = 90-fPhi1;
   if (ddp<0) ddp+= 360;
   if (ddp>360) ddp-=360;
   if (ddp<=fDphi) ymax = fR+fRmax;
   ddp = 180-fPhi1;
   if (ddp<0) ddp+= 360;
   if (ddp>360) ddp-=360;
   if (ddp<=fDphi) xmin = -(fR+fRmax);
   ddp = 270-fPhi1;
   if (ddp<0) ddp+= 360;
   if (ddp>360) ddp-=360;
   if (ddp<=fDphi) ymin = -(fR+fRmax);
   fOrigin[0] = (xmax+xmin)/2;
   fOrigin[1] = (ymax+ymin)/2;
   fOrigin[2] = 0;
   fDX = (xmax-xmin)/2;
   fDY = (ymax-ymin)/2;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute normal to closest surface from POINT.

void TGeoTorus::ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm)
{
   Double_t phi = TMath::ATan2(point[1],point[0]);
   if (fDphi<360) {
      Double_t phi1 = fPhi1*TMath::DegToRad();
      Double_t phi2 = (fPhi1+fDphi)*TMath::DegToRad();
      Double_t c1 = TMath::Cos(phi1);
      Double_t s1 = TMath::Sin(phi1);
      Double_t c2 = TMath::Cos(phi2);
      Double_t s2 = TMath::Sin(phi2);

      Double_t daxis = Daxis(point,dir,0);
      if ((fRmax-daxis)>1E-5) {
         if (TGeoShape::IsSameWithinTolerance(fRmin,0) || (daxis-fRmin)>1E-5) {
            TGeoShape::NormalPhi(point,dir,norm,c1,s1,c2,s2);
            return;
         }
      }
   }
   Double_t r0[3];
   r0[0] = fR*TMath::Cos(phi);
   r0[1] = fR*TMath::Sin(phi);
   r0[2] = 0;
   Double_t normsq = 0;
   for (Int_t i=0; i<3; i++) {
      norm[i] = point[i] - r0[i];
      normsq += norm[i]*norm[i];
   }

   normsq = TMath::Sqrt(normsq);
   norm[0] /= normsq;
   norm[1] /= normsq;
   norm[2] /= normsq;
   if (dir[0]*norm[0]+dir[1]*norm[1]+dir[2]*norm[2] < 0) {
      norm[0] = -norm[0];
      norm[1] = -norm[1];
      norm[2] = -norm[2];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Test if point is inside the torus.
/// check phi range

Bool_t TGeoTorus::Contains(const Double_t *point) const
{
   if (!TGeoShape::IsSameWithinTolerance(fDphi,360)) {
      Double_t phi = TMath::ATan2(point[1], point[0]) * TMath::RadToDeg();
      if (phi < 0) phi+=360.0;
      Double_t ddp = phi-fPhi1;
      if (ddp<0) ddp+=360.;
      if (ddp>fDphi) return kFALSE;
   }
   //check radius
   Double_t rxy = TMath::Sqrt(point[0]*point[0]+point[1]*point[1]);
   Double_t radsq = (rxy-fR)*(rxy-fR) + point[2]*point[2];
   if (radsq<fRmin*fRmin) return kFALSE;
   if (radsq>fRmax*fRmax) return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute closest distance from point px,py to each vertex.

Int_t TGeoTorus::DistancetoPrimitive(Int_t px, Int_t py)
{
   Int_t n = gGeoManager->GetNsegments()+1;
   Int_t numPoints = n*(n-1);
   if (fRmin>0) numPoints *= 2;
   else if (fDphi<360) numPoints += 2;
   return ShapeDistancetoPrimitive(numPoints, px, py);
}

////////////////////////////////////////////////////////////////////////////////
/// Computes distance to axis of the torus from point pt + t*dir;

Double_t TGeoTorus::Daxis(const Double_t *pt, const Double_t *dir, Double_t t) const
{
   Double_t p[3];
   for (Int_t i=0; i<3; i++) p[i] = pt[i]+t*dir[i];
   Double_t rxy = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]);
   return TMath::Sqrt((rxy-fR)*(rxy-fR)+p[2]*p[2]);
}

////////////////////////////////////////////////////////////////////////////////
/// Computes derivative w.r.t. t of the distance to axis of the torus from point pt + t*dir;

Double_t TGeoTorus::DDaxis(const Double_t *pt, const Double_t *dir, Double_t t) const
{
   Double_t p[3];
   for (Int_t i=0; i<3; i++) p[i] = pt[i]+t*dir[i];
   Double_t rxy = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]);
   if (rxy<1E-4) return ((p[2]*dir[2]-fR*TMath::Sqrt(dir[0]*dir[0]+dir[1]*dir[1]))/TMath::Sqrt(fR*fR+p[2]*p[2]));
   Double_t d = TMath::Sqrt((rxy-fR)*(rxy-fR)+p[2]*p[2]);
   if (TGeoShape::IsSameWithinTolerance(d,0)) return 0.;
   Double_t dd = (p[0]*dir[0]+p[1]*dir[1]+p[2]*dir[2] - (p[0]*dir[0]+p[1]*dir[1])*fR/rxy)/d;
   return dd;
}

////////////////////////////////////////////////////////////////////////////////
/// Second derivative of distance to torus axis w.r.t t.

Double_t TGeoTorus::DDDaxis(const Double_t *pt, const Double_t *dir, Double_t t) const
{
   Double_t p[3];
   for (Int_t i=0; i<3; i++) p[i] = pt[i]+t*dir[i];
   Double_t rxy = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]);
   if (rxy<1E-6) return 0;
   Double_t daxis = TMath::Sqrt((rxy-fR)*(rxy-fR)+p[2]*p[2]);
   if (TGeoShape::IsSameWithinTolerance(daxis,0)) return 0;
   Double_t ddaxis = (p[0]*dir[0]+p[1]*dir[1]+p[2]*dir[2] - (p[0]*dir[0]+p[1]*dir[1])*fR/rxy)/daxis;
   Double_t dddaxis = 1 - ddaxis*ddaxis - (1-dir[2]*dir[2])*fR/rxy +
                      fR*(p[0]*dir[0]+p[1]*dir[1])*(p[0]*dir[0]+p[1]*dir[1])/(rxy*rxy*rxy);
   dddaxis /= daxis;
   return dddaxis;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from inside point to surface of the torus.

Double_t TGeoTorus::DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
{
   if (iact<3 && safe) {
      *safe = Safety(point, kTRUE);
      if (iact==0) return TGeoShape::Big();
      if ((iact==1) && (step<=*safe)) return TGeoShape::Big();
   }
   Double_t snext = TGeoShape::Big();
   Bool_t hasphi = (fDphi<360)?kTRUE:kFALSE;
   Bool_t hasrmin = (fRmin>0)?kTRUE:kFALSE;
   Double_t dout = ToBoundary(point,dir,fRmax,kTRUE);
//   Double_t dax = Daxis(point,dir,dout);
   Double_t din = (hasrmin)?ToBoundary(point,dir,fRmin,kTRUE):TGeoShape::Big();
   snext = TMath::Min(dout,din);
   if (snext>1E10) return TGeoShape::Tolerance();
   Double_t dphi = TGeoShape::Big();
   if (hasphi) {
      // Torus segment case.
      Double_t c1,s1,c2,s2,cm,sm,cdfi;
      Double_t phi1=fPhi1*TMath::DegToRad();
      Double_t phi2=(fPhi1+fDphi)*TMath::DegToRad();
      c1=TMath::Cos(phi1);
      s1=TMath::Sin(phi1);
      c2=TMath::Cos(phi2);
      s2=TMath::Sin(phi2);
      Double_t fio=0.5*(phi1+phi2);
      cm=TMath::Cos(fio);
      sm=TMath::Sin(fio);
      cdfi = TMath::Cos(0.5*(phi2-phi1));
      dphi = TGeoTubeSeg::DistFromInsideS(point,dir,fR-fRmax,fR+fRmax, fRmax, c1,s1,c2,s2,cm,sm,cdfi);
      Double_t daxis = Daxis(point,dir,dphi);
      if (daxis>=fRmin+1.E-8 && daxis<=fRmax-1.E-8) snext=TMath::Min(snext,dphi);
   }
   return snext;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from outside point to surface of the torus.

Double_t TGeoTorus::DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
{
   if (iact<3 && safe) {
      *safe = Safety(point, kFALSE);
      if (iact==0) return TGeoShape::Big();
      if ((iact==1) && (step<=*safe)) return TGeoShape::Big();
   }
// Check if the bounding box is crossed within the requested distance
   Double_t sdist = TGeoBBox::DistFromOutside(point,dir, fDX, fDY, fDZ, fOrigin, step);
   if (sdist>=step) return TGeoShape::Big();
   Double_t daxis;
   Bool_t hasphi = (fDphi<360)?kTRUE:kFALSE;
//   Bool_t hasrmin = (fRmin>0)?kTRUE:kFALSE;
   Double_t c1=0,s1=0,c2=0,s2=0,cm=0,sm=0,cdfi=0;
   Bool_t inphi = kFALSE;
   Double_t phi, ddp, phi1,phi2,fio;
   Double_t rxy2,dd;
   Double_t snext;
   Double_t pt[3];
   Int_t i;

   if (hasphi) {
      // Torus segment case.
      phi=TMath::ATan2(point[1], point[0])*TMath::RadToDeg();;
      if (phi<0) phi+=360;
      ddp = phi-fPhi1;
      if (ddp<0) ddp+=360;;
      if (ddp<=fDphi) inphi=kTRUE;
      phi1=fPhi1*TMath::DegToRad();
      phi2=(fPhi1+fDphi)*TMath::DegToRad();
      c1=TMath::Cos(phi1);
      s1=TMath::Sin(phi1);
      c2=TMath::Cos(phi2);
      s2=TMath::Sin(phi2);
      fio=0.5*(phi1+phi2);
      cm=TMath::Cos(fio);
      sm=TMath::Sin(fio);
      cdfi=TMath::Cos(0.5*(phi2-phi1));
   }
   // Check if we are inside or outside the bounding ring.
   Bool_t inbring = kFALSE;
   if (TMath::Abs(point[2]) <= fRmax) {
      rxy2 = point[0]*point[0]+point[1]*point[1];
      if ((rxy2>=(fR-fRmax)*(fR-fRmax)) && (rxy2<=(fR+fRmax)*(fR+fRmax))) {
         if (!hasphi || inphi) inbring=kTRUE;
      }
   }

   // If outside the ring, compute distance to it.
   Double_t dring = TGeoShape::Big();
   Double_t eps = 1.E-8;
   snext = 0;
   daxis = -1;
   memcpy(pt,point,3*sizeof(Double_t));
   if (!inbring) {
      if (hasphi) dring = TGeoTubeSeg::DistFromOutsideS(point,dir,TMath::Max(0.,fR-fRmax-eps),fR+fRmax+eps, fRmax+eps, c1,s1,c2,s2,cm,sm,cdfi);
      else        dring = TGeoTube::DistFromOutsideS(point,dir,TMath::Max(0.,fR-fRmax-eps),fR+fRmax+eps, fRmax+eps);
      // If not crossing it, return BIG.
      if (dring>1E10) return TGeoShape::Big();
      snext = dring;
      // Check if the crossing is due to phi.
      daxis = Daxis(point,dir,snext);
      if (daxis>=fRmin && daxis<fRmax) return snext;
      // Not a phi crossing -> propagate until we cross the ring.
      for (i=0; i<3; i++) pt[i] = point[i]+snext*dir[i];
   }
   // Point pt is inside the bounding ring, no phi crossing so far.
   // Check if we are in the hole.
   if (daxis<0) daxis = Daxis(pt,dir,0);
   if (daxis<fRmin+1.E-8) {
      // We are in the hole. Check if we came from outside.
      if (snext>0) {
         // we can cross either the inner torus or exit the other hole.
         snext += 0.1*eps;
         for (i=0; i<3; i++) pt[i] += 0.1*eps*dir[i];
      }
      // We are in the hole from the beginning.
      // find first crossing with inner torus
      dd = ToBoundary(pt,dir, fRmin,kFALSE);
      // find exit distance from inner bounding ring
      if (hasphi) dring = TGeoTubeSeg::DistFromInsideS(pt,dir,fR-fRmin,fR+fRmin, fRmin, c1,s1,c2,s2,cm,sm,cdfi);
      else        dring = TGeoTube::DistFromInsideS(pt,dir,fR-fRmin,fR+fRmin, fRmin);
      if (dd<dring) return (snext+dd);
      // we were exiting a hole inside phi hole
      snext += dring+ eps;
      for (i=0; i<3; i++) pt[i] = point[i] + snext*dir[i];
      snext += DistFromOutside(pt,dir,3);
      return snext;
   }
   // We are inside the outer ring, having daxis>fRmax
   // Compute distance to exit the bounding ring (again)
   if (snext>0) {
      // we can cross either the inner torus or exit the other hole.
      snext += 0.1*eps;
      for (i=0; i<3; i++) pt[i] += 0.1*eps*dir[i];
   }
   // Check intersection with outer torus
   dd = ToBoundary(pt, dir, fRmax, kFALSE);
   if (hasphi) dring = TGeoTubeSeg::DistFromInsideS(pt,dir,TMath::Max(0.,fR-fRmax-eps),fR+fRmax+eps, fRmax+eps, c1,s1,c2,s2,cm,sm,cdfi);
   else        dring = TGeoTube::DistFromInsideS(pt,dir,TMath::Max(0.,fR-fRmax-eps),fR+fRmax+eps, fRmax+eps);
   if (dd<dring) {
      snext += dd;
      return snext;
   }
   // We are exiting the bounding ring before crossing the torus -> propagate
   snext += dring+eps;
   for (i=0; i<3; i++) pt[i] = point[i] + snext*dir[i];
   snext += DistFromOutside(pt,dir,3);
   return snext;
}

////////////////////////////////////////////////////////////////////////////////
/// Divide this torus shape belonging to volume "voldiv" into ndiv volumes
/// called divname, from start position with the given step.

TGeoVolume *TGeoTorus::Divide(TGeoVolume * /*voldiv*/, const char * /*divname*/, Int_t /*iaxis*/, Int_t /*ndiv*/,
                              Double_t /*start*/, Double_t /*step*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns name of axis IAXIS.

const char *TGeoTorus::GetAxisName(Int_t iaxis) const
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

Double_t TGeoTorus::GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const
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
         xhi = fPhi1+fDphi;
         dx = fDphi;
         return dx;
      case 3:
         dx = 0;
         return dx;
   }
   return dx;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill vector param[4] with the bounding cylinder parameters. The order
/// is the following : Rmin, Rmax, Phi1, Phi2, dZ

void TGeoTorus::GetBoundingCylinder(Double_t *param) const
{
   param[0] = (fR-fRmax); // Rmin
   param[1] = (fR+fRmax); // Rmax
   param[2] = fPhi1;    // Phi1
   param[3] = fPhi1+fDphi;  // Phi2
}

////////////////////////////////////////////////////////////////////////////////
/// Create a shape fitting the mother.

TGeoShape *TGeoTorus::GetMakeRuntimeShape(TGeoShape * /*mother*/, TGeoMatrix * /*mat*/) const
{
   if (!TestShapeBit(kGeoRunTimeShape)) return 0;
   Error("GetMakeRuntimeShape", "parametrized toruses not supported");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// print shape parameters

void TGeoTorus::InspectShape() const
{
   printf("*** Shape %s: TGeoTorus ***\n", GetName());
   printf("    R    = %11.5f\n", fR);
   printf("    Rmin = %11.5f\n", fRmin);
   printf("    Rmax = %11.5f\n", fRmax);
   printf("    Phi1 = %11.5f\n", fPhi1);
   printf("    Dphi = %11.5f\n", fDphi);
   printf(" Bounding box:\n");
   TGeoBBox::InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a TBuffer3D describing *this* shape.
/// Coordinates are in local reference frame.

TBuffer3D *TGeoTorus::MakeBuffer3D() const
{
   Int_t n = gGeoManager->GetNsegments()+1;
   Int_t nbPnts = n*(n-1);
   Bool_t hasrmin = (GetRmin()>0)?kTRUE:kFALSE;
   Bool_t hasphi  = (GetDphi()<360)?kTRUE:kFALSE;
   if (hasrmin) nbPnts *= 2;
   else if (hasphi) nbPnts += 2;

   Int_t nbSegs = (2*n-1)*(n-1);
   Int_t nbPols = (n-1)*(n-1);
   if (hasrmin) {
      nbSegs += (2*n-1)*(n-1);
      nbPols += (n-1)*(n-1);
   }
   if (hasphi) {
      nbSegs += 2*(n-1);
      nbPols += 2*(n-1);
   }

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

void TGeoTorus::SetSegsAndPols(TBuffer3D &buff) const
{
   Int_t i, j;
   Int_t n = gGeoManager->GetNsegments()+1;
   Int_t nbPnts = n*(n-1);
   Int_t indx, indp, startcap=0;
   Bool_t hasrmin = (GetRmin()>0)?kTRUE:kFALSE;
   Bool_t hasphi  = (GetDphi()<360)?kTRUE:kFALSE;
   if (hasrmin) nbPnts *= 2;
   else if (hasphi) nbPnts += 2;
   Int_t c = GetBasicColor();

   indp = n*(n-1); // start index for points on inner surface
   memset(buff.fSegs, 0, buff.NbSegs()*3*sizeof(Int_t));

   // outer surface phi circles = n*(n-1) -> [0, n*(n-1) -1]
   // connect point j with point j+1 on same row
   indx = 0;
   for (i = 0; i < n; i++) { // rows [0,n-1]
      for (j = 0; j < n-1; j++) {  // points on a row [0, n-2]
         buff.fSegs[indx+(i*(n-1)+j)*3] = c;
         buff.fSegs[indx+(i*(n-1)+j)*3+1] = i*(n-1)+j;   // j on row i
         buff.fSegs[indx+(i*(n-1)+j)*3+2] = i*(n-1)+((j+1)%(n-1)); // j+1 on row i
      }
   }
   indx += 3*n*(n-1);
   // outer surface generators = (n-1)*(n-1) -> [n*(n-1), (2*n-1)*(n-1) -1]
   // connect point j on row i with point j on row i+1
   for (i = 0; i < n-1; i++) { // rows [0, n-2]
      for (j = 0; j < n-1; j++) {  // points on a row [0, n-2]
         buff.fSegs[indx+(i*(n-1)+j)*3] = c;
         buff.fSegs[indx+(i*(n-1)+j)*3+1] = i*(n-1)+j;     // j on row i
         buff.fSegs[indx+(i*(n-1)+j)*3+2] = (i+1)*(n-1)+j; // j on row i+1
      }
   }
   indx += 3*(n-1)*(n-1);
   startcap = (2*n-1)*(n-1);

   if (hasrmin) {
      // inner surface phi circles = n*(n-1) -> [(2*n-1)*(n-1), (3*n-1)*(n-1) -1]
      // connect point j with point j+1 on same row
      for (i = 0; i < n; i++) { // rows [0, n-1]
         for (j = 0; j < n-1; j++) {  // points on a row [0, n-2]
            buff.fSegs[indx+(i*(n-1)+j)*3] = c;              // lighter color
            buff.fSegs[indx+(i*(n-1)+j)*3+1] = indp + i*(n-1)+j;   // j on row i
            buff.fSegs[indx+(i*(n-1)+j)*3+2] = indp + i*(n-1)+((j+1)%(n-1)); // j+1 on row i
         }
      }
      indx += 3*n*(n-1);
      // inner surface generators = (n-1)*n -> [(3*n-1)*(n-1), (4*n-2)*(n-1) -1]
      // connect point j on row i with point j on row i+1
      for (i = 0; i < n-1; i++) { // rows [0, n-2]
         for (j = 0; j < n-1; j++) {  // points on a row [0, n-2]
            buff.fSegs[indx+(i*(n-1)+j)*3] = c;                // lighter color
            buff.fSegs[indx+(i*(n-1)+j)*3+1] = indp + i*(n-1)+j;     // j on row i
            buff.fSegs[indx+(i*(n-1)+j)*3+2] = indp + (i+1)*(n-1)+j; // j on row i+1
         }
      }
      indx += 3*(n-1)*(n-1);
      startcap = (4*n-2)*(n-1);
   }

   if (hasphi) {
      if (hasrmin) {
         // endcaps = 2*(n-1) -> [(4*n-2)*(n-1), 4*n*(n-1)-1]
         i = 0;
         for (j = 0; j < n-1; j++) {
            buff.fSegs[indx+j*3] = c+1;
            buff.fSegs[indx+j*3+1] = (n-1)*i+j;     // outer j on row 0
            buff.fSegs[indx+j*3+2] = indp+(n-1)*i+j; // inner j on row 0
         }
         indx += 3*(n-1);
         i = n-1;
         for (j = 0; j < n-1; j++) {
            buff.fSegs[indx+j*3] = c+1;
            buff.fSegs[indx+j*3+1] = (n-1)*i+j;     // outer j on row n-1
            buff.fSegs[indx+j*3+2] = indp+(n-1)*i+j; // inner j on row n-1
         }
         indx += 3*(n-1);
      } else {
         i = 0;
         for (j = 0; j < n-1; j++) {
            buff.fSegs[indx+j*3] = c+1;
            buff.fSegs[indx+j*3+1] = (n-1)*i+j;     // outer j on row 0
            buff.fSegs[indx+j*3+2] = n*(n-1);       // center of first endcap
         }
         indx += 3*(n-1);
         i = n-1;
         for (j = 0; j < n-1; j++) {
            buff.fSegs[indx+j*3] = c+1;
            buff.fSegs[indx+j*3+1] = (n-1)*i+j;     // outer j on row n-1
            buff.fSegs[indx+j*3+2] = n*(n-1)+1;     // center of second endcap
         }
         indx += 3*(n-1);
      }
   }

   indx = 0;
   memset(buff.fPols, 0, buff.NbPols()*6*sizeof(Int_t));

   // outer surface = (n-1)*(n-1) -> [0, (n-1)*(n-1)-1]
   // normal pointing out
   for (i=0; i<n-1; i++) {
      for (j=0; j<n-1; j++) {
         buff.fPols[indx++] = c;
         buff.fPols[indx++] = 4;
         buff.fPols[indx++] = n*(n-1)+(n-1)*i+((j+1)%(n-1)); // generator j+1 on outer row i
         buff.fPols[indx++] = (n-1)*(i+1)+j; // seg j on outer row i+1
         buff.fPols[indx++] = n*(n-1)+(n-1)*i+j; // generator j on outer row i
         buff.fPols[indx++] = (n-1)*i+j; // seg j on outer row i
      }
   }
   if (hasrmin) {
      indp = (2*n-1)*(n-1); // start index of inner segments
      // inner surface = (n-1)*(n-1) -> [(n-1)*(n-1), 2*(n-1)*(n-1)-1]
      // normal pointing out
      for (i=0; i<n-1; i++) {
         for (j=0; j<n-1; j++) {
            buff.fPols[indx++] = c;
            buff.fPols[indx++] = 4;
            buff.fPols[indx++] = indp+n*(n-1)+(n-1)*i+j; // generator j on inner row i
            buff.fPols[indx++] = indp+(n-1)*(i+1)+j; // seg j on inner row i+1
            buff.fPols[indx++] = indp+n*(n-1)+(n-1)*i+((j+1)%(n-1)); // generator j+1 on inner r>
            buff.fPols[indx++] = indp+(n-1)*i+j; // seg j on inner row i
         }
      }
   }
   if (hasphi) {
      // endcaps = 2*(n-1) -> [2*(n-1)*(n-1), 2*n*(n-1)-1]
      i=0; // row 0
      Int_t np = (hasrmin)?4:3;
      for (j=0; j<n-1; j++) {
         buff.fPols[indx++] = c+1;
         buff.fPols[indx++] = np;
         buff.fPols[indx++] = j;         // seg j on outer row 0  a
         buff.fPols[indx++] = startcap+j;        // endcap j on row 0  d
         if(hasrmin)
            buff.fPols[indx++] = indp+j; // seg j on inner row 0  c
         buff.fPols[indx++] = startcap+((j+1)%(n-1)); // endcap j+1 on row 0  b
      }

      i=n-1; // row n-1
      for (j=0; j<n-1; j++) {
         buff.fPols[indx++] = c+1;
         buff.fPols[indx++] = np;
         buff.fPols[indx++] = (n-1)*i+j;         // seg j on outer row n-1 a
         buff.fPols[indx++] = startcap+(n-1)+((j+1)%(n-1));    // endcap j+1 on row n-1 d
         if (hasrmin)
            buff.fPols[indx++] = indp+(n-1)*i+j; // seg j on inner row n-1 c
         buff.fPols[indx++] = startcap+(n-1)+j;      // endcap j on row n-1 b
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// computes the closest distance from given point to this shape, according
/// to option. The matching point on the shape is stored in spoint.

Double_t TGeoTorus::Safety(const Double_t *point, Bool_t in) const
{
   Double_t saf[2];
   Int_t i;
   Double_t rxy = TMath::Sqrt(point[0]*point[0]+point[1]*point[1]);
   Double_t rad = TMath::Sqrt((rxy-fR)*(rxy-fR) + point[2]*point[2]);
   saf[0] = rad-fRmin;
   saf[1] = fRmax-rad;
   if (TGeoShape::IsSameWithinTolerance(fDphi,360)) {
      if (in) return TMath::Min(saf[0],saf[1]);
      for (i=0; i<2; i++) saf[i]=-saf[i];
      return TMath::Max(saf[0], saf[1]);
   }

   Double_t safphi = TGeoShape::SafetyPhi(point,in,fPhi1, fPhi1+fDphi);
   Double_t safe = TGeoShape::Big();
   if (in) {
      safe = TMath::Min(saf[0], saf[1]);
      return TMath::Min(safe, safphi);
   }
   for (i=0; i<2; i++) saf[i]=-saf[i];
   safe = TMath::Max(saf[0], saf[1]);
   return TMath::Max(safe, safphi);
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoTorus::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive)) return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   r    = " << fR << ";" << std::endl;
   out << "   rmin = " << fRmin << ";" << std::endl;
   out << "   rmax = " << fRmax << ";" << std::endl;
   out << "   phi1 = " << fPhi1 << ";" << std::endl;
   out << "   dphi = " << fDphi << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoTorus(\"" << GetName() << "\",r,rmin,rmax,phi1,dphi);" << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Set torus dimensions.

void TGeoTorus::SetTorusDimensions(Double_t r, Double_t rmin, Double_t rmax,
                          Double_t phi1, Double_t dphi)
{
   fR = r;
   fRmin = rmin;
   fRmax = rmax;
   fPhi1 = phi1;
   if (fPhi1<0) fPhi1+=360.;
   fDphi = dphi;
}

////////////////////////////////////////////////////////////////////////////////
/// Set torus dimensions starting from a list.

void TGeoTorus::SetDimensions(Double_t *param)
{
   SetTorusDimensions(param[0], param[1], param[2], param[3], param[4]);
}

////////////////////////////////////////////////////////////////////////////////
/// Create torus mesh points

void TGeoTorus::SetPoints(Double_t *points) const
{
   if (!points) return;
   Int_t n = gGeoManager->GetNsegments()+1;
   Double_t phin, phout;
   Double_t dpin = 360./(n-1);
   Double_t dpout = fDphi/(n-1);
   Double_t co,so,ci,si;
   Bool_t havermin = (fRmin<TGeoShape::Tolerance())?kFALSE:kTRUE;
   Int_t i,j;
   Int_t indx = 0;
   // loop outer mesh -> n*(n-1) points [0, 3*n*(n-1)-1]
   for (i=0; i<n; i++) {
      phout = (fPhi1+i*dpout)*TMath::DegToRad();
      co = TMath::Cos(phout);
      so = TMath::Sin(phout);
      for (j=0; j<n-1; j++) {
         phin = j*dpin*TMath::DegToRad();
         ci = TMath::Cos(phin);
         si = TMath::Sin(phin);
         points[indx++] = (fR+fRmax*ci)*co;
         points[indx++] = (fR+fRmax*ci)*so;
         points[indx++] = fRmax*si;
      }
   }

   if (havermin) {
    // loop inner mesh -> n*(n-1) points [3*n*(n-1), 6*n*(n-1)]
      for (i=0; i<n; i++) {
         phout = (fPhi1+i*dpout)*TMath::DegToRad();
         co = TMath::Cos(phout);
         so = TMath::Sin(phout);
         for (j=0; j<n-1; j++) {
            phin = j*dpin*TMath::DegToRad();
            ci = TMath::Cos(phin);
            si = TMath::Sin(phin);
            points[indx++] = (fR+fRmin*ci)*co;
            points[indx++] = (fR+fRmin*ci)*so;
            points[indx++] = fRmin*si;
         }
      }
   } else {
      if (fDphi<360.) {
      // just add extra 2 points on the centers of the 2 phi cuts [3*n*n, 3*n*n+1]
         points[indx++] = fR*TMath::Cos(fPhi1*TMath::DegToRad());
         points[indx++] = fR*TMath::Sin(fPhi1*TMath::DegToRad());
         points[indx++] = 0;
         points[indx++] = fR*TMath::Cos((fPhi1+fDphi)*TMath::DegToRad());
         points[indx++] = fR*TMath::Sin((fPhi1+fDphi)*TMath::DegToRad());
         points[indx++] = 0;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create torus mesh points

void TGeoTorus::SetPoints(Float_t *points) const
{
   if (!points) return;
   Int_t n = gGeoManager->GetNsegments()+1;
   Double_t phin, phout;
   Double_t dpin = 360./(n-1);
   Double_t dpout = fDphi/(n-1);
   Double_t co,so,ci,si;
   Bool_t havermin = (fRmin<TGeoShape::Tolerance())?kFALSE:kTRUE;
   Int_t i,j;
   Int_t indx = 0;
   // loop outer mesh -> n*(n-1) points [0, 3*n*(n-1)-1]
   // plane i = 0, n-1  point j = 0, n-1  ipoint = n*i + j
   for (i=0; i<n; i++) {
      phout = (fPhi1+i*dpout)*TMath::DegToRad();
      co = TMath::Cos(phout);
      so = TMath::Sin(phout);
      for (j=0; j<n-1; j++) {
         phin = j*dpin*TMath::DegToRad();
         ci = TMath::Cos(phin);
         si = TMath::Sin(phin);
         points[indx++] = (fR+fRmax*ci)*co;
         points[indx++] = (fR+fRmax*ci)*so;
         points[indx++] = fRmax*si;
      }
   }

   if (havermin) {
    // loop inner mesh -> n*(n-1) points [3*n*(n-1), 6*n*(n-1)]
      // plane i = 0, n-1  point j = 0, n-1  ipoint = n*n + n*i + j
      for (i=0; i<n; i++) {
         phout = (fPhi1+i*dpout)*TMath::DegToRad();
         co = TMath::Cos(phout);
         so = TMath::Sin(phout);
         for (j=0; j<n-1; j++) {
            phin = j*dpin*TMath::DegToRad();
            ci = TMath::Cos(phin);
            si = TMath::Sin(phin);
            points[indx++] = (fR+fRmin*ci)*co;
            points[indx++] = (fR+fRmin*ci)*so;
            points[indx++] = fRmin*si;
         }
      }
   } else {
      if (fDphi<360.) {
      // just add extra 2 points on the centers of the 2 phi cuts [n*n, n*n+1]
      // ip1 = n*(n-1) + 0;
      // ip2 = n*(n-1) + 1
         points[indx++] = fR*TMath::Cos(fPhi1*TMath::DegToRad());
         points[indx++] = fR*TMath::Sin(fPhi1*TMath::DegToRad());
         points[indx++] = 0;
         points[indx++] = fR*TMath::Cos((fPhi1+fDphi)*TMath::DegToRad());
         points[indx++] = fR*TMath::Sin((fPhi1+fDphi)*TMath::DegToRad());
         points[indx++] = 0;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of vertices of the mesh representation

Int_t TGeoTorus::GetNmeshVertices() const
{
   Int_t n = gGeoManager->GetNsegments()+1;
   Int_t numPoints = n*(n-1);
   if (fRmin>TGeoShape::Tolerance()) numPoints *= 2;
   else if (fDphi<360.)              numPoints += 2;
   return numPoints;
}

////////////////////////////////////////////////////////////////////////////////
/// fill size of this 3-D object

void TGeoTorus::Sizeof3D() const
{
}

////////////////////////////////////////////////////////////////////////////////
/// Find real solutions of the cubic equation : x^3 + a*x^2 + b*x + c = 0
/// Input: a,b,c
/// Output: x[3] real solutions
/// Returns number of real solutions (1 or 3)

Int_t TGeoTorus::SolveCubic(Double_t a, Double_t b, Double_t c, Double_t *x) const
{
   const Double_t ott = 1./3.;
   const Double_t sq3 = TMath::Sqrt(3.);
   Int_t ireal = 1;
   Double_t p = b-a*a*ott;
   Double_t q = c-a*b*ott+2.*a*a*a*ott*ott*ott;
   Double_t delta = 4*p*p*p+27*q*q;
//   Double_t y1r, y1i, y2r, y2i;
   Double_t t,u;
   if (delta>=0) {
      delta = TMath::Sqrt(delta);
      t = (-3*q*sq3+delta)/(6*sq3);
      u = (3*q*sq3+delta)/(6*sq3);
      x[0] = TMath::Sign(1.,t)*TMath::Power(TMath::Abs(t),ott)-
             TMath::Sign(1.,u)*TMath::Power(TMath::Abs(u),ott)-a*ott;
   } else {
      delta = TMath::Sqrt(-delta);
      t = -0.5*q;
      u = delta/(6*sq3);
      x[0] = 2.*TMath::Power(t*t+u*u,0.5*ott) * TMath::Cos(ott*TMath::ATan2(u,t));
      x[0] -= a*ott;
   }

   t = x[0]*x[0]+a*x[0]+b;
   u = a+x[0];
   delta = u*u-4.*t;
   if (delta>=0) {
      ireal = 3;
      delta = TMath::Sqrt(delta);
      x[1] = 0.5*(-u-delta);
      x[2] = 0.5*(-u+delta);
   }
   return ireal;
}

////////////////////////////////////////////////////////////////////////////////
/// Find real solutions of the quartic equation : x^4 + a*x^3 + b*x^2 + c*x + d = 0
/// Input: a,b,c,d
/// Output: x[4] - real solutions
/// Returns number of real solutions (0 to 3)

Int_t TGeoTorus::SolveQuartic(Double_t a, Double_t b, Double_t c, Double_t d, Double_t *x) const
{
   Double_t e = b-3.*a*a/8.;
   Double_t f = c+a*a*a/8.-0.5*a*b;
   Double_t g = d-3.*a*a*a*a/256. + a*a*b/16. - a*c/4.;
   Double_t xx[4];
   Int_t    ind[4];
   Double_t delta;
   Double_t h=0;
   Int_t ireal = 0;
   Int_t i;
   if (TGeoShape::IsSameWithinTolerance(f,0)) {
      delta = e*e-4.*g;
      if (delta<0) return 0;
      delta = TMath::Sqrt(delta);
      h = 0.5*(-e-delta);
      if (h>=0) {
         h = TMath::Sqrt(h);
         x[ireal++] = -h-0.25*a;
         x[ireal++] = h-0.25*a;
      }
      h = 0.5*(-e+delta);
      if (h>=0) {
         h = TMath::Sqrt(h);
         x[ireal++] = -h-0.25*a;
         x[ireal++] = h-0.25*a;
      }
      if (ireal>0) {
         TMath::Sort(ireal, x, ind,kFALSE);
         for (i=0; i<ireal; i++) xx[i] = x[ind[i]];
         memcpy(x,xx,ireal*sizeof(Double_t));
      }
      return ireal;
   }

   if (TGeoShape::IsSameWithinTolerance(g,0)) {
      x[ireal++] = -0.25*a;
      ind[0] = SolveCubic(0,e,f,xx);
      for (i=0; i<ind[0]; i++) x[ireal++] = xx[i]-0.25*a;
      if (ireal>0) {
         TMath::Sort(ireal, x, ind,kFALSE);
         for (i=0; i<ireal; i++) xx[i] = x[ind[i]];
         memcpy(x,xx,ireal*sizeof(Double_t));
      }
      return ireal;
   }


   ireal = SolveCubic(2.*e, e*e-4.*g, -f*f, xx);
   if (ireal==1) {
      if (xx[0]<=0) return 0;
      h = TMath::Sqrt(xx[0]);
   } else {
      // 3 real solutions of the cubic
      for (i=0; i<3; i++) {
         h = xx[i];
         if (h>=0) break;
      }
      if (h<=0) return 0;
      h = TMath::Sqrt(h);
   }
   Double_t j = 0.5*(e+h*h-f/h);
   ireal = 0;
   delta = h*h-4.*j;
   if (delta>=0) {
      delta = TMath::Sqrt(delta);
      x[ireal++] = 0.5*(-h-delta)-0.25*a;
      x[ireal++] = 0.5*(-h+delta)-0.25*a;
   }
   delta = h*h-4.*g/j;
   if (delta>=0) {
      delta = TMath::Sqrt(delta);
      x[ireal++] = 0.5*(h-delta)-0.25*a;
      x[ireal++] = 0.5*(h+delta)-0.25*a;
   }
   if (ireal>0) {
      TMath::Sort(ireal, x, ind,kFALSE);
      for (i=0; i<ireal; i++) xx[i] = x[ind[i]];
      memcpy(x,xx,ireal*sizeof(Double_t));
   }
   return ireal;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns distance to the surface or the torus (fR,r) from a point, along
/// a direction. Point is close enough to the boundary so that the distance
/// to the torus is decreasing while moving along the given direction.

Double_t TGeoTorus::ToBoundary(const Double_t *pt, const Double_t *dir, Double_t r, Bool_t in) const
{
   // Compute coefficients of the quartic
   Double_t s = TGeoShape::Big();
   Double_t tol = TGeoShape::Tolerance();
   Double_t r0sq  = pt[0]*pt[0]+pt[1]*pt[1]+pt[2]*pt[2];
   Double_t rdotn = pt[0]*dir[0]+pt[1]*dir[1]+pt[2]*dir[2];
   Double_t rsumsq = fR*fR+r*r;
   Double_t a = 4.*rdotn;
   Double_t b = 2.*(r0sq+2.*rdotn*rdotn-rsumsq+2.*fR*fR*dir[2]*dir[2]);
   Double_t c = 4.*(r0sq*rdotn-rsumsq*rdotn+2.*fR*fR*pt[2]*dir[2]);
   Double_t d = r0sq*r0sq-2.*r0sq*rsumsq+4.*fR*fR*pt[2]*pt[2]+(fR*fR-r*r)*(fR*fR-r*r);

   Double_t x[4],y[4];
   Int_t nsol = 0;

   if (TMath::Abs(dir[2])<1E-3 && TMath::Abs(pt[2])<0.1*r) {
      Double_t r0 = fR - TMath::Sqrt((r-pt[2])*(r+pt[2]));
      Double_t b0 = (pt[0]*dir[0]+pt[1]*dir[1])/(dir[0]*dir[0]+dir[1]*dir[1]);
      Double_t c0 = (pt[0]*pt[0] + (pt[1]-r0)*(pt[1]+r0))/(dir[0]*dir[0]+dir[1]*dir[1]);
      Double_t delta = b0*b0-c0;
      if (delta>0) {
         y[nsol] = -b0-TMath::Sqrt(delta);
         if (y[nsol]>-tol) nsol++;
         y[nsol] = -b0+TMath::Sqrt(delta);
         if (y[nsol]>-tol) nsol++;
      }
      r0 = fR + TMath::Sqrt((r-pt[2])*(r+pt[2]));
      c0 = (pt[0]*pt[0] + (pt[1]-r0)*(pt[1]+r0))/(dir[0]*dir[0]+dir[1]*dir[1]);
      delta = b0*b0-c0;
      if (delta>0) {
         y[nsol] = -b0-TMath::Sqrt(delta);
         if (y[nsol]>-tol) nsol++;
         y[nsol] = -b0+TMath::Sqrt(delta);
         if (y[nsol]>-tol) nsol++;
      }
      if (nsol) {
         // Sort solutions
         Int_t ind[4];
         TMath::Sort(nsol, y, ind,kFALSE);
         for (Int_t j=0; j<nsol; j++) x[j] = y[ind[j]];
      }
   } else {
      nsol = SolveQuartic(a,b,c,d,x);
   }
   if (!nsol) return TGeoShape::Big();
   // look for first positive solution
   Double_t phi, ndotd;
   Double_t r0[3], norm[3];
   Bool_t inner = (TMath::Abs(r-fRmin)<TGeoShape::Tolerance())?kTRUE:kFALSE;
   for (Int_t i=0; i<nsol; i++) {
      if (x[i]<-10) continue;
      phi = TMath::ATan2(pt[1]+x[i]*dir[1],pt[0]+x[i]*dir[0]);
      r0[0] = fR*TMath::Cos(phi);
      r0[1] = fR*TMath::Sin(phi);
      r0[2] = 0;
      for (Int_t ipt=0; ipt<3; ipt++) norm[ipt] = pt[ipt]+x[i]*dir[ipt] - r0[ipt];
      ndotd = norm[0]*dir[0]+norm[1]*dir[1]+norm[2]*dir[2];
      if (inner^in) {
         if (ndotd<0) continue;
      } else {
         if (ndotd>0) continue;
      }
      s = x[i];
      Double_t eps = TGeoShape::Big();
      Double_t delta = s*s*s*s + a*s*s*s + b*s*s + c*s + d;
      Double_t eps0 = -delta/(4.*s*s*s + 3.*a*s*s + 2.*b*s + c);
      while (TMath::Abs(eps)>TGeoShape::Tolerance()) {
         if (TMath::Abs(eps0)>100) break;
         s += eps0;
         if (TMath::Abs(s+eps0)<TGeoShape::Tolerance()) break;
         delta = s*s*s*s + a*s*s*s + b*s*s + c*s + d;
         eps = -delta/(4.*s*s*s + 3.*a*s*s + 2.*b*s + c);
         if (TMath::Abs(eps)>TMath::Abs(eps0)) break;
         eps0 = eps;
      }
      if (s<-TGeoShape::Tolerance()) continue;
      return TMath::Max(0.,s);
   }
   return TGeoShape::Big();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoTorus::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   Int_t n = gGeoManager->GetNsegments()+1;
   nvert = n*(n-1);
   Bool_t hasrmin = (GetRmin()>0)?kTRUE:kFALSE;
   Bool_t hasphi  = (GetDphi()<360)?kTRUE:kFALSE;
   if (hasrmin) nvert *= 2;
   else if (hasphi) nvert += 2;
   nsegs = (2*n-1)*(n-1);
   npols = (n-1)*(n-1);
   if (hasrmin) {
      nsegs += (2*n-1)*(n-1);
      npols += (n-1)*(n-1);
   }
   if (hasphi) {
      nsegs += 2*(n-1);
      npols += 2*(n-1);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fills a static 3D buffer and returns a reference.

const TBuffer3D & TGeoTorus::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
{
   static TBuffer3D buffer(TBuffer3DTypes::kGeneric);

   TGeoBBox::FillBuffer3D(buffer, reqSections, localFrame);

   if (reqSections & TBuffer3D::kRawSizes) {
      Int_t n = gGeoManager->GetNsegments()+1;
      Int_t nbPnts = n*(n-1);
      Bool_t hasrmin = (GetRmin()>0)?kTRUE:kFALSE;
      Bool_t hasphi  = (GetDphi()<360)?kTRUE:kFALSE;
      if (hasrmin) nbPnts *= 2;
      else if (hasphi) nbPnts += 2;

      Int_t nbSegs = (2*n-1)*(n-1);
      Int_t nbPols = (n-1)*(n-1);
      if (hasrmin) {
         nbSegs += (2*n-1)*(n-1);
         nbPols += (n-1)*(n-1);
      }
      if (hasphi) {
         nbSegs += 2*(n-1);
         nbPols += 2*(n-1);
      }

      if (buffer.SetRawSizes(nbPnts, 3*nbPnts, nbSegs, 3*nbSegs, nbPols, 6*nbPols)) {
         buffer.SetSectionsValid(TBuffer3D::kRawSizes);
      }
   }
   // TODO: Push down to TGeoShape?? But would have to do raw sizes set first..
   // can rest of TGeoShape be deferred until after
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

void TGeoTorus::Contains_v(const Double_t *points, Bool_t *inside, Int_t vecsize) const
{
   for (Int_t i=0; i<vecsize; i++) inside[i] = Contains(&points[3*i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the normal for an array o points so that norm.dot.dir is positive
/// Input: Arrays of point coordinates and directions + vector size
/// Output: Array of normal directions

void TGeoTorus::ComputeNormal_v(const Double_t *points, const Double_t *dirs, Double_t *norms, Int_t vecsize)
{
   for (Int_t i=0; i<vecsize; i++) ComputeNormal(&points[3*i], &dirs[3*i], &norms[3*i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoTorus::DistFromInside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromInside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoTorus::DistFromOutside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromOutside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute safe distance from each of the points in the input array.
/// Input: Array of point coordinates, array of statuses for these points, size of the arrays
/// Output: Safety values

void TGeoTorus::Safety_v(const Double_t *points, const Bool_t *inside, Double_t *safe, Int_t vecsize) const
{
   for (Int_t i=0; i<vecsize; i++) safe[i] = Safety(&points[3*i], inside[i]);
}
