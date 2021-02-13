// @(#)root/geom:$Id$
// Author: Mihaela Gheata   20/11/04

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#include <iostream>

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TVirtualGeoPainter.h"
#include "TGeoHype.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"

/** \class TGeoHype
\ingroup Geometry_classes

Hyperboloid class defined by 5 parameters. Bounded by:
  - Two z planes at z=+/-dz
  - Inner and outer lateral surfaces. These represent the surfaces
    described by the revolution of 2 hyperbolas about the Z axis:
     r^2 - (t*z)^2 = a^2

  - r = distance between hyperbola and Z axis at coordinate z
  - t = tangent of the stereo angle (angle made by hyperbola
        asymptotic lines and Z axis). t=0 means cylindrical surface.
  - a = distance between hyperbola and Z axis at z=0

The inner hyperbolic surface is described by:
    r^2 - (tin*z)^2 = rin^2
  - absence of the inner surface (filled hyperboloid can be forced
    by rin=0 and sin=0
The outer hyperbolic surface is described by:
    r^2 - (tout*z)^2 = rout^2
TGeoHype parameters: dz[cm], rin[cm], sin[deg], rout[cm], sout[deg].
MANDATORY conditions:

  - rin < rout
  - rout > 0
  - rin^2 + (tin*dz)^2 > rout^2 + (tout*dz)^2

SUPPORTED CASES:

  - rin = 0, tin != 0     => inner surface conical
  - tin=0 AND/OR tout=0   => corresponding surface(s) cylindrical
    e.g. tin=0 AND tout=0 => shape becomes a tube with: rmin,rmax,dz
*/

ClassImp(TGeoHype);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoHype::TGeoHype()
{
   SetShapeBit(TGeoShape::kGeoHype);
   fStIn = 0.;
   fStOut = 0.;
   fTin = 0.;
   fTinsq = 0.;
   fTout = 0.;
   fToutsq = 0.;
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor specifying hyperboloid parameters.

TGeoHype::TGeoHype(Double_t rin, Double_t stin, Double_t rout, Double_t stout, Double_t dz)
         :TGeoTube(rin, rout, dz)
{
   SetShapeBit(TGeoShape::kGeoHype);
   SetHypeDimensions(rin, stin, rout, stout, dz);
   // dz<0 can be used to force dz of hyperboloid fit the container volume
   if (fDz<0) SetShapeBit(kGeoRunTimeShape);
   ComputeBBox();
}
////////////////////////////////////////////////////////////////////////////////
/// Constructor specifying parameters and name.

TGeoHype::TGeoHype(const char *name,Double_t rin, Double_t stin, Double_t rout, Double_t stout, Double_t dz)
         :TGeoTube(name, rin, rout, dz)
{
   SetShapeBit(TGeoShape::kGeoHype);
   SetHypeDimensions(rin, stin, rout, stout, dz);
   // dz<0 can be used to force dz of hyperboloid fit the container volume
   if (fDz<0) SetShapeBit(kGeoRunTimeShape);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying a list of parameters
///  - param[0] = dz
///  - param[1] = rin
///  - param[2] = stin
///  - param[3] = rout
///  - param[4] = stout

TGeoHype::TGeoHype(Double_t *param)
         :TGeoTube(param[1],param[3],param[0])
{
   SetShapeBit(TGeoShape::kGeoHype);
   SetDimensions(param);
   // dz<0 can be used to force dz of hyperboloid fit the container volume
   if (fDz<0) SetShapeBit(kGeoRunTimeShape);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoHype::~TGeoHype()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Computes capacity of the shape in [length^3]

Double_t TGeoHype::Capacity() const
{
   Double_t capacity = 2.*TMath::Pi()*fDz*(fRmax*fRmax-fRmin*fRmin) +
                       (2.*TMath::Pi()/3.)*fDz*fDz*fDz*(fToutsq-fTinsq);
   return capacity;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute bounding box of the hyperboloid

void TGeoHype::ComputeBBox()
{
   if (fRmin<0.) {
      Warning("ComputeBBox", "Shape %s has invalid rmin=%g ! SET TO 0.", GetName(),fRmin);
      fRmin = 0.;
   }
   if ((fRmin>fRmax) || (fRmin*fRmin+fTinsq*fDz*fDz > fRmax*fRmax+fToutsq*fDz*fDz)) {
      SetShapeBit(kGeoInvalidShape);
      Error("ComputeBBox", "Shape %s hyperbolic surfaces are malformed: rin=%g, stin=%g, rout=%g, stout=%g",
             GetName(), fRmin, fStIn, fRmax, fStOut);
      return;
   }

   fDX = fDY = TMath::Sqrt(RadiusHypeSq(fDz, kFALSE));
   fDZ = fDz;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute normal to closest surface from POINT.

void TGeoHype::ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm)
{
   Double_t saf[3];
   Double_t rsq = point[0]*point[0]+point[1]*point[1];
   Double_t r = TMath::Sqrt(rsq);
   Double_t rin = (HasInner())?(TMath::Sqrt(RadiusHypeSq(point[2],kTRUE))):0.;
   Double_t rout = TMath::Sqrt(RadiusHypeSq(point[2],kFALSE));
   saf[0] = TMath::Abs(fDz-TMath::Abs(point[2]));
   saf[1] = (HasInner())?TMath::Abs(rin-r):TGeoShape::Big();
   saf[2] = TMath::Abs(rout-r);
   Int_t i = TMath::LocMin(3,saf);
   if (i==0 || r<1.E-10) {
      norm[0] = norm[1] = 0.;
      norm[2] = TMath::Sign(1.,dir[2]);
      return;
   }
   Double_t t = (i==1)?fTinsq:fToutsq;;
   t *= -point[2]/r;
   Double_t ct = TMath::Sqrt(1./(1.+t*t));
   Double_t st = t * ct;
   Double_t phi = TMath::ATan2(point[1], point[0]);
   Double_t cphi = TMath::Cos(phi);
   Double_t sphi = TMath::Sin(phi);

   norm[0] = ct*cphi;
   norm[1] = ct*sphi;
   norm[2] = st;
   if (norm[0]*dir[0]+norm[1]*dir[1]+norm[2]*dir[2]<0) {
      norm[0] = -norm[0];
      norm[1] = -norm[1];
      norm[2] = -norm[2];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// test if point is inside this tube

Bool_t TGeoHype::Contains(const Double_t *point) const
{
   if (TMath::Abs(point[2]) > fDz) return kFALSE;
   Double_t r2 = point[0]*point[0]+point[1]*point[1];
   Double_t routsq = RadiusHypeSq(point[2], kFALSE);
   if (r2>routsq) return kFALSE;
   if (!HasInner()) return kTRUE;
   Double_t rinsq = RadiusHypeSq(point[2], kTRUE);
   if (r2<rinsq) return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// compute closest distance from point px,py to each corner

Int_t TGeoHype::DistancetoPrimitive(Int_t px, Int_t py)
{
   Int_t numPoints = GetNmeshVertices();
   return ShapeDistancetoPrimitive(numPoints, px, py);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from inside point to surface of the hyperboloid.

Double_t TGeoHype::DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
{
   if (iact<3 && safe) {
      *safe = Safety(point, kTRUE);
      if (iact==0) return TGeoShape::Big();
      if ((iact==1) && (*safe>step)) return TGeoShape::Big();
   }
   // compute distance to surface
   // Do Z
   Double_t sz = TGeoShape::Big();
   if (dir[2]>0) {
      sz = (fDz-point[2])/dir[2];
      if (sz<=0.) return 0.;
   } else {
      if (dir[2]<0) {
         sz = -(fDz+point[2])/dir[2];
         if (sz<=0.) return 0.;
      }
   }


   // Do R
   Double_t srin = TGeoShape::Big();
   Double_t srout = TGeoShape::Big();
   Double_t sr;
   // inner and outer surfaces
   Double_t s[2];
   Int_t npos;
   npos = DistToHype(point, dir, s, kTRUE, kTRUE);
   if (npos) srin = s[0];
   npos = DistToHype(point, dir, s, kFALSE, kTRUE);
   if (npos) srout = s[0];
   sr = TMath::Min(srin, srout);
   return TMath::Min(sz,sr);
}


////////////////////////////////////////////////////////////////////////////////
/// compute distance from outside point to surface of the hyperboloid.

Double_t TGeoHype::DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
{
   if (iact<3 && safe) {
      *safe = Safety(point, kFALSE);
      if (iact==0) return TGeoShape::Big();
      if ((iact==1) && (step<=*safe)) return TGeoShape::Big();
   }
// Check if the bounding box is crossed within the requested distance
   Double_t sdist = TGeoBBox::DistFromOutside(point,dir, fDX, fDY, fDZ, fOrigin, step);
   if (sdist>=step) return TGeoShape::Big();
   // find distance to shape
   // Do Z
   Double_t xi, yi, zi;
   Double_t sz = TGeoShape::Big();
   if (TMath::Abs(point[2])>=fDz) {
      // We might find Z plane crossing
      if ((point[2]*dir[2]) < 0) {
         // Compute distance to Z (always positive)
         sz = (TMath::Abs(point[2])-fDz)/TMath::Abs(dir[2]);
         // Extrapolate
         xi = point[0]+sz*dir[0];
         yi = point[1]+sz*dir[1];
         Double_t r2 = xi*xi + yi*yi;
         Double_t rmin2 = RadiusHypeSq(fDz, kTRUE);
         if (r2 >= rmin2) {
            Double_t rmax2 = RadiusHypeSq(fDz, kFALSE);
            if (r2 <= rmax2) return sz;
         }
      }
   }
   // We do not cross Z planes.
   Double_t sin = TGeoShape::Big();
   Double_t sout = TGeoShape::Big();
   Double_t s[2];
   Int_t npos;
   npos = DistToHype(point, dir, s, kTRUE, kFALSE);
   if (npos) {
      zi = point[2] + s[0]*dir[2];
      if (TMath::Abs(zi) <= fDz) sin = s[0];
      else if (npos==2) {
         zi = point[2] + s[1]*dir[2];
         if (TMath::Abs(zi) <= fDz) sin = s[1];
      }
   }
   npos = DistToHype(point, dir, s, kFALSE, kFALSE);
   if (npos) {
      zi = point[2] + s[0]*dir[2];
      if (TMath::Abs(zi) <= fDz) sout = s[0];
      else if (npos==2) {
         zi = point[2] + s[1]*dir[2];
         if (TMath::Abs(zi) <= fDz) sout = s[1];
      }
   }
   return TMath::Min(sin, sout);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from an arbitrary point to inner/outer surface of hyperboloid.
/// Returns number of positive solutions. S[2] contains the solutions.

Int_t TGeoHype::DistToHype(const Double_t *point, const Double_t *dir, Double_t *s, Bool_t inner, Bool_t in) const
{
   Double_t r0, t0, snext;
   if (inner) {
      if (!HasInner()) return 0;
      r0 = fRmin;
      t0 = fTinsq;
   } else {
      r0 = fRmax;
      t0 = fToutsq;
   }
   Double_t a = dir[0]*dir[0] + dir[1]*dir[1] - t0*dir[2]*dir[2];
   Double_t b = t0*point[2]*dir[2] - point[0]*dir[0] - point[1]*dir[1];
   Double_t c = point[0]*point[0] + point[1]*point[1] - t0*point[2]*point[2] - r0*r0;

   if (TMath::Abs(a) < TGeoShape::Tolerance()) {
      if (TMath::Abs(b) < TGeoShape::Tolerance()) return 0;
      snext = 0.5*c/b;
      if (snext < 0.) return 0;
      s[0] = snext;
      return 1;
   }

   Double_t delta = b*b - a*c;
   Double_t ainv = 1./a;
   Int_t npos = 0;
   if (delta < 0.) return 0;
   delta = TMath::Sqrt(delta);
   Double_t sone = TMath::Sign(1.,ainv);
   Int_t i = -1;
   while (i<2) {
      snext = (b + i*sone*delta)*ainv;
      i += 2;
      if (snext<0) continue;
      if (snext<1.E-8) {
         Double_t r = TMath::Sqrt(point[0]*point[0]+point[1]*point[1]);
         Double_t t = (inner)?fTinsq:fToutsq;
         t *= -point[2]/r;
         Double_t phi = TMath::ATan2(point[1], point[0]);
         Double_t ndotd = TMath::Cos(phi)*dir[0]+TMath::Sin(phi)*dir[1]+t*dir[2];
         if (inner) ndotd *= -1;
         if (in) ndotd *= -1;
         if (ndotd<0) s[npos++] = snext;
      } else          s[npos++] = snext;
   }
   return npos;
}

////////////////////////////////////////////////////////////////////////////////
/// Cannot divide hyperboloids.

TGeoVolume *TGeoHype::Divide(TGeoVolume * /*voldiv*/, const char *divname, Int_t /*iaxis*/, Int_t /*ndiv*/,
                             Double_t /*start*/, Double_t /*step*/)
{
   Error("Divide", "Hyperboloids cannot be divided. Division volume %s not created", divname);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get range of shape for a given axis.

Double_t TGeoHype::GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const
{
   xlo = 0;
   xhi = 0;
   Double_t dx = 0;
   switch (iaxis) {
      case 1: // R
         xlo = fRmin;
         xhi = TMath::Sqrt(RadiusHypeSq(fDz, kFALSE));
         dx = xhi-xlo;
         return dx;
      case 2: // Phi
         xlo = 0;
         xhi = 360;
         dx = 360;
         return dx;
      case 3: // Z
         xlo = -fDz;
         xhi = fDz;
         dx = xhi-xlo;
         return dx;
   }
   return dx;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill vector param[4] with the bounding cylinder parameters. The order
/// is the following : Rmin, Rmax, Phi1, Phi2, dZ

void TGeoHype::GetBoundingCylinder(Double_t *param) const
{
   param[0] = fRmin; // Rmin
   param[0] *= param[0];
   param[1] = TMath::Sqrt(RadiusHypeSq(fDz, kFALSE)); // Rmax
   param[1] *= param[1];
   param[2] = 0.;    // Phi1
   param[3] = 360.;  // Phi2
}

////////////////////////////////////////////////////////////////////////////////
/// in case shape has some negative parameters, these has to be computed
/// in order to fit the mother

TGeoShape *TGeoHype::GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix * /*mat*/) const
{
   if (!TestShapeBit(kGeoRunTimeShape)) return 0;
   Double_t dz;
   Double_t zmin,zmax;
   dz = fDz;
   if (fDz<0) {
      mother->GetAxisRange(3,zmin,zmax);
      if (zmax<0) return 0;
      dz=zmax;
   } else {
      Error("GetMakeRuntimeShape", "Shape %s does not have negative Z range", GetName());
      return 0;
   }
   TGeoShape *hype = new TGeoHype(GetName(), dz, fRmax, fStOut, fRmin, fStIn);
   return hype;
}

////////////////////////////////////////////////////////////////////////////////
/// print shape parameters

void TGeoHype::InspectShape() const
{
   printf("*** Shape %s: TGeoHype ***\n", GetName());
   printf("    Rin  = %11.5f\n", fRmin);
   printf("    sin  = %11.5f\n", fStIn);
   printf("    Rout = %11.5f\n", fRmax);
   printf("    sout = %11.5f\n", fStOut);
   printf("    dz   = %11.5f\n", fDz);

   printf(" Bounding box:\n");
   TGeoBBox::InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a TBuffer3D describing *this* shape.
/// Coordinates are in local reference frame.

TBuffer3D *TGeoHype::MakeBuffer3D() const
{
   Int_t n = gGeoManager->GetNsegments();
   Bool_t hasRmin = HasInner();
   Int_t nbPnts = (hasRmin)?(2*n*n):(n*n+2);
   Int_t nbSegs = (hasRmin)?(4*n*n):(n*(2*n+1));
   Int_t nbPols = (hasRmin)?(2*n*n):(n*(n+1));

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

void TGeoHype::SetSegsAndPols(TBuffer3D &buff) const
{
   Int_t c = GetBasicColor();
   Int_t i, j, n;
   n = gGeoManager->GetNsegments();
   Bool_t hasRmin = HasInner();
   Int_t irin = 0;
   Int_t irout = (hasRmin)?(n*n):2;
   // Fill segments
   // Case hasRmin:
   //   Inner circles:  [isin = 0], n (per circle) * n ( circles)
   //        iseg = isin+n*i+j , i = 0, n-1   , j = 0, n-1
   //        seg(i=1,n; j=1,n) = [irin+n*i+j] and [irin+n*i+(j+1)%n]
   //   Inner generators: [isgenin = isin+n*n], n (per circle) *(n-1) (slices)
   //        iseg = isgenin + i*n + j, i=0,n-2,  j=0,n-1
   //        seg(i,j) = [irin+n*i+j] and [irin+n*(i+1)+j]
   //   Outer circles:  [isout = isgenin+n*(n-1)], n (per circle) * n ( circles)
   //        iseg = isout + i*n + j , iz = 0, n-1   , j = 0, n-1
   //        seg(i=1,n; j=1,n) = [irout+n*i+j] and [irout+n*i+(j+1)%n]
   //   Outer generators: [isgenout = isout+n*n], n (per circle) *(n-1) (slices)
   //        iseg = isgenout + i*n + j, i=0,n-2,  j=0,n-1
   //        seg(i,j) = [irout+n*i+j] and [irout+n*(i+1)+j]
   //   Lower cap : [islow = isgenout + n*(n-1)], n radial segments
   //        iseg = islow + j,  j=0,n-1
   //        seg(j) = [irin + j] and [irout+j]
   //   Upper cap: [ishi = islow + n], nradial segments
   //        iseg = ishi + j, j=0,n-1
   //        seg[j] = [irin + n*(n-1) + j] and [irout+n*(n-1) + j]
   //
   // Case !hasRmin:
   //   Outer circles: [isout=0], same outer circles (n*n)
   // Outer generators: isgenout = isout + n*n
   //   Lower cap: [islow = isgenout+n*(n-1)], n seg.
   //        iseg = islow + j, j=0,n-1
   //        seg[j] = [irin] and [irout+j]
   //   Upper cap: [ishi = islow +n]
   //        iseg = ishi + j, j=0,n-1
   //        seg[j] = [irin+1] and [irout+n*(n-1) + j]

   Int_t isin = 0;
   Int_t isgenin = (hasRmin)?(isin+n*n):0;
   Int_t isout = (hasRmin)?(isgenin+n*(n-1)):0;
   Int_t isgenout  = isout+n*n;
   Int_t islo = isgenout+n*(n-1);
   Int_t ishi = islo + n;

   Int_t npt = 0;
   // Fill inner circle segments (n*n)
   if (hasRmin) {
      for (i=0; i<n; i++) {
         for (j=0; j<n; j++) {
            npt = 3*(isin+n*i+j);
            buff.fSegs[npt]   = c;
            buff.fSegs[npt+1] = irin+n*i+j;
            buff.fSegs[npt+2] = irin+n*i+((j+1)%n);
         }
      }
      // Fill inner generators (n*(n-1))
      for (i=0; i<n-1; i++) {
         for (j=0; j<n; j++) {
            npt = 3*(isgenin+n*i+j);
            buff.fSegs[npt]   = c;
            buff.fSegs[npt+1] = irin+n*i+j;
            buff.fSegs[npt+2] = irin+n*(i+1)+j;
         }
      }
   }
   // Fill outer circle segments (n*n)
   for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
         npt = 3*(isout + n*i+j);
         buff.fSegs[npt]   = c;
         buff.fSegs[npt+1] = irout+n*i+j;
         buff.fSegs[npt+2] = irout+n*i+((j+1)%n);
      }
   }
   // Fill outer generators (n*(n-1))
   for (i=0; i<n-1; i++) {
      for (j=0; j<n; j++) {
         npt = 3*(isgenout+n*i+j);
         buff.fSegs[npt]   = c;
         buff.fSegs[npt+1] = irout+n*i+j;
         buff.fSegs[npt+2] = irout+n*(i+1)+j;
      }
   }
   // Fill lower cap (n)
   for (j=0; j<n; j++) {
      npt = 3*(islo+j);
      buff.fSegs[npt]   = c;
      buff.fSegs[npt+1] = irin;
      if (hasRmin) buff.fSegs[npt+1] += j;
      buff.fSegs[npt+2] = irout + j;
   }
   // Fill upper cap (n)
   for (j=0; j<n; j++) {
      npt = 3*(ishi+j);
      buff.fSegs[npt]   = c;
      buff.fSegs[npt+1] = irin+1;
      if (hasRmin) buff.fSegs[npt+1] += n*(n-1)+j-1;
      buff.fSegs[npt+2] = irout + n*(n-1)+j;
   }

   // Fill polygons
   // Inner polygons: [ipin = 0] (n-1) slices * n (edges)
   //   ipoly = ipin + n*i + j;  i=0,n-2   j=0,n-1
   //   poly[i,j] = [isin+n*i+j]  [isgenin+i*n+(j+1)%n]  [isin+n*(i+1)+j]  [isgenin+i*n+j]
   // Outer polygons: [ipout = ipin+n*(n-1)]  also (n-1)*n
   //   ipoly = ipout + n*i + j; i=0,n-2   j=0,n-1
   //   poly[i,j] = [isout+n*i+j]  [isgenout+i*n+j]  [isout+n*(i+1)+j]  [isgenout+i*n+(j+1)%n]
   // Lower cap: [iplow = ipout+n*(n-1):  n polygons
   //   ipoly = iplow + j;  j=0,n-1
   //   poly[i=0,j] = [isin+j] [islow+j] [isout+j] [islow+(j+1)%n]
   // Upper cap: [ipup = iplow+n] : n polygons
   //   ipoly = ipup + j;  j=0,n-1
   //   poly[i=n-1, j] = [isin+n*(n-1)+j] [ishi+(j+1)%n] [isout+n*(n-1)+j] [ishi+j]
   //
   // Case !hasRmin:
   // ipin = 0 no inner polygons
   // ipout = 0 same outer polygons
   // Lower cap: iplow = ipout+n*(n-1):  n polygons with 3 segments
   //   poly[i=0,j] = [isout+j] [islow+(j+1)%n] [islow+j]
   // Upper cap: ipup = iplow+n;
   //   poly[i=n-1,j] = [isout+n*(n-1)+j] [ishi+j] [ishi+(j+1)%n]

   Int_t ipin = 0;
   Int_t ipout = (hasRmin)?(ipin+n*(n-1)):0;
   Int_t iplo = ipout+n*(n-1);
   Int_t ipup = iplo+n;
   // Inner polygons n*(n-1)
   if (hasRmin) {
      for (i=0; i<n-1; i++) {
         for (j=0; j<n; j++) {
            npt = 6*(ipin+n*i+j);
            buff.fPols[npt]   = c;
            buff.fPols[npt+1] = 4;
            buff.fPols[npt+2] = isin+n*i+j;
            buff.fPols[npt+3] = isgenin+i*n+((j+1)%n);
            buff.fPols[npt+4] = isin+n*(i+1)+j;
            buff.fPols[npt+5] = isgenin+i*n+j;
         }
      }
   }
   // Outer polygons n*(n-1)
   for (i=0; i<n-1; i++) {
      for (j=0; j<n; j++) {
         npt = 6*(ipout+n*i+j);
         buff.fPols[npt]   = c;
         buff.fPols[npt+1] = 4;
         buff.fPols[npt+2] = isout+n*i+j;
         buff.fPols[npt+3] = isgenout+i*n+j;
         buff.fPols[npt+4] = isout+n*(i+1)+j;
         buff.fPols[npt+5] = isgenout+i*n+((j+1)%n);
      }
   }
   // End caps
   if (hasRmin) {
      for (j=0; j<n; j++) {
         npt = 6*(iplo+j);
         buff.fPols[npt]   = c+1;
         buff.fPols[npt+1] = 4;
         buff.fPols[npt+2] = isin+j;
         buff.fPols[npt+3] = islo+j;
         buff.fPols[npt+4] = isout+j;
         buff.fPols[npt+5] = islo+((j+1)%n);
      }
      for (j=0; j<n; j++) {
         npt = 6*(ipup+j);
         buff.fPols[npt]   = c+2;
         buff.fPols[npt+1] = 4;
         buff.fPols[npt+2] = isin+n*(n-1)+j;
         buff.fPols[npt+3] = ishi+((j+1)%n);
         buff.fPols[npt+4] = isout+n*(n-1)+j;
         buff.fPols[npt+5] = ishi+j;
      }
   } else {
      for (j=0; j<n; j++) {
         npt = 6*iplo+5*j;
         buff.fPols[npt]   = c+1;
         buff.fPols[npt+1] = 3;
         buff.fPols[npt+2] = isout+j;
         buff.fPols[npt+3] = islo+((j+1)%n);
         buff.fPols[npt+4] = islo+j;
      }
      for (j=0; j<n; j++) {
         npt = 6*iplo+5*(n+j);
         buff.fPols[npt]   = c+2;
         buff.fPols[npt+1] = 3;
         buff.fPols[npt+2] = isout+n*(n-1)+j;
         buff.fPols[npt+3] = ishi+j;
         buff.fPols[npt+4] = ishi+((j+1)%n);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute r^2 = x^2 + y^2 at a given z coordinate, for either inner or outer hyperbolas.

Double_t TGeoHype::RadiusHypeSq(Double_t z, Bool_t inner) const
{
   Double_t r0, tsq;
   if (inner) {
      r0 = fRmin;
      tsq = fTinsq;
   } else {
      r0 = fRmax;
      tsq = fToutsq;
   }
   return (r0*r0+tsq*z*z);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute z^2 at a given  r^2, for either inner or outer hyperbolas.

Double_t TGeoHype::ZHypeSq(Double_t r, Bool_t inner) const
{
   Double_t r0, tsq;
   if (inner) {
      r0 = fRmin;
      tsq = fTinsq;
   } else {
      r0 = fRmax;
      tsq = fToutsq;
   }
   if (TMath::Abs(tsq) < TGeoShape::Tolerance()) return TGeoShape::Big();
   return ((r*r-r0*r0)/tsq);
}

////////////////////////////////////////////////////////////////////////////////
/// computes the closest distance from given point to this shape, according
/// to option. The matching point on the shape is stored in spoint.

Double_t TGeoHype::Safety(const Double_t *point, Bool_t in) const
{
   Double_t safe, safrmin, safrmax;
   if (in) {
      safe    = fDz-TMath::Abs(point[2]);
      safrmin = SafetyToHype(point, kTRUE, in);
      if (safrmin < safe) safe = safrmin;
      safrmax = SafetyToHype(point, kFALSE,in);
      if (safrmax < safe) safe = safrmax;
   } else {
      safe    = -fDz+TMath::Abs(point[2]);
      safrmin = SafetyToHype(point, kTRUE, in);
      if (safrmin > safe) safe = safrmin;
      safrmax = SafetyToHype(point, kFALSE,in);
      if (safrmax > safe) safe = safrmax;
   }
   return safe;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute an underestimate of the closest distance from a point to inner or
/// outer infinite hyperbolas.

Double_t TGeoHype::SafetyToHype(const Double_t *point, Bool_t inner, Bool_t in) const
{
   Double_t r, rsq, rhsq, rh, dr, tsq, saf;
   if (inner && !HasInner()) return (in)?TGeoShape::Big():-TGeoShape::Big();
   rsq = point[0]*point[0]+point[1]*point[1];
   r = TMath::Sqrt(rsq);
   rhsq = RadiusHypeSq(point[2], inner);
   rh = TMath::Sqrt(rhsq);
   dr = r - rh;
   if (inner) {
      if (!in && dr>0) return -TGeoShape::Big();
      if (TMath::Abs(fStIn) < TGeoShape::Tolerance()) return TMath::Abs(dr);
      if (fRmin<TGeoShape::Tolerance()) return TMath::Abs(dr/TMath::Sqrt(1.+ fTinsq));
      tsq = fTinsq;
   } else {
      if (!in && dr<0) return -TGeoShape::Big();
      if (TMath::Abs(fStOut) < TGeoShape::Tolerance()) return TMath::Abs(dr);
      tsq = fToutsq;
   }
   if (TMath::Abs(dr)<TGeoShape::Tolerance()) return 0.;
   // 1. dr<0 => approximate safety with distance to tangent to hyperbola in z = |point[2]|
   Double_t m;
   if (dr<0) {
      m = rh/(tsq*TMath::Abs(point[2]));
      saf = -m*dr/TMath::Sqrt(1.+m*m);
      return saf;
   }
   // 2. dr>0 => approximate safety with distance from point to segment P1(r(z0),z0) and P2(r0, z(r0))
   m = (TMath::Sqrt(ZHypeSq(r,inner)) - TMath::Abs(point[2]))/dr;
   saf = m*dr/TMath::Sqrt(1.+m*m);
   return saf;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoHype::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive)) return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   rin   = " << fRmin << ";" << std::endl;
   out << "   stin  = " << fStIn << ";" << std::endl;
   out << "   rout  = " << fRmax << ";" << std::endl;
   out << "   stout = " << fStOut << ";" << std::endl;
   out << "   dz    = " << fDz << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoHype(\"" << GetName() << "\",rin,stin,rout,stout,dz);" << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Set dimensions of the hyperboloid.

void TGeoHype::SetHypeDimensions(Double_t rin, Double_t stin, Double_t rout, Double_t stout, Double_t dz)
{
   fRmin = rin;
   fRmax = rout;
   fDz   = dz;
   fStIn = stin;
   fStOut = stout;
   fTin = TMath::Tan(fStIn*TMath::DegToRad());
   fTinsq = fTin*fTin;
   fTout = TMath::Tan(fStOut*TMath::DegToRad());
   fToutsq = fTout*fTout;
   if ((fRmin==0) && (fStIn==0)) SetShapeBit(kGeoRSeg, kTRUE);
   else                          SetShapeBit(kGeoRSeg, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Set dimensions of the hyperboloid starting from an array.
///  - param[0] = dz
///  - param[1] = rin
///  - param[2] = stin
///  - param[3] = rout
///  - param[4] = stout

void TGeoHype::SetDimensions(Double_t *param)
{
   Double_t dz = param[0];
   Double_t rin = param[1];
   Double_t stin = param[2];
   Double_t rout = param[3];
   Double_t stout = param[4];
   SetHypeDimensions(rin, stin, rout, stout, dz);
}

////////////////////////////////////////////////////////////////////////////////
/// create tube mesh points

void TGeoHype::SetPoints(Double_t *points) const
{
   Double_t z,dz,r;
   Int_t i,j, n;
   if (!points) return;
   n = gGeoManager->GetNsegments();
   Double_t dphi = 360./n;
   Double_t phi = 0;
   dz = 2.*fDz/(n-1);

   Int_t indx = 0;

   if (HasInner()) {
      // Inner surface points
      for (i=0; i<n; i++) {
         z = -fDz+i*dz;
         r = TMath::Sqrt(RadiusHypeSq(z, kTRUE));
         for (j=0; j<n; j++) {
            phi = j*dphi*TMath::DegToRad();
            points[indx++] = r * TMath::Cos(phi);
            points[indx++] = r * TMath::Sin(phi);
            points[indx++] = z;
         }
      }
   } else {
      points[indx++] = 0.;
      points[indx++] = 0.;
      points[indx++] = -fDz;
      points[indx++] = 0.;
      points[indx++] = 0.;
      points[indx++] = fDz;
   }
   // Outer surface points
   for (i=0; i<n; i++) {
      z = -fDz + i*dz;
      r = TMath::Sqrt(RadiusHypeSq(z, kFALSE));
      for (j=0; j<n; j++) {
         phi = j*dphi*TMath::DegToRad();
         points[indx++] = r * TMath::Cos(phi);
         points[indx++] = r * TMath::Sin(phi);
         points[indx++] = z;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// create tube mesh points

void TGeoHype::SetPoints(Float_t *points) const
{
   Double_t z,dz,r;
   Int_t i,j, n;
   if (!points) return;
   n = gGeoManager->GetNsegments();
   Double_t dphi = 360./n;
   Double_t phi = 0;
   dz = 2.*fDz/(n-1);

   Int_t indx = 0;

   if (HasInner()) {
      // Inner surface points
      for (i=0; i<n; i++) {
         z = -fDz+i*dz;
         r = TMath::Sqrt(RadiusHypeSq(z, kTRUE));
         for (j=0; j<n; j++) {
            phi = j*dphi*TMath::DegToRad();
            points[indx++] = r * TMath::Cos(phi);
            points[indx++] = r * TMath::Sin(phi);
            points[indx++] = z;
         }
      }
   } else {
      points[indx++] = 0.;
      points[indx++] = 0.;
      points[indx++] = -fDz;
      points[indx++] = 0.;
      points[indx++] = 0.;
      points[indx++] = fDz;
   }
   // Outer surface points
   for (i=0; i<n; i++) {
      z = -fDz + i*dz;
      r = TMath::Sqrt(RadiusHypeSq(z, kFALSE));
      for (j=0; j<n; j++) {
         phi = j*dphi*TMath::DegToRad();
         points[indx++] = r * TMath::Cos(phi);
         points[indx++] = r * TMath::Sin(phi);
         points[indx++] = z;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoHype::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   Int_t n = gGeoManager->GetNsegments();
   Bool_t hasRmin = HasInner();
   nvert = (hasRmin)?(2*n*n):(n*n+2);
   nsegs = (hasRmin)?(4*n*n):(n*(2*n+1));
   npols = (hasRmin)?(2*n*n):(n*(n+1));
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of vertices of the mesh representation

Int_t TGeoHype::GetNmeshVertices() const
{
   Int_t n = gGeoManager->GetNsegments();
   Int_t numPoints = (HasRmin())?(2*n*n):(n*n+2);
   return numPoints;
}

////////////////////////////////////////////////////////////////////////////////
/// fill size of this 3-D object

void TGeoHype::Sizeof3D() const
{
}

////////////////////////////////////////////////////////////////////////////////
/// Fills a static 3D buffer and returns a reference.

const TBuffer3D & TGeoHype::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
{
   static TBuffer3D buffer(TBuffer3DTypes::kGeneric);

   TGeoBBox::FillBuffer3D(buffer, reqSections, localFrame);

   if (reqSections & TBuffer3D::kRawSizes) {
      Int_t n = gGeoManager->GetNsegments();
      Bool_t hasRmin = HasInner();
      Int_t nbPnts = (hasRmin)?(2*n*n):(n*n+2);
      Int_t nbSegs = (hasRmin)?(4*n*n):(n*(2*n+1));
      Int_t nbPols = (hasRmin)?(2*n*n):(n*(n+1));
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

void TGeoHype::Contains_v(const Double_t *points, Bool_t *inside, Int_t vecsize) const
{
   for (Int_t i=0; i<vecsize; i++) inside[i] = Contains(&points[3*i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the normal for an array o points so that norm.dot.dir is positive
/// Input: Arrays of point coordinates and directions + vector size
/// Output: Array of normal directions

void TGeoHype::ComputeNormal_v(const Double_t *points, const Double_t *dirs, Double_t *norms, Int_t vecsize)
{
   for (Int_t i=0; i<vecsize; i++) ComputeNormal(&points[3*i], &dirs[3*i], &norms[3*i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoHype::DistFromInside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromInside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoHype::DistFromOutside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromOutside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute safe distance from each of the points in the input array.
/// Input: Array of point coordinates, array of statuses for these points, size of the arrays
/// Output: Safety values

void TGeoHype::Safety_v(const Double_t *points, const Bool_t *inside, Double_t *safe, Int_t vecsize) const
{
   for (Int_t i=0; i<vecsize; i++) safe[i] = Safety(&points[3*i], inside[i]);
}
