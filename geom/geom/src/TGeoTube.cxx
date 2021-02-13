// @(#)root/geom:$Id$
// Author: Andrei Gheata   24/10/01
// TGeoTube::Contains() and DistFromInside/In() implemented by Mihaela Gheata

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoTube
\ingroup Geometry_classes

Cylindrical tube  class. A tube has 3 parameters :
  - Rmin - minimum radius
  - Rmax - maximum radius
  - dz - half length
Begin_Macro(source)
{
   TCanvas *c = new TCanvas("c", "c",0,0,600,600);
   new TGeoManager("tube", "poza2");
   TGeoMaterial *mat = new TGeoMaterial("Al", 26.98,13,2.7);
   TGeoMedium *med = new TGeoMedium("MED",1,mat);
   TGeoVolume *top = gGeoManager->MakeBox("TOP",med,100,100,100);
   gGeoManager->SetTopVolume(top);
   TGeoVolume *vol = gGeoManager->MakeTube("TUBE",med, 20,30,40);
   vol->SetLineWidth(2);
   top->AddNode(vol,1);
   gGeoManager->CloseGeometry();
   gGeoManager->SetNsegments(80);
   top->Draw();
   TView *view = gPad->GetView();
   view->ShowAxis();
}
End_Macro
*/

/** \class TGeoTubeSeg
\ingroup Geometry_classes

A phi segment of a tube. Has 5 parameters :
  - the same 3 as a tube;
  - first phi limit (in degrees)
  - second phi limit
The segment will be be placed from the first angle (first phi limit)
to the second angle (second phi limit)
Begin_Macro(source)
{
   TCanvas *c = new TCanvas("c", "c",0,0,600,600);
   new TGeoManager("tubeseg", "poza3");
   TGeoMaterial *mat = new TGeoMaterial("Al", 26.98,13,2.7);
   TGeoMedium *med = new TGeoMedium("MED",1,mat);
   TGeoVolume *top = gGeoManager->MakeBox("TOP",med,100,100,100);
   gGeoManager->SetTopVolume(top);
   TGeoVolume *vol = gGeoManager->MakeTubs("TUBESEG",med, 20,30,40,-30,270);
   vol->SetLineWidth(2);
   top->AddNode(vol,1);
   gGeoManager->CloseGeometry();
   gGeoManager->SetNsegments(80);
   top->Draw();
   TView *view = gPad->GetView();
   view->ShowAxis();
}

End_Macro
*/

/** \class TGeoCtub
\ingroup Geometry_classes

A tube segment cut with 2 planes. Has 11 parameters :
  - the same 5 as a tube segment;
  - x, y, z components of the normal to the -dZ cut plane in
    point (0, 0, -dZ);
  - x, y, z components of the normal to the +dZ cut plane in
    point (0, 0, dZ);

Begin_Macro(source)
{
   TCanvas *c = new TCanvas("c", "c",0,0,600,600);
   new TGeoManager("ctub", "poza3");
   TGeoMaterial *mat = new TGeoMaterial("Al", 26.98,13,2.7);
   TGeoMedium *med = new TGeoMedium("MED",1,mat);
   TGeoVolume *top = gGeoManager->MakeBox("TOP",med,100,100,100);
   gGeoManager->SetTopVolume(top);
   Double_t theta = 160.*TMath::Pi()/180.;
   Double_t phi   = 30.*TMath::Pi()/180.;
   Double_t nlow[3];
   nlow[0] = TMath::Sin(theta)*TMath::Cos(phi);
   nlow[1] = TMath::Sin(theta)*TMath::Sin(phi);
   nlow[2] = TMath::Cos(theta);
   theta = 20.*TMath::Pi()/180.;
   phi   = 60.*TMath::Pi()/180.;
   Double_t nhi[3];
   nhi[0] = TMath::Sin(theta)*TMath::Cos(phi);
   nhi[1] = TMath::Sin(theta)*TMath::Sin(phi);
   nhi[2] = TMath::Cos(theta);
   TGeoVolume *vol = gGeoManager->MakeCtub("CTUB",med, 20,30,40,-30,250, nlow[0], nlow[1], nlow[2], nhi[0],nhi[1],nhi[2]);
   vol->SetLineWidth(2);
   top->AddNode(vol,1);
   gGeoManager->CloseGeometry();
   gGeoManager->SetNsegments(80);
   top->Draw();
   TView *view = gPad->GetView();
   view->ShowAxis();
}
End_Macro
*/

#include <iostream>

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TVirtualGeoPainter.h"
#include "TGeoTube.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"

ClassImp(TGeoTube);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoTube::TGeoTube()
{
   SetShapeBit(TGeoShape::kGeoTube);
   fRmin = 0.0;
   fRmax = 0.0;
   fDz   = 0.0;
}


////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius

TGeoTube::TGeoTube(Double_t rmin, Double_t rmax, Double_t dz)
           :TGeoBBox(0, 0, 0)
{
   SetShapeBit(TGeoShape::kGeoTube);
   SetTubeDimensions(rmin, rmax, dz);
   if ((fDz<0) || (fRmin<0) || (fRmax<0)) {
      SetShapeBit(kGeoRunTimeShape);
//      if (fRmax<=fRmin) SetShapeBit(kGeoInvalidShape);
//      printf("tube : dz=%f rmin=%f rmax=%f\n", dz, rmin, rmax);
   }
   ComputeBBox();
}
////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius

TGeoTube::TGeoTube(const char *name, Double_t rmin, Double_t rmax, Double_t dz)
           :TGeoBBox(name, 0, 0, 0)
{
   SetShapeBit(TGeoShape::kGeoTube);
   SetTubeDimensions(rmin, rmax, dz);
   if ((fDz<0) || (fRmin<0) || (fRmax<0)) {
      SetShapeBit(kGeoRunTimeShape);
//      if (fRmax<=fRmin) SetShapeBit(kGeoInvalidShape);
//      printf("tube : dz=%f rmin=%f rmax=%f\n", dz, rmin, rmax);
   }
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius
///  - param[0] = Rmin
///  - param[1] = Rmax
///  - param[2] = dz

TGeoTube::TGeoTube(Double_t *param)
         :TGeoBBox(0, 0, 0)
{
   SetShapeBit(TGeoShape::kGeoTube);
   SetDimensions(param);
   if ((fDz<0) || (fRmin<0) || (fRmax<0)) SetShapeBit(kGeoRunTimeShape);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoTube::~TGeoTube()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Computes capacity of the shape in [length^3]

Double_t TGeoTube::Capacity() const
{
   return TGeoTube::Capacity(fRmin,fRmax, fDz);
}

////////////////////////////////////////////////////////////////////////////////
/// Computes capacity of the shape in [length^3]

Double_t TGeoTube::Capacity(Double_t rmin, Double_t rmax, Double_t dz)
{
   Double_t capacity = 2.*TMath::Pi()*(rmax*rmax-rmin*rmin)*dz;
   return capacity;
}

////////////////////////////////////////////////////////////////////////////////
/// compute bounding box of the tube

void TGeoTube::ComputeBBox()
{
   fDX = fDY = fRmax;
   fDZ = fDz;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute normal to closest surface from POINT.

void TGeoTube::ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm)
{
   Double_t saf[3];
   Double_t rsq = point[0]*point[0]+point[1]*point[1];
   Double_t r = TMath::Sqrt(rsq);
   saf[0] = TMath::Abs(fDz-TMath::Abs(point[2]));
   saf[1] = (fRmin>1E-10)?TMath::Abs(r-fRmin):TGeoShape::Big();
   saf[2] = TMath::Abs(fRmax-r);
   Int_t i = TMath::LocMin(3,saf);
   if (i==0) {
      norm[0] = norm[1] = 0.;
      norm[2] = TMath::Sign(1.,dir[2]);
      return;
   }
   norm[2] = 0;
   Double_t phi = TMath::ATan2(point[1], point[0]);
   norm[0] = TMath::Cos(phi);
   norm[1] = TMath::Sin(phi);
   if (norm[0]*dir[0]+norm[1]*dir[1]<0) {
      norm[0] = -norm[0];
      norm[1] = -norm[1];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute normal to closest surface from POINT.

void TGeoTube::ComputeNormalS(const Double_t *point, const Double_t *dir, Double_t *norm,
                              Double_t /*rmin*/, Double_t /*rmax*/, Double_t /*dz*/)
{
   norm[2] = 0;
   Double_t phi = TMath::ATan2(point[1], point[0]);
   norm[0] = TMath::Cos(phi);
   norm[1] = TMath::Sin(phi);
   if (norm[0]*dir[0]+norm[1]*dir[1]<0) {
      norm[0] = -norm[0];
      norm[1] = -norm[1];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// test if point is inside this tube

Bool_t TGeoTube::Contains(const Double_t *point) const
{
   if (TMath::Abs(point[2]) > fDz) return kFALSE;
   Double_t r2 = point[0]*point[0]+point[1]*point[1];
   if ((r2<fRmin*fRmin) || (r2>fRmax*fRmax)) return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// compute closest distance from point px,py to each corner

Int_t TGeoTube::DistancetoPrimitive(Int_t px, Int_t py)
{
   Int_t n = gGeoManager->GetNsegments();
   Int_t numPoints = 4*n;
   if (!HasRmin()) numPoints = 2*(n+1);
   return ShapeDistancetoPrimitive(numPoints, px, py);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from inside point to surface of the tube (static)
/// Boundary safe algorithm.
/// compute distance to surface
/// Do Z

Double_t TGeoTube::DistFromInsideS(const Double_t *point, const Double_t *dir, Double_t rmin, Double_t rmax, Double_t dz)
{
   Double_t sz = TGeoShape::Big();
   if (dir[2]) {
      sz = (TMath::Sign(dz, dir[2])-point[2])/dir[2];
      if (sz<=0) return 0.0;
   }
   // Do R
   Double_t nsq=dir[0]*dir[0]+dir[1]*dir[1];
   if (TMath::Abs(nsq)<TGeoShape::Tolerance()) return sz;
   Double_t rsq=point[0]*point[0]+point[1]*point[1];
   Double_t rdotn=point[0]*dir[0]+point[1]*dir[1];
   Double_t b,d;
   Double_t sr = TGeoShape::Big();
   // inner cylinder
   if (rmin>0) {
      // Protection in case point is actually outside the tube
      if (rsq <= rmin*rmin+TGeoShape::Tolerance()) {
         if (rdotn<0) return 0.0;
      } else {
         if (rdotn<0) {
            DistToTube(rsq,nsq,rdotn,rmin,b,d);
            if (d>0) {
               sr=-b-d;
               if (sr>0) return TMath::Min(sz,sr);
            }
         }
      }
   }
   // outer cylinder
   if (rsq >= rmax*rmax-TGeoShape::Tolerance()) {
      if (rdotn>=0) return 0.0;
   }
   DistToTube(rsq,nsq,rdotn,rmax,b,d);
   if (d>0) {
      sr=-b+d;
      if (sr>0) return TMath::Min(sz,sr);
   }
   return 0.;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from inside point to surface of the tube
/// Boundary safe algorithm.

Double_t TGeoTube::DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
{
   if (iact<3 && safe) {
      *safe = Safety(point, kTRUE);
      if (iact==0) return TGeoShape::Big();
      if ((iact==1) && (*safe>step)) return TGeoShape::Big();
   }
   // compute distance to surface
   return DistFromInsideS(point, dir, fRmin, fRmax, fDz);
}

////////////////////////////////////////////////////////////////////////////////
/// Static method to compute distance from outside point to a tube with given parameters
/// Boundary safe algorithm.
/// check Z planes

Double_t TGeoTube::DistFromOutsideS(const Double_t *point, const Double_t *dir, Double_t rmin, Double_t rmax, Double_t dz)
{
   Double_t xi,yi,zi;
   Double_t rmaxsq = rmax*rmax;
   Double_t rminsq = rmin*rmin;
   zi = dz - TMath::Abs(point[2]);
   Double_t s = TGeoShape::Big();
   Bool_t in = kFALSE;
   Bool_t inz = (zi<0)?kFALSE:kTRUE;
   if (!inz) {
      if (point[2]*dir[2]>=0) return TGeoShape::Big();
      s  = -zi/TMath::Abs(dir[2]);
      xi = point[0]+s*dir[0];
      yi = point[1]+s*dir[1];
      Double_t r2=xi*xi+yi*yi;
      if ((rminsq<=r2) && (r2<=rmaxsq)) return s;
   }

   Double_t rsq = point[0]*point[0]+point[1]*point[1];
   // check outer cyl. surface
   Double_t nsq=dir[0]*dir[0]+dir[1]*dir[1];
   Double_t rdotn=point[0]*dir[0]+point[1]*dir[1];
   Double_t b,d;
   Bool_t inrmax = kFALSE;
   Bool_t inrmin = kFALSE;
   if (rsq<=rmaxsq+TGeoShape::Tolerance()) inrmax = kTRUE;
   if (rsq>=rminsq-TGeoShape::Tolerance()) inrmin = kTRUE;
   in = inz & inrmin & inrmax;
   // If inside, we are most likely on a boundary within machine precision.
   if (in) {
      Bool_t checkout = kFALSE;
      Double_t r = TMath::Sqrt(rsq);
      if (zi<rmax-r) {
         if ((TGeoShape::IsSameWithinTolerance(rmin,0)) || (zi<r-rmin)) {
            if (point[2]*dir[2]<0) return 0.0;
            return TGeoShape::Big();
         }
      }
      if ((rmaxsq-rsq) < (rsq-rminsq)) checkout = kTRUE;
      if (checkout) {
         if (rdotn>=0) return TGeoShape::Big();
         return 0.0;
      }
      if (TGeoShape::IsSameWithinTolerance(rmin,0)) return 0.0;
      if (rdotn>=0) return 0.0;
      // Ray exiting rmin -> check (+) solution for inner tube
      if (TMath::Abs(nsq)<TGeoShape::Tolerance()) return TGeoShape::Big();
      DistToTube(rsq, nsq, rdotn, rmin, b, d);
      if (d>0) {
         s=-b+d;
         if (s>0) {
            zi=point[2]+s*dir[2];
            if (TMath::Abs(zi)<=dz) return s;
         }
      }
      return TGeoShape::Big();
   }
   // Check outer cylinder (only r>rmax has to be considered)
   if (TMath::Abs(nsq)<TGeoShape::Tolerance()) return TGeoShape::Big();
   if (!inrmax) {
      DistToTube(rsq, nsq, rdotn, rmax, b, d);
      if (d>0) {
         s=-b-d;
         if (s>0) {
            zi=point[2]+s*dir[2];
            if (TMath::Abs(zi)<=dz) return s;
         }
      }
   }
   // check inner cylinder
   if (rmin>0) {
      DistToTube(rsq, nsq, rdotn, rmin, b, d);
      if (d>0) {
         s=-b+d;
         if (s>0) {
            zi=point[2]+s*dir[2];
            if (TMath::Abs(zi)<=dz) return s;
         }
      }
   }
   return TGeoShape::Big();
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from outside point to surface of the tube and safe distance
/// Boundary safe algorithm.
/// fist localize point w.r.t tube

Double_t TGeoTube::DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
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
   return DistFromOutsideS(point, dir, fRmin, fRmax, fDz);
}

////////////////////////////////////////////////////////////////////////////////
/// Static method computing the distance to a tube with given radius, starting from
/// POINT along DIR director cosines. The distance is computed as :
///    RSQ   = point[0]*point[0]+point[1]*point[1]
///    NSQ   = dir[0]*dir[0]+dir[1]*dir[1]  ---> should NOT be 0 !!!
///    RDOTN = point[0]*dir[0]+point[1]*dir[1]
/// The distance can be computed as :
///    D = -B +/- DELTA
/// where DELTA.GT.0 and D.GT.0

void TGeoTube::DistToTube(Double_t rsq, Double_t nsq, Double_t rdotn, Double_t radius, Double_t &b, Double_t &delta)
{
   Double_t t1 = 1./nsq;
   Double_t t3=rsq-(radius*radius);
   b          = t1*rdotn;
   Double_t c =t1*t3;
   delta = b*b-c;
   if (delta>0) {
      delta=TMath::Sqrt(delta);
   } else {
      delta = -1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Divide this tube shape belonging to volume "voldiv" into ndiv volumes
/// called divname, from start position with the given step. Returns pointer
/// to created division cell volume in case of Z divisions. For radial division
/// creates all volumes with different shapes and returns pointer to volume that
/// was divided. In case a wrong division axis is supplied, returns pointer to
/// volume that was divided.

TGeoVolume *TGeoTube::Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv,
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
         finder = new TGeoPatternCylR(voldiv, ndiv, start, end);
         vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
         voldiv->SetFinder(finder);
         finder->SetDivIndex(voldiv->GetNdaughters());
         for (id=0; id<ndiv; id++) {
            shape = new TGeoTube(start+id*step, start+(id+1)*step, fDz);
            vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
            vmulti->AddVolume(vol);
            opt = "R";
            voldiv->AddNodeOffset(vol, id, 0, opt.Data());
            ((TGeoNodeOffset*)voldiv->GetNodes()->At(voldiv->GetNdaughters()-1))->SetFinder(finder);
         }
         return vmulti;
      case 2:  //---                Phi division
         finder = new TGeoPatternCylPhi(voldiv, ndiv, start, end);
         voldiv->SetFinder(finder);
         finder->SetDivIndex(voldiv->GetNdaughters());
         shape = new TGeoTubeSeg(fRmin, fRmax, fDz, -step/2, step/2);
         vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
         vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
         vmulti->AddVolume(vol);
         opt = "Phi";
         for (id=0; id<ndiv; id++) {
            voldiv->AddNodeOffset(vol, id, start+id*step+step/2, opt.Data());
            ((TGeoNodeOffset*)voldiv->GetNodes()->At(voldiv->GetNdaughters()-1))->SetFinder(finder);
         }
         return vmulti;
      case 3: //---                  Z division
         finder = new TGeoPatternZ(voldiv, ndiv, start, start+ndiv*step);
         voldiv->SetFinder(finder);
         finder->SetDivIndex(voldiv->GetNdaughters());
         shape = new TGeoTube(fRmin, fRmax, step/2);
         vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
         vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
         vmulti->AddVolume(vol);
         opt = "Z";
         for (id=0; id<ndiv; id++) {
            voldiv->AddNodeOffset(vol, id, start+step/2+id*step, opt.Data());
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

const char *TGeoTube::GetAxisName(Int_t iaxis) const
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

Double_t TGeoTube::GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const
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
         xlo = 0;
         xhi = 360;
         dx = 360;
         return dx;
      case 3:
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

void TGeoTube::GetBoundingCylinder(Double_t *param) const
{
   param[0] = fRmin; // Rmin
   param[0] *= param[0];
   param[1] = fRmax; // Rmax
   param[1] *= param[1];
   param[2] = 0.;    // Phi1
   param[3] = 360.;  // Phi2
}

////////////////////////////////////////////////////////////////////////////////
/// in case shape has some negative parameters, these has to be computed
/// in order to fit the mother

TGeoShape *TGeoTube::GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix * /*mat*/) const
{
   if (!TestShapeBit(kGeoRunTimeShape)) return 0;
   Double_t rmin, rmax, dz;
   Double_t xmin,xmax;
   rmin = fRmin;
   rmax = fRmax;
   dz = fDz;
   if (fDz<0) {
      mother->GetAxisRange(3,xmin,xmax);
      if (xmax<0) return 0;
      dz=xmax;
   }
   mother->GetAxisRange(1,xmin,xmax);
   if (fRmin<0) {
      if (xmin<0) return 0;
      rmin = xmin;
   }
   if (fRmax<0) {
      if (xmax<=0) return 0;
      rmax = xmax;
   }

   return (new TGeoTube(GetName(), rmin, rmax, dz));
}

////////////////////////////////////////////////////////////////////////////////
/// print shape parameters

void TGeoTube::InspectShape() const
{
   printf("*** Shape %s: TGeoTube ***\n", GetName());
   printf("    Rmin = %11.5f\n", fRmin);
   printf("    Rmax = %11.5f\n", fRmax);
   printf("    dz   = %11.5f\n", fDz);
   printf(" Bounding box:\n");
   TGeoBBox::InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a TBuffer3D describing *this* shape.
/// Coordinates are in local reference frame.

TBuffer3D *TGeoTube::MakeBuffer3D() const
{
   Int_t n = gGeoManager->GetNsegments();
   Int_t nbPnts = 4*n;
   Int_t nbSegs = 8*n;
   Int_t nbPols = 4*n;
   if (!HasRmin()) {
      nbPnts = 2*(n+1);
      nbSegs = 5*n;
      nbPols = 3*n;
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

void TGeoTube::SetSegsAndPols(TBuffer3D &buffer) const
{
   Int_t i, j,indx;
   Int_t n = gGeoManager->GetNsegments();
   Int_t c = (((buffer.fColor) %8) -1) * 4;
   if (c < 0) c = 0;

   if (HasRmin()) {
      // circle segments:
      // lower rmin circle: i=0, (0, n-1)
      // lower rmax circle: i=1, (n, 2n-1)
      // upper rmin circle: i=2, (2n, 3n-1)
      // upper rmax circle: i=1, (3n, 4n-1)
      for (i = 0; i < 4; i++) {
         for (j = 0; j < n; j++) {
            indx = 3*(i*n+j);
            buffer.fSegs[indx  ] = c;
            buffer.fSegs[indx+1] = i*n+j;
            buffer.fSegs[indx+2] = i*n+(j+1)%n;
         }
      }
      // Z-parallel segments
      // inner: i=4, (4n, 5n-1)
      // outer: i=5, (5n, 6n-1)
      for (i = 4; i < 6; i++) {
         for (j = 0; j < n; j++) {
            indx = 3*(i*n+j);
            buffer.fSegs[indx  ] = c+1;
            buffer.fSegs[indx+1] = (i-4)*n+j;
            buffer.fSegs[indx+2] = (i-2)*n+j;
         }
      }
      // Radial segments
      // lower: i=6, (6n, 7n-1)
      // upper: i=7, (7n, 8n-1)
      for (i = 6; i < 8; i++) {
         for (j = 0; j < n; j++) {
            indx = 3*(i*n+j);
            buffer.fSegs[indx  ] = c;
            buffer.fSegs[indx+1] = 2*(i-6)*n+j;
            buffer.fSegs[indx+2] = (2*(i-6)+1)*n+j;
         }
      }
      // Polygons
      i=0;
      // Inner lateral (0, n-1)
      for (j = 0; j < n; j++) {
         indx = 6*(i*n+j);
         buffer.fPols[indx  ] = c;
         buffer.fPols[indx+1] = 4;
         buffer.fPols[indx+2] = j;
         buffer.fPols[indx+3] = 4*n+(j+1)%n;
         buffer.fPols[indx+4] = 2*n+j;
         buffer.fPols[indx+5] = 4*n+j;
      }
      i=1;
      // Outer lateral (n,2n-1)
      for (j = 0; j < n; j++) {
         indx = 6*(i*n+j);
         buffer.fPols[indx  ] = c+1;
         buffer.fPols[indx+1] = 4;
         buffer.fPols[indx+2] = n+j;
         buffer.fPols[indx+3] = 5*n+j;
         buffer.fPols[indx+4] = 3*n+j;
         buffer.fPols[indx+5] = 5*n+(j+1)%n;
      }
      i=2;
      // lower disc (2n, 3n-1)
      for (j = 0; j < n; j++) {
         indx = 6*(i*n+j);
         buffer.fPols[indx  ] = c;
         buffer.fPols[indx+1] = 4;
         buffer.fPols[indx+2] = j;
         buffer.fPols[indx+3] = 6*n+j;
         buffer.fPols[indx+4] = n+j;
         buffer.fPols[indx+5] = 6*n+(j+1)%n;
      }
      i=3;
      // upper disc (3n, 4n-1)
      for (j = 0; j < n; j++) {
         indx = 6*(i*n+j);
         buffer.fPols[indx  ] = c;
         buffer.fPols[indx+1] = 4;
         buffer.fPols[indx+2] = 2*n+j;
         buffer.fPols[indx+3] = 7*n+(j+1)%n;
         buffer.fPols[indx+4] = 3*n+j;
         buffer.fPols[indx+5] = 7*n+j;
      }
      return;
   }
   // Rmin=0 tubes
   // circle segments
   // lower rmax circle: i=0, (0, n-1)
   // upper rmax circle: i=1, (n, 2n-1)
   for (i = 0; i < 2; i++) {
      for (j = 0; j < n; j++) {
         indx = 3*(i*n+j);
         buffer.fSegs[indx  ] = c;
         buffer.fSegs[indx+1] = 2+i*n+j;
         buffer.fSegs[indx+2] = 2+i*n+(j+1)%n;
      }
   }
   // Z-parallel segments (2n,3n-1)
   for (j = 0; j < n; j++) {
      indx = 3*(2*n+j);
      buffer.fSegs[indx  ] = c+1;
      buffer.fSegs[indx+1] = 2+j;
      buffer.fSegs[indx+2] = 2+n+j;
   }
   // Radial segments
   // Lower circle: i=3, (3n,4n-1)
   // Upper circle: i=4, (4n,5n-1)
   for (i = 3; i < 5; i++) {
      for (j = 0; j < n; j++) {
         indx = 3*(i*n+j);
         buffer.fSegs[indx  ] = c;
         buffer.fSegs[indx+1] = i-3;
         buffer.fSegs[indx+2] = 2+(i-3)*n+j;
      }
   }
   // Polygons
   // lateral (0,n-1)
   for (j = 0; j < n; j++) {
      indx = 6*j;
      buffer.fPols[indx  ] = c+1;
      buffer.fPols[indx+1] = 4;
      buffer.fPols[indx+2] = j;
      buffer.fPols[indx+3] = 2*n+j;
      buffer.fPols[indx+4] = n+j;
      buffer.fPols[indx+5] = 2*n+(j+1)%n;
   }
   // bottom triangles (n,2n-1)
   for (j = 0; j < n; j++) {
      indx = 6*n + 5*j;
      buffer.fPols[indx  ] = c;
      buffer.fPols[indx+1] = 3;
      buffer.fPols[indx+2] = j;
      buffer.fPols[indx+3] = 3*n+(j+1)%n;
      buffer.fPols[indx+4] = 3*n+j;
   }
   // top triangles (2n,3n-1)
   for (j = 0; j < n; j++) {
      indx = 6*n + 5*n + 5*j;
      buffer.fPols[indx  ] = c;
      buffer.fPols[indx+1] = 3;
      buffer.fPols[indx+2] = n+j;
      buffer.fPols[indx+3] = 4*n+j;
      buffer.fPols[indx+4] = 4*n+(j+1)%n;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// computes the closest distance from given point to this shape, according
/// to option. The matching point on the shape is stored in spoint.

Double_t TGeoTube::Safety(const Double_t *point, Bool_t in) const
{
#ifndef NEVER
   Double_t r = TMath::Sqrt(point[0]*point[0]+point[1]*point[1]);
   Double_t safe, safrmin, safrmax;
   if (in) {
      safe    = fDz-TMath::Abs(point[2]); // positive if inside
      if (fRmin>1E-10) {
         safrmin = r-fRmin;
         if (safrmin < safe) safe = safrmin;
      }
      safrmax = fRmax-r;
      if (safrmax < safe) safe = safrmax;
   } else {
      safe    = -fDz+TMath::Abs(point[2]);
      if (fRmin>1E-10) {
         safrmin = -r+fRmin;
         if (safrmin > safe) safe = safrmin;
      }
      safrmax = -fRmax+r;
      if (safrmax > safe) safe = safrmax;
   }
   return safe;
#else
   Double_t saf[3];
   Double_t rsq = point[0]*point[0]+point[1]*point[1];
   Double_t r = TMath::Sqrt(rsq);
   saf[0] = fDz-TMath::Abs(point[2]); // positive if inside
   saf[1] = (fRmin>1E-10)?(r-fRmin):TGeoShape::Big();
   saf[2] = fRmax-r;
   if (in) return saf[TMath::LocMin(3,saf)];
   for (Int_t i=0; i<3; i++) saf[i]=-saf[i];
   return saf[TMath::LocMax(3,saf)];
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// computes the closest distance from given point to this shape, according
/// to option. The matching point on the shape is stored in spoint.

Double_t TGeoTube::SafetyS(const Double_t *point, Bool_t in, Double_t rmin, Double_t rmax, Double_t dz, Int_t skipz)
{
   Double_t saf[3];
   Double_t rsq = point[0]*point[0]+point[1]*point[1];
   Double_t r = TMath::Sqrt(rsq);
   switch (skipz) {
      case 1: // skip lower Z plane
         saf[0] = dz - point[2];
         break;
      case 2: // skip upper Z plane
         saf[0] = dz + point[2];
         break;
      case 3: // skip both
         saf[0] = TGeoShape::Big();
         break;
      default:
         saf[0] = dz-TMath::Abs(point[2]);
   }
   saf[1] = (rmin>1E-10)?(r-rmin):TGeoShape::Big();
   saf[2] = rmax-r;
//   printf("saf0=%g saf1=%g saf2=%g in=%d skipz=%d\n", saf[0],saf[1],saf[2],in,skipz);
   if (in) return saf[TMath::LocMin(3,saf)];
   for (Int_t i=0; i<3; i++) saf[i]=-saf[i];
   return saf[TMath::LocMax(3,saf)];
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoTube::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive)) return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   rmin = " << fRmin << ";" << std::endl;
   out << "   rmax = " << fRmax << ";" << std::endl;
   out << "   dz   = " << fDz << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoTube(\"" << GetName() << "\",rmin,rmax,dz);" << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Set tube dimensions.

void TGeoTube::SetTubeDimensions(Double_t rmin, Double_t rmax, Double_t dz)
{
   fRmin = rmin;
   fRmax = rmax;
   fDz   = dz;
   if (fRmin>0 && fRmax>0 && fRmin>=fRmax)
      Error("SetTubeDimensions", "In shape %s wrong rmin=%g rmax=%g", GetName(), rmin,rmax);
}

////////////////////////////////////////////////////////////////////////////////
/// Set tube dimensions starting from a list.

void TGeoTube::SetDimensions(Double_t *param)
{
   Double_t rmin = param[0];
   Double_t rmax = param[1];
   Double_t dz   = param[2];
   SetTubeDimensions(rmin, rmax, dz);
}

////////////////////////////////////////////////////////////////////////////////
/// Fills array with n random points located on the line segments of the shape mesh.
/// The output array must be provided with a length of minimum 3*npoints. Returns
/// true if operation is implemented.

Bool_t TGeoTube::GetPointsOnSegments(Int_t npoints, Double_t *array) const
{
   if (npoints > (npoints/2)*2) {
      Error("GetPointsOnSegments","Npoints must be even number");
      return kFALSE;
   }
   Int_t nc = 0;
   if (HasRmin()) nc = (Int_t)TMath::Sqrt(0.5*npoints);
   else           nc = (Int_t)TMath::Sqrt(1.*npoints);
   Double_t dphi = TMath::TwoPi()/nc;
   Double_t phi = 0;
   Int_t ntop = 0;
   if (HasRmin()) ntop = npoints/2 - nc*(nc-1);
   else           ntop = npoints - nc*(nc-1);
   Double_t dz = 2*fDz/(nc-1);
   Double_t z = 0;
   Int_t icrt = 0;
   Int_t nphi = nc;
   // loop z sections
   for (Int_t i=0; i<nc; i++) {
      if (i == (nc-1)) nphi = ntop;
      z = -fDz + i*dz;
      // loop points on circle sections
      for (Int_t j=0; j<nphi; j++) {
         phi = j*dphi;
         if (HasRmin()) {
            array[icrt++] = fRmin * TMath::Cos(phi);
            array[icrt++] = fRmin * TMath::Sin(phi);
            array[icrt++] = z;
         }
         array[icrt++] = fRmax * TMath::Cos(phi);
         array[icrt++] = fRmax * TMath::Sin(phi);
         array[icrt++] = z;
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// create tube mesh points

void TGeoTube::SetPoints(Double_t *points) const
{
   Double_t dz;
   Int_t j, n;
   n = gGeoManager->GetNsegments();
   Double_t dphi = 360./n;
   Double_t phi = 0;
   dz = fDz;
   Int_t indx = 0;
   if (points) {
      if (HasRmin()) {
         // 4*n points
         // (0,n-1) lower rmin circle
         // (2n, 3n-1) upper rmin circle
         for (j = 0; j < n; j++) {
            phi = j*dphi*TMath::DegToRad();
            points[indx+6*n] = points[indx] = fRmin * TMath::Cos(phi);
            indx++;
            points[indx+6*n] = points[indx] = fRmin * TMath::Sin(phi);
            indx++;
            points[indx+6*n] = dz;
            points[indx]     =-dz;
            indx++;
         }
         // (n, 2n-1) lower rmax circle
         // (3n, 4n-1) upper rmax circle
         for (j = 0; j < n; j++) {
            phi = j*dphi*TMath::DegToRad();
            points[indx+6*n] = points[indx] = fRmax * TMath::Cos(phi);
            indx++;
            points[indx+6*n] = points[indx] = fRmax * TMath::Sin(phi);
            indx++;
            points[indx+6*n]= dz;
            points[indx]    =-dz;
            indx++;
         }
      } else {
         // centers of lower/upper circles (0,1)
         points[indx++] = 0.;
         points[indx++] = 0.;
         points[indx++] = -dz;
         points[indx++] = 0.;
         points[indx++] = 0.;
         points[indx++] = dz;
         // lower rmax circle (2, 2+n-1)
         // upper rmax circle (2+n, 2+2n-1)
         for (j = 0; j < n; j++) {
            phi = j*dphi*TMath::DegToRad();
            points[indx+3*n] = points[indx] = fRmax * TMath::Cos(phi);
            indx++;
            points[indx+3*n] = points[indx] = fRmax * TMath::Sin(phi);
            indx++;
            points[indx+3*n]= dz;
            points[indx]    =-dz;
            indx++;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// create tube mesh points

void TGeoTube::SetPoints(Float_t *points) const
{
   Double_t dz;
   Int_t j, n;
   n = gGeoManager->GetNsegments();
   Double_t dphi = 360./n;
   Double_t phi = 0;
   dz = fDz;
   Int_t indx = 0;
   if (points) {
      if (HasRmin()) {
         // 4*n points
         // (0,n-1) lower rmin circle
         // (2n, 3n-1) upper rmin circle
         for (j = 0; j < n; j++) {
            phi = j*dphi*TMath::DegToRad();
            points[indx+6*n] = points[indx] = fRmin * TMath::Cos(phi);
            indx++;
            points[indx+6*n] = points[indx] = fRmin * TMath::Sin(phi);
            indx++;
            points[indx+6*n] = dz;
            points[indx]     =-dz;
            indx++;
         }
         // (n, 2n-1) lower rmax circle
         // (3n, 4n-1) upper rmax circle
         for (j = 0; j < n; j++) {
            phi = j*dphi*TMath::DegToRad();
            points[indx+6*n] = points[indx] = fRmax * TMath::Cos(phi);
            indx++;
            points[indx+6*n] = points[indx] = fRmax * TMath::Sin(phi);
            indx++;
            points[indx+6*n]= dz;
            points[indx]    =-dz;
            indx++;
         }
      } else {
         // centers of lower/upper circles (0,1)
         points[indx++] = 0.;
         points[indx++] = 0.;
         points[indx++] = -dz;
         points[indx++] = 0.;
         points[indx++] = 0.;
         points[indx++] = dz;
         // lower rmax circle (2, 2+n-1)
         // upper rmax circle (2+n, 2+2n-1)
         for (j = 0; j < n; j++) {
            phi = j*dphi*TMath::DegToRad();
            points[indx+3*n] = points[indx] = fRmax * TMath::Cos(phi);
            indx++;
            points[indx+3*n] = points[indx] = fRmax * TMath::Sin(phi);
            indx++;
            points[indx+3*n]= dz;
            points[indx]    =-dz;
            indx++;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of vertices of the mesh representation

Int_t TGeoTube::GetNmeshVertices() const
{
   Int_t n = gGeoManager->GetNsegments();
   Int_t numPoints = n*4;
   if (!HasRmin()) numPoints = 2*(n+1);
   return numPoints;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoTube::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   Int_t n = gGeoManager->GetNsegments();
   nvert = n*4;
   nsegs = n*8;
   npols = n*4;
   if (!HasRmin()) {
      nvert = 2*(n+1);
      nsegs = 5*n;
      npols = 3*n;
   } else {
      nvert = n*4;
      nsegs = n*8;
      npols = n*4;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// fill size of this 3-D object

void TGeoTube::Sizeof3D() const
{
}

////////////////////////////////////////////////////////////////////////////////
/// Fills a static 3D buffer and returns a reference.

const TBuffer3D & TGeoTube::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
{
   static TBuffer3DTube buffer;
   TGeoBBox::FillBuffer3D(buffer, reqSections, localFrame);

   if (reqSections & TBuffer3D::kShapeSpecific) {
      buffer.fRadiusInner  = fRmin;
      buffer.fRadiusOuter  = fRmax;
      buffer.fHalfLength   = fDz;
      buffer.SetSectionsValid(TBuffer3D::kShapeSpecific);
   }
   if (reqSections & TBuffer3D::kRawSizes) {
      Int_t n = gGeoManager->GetNsegments();
      Int_t nbPnts = 4*n;
      Int_t nbSegs = 8*n;
      Int_t nbPols = 4*n;
      if (!HasRmin()) {
         nbPnts = 2*(n+1);
         nbSegs = 5*n;
         nbPols = 3*n;
      }
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

void TGeoTube::Contains_v(const Double_t *points, Bool_t *inside, Int_t vecsize) const
{
   for (Int_t i=0; i<vecsize; i++) inside[i] = Contains(&points[3*i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the normal for an array o points so that norm.dot.dir is positive
/// Input: Arrays of point coordinates and directions + vector size
/// Output: Array of normal directions

void TGeoTube::ComputeNormal_v(const Double_t *points, const Double_t *dirs, Double_t *norms, Int_t vecsize)
{
   for (Int_t i=0; i<vecsize; i++) ComputeNormal(&points[3*i], &dirs[3*i], &norms[3*i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoTube::DistFromInside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromInside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoTube::DistFromOutside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromOutside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute safe distance from each of the points in the input array.
/// Input: Array of point coordinates, array of statuses for these points, size of the arrays
/// Output: Safety values

void TGeoTube::Safety_v(const Double_t *points, const Bool_t *inside, Double_t *safe, Int_t vecsize) const
{
   for (Int_t i=0; i<vecsize; i++) safe[i] = Safety(&points[3*i], inside[i]);
}

ClassImp(TGeoTubeSeg);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoTubeSeg::TGeoTubeSeg()
            :TGeoTube(),
             fPhi1(0.), fPhi2(0.), fS1(0.), fC1(0.), fS2(0.), fC2(0.), fSm(0.), fCm(0.), fCdfi(0.)
{
   SetShapeBit(TGeoShape::kGeoTubeSeg);
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius.
/// The segment will be from phiStart to phiEnd expressed in degree.

TGeoTubeSeg::TGeoTubeSeg(Double_t rmin, Double_t rmax, Double_t dz,
                          Double_t phiStart, Double_t phiEnd)
            :TGeoTube(rmin, rmax, dz),
             fPhi1(0.), fPhi2(0.), fS1(0.), fC1(0.), fS2(0.), fC2(0.), fSm(0.), fCm(0.), fCdfi(0.)
{
   SetShapeBit(TGeoShape::kGeoTubeSeg);
   SetTubsDimensions(rmin, rmax, dz, phiStart, phiEnd);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius
/// The segment will be from phiStart to phiEnd expressed in degree.

TGeoTubeSeg::TGeoTubeSeg(const char *name, Double_t rmin, Double_t rmax, Double_t dz,
                          Double_t phiStart, Double_t phiEnd)
            :TGeoTube(name, rmin, rmax, dz)
{
   SetShapeBit(TGeoShape::kGeoTubeSeg);
   SetTubsDimensions(rmin, rmax, dz, phiStart, phiEnd);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius
///  - param[0] = Rmin
///  - param[1] = Rmax
///  - param[2] = dz
///  - param[3] = phi1
///  - param[4] = phi2

TGeoTubeSeg::TGeoTubeSeg(Double_t *param)
            :TGeoTube(0, 0, 0)
{
   SetShapeBit(TGeoShape::kGeoTubeSeg);
   SetDimensions(param);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoTubeSeg::~TGeoTubeSeg()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Function called after streaming an object of this class.

void TGeoTubeSeg::AfterStreamer()
{
   InitTrigonometry();
}

////////////////////////////////////////////////////////////////////////////////
/// Init frequently used trigonometric values

void TGeoTubeSeg::InitTrigonometry()
{
   Double_t phi1 = fPhi1*TMath::DegToRad();
   Double_t phi2 = fPhi2*TMath::DegToRad();
   fC1 = TMath::Cos(phi1);
   fS1 = TMath::Sin(phi1);
   fC2 = TMath::Cos(phi2);
   fS2 = TMath::Sin(phi2);
   Double_t fio = 0.5*(phi1+phi2);
   fCm = TMath::Cos(fio);
   fSm = TMath::Sin(fio);
   Double_t dfi = 0.5*(phi2-phi1);
   fCdfi = TMath::Cos(dfi);
}

////////////////////////////////////////////////////////////////////////////////
/// Computes capacity of the shape in [length^3]

Double_t TGeoTubeSeg::Capacity() const
{
   return TGeoTubeSeg::Capacity(fRmin,fRmax,fDz,fPhi1,fPhi2);
}

////////////////////////////////////////////////////////////////////////////////
/// Computes capacity of the shape in [length^3]

Double_t TGeoTubeSeg::Capacity(Double_t rmin, Double_t rmax, Double_t dz, Double_t phiStart, Double_t phiEnd)
{
   Double_t capacity = TMath::Abs(phiEnd-phiStart)*TMath::DegToRad()*(rmax*rmax-rmin*rmin)*dz;
   return capacity;
}

////////////////////////////////////////////////////////////////////////////////
/// compute bounding box of the tube segment

void TGeoTubeSeg::ComputeBBox()
{
   Double_t xc[4];
   Double_t yc[4];
   xc[0] = fRmax*fC1;
   yc[0] = fRmax*fS1;
   xc[1] = fRmax*fC2;
   yc[1] = fRmax*fS2;
   xc[2] = fRmin*fC1;
   yc[2] = fRmin*fS1;
   xc[3] = fRmin*fC2;
   yc[3] = fRmin*fS2;

   Double_t xmin = xc[TMath::LocMin(4, &xc[0])];
   Double_t xmax = xc[TMath::LocMax(4, &xc[0])];
   Double_t ymin = yc[TMath::LocMin(4, &yc[0])];
   Double_t ymax = yc[TMath::LocMax(4, &yc[0])];

   Double_t dp = fPhi2-fPhi1;
   if (dp<0) dp+=360;
   Double_t ddp = -fPhi1;
   if (ddp<0) ddp+= 360;
   if (ddp>360) ddp-=360;
   if (ddp<=dp) xmax = fRmax;
   ddp = 90-fPhi1;
   if (ddp<0) ddp+= 360;
   if (ddp>360) ddp-=360;
   if (ddp<=dp) ymax = fRmax;
   ddp = 180-fPhi1;
   if (ddp<0) ddp+= 360;
   if (ddp>360) ddp-=360;
   if (ddp<=dp) xmin = -fRmax;
   ddp = 270-fPhi1;
   if (ddp<0) ddp+= 360;
   if (ddp>360) ddp-=360;
   if (ddp<=dp) ymin = -fRmax;
   fOrigin[0] = (xmax+xmin)/2;
   fOrigin[1] = (ymax+ymin)/2;
   fOrigin[2] = 0;
   fDX = (xmax-xmin)/2;
   fDY = (ymax-ymin)/2;
   fDZ = fDz;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute normal to closest surface from POINT.

void TGeoTubeSeg::ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm)
{
   Double_t saf[3];
   Double_t rsq = point[0]*point[0]+point[1]*point[1];
   Double_t r = TMath::Sqrt(rsq);
   saf[0] = TMath::Abs(fDz-TMath::Abs(point[2]));
   saf[1] = (fRmin>1E-10)?TMath::Abs(r-fRmin):TGeoShape::Big();
   saf[2] = TMath::Abs(fRmax-r);
   Int_t i = TMath::LocMin(3,saf);
   if (((fPhi2-fPhi1)<360.) && TGeoShape::IsCloseToPhi(saf[i], point,fC1,fS1,fC2,fS2)) {
      TGeoShape::NormalPhi(point,dir,norm,fC1,fS1,fC2,fS2);
      return;
   }
   if (i==0) {
      norm[0] = norm[1] = 0.;
      norm[2] = TMath::Sign(1.,dir[2]);
      return;
   };
   norm[2] = 0;
   Double_t phi = TMath::ATan2(point[1], point[0]);
   norm[0] = TMath::Cos(phi);
   norm[1] = TMath::Sin(phi);
   if (norm[0]*dir[0]+norm[1]*dir[1]<0) {
      norm[0] = -norm[0];
      norm[1] = -norm[1];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute normal to closest surface from POINT.

void TGeoTubeSeg::ComputeNormalS(const Double_t *point, const Double_t *dir, Double_t *norm,
                                 Double_t rmin, Double_t rmax, Double_t /*dz*/,
                                 Double_t c1, Double_t s1, Double_t c2, Double_t s2)
{
   Double_t saf[2];
   Double_t rsq = point[0]*point[0]+point[1]*point[1];
   Double_t r = TMath::Sqrt(rsq);
   saf[0] = (rmin>1E-10)?TMath::Abs(r-rmin):TGeoShape::Big();
   saf[1] = TMath::Abs(rmax-r);
   Int_t i = TMath::LocMin(2,saf);
   if (TGeoShape::IsCloseToPhi(saf[i], point,c1,s1,c2,s2)) {
      TGeoShape::NormalPhi(point,dir,norm,c1,s1,c2,s2);
      return;
   }
   norm[2] = 0;
   Double_t phi = TMath::ATan2(point[1], point[0]);
   norm[0] = TMath::Cos(phi);
   norm[1] = TMath::Sin(phi);
   if (norm[0]*dir[0]+norm[1]*dir[1]<0) {
      norm[0] = -norm[0];
      norm[1] = -norm[1];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// test if point is inside this tube segment
/// first check if point is inside the tube

Bool_t TGeoTubeSeg::Contains(const Double_t *point) const
{
   if (!TGeoTube::Contains(point)) return kFALSE;
   return IsInPhiRange(point, fPhi1, fPhi2);
}

////////////////////////////////////////////////////////////////////////////////
/// compute closest distance from point px,py to each corner

Int_t TGeoTubeSeg::DistancetoPrimitive(Int_t px, Int_t py)
{
   Int_t n = gGeoManager->GetNsegments()+1;
   const Int_t numPoints = 4*n;
   return ShapeDistancetoPrimitive(numPoints, px, py);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from inside point to surface of the tube segment (static)
/// Boundary safe algorithm.
/// Do Z

Double_t TGeoTubeSeg::DistFromInsideS(const Double_t *point, const Double_t *dir, Double_t rmin, Double_t rmax, Double_t dz,
                                 Double_t c1, Double_t s1, Double_t c2, Double_t s2, Double_t cm, Double_t sm, Double_t cdfi)
{
   Double_t stube = TGeoTube::DistFromInsideS(point,dir,rmin,rmax,dz);
   if (stube<=0) return 0.0;
   Double_t sfmin = TGeoShape::Big();
   Double_t rsq = point[0]*point[0]+point[1]*point[1];
   Double_t r = TMath::Sqrt(rsq);
   Double_t cpsi=point[0]*cm+point[1]*sm;
   if (cpsi>r*cdfi+TGeoShape::Tolerance())  {
      sfmin = TGeoShape::DistToPhiMin(point, dir, s1, c1, s2, c2, sm, cm);
      return TMath::Min(stube,sfmin);
   }
   // Point on the phi boundary or outside
   // which one: phi1 or phi2
   Double_t ddotn, xi, yi;
   if (TMath::Abs(point[1]-s1*r) < TMath::Abs(point[1]-s2*r)) {
      ddotn = s1*dir[0]-c1*dir[1];
      if (ddotn>=0) return 0.0;
      ddotn = -s2*dir[0]+c2*dir[1];
      if (ddotn<=0) return stube;
      sfmin = s2*point[0]-c2*point[1];
      if (sfmin<=0) return stube;
      sfmin /= ddotn;
      if (sfmin >= stube) return stube;
      xi = point[0]+sfmin*dir[0];
      yi = point[1]+sfmin*dir[1];
      if (yi*cm-xi*sm<0) return stube;
      return sfmin;
   }
   ddotn = -s2*dir[0]+c2*dir[1];
   if (ddotn>=0) return 0.0;
   ddotn = s1*dir[0]-c1*dir[1];
   if (ddotn<=0) return stube;
   sfmin = -s1*point[0]+c1*point[1];
   if (sfmin<=0) return stube;
   sfmin /= ddotn;
   if (sfmin >= stube) return stube;
   xi = point[0]+sfmin*dir[0];
   yi = point[1]+sfmin*dir[1];
   if (yi*cm-xi*sm>0) return stube;
   return sfmin;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from inside point to surface of the tube segment
/// Boundary safe algorithm.

Double_t TGeoTubeSeg::DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
{
   if (iact<3 && safe) {
      *safe = SafetyS(point, kTRUE, fRmin, fRmax, fDz, fPhi1, fPhi2);
      if (iact==0) return TGeoShape::Big();
      if ((iact==1) && (*safe>step)) return TGeoShape::Big();
   }
   if ((fPhi2-fPhi1)>=360.) return TGeoTube::DistFromInsideS(point,dir,fRmin,fRmax,fDz);

   // compute distance to surface
   return TGeoTubeSeg::DistFromInsideS(point,dir,fRmin,fRmax,fDz,fC1,fS1,fC2,fS2,fCm,fSm,fCdfi);
}

////////////////////////////////////////////////////////////////////////////////
/// Static method to compute distance to arbitrary tube segment from outside point
/// Boundary safe algorithm.

Double_t TGeoTubeSeg::DistFromOutsideS(const Double_t *point, const Double_t *dir, Double_t rmin, Double_t rmax,
                                Double_t dz, Double_t c1, Double_t s1, Double_t c2, Double_t s2,
                                Double_t cm, Double_t sm, Double_t cdfi)
{
   Double_t r2, cpsi;
   // check Z planes
   Double_t xi, yi, zi;
   zi = dz - TMath::Abs(point[2]);
   Double_t rmaxsq = rmax*rmax;
   Double_t rminsq = rmin*rmin;
   Double_t s = TGeoShape::Big();
   Double_t snxt=TGeoShape::Big();
   Bool_t in = kFALSE;
   Bool_t inz = (zi<0)?kFALSE:kTRUE;
   if (!inz) {
      if (point[2]*dir[2]>=0) return TGeoShape::Big();
      s = -zi/TMath::Abs(dir[2]);
      xi = point[0]+s*dir[0];
      yi = point[1]+s*dir[1];
      r2=xi*xi+yi*yi;
      if ((rminsq<=r2) && (r2<=rmaxsq)) {
         cpsi=(xi*cm+yi*sm)/TMath::Sqrt(r2);
         if (cpsi>=cdfi) return s;
      }
   }

   // check outer cyl. surface
   Double_t rsq = point[0]*point[0]+point[1]*point[1];
   Double_t r = TMath::Sqrt(rsq);
   Double_t nsq=dir[0]*dir[0]+dir[1]*dir[1];
   Double_t rdotn=point[0]*dir[0]+point[1]*dir[1];
   Double_t b,d;
   Bool_t inrmax = kFALSE;
   Bool_t inrmin = kFALSE;
   Bool_t inphi  = kFALSE;
   if (rsq<=rmaxsq+TGeoShape::Tolerance()) inrmax = kTRUE;
   if (rsq>=rminsq-TGeoShape::Tolerance()) inrmin = kTRUE;
   cpsi=point[0]*cm+point[1]*sm;
   if (cpsi>r*cdfi-TGeoShape::Tolerance())  inphi = kTRUE;
   in = inz & inrmin & inrmax & inphi;
   // If inside, we are most likely on a boundary within machine precision.
   if (in) {
      Bool_t checkout = kFALSE;
      Double_t safphi = (cpsi-r*cdfi)*TMath::Sqrt(1.-cdfi*cdfi);
//      Double_t sch, cch;
      // check if on Z boundaries
      if (zi<rmax-r) {
         if (TGeoShape::IsSameWithinTolerance(rmin,0) || (zi<r-rmin)) {
            if (zi<safphi) {
               if (point[2]*dir[2]<0) return 0.0;
               return TGeoShape::Big();
            }
         }
      }
      if ((rmaxsq-rsq) < (rsq-rminsq)) checkout = kTRUE;
      // check if on Rmax boundary
      if (checkout && (rmax-r<safphi)) {
         if (rdotn>=0) return TGeoShape::Big();
         return 0.0;
      }
      if (TMath::Abs(nsq)<TGeoShape::Tolerance()) return TGeoShape::Big();
      // check if on phi boundary
      if (TGeoShape::IsSameWithinTolerance(rmin,0) || (safphi<r-rmin)) {
         // We may cross again a phi of rmin boundary
         // check first if we are on phi1 or phi2
         Double_t un;
         if (point[0]*c1 + point[1]*s1 > point[0]*c2 + point[1]*s2) {
            un = dir[0]*s1-dir[1]*c1;
            if (un < 0) return 0.0;
            if (cdfi>=0) return TGeoShape::Big();
            un = -dir[0]*s2+dir[1]*c2;
            if (un<0) {
               s = -point[0]*s2+point[1]*c2;
               if (s>0) {
                  s /= (-un);
                  zi = point[2]+s*dir[2];
                  if (TMath::Abs(zi)<=dz) {
                     xi = point[0]+s*dir[0];
                     yi = point[1]+s*dir[1];
                     r2=xi*xi+yi*yi;
                     if ((rminsq<=r2) && (r2<=rmaxsq)) {
                        if ((yi*cm-xi*sm)>0) return s;
                     }
                  }
               }
            }
         } else {
            un = -dir[0]*s2+dir[1]*c2;
            if (un < 0) return 0.0;
            if (cdfi>=0) return TGeoShape::Big();
            un = dir[0]*s1-dir[1]*c1;
            if (un<0) {
               s = point[0]*s1-point[1]*c1;
               if (s>0) {
                  s /= (-un);
                  zi = point[2]+s*dir[2];
                  if (TMath::Abs(zi)<=dz) {
                     xi = point[0]+s*dir[0];
                     yi = point[1]+s*dir[1];
                     r2=xi*xi+yi*yi;
                     if ((rminsq<=r2) && (r2<=rmaxsq)) {
                        if ((yi*cm-xi*sm)<0) return s;
                     }
                  }
               }
            }
         }
         // We may also cross rmin, (+) solution
         if (rdotn>=0) return TGeoShape::Big();
         if (cdfi>=0) return TGeoShape::Big();
         DistToTube(rsq, nsq, rdotn, rmin, b, d);
         if (d>0) {
            s=-b+d;
            if (s>0) {
               zi=point[2]+s*dir[2];
               if (TMath::Abs(zi)<=dz) {
                  xi=point[0]+s*dir[0];
                  yi=point[1]+s*dir[1];
                  if ((xi*cm+yi*sm) >= rmin*cdfi) return s;
               }
            }
         }
         return TGeoShape::Big();
      }
      // we are on rmin boundary: we may cross again rmin or a phi facette
      if (rdotn>=0) return 0.0;
      DistToTube(rsq, nsq, rdotn, rmin, b, d);
      if (d>0) {
         s=-b+d;
         if (s>0) {
            zi=point[2]+s*dir[2];
            if (TMath::Abs(zi)<=dz) {
               // now check phi range
               xi=point[0]+s*dir[0];
               yi=point[1]+s*dir[1];
               r2=xi*xi+yi*yi;
               if ((xi*cm+yi*sm) >= rmin*cdfi) return s;
               // now we really have to check any phi crossing
               Double_t un=-dir[0]*s1+dir[1]*c1;
               if (un > 0) {
                  s=point[0]*s1-point[1]*c1;
                  if (s>=0) {
                     s /= un;
                     zi=point[2]+s*dir[2];
                     if (TMath::Abs(zi)<=dz) {
                        xi=point[0]+s*dir[0];
                        yi=point[1]+s*dir[1];
                        r2=xi*xi+yi*yi;
                        if ((rminsq<=r2) && (r2<=rmaxsq)) {
                           if ((yi*cm-xi*sm)<=0) {
                              if (s<snxt) snxt=s;
                           }
                        }
                     }
                  }
               }
               un=dir[0]*s2-dir[1]*c2;
               if (un > 0) {
                  s=(point[1]*c2-point[0]*s2)/un;
                  if (s>=0 && s<snxt) {
                     zi=point[2]+s*dir[2];
                     if (TMath::Abs(zi)<=dz) {
                        xi=point[0]+s*dir[0];
                        yi=point[1]+s*dir[1];
                        r2=xi*xi+yi*yi;
                        if ((rminsq<=r2) && (r2<=rmaxsq)) {
                           if ((yi*cm-xi*sm)>=0) {
                              return s;
                           }
                        }
                     }
                  }
               }
               return snxt;
            }
         }
      }
      return TGeoShape::Big();
   }
   // only r>rmax has to be considered
   if (TMath::Abs(nsq)<TGeoShape::Tolerance()) return TGeoShape::Big();
   if (rsq>=rmax*rmax) {
      if (rdotn>=0) return TGeoShape::Big();
      TGeoTube::DistToTube(rsq, nsq, rdotn, rmax, b, d);
      if (d>0) {
         s=-b-d;
         if (s>0) {
            zi=point[2]+s*dir[2];
            if (TMath::Abs(zi)<=dz) {
               xi=point[0]+s*dir[0];
               yi=point[1]+s*dir[1];
               cpsi = xi*cm+yi*sm;
               if (cpsi>=rmax*cdfi) return s;
            }
         }
      }
   }
   // check inner cylinder
   if (rmin>0) {
      TGeoTube::DistToTube(rsq, nsq, rdotn, rmin, b, d);
      if (d>0) {
         s=-b+d;
         if (s>0) {
            zi=point[2]+s*dir[2];
            if (TMath::Abs(zi)<=dz) {
               xi=point[0]+s*dir[0];
               yi=point[1]+s*dir[1];
               cpsi = xi*cm+yi*sm;
               if (cpsi>=rmin*cdfi) snxt=s;
            }
         }
      }
   }
   // check phi planes
   Double_t un=-dir[0]*s1+dir[1]*c1;
   if (un > 0) {
      s=point[0]*s1-point[1]*c1;
      if (s>=0) {
         s /= un;
         zi=point[2]+s*dir[2];
         if (TMath::Abs(zi)<=dz) {
            xi=point[0]+s*dir[0];
            yi=point[1]+s*dir[1];
            r2=xi*xi+yi*yi;
            if ((rminsq<=r2) && (r2<=rmaxsq)) {
               if ((yi*cm-xi*sm)<=0) {
                  if (s<snxt) snxt=s;
               }
            }
         }
      }
   }
   un=dir[0]*s2-dir[1]*c2;
   if (un > 0) {
      s=point[1]*c2-point[0]*s2;
      if (s>=0) {
         s /= un;
         zi=point[2]+s*dir[2];
         if (TMath::Abs(zi)<=dz) {
            xi=point[0]+s*dir[0];
            yi=point[1]+s*dir[1];
            r2=xi*xi+yi*yi;
            if ((rminsq<=r2) && (r2<=rmaxsq)) {
               if ((yi*cm-xi*sm)>=0) {
                  if (s<snxt) snxt=s;
               }
            }
         }
      }
   }
   return snxt;
}

////////////////////////////////////////////////////////////////////////////////
/// compute distance from outside point to surface of the tube segment
/// fist localize point w.r.t tube

Double_t TGeoTubeSeg::DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
{
   if (iact<3 && safe) {
      *safe = SafetyS(point, kFALSE, fRmin, fRmax, fDz, fPhi1, fPhi2);
      if (iact==0) return TGeoShape::Big();
      if ((iact==1) && (step<=*safe)) return TGeoShape::Big();
   }
// Check if the bounding box is crossed within the requested distance
   Double_t sdist = TGeoBBox::DistFromOutside(point,dir, fDX, fDY, fDZ, fOrigin, step);
   if (sdist>=step) return TGeoShape::Big();
   if ((fPhi2-fPhi1)>=360.) return TGeoTube::DistFromOutsideS(point,dir,fRmin,fRmax,fDz);

   // find distance to shape
   return TGeoTubeSeg::DistFromOutsideS(point, dir, fRmin, fRmax, fDz, fC1, fS1, fC2, fS2, fCm, fSm, fCdfi);
}

////////////////////////////////////////////////////////////////////////////////
/// Divide this tube segment shape belonging to volume "voldiv" into ndiv volumes
/// called divname, from start position with the given step. Returns pointer
/// to created division cell volume in case of Z divisions. For radialdivision
/// creates all volumes with different shapes and returns pointer to volume that
/// was divided. In case a wrong division axis is supplied, returns pointer to
/// volume that was divided.

TGeoVolume *TGeoTubeSeg::Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv,
                             Double_t start, Double_t step)
{
   TGeoShape *shape;           //--- shape to be created
   TGeoVolume *vol;            //--- division volume to be created
   TGeoVolumeMulti *vmulti;    //--- generic divided volume
   TGeoPatternFinder *finder;  //--- finder to be attached
   TString opt = "";           //--- option to be attached
   Double_t dphi;
   Int_t id;
   Double_t end = start+ndiv*step;
   switch (iaxis) {
      case 1:  //---                 R division
         finder = new TGeoPatternCylR(voldiv, ndiv, start, end);
         vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
         voldiv->SetFinder(finder);
         finder->SetDivIndex(voldiv->GetNdaughters());
         for (id=0; id<ndiv; id++) {
            shape = new TGeoTubeSeg(start+id*step, start+(id+1)*step, fDz, fPhi1, fPhi2);
            vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
            vmulti->AddVolume(vol);
            opt = "R";
            voldiv->AddNodeOffset(vol, id, 0, opt.Data());
            ((TGeoNodeOffset*)voldiv->GetNodes()->At(voldiv->GetNdaughters()-1))->SetFinder(finder);
         }
         return vmulti;
      case 2:  //---                 Phi division
         dphi = fPhi2-fPhi1;
         if (dphi<0) dphi+=360.;
         if (step<=0) {step=dphi/ndiv; start=fPhi1; end=fPhi2;}
         finder = new TGeoPatternCylPhi(voldiv, ndiv, start, end);
         voldiv->SetFinder(finder);
         finder->SetDivIndex(voldiv->GetNdaughters());
         shape = new TGeoTubeSeg(fRmin, fRmax, fDz, -step/2, step/2);
         vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
         vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
         vmulti->AddVolume(vol);
         opt = "Phi";
         for (id=0; id<ndiv; id++) {
            voldiv->AddNodeOffset(vol, id, start+id*step+step/2, opt.Data());
            ((TGeoNodeOffset*)voldiv->GetNodes()->At(voldiv->GetNdaughters()-1))->SetFinder(finder);
         }
         return vmulti;
      case 3: //---                  Z division
         finder = new TGeoPatternZ(voldiv, ndiv, start, end);
         voldiv->SetFinder(finder);
         finder->SetDivIndex(voldiv->GetNdaughters());
         shape = new TGeoTubeSeg(fRmin, fRmax, step/2, fPhi1, fPhi2);
         vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
         vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
         vmulti->AddVolume(vol);
         opt = "Z";
         for (id=0; id<ndiv; id++) {
            voldiv->AddNodeOffset(vol, id, start+step/2+id*step, opt.Data());
            ((TGeoNodeOffset*)voldiv->GetNodes()->At(voldiv->GetNdaughters()-1))->SetFinder(finder);
         }
         return vmulti;
      default:
         Error("Divide", "In shape %s wrong axis type for division", GetName());
         return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get range of shape for a given axis.

Double_t TGeoTubeSeg::GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const
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
         xlo = -fDz;
         xhi = fDz;
         dx = xhi-xlo;
         return dx;
   }
   return dx;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill vector param[4] with the bounding cylinder parameters. The order
/// is the following : Rmin, Rmax, Phi1, Phi2

void TGeoTubeSeg::GetBoundingCylinder(Double_t *param) const
{
   param[0] = fRmin;
   param[0] *= param[0];
   param[1] = fRmax;
   param[1] *= param[1];
   param[2] = fPhi1;
   param[3] = fPhi2;
}

////////////////////////////////////////////////////////////////////////////////
/// in case shape has some negative parameters, these has to be computed
/// in order to fit the mother

TGeoShape *TGeoTubeSeg::GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix * /*mat*/) const
{
   if (!TestShapeBit(kGeoRunTimeShape)) return 0;
   if (!mother->TestShapeBit(kGeoTube)) {
      Error("GetMakeRuntimeShape", "Invalid mother for shape %s", GetName());
      return 0;
   }
   Double_t rmin, rmax, dz;
   rmin = fRmin;
   rmax = fRmax;
   dz = fDz;
   if (fDz<0) dz=((TGeoTube*)mother)->GetDz();
   if (fRmin<0)
      rmin = ((TGeoTube*)mother)->GetRmin();
   if ((fRmax<0) || (fRmax<=fRmin))
      rmax = ((TGeoTube*)mother)->GetRmax();

   return (new TGeoTubeSeg(GetName(),rmin, rmax, dz, fPhi1, fPhi2));
}

////////////////////////////////////////////////////////////////////////////////
/// print shape parameters

void TGeoTubeSeg::InspectShape() const
{
   printf("*** Shape %s: TGeoTubeSeg ***\n", GetName());
   printf("    Rmin = %11.5f\n", fRmin);
   printf("    Rmax = %11.5f\n", fRmax);
   printf("    dz   = %11.5f\n", fDz);
   printf("    phi1 = %11.5f\n", fPhi1);
   printf("    phi2 = %11.5f\n", fPhi2);
   printf(" Bounding box:\n");
   TGeoBBox::InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a TBuffer3D describing *this* shape.
/// Coordinates are in local reference frame.

TBuffer3D *TGeoTubeSeg::MakeBuffer3D() const
{
   Int_t n = gGeoManager->GetNsegments()+1;
   Int_t nbPnts = 4*n;
   Int_t nbSegs = 2*nbPnts;
   Int_t nbPols = nbPnts-2;

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

void TGeoTubeSeg::SetSegsAndPols(TBuffer3D &buff) const
{
   Int_t i, j;
   Int_t n = gGeoManager->GetNsegments()+1;
   Int_t c = GetBasicColor();

   memset(buff.fSegs, 0, buff.NbSegs()*3*sizeof(Int_t));
   for (i = 0; i < 4; i++) {
      for (j = 1; j < n; j++) {
         buff.fSegs[(i*n+j-1)*3  ] = c;
         buff.fSegs[(i*n+j-1)*3+1] = i*n+j-1;
         buff.fSegs[(i*n+j-1)*3+2] = i*n+j;
      }
   }
   for (i = 4; i < 6; i++) {
      for (j = 0; j < n; j++) {
         buff.fSegs[(i*n+j)*3  ] = c+1;
         buff.fSegs[(i*n+j)*3+1] = (i-4)*n+j;
         buff.fSegs[(i*n+j)*3+2] = (i-2)*n+j;
      }
   }
   for (i = 6; i < 8; i++) {
      for (j = 0; j < n; j++) {
         buff.fSegs[(i*n+j)*3  ] = c;
         buff.fSegs[(i*n+j)*3+1] = 2*(i-6)*n+j;
         buff.fSegs[(i*n+j)*3+2] = (2*(i-6)+1)*n+j;
      }
   }

   Int_t indx = 0;
   memset(buff.fPols, 0, buff.NbPols()*6*sizeof(Int_t));
   i = 0;
   for (j = 0; j < n-1; j++) {
      buff.fPols[indx++] = c;
      buff.fPols[indx++] = 4;
      buff.fPols[indx++] = (4+i)*n+j+1;
      buff.fPols[indx++] = (2+i)*n+j;
      buff.fPols[indx++] = (4+i)*n+j;
      buff.fPols[indx++] = i*n+j;
   }
   i = 1;
   for (j = 0; j < n-1; j++) {
      buff.fPols[indx++] = c;
      buff.fPols[indx++] = 4;
      buff.fPols[indx++] = i*n+j;
      buff.fPols[indx++] = (4+i)*n+j;
      buff.fPols[indx++] = (2+i)*n+j;
      buff.fPols[indx++] = (4+i)*n+j+1;
   }
   i = 2;
   for (j = 0; j < n-1; j++) {
      buff.fPols[indx++] = c+i;
      buff.fPols[indx++] = 4;
      buff.fPols[indx++] = (i-2)*2*n+j;
      buff.fPols[indx++] = (4+i)*n+j;
      buff.fPols[indx++] = ((i-2)*2+1)*n+j;
      buff.fPols[indx++] = (4+i)*n+j+1;
   }
   i = 3;
   for (j = 0; j < n-1; j++) {
      buff.fPols[indx++] = c+i;
      buff.fPols[indx++] = 4;
      buff.fPols[indx++] = (4+i)*n+j+1;
      buff.fPols[indx++] = ((i-2)*2+1)*n+j;
      buff.fPols[indx++] = (4+i)*n+j;
      buff.fPols[indx++] = (i-2)*2*n+j;
   }
   buff.fPols[indx++] = c+2;
   buff.fPols[indx++] = 4;
   buff.fPols[indx++] = 6*n;
   buff.fPols[indx++] = 4*n;
   buff.fPols[indx++] = 7*n;
   buff.fPols[indx++] = 5*n;
   buff.fPols[indx++] = c+2;
   buff.fPols[indx++] = 4;
   buff.fPols[indx++] = 6*n-1;
   buff.fPols[indx++] = 8*n-1;
   buff.fPols[indx++] = 5*n-1;
   buff.fPols[indx++] = 7*n-1;
}

////////////////////////////////////////////////////////////////////////////////
/// computes the closest distance from given point InitTrigonometry();to this shape, according
/// to option. The matching point on the shape is stored in spoint.

Double_t TGeoTubeSeg::Safety(const Double_t *point, Bool_t in) const
{
   Double_t safe = TGeoShape::Big();
   Double_t saf[3];
   Double_t rsq = point[0]*point[0]+point[1]*point[1];
   Double_t r = TMath::Sqrt(rsq);
   if (in) {
      saf[0] = fDz-TMath::Abs(point[2]);
      saf[1] = r-fRmin;
      saf[2] = fRmax-r;
      safe   = saf[TMath::LocMin(3,saf)];
      if ((fPhi2-fPhi1)>=360.) return safe;
      Double_t safphi = TGeoShape::SafetyPhi(point,in,fPhi1,fPhi2);
      return TMath::Min(safe, safphi);
   }
   // Point expected to be outside
   Bool_t inphi  = kFALSE;
   Double_t cpsi=point[0]*fCm+point[1]*fSm;
   saf[0] = TMath::Abs(point[2])-fDz;
   if (cpsi>r*fCdfi-TGeoShape::Tolerance())  inphi = kTRUE;
   if (inphi) {
      saf[1] = fRmin-r;
      saf[2] = r-fRmax;
      safe = saf[TMath::LocMax(3,saf)];
      safe = TMath::Max(0., safe);
      return safe;
   }
   // Point outside the phi range
   // Compute projected radius of the (r,phi) position vector onto
   // phi1 and phi2 edges and take the maximum for choosing the side.
   Double_t rproj = TMath::Max(point[0]*fC1+point[1]*fS1, point[0]*fC2+point[1]*fS2);
   saf[1] = fRmin-rproj;
   saf[2] = rproj-fRmax;
   safe = TMath::Max(saf[1], saf[2]);
   if ((fPhi2-fPhi1)>=360.) return TMath::Max(safe,saf[0]);
   if (safe>0) {
      // rproj not within (rmin,rmax) - > no need to calculate safphi
      safe = TMath::Sqrt(rsq-rproj*rproj+safe*safe);
      return (saf[0]<0) ? safe : TMath::Sqrt(safe*safe+saf[0]*saf[0]);
   }
   Double_t safphi = TGeoShape::SafetyPhi(point,in,fPhi1,fPhi2);
   return (saf[0]<0) ? safphi : TMath::Sqrt(saf[0]*saf[0]+safphi*safphi);
}

////////////////////////////////////////////////////////////////////////////////
/// Static method to compute the closest distance from given point to this shape.

Double_t TGeoTubeSeg::SafetyS(const Double_t *point, Bool_t in, Double_t rmin, Double_t rmax, Double_t dz,
                              Double_t phi1d, Double_t phi2d, Int_t skipz)
{
   Double_t safe = TGeoShape::Big();
   Double_t saf[3];
   Double_t rsq = point[0]*point[0]+point[1]*point[1];
   Double_t r = TMath::Sqrt(rsq);

   switch (skipz) {
      case 1: // skip lower Z plane
         saf[0] = dz - point[2];
         break;
      case 2: // skip upper Z plane
         saf[0] = dz + point[2];
         break;
      case 3: // skip both
         saf[0] = TGeoShape::Big();
         break;
      default:
         saf[0] = dz-TMath::Abs(point[2]);
   }

   if (in) {
      saf[1] = r-rmin;
      saf[2] = rmax-r;
      safe   = saf[TMath::LocMin(3,saf)];
      if ((phi2d-phi1d)>=360.) return safe;
      Double_t safphi = TGeoShape::SafetyPhi(point,in,phi1d,phi2d);
      return TMath::Min(safe, safphi);
   }
   // Point expected to be outside
   saf[0] = -saf[0];
   Bool_t inphi  = kFALSE;
   Double_t phi1 = phi1d*TMath::DegToRad();
   Double_t phi2 = phi2d*TMath::DegToRad();

   Double_t fio = 0.5*(phi1+phi2);
   Double_t cm = TMath::Cos(fio);
   Double_t sm = TMath::Sin(fio);
   Double_t cpsi=point[0]*cm+point[1]*sm;
   Double_t dfi = 0.5*(phi2-phi1);
   Double_t cdfi = TMath::Cos(dfi);
   if (cpsi>r*cdfi-TGeoShape::Tolerance())  inphi = kTRUE;
   if (inphi) {
      saf[1] = rmin-r;
      saf[2] = r-rmax;
      safe = saf[TMath::LocMax(3,saf)];
      safe = TMath::Max(0., safe);
      return safe;
   }
   // Point outside the phi range
   // Compute projected radius of the (r,phi) position vector onto
   // phi1 and phi2 edges and take the maximum for choosing the side.
   Double_t c1 = TMath::Cos(phi1);
   Double_t s1 = TMath::Sin(phi1);
   Double_t c2 = TMath::Cos(phi2);
   Double_t s2 = TMath::Sin(phi2);

   Double_t rproj = TMath::Max(point[0]*c1+point[1]*s1, point[0]*c2+point[1]*s2);
   saf[1] = rmin-rproj;
   saf[2] = rproj-rmax;
   safe   = TMath::Max(saf[1], saf[2]);
   if ((phi2d-phi1d)>=360.) return TMath::Max(safe,saf[0]);
   if (safe>0) {
      // rproj not within (rmin,rmax) - > no need to calculate safphi
      safe = TMath::Sqrt(rsq-rproj*rproj+safe*safe);
      return (saf[0]<0) ? safe : TMath::Sqrt(safe*safe+saf[0]*saf[0]);
   }
   Double_t safphi = TGeoShape::SafetyPhi(point,in,phi1d,phi2d);
   return (saf[0]<0) ? safphi : TMath::Sqrt(saf[0]*saf[0]+safphi*safphi);
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoTubeSeg::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive)) return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   rmin = " << fRmin << ";" << std::endl;
   out << "   rmax = " << fRmax << ";" << std::endl;
   out << "   dz   = " << fDz << ";" << std::endl;
   out << "   phi1 = " << fPhi1 << ";" << std::endl;
   out << "   phi2 = " << fPhi2 << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoTubeSeg(\"" << GetName() << "\",rmin,rmax,dz,phi1,phi2);" << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Set dimensions of the tube segment.
/// The segment will be from phiStart to phiEnd expressed in degree.

void TGeoTubeSeg::SetTubsDimensions(Double_t rmin, Double_t rmax, Double_t dz,
                                    Double_t phiStart, Double_t phiEnd)
{
   fRmin = rmin;
   fRmax = rmax;
   fDz   = dz;
   fPhi1 = phiStart;
   if (fPhi1 < 0) fPhi1 += 360.;
   fPhi2 = phiEnd;
   while (fPhi2<=fPhi1) fPhi2+=360.;
   if (TGeoShape::IsSameWithinTolerance(fPhi1,fPhi2)) Fatal("SetTubsDimensions", "In shape %s invalid phi1=%g, phi2=%g\n", GetName(), fPhi1, fPhi2);
   InitTrigonometry();
}

////////////////////////////////////////////////////////////////////////////////
/// Set dimensions of the tube segment starting from a list.

void TGeoTubeSeg::SetDimensions(Double_t *param)
{
   Double_t rmin = param[0];
   Double_t rmax = param[1];
   Double_t dz   = param[2];
   Double_t phi1 = param[3];
   Double_t phi2 = param[4];
   SetTubsDimensions(rmin, rmax, dz, phi1, phi2);
}

////////////////////////////////////////////////////////////////////////////////
/// Fills array with n random points located on the line segments of the shape mesh.
/// The output array must be provided with a length of minimum 3*npoints. Returns
/// true if operation is implemented.

Bool_t TGeoTubeSeg::GetPointsOnSegments(Int_t npoints, Double_t *array) const
{
   if (npoints > (npoints/2)*2) {
      Error("GetPointsOnSegments","Npoints must be even number");
      return kFALSE;
   }
   Int_t nc = (Int_t)TMath::Sqrt(0.5*npoints);
   Double_t dphi = (fPhi2-fPhi1)*TMath::DegToRad()/(nc-1);
   Double_t phi = 0;
   Double_t phi1 = fPhi1 * TMath::DegToRad();
   Int_t ntop = npoints/2 - nc*(nc-1);
   Double_t dz = 2*fDz/(nc-1);
   Double_t z = 0;
   Int_t icrt = 0;
   Int_t nphi = nc;
   // loop z sections
   for (Int_t i=0; i<nc; i++) {
      if (i == (nc-1)) {
         nphi = ntop;
         dphi = (fPhi2-fPhi1)*TMath::DegToRad()/(nphi-1);
      }
      z = -fDz + i*dz;
      // loop points on circle sections
      for (Int_t j=0; j<nphi; j++) {
         phi = phi1 + j*dphi;
         array[icrt++] = fRmin * TMath::Cos(phi);
         array[icrt++] = fRmin * TMath::Sin(phi);
         array[icrt++] = z;
         array[icrt++] = fRmax * TMath::Cos(phi);
         array[icrt++] = fRmax * TMath::Sin(phi);
         array[icrt++] = z;
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create tube segment mesh points.

void TGeoTubeSeg::SetPoints(Double_t *points) const
{
   Double_t dz;
   Int_t j, n;
   Double_t phi, phi1, phi2, dphi;
   phi1 = fPhi1;
   phi2 = fPhi2;
   if (phi2<phi1) phi2+=360.;
   n = gGeoManager->GetNsegments()+1;

   dphi = (phi2-phi1)/(n-1);
   dz   = fDz;

   if (points) {
      Int_t indx = 0;

      for (j = 0; j < n; j++) {
         phi = (phi1+j*dphi)*TMath::DegToRad();
         points[indx+6*n] = points[indx] = fRmin * TMath::Cos(phi);
         indx++;
         points[indx+6*n] = points[indx] = fRmin * TMath::Sin(phi);
         indx++;
         points[indx+6*n] = dz;
         points[indx]     =-dz;
         indx++;
      }
      for (j = 0; j < n; j++) {
         phi = (phi1+j*dphi)*TMath::DegToRad();
         points[indx+6*n] = points[indx] = fRmax * TMath::Cos(phi);
         indx++;
         points[indx+6*n] = points[indx] = fRmax * TMath::Sin(phi);
         indx++;
         points[indx+6*n]= dz;
         points[indx]    =-dz;
         indx++;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create tube segment mesh points.

void TGeoTubeSeg::SetPoints(Float_t *points) const
{
   Double_t dz;
   Int_t j, n;
   Double_t phi, phi1, phi2, dphi;
   phi1 = fPhi1;
   phi2 = fPhi2;
   if (phi2<phi1) phi2+=360.;
   n = gGeoManager->GetNsegments()+1;

   dphi = (phi2-phi1)/(n-1);
   dz   = fDz;

   if (points) {
      Int_t indx = 0;

      for (j = 0; j < n; j++) {
         phi = (phi1+j*dphi)*TMath::DegToRad();
         points[indx+6*n] = points[indx] = fRmin * TMath::Cos(phi);
         indx++;
         points[indx+6*n] = points[indx] = fRmin * TMath::Sin(phi);
         indx++;
         points[indx+6*n] = dz;
         points[indx]     =-dz;
         indx++;
      }
      for (j = 0; j < n; j++) {
         phi = (phi1+j*dphi)*TMath::DegToRad();
         points[indx+6*n] = points[indx] = fRmax * TMath::Cos(phi);
         indx++;
         points[indx+6*n] = points[indx] = fRmax * TMath::Sin(phi);
         indx++;
         points[indx+6*n]= dz;
         points[indx]    =-dz;
         indx++;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoTubeSeg::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   Int_t n = gGeoManager->GetNsegments()+1;
   nvert = n*4;
   nsegs = n*8;
   npols = n*4 - 2;
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of vertices of the mesh representation

Int_t TGeoTubeSeg::GetNmeshVertices() const
{
   Int_t n = gGeoManager->GetNsegments()+1;
   Int_t numPoints = n*4;
   return numPoints;
}

////////////////////////////////////////////////////////////////////////////////
/// fill size of this 3-D object

void TGeoTubeSeg::Sizeof3D() const
{
}

////////////////////////////////////////////////////////////////////////////////
/// Fills a static 3D buffer and returns a reference.

const TBuffer3D & TGeoTubeSeg::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
{
   static TBuffer3DTubeSeg buffer;
   TGeoBBox::FillBuffer3D(buffer, reqSections, localFrame);

   if (reqSections & TBuffer3D::kShapeSpecific) {
      // These from TBuffer3DTube / TGeoTube
      buffer.fRadiusInner  = fRmin;
      buffer.fRadiusOuter  = fRmax;
      buffer.fHalfLength   = fDz;
      buffer.fPhiMin       = fPhi1;
      buffer.fPhiMax       = fPhi2;
      buffer.SetSectionsValid(TBuffer3D::kShapeSpecific);
   }
   if (reqSections & TBuffer3D::kRawSizes) {
      Int_t n = gGeoManager->GetNsegments()+1;
      Int_t nbPnts = 4*n;
      Int_t nbSegs = 2*nbPnts;
      Int_t nbPols = nbPnts-2;
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

void TGeoTubeSeg::Contains_v(const Double_t *points, Bool_t *inside, Int_t vecsize) const
{
   for (Int_t i=0; i<vecsize; i++) inside[i] = Contains(&points[3*i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the normal for an array o points so that norm.dot.dir is positive
/// Input: Arrays of point coordinates and directions + vector size
/// Output: Array of normal directions

void TGeoTubeSeg::ComputeNormal_v(const Double_t *points, const Double_t *dirs, Double_t *norms, Int_t vecsize)
{
   for (Int_t i=0; i<vecsize; i++) ComputeNormal(&points[3*i], &dirs[3*i], &norms[3*i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoTubeSeg::DistFromInside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromInside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoTubeSeg::DistFromOutside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromOutside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute safe distance from each of the points in the input array.
/// Input: Array of point coordinates, array of statuses for these points, size of the arrays
/// Output: Safety values

void TGeoTubeSeg::Safety_v(const Double_t *points, const Bool_t *inside, Double_t *safe, Int_t vecsize) const
{
   for (Int_t i=0; i<vecsize; i++) safe[i] = Safety(&points[3*i], inside[i]);
}

ClassImp(TGeoCtub);

TGeoCtub::TGeoCtub()
{
// default ctor
   fNlow[0] = fNlow[1] = fNhigh[0] = fNhigh[1] = 0.;
   fNlow[2] = -1;
   fNhigh[2] = 1;
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoCtub::TGeoCtub(Double_t rmin, Double_t rmax, Double_t dz, Double_t phi1, Double_t phi2,
                   Double_t lx, Double_t ly, Double_t lz, Double_t tx, Double_t ty, Double_t tz)
         :TGeoTubeSeg(rmin, rmax, dz, phi1, phi2)
{
   fNlow[0] = lx;
   fNlow[1] = ly;
   fNlow[2] = lz;
   fNhigh[0] = tx;
   fNhigh[1] = ty;
   fNhigh[2] = tz;
   SetShapeBit(kGeoCtub);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoCtub::TGeoCtub(const char *name, Double_t rmin, Double_t rmax, Double_t dz, Double_t phi1, Double_t phi2,
                   Double_t lx, Double_t ly, Double_t lz, Double_t tx, Double_t ty, Double_t tz)
         :TGeoTubeSeg(name, rmin, rmax, dz, phi1, phi2)
{
   fNlow[0] = lx;
   fNlow[1] = ly;
   fNlow[2] = lz;
   fNhigh[0] = tx;
   fNhigh[1] = ty;
   fNhigh[2] = tz;
   SetShapeBit(kGeoCtub);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// ctor with parameters

TGeoCtub::TGeoCtub(Double_t *params)
         :TGeoTubeSeg(0,0,0,0,0)
{
   SetCtubDimensions(params[0], params[1], params[2], params[3], params[4], params[5],
                     params[6], params[7], params[8], params[9], params[10]);
   SetShapeBit(kGeoCtub);
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoCtub::~TGeoCtub()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Computes capacity of the shape in [length^3]

Double_t TGeoCtub::Capacity() const
{
   Double_t capacity = TGeoTubeSeg::Capacity();
   return capacity;
}

////////////////////////////////////////////////////////////////////////////////
/// compute minimum bounding box of the ctub

void TGeoCtub::ComputeBBox()
{
   TGeoTubeSeg::ComputeBBox();
   if ((fNlow[2]>-(1E-10)) || (fNhigh[2]<1E-10)) {
      Error("ComputeBBox", "In shape %s wrong definition of cut planes", GetName());
      return;
   }
   Double_t xc=0, yc=0;
   Double_t zmin=0, zmax=0;
   Double_t z1;
   Double_t z[8];
   // check if nxy is in the phi range
   Double_t phi_low = TMath::ATan2(fNlow[1], fNlow[0]) *TMath::RadToDeg();
   Double_t phi_hi = TMath::ATan2(fNhigh[1], fNhigh[0]) *TMath::RadToDeg();
   Bool_t in_range_low = kFALSE;
   Bool_t in_range_hi = kFALSE;

   Int_t i;
   for (i=0; i<2; i++) {
      if (phi_low<0) phi_low+=360.;
      Double_t dphi = fPhi2 -fPhi1;
      if (dphi < 0) dphi+=360.;
      Double_t ddp = phi_low-fPhi1;
      if (ddp<0) ddp += 360.;
      if (ddp <= dphi) {
         xc = fRmin*TMath::Cos(phi_low*TMath::DegToRad());
         yc = fRmin*TMath::Sin(phi_low*TMath::DegToRad());
         z1 = GetZcoord(xc, yc, -fDz);
         xc = fRmax*TMath::Cos(phi_low*TMath::DegToRad());
         yc = fRmax*TMath::Sin(phi_low*TMath::DegToRad());
         z1 = TMath::Min(z1, GetZcoord(xc, yc, -fDz));
         if (in_range_low)
            zmin = TMath::Min(zmin, z1);
         else
            zmin = z1;
         in_range_low = kTRUE;
      }
      phi_low += 180;
      if (phi_low>360) phi_low-=360.;
   }

   for (i=0; i<2; i++) {
      if (phi_hi<0) phi_hi+=360.;
      Double_t dphi = fPhi2 -fPhi1;
      if (dphi < 0) dphi+=360.;
      Double_t ddp = phi_hi-fPhi1;
      if (ddp<0) ddp += 360.;
      if (ddp <= dphi) {
         xc = fRmin*TMath::Cos(phi_hi*TMath::DegToRad());
         yc = fRmin*TMath::Sin(phi_hi*TMath::DegToRad());
         z1 = GetZcoord(xc, yc, fDz);
         xc = fRmax*TMath::Cos(phi_hi*TMath::DegToRad());
         yc = fRmax*TMath::Sin(phi_hi*TMath::DegToRad());
         z1 = TMath::Max(z1, GetZcoord(xc, yc, fDz));
         if (in_range_hi)
            zmax = TMath::Max(zmax, z1);
         else
            zmax = z1;
         in_range_hi = kTRUE;
      }
      phi_hi += 180;
      if (phi_hi>360) phi_hi-=360.;
   }


   xc = fRmin*fC1;
   yc = fRmin*fS1;
   z[0] = GetZcoord(xc, yc, -fDz);
   z[4] = GetZcoord(xc, yc, fDz);

   xc = fRmin*fC2;
   yc = fRmin*fS2;
   z[1] = GetZcoord(xc, yc, -fDz);
   z[5] = GetZcoord(xc, yc, fDz);

   xc = fRmax*fC1;
   yc = fRmax*fS1;
   z[2] = GetZcoord(xc, yc, -fDz);
   z[6] = GetZcoord(xc, yc, fDz);

   xc = fRmax*fC2;
   yc = fRmax*fS2;
   z[3] = GetZcoord(xc, yc, -fDz);
   z[7] = GetZcoord(xc, yc, fDz);

   z1 = z[TMath::LocMin(4, &z[0])];
   if (in_range_low)
      zmin = TMath::Min(zmin, z1);
   else
      zmin = z1;

   z1 = z[TMath::LocMax(4, &z[4])+4];
   if (in_range_hi)
      zmax = TMath::Max(zmax, z1);
   else
      zmax = z1;

   fDZ = 0.5*(zmax-zmin);
   fOrigin[2] = 0.5*(zmax+zmin);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute normal to closest surface from POINT.

void TGeoCtub::ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm)
{
   Double_t saf[4];
   Bool_t isseg = kTRUE;
   if (TMath::Abs(fPhi2-fPhi1-360.)<1E-8) isseg=kFALSE;
   Double_t rsq = point[0]*point[0]+point[1]*point[1];
   Double_t r = TMath::Sqrt(rsq);

   saf[0] = TMath::Abs(point[0]*fNlow[0] + point[1]*fNlow[1] + (fDz+point[2])*fNlow[2]);
   saf[1] = TMath::Abs(point[0]*fNhigh[0] + point[1]*fNhigh[1] - (fDz-point[2])*fNhigh[2]);
   saf[2] = (fRmin>1E-10)?TMath::Abs(r-fRmin):TGeoShape::Big();
   saf[3] = TMath::Abs(fRmax-r);
   Int_t i = TMath::LocMin(4,saf);
   if (isseg) {
      if (TGeoShape::IsCloseToPhi(saf[i], point,fC1,fS1,fC2,fS2)) {
         TGeoShape::NormalPhi(point,dir,norm,fC1,fS1,fC2,fS2);
         return;
      }
   }
   if (i==0) {
      memcpy(norm, fNlow, 3*sizeof(Double_t));
      if (norm[0]*dir[0]+norm[1]*dir[1]+norm[2]*dir[2]<0) {
         norm[0] = -norm[0];
         norm[1] = -norm[1];
         norm[2] = -norm[2];
      }
      return;
   }
   if (i==1) {
      memcpy(norm, fNhigh, 3*sizeof(Double_t));
      if (norm[0]*dir[0]+norm[1]*dir[1]+norm[2]*dir[2]<0) {
         norm[0] = -norm[0];
         norm[1] = -norm[1];
         norm[2] = -norm[2];
      }
      return;
   }

   norm[2] = 0;
   Double_t phi = TMath::ATan2(point[1], point[0]);
   norm[0] = TMath::Cos(phi);
   norm[1] = TMath::Sin(phi);
   if (norm[0]*dir[0]+norm[1]*dir[1]<0) {
      norm[0] = -norm[0];
      norm[1] = -norm[1];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// check if point is contained in the cut tube
/// check the lower cut plane

Bool_t TGeoCtub::Contains(const Double_t *point) const
{
   Double_t zin = point[0]*fNlow[0]+point[1]*fNlow[1]+(point[2]+fDz)*fNlow[2];
   if (zin>0) return kFALSE;
   // check the higher cut plane
   zin = point[0]*fNhigh[0]+point[1]*fNhigh[1]+(point[2]-fDz)*fNhigh[2];
   if (zin>0) return kFALSE;
   // check radius
   Double_t r2 = point[0]*point[0]+point[1]*point[1];
   if ((r2<fRmin*fRmin) || (r2>fRmax*fRmax)) return kFALSE;
   // check phi
   Double_t phi = TMath::ATan2(point[1], point[0]) * TMath::RadToDeg();
   if (phi < 0 ) phi+=360.;
   Double_t dphi = fPhi2 -fPhi1;
   Double_t ddp = phi-fPhi1;
   if (ddp<0) ddp += 360.;
//   if (ddp>360) ddp-=360;
   if (ddp > dphi) return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Get range of shape for a given axis.

Double_t TGeoCtub::GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const
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
   }
   return dx;
}

////////////////////////////////////////////////////////////////////////////////
/// compute real Z coordinate of a point belonging to either lower or
/// higher caps (z should be either +fDz or -fDz)

Double_t TGeoCtub::GetZcoord(Double_t xc, Double_t yc, Double_t zc) const
{
   Double_t newz = 0;
   if (zc<0) newz =  -fDz-(xc*fNlow[0]+yc*fNlow[1])/fNlow[2];
   else      newz = fDz-(xc*fNhigh[0]+yc*fNhigh[1])/fNhigh[2];
   return newz;
}

////////////////////////////////////////////////////////////////////////////////
/// compute distance from outside point to surface of the cut tube

Double_t TGeoCtub::DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
{
   if (iact<3 && safe) {
      *safe = Safety(point, kFALSE);
      if (iact==0) return TGeoShape::Big();
      if ((iact==1) && (step<=*safe)) return TGeoShape::Big();
   }
// Check if the bounding box is crossed within the requested distance
   Double_t sdist = TGeoBBox::DistFromOutside(point,dir, fDX, fDY, fDZ, fOrigin, step);
   if (sdist>=step) return TGeoShape::Big();
   Double_t saf[2];
   saf[0] = point[0]*fNlow[0] + point[1]*fNlow[1] + (fDz+point[2])*fNlow[2];
   saf[1] = point[0]*fNhigh[0] + point[1]*fNhigh[1] + (point[2]-fDz)*fNhigh[2];
   Double_t rsq = point[0]*point[0]+point[1]*point[1];
   Double_t r = TMath::Sqrt(rsq);
   Double_t cpsi=0;
   Bool_t tub = kFALSE;
   if (TMath::Abs(fPhi2-fPhi1-360.)<1E-8) tub = kTRUE;

   // find distance to shape
   Double_t r2;
   Double_t calf = dir[0]*fNlow[0]+dir[1]*fNlow[1]+dir[2]*fNlow[2];
   // check Z planes
   Double_t xi, yi, zi;
   Double_t s = TGeoShape::Big();
   if (saf[0]>0) {
      if (calf<0) {
         s = -saf[0]/calf;
         xi = point[0]+s*dir[0];
         yi = point[1]+s*dir[1];
         r2=xi*xi+yi*yi;
         if (((fRmin*fRmin)<=r2) && (r2<=(fRmax*fRmax))) {
            if (tub) return s;
            cpsi=(xi*fCm+yi*fSm)/TMath::Sqrt(r2);
            if (cpsi>=fCdfi) return s;
         }
      }
   }
   calf = dir[0]*fNhigh[0]+dir[1]*fNhigh[1]+dir[2]*fNhigh[2];
   if (saf[1]>0) {
      if (calf<0) {
         s = -saf[1]/calf;
         xi = point[0]+s*dir[0];
         yi = point[1]+s*dir[1];
         r2=xi*xi+yi*yi;
         if (((fRmin*fRmin)<=r2) && (r2<=(fRmax*fRmax))) {
            if (tub) return s;
            cpsi=(xi*fCm+yi*fSm)/TMath::Sqrt(r2);
            if (cpsi>=fCdfi) return s;
         }
      }
   }

   // check outer cyl. surface
   Double_t nsq=dir[0]*dir[0]+dir[1]*dir[1];
   if (TMath::Abs(nsq)<1E-10) return TGeoShape::Big();
   Double_t rdotn=point[0]*dir[0]+point[1]*dir[1];
   Double_t b,d;
   // only r>fRmax coming inwards has to be considered
   if (r>fRmax && rdotn<0) {
      TGeoTube::DistToTube(rsq, nsq, rdotn, fRmax, b, d);
      if (d>0) {
         s=-b-d;
         if (s>0) {
            xi=point[0]+s*dir[0];
            yi=point[1]+s*dir[1];
            zi=point[2]+s*dir[2];
            if ((-xi*fNlow[0]-yi*fNlow[1]-(zi+fDz)*fNlow[2])>0) {
               if ((-xi*fNhigh[0]-yi*fNhigh[1]+(fDz-zi)*fNhigh[2])>0) {
                  if (tub) return s;
                  cpsi=(xi*fCm+yi*fSm)/fRmax;
                  if (cpsi>=fCdfi) return s;
               }
            }
         }
      }
   }
   // check inner cylinder
   Double_t snxt=TGeoShape::Big();
   if (fRmin>0) {
      TGeoTube::DistToTube(rsq, nsq, rdotn, fRmin, b, d);
      if (d>0) {
         s=-b+d;
         if (s>0) {
            xi=point[0]+s*dir[0];
            yi=point[1]+s*dir[1];
            zi=point[2]+s*dir[2];
            if ((-xi*fNlow[0]-yi*fNlow[1]-(zi+fDz)*fNlow[2])>0) {
               if ((-xi*fNhigh[0]-yi*fNhigh[1]+(fDz-zi)*fNhigh[2])>0) {
                  if (tub) return s;
                  cpsi=(xi*fCm+yi*fSm)/fRmin;
                  if (cpsi>=fCdfi) snxt=s;
               }
            }
         }
      }
   }
   // check phi planes
   if (tub) return snxt;
   Double_t un=dir[0]*fS1-dir[1]*fC1;
   if (un<-TGeoShape::Tolerance()) {
      s=(point[1]*fC1-point[0]*fS1)/un;
      if (s>=0) {
         xi=point[0]+s*dir[0];
         yi=point[1]+s*dir[1];
         zi=point[2]+s*dir[2];
         if ((-xi*fNlow[0]-yi*fNlow[1]-(zi+fDz)*fNlow[2])>0) {
            if ((-xi*fNhigh[0]-yi*fNhigh[1]+(fDz-zi)*fNhigh[2])>0) {
               r2=xi*xi+yi*yi;
               if ((fRmin*fRmin<=r2) && (r2<=fRmax*fRmax)) {
                  if ((yi*fCm-xi*fSm)<=0) {
                     if (s<snxt) snxt=s;
                  }
               }
            }
         }
      }
   }
   un=dir[0]*fS2-dir[1]*fC2;
   if (un>TGeoShape::Tolerance()) {
      s=(point[1]*fC2-point[0]*fS2)/un;
      if (s>=0) {
         xi=point[0]+s*dir[0];
         yi=point[1]+s*dir[1];
         zi=point[2]+s*dir[2];
         if ((-xi*fNlow[0]-yi*fNlow[1]-(zi+fDz)*fNlow[2])>0) {
            if ((-xi*fNhigh[0]-yi*fNhigh[1]+(fDz-zi)*fNhigh[2])>0) {
               r2=xi*xi+yi*yi;
               if ((fRmin*fRmin<=r2) && (r2<=fRmax*fRmax)) {
                  if ((yi*fCm-xi*fSm)>=0) {
                     if (s<snxt) snxt=s;
                  }
               }
            }
         }
      }
   }
   return snxt;
}

////////////////////////////////////////////////////////////////////////////////
/// compute distance from inside point to surface of the cut tube

Double_t TGeoCtub::DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
{
   if (iact<3 && safe) *safe = Safety(point, kTRUE);
   if (iact==0) return TGeoShape::Big();
   if ((iact==1) && (*safe>step)) return TGeoShape::Big();
   Double_t rsq = point[0]*point[0]+point[1]*point[1];
   Bool_t tub = kFALSE;
   if (TMath::Abs(fPhi2-fPhi1-360.)<1E-8) tub = kTRUE;
   // compute distance to surface
   // Do Z
   Double_t sz = TGeoShape::Big();
   Double_t saf[2];
   saf[0] = -point[0]*fNlow[0] - point[1]*fNlow[1] - (fDz+point[2])*fNlow[2];
   saf[1] = -point[0]*fNhigh[0] - point[1]*fNhigh[1] + (fDz-point[2])*fNhigh[2];
   Double_t calf = dir[0]*fNlow[0]+dir[1]*fNlow[1]+dir[2]*fNlow[2];
   if (calf>0) sz = saf[0]/calf;

   Double_t sz1=TGeoShape::Big();
   calf = dir[0]*fNhigh[0]+dir[1]*fNhigh[1]+dir[2]*fNhigh[2];
   if (calf>0) {
      sz1 = saf[1]/calf;
      if (sz1<sz) sz = sz1;
   }

   // Do R
   Double_t nsq=dir[0]*dir[0]+dir[1]*dir[1];
   // track parallel to Z
   if (TMath::Abs(nsq)<1E-10) return sz;
   Double_t rdotn=point[0]*dir[0]+point[1]*dir[1];
   Double_t sr=TGeoShape::Big();
   Double_t b, d;
   Bool_t skip_outer = kFALSE;
   // inner cylinder
   if (fRmin>1E-10) {
      TGeoTube::DistToTube(rsq, nsq, rdotn, fRmin, b, d);
      if (d>0) {
         sr=-b-d;
         if (sr>0) skip_outer = kTRUE;
      }
   }
   // outer cylinder
   if (!skip_outer) {
      TGeoTube::DistToTube(rsq, nsq, rdotn, fRmax, b, d);
      if (d>0) {
         sr=-b+d;
         if (sr<0) sr=TGeoShape::Big();
      } else {
         return 0.; // already outside
      }
   }
   // phi planes
   Double_t sfmin = TGeoShape::Big();
   if (!tub) sfmin=TGeoShape::DistToPhiMin(point, dir, fS1, fC1, fS2, fC2, fSm, fCm);
   return TMath::Min(TMath::Min(sz,sr), sfmin);
}

////////////////////////////////////////////////////////////////////////////////
/// Divide the tube along one axis.

TGeoVolume *TGeoCtub::Divide(TGeoVolume * /*voldiv*/, const char * /*divname*/, Int_t /*iaxis*/, Int_t /*ndiv*/,
                             Double_t /*start*/, Double_t /*step*/)
{
   Warning("Divide", "In shape %s division of a cut tube not implemented", GetName());
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// in case shape has some negative parameters, these has to be computed
/// in order to fit the mother

TGeoShape *TGeoCtub::GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix * /*mat*/) const
{
   if (!TestShapeBit(kGeoRunTimeShape)) return 0;
   if (!mother->TestShapeBit(kGeoTube)) {
      Error("GetMakeRuntimeShape", "Invalid mother for shape %s", GetName());
      return 0;
   }
   Double_t rmin, rmax, dz;
   rmin = fRmin;
   rmax = fRmax;
   dz = fDz;
   if (fDz<0) dz=((TGeoTube*)mother)->GetDz();
   if (fRmin<0)
      rmin = ((TGeoTube*)mother)->GetRmin();
   if ((fRmax<0) || (fRmax<=fRmin))
      rmax = ((TGeoTube*)mother)->GetRmax();

   return (new TGeoCtub(rmin, rmax, dz, fPhi1, fPhi2, fNlow[0], fNlow[1], fNlow[2],
                        fNhigh[0], fNhigh[1], fNhigh[2]));
}

////////////////////////////////////////////////////////////////////////////////
/// print shape parameters

void TGeoCtub::InspectShape() const
{
   printf("*** Shape %s: TGeoCtub ***\n", GetName());
   printf("    lx = %11.5f\n", fNlow[0]);
   printf("    ly = %11.5f\n", fNlow[1]);
   printf("    lz = %11.5f\n", fNlow[2]);
   printf("    tx = %11.5f\n", fNhigh[0]);
   printf("    ty = %11.5f\n", fNhigh[1]);
   printf("    tz = %11.5f\n", fNhigh[2]);
   TGeoTubeSeg::InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// computes the closest distance from given point to this shape, according
/// to option. The matching point on the shape is stored in spoint.

Double_t TGeoCtub::Safety(const Double_t *point, Bool_t in) const
{
   Double_t saf[4];
   Double_t rsq = point[0]*point[0]+point[1]*point[1];
   Double_t r = TMath::Sqrt(rsq);
   Bool_t isseg = kTRUE;
   if (TMath::Abs(fPhi2-fPhi1-360.)<1E-8) isseg=kFALSE;

   saf[0] = -point[0]*fNlow[0] - point[1]*fNlow[1] - (fDz+point[2])*fNlow[2];
   saf[1] = -point[0]*fNhigh[0] - point[1]*fNhigh[1] + (fDz-point[2])*fNhigh[2];
   saf[2] = (fRmin<1E-10 && !isseg)?TGeoShape::Big():(r-fRmin);
   saf[3] = fRmax-r;
   Double_t safphi = TGeoShape::Big();
   Double_t safe = TGeoShape::Big();
   if (isseg) safphi =  TGeoShape::SafetyPhi(point, in, fPhi1, fPhi2);

   if (in) {
      safe = saf[TMath::LocMin(4,saf)];
      return TMath::Min(safe, safphi);
   }
   for (Int_t i=0; i<4; i++) saf[i]=-saf[i];
   safe = saf[TMath::LocMax(4,saf)];
   if (isseg) return TMath::Max(safe, safphi);
   return safe;
}

////////////////////////////////////////////////////////////////////////////////
/// set dimensions of a cut tube

void TGeoCtub::SetCtubDimensions(Double_t rmin, Double_t rmax, Double_t dz, Double_t phi1, Double_t phi2,
                   Double_t lx, Double_t ly, Double_t lz, Double_t tx, Double_t ty, Double_t tz)
{
   SetTubsDimensions(rmin, rmax, dz, phi1, phi2);
   fNlow[0] = lx;
   fNlow[1] = ly;
   fNlow[2] = lz;
   fNhigh[0] = tx;
   fNhigh[1] = ty;
   fNhigh[2] = tz;
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoCtub::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive)) return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   rmin = " << fRmin << ";" << std::endl;
   out << "   rmax = " << fRmax << ";" << std::endl;
   out << "   dz   = " << fDz << ";" << std::endl;
   out << "   phi1 = " << fPhi1 << ";" << std::endl;
   out << "   phi2 = " << fPhi2 << ";" << std::endl;
   out << "   lx   = " << fNlow[0] << ";" << std::endl;
   out << "   ly   = " << fNlow[1] << ";" << std::endl;
   out << "   lz   = " << fNlow[2] << ";" << std::endl;
   out << "   tx   = " << fNhigh[0] << ";" << std::endl;
   out << "   ty   = " << fNhigh[1] << ";" << std::endl;
   out << "   tz   = " << fNhigh[2] << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoCtub(\"" << GetName() << "\",rmin,rmax,dz,phi1,phi2,lx,ly,lz,tx,ty,tz);" << std::endl;   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Set dimensions of the cut tube starting from a list.

void TGeoCtub::SetDimensions(Double_t *param)
{
   SetCtubDimensions(param[0], param[1], param[2], param[3], param[4], param[5],
                     param[6], param[7], param[8], param[9], param[10]);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Fills array with n random points located on the line segments of the shape mesh.
/// The output array must be provided with a length of minimum 3*npoints. Returns
/// true if operation is implemented.

Bool_t TGeoCtub::GetPointsOnSegments(Int_t /*npoints*/, Double_t * /*array*/) const
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create mesh points for the cut tube.

void TGeoCtub::SetPoints(Double_t *points) const
{
   Double_t dz;
   Int_t j, n;
   Double_t phi, phi1, phi2, dphi;
   phi1 = fPhi1;
   phi2 = fPhi2;
   if (phi2<phi1) phi2+=360.;
   n = gGeoManager->GetNsegments()+1;

   dphi = (phi2-phi1)/(n-1);
   dz   = fDz;

   if (points) {
      Int_t indx = 0;

      for (j = 0; j < n; j++) {
         phi = (phi1+j*dphi)*TMath::DegToRad();
         points[indx+6*n] = points[indx] = fRmin * TMath::Cos(phi);
         indx++;
         points[indx+6*n] = points[indx] = fRmin * TMath::Sin(phi);
         indx++;
         points[indx+6*n] = GetZcoord(points[indx-2], points[indx-1], dz);
         points[indx]     = GetZcoord(points[indx-2], points[indx-1], -dz);
         indx++;
      }
      for (j = 0; j < n; j++) {
         phi = (phi1+j*dphi)*TMath::DegToRad();
         points[indx+6*n] = points[indx] = fRmax * TMath::Cos(phi);
         indx++;
         points[indx+6*n] = points[indx] = fRmax * TMath::Sin(phi);
         indx++;
         points[indx+6*n]= GetZcoord(points[indx-2], points[indx-1], dz);
         points[indx]    = GetZcoord(points[indx-2], points[indx-1], -dz);
         indx++;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create mesh points for the cut tube.

void TGeoCtub::SetPoints(Float_t *points) const
{
   Double_t dz;
   Int_t j, n;
   Double_t phi, phi1, phi2, dphi;
   phi1 = fPhi1;
   phi2 = fPhi2;
   if (phi2<phi1) phi2+=360.;
   n = gGeoManager->GetNsegments()+1;

   dphi = (phi2-phi1)/(n-1);
   dz   = fDz;

   if (points) {
      Int_t indx = 0;

      for (j = 0; j < n; j++) {
         phi = (phi1+j*dphi)*TMath::DegToRad();
         points[indx+6*n] = points[indx] = fRmin * TMath::Cos(phi);
         indx++;
         points[indx+6*n] = points[indx] = fRmin * TMath::Sin(phi);
         indx++;
         points[indx+6*n] = GetZcoord(points[indx-2], points[indx-1], dz);
         points[indx]     = GetZcoord(points[indx-2], points[indx-1], -dz);
         indx++;
      }
      for (j = 0; j < n; j++) {
         phi = (phi1+j*dphi)*TMath::DegToRad();
         points[indx+6*n] = points[indx] = fRmax * TMath::Cos(phi);
         indx++;
         points[indx+6*n] = points[indx] = fRmax * TMath::Sin(phi);
         indx++;
         points[indx+6*n]= GetZcoord(points[indx-2], points[indx-1], dz);
         points[indx]    = GetZcoord(points[indx-2], points[indx-1], -dz);
         indx++;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoCtub::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   TGeoTubeSeg::GetMeshNumbers(nvert,nsegs,npols);
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of vertices of the mesh representation

Int_t TGeoCtub::GetNmeshVertices() const
{
   Int_t n = gGeoManager->GetNsegments()+1;
   Int_t numPoints = n*4;
   return numPoints;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills a static 3D buffer and returns a reference.

const TBuffer3D & TGeoCtub::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
{
   static TBuffer3DCutTube buffer;

   TGeoBBox::FillBuffer3D(buffer, reqSections, localFrame);

   if (reqSections & TBuffer3D::kShapeSpecific) {
      // These from TBuffer3DCutTube / TGeoCtub
      buffer.fRadiusInner  = fRmin;
      buffer.fRadiusOuter  = fRmax;
      buffer.fHalfLength   = fDz;
      buffer.fPhiMin       = fPhi1;
      buffer.fPhiMax       = fPhi2;

      for (UInt_t i = 0; i < 3; i++ ) {
         buffer.fLowPlaneNorm[i] = fNlow[i];
         buffer.fHighPlaneNorm[i] = fNhigh[i];
      }
      buffer.SetSectionsValid(TBuffer3D::kShapeSpecific);
   }
   if (reqSections & TBuffer3D::kRawSizes) {
      Int_t n = gGeoManager->GetNsegments()+1;
      Int_t nbPnts = 4*n;
      Int_t nbSegs = 2*nbPnts;
      Int_t nbPols = nbPnts-2;
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

void TGeoCtub::Contains_v(const Double_t *points, Bool_t *inside, Int_t vecsize) const
{
   for (Int_t i=0; i<vecsize; i++) inside[i] = Contains(&points[3*i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the normal for an array o points so that norm.dot.dir is positive
/// Input: Arrays of point coordinates and directions + vector size
/// Output: Array of normal directions

void TGeoCtub::ComputeNormal_v(const Double_t *points, const Double_t *dirs, Double_t *norms, Int_t vecsize)
{
   for (Int_t i=0; i<vecsize; i++) ComputeNormal(&points[3*i], &dirs[3*i], &norms[3*i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoCtub::DistFromInside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromInside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoCtub::DistFromOutside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromOutside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute safe distance from each of the points in the input array.
/// Input: Array of point coordinates, array of statuses for these points, size of the arrays
/// Output: Safety values

void TGeoCtub::Safety_v(const Double_t *points, const Bool_t *inside, Double_t *safe, Int_t vecsize) const
{
   for (Int_t i=0; i<vecsize; i++) safe[i] = Safety(&points[3*i], inside[i]);
}
