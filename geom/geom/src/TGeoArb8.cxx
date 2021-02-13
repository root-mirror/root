// @(#)root/geom:$Id$
// Author: Andrei Gheata   31/01/02

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGeoArb8.h"

#include <iostream>
#include "TBuffer.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TMath.h"

ClassImp(TGeoArb8);

/** \class TGeoArb8
\ingroup Geometry_classes

An arbitrary trapezoid with less than 8 vertices standing on
two parallel planes perpendicular to Z axis. Parameters :
  - dz - half length in Z;
  - xy[8][2] - vector of (x,y) coordinates of vertices
     - first four points (xy[i][j], i<4, j<2) are the (x,y)
       coordinates of the vertices sitting on the -dz plane;
     - last four points (xy[i][j], i>=4, j<2) are the (x,y)
       coordinates of the vertices sitting on the +dz plane;

The order of defining the vertices of an arb8 is the following :
  - point 0 is connected with points 1,3,4
  - point 1 is connected with points 0,2,5
  - point 2 is connected with points 1,3,6
  - point 3 is connected with points 0,2,7
  - point 4 is connected with points 0,5,7
  - point 5 is connected with points 1,4,6
  - point 6 is connected with points 2,5,7
  - point 7 is connected with points 3,4,6

Points can be identical in order to create shapes with less than
8 vertices.

Begin_Macro(source)
{
   TCanvas *c = new TCanvas("c", "c",0,0,600,600);
   new TGeoManager("arb8", "poza12");
   TGeoMaterial *mat = new TGeoMaterial("Al", 26.98,13,2.7);
   TGeoMedium *med = new TGeoMedium("MED",1,mat);
   TGeoVolume *top = gGeoManager->MakeBox("TOP",med,100,100,100);
   gGeoManager->SetTopVolume(top);
   TGeoArb8 *arb = new TGeoArb8(20);
   arb->SetVertex(0,-30,-25);
   arb->SetVertex(1,-25,25);
   arb->SetVertex(2,5,25);
   arb->SetVertex(3,25,-25);
   arb->SetVertex(4,-28,-23);
   arb->SetVertex(5,-23,27);
   arb->SetVertex(6,-23,27);
   arb->SetVertex(7,13,-27);
   TGeoVolume *vol = new TGeoVolume("ARB8",arb,med);
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

/** \class TGeoGtra
\ingroup Geometry_classes

Gtra is a twisted trapezoid.
i.e. one for which the faces perpendicular
to z are trapezia and their centres are not the same x, y. It has 12
parameters: the half length in z, the polar angles from the centre of
the face at low z to that at high z, twist, H1 the half length in y at low z,
LB1 the half length in x at low z and y low edge, LB2 the half length
in x at low z and y high edge, TH1 the angle w.r.t. the y axis from the
centre of low y edge to the centre of the high y edge, and H2, LB2,
LH2, TH2, the corresponding quantities at high z.

Begin_Macro(source)
{
   TCanvas *c = new TCanvas("c", "c",0,0,600,600);
   new TGeoManager("gtra", "poza11");
   TGeoMaterial *mat = new TGeoMaterial("Al", 26.98,13,2.7);
   TGeoMedium *med = new TGeoMedium("MED",1,mat);
   TGeoVolume *top = gGeoManager->MakeBox("TOP",med,100,100,100);
   gGeoManager->SetTopVolume(top);
   TGeoVolume *vol = gGeoManager->MakeGtra("Gtra",med, 30,15,30,30,20,10,15,0,20,10,15,0);
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

/** \class TGeoTrap
\ingroup Geometry_classes

TRAP is a general trapezoid, i.e. one for which the faces perpendicular
to z are trapezia and their centres are not the same x, y. It has 11
parameters: the half length in z, the polar angles from the centre of
the face at low z to that at high z, H1 the half length in y at low z,
LB1 the half length in x at low z and y low edge, LB2 the half length
in x at low z and y high edge, TH1 the angle w.r.t. the y axis from the
centre of low y edge to the centre of the high y edge, and H2, LB2,
LH2, TH2, the corresponding quantities at high z.

Begin_Macro(source)
{
   TCanvas *c = new TCanvas("c", "c",0,0,600,600);
   new TGeoManager("trap", "poza10");
   TGeoMaterial *mat = new TGeoMaterial("Al", 26.98,13,2.7);
   TGeoMedium *med = new TGeoMedium("MED",1,mat);
   TGeoVolume *top = gGeoManager->MakeBox("TOP",med,100,100,100);
   gGeoManager->SetTopVolume(top);
   TGeoVolume *vol = gGeoManager->MakeTrap("Trap",med, 30,15,30,20,10,15,0,20,10,15,0);
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

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TGeoArb8::TGeoArb8()
{
   fDz = 0;
   for (Int_t i=0; i<8; i++) {
      fXY[i][0] = 0.0;
      fXY[i][1] = 0.0;
   }
   SetShapeBit(kGeoArb8);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor. If the array of vertices is not null, this should be
/// in the format : (x0, y0, x1, y1, ... , x7, y7)

TGeoArb8::TGeoArb8(Double_t dz, Double_t *vertices)
         :TGeoBBox(0,0,0)
{
   fDz = dz;
   SetShapeBit(kGeoArb8);
   if (vertices) {
      for (Int_t i=0; i<8; i++) {
         fXY[i][0] = vertices[2*i];
         fXY[i][1] = vertices[2*i+1];
      }
      ComputeTwist();
      ComputeBBox();
   } else {
      for (Int_t i=0; i<8; i++) {
         fXY[i][0] = 0.0;
         fXY[i][1] = 0.0;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Named constructor. If the array of vertices is not null, this should be
/// in the format : (x0, y0, x1, y1, ... , x7, y7)

TGeoArb8::TGeoArb8(const char *name, Double_t dz, Double_t *vertices)
         :TGeoBBox(name, 0,0,0)
{
   fDz = dz;
   SetShapeBit(kGeoArb8);
   if (vertices) {
      for (Int_t i=0; i<8; i++) {
         fXY[i][0] = vertices[2*i];
         fXY[i][1] = vertices[2*i+1];
      }
      ComputeTwist();
      ComputeBBox();
   } else {
      for (Int_t i=0; i<8; i++) {
         fXY[i][0] = 0.0;
         fXY[i][1] = 0.0;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TGeoArb8::TGeoArb8(const TGeoArb8& ga8) :
  TGeoBBox(ga8),
  fDz(ga8.fDz)
{
   for(Int_t i=0; i<8; i++) {
      fXY[i][0]=ga8.fXY[i][0];
      fXY[i][1]=ga8.fXY[i][1];
   }
   CopyTwist(ga8.fTwist);
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TGeoArb8& TGeoArb8::operator=(const TGeoArb8& ga8)
{
   if(this!=&ga8) {
      TGeoBBox::operator=(ga8);
      fDz=ga8.fDz;
      CopyTwist(ga8.fTwist);
      for(Int_t i=0; i<8; i++) {
         fXY[i][0]=ga8.fXY[i][0];
         fXY[i][1]=ga8.fXY[i][1];
      }
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGeoArb8::~TGeoArb8()
{
   if (fTwist) delete [] fTwist;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy twist values from source array

void TGeoArb8::CopyTwist(Double_t *twist)
{
   if (twist) {
      if (!fTwist) fTwist = new Double_t[4];
      memcpy(fTwist, twist, 4*sizeof(Double_t));
   } else if (fTwist) {
      delete [] fTwist;
      fTwist = nullptr;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Computes capacity of the shape in [length^3].

Double_t TGeoArb8::Capacity() const
{
   Int_t i,j;
   Double_t capacity = 0;
   for (i=0; i<4; i++) {
      j = (i+1)%4;
      capacity += 0.25*fDz*((fXY[i][0]+fXY[i+4][0])*(fXY[j][1]+fXY[j+4][1]) -
                            (fXY[j][0]+fXY[j+4][0])*(fXY[i][1]+fXY[i+4][1]) +
                    (1./3)*((fXY[i+4][0]-fXY[i][0])*(fXY[j+4][1]-fXY[j][1]) -
                            (fXY[j][0]-fXY[j+4][0])*(fXY[i][1]-fXY[i+4][1])));
   }
   return TMath::Abs(capacity);
}

////////////////////////////////////////////////////////////////////////////////
/// Computes bounding box for an Arb8 shape.

void TGeoArb8::ComputeBBox()
{
   Double_t xmin, xmax, ymin, ymax;
   xmin = xmax = fXY[0][0];
   ymin = ymax = fXY[0][1];

   for (Int_t i=1; i<8; i++) {
      if (xmin>fXY[i][0]) xmin=fXY[i][0];
      if (xmax<fXY[i][0]) xmax=fXY[i][0];
      if (ymin>fXY[i][1]) ymin=fXY[i][1];
      if (ymax<fXY[i][1]) ymax=fXY[i][1];
   }
   fDX = 0.5*(xmax-xmin);
   fDY = 0.5*(ymax-ymin);
   fDZ = fDz;
   fOrigin[0] = 0.5*(xmax+xmin);
   fOrigin[1] = 0.5*(ymax+ymin);
   fOrigin[2] = 0;
   SetShapeBit(kGeoClosedShape);
}

////////////////////////////////////////////////////////////////////////////////
/// Computes tangents of twist angles (angles between projections on XY plane
/// of corresponding -dz +dz edges). Computes also if the vertices are defined
/// clockwise or anti-clockwise.

void TGeoArb8::ComputeTwist()
{
   Double_t twist[4];
   Bool_t twisted = kFALSE;
   Double_t dx1, dy1, dx2, dy2;
   Bool_t singleBottom = kTRUE;
   Bool_t singleTop = kTRUE;
   Int_t i;
   for (i=0; i<4; i++) {
      dx1 = fXY[(i+1)%4][0]-fXY[i][0];
      dy1 = fXY[(i+1)%4][1]-fXY[i][1];
      if (TMath::Abs(dx1)<TGeoShape::Tolerance() && TMath::Abs(dy1)<TGeoShape::Tolerance()) {
         twist[i] = 0;
         continue;
      }
      singleBottom = kFALSE;
      dx2 = fXY[4+(i+1)%4][0]-fXY[4+i][0];
      dy2 = fXY[4+(i+1)%4][1]-fXY[4+i][1];
      if (TMath::Abs(dx2)<TGeoShape::Tolerance() && TMath::Abs(dy2)<TGeoShape::Tolerance()) {
         twist[i] = 0;
         continue;
      }
      singleTop = kFALSE;
      twist[i] = dy1*dx2 - dx1*dy2;
      if (TMath::Abs(twist[i])<TGeoShape::Tolerance()) {
         twist[i] = 0;
         continue;
      }
      twist[i] = TMath::Sign(1.,twist[i]);
      twisted = kTRUE;
   }

   CopyTwist(twisted ? twist : nullptr);

   if (singleBottom) {
      for (i=0; i<4; i++) {
         fXY[i][0] += 1.E-8*fXY[i+4][0];
         fXY[i][1] += 1.E-8*fXY[i+4][1];
      }
   }
   if (singleTop) {
      for (i=0; i<4; i++) {
         fXY[i+4][0] += 1.E-8*fXY[i][0];
         fXY[i+4][1] += 1.E-8*fXY[i][1];
      }
   }
   Double_t sum1 = 0.;
   Double_t sum2 = 0.;
   Int_t j;
   for (i=0; i<4; i++) {
      j = (i+1)%4;
      sum1 += fXY[i][0]*fXY[j][1]-fXY[j][0]*fXY[i][1];
      sum2 += fXY[i+4][0]*fXY[j+4][1]-fXY[j+4][0]*fXY[i+4][1];
   }
   if (sum1*sum2 < -TGeoShape::Tolerance()) {
      Fatal("ComputeTwist", "Shape %s type Arb8: Lower/upper faces defined with opposite clockwise", GetName());
      return;
   }
   if (sum1>TGeoShape::Tolerance()) {
      Error("ComputeTwist", "Shape %s type Arb8: Vertices must be defined clockwise in XY planes. Re-ordering...", GetName());
      Double_t xtemp, ytemp;
      xtemp = fXY[1][0];
      ytemp = fXY[1][1];
      fXY[1][0] = fXY[3][0];
      fXY[1][1] = fXY[3][1];
      fXY[3][0] = xtemp;
      fXY[3][1] = ytemp;
      xtemp = fXY[5][0];
      ytemp = fXY[5][1];
      fXY[5][0] = fXY[7][0];
      fXY[5][1] = fXY[7][1];
      fXY[7][0] = xtemp;
      fXY[7][1] = ytemp;
   }
   // Check for illegal crossings.
   Bool_t illegal_cross = kFALSE;
   illegal_cross = TGeoShape::IsSegCrossing(fXY[0][0],fXY[0][1],fXY[1][0],fXY[1][1],
                                            fXY[2][0],fXY[2][1],fXY[3][0],fXY[3][1]);
   if (!illegal_cross)
      illegal_cross = TGeoShape::IsSegCrossing(fXY[4][0],fXY[4][1],fXY[5][0],fXY[5][1],
                                               fXY[6][0],fXY[6][1],fXY[7][0],fXY[7][1]);
   if (illegal_cross) {
      Error("ComputeTwist", "Shape %s type Arb8: Malformed polygon with crossing opposite segments", GetName());
      InspectShape();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get twist for segment I in range [0,3]

Double_t TGeoArb8::GetTwist(Int_t iseg) const
{
   return (!fTwist || iseg<0 || iseg>3) ? 0. : fTwist[iseg];
}

////////////////////////////////////////////////////////////////////////////////
/// Get index of the edge of the quadrilater represented by vert closest to point.
/// If [P1,P2] is the closest segment and P is the point, the function returns the fraction of the
/// projection of (P1P) over (P1P2). If projection of P is not in range [P1,P2] return -1.

Double_t TGeoArb8::GetClosestEdge(const Double_t *point, Double_t *vert, Int_t &isegment) const
{
   isegment = 0;
   Int_t isegmin = 0;
   Int_t i1, i2;
   Double_t p1[2], p2[2];
   Double_t lsq, ssq, dx, dy, dpx, dpy, u;
   Double_t umin = -1.;
   Double_t safe=1E30;
   for (i1=0; i1<4; i1++) {
      if (TGeoShape::IsSameWithinTolerance(safe,0)) {
         isegment = isegmin;
         return umin;
      }
      i2 = (i1+1)%4;
      p1[0] = vert[2*i1];
      p1[1] = vert[2*i1+1];
      p2[0] = vert[2*i2];
      p2[1] = vert[2*i2+1];
      dx = p2[0] - p1[0];
      dy = p2[1] - p1[1];
      dpx = point[0] - p1[0];
      dpy = point[1] - p1[1];
      lsq = dx*dx + dy*dy;
      // Current segment collapsed to a point
      if (TGeoShape::IsSameWithinTolerance(lsq,0)) {
         ssq = dpx*dpx + dpy*dpy;
         if (ssq < safe) {
            safe = ssq;
            isegmin = i1;
            umin = -1;
         }
         continue;
      }
      // Projection fraction
      u = (dpx*dx + dpy*dy)/lsq;
      // If projection of P is outside P1P2 limits, take the distance from P to the closest of P1 and P2
      if (u>1) {
         // Outside, closer to P2
         dpx = point[0]-p2[0];
         dpy = point[1]-p2[1];
         u = -1.;
      } else {
         if (u>=0) {
            // Projection inside
            dpx -= u*dx;
            dpy -= u*dy;
         } else {
            // Outside, closer to P1
            u = -1.;
         }
      }
      ssq = dpx*dpx + dpy*dpy;
      if (ssq < safe) {
         safe = ssq;
         isegmin = i1;
         umin = u;
      }
   }
   isegment = isegmin;
   // safe = TMath::Sqrt(safe);
   return umin;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute normal to closest surface from POINT.

void TGeoArb8::ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm)
{
   Double_t safc;
   Double_t x0, y0, z0, x1, y1, z1, x2, y2, z2;
   Double_t ax, ay, az, bx, by, bz;
   Double_t fn;
   safc = fDz-TMath::Abs(point[2]);
   if (safc<10.*TGeoShape::Tolerance()) {
      memset(norm,0,3*sizeof(Double_t));
      norm[2] = (dir[2]>0)?1:(-1);
      return;
   }
   Double_t vert[8];
   SetPlaneVertices(point[2], vert);
   // Get the closest edge (point should be on this edge within tolerance)
   Int_t iseg;
   Double_t frac = GetClosestEdge(point, vert, iseg);
   if (frac<0) frac = 0.;
   Int_t jseg = (iseg+1)%4;
   x0 = vert[2*iseg];
   y0 = vert[2*iseg+1];
   z0 = point[2];
   x2 = vert[2*jseg];
   y2 = vert[2*jseg+1];
   z2 = point[2];
   x0 += frac*(x2-x0);
   y0 += frac*(y2-y0);
   x1 = fXY[iseg+4][0];
   y1 = fXY[iseg+4][1];
   z1 = fDz;
   x1 += frac*(fXY[jseg+4][0]-x1);
   y1 += frac*(fXY[jseg+4][1]-y1);
   ax = x1-x0;
   ay = y1-y0;
   az = z1-z0;
   bx = x2-x0;
   by = y2-y0;
   bz = z2-z0;
   // Cross product of the vector given by the section segment (that contains the point) at z=point[2]
   // and the vector connecting the point projection to its correspondent on the top edge (hard to swallow, isn't it?)
   norm[0] = ay*bz-az*by;
   norm[1] = az*bx-ax*bz;
   norm[2] = ax*by-ay*bx;
   fn = TMath::Sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
   // If point is on one edge, fn may be very small and the normal does not make sense - avoid divzero
   if (fn<1E-10) {
      norm[0] = 1.;
      norm[1] = 0.;
      norm[2] = 0.;
   } else {
      norm[0] /= fn;
      norm[1] /= fn;
      norm[2] /= fn;
   }
   // Make sure the dot product of the normal and the direction is positive.
   if (dir[0]>-2. && dir[0]*norm[0]+dir[1]*norm[1]+dir[2]*norm[2] < 0) {
      norm[0] = -norm[0];
      norm[1] = -norm[1];
      norm[2] = -norm[2];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Test if point is inside this shape.

Bool_t TGeoArb8::Contains(const Double_t *point) const
{
   // first check Z range
   if (TMath::Abs(point[2]) > fDz) return kFALSE;
   // compute intersection between Z plane containing point and the arb8
   Double_t poly[8];
   // memset(&poly[0], 0, 8*sizeof(Double_t));
   // SetPlaneVertices(point[2], &poly[0]);
   Double_t cf = 0.5*(fDz-point[2])/fDz;
   Int_t i;
   for (i=0; i<4; i++) {
      poly[2*i]   = fXY[i+4][0]+cf*(fXY[i][0]-fXY[i+4][0]);
      poly[2*i+1] = fXY[i+4][1]+cf*(fXY[i][1]-fXY[i+4][1]);
   }
   return InsidePolygon(point[0],point[1],poly);
}

////////////////////////////////////////////////////////////////////////////////
/// Computes distance to plane ipl :
///  - ipl=0 : points 0,4,1,5
///  - ipl=1 : points 1,5,2,6
///  - ipl=2 : points 2,6,3,7
///  - ipl=3 : points 3,7,0,4

Double_t TGeoArb8::DistToPlane(const Double_t *point, const Double_t *dir, Int_t ipl, Bool_t in) const
{
   Double_t xa,xb,xc,xd;
   Double_t ya,yb,yc,yd;
   Double_t eps = 10.*TGeoShape::Tolerance();
   Double_t norm[3];
   Double_t dirp[3];
   Double_t ndotd = 0;
   Int_t j = (ipl+1)%4;
   xa=fXY[ipl][0];
   ya=fXY[ipl][1];
   xb=fXY[ipl+4][0];
   yb=fXY[ipl+4][1];
   xc=fXY[j][0];
   yc=fXY[j][1];
   xd=fXY[4+j][0];
   yd=fXY[4+j][1];
   Double_t dz2 =0.5/fDz;
   Double_t tx1 =dz2*(xb-xa);
   Double_t ty1 =dz2*(yb-ya);
   Double_t tx2 =dz2*(xd-xc);
   Double_t ty2 =dz2*(yd-yc);
   Double_t dzp =fDz+point[2];
   Double_t xs1 =xa+tx1*dzp;
   Double_t ys1 =ya+ty1*dzp;
   Double_t xs2 =xc+tx2*dzp;
   Double_t ys2 =yc+ty2*dzp;
   Double_t dxs =xs2-xs1;
   Double_t dys =ys2-ys1;
   Double_t dtx =tx2-tx1;
   Double_t dty =ty2-ty1;
   Double_t a=(dtx*dir[1]-dty*dir[0]+(tx1*ty2-tx2*ty1)*dir[2])*dir[2];
   Double_t signa = TMath::Sign(1.,a);
   Double_t b=dxs*dir[1]-dys*dir[0]+(dtx*point[1]-dty*point[0]+ty2*xs1-ty1*xs2
              +tx1*ys2-tx2*ys1)*dir[2];
   Double_t c=dxs*point[1]-dys*point[0]+xs1*ys2-xs2*ys1;
   Double_t s=TGeoShape::Big();
   Double_t x1,x2,y1,y2,xp,yp,zi;
   if (TMath::Abs(a)<eps) {
      // Surface is planar
      if (TMath::Abs(b)<eps) return TGeoShape::Big(); // Track parallel to surface
      s=-c/b;
      if (TMath::Abs(s)<1.E-6 && TMath::Abs(TMath::Abs(point[2])-fDz)>eps) {
         memcpy(dirp,dir,3*sizeof(Double_t));
         dirp[0] = -3;
         // Compute normal pointing outside
         ((TGeoArb8*)this)->ComputeNormal(point,dirp,norm);
         ndotd = dir[0]*norm[0]+dir[1]*norm[1]+dir[2]*norm[2];
         if (!in) ndotd*=-1.;
         if (ndotd>0) {
            s = TMath::Max(0.,s);
            zi = (point[0]-xs1)*(point[0]-xs2)+(point[1]-ys1)*(point[1]-ys2);
            if (zi<=0) return s;
         }
         return TGeoShape::Big();
      }
      if (s<0) return TGeoShape::Big();
   } else {
      Double_t d=b*b-4*a*c;
      if (d<0) return TGeoShape::Big();
      Double_t smin=0.5*(-b-signa*TMath::Sqrt(d))/a;
      Double_t smax=0.5*(-b+signa*TMath::Sqrt(d))/a;
      s = smin;
      if (TMath::Abs(s)<1.E-6 && TMath::Abs(TMath::Abs(point[2])-fDz)>eps) {
         memcpy(dirp,dir,3*sizeof(Double_t));
         dirp[0] = -3;
         // Compute normal pointing outside
         ((TGeoArb8*)this)->ComputeNormal(point,dirp,norm);
         ndotd = dir[0]*norm[0]+dir[1]*norm[1]+dir[2]*norm[2];
         if (!in) ndotd*=-1.;
         if (ndotd>0) return TMath::Max(0.,s);
         s = 0.; // ignore
      }
      if (s>eps) {
         // Check smin
         zi=point[2]+s*dir[2];
         if (TMath::Abs(zi)<fDz) {
            x1=xs1+tx1*dir[2]*s;
            x2=xs2+tx2*dir[2]*s;
            xp=point[0]+s*dir[0];
            y1=ys1+ty1*dir[2]*s;
            y2=ys2+ty2*dir[2]*s;
            yp=point[1]+s*dir[1];
            zi = (xp-x1)*(xp-x2)+(yp-y1)*(yp-y2);
            if (zi<=0) return s;
         }
      }
      // Smin failed, try smax
      s=smax;
      if (TMath::Abs(s)<1.E-6 && TMath::Abs(TMath::Abs(point[2])-fDz)>eps) {
         memcpy(dirp,dir,3*sizeof(Double_t));
         dirp[0] = -3;
         // Compute normal pointing outside
         ((TGeoArb8*)this)->ComputeNormal(point,dirp,norm);
         ndotd = dir[0]*norm[0]+dir[1]*norm[1]+dir[2]*norm[2];
         if (!in) ndotd*=-1.;
         if (ndotd>0) s = TMath::Max(0.,s);
         else         s = TGeoShape::Big();
         return s;
      }
   }
   if (s>eps) {
      // Check smin
      zi=point[2]+s*dir[2];
      if (TMath::Abs(zi)<fDz) {
         x1=xs1+tx1*dir[2]*s;
         x2=xs2+tx2*dir[2]*s;
         xp=point[0]+s*dir[0];
         y1=ys1+ty1*dir[2]*s;
         y2=ys2+ty2*dir[2]*s;
         yp=point[1]+s*dir[1];
         zi = (xp-x1)*(xp-x2)+(yp-y1)*(yp-y2);
         if (zi<=0) return s;
      }
   }
   return TGeoShape::Big();
}

////////////////////////////////////////////////////////////////////////////////
/// Computes distance from outside point to surface of the shape.

Double_t TGeoArb8::DistFromOutside(const Double_t *point, const Double_t *dir, Int_t /*iact*/, Double_t step, Double_t * /*safe*/) const
{
   Double_t sdist = TGeoBBox::DistFromOutside(point,dir, fDX, fDY, fDZ, fOrigin, step);
   if (sdist>=step) return TGeoShape::Big();
   Double_t snext;
   // check Z planes
   if (TMath::Abs(point[2])>fDz-1.E-8) {
      Double_t pt[3];
      if (point[2]*dir[2]<0) {
         pt[2]=fDz*TMath::Sign(1.,point[2]);
         snext = TMath::Max((pt[2]-point[2])/dir[2],0.);
         for (Int_t j=0; j<2; j++) pt[j]=point[j]+snext*dir[j];
         if (Contains(&pt[0])) return snext;
      }
   }
   // check lateral faces
   Double_t dist;
   snext = TGeoShape::Big();
   for (Int_t i=0; i<4; i++) {
      dist = DistToPlane(point, dir, i, kFALSE);
      if (dist<snext) snext = dist;
   }
   return snext;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from inside point to surface of the shape.

Double_t TGeoArb8::DistFromInside(const Double_t *point, const Double_t *dir, Int_t /*iact*/, Double_t /*step*/, Double_t * /*safe*/) const
{
   Int_t i;
   Double_t distz = TGeoShape::Big();
   Double_t distl = TGeoShape::Big();
   Double_t dist;
   Double_t pt[3] = {0.,0.,0.};
   if (dir[2]<0) {
      distz=(-fDz-point[2])/dir[2];
      pt[2] = -fDz;
   } else {
      if (dir[2]>0) distz=(fDz-point[2])/dir[2];
      pt[2] = fDz;
   }
   for (i=0; i<4; i++) {
      dist=DistToPlane(point, dir, i, kTRUE);
      if (dist<distl) distl = dist;
   }
   if (distz<distl) {
      pt[0] = point[0]+distz*dir[0];
      pt[1] = point[1]+distz*dir[1];
      if (!Contains(pt)) distz = TGeoShape::Big();
   }
   dist = TMath::Min(distz, distl);
   if (dist<0 || dist>1.E10) return 0.;
   return dist;
#ifdef OLDALGORITHM
//#else
// compute distance to plane ipl :
// ipl=0 : points 0,4,1,5
// ipl=1 : points 1,5,2,6
// ipl=2 : points 2,6,3,7
// ipl=3 : points 3,7,0,4
   Double_t distmin;
   Bool_t lateral_cross = kFALSE;
   if (dir[2]<0) {
      distmin=(-fDz-point[2])/dir[2];
   } else {
      if (dir[2]>0) distmin =(fDz-point[2])/dir[2];
      else          distmin = TGeoShape::Big();
   }
   Double_t dz2 =0.5/fDz;
   Double_t xa,xb,xc,xd;
   Double_t ya,yb,yc,yd;
   Double_t eps = 100.*TGeoShape::Tolerance();
   for (Int_t ipl=0;ipl<4;ipl++) {
      Int_t j = (ipl+1)%4;
      xa=fXY[ipl][0];
      ya=fXY[ipl][1];
      xb=fXY[ipl+4][0];
      yb=fXY[ipl+4][1];
      xc=fXY[j][0];
      yc=fXY[j][1];
      xd=fXY[4+j][0];
      yd=fXY[4+j][1];

      Double_t tx1 =dz2*(xb-xa);
      Double_t ty1 =dz2*(yb-ya);
      Double_t tx2 =dz2*(xd-xc);
      Double_t ty2 =dz2*(yd-yc);
      Double_t dzp =fDz+point[2];
      Double_t xs1 =xa+tx1*dzp;
      Double_t ys1 =ya+ty1*dzp;
      Double_t xs2 =xc+tx2*dzp;
      Double_t ys2 =yc+ty2*dzp;
      Double_t dxs =xs2-xs1;
      Double_t dys =ys2-ys1;
      Double_t dtx =tx2-tx1;
      Double_t dty =ty2-ty1;
      Double_t a=(dtx*dir[1]-dty*dir[0]+(tx1*ty2-tx2*ty1)*dir[2])*dir[2];
      Double_t b=dxs*dir[1]-dys*dir[0]+(dtx*point[1]-dty*point[0]+ty2*xs1-ty1*xs2
              +tx1*ys2-tx2*ys1)*dir[2];
      Double_t c=dxs*point[1]-dys*point[0]+xs1*ys2-xs2*ys1;
      Double_t s=TGeoShape::Big();
      if (TMath::Abs(a)<eps) {
         if (TMath::Abs(b)<eps) continue;
         s=-c/b;
         if (s>eps && s < distmin) {
            distmin =s;
            lateral_cross=kTRUE;
         }
         continue;
      }
      Double_t d=b*b-4*a*c;
      if (d>=0.) {
         if (a>0) s=0.5*(-b-TMath::Sqrt(d))/a;
         else     s=0.5*(-b+TMath::Sqrt(d))/a;
         if (s>eps) {
            if (s < distmin) {
               distmin = s;
               lateral_cross = kTRUE;
            }
         } else {
            if (a>0) s=0.5*(-b+TMath::Sqrt(d))/a;
            else     s=0.5*(-b-TMath::Sqrt(d))/a;
            if (s>eps && s < distmin) {
               distmin =s;
               lateral_cross = kTRUE;
            }
         }
      }
   }

   if (!lateral_cross) {
      // We have to make sure that track crosses the top or bottom.
      if (distmin > 1.E10) return TGeoShape::Tolerance();
      Double_t pt[2];
      pt[0] = point[0]+distmin*dir[0];
      pt[1] = point[1]+distmin*dir[1];
      // Check if propagated point is in the polygon
      Double_t poly[8];
      Int_t i = 0;
      if (dir[2]>0.) i=4;
      for (Int_t j=0; j<4; j++) {
         poly[2*j]   = fXY[j+i][0];
         poly[2*j+1] = fXY[j+i][1];
      }
      if (!InsidePolygon(pt[0],pt[1],poly)) return TGeoShape::Tolerance();
   }
   return distmin;
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Divide this shape along one axis.

TGeoVolume *TGeoArb8::Divide(TGeoVolume *voldiv, const char * /*divname*/, Int_t /*iaxis*/, Int_t /*ndiv*/,
                             Double_t /*start*/, Double_t /*step*/)
{
   Error("Divide", "Division of an arbitrary trapezoid not implemented");
   return voldiv;
}

////////////////////////////////////////////////////////////////////////////////
/// Get shape range on a given axis.

Double_t TGeoArb8::GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const
{
   xlo = 0;
   xhi = 0;
   Double_t dx = 0;
   if (iaxis==3) {
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

void TGeoArb8::GetBoundingCylinder(Double_t *param) const
{
   // first compute rmin/rmax
   Double_t rmaxsq = 0;
   Double_t rsq;
   Int_t i;
   for (i=0; i<8; i++) {
      rsq = fXY[i][0]*fXY[i][0] + fXY[i][1]*fXY[i][1];
      rmaxsq = TMath::Max(rsq, rmaxsq);
   }
   param[0] = 0.;                  // Rmin
   param[1] = rmaxsq;              // Rmax
   param[2] = 0.;                  // Phi1
   param[3] = 360.;                // Phi2
}

////////////////////////////////////////////////////////////////////////////////
/// Fills real parameters of a positioned box inside this arb8. Returns 0 if successful.

Int_t TGeoArb8::GetFittingBox(const TGeoBBox *parambox, TGeoMatrix *mat, Double_t &dx, Double_t &dy, Double_t &dz) const
{
   dx=dy=dz=0;
   if (mat->IsRotation()) {
      Error("GetFittingBox", "cannot handle parametrized rotated volumes");
      return 1; // ### rotation not accepted ###
   }
   //--> translate the origin of the parametrized box to the frame of this box.
   Double_t origin[3];
   mat->LocalToMaster(parambox->GetOrigin(), origin);
   if (!Contains(origin)) {
      Error("GetFittingBox", "wrong matrix - parametrized box is outside this");
      return 1; // ### wrong matrix ###
   }
   //--> now we have to get the valid range for all parametrized axis
   Double_t dd[3];
   dd[0] = parambox->GetDX();
   dd[1] = parambox->GetDY();
   dd[2] = parambox->GetDZ();
   //-> check if Z range is fixed
   if (dd[2]<0) {
      dd[2] = TMath::Min(origin[2]+fDz, fDz-origin[2]);
      if (dd[2]<0) {
         Error("GetFittingBox", "wrong matrix");
         return 1;
      }
   }
   if (dd[0]>=0 && dd[1]>=0) {
      dx = dd[0];
      dy = dd[1];
      dz = dd[2];
      return 0;
   }
   //-> check now vertices at Z = origin[2] +/- dd[2]
   Double_t upper[8];
   Double_t lower[8];
   SetPlaneVertices(origin[2]-dd[2], lower);
   SetPlaneVertices(origin[2]+dd[2], upper);
   Double_t ddmin=TGeoShape::Big();
   for (Int_t iaxis=0; iaxis<2; iaxis++) {
      if (dd[iaxis]>=0) continue;
      ddmin=TGeoShape::Big();
      for (Int_t ivert=0; ivert<4; ivert++) {
         ddmin = TMath::Min(ddmin, TMath::Abs(origin[iaxis]-lower[2*ivert+iaxis]));
         ddmin = TMath::Min(ddmin, TMath::Abs(origin[iaxis]-upper[2*ivert+iaxis]));
      }
      dd[iaxis] = ddmin;
   }
   dx = dd[0];
   dy = dd[1];
   dz = dd[2];
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Computes normal to plane defined by P1, P2 and P3

void TGeoArb8::GetPlaneNormal(Double_t *p1, Double_t *p2, Double_t *p3, Double_t *norm)
{
   Double_t cross = 0.;
   Double_t v1[3], v2[3];
   Int_t i;
   for (i=0; i<3; i++) {
      v1[i] = p2[i] - p1[i];
      v2[i] = p3[i] - p1[i];
   }
   norm[0] = v1[1]*v2[2]-v1[2]*v2[1];
   cross += norm[0]*norm[0];
   norm[1] = v1[2]*v2[0]-v1[0]*v2[2];
   cross += norm[1]*norm[1];
   norm[2] = v1[0]*v2[1]-v1[1]*v2[0];
   cross += norm[2]*norm[2];
   if (TMath::Abs(cross) < TGeoShape::Tolerance()) return;
   cross = 1./TMath::Sqrt(cross);
   for (i=0; i<3; i++) norm[i] *= cross;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills array with n random points located on the surface of indexed facet.
/// The output array must be provided with a length of minimum 3*npoints. Returns
/// true if operation succeeded.
/// Possible index values:
///  - 0 - all facets together
///  - 1 to 6 - facet index from bottom to top Z

Bool_t TGeoArb8::GetPointsOnFacet(Int_t /*index*/, Int_t /*npoints*/, Double_t * /* array */) const
{
   return kFALSE;
/*
   if (index<0 || index>6) return kFALSE;
   if (index==0) {
      // Just generate same number of points on each facet
      Int_t npts = npoints/6.;
      Int_t count = 0;
      for (Int_t ifacet=0; ifacet<6; ifacet++) {
         if (GetPointsOnFacet(ifacet+1, npts, &array[3*count])) count += npts;
         if (ifacet<5) npts = (npoints-count)/(5.-ifacet);
      }
      if (count>0) return kTRUE;
      return kFALSE;
   }
   Double_t z, cf;
   Double_t xmin=TGeoShape::Big();
   Double_t xmax=-xmin;
   Double_t ymin=TGeoShape::Big();
   Double_t ymax=-ymin;
   Double_t dy=0.;
   Double_t poly[8];
   Double_t point[2];
   Int_t i;
   if (index==1 || index==6) {
      z = (index==1)?-fDz:fDz;
      cf = 0.5*(fDz-z)/fDz;
      for (i=0; i<4; i++) {
         poly[2*i]   = fXY[i+4][0]+cf*(fXY[i][0]-fXY[i+4][0]);
         poly[2*i+1] = fXY[i+4][1]+cf*(fXY[i][1]-fXY[i+4][1]);
         xmin = TMath::Min(xmin, poly[2*i]);
         xmax = TMath::Max(xmax, poly[2*i]);
         ymin = TMath::Min(ymin, poly[2*i]);
         ymax = TMath::Max(ymax, poly[2*i]);
      }
   }
   Int_t nshoot = 0;
   Int_t nmiss = 0;
   for (i=0; i<npoints; i++) {
      Double_t *point = &array[3*i];
      switch (surfindex) {
         case 1:
         case 6:
            while (nmiss<1000) {
               point[0] = xmin + (xmax-xmin)*gRandom->Rndm();
               point[1] = ymin + (ymax-ymin)*gRandom->Rndm();
            }

   return InsidePolygon(point[0],point[1],poly);
*/
}

////////////////////////////////////////////////////////////////////////////////
/// Finds if a point in XY plane is inside the polygon defines by PTS.

Bool_t TGeoArb8::InsidePolygon(Double_t x, Double_t y, Double_t *pts)
{
   Int_t i,j;
   Double_t x1,y1,x2,y2;
   Double_t cross;
   for (i=0; i<4; i++) {
      j = (i+1)%4;
      x1 = pts[i<<1];
      y1 = pts[(i<<1)+1];
      x2 = pts[j<<1];
      y2 = pts[(j<<1)+1];
      cross = (x-x1)*(y2-y1)-(y-y1)*(x2-x1);
      if (cross<0) return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Prints shape parameters

void TGeoArb8::InspectShape() const
{
   printf("*** Shape %s: TGeoArb8 ***\n", GetName());
   if (IsTwisted()) printf("  = TWISTED\n");
   for (Int_t ip=0; ip<8; ip++) {
      printf("    point #%i : x=%11.5f y=%11.5f z=%11.5f\n",
             ip, fXY[ip][0], fXY[ip][1], fDz*((ip<4)?-1:1));
   }
   printf(" Bounding box:\n");
   TGeoBBox::InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Computes the closest distance from given point to this shape.

Double_t TGeoArb8::Safety(const Double_t *point, Bool_t in) const
{
   Double_t safz = fDz-TMath::Abs(point[2]);
   if (!in) safz = -safz;
   Int_t iseg;
   Double_t safe = TGeoShape::Big();
   Double_t lsq, ssq, dx, dy, dpx, dpy, u;
   if (IsTwisted()) {
      if (!in) {
         if (!TGeoBBox::Contains(point)) return TGeoBBox::Safety(point,kFALSE);
      }
      // Point is also in the bounding box ;-(
      // Compute closest distance to any segment
      Double_t vert[8];
      Double_t *p1, *p2;
      Int_t isegmin=0;
      Double_t umin = 0.;
      SetPlaneVertices (point[2], vert);
      for (iseg=0; iseg<4; iseg++) {
         if (safe<TGeoShape::Tolerance()) return 0.;
         p1 = &vert[2*iseg];
         p2 = &vert[2*((iseg+1)%4)];
         dx = p2[0] - p1[0];
         dy = p2[1] - p1[1];
         dpx = point[0] - p1[0];
         dpy = point[1] - p1[1];

         lsq = dx*dx + dy*dy;
         u = (dpx*dx + dpy*dy)/lsq;
         if (u>1) {
            dpx = point[0]-p2[0];
            dpy = point[1]-p2[1];
         } else {
            if (u>=0) {
               dpx -= u*dx;
               dpy -= u*dy;
            }
         }
         ssq = dpx*dpx + dpy*dpy;
         if (ssq < safe) {
            isegmin = iseg;
            umin = u;
            safe = ssq;
         }
      }
      if (umin<0) umin = 0.;
      if (umin>1) {
         isegmin = (isegmin+1)%4;
         umin = 0.;
      }
      Int_t i1 = isegmin;
      Int_t i2 = (isegmin+1)%4;
      Double_t dx1 = fXY[i2][0]-fXY[i1][0];
      Double_t dx2 = fXY[i2+4][0]-fXY[i1+4][0];
      Double_t dy1 = fXY[i2][1]-fXY[i1][1];
      Double_t dy2 = fXY[i2+4][1]-fXY[i1+4][1];
      dx = dx1 + umin*(dx2-dx1);
      dy = dy1 + umin*(dy2-dy1);
      safe *= 1.- 4.*fDz*fDz/(dx*dx+dy*dy+4.*fDz*fDz);
      safe = TMath::Sqrt(safe);
      if (in) return TMath::Min(safz,safe);
      return TMath::Max(safz,safe);
   }

   Double_t saf[5];
   saf[0] = safz;

   for (iseg=0; iseg<4; iseg++) saf[iseg+1] = SafetyToFace(point,iseg,in);
   if (in) safe = saf[TMath::LocMin(5, saf)];
   else    safe = saf[TMath::LocMax(5, saf)];
   if (safe<0) return 0.;
   return safe;
}

////////////////////////////////////////////////////////////////////////////////
/// Estimate safety to lateral plane defined by segment iseg in range [0,3]
/// Might be negative: plane seen only from inside.

Double_t TGeoArb8::SafetyToFace(const Double_t *point, Int_t iseg, Bool_t in) const
{
   Double_t vertices[12];
   Int_t ipln = (iseg+1)%4;
   // point 1
   vertices[0] = fXY[iseg][0];
   vertices[1] = fXY[iseg][1];
   vertices[2] = -fDz;
   // point 2
   vertices[3] = fXY[ipln][0];
   vertices[4] = fXY[ipln][1];
   vertices[5] = -fDz;
   // point 3
   vertices[6] = fXY[ipln+4][0];
   vertices[7] = fXY[ipln+4][1];
   vertices[8] = fDz;
   // point 4
   vertices[9] = fXY[iseg+4][0];
   vertices[10] = fXY[iseg+4][1];
   vertices[11] = fDz;
   Double_t safe;
   Double_t norm[3];
   Double_t *p1, *p2, *p3;
   p1 = &vertices[0];
   p2 = &vertices[9];
   p3 = &vertices[6];
   if (IsSamePoint(p2,p3)) {
      p3 = &vertices[3];
      if (IsSamePoint(p1,p3)) return -TGeoShape::Big(); // skip single segment
   }
   GetPlaneNormal(p1,p2,p3,norm);
   safe = (point[0]-p1[0])*norm[0]+(point[1]-p1[1])*norm[1]+(point[2]-p1[2])*norm[2];
   if (in) return (-safe);
   return safe;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoArb8::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive)) return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   dz       = " << fDz << ";" << std::endl;
   out << "   vert[0]  = " << fXY[0][0] << ";" << std::endl;
   out << "   vert[1]  = " << fXY[0][1] << ";" << std::endl;
   out << "   vert[2]  = " << fXY[1][0] << ";" << std::endl;
   out << "   vert[3]  = " << fXY[1][1] << ";" << std::endl;
   out << "   vert[4]  = " << fXY[2][0] << ";" << std::endl;
   out << "   vert[5]  = " << fXY[2][1] << ";" << std::endl;
   out << "   vert[6]  = " << fXY[3][0] << ";" << std::endl;
   out << "   vert[7]  = " << fXY[3][1] << ";" << std::endl;
   out << "   vert[8]  = " << fXY[4][0] << ";" << std::endl;
   out << "   vert[9]  = " << fXY[4][1] << ";" << std::endl;
   out << "   vert[10] = " << fXY[5][0] << ";" << std::endl;
   out << "   vert[11] = " << fXY[5][1] << ";" << std::endl;
   out << "   vert[12] = " << fXY[6][0] << ";" << std::endl;
   out << "   vert[13] = " << fXY[6][1] << ";" << std::endl;
   out << "   vert[14] = " << fXY[7][0] << ";" << std::endl;
   out << "   vert[15] = " << fXY[7][1] << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoArb8(\"" << GetName() << "\", dz,vert);" << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Computes intersection points between plane at zpl and non-horizontal edges.

void TGeoArb8::SetPlaneVertices(Double_t zpl, Double_t *vertices) const
{
   Double_t cf = 0.5*(fDz-zpl)/fDz;
   for (Int_t i=0; i<4; i++) {
      vertices[2*i]   = fXY[i+4][0]+cf*(fXY[i][0]-fXY[i+4][0]);
      vertices[2*i+1] = fXY[i+4][1]+cf*(fXY[i][1]-fXY[i+4][1]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set all arb8 params in one step.
/// param[0] = dz
/// param[1] = x0
/// param[2] = y0
/// ...

void TGeoArb8::SetDimensions(Double_t *param)
{
   fDz      = param[0];
   for (Int_t i=0; i<8; i++) {
      fXY[i][0] = param[2*i+1];
      fXY[i][1] = param[2*i+2];
   }
   ComputeTwist();
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates arb8 mesh points

void TGeoArb8::SetPoints(Double_t *points) const
{
   for (Int_t i=0; i<8; i++) {
      points[3*i] = fXY[i][0];
      points[3*i+1] = fXY[i][1];
      points[3*i+2] = (i<4)?-fDz:fDz;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Creates arb8 mesh points

void TGeoArb8::SetPoints(Float_t *points) const
{
   for (Int_t i=0; i<8; i++) {
      points[3*i] = fXY[i][0];
      points[3*i+1] = fXY[i][1];
      points[3*i+2] = (i<4)?-fDz:fDz;
   }
}

////////////////////////////////////////////////////////////////////////////////
///  Set values for a given vertex.

void TGeoArb8::SetVertex(Int_t vnum, Double_t x, Double_t y)
{
   if (vnum<0 || vnum >7) {
      Error("SetVertex", "Invalid vertex number");
      return;
   }
   fXY[vnum][0] = x;
   fXY[vnum][1] = y;
   if (vnum == 7) {
      ComputeTwist();
      ComputeBBox();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill size of this 3-D object

void TGeoArb8::Sizeof3D() const
{
   TGeoBBox::Sizeof3D();
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TGeoManager.

void TGeoArb8::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TGeoArb8::Class(), this);
      ComputeTwist();
   } else {
      R__b.WriteClassBuffer(TGeoArb8::Class(), this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check the inside status for each of the points in the array.
/// Input: Array of point coordinates + vector size
/// Output: Array of Booleans for the inside of each point

void TGeoArb8::Contains_v(const Double_t *points, Bool_t *inside, Int_t vecsize) const
{
   for (Int_t i=0; i<vecsize; i++) inside[i] = Contains(&points[3*i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the normal for an array o points so that norm.dot.dir is positive
/// Input: Arrays of point coordinates and directions + vector size
/// Output: Array of normal directions

void TGeoArb8::ComputeNormal_v(const Double_t *points, const Double_t *dirs, Double_t *norms, Int_t vecsize)
{
   for (Int_t i=0; i<vecsize; i++) ComputeNormal(&points[3*i], &dirs[3*i], &norms[3*i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoArb8::DistFromInside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromInside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoArb8::DistFromOutside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromOutside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute safe distance from each of the points in the input array.
/// Input: Array of point coordinates, array of statuses for these points, size of the arrays
/// Output: Safety values

void TGeoArb8::Safety_v(const Double_t *points, const Bool_t *inside, Double_t *safe, Int_t vecsize) const
{
   for (Int_t i=0; i<vecsize; i++) safe[i] = Safety(&points[3*i], inside[i]);
}

ClassImp(TGeoTrap);

////////////////////////////////////////////////////////////////////////////////
/// Default ctor

TGeoTrap::TGeoTrap()
{
   fDz = 0;
   fTheta = 0;
   fPhi = 0;
   fH1 = fH2 = fBl1 = fBl2 = fTl1 = fTl2 = fAlpha1 = fAlpha2 = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor providing just a range in Z, theta and phi.

TGeoTrap::TGeoTrap(Double_t dz, Double_t theta, Double_t phi)
         :TGeoArb8("", 0, 0)
{
   fDz = dz;
   fTheta = theta;
   fPhi = phi;
   fH1 = fH2 = fBl1 = fBl2 = fTl1 = fTl2 = fAlpha1 = fAlpha2 = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor.

TGeoTrap::TGeoTrap(Double_t dz, Double_t theta, Double_t phi, Double_t h1,
              Double_t bl1, Double_t tl1, Double_t alpha1, Double_t h2, Double_t bl2,
              Double_t tl2, Double_t alpha2)
         :TGeoArb8("", 0, 0)
{
   fDz = dz;
   fTheta = theta;
   fPhi = phi;
   fH1 = h1;
   fH2 = h2;
   fBl1 = bl1;
   fBl2 = bl2;
   fTl1 = tl1;
   fTl2 = tl2;
   fAlpha1 = alpha1;
   fAlpha2 = alpha2;
   Double_t tx = TMath::Tan(theta*TMath::DegToRad())*TMath::Cos(phi*TMath::DegToRad());
   Double_t ty = TMath::Tan(theta*TMath::DegToRad())*TMath::Sin(phi*TMath::DegToRad());
   Double_t ta1 = TMath::Tan(alpha1*TMath::DegToRad());
   Double_t ta2 = TMath::Tan(alpha2*TMath::DegToRad());
   fXY[0][0] = -dz*tx-h1*ta1-bl1;    fXY[0][1] = -dz*ty-h1;
   fXY[1][0] = -dz*tx+h1*ta1-tl1;    fXY[1][1] = -dz*ty+h1;
   fXY[2][0] = -dz*tx+h1*ta1+tl1;    fXY[2][1] = -dz*ty+h1;
   fXY[3][0] = -dz*tx-h1*ta1+bl1;    fXY[3][1] = -dz*ty-h1;
   fXY[4][0] = dz*tx-h2*ta2-bl2;    fXY[4][1] = dz*ty-h2;
   fXY[5][0] = dz*tx+h2*ta2-tl2;    fXY[5][1] = dz*ty+h2;
   fXY[6][0] = dz*tx+h2*ta2+tl2;    fXY[6][1] = dz*ty+h2;
   fXY[7][0] = dz*tx-h2*ta2+bl2;    fXY[7][1] = dz*ty-h2;
   ComputeTwist();
   if ((dz<0) || (h1<0) || (bl1<0) || (tl1<0) ||
       (h2<0) || (bl2<0) || (tl2<0)) {
      SetShapeBit(kGeoRunTimeShape);
   }
   else TGeoArb8::ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor with name.

TGeoTrap::TGeoTrap(const char *name, Double_t dz, Double_t theta, Double_t phi, Double_t h1,
              Double_t bl1, Double_t tl1, Double_t alpha1, Double_t h2, Double_t bl2,
              Double_t tl2, Double_t alpha2)
         :TGeoArb8(name, 0, 0)
{
   SetName(name);
   fDz = dz;
   fTheta = theta;
   fPhi = phi;
   fH1 = h1;
   fH2 = h2;
   fBl1 = bl1;
   fBl2 = bl2;
   fTl1 = tl1;
   fTl2 = tl2;
   fAlpha1 = alpha1;
   fAlpha2 = alpha2;
   for (Int_t i=0; i<8; i++) {
      fXY[i][0] = 0.0;
      fXY[i][1] = 0.0;
   }
   Double_t tx = TMath::Tan(theta*TMath::DegToRad())*TMath::Cos(phi*TMath::DegToRad());
   Double_t ty = TMath::Tan(theta*TMath::DegToRad())*TMath::Sin(phi*TMath::DegToRad());
   Double_t ta1 = TMath::Tan(alpha1*TMath::DegToRad());
   Double_t ta2 = TMath::Tan(alpha2*TMath::DegToRad());
   fXY[0][0] = -dz*tx-h1*ta1-bl1;    fXY[0][1] = -dz*ty-h1;
   fXY[1][0] = -dz*tx+h1*ta1-tl1;    fXY[1][1] = -dz*ty+h1;
   fXY[2][0] = -dz*tx+h1*ta1+tl1;    fXY[2][1] = -dz*ty+h1;
   fXY[3][0] = -dz*tx-h1*ta1+bl1;    fXY[3][1] = -dz*ty-h1;
   fXY[4][0] = dz*tx-h2*ta2-bl2;    fXY[4][1] = dz*ty-h2;
   fXY[5][0] = dz*tx+h2*ta2-tl2;    fXY[5][1] = dz*ty+h2;
   fXY[6][0] = dz*tx+h2*ta2+tl2;    fXY[6][1] = dz*ty+h2;
   fXY[7][0] = dz*tx-h2*ta2+bl2;    fXY[7][1] = dz*ty-h2;
   ComputeTwist();
   if ((dz<0) || (h1<0) || (bl1<0) || (tl1<0) ||
       (h2<0) || (bl2<0) || (tl2<0)) {
      SetShapeBit(kGeoRunTimeShape);
   }
   else TGeoArb8::ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGeoTrap::~TGeoTrap()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from inside point to surface of the trapezoid

Double_t TGeoTrap::DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
{
   if (iact<3 && safe) {
      // compute safe distance
      *safe = Safety(point, kTRUE);
      if (iact==0) return TGeoShape::Big();
      if (iact==1 && step<*safe) return TGeoShape::Big();
   }
   // compute distance to get outside this shape
   // return TGeoArb8::DistFromInside(point, dir, iact, step, safe);

// compute distance to plane ipl :
// ipl=0 : points 0,4,1,5
// ipl=1 : points 1,5,2,6
// ipl=2 : points 2,6,3,7
// ipl=3 : points 3,7,0,4
   Double_t distmin;
   if (dir[2]<0) {
      distmin=(-fDz-point[2])/dir[2];
   } else {
      if (dir[2]>0) distmin =(fDz-point[2])/dir[2];
      else          distmin = TGeoShape::Big();
   }
   Double_t xa,xb,xc;
   Double_t ya,yb,yc;
   for (Int_t ipl=0;ipl<4;ipl++) {
      Int_t j = (ipl+1)%4;
      xa=fXY[ipl][0];
      ya=fXY[ipl][1];
      xb=fXY[ipl+4][0];
      yb=fXY[ipl+4][1];
      xc=fXY[j][0];
      yc=fXY[j][1];
      Double_t ax,ay,az;
      ax = xb-xa;
      ay = yb-ya;
      az = 2.*fDz;
      Double_t bx,by;
      bx = xc-xa;
      by = yc-ya;
      Double_t ddotn = -dir[0]*az*by + dir[1]*az*bx+dir[2]*(ax*by-ay*bx);
      if (ddotn<=0) continue; // entering
      Double_t saf = -(point[0]-xa)*az*by + (point[1]-ya)*az*bx + (point[2]+fDz)*(ax*by-ay*bx);
      if (saf>=0.0) return 0.0;
      Double_t s = -saf/ddotn;
      if (s<distmin) distmin=s;
   }
   return distmin;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from outside point to surface of the trapezoid

Double_t TGeoTrap::DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
{
   if (iact<3 && safe) {
      // compute safe distance
      *safe = Safety(point, kFALSE);
      if (iact==0) return TGeoShape::Big();
      if (iact==1 && step<*safe) return TGeoShape::Big();
   }
   // Check if the bounding box is crossed within the requested distance
   Double_t sdist = TGeoBBox::DistFromOutside(point,dir, fDX, fDY, fDZ, fOrigin, step);
   if (sdist>=step) return TGeoShape::Big();
   // compute distance to get outside this shape
   Bool_t in = kTRUE;
   Double_t pts[8];
   Double_t snxt;
   Double_t xnew, ynew, znew;
   Int_t i,j;
   if (point[2]<-fDz+TGeoShape::Tolerance()) {
      if (dir[2] < TGeoShape::Tolerance()) return TGeoShape::Big();
      in = kFALSE;
      snxt = -(fDz+point[2])/dir[2];
      xnew = point[0] + snxt*dir[0];
      ynew = point[1] + snxt*dir[1];
      for (i=0;i<4;i++) {
         j = i<<1;
         pts[j] = fXY[i][0];
         pts[j+1] = fXY[i][1];
      }
      if (InsidePolygon(xnew,ynew,pts)) return snxt;
   } else if (point[2]>fDz-TGeoShape::Tolerance()) {
      if (dir[2] > -TGeoShape::Tolerance()) return TGeoShape::Big();
      in = kFALSE;
      snxt = (fDz-point[2])/dir[2];
      xnew = point[0] + snxt*dir[0];
      ynew = point[1] + snxt*dir[1];
      for (i=0;i<4;i++) {
         j = i<<1;
         pts[j] = fXY[i+4][0];
         pts[j+1] = fXY[i+4][1];
      }
      if (InsidePolygon(xnew,ynew,pts)) return snxt;
   }
   snxt = TGeoShape::Big();


   // check lateral faces
   Double_t dz2 =0.5/fDz;
   Double_t xa,xb,xc,xd;
   Double_t ya,yb,yc,yd;
   Double_t ax,ay,az;
   Double_t bx,by;
   Double_t ddotn, saf;
   Double_t safmin = TGeoShape::Big();
   Bool_t exiting = kFALSE;
   for (i=0; i<4; i++) {
      j = (i+1)%4;
      xa=fXY[i][0];
      ya=fXY[i][1];
      xb=fXY[i+4][0];
      yb=fXY[i+4][1];
      xc=fXY[j][0];
      yc=fXY[j][1];
      xd=fXY[4+j][0];
      yd=fXY[4+j][1];
      ax = xb-xa;
      ay = yb-ya;
      az = 2.*fDz;
      bx = xc-xa;
      by = yc-ya;
      ddotn = -dir[0]*az*by + dir[1]*az*bx+dir[2]*(ax*by-ay*bx);
      saf = (point[0]-xa)*az*by - (point[1]-ya)*az*bx - (point[2]+fDz)*(ax*by-ay*bx);

      if (saf<=0) {
         // face visible from point outside
         in = kFALSE;
         if (ddotn>=0) return TGeoShape::Big();
         snxt = saf/ddotn;
         znew = point[2]+snxt*dir[2];
         if (TMath::Abs(znew)<=fDz) {
            xnew = point[0]+snxt*dir[0];
            ynew = point[1]+snxt*dir[1];
            Double_t tx1 =dz2*(xb-xa);
            Double_t ty1 =dz2*(yb-ya);
            Double_t tx2 =dz2*(xd-xc);
            Double_t ty2 =dz2*(yd-yc);
            Double_t dzp =fDz+znew;
            Double_t xs1 =xa+tx1*dzp;
            Double_t ys1 =ya+ty1*dzp;
            Double_t xs2 =xc+tx2*dzp;
            Double_t ys2 =yc+ty2*dzp;
            if (TMath::Abs(xs1-xs2)>TMath::Abs(ys1-ys2)) {
               if ((xnew-xs1)*(xs2-xnew)>=0) return snxt;
            } else {
               if ((ynew-ys1)*(ys2-ynew)>=0) return snxt;
            }
         }
      } else {
         if (saf<safmin) {
            safmin = saf;
            if (ddotn>=0) exiting = kTRUE;
            else exiting = kFALSE;
         }
      }
   }
   // Check also Z boundaries (point may be inside and close to Z) - Christian Hammann
   saf = fDz - TMath::Abs(point[2]);
   if (saf>0 && saf<safmin) exiting = (point[2]*dir[2] > 0)?kTRUE:kFALSE;
   if (!in) return TGeoShape::Big();
   if (exiting) return TGeoShape::Big();
   return 0.0;
}

////////////////////////////////////////////////////////////////////////////////
/// Divide this trapezoid shape belonging to volume "voldiv" into ndiv volumes
/// called divname, from start position with the given step. Only Z divisions
/// are supported. For Z divisions just return the pointer to the volume to be
/// divided. In case a wrong division axis is supplied, returns pointer to
/// volume that was divided.

TGeoVolume *TGeoTrap::Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv,
                             Double_t start, Double_t step)
{
   TGeoShape *shape;           //--- shape to be created
   TGeoVolume *vol;            //--- division volume to be created
   TGeoVolumeMulti *vmulti;    //--- generic divided volume
   TGeoPatternFinder *finder;  //--- finder to be attached
   TString opt = "";           //--- option to be attached
   if (iaxis!=3) {
      Error("Divide", "cannot divide trapezoids on other axis than Z");
      return 0;
   }
   Double_t end = start+ndiv*step;
   Double_t points_lo[8];
   Double_t points_hi[8];
   finder = new TGeoPatternTrapZ(voldiv, ndiv, start, end);
   voldiv->SetFinder(finder);
   finder->SetDivIndex(voldiv->GetNdaughters());
   opt = "Z";
   vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
   Double_t txz = ((TGeoPatternTrapZ*)finder)->GetTxz();
   Double_t tyz = ((TGeoPatternTrapZ*)finder)->GetTyz();
   Double_t zmin, zmax, ox,oy,oz;
   for (Int_t idiv=0; idiv<ndiv; idiv++) {
      zmin = start+idiv*step;
      zmax = start+(idiv+1)*step;
      oz = start+idiv*step+step/2;
      ox = oz*txz;
      oy = oz*tyz;
      SetPlaneVertices(zmin, &points_lo[0]);
      SetPlaneVertices(zmax, &points_hi[0]);
      shape = new TGeoTrap(step/2, fTheta, fPhi);
      for (Int_t vert1=0; vert1<4; vert1++)
         ((TGeoArb8*)shape)->SetVertex(vert1, points_lo[2*vert1]-ox, points_lo[2*vert1+1]-oy);
      for (Int_t vert2=0; vert2<4; vert2++)
         ((TGeoArb8*)shape)->SetVertex(vert2+4, points_hi[2*vert2]-ox, points_hi[2*vert2+1]-oy);
      vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
      vmulti->AddVolume(vol);
      voldiv->AddNodeOffset(vol, idiv, oz, opt.Data());
      ((TGeoNodeOffset*)voldiv->GetNodes()->At(voldiv->GetNdaughters()-1))->SetFinder(finder);
   }
   return vmulti;
}

////////////////////////////////////////////////////////////////////////////////
/// In case shape has some negative parameters, these have to be computed
/// in order to fit the mother.

TGeoShape *TGeoTrap::GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix * /*mat*/) const
{
   if (!TestShapeBit(kGeoRunTimeShape)) return 0;
   if (mother->IsRunTimeShape()) {
      Error("GetMakeRuntimeShape", "invalid mother");
      return 0;
   }
   Double_t dz, h1, bl1, tl1, h2, bl2, tl2;
   if (fDz<0) dz=((TGeoTrap*)mother)->GetDz();
   else dz=fDz;

   if (fH1<0)  h1 = ((TGeoTrap*)mother)->GetH1();
   else        h1 = fH1;

   if (fH2<0)  h2 = ((TGeoTrap*)mother)->GetH2();
   else        h2 = fH2;

   if (fBl1<0) bl1 = ((TGeoTrap*)mother)->GetBl1();
   else        bl1 = fBl1;

   if (fBl2<0) bl2 = ((TGeoTrap*)mother)->GetBl2();
   else        bl2 = fBl2;

   if (fTl1<0) tl1 = ((TGeoTrap*)mother)->GetTl1();
   else        tl1 = fTl1;

   if (fTl2<0) tl2 = ((TGeoTrap*)mother)->GetTl2();
   else        tl2 = fTl2;

   return (new TGeoTrap(dz, fTheta, fPhi, h1, bl1, tl1, fAlpha1, h2, bl2, tl2, fAlpha2));
}

////////////////////////////////////////////////////////////////////////////////
/// Computes the closest distance from given point to this shape.

Double_t TGeoTrap::Safety(const Double_t *point, Bool_t in) const
{
   Double_t safe = TGeoShape::Big();
   Double_t saf[5];
   Double_t norm[3]; // normal to current facette
   Int_t i, j;       // current facette index
   Double_t x0, y0, z0=-fDz, x1, y1, z1=fDz, x2, y2;
   Double_t ax, ay, az=z1-z0, bx, by;
   Double_t fn;
   //---> compute safety for lateral planes
   for (i=0; i<4; i++) {
      if (in) saf[i] = TGeoShape::Big();
      else    saf[i] = 0.;
      x0 = fXY[i][0];
      y0 = fXY[i][1];
      x1 = fXY[i+4][0];
      y1 = fXY[i+4][1];
      ax = x1-x0;
      ay = y1-y0;
      az = z1-z0;
      j  = (i+1)%4;
      x2 = fXY[j][0];
      y2 = fXY[j][1];
      bx = x2-x0;
      by = y2-y0;
      if (TMath::Abs(bx)<TGeoShape::Tolerance() && TMath::Abs(by)<TGeoShape::Tolerance()) {
         x2 = fXY[4+j][0];
         y2 = fXY[4+j][1];
         bx = x2-x1;
         by = y2-y1;
         if (TMath::Abs(bx)<TGeoShape::Tolerance() && TMath::Abs(by)<TGeoShape::Tolerance()) continue;
      }
      norm[0] = -az*by;
      norm[1] = az*bx;
      norm[2] = ax*by-ay*bx;
      fn = TMath::Sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
      if (fn<1E-10) continue;
      saf[i] = (x0-point[0])*norm[0]+(y0-point[1])*norm[1]+(-fDz-point[2])*norm[2];
      if (in) {
         saf[i]=TMath::Abs(saf[i])/fn; // they should be all positive anyway
      } else {
         saf[i] = -saf[i]/fn;   // only negative values are interesting
      }
   }
   saf[4] = fDz-TMath::Abs(point[2]);
   if (in) {
      safe = saf[0];
      for (j=1;j<5;j++) if (saf[j] <safe) safe = saf[j];
   } else {
      saf[4]=-saf[4];
      safe = saf[0];
      for (j=1;j<5;j++) if (saf[j] >safe) safe = saf[j];
   }
   return safe;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoTrap::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive)) return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   dz     = " << fDz << ";" << std::endl;
   out << "   theta  = " << fTheta << ";" << std::endl;
   out << "   phi    = " << fPhi << ";" << std::endl;
   out << "   h1     = " << fH1<< ";" << std::endl;
   out << "   bl1    = " << fBl1<< ";" << std::endl;
   out << "   tl1    = " << fTl1<< ";" << std::endl;
   out << "   alpha1 = " << fAlpha1 << ";" << std::endl;
   out << "   h2     = " << fH2 << ";" << std::endl;
   out << "   bl2    = " << fBl2<< ";" << std::endl;
   out << "   tl2    = " << fTl2<< ";" << std::endl;
   out << "   alpha2 = " << fAlpha2 << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoTrap(\"" << GetName() << "\", dz,theta,phi,h1,bl1,tl1,alpha1,h2,bl2,tl2,alpha2);" << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Set all arb8 params in one step.
///  - param[0] = dz
///  - param[1] = theta
///  - param[2] = phi
///  - param[3] = h1
///  - param[4] = bl1
///  - param[5] = tl1
///  - param[6] = alpha1
///  - param[7] = h2
///  - param[8] = bl2
///  - param[9] = tl2
///  - param[10] = alpha2

void TGeoTrap::SetDimensions(Double_t *param)
{
   fDz      = param[0];
   fTheta = param[1];
   fPhi   = param[2];
   fH1 = param[3];
   fH2 = param[7];
   fBl1 = param[4];
   fBl2 = param[8];
   fTl1 = param[5];
   fTl2 = param[9];
   fAlpha1 = param[6];
   fAlpha2 = param[10];
   Double_t tx = TMath::Tan(fTheta*TMath::DegToRad())*TMath::Cos(fPhi*TMath::DegToRad());
   Double_t ty = TMath::Tan(fTheta*TMath::DegToRad())*TMath::Sin(fPhi*TMath::DegToRad());
   Double_t ta1 = TMath::Tan(fAlpha1*TMath::DegToRad());
   Double_t ta2 = TMath::Tan(fAlpha2*TMath::DegToRad());
   fXY[0][0] = -fDz*tx-fH1*ta1-fBl1;    fXY[0][1] = -fDz*ty-fH1;
   fXY[1][0] = -fDz*tx+fH1*ta1-fTl1;    fXY[1][1] = -fDz*ty+fH1;
   fXY[2][0] = -fDz*tx+fH1*ta1+fTl1;    fXY[2][1] = -fDz*ty+fH1;
   fXY[3][0] = -fDz*tx-fH1*ta1+fBl1;    fXY[3][1] = -fDz*ty-fH1;
   fXY[4][0] = fDz*tx-fH2*ta2-fBl2;    fXY[4][1] = fDz*ty-fH2;
   fXY[5][0] = fDz*tx+fH2*ta2-fTl2;    fXY[5][1] = fDz*ty+fH2;
   fXY[6][0] = fDz*tx+fH2*ta2+fTl2;    fXY[6][1] = fDz*ty+fH2;
   fXY[7][0] = fDz*tx-fH2*ta2+fBl2;    fXY[7][1] = fDz*ty-fH2;
   ComputeTwist();
   if ((fDz<0) || (fH1<0) || (fBl1<0) || (fTl1<0) ||
       (fH2<0) || (fBl2<0) || (fTl2<0)) {
      SetShapeBit(kGeoRunTimeShape);
   }
   else TGeoArb8::ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoTrap::DistFromInside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromInside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoTrap::DistFromOutside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromOutside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute safe distance from each of the points in the input array.
/// Input: Array of point coordinates, array of statuses for these points, size of the arrays
/// Output: Safety values

void TGeoTrap::Safety_v(const Double_t *points, const Bool_t *inside, Double_t *safe, Int_t vecsize) const
{
   for (Int_t i=0; i<vecsize; i++) safe[i] = Safety(&points[3*i], inside[i]);
}

ClassImp(TGeoGtra);

////////////////////////////////////////////////////////////////////////////////
/// Default ctor

TGeoGtra::TGeoGtra()
{
   fTwistAngle = 0;

}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGeoGtra::TGeoGtra(Double_t dz, Double_t theta, Double_t phi, Double_t twist, Double_t h1,
              Double_t bl1, Double_t tl1, Double_t alpha1, Double_t h2, Double_t bl2,
              Double_t tl2, Double_t alpha2)
         :TGeoTrap(dz, theta, phi, h1, bl1, tl1, alpha1, h2, bl2, tl2, alpha2)
{
   fTwistAngle = twist;
   Double_t x,y;
   Double_t th = theta*TMath::DegToRad();
   Double_t ph = phi*TMath::DegToRad();
   // Coordinates of the center of the bottom face
   Double_t xc = -dz*TMath::Sin(th)*TMath::Cos(ph);
   Double_t yc = -dz*TMath::Sin(th)*TMath::Sin(ph);

   Int_t i;

   for (i=0; i<4; i++) {
      x = fXY[i][0] - xc;
      y = fXY[i][1] - yc;
      fXY[i][0] = x*TMath::Cos(-0.5*twist*TMath::DegToRad()) + y*TMath::Sin(-0.5*twist*TMath::DegToRad()) + xc;
      fXY[i][1] = -x*TMath::Sin(-0.5*twist*TMath::DegToRad()) + y*TMath::Cos(-0.5*twist*TMath::DegToRad()) + yc;
   }
   // Coordinates of the center of the top face
   xc = -xc;
   yc = -yc;
   for (i=4; i<8; i++) {
      x = fXY[i][0] - xc;
      y = fXY[i][1] - yc;
      fXY[i][0] = x*TMath::Cos(0.5*twist*TMath::DegToRad()) + y*TMath::Sin(0.5*twist*TMath::DegToRad()) + xc;
      fXY[i][1] = -x*TMath::Sin(0.5*twist*TMath::DegToRad()) + y*TMath::Cos(0.5*twist*TMath::DegToRad()) + yc;
   }
   ComputeTwist();
   if ((dz<0) || (h1<0) || (bl1<0) || (tl1<0) ||
       (h2<0) || (bl2<0) || (tl2<0)) SetShapeBit(kGeoRunTimeShape);
   else TGeoArb8::ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor providing the name of the shape.

TGeoGtra::TGeoGtra(const char *name, Double_t dz, Double_t theta, Double_t phi, Double_t twist, Double_t h1,
              Double_t bl1, Double_t tl1, Double_t alpha1, Double_t h2, Double_t bl2,
              Double_t tl2, Double_t alpha2)
         :TGeoTrap(name, dz, theta, phi, h1, bl1, tl1, alpha1, h2, bl2, tl2, alpha2)
{
   fTwistAngle = twist;
   Double_t x,y;
   Double_t th = theta*TMath::DegToRad();
   Double_t ph = phi*TMath::DegToRad();
   // Coordinates of the center of the bottom face
   Double_t xc = -dz*TMath::Sin(th)*TMath::Cos(ph);
   Double_t yc = -dz*TMath::Sin(th)*TMath::Sin(ph);

   Int_t i;

   for (i=0; i<4; i++) {
      x = fXY[i][0] - xc;
      y = fXY[i][1] - yc;
      fXY[i][0] = x*TMath::Cos(-0.5*twist*TMath::DegToRad()) + y*TMath::Sin(-0.5*twist*TMath::DegToRad()) + xc;
      fXY[i][1] = -x*TMath::Sin(-0.5*twist*TMath::DegToRad()) + y*TMath::Cos(-0.5*twist*TMath::DegToRad()) + yc;
   }
   // Coordinates of the center of the top face
   xc = -xc;
   yc = -yc;
   for (i=4; i<8; i++) {
      x = fXY[i][0] - xc;
      y = fXY[i][1] - yc;
      fXY[i][0] = x*TMath::Cos(0.5*twist*TMath::DegToRad()) + y*TMath::Sin(0.5*twist*TMath::DegToRad()) + xc;
      fXY[i][1] = -x*TMath::Sin(0.5*twist*TMath::DegToRad()) + y*TMath::Cos(0.5*twist*TMath::DegToRad()) + yc;
   }
   ComputeTwist();
   if ((dz<0) || (h1<0) || (bl1<0) || (tl1<0) ||
       (h2<0) || (bl2<0) || (tl2<0)) SetShapeBit(kGeoRunTimeShape);
   else TGeoArb8::ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGeoGtra::~TGeoGtra()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from inside point to surface of the shape.

Double_t TGeoGtra::DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
{
   if (iact<3 && safe) {
      // compute safe distance
      *safe = Safety(point, kTRUE);
      if (iact==0) return TGeoShape::Big();
      if (iact==1 && step<*safe) return TGeoShape::Big();
   }
   // compute distance to get outside this shape
   return TGeoArb8::DistFromInside(point, dir, iact, step, safe);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from inside point to surface of the shape.

Double_t TGeoGtra::DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
{
   if (iact<3 && safe) {
      // compute safe distance
      *safe = Safety(point, kTRUE);
      if (iact==0) return TGeoShape::Big();
      if (iact==1 && step<*safe) return TGeoShape::Big();
   }
   // compute distance to get outside this shape
   return TGeoArb8::DistFromOutside(point, dir, iact, step, safe);
}

////////////////////////////////////////////////////////////////////////////////
/// In case shape has some negative parameters, these has to be computed
/// in order to fit the mother

TGeoShape *TGeoGtra::GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix * /*mat*/) const
{
   if (!TestShapeBit(kGeoRunTimeShape)) return 0;
   if (mother->IsRunTimeShape()) {
      Error("GetMakeRuntimeShape", "invalid mother");
      return 0;
   }
   Double_t dz, h1, bl1, tl1, h2, bl2, tl2;
   if (fDz<0) dz=((TGeoTrap*)mother)->GetDz();
   else dz=fDz;
   if (fH1<0)
      h1 = ((TGeoTrap*)mother)->GetH1();
   else
      h1 = fH1;
   if (fH2<0)
      h2 = ((TGeoTrap*)mother)->GetH2();
   else
      h2 = fH2;
   if (fBl1<0)
      bl1 = ((TGeoTrap*)mother)->GetBl1();
   else
      bl1 = fBl1;
   if (fBl2<0)
      bl2 = ((TGeoTrap*)mother)->GetBl2();
   else
      bl2 = fBl2;
   if (fTl1<0)
      tl1 = ((TGeoTrap*)mother)->GetTl1();
   else
      tl1 = fTl1;
   if (fTl2<0)
      tl2 = ((TGeoTrap*)mother)->GetTl2();
   else
      tl2 = fTl2;
   return (new TGeoGtra(dz, fTheta, fPhi, fTwistAngle ,h1, bl1, tl1, fAlpha1, h2, bl2, tl2, fAlpha2));
}

////////////////////////////////////////////////////////////////////////////////
/// Computes the closest distance from given point to this shape.

Double_t TGeoGtra::Safety(const Double_t *point, Bool_t in) const
{
   return TGeoArb8::Safety(point,in);
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoGtra::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive)) return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   dz     = " << fDz << ";" << std::endl;
   out << "   theta  = " << fTheta << ";" << std::endl;
   out << "   phi    = " << fPhi << ";" << std::endl;
   out << "   twist  = " << fTwistAngle << ";" << std::endl;
   out << "   h1     = " << fH1<< ";" << std::endl;
   out << "   bl1    = " << fBl1<< ";" << std::endl;
   out << "   tl1    = " << fTl1<< ";" << std::endl;
   out << "   alpha1 = " << fAlpha1 << ";" << std::endl;
   out << "   h2     = " << fH2 << ";" << std::endl;
   out << "   bl2    = " << fBl2<< ";" << std::endl;
   out << "   tl2    = " << fTl2<< ";" << std::endl;
   out << "   alpha2 = " << fAlpha2 << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoGtra(\"" << GetName() << "\", dz,theta,phi,twist,h1,bl1,tl1,alpha1,h2,bl2,tl2,alpha2);" << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Set all arb8 params in one step.
///  - param[0] = dz
///  - param[1] = theta
///  - param[2] = phi
///  - param[3] = h1
///  - param[4] = bl1
///  - param[5] = tl1
///  - param[6] = alpha1
///  - param[7] = h2
///  - param[8] = bl2
///  - param[9] = tl2
///  - param[10] = alpha2
///  - param[11] = twist

void TGeoGtra::SetDimensions(Double_t *param)
{
   TGeoTrap::SetDimensions(param);
   fTwistAngle = param[11];
   Double_t x,y;
   Double_t twist = fTwistAngle;
   Double_t th = fTheta*TMath::DegToRad();
   Double_t ph = fPhi*TMath::DegToRad();
   // Coordinates of the center of the bottom face
   Double_t xc = -fDz*TMath::Sin(th)*TMath::Cos(ph);
   Double_t yc = -fDz*TMath::Sin(th)*TMath::Sin(ph);

   Int_t i;

   for (i=0; i<4; i++) {
      x = fXY[i][0] - xc;
      y = fXY[i][1] - yc;
      fXY[i][0] = x*TMath::Cos(-0.5*twist*TMath::DegToRad()) + y*TMath::Sin(-0.5*twist*TMath::DegToRad()) + xc;
      fXY[i][1] = -x*TMath::Sin(-0.5*twist*TMath::DegToRad()) + y*TMath::Cos(-0.5*twist*TMath::DegToRad()) + yc;
   }
   // Coordinates of the center of the top face
   xc = -xc;
   yc = -yc;
   for (i=4; i<8; i++) {
      x = fXY[i][0] - xc;
      y = fXY[i][1] - yc;
      fXY[i][0] = x*TMath::Cos(0.5*twist*TMath::DegToRad()) + y*TMath::Sin(0.5*twist*TMath::DegToRad()) + xc;
      fXY[i][1] = -x*TMath::Sin(0.5*twist*TMath::DegToRad()) + y*TMath::Cos(0.5*twist*TMath::DegToRad()) + yc;
   }
   ComputeTwist();
   if ((fDz<0) || (fH1<0) || (fBl1<0) || (fTl1<0) ||
       (fH2<0) || (fBl2<0) || (fTl2<0)) SetShapeBit(kGeoRunTimeShape);
   else TGeoArb8::ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoGtra::DistFromInside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromInside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from array of input points having directions specified by dirs. Store output in dists

void TGeoGtra::DistFromOutside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize, Double_t* step) const
{
   for (Int_t i=0; i<vecsize; i++) dists[i] = DistFromOutside(&points[3*i], &dirs[3*i], 3, step[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute safe distance from each of the points in the input array.
/// Input: Array of point coordinates, array of statuses for these points, size of the arrays
/// Output: Safety values

void TGeoGtra::Safety_v(const Double_t *points, const Bool_t *inside, Double_t *safe, Int_t vecsize) const
{
   for (Int_t i=0; i<vecsize; i++) safe[i] = Safety(&points[3*i], inside[i]);
}
