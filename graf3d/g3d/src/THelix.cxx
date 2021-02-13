// @(#)root/g3d:$Id$
// Author: Ping Yeh   19/12/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class THelix
\ingroup g3d
THelix has two different constructors.

If a particle with charge q passes through a point (x,y,z)
with momentum (px,py,pz) with magnetic field B along an axis (nx,ny,nz),
this helix can be constructed like:

~~~ {.cpp}
      THelix p(x,y,z, px,py,pz, q*B, nx,ny,nz);
      (nx,ny,nz) defaults to (0,0,1).
~~~

A helix in its own frame can be defined with a pivotal point
(x0,y0,z0), the velocity at that point (vx0,vy0,vz0), and
an angular frequency w.  Combining vx0 and vy0 to a transverse
velocity vt0 one can parametrize the helix as:

~~~ {.cpp}
    x(t) = x0 - vt0 / w * sin(-w * t + phi0)
    y(t) = y0 + vt0 / w * cos(-w * t + phi0)
    z(t) = z0 + vz0 * t
~~~

The second constructor has 6 parameters,

Example:

~~~ {.cpp}
      THelix pl1(xyz, v, w, range, rtype, axis);
~~~

where:

  - xyz  : array of initial position
  - v    : array of initial velocity
  - w    : angular frequency
  - range: helix range
  - rtype: kHelixZ specifies allowed drawing range in helix Z direction, i.e., along B field.
           kLabZ specifies drawing range in lab frame.
           kHelixX, kHelixY, kLabX, kLabY, kUnchanged ... etc can also be specified
  - axis : helix axis

Example constructing a helix with several default values and drawing it:

Begin_Macro(source)
{
   TCanvas* helix_example_c1 = new TCanvas("helix_example_c1");
   TView *view = TView::CreateView(1);
   view->SetRange(-1,-1,-1,1,1,1);
   THelix *helix = new THelix(0., 0., 0., 1., 0., 0.3, 10.);
   helix->Draw();
}
End_Macro

This initializes a helix with its axis in Z direction (rtype=kHelixZ).
*/

#include <iostream>
#include "TBuffer.h"
#include "TROOT.h"
#include "THelix.h"
#include "TMath.h"

Int_t THelix::fgMinNSeg=5;        // at least 5 line segments in TPolyLine3D

ClassImp(THelix);

////////////////////////////////////////////////////////////////////////////////
/// Set all helix parameters.

void  THelix::SetHelix(Double_t const* xyz,  Double_t const* v,  Double_t w,
                       Double_t const* range, EHelixRangeType rType,
                       Double_t const* axis )
{
   // Define the helix frame by setting the helix axis and rotation matrix
   SetAxis(axis);

   // Calculate initial position and velocity in helix frame
   fW    = w;
   Double_t * m = fRotMat->GetMatrix();
   Double_t vx0, vy0, vz0;
   vx0   = v[0] * m[0] + v[1] * m[1] + v[2] * m[2];
   vy0   = v[0] * m[3] + v[1] * m[4] + v[2] * m[5];
   vz0   = v[0] * m[6] + v[1] * m[7] + v[2] * m[8];
   fVt   = TMath::Sqrt(vx0*vx0 + vy0*vy0);
   fPhi0 = TMath::ATan2(vy0,vx0);
   fVz   = vz0;
   fX0   = xyz[0] * m[0] +  xyz[1] * m[1] +  xyz[2] * m[2];
   fY0   = xyz[0] * m[3] +  xyz[1] * m[4] +  xyz[2] * m[5];
   fZ0   = xyz[0] * m[6] +  xyz[1] * m[7] +  xyz[2] * m[8];
   if (fW != 0) {
      fX0 += fVt / fW * TMath::Sin(fPhi0);
      fY0 -= fVt / fW * TMath::Cos(fPhi0);
   }

   // Then calculate the range in t and set the polyline representation
   Double_t r1 = 0;
   Double_t r2 = 1;
   if (range) {r1 = range[0]; r2 = range[1];}
   if (rType != kUnchanged) {
      fRange[0] = 0.0;   fRange[1] = TMath::Pi();   // initialize to half round
      SetRange(r1,r2,rType);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Helix default constructor.

THelix::THelix()
{
   fX0 = fY0 = fZ0 = fVt = fPhi0 = fVz = fAxis[0] = fAxis[1] = 0.0;
   fAxis[2]  = 1.0;
   fW        = 1.5E7;   // roughly the cyclon frequency of proton in AMS
   fRange[0] = 0.0;
   fRange[1] = 1.0;
   fRotMat   = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Helix normal constructor.

THelix::THelix(Double_t x,  Double_t y,  Double_t z,
               Double_t vx, Double_t vy, Double_t vz,
               Double_t w)
        : TPolyLine3D()
{
   Double_t p[3], v[3];
   p[0] = x;
   p[1] = y;
   p[2] = z;
   v[0] = vx;
   v[1] = vy;
   v[2] = vz;
   Double_t *range = 0;
   fRotMat   = 0;

   SetHelix(p, v, w, range, kHelixZ);
   fOption = "";
}

////////////////////////////////////////////////////////////////////////////////
/// Helix normal constructor.

THelix::THelix(Double_t const* xyz, Double_t const* v, Double_t w,
               Double_t const* range, EHelixRangeType rType, Double_t const* axis)
        : TPolyLine3D()
{
   Double_t r[2];
   if ( range ) {
      r[0] = range[0];   r[1] = range[1];
   } else {
      r[0] = 0.0;        r[1] = 1.0;
   }

   fRotMat   = 0;
   if ( axis ) {                        // specify axis
      SetHelix(xyz, v, w, r, rType, axis);
   } else {                             // default axis
      SetHelix(xyz, v, w, r, rType);
   }

   fOption = "";
}


#if 0
////////////////////////////////////////////////////////////////////////////////
/// Helix copy constructor.

THelix::THelix(const THelix &h) : TPolyLine3D()
{
   fX0   = h.fX0;
   fY0   = h.fY0;
   fZ0   = h.fZ0;
   fVt   = h.fVt;
   fPhi0 = h.fPhi0;
   fVz   = h.fVz;
   fW    = h.fW;
   for (Int_t i=0; i<3; i++) fAxis[i] = h.fAxis[i];
   fRotMat = new TRotMatrix(*(h.fRotMat));
   fRange[0] = h.fRange[0];
   fRange[1] = h.fRange[1];

   fOption = h.fOption;
}
#endif

////////////////////////////////////////////////////////////////////////////////
/// assignment operator

THelix& THelix::operator=(const THelix& hx)
{
   if(this!=&hx) {
      TPolyLine3D::operator=(hx);
      fX0=hx.fX0;
      fY0=hx.fY0;
      fZ0=hx.fZ0;
      fVt=hx.fVt;
      fPhi0=hx.fPhi0;
      fVz=hx.fVz;
      fW=hx.fW;
      for(Int_t i=0; i<3; i++)
         fAxis[i]=hx.fAxis[i];
      fRotMat=hx.fRotMat;
      for(Int_t i=0; i<2; i++)
         fRange[i]=hx.fRange[i];
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Helix destructor.

THelix::~THelix()
{
   if (fRotMat) delete fRotMat;
}


////////////////////////////////////////////////////////////////////////////////
/// Helix copy constructor.

THelix::THelix(const THelix &helix) : TPolyLine3D(helix)
{
   fRotMat=0;
   ((THelix&)helix).THelix::Copy(*this);
}


////////////////////////////////////////////////////////////////////////////////
/// Copy this helix to obj.

void THelix::Copy(TObject &obj) const
{
   TObject::Copy(obj);
   TAttLine::Copy(((THelix&)obj));

   ((THelix&)obj).fX0        = fX0;
   ((THelix&)obj).fY0        = fY0;
   ((THelix&)obj).fZ0        = fZ0;
   ((THelix&)obj).fVt        = fVt;
   ((THelix&)obj).fPhi0      = fPhi0;
   ((THelix&)obj).fVz        = fVz;
   ((THelix&)obj).fW         = fW;
   for (Int_t i=0; i<3; i++)
      ((THelix&)obj).fAxis[i] = fAxis[i];

   if (((THelix&)obj).fRotMat)
      delete ((THelix&)obj).fRotMat;
   ((THelix&)obj).fRotMat    = new TRotMatrix(*fRotMat);

   ((THelix&)obj).fRange[0]  = fRange[0];
   ((THelix&)obj).fRange[1]  = fRange[1];

   ((THelix&)obj).fOption    = fOption;

   //
   // Set range and make the graphic representation
   //
   ((THelix&)obj).SetRange(fRange[0], fRange[1], kHelixT);
}


////////////////////////////////////////////////////////////////////////////////
/// Draw this helix with its current attributes.

void THelix::Draw(Option_t *option)
{
   AppendPad(option);
}


////////////////////////////////////////////////////////////////////////////////
/// Dump this helix with its attributes.

void THelix::Print(Option_t *option) const
{
   std::cout <<"    THelix Printing N=" <<fN<<" Option="<<option<<std::endl;
}


////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out.

void THelix::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   char quote = '"';
   out<<"   "<<std::endl;
   if (gROOT->ClassSaved(THelix::Class())) {
      out<<"   ";
   } else {
      out<<"   THelix *";
   }
   out<<"helix = new THelix("<<fX0<<","<<fY0<<","<<fZ0<<","
      <<fVt*TMath::Cos(fPhi0)<<","<<fVt*TMath::Sin(fPhi0)<<","<<fVz<<","
      <<fW<<","<<fRange[0]<<","<<fRange[1]<<","<<(Int_t)kHelixT<<","
      <<fAxis[0]<<","<<fAxis[1]<<","<<fAxis[2]<<","
      <<quote<<fOption<<quote<<");"<<std::endl;

   SaveLineAttributes(out,"helix",1,1,1);

   out<<"   helix->Draw();"<<std::endl;
}


////////////////////////////////////////////////////////////////////////////////
/// Set a new axis for the helix.  This will make a new rotation matrix.

void THelix::SetAxis(Double_t const* axis)
{
   if (axis) {
      Double_t len = TMath::Sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
      if (len <= 0) {
         Error("SetAxis()", "Impossible! axis length %lf <= 0!", len);
         return;
      }
      fAxis[0] = axis[0]/len;
      fAxis[1] = axis[1]/len;
      fAxis[2] = axis[2]/len;
   } else {
      fAxis[0] = 0;
      fAxis[1] = 0;
      fAxis[2] = 1;
   }

   // Construct the rotational matrix from the axis
   SetRotMatrix();
}


////////////////////////////////////////////////////////////////////////////////
/// Set axis.

void THelix::SetAxis(Double_t x, Double_t y, Double_t z)
{
   Double_t axis[3];    axis[0] = x;    axis[1] = y;    axis[2] = z;
   SetAxis(axis);
}


////////////////////////////////////////////////////////////////////////////////
/// Set a new range for the helix.  This will remake the polyline.

void THelix::SetRange(Double_t * range, EHelixRangeType rType)
{
   Double_t a[2];
   Double_t halfpi = TMath::Pi()/2.0;
   Int_t i;
   Double_t vx = fVt * TMath::Cos(fPhi0);
   Double_t vy = fVt * TMath::Sin(fPhi0);
   Double_t phase;

   if ( fW != 0 && fVz != 0 ) {         // general case
      switch ( rType ) {
         case kHelixT :
            fRange[0] = range[0];  fRange[1] = range[1];  break;

         case kHelixX :
            for (i=0; i<2; i++ ) {
               a[i] = fW / fVt * (range[i] - fX0);
               if ( a[i] < -1 || a[i] > 1 ) {
                  Error("SetRange()",
                        "range out of bound (%lf:%lf): %lf.  Default used: %lf",
                        fX0-fVt/fW, fX0+fVt/fW, range[i], fRange[i]);
                  return;
               }
               phase = FindClosestPhase(fPhi0+halfpi, a[i]);
               fRange[i] = ( fPhi0 + halfpi - phase ) / fW;
            }
            break;

         case kHelixY :
            for (i=0; i<2; i++ ) {
               a[i] = fW / fVt * (range[i] - fY0);
               if ( a[i] < -1 || a[i] > 1 ) {
                  Error("SetRange()",
                        "range out of bound (%lf:%lf): %lf.  Default used: %lf",
                         fY0-fVt/fW, fY0+fVt/fW, range[i], fRange[i]);
                  return;
               }
               phase = FindClosestPhase(fPhi0, a[i]);
               fRange[i] = ( fPhi0 - phase ) / fW;
            }
            break;

         case kHelixZ :
            if ( fVz != 0 ) {
               for (i=0; i<2; i++ ) {
                  fRange[i] = (range[i] - fZ0) / fVz;
               }
            } else {                // fVz = 0, z = constant = fZ0
               Error("SetRange()",
                     "Vz = 0 and attempts to set range along helix axis!");
               return;
            }
            break;

         case kLabX :
         case kLabY :
         case kLabZ :
            printf("setting range in lab axes is not implemented yet\n");
            break;
         default:
            Error("SetRange()","unknown range type %d", rType);
            break;
      }
   } else if ( fW == 0 ) {                // straight line: x = x0 + vx * t
      switch ( rType ) {
         case kHelixT :
            fRange[0] = range[0];  fRange[1] = range[1];
            break;
         case kHelixX :
            if ( vx != 0 ) {
               fRange[0] = (range[0] - fX0) / vx;
               fRange[1] = (range[1] - fX0) / vx;
            } else {
               Error("SetRange()",
                     "Vx = 0 and attempts to set range on helix x axis!");
               return;
            }
            break;
         case kHelixY :
            if ( vy != 0 ) {
               fRange[0] = (range[0] - fY0) / vy;
               fRange[1] = (range[1] - fY0) / vy;
            } else {
               Error("SetRange()",
                     "Vy = 0 and attempts to set range on helix y axis!");
               return;
            }
            break;
         case kHelixZ :
            if ( fVz != 0 ) {
               fRange[0] = (range[0] - fZ0) / fVz;
               fRange[1] = (range[1] - fZ0) / fVz;
            } else {
               Error("SetRange()",
                     "Vz = 0 and attempts to set range on helix z axis!");
               return;
            }
            break;
         case kLabX   :
         case kLabY   :
         case kLabZ   :
            printf("setting range in lab axes is not implemented yet\n");
            break;
         default      :
            Error("SetRange()","unknown range type %d", rType);
            break;
      }
   } else if ( fVz == 0 ) {               // a circle, not fully implemented yet
      switch ( rType ) {
         case kHelixT :
            fRange[0] = range[0];  fRange[1] = range[1];  break;
         case kHelixX :
            if ( vx != 0 ) {
               fRange[0] = (range[0] - fX0) / vx;
               fRange[1] = (range[1] - fX0) / vx;
            } else {
               Error("SetRange()",
                     "Vx = 0 and attempts to set range on helix x axis!");
               return;
            }
            break;
         case kHelixY :
            if ( vy != 0 ) {
               fRange[0] = (range[0] - fY0) / vy;
               fRange[1] = (range[1] - fY0) / vy;
            } else {
               Error("SetRange()",
                     "Vy = 0 and attempts to set range on helix y axis!");
               return;
            }
            break;
         case kHelixZ :
            Error("SetRange()",
                  "Vz = 0 and attempts to set range on helix z axis!");
            return;
         case kLabX   :
         case kLabY   :
         case kLabZ   :
            printf("setting range in lab axes is not implemented yet\n");
            break;
         default      :
            Error("SetRange()","unknown range type %d", rType);
            break;
      }
   }

   if ( fRange[0] > fRange[1] ) {
      Double_t temp = fRange[1];   fRange[1] = fRange[0];  fRange[0] = temp;
   }

   // Set the polylines in global coordinates
   Double_t degrad  = TMath::Pi() / 180.0;
   Double_t segment = 5.0 * degrad;             // 5 degree segments
   Double_t dt      = segment / TMath::Abs(fW); // parameter span on each segm.

   Int_t    nSeg    = Int_t((fRange[1]-fRange[0]) / dt) + 1;
   if (nSeg < THelix::fgMinNSeg) {
      nSeg = THelix::fgMinNSeg;
      dt = (fRange[1]-fRange[0])/nSeg;
   }

   Double_t * xl    = new Double_t[nSeg+1];     // polyline in local coordinates
   Double_t * yl    = new Double_t[nSeg+1];
   Double_t * zl    = new Double_t[nSeg+1];

   for (i=0; i<=nSeg; i++) {                    // calculate xl[], yl[], zl[];
      Double_t t, phase2;
      if (i==nSeg) t = fRange[1];                // the last point
      else         t = fRange[0] + dt * i;
      phase2 = -fW * t + fPhi0;
      xl[i] = fX0 - fVt/fW * TMath::Sin(phase2);
      yl[i] = fY0 + fVt/fW * TMath::Cos(phase2);
      zl[i] = fZ0 + fVz * t;
   }

   Float_t xg, yg,zg;     // global coordinates
                          // must be Float_t to call TPolyLine3D::SetPoint()
   Double_t * m = fRotMat->GetMatrix();
   TPolyLine3D::SetPolyLine(nSeg+1);
   for (i=0; i<=nSeg; i++) {                    // m^{-1} = transpose of m
      xg =  xl[i] * m[0]  +  yl[i] * m[3]  +  zl[i] * m[6] ;
      yg =  xl[i] * m[1]  +  yl[i] * m[4]  +  zl[i] * m[7] ;
      zg =  xl[i] * m[2]  +  yl[i] * m[5]  +  zl[i] * m[8] ;
      TPolyLine3D::SetPoint(i,xg,yg,zg);
   }

   delete[] xl;  delete[] yl;    delete[] zl;
}


////////////////////////////////////////////////////////////////////////////////
/// Set range.

void THelix::SetRange(Double_t r1, Double_t r2, EHelixRangeType rType)
{
   Double_t range[2];
   range[0] = r1;       range[1] = r2;
   SetRange(range, rType);
}


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                   Protected     Member     Functions                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
/// Set the rotational matrix according to the helix axis.

void THelix::SetRotMatrix()
{
   // Calculate all 6 angles.
   // Note that TRotMatrix::TRotMatrix() expects angles in degrees.
   Double_t raddeg = 180.0 / TMath::Pi();
   Double_t halfpi = TMath::Pi()/2.0 * raddeg;
                                 // (theta3,phi3) is the helix axis
   Double_t theta3 = TMath::ACos(fAxis[2]) * raddeg;
   Double_t phi3   = TMath::ATan2(fAxis[1], fAxis[0]) * raddeg;
                                 //  (theta1,phi1) is the x-axis in helix frame
   Double_t theta1 = theta3 + halfpi;
   Double_t phi1   = phi3;
                                 //  (theta2,phi2) is the y-axis in helix frame
   Double_t theta2 = halfpi;
   Double_t phi2   = phi1 + halfpi;

   // Delete the old rotation matrix
   if (fRotMat) delete fRotMat;

   // Make a new rotation matrix
   fRotMat = new TRotMatrix("HelixRotMat", "Master frame -> Helix frame",
                            theta1, phi1,  theta2, phi2,  theta3, phi3 );
   return;
}


////////////////////////////////////////////////////////////////////////////////
/// Finds the closest phase to phi0 that gives cos(phase) = cosine

Double_t  THelix::FindClosestPhase(Double_t phi0,  Double_t cosine)
{
   const Double_t pi    = TMath::Pi();
   const Double_t twopi = TMath::Pi() * 2.0;
   Double_t phi1 = TMath::ACos(cosine);
   Double_t phi2 = - phi1;

   while ( phi1 - phi0 >  pi )   phi1 -= twopi;
   while ( phi1 - phi0 < -pi )   phi1 += twopi;

   while ( phi2 - phi0 >  pi )   phi2 -= twopi;
   while ( phi2 - phi0 < -pi )   phi2 += twopi;

   // Now phi1, phi2 and phi0 are within the same 2pi range
   // and cos(phi1) = cos(phi2) = cosine
   if ( TMath::Abs(phi1-phi0) < TMath::Abs(phi2-phi0) )  return phi1;
   else                                                  return phi2;
}


////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class THelix.

void THelix::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      if (R__v > 1) {
         R__b.ReadClassBuffer(THelix::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TPolyLine3D::Streamer(R__b);
      R__b >> fX0;
      R__b >> fY0;
      R__b >> fZ0;
      R__b >> fVt;
      R__b >> fPhi0;
      R__b >> fVz;
      R__b >> fW;
      R__b.ReadStaticArray(fAxis);
      R__b >> fRotMat;
      R__b.ReadStaticArray(fRange);
      R__b.CheckByteCount(R__s, R__c, THelix::IsA());
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(THelix::Class(),this);
   }
}
