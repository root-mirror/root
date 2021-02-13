// @(#)root/g3d:$Id$
// Author: Rene Brun   13/06/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TSPHE.h"
#include "TNode.h"
#include "TBuffer.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TGeometry.h"
#include "TMath.h"

ClassImp(TSPHE);

/** \class TSPHE
\ingroup g3d
A Sphere.

It has 9 parameters:

  - name:       name of the shape
  - title:      shape's title
  - material:  (see TMaterial)
  - rmin:       minimum radius
  - rmax:       maximum radius
  - themin:     theta min
  - themax:     theta max
  - phimin:     phi min
  - phimax:     phi max

ROOT color indx = max(i-i0,j-j0);
*/

////////////////////////////////////////////////////////////////////////////////
/// SPHE shape default constructor

TSPHE::TSPHE()
{
   fRmin       = 0;
   fRmax       = 0;
   fThemin     = 0;
   fThemax     = 0;
   fPhimin     = 0;
   fPhimax     = 0;
   fSiTab      = 0;
   fCoTab      = 0;
   fCoThetaTab = 0;
   fNdiv       = 0;
   fAspectRatio=1.0;
   faX = faY = faZ = 1.0;      // Coeff along Ox
   fNz = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// SPHE shape normal constructor

TSPHE::TSPHE(const char *name, const char *title, const char *material, Float_t rmin, Float_t rmax, Float_t themin,
             Float_t themax, Float_t phimin, Float_t phimax)
     : TShape(name, title,material)
{
   fRmin   = rmin;
   fRmax   = rmax;
   fThemin = themin;
   fThemax = themax;
   fPhimin = phimin;
   fPhimax = phimax;

   fSiTab      = 0;
   fCoTab      = 0;
   fCoThetaTab = 0;
   fNdiv       = 0;

   fAspectRatio=1.0;
   faX = faY = faZ = 1.0;      // Coeff along Ox

   SetNumberOfDivisions (20);
}

////////////////////////////////////////////////////////////////////////////////
/// SPHE shape "simplified" constructor

TSPHE::TSPHE(const char *name, const char *title, const char *material, Float_t rmax)
     : TShape(name, title,material)
{
   fRmin   = 0;
   fRmax   = rmax;
   fThemin = 0;
   fThemax = 180;
   fPhimin = 0;
   fPhimax = 360;

   fSiTab      = 0;
   fCoTab      = 0;
   fCoThetaTab = 0;
   fNdiv       = 0;

   fAspectRatio=1.0;
   faX = faY = faZ = 1.0;      // Coeff along Ox

   SetNumberOfDivisions (20);
}

////////////////////////////////////////////////////////////////////////////////
/// SPHE shape default destructor

TSPHE::~TSPHE()
{
   if (fCoThetaTab) delete [] fCoThetaTab;
   if (fSiTab) delete [] fSiTab;
   if (fCoTab) delete [] fCoTab;

   fCoTab = 0;
   fSiTab = 0;
   fCoThetaTab=0;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from point px,py to a PSPHE
///
/// Compute the closest distance of approach from point px,py to each
/// computed outline point of the PSPHE (stolen from PCON).

Int_t TSPHE::DistancetoPrimitive(Int_t px, Int_t py)
{
   Int_t n = GetNumberOfDivisions()+1;
   Int_t numPoints = 2*n*(fNz+1);
   return ShapeDistancetoPrimitive(numPoints,px,py);
}

////////////////////////////////////////////////////////////////////////////////
/// Set ellipse.

void TSPHE::SetEllipse(const Float_t *factors)
{
   if (factors[0] > 0) faX = factors[0];
   if (factors[1] > 0) faY = factors[1];
   if (factors[2] > 0) faZ = factors[2];
}

////////////////////////////////////////////////////////////////////////////////
/// Set number of divisions.

void TSPHE::SetNumberOfDivisions (Int_t p)
{
   if (GetNumberOfDivisions () == p) return;
   fNdiv=p;
   fNz = Int_t(fAspectRatio*fNdiv*(fThemax - fThemin )/(fPhimax - fPhimin )) + 1;
   MakeTableOfCoSin();
}

////////////////////////////////////////////////////////////////////////////////
/// Create SPHE points

void TSPHE::SetPoints(Double_t *points) const
{
   Int_t i, j, n;
   Int_t indx = 0;

   n = GetNumberOfDivisions()+1;

   if (points) {
      if (!fCoTab)   MakeTableOfCoSin();
      Float_t z;
      for (i = 0; i < fNz+1; i++) {
         z = fRmin * fCoThetaTab[i]; // fSinPhiTab[i];
         Float_t sithet = TMath::Sqrt(TMath::Abs(1-fCoThetaTab[i]*fCoThetaTab[i]));
         Float_t zi = fRmin*sithet;
         for (j = 0; j < n; j++) {
            points[indx++] = faX*zi * fCoTab[j];
            points[indx++] = faY*zi * fSiTab[j];
            points[indx++] = faZ*z;
         }
         z = fRmax * fCoThetaTab[i];
         zi = fRmax*sithet;
         for (j = 0; j < n; j++) {
            points[indx++] = faX*zi * fCoTab[j];
            points[indx++] = faY*zi * fSiTab[j];
            points[indx++] = faZ*z;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return total X3D needed by TNode::ls (when called with option "x")

void TSPHE::Sizeof3D() const
{
   Int_t n;

   n = GetNumberOfDivisions()+1;
   Int_t nz = fNz+1;
   Bool_t specialCase = kFALSE;

   if (TMath::Abs(TMath::Sin(2*(fPhimax - fPhimin))) <= 0.01)  //mark this as a very special case, when
         specialCase = kTRUE;                                  //we have to draw this PCON like a TUBE

   gSize3D.numPoints += 2*n*nz;
   gSize3D.numSegs   += 4*(nz*n-1+(specialCase == kTRUE));
   gSize3D.numPolys  += 2*(nz*n-1+(specialCase == kTRUE));
}

////////////////////////////////////////////////////////////////////////////////
/// Make table of sine and cosine.

void TSPHE::MakeTableOfCoSin() const
{
   const Double_t pi  = TMath::ATan(1) * 4.0;
   const Double_t ragrad  = pi/180.0;

   Float_t dphi = fPhimax - fPhimin;
   while (dphi > 360) dphi -= 360;

   Float_t dtet = fThemax - fThemin;
   while (dtet > 180) dtet -= 180;

   Int_t j;
   Int_t n = GetNumberOfDivisions () + 1;
   if (fCoTab)
      delete [] fCoTab; // Delete the old tab if any
   fCoTab = new Double_t [n];
   if (!fCoTab ) return;

   if (fSiTab)
      delete [] fSiTab; // Delete the old tab if any
   fSiTab = new Double_t [n];
   if (!fSiTab ) return;

   Double_t range   = Double_t(dphi * ragrad);
   Double_t phi1    = Double_t(fPhimin  * ragrad);
   Double_t angstep = range/(n-1);

   Double_t ph = phi1;
   for (j = 0; j < n; j++)
   {
      ph = phi1 + j*angstep;
      fCoTab[j] = TMath::Cos(ph);
      fSiTab[j] = TMath::Sin(ph);
   }

   n  = fNz + 1;

   if (fCoThetaTab)
      delete [] fCoThetaTab; // Delete the old tab if any
   fCoThetaTab = new Double_t [n];
   if (!fCoThetaTab ) return;

   range   = Double_t(dtet * ragrad);
   phi1    = Double_t(fThemin  * ragrad);
   angstep = range/(n-1);

   ph = phi1;
   for (j = 0; j < n; j++)
   {
      fCoThetaTab[n-j-1] = TMath::Cos(ph);
      ph += angstep;
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Stream a class object

void TSPHE::Streamer(TBuffer &b)
{
   if (b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         b.ReadClassBuffer(TSPHE::Class(), this, R__v, R__s, R__c);
         Int_t ndiv = fNdiv;
         fNdiv = 0;
         SetNumberOfDivisions (ndiv);
         return;
      }
      //====process old versions before automatic schema evolution
      TShape::Streamer(b);
      b >> fRmin;    // minimum radius
      b >> fRmax;    // maximum radius
      b >> fThemin;  // minimum theta
      b >> fThemax;  // maximum theta
      b >> fPhimin;  // minimum phi
      b >> fPhimax;  // maximum phi
      Int_t tNdiv;   // XXX added by RvdE XXX (fNdiv is set by SetNumberOfDivisions)
      b >> tNdiv;
      if (R__v > 1) {
         b >> faX;
         b >> faY;
         b >> faZ;
      }
      SetNumberOfDivisions (tNdiv); // XXX added by RvdE
      b.CheckByteCount(R__s, R__c, TSPHE::IsA());
      //====end of old versions

   } else {
      b.WriteClassBuffer(TSPHE::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get buffer 3d.

const TBuffer3D & TSPHE::GetBuffer3D(Int_t reqSections) const
{
   static TBuffer3D buffer(TBuffer3DTypes::kGeneric);

   TShape::FillBuffer3D(buffer, reqSections);

   // Needed by kRawSizes / kRaw
   const Int_t n = GetNumberOfDivisions()+1;
   const Int_t nz = fNz+1;
   Bool_t specialCase = (TMath::Abs(TMath::Sin(2*(fPhimax - fPhimin))) <= 0.01);

   if (reqSections & TBuffer3D::kRawSizes) {
      Int_t nbPnts = 2*n*nz;
      Int_t nbSegs = 4*(nz*n-1+(specialCase == kTRUE));
      Int_t nbPols = 2*(nz*n-1+(specialCase == kTRUE));
      if (buffer.SetRawSizes(nbPnts, 3*nbPnts, nbSegs, 3*nbSegs, nbPols, 6*nbPols)) {
         buffer.SetSectionsValid(TBuffer3D::kRawSizes);
      }
   }
   if ((reqSections & TBuffer3D::kRaw) && buffer.SectionsValid(TBuffer3D::kRawSizes)) {
      // Points
      SetPoints(buffer.fPnts);
      if (!buffer.fLocalFrame) {
         TransformPoints(buffer.fPnts, buffer.NbPnts());
      }

      Int_t c = GetBasicColor();

      // Segments
      Int_t indx = 0;
      Int_t indx2 = 0;
      Int_t i, j, k;
      //inside & outside spheres, number of segments: 2*nz*(n-1)
      //             special case number of segments: 2*nz*n
      for (i = 0; i < nz*2; i++) {
         indx2 = i*n;
         for (j = 1; j < n; j++) {
            buffer.fSegs[indx++] = c;
            buffer.fSegs[indx++] = indx2+j-1;
            buffer.fSegs[indx++] = indx2+j;
         }
         if (specialCase) {
            buffer.fSegs[indx++] = c;
            buffer.fSegs[indx++] = indx2+j-1;
            buffer.fSegs[indx++] = indx2;
         }
      }

      //bottom & top lines, number of segments: 2*n
      for (i = 0; i < 2; i++) {
         indx2 = i*(nz-1)*2*n;
         for (j = 0; j < n; j++) {
            buffer.fSegs[indx++] = c;
            buffer.fSegs[indx++] = indx2+j;
            buffer.fSegs[indx++] = indx2+n+j;
         }
      }

      //inside & outside spheres, number of segments: 2*(nz-1)*n
      for (i = 0; i < (nz-1); i++) {

         //inside sphere
         indx2 = i*n*2;
         for (j = 0; j < n; j++) {
            buffer.fSegs[indx++] = c+2;
            buffer.fSegs[indx++] = indx2+j;
            buffer.fSegs[indx++] = indx2+n*2+j;
         }
         //outside sphere
         indx2 = i*n*2+n;
         for (j = 0; j < n; j++) {
            buffer.fSegs[indx++] = c+3;
            buffer.fSegs[indx++] = indx2+j;
            buffer.fSegs[indx++] = indx2+n*2+j;
         }
      }

      //left & right sections, number of segments: 2*(nz-2)
      //          special case number of segments: 0
      if (!specialCase) {
         for (i = 1; i < (nz-1); i++) {
            for (j = 0; j < 2; j++) {
               buffer.fSegs[indx++] = c;
               buffer.fSegs[indx++] =  2*i    * n + j*(n-1);
               buffer.fSegs[indx++] = (2*i+1) * n + j*(n-1);
            }
         }
      }

      // Polygons
      Int_t m = n - 1 + (specialCase == kTRUE);
      indx = 0;

      //bottom & top, number of polygons: 2*(n-1)
      // special case number of polygons: 2*n
      for (j = 0; j < n-1; j++) {
         buffer.fPols[indx++] = c+3;
         buffer.fPols[indx++] = 4;
         buffer.fPols[indx++] = 2*nz*m+j;
         buffer.fPols[indx++] = m+j;
         buffer.fPols[indx++] = 2*nz*m+j+1;
         buffer.fPols[indx++] = j;
      }
      for (j = 0; j < n-1; j++) {
         buffer.fPols[indx++] = c+3;
         buffer.fPols[indx++] = 4;
         buffer.fPols[indx++] = 2*nz*m+n+j;
         buffer.fPols[indx++] = (nz*2-2)*m+j;
         buffer.fPols[indx++] = 2*nz*m+n+j+1;
         buffer.fPols[indx++] = (nz*2-2)*m+m+j;
      }
      if (specialCase) {
         buffer.fPols[indx++] = c+3;
         buffer.fPols[indx++] = 4;
         buffer.fPols[indx++] = 2*nz*m+j;
         buffer.fPols[indx++] = m+j;
         buffer.fPols[indx++] = 2*nz*m;
         buffer.fPols[indx++] = j;

         buffer.fPols[indx++] = c+3;
         buffer.fPols[indx++] = 4;
         buffer.fPols[indx++] = 2*nz*m+n+j;
         buffer.fPols[indx++] = (nz*2-2)*m+j;
         buffer.fPols[indx++] = 2*nz*m+n;
         buffer.fPols[indx++] = (nz*2-2)*m+m+j;
      }

      //inside & outside, number of polygons: (nz-1)*2*(n-1)
      for (k = 0; k < (nz-1); k++) {
         for (j = 0; j < n-1; j++) {
            buffer.fPols[indx++] = c;
            buffer.fPols[indx++] = 4;
            buffer.fPols[indx++] = 2*k*m+j;
            buffer.fPols[indx++] = nz*2*m+(2*k+2)*n+j+1;
            buffer.fPols[indx++] = (2*k+2)*m+j;
            buffer.fPols[indx++] = nz*2*m+(2*k+2)*n+j;
         }
         for (j = 0; j < n-1; j++) {
            buffer.fPols[indx++] = c+1;
            buffer.fPols[indx++] = 4;
            buffer.fPols[indx++] = (2*k+1)*m+j;
            buffer.fPols[indx++] = nz*2*m+(2*k + 3)*n+j;
            buffer.fPols[indx++] = (2*k+ 3)*m+j;
            buffer.fPols[indx++] = nz*2*m+(2*k+3)*n+j+1;
         }

         if (specialCase) {
            buffer.fPols[indx++] = c;
            buffer.fPols[indx++] = 4;
            buffer.fPols[indx++] = 2*k*m+j;
            buffer.fPols[indx++] = nz*2*m+(2*k+2)*n+j;
            buffer.fPols[indx++] = (2*k+2)*m+j;
            buffer.fPols[indx++] = nz*2*m+(2*k+2)*n;

            buffer.fPols[indx++] = c+1;
            buffer.fPols[indx++] = 4;
            buffer.fPols[indx++] = (2*k+1)*m+j;
            buffer.fPols[indx++] = nz*2*m+(2*k+3)*n+j;
            buffer.fPols[indx++] = (2*k+3)*m+j;
            buffer.fPols[indx++] = nz*2*m+(2*k+3)*n;
         }
      }

      //left & right sections, number of polygons: 2*(nz-1)
      //          special case number of polygons: 0
      if (!specialCase) {
         indx2 = nz*2*(n-1);
         for (k = 0; k < (nz-1); k++) {
            buffer.fPols[indx++] = c+2;
            buffer.fPols[indx++] = 4;
            buffer.fPols[indx++] = k==0 ? indx2 : indx2+2*nz*n+2*(k-1);
            buffer.fPols[indx++] = indx2+2*(k+1)*n;
            buffer.fPols[indx++] = indx2+2*nz*n+2*k;
            buffer.fPols[indx++] = indx2+(2*k+3)*n;

            buffer.fPols[indx++] = c+2;
            buffer.fPols[indx++] = 4;
            buffer.fPols[indx++] = k==0 ? indx2+n-1 : indx2+2*nz*n+2*(k-1)+1;
            buffer.fPols[indx++] = indx2+(2*k+3)*n+n-1;
            buffer.fPols[indx++] = indx2+2*nz*n+2*k+1;
            buffer.fPols[indx++] = indx2+2*(k+1)*n+n-1;
         }

         buffer.fPols[indx-8] = indx2+n;
         buffer.fPols[indx-2] = indx2+2*n-1;
      }

      buffer.SetSectionsValid(TBuffer3D::kRaw);
   }

   return buffer;
}
