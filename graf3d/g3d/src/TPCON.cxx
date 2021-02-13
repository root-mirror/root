// @(#)root/g3d:$Id$
// Author: Nenad Buncic   29/09/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TPCON.h"
#include "TNode.h"
#include "TMath.h"
#include "TBuffer.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TGeometry.h"

ClassImp(TPCON);

/** \class TPCON
\ingroup g3d
A polycone

\image html g3d_pcon.png
It has the following parameters:

  - name:       name of the shape
  - title:      shape's title
  - material:   (see TMaterial)
  - phi1:       the azimuthal angle phi at which the volume begins (angles
                are counted counterclockwise)
  - dphi:       opening angle of the volume, which extends from
                phi1 to phi1+dphi
  - nz:         number of planes perpendicular to the z axis where
                the dimension of the section is given -- this number
                should be at least 2
  - rmin:       array of dimension nz with minimum radius at a given plane
  - rmax:       array of dimension nz with maximum radius at a given plane
  - z:          array of dimension nz with z position of given plane
*/

////////////////////////////////////////////////////////////////////////////////
/// PCON shape default constructor

TPCON::TPCON()
{
   fRmin  = 0;
   fRmax  = 0;
   fDz    = 0;
   fCoTab = 0;
   fSiTab = 0;
   fPhi1  = 0.;
   fDphi1 = 0.;
   fNz    = 0;
   fNdiv  = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// PCON shape normal constructor
///
/// Parameters of the nz positions must be entered via TPCON::DefineSection.

TPCON::TPCON(const char *name, const char *title, const char *material, Float_t phi1, Float_t dphi1, Int_t nz)
      : TShape(name, title,material)
{
   if (nz < 2 ) {
      Error(name, "number of z planes for %s must be at least two !", name);
      return;
   }
   fPhi1  = phi1;
   fDphi1 = dphi1;
   fNz    = nz;
   fNdiv  = 0;
   fRmin  = new Float_t [nz+1];
   fRmax  = new Float_t [nz+1];
   fDz    = new Float_t [nz+1];

   fCoTab = 0;
   fSiTab = 0;

   while (fDphi1 > 360) fDphi1 -= 360;

   MakeTableOfCoSin();
}

////////////////////////////////////////////////////////////////////////////////
/// copy constructor

TPCON::TPCON(const TPCON& pc) :
  TShape(pc),
  fSiTab(pc.fSiTab),
  fCoTab(pc.fCoTab),
  fPhi1(pc.fPhi1),
  fDphi1(pc.fDphi1),
  fNdiv(pc.fNdiv),
  fNz(pc.fNz),
  fRmin(pc.fRmin),
  fRmax(pc.fRmax),
  fDz(pc.fDz)
{
}

////////////////////////////////////////////////////////////////////////////////
/// assignment operator

TPCON& TPCON::operator=(const TPCON& pc)
{
   if(this!=&pc) {
      TShape::operator=(pc);
      fSiTab=pc.fSiTab;
      fCoTab=pc.fCoTab;
      fPhi1=pc.fPhi1;
      fDphi1=pc.fDphi1;
      fNdiv=pc.fNdiv;
      fNz=pc.fNz;
      fRmin=pc.fRmin;
      fRmax=pc.fRmax;
      fDz=pc.fDz;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Make table of cosine and sine

void TPCON::MakeTableOfCoSin() const
{
   const Double_t pi  = TMath::ATan(1) * 4.0;
   const Double_t ragrad  = pi/180.0;

   Int_t n = GetNumberOfDivisions () + 1;
   if (fCoTab) delete [] fCoTab; // Delete the old tab if any
   fCoTab = new Double_t [n];
   if (!fCoTab ) return;

   if (fSiTab) delete [] fSiTab; // Delete the old tab if any
   fSiTab = new Double_t [n];
   if (!fSiTab ) return;

   Double_t range   = Double_t(fDphi1 * ragrad);
   Double_t phi1    = Double_t(fPhi1  * ragrad);
   Double_t angstep = range/(n-1);

   FillTableOfCoSin(phi1,angstep,n);
}

////////////////////////////////////////////////////////////////////////////////
/// PCON shape default destructor

TPCON::~TPCON()
{
   if (fRmin) delete [] fRmin;
   if (fRmax) delete [] fRmax;
   if (fDz)   delete [] fDz;
   if (fSiTab) delete [] fSiTab;
   if (fCoTab) delete [] fCoTab;

   fRmin = 0;
   fRmax = 0;
   fDz   = 0;
   fCoTab = 0;
   fSiTab = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Defines section secNum of the polycone
///
///  - rmin  radius of the inner circle in the cross-section
///  - rmax  radius of the outer circle in the cross-section
///  - z     z coordinate of the section

void TPCON::DefineSection(Int_t secNum, Float_t z, Float_t rmin, Float_t rmax)
{
   if ((secNum < 0) || (secNum >= fNz)) return;

   fRmin[secNum] = rmin;
   fRmax[secNum] = rmax;
   fDz[secNum]   = z;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from point px,py to a PCON
///
/// Compute the closest distance of approach from point px,py to each
/// computed outline point of the PCON.

Int_t TPCON::DistancetoPrimitive(Int_t px, Int_t py)
{
   Int_t n = GetNumberOfDivisions()+1;
   Int_t numPoints = fNz*2*n;
   return ShapeDistancetoPrimitive(numPoints,px,py);
}

////////////////////////////////////////////////////////////////////////////////
/// Fill the table of cos and sin to prepare drawing

void  TPCON::FillTableOfCoSin(Double_t phi, Double_t angstep,Int_t n) const
{
   Double_t ph = phi-angstep;
   for (Int_t j = 0; j < n; j++) {
      ph += angstep;
      fCoTab[j] = TMath::Cos(ph);
      fSiTab[j] = TMath::Sin(ph);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set number of divisions.

void TPCON::SetNumberOfDivisions (Int_t p)
{
   if (GetNumberOfDivisions () == p) return;
   fNdiv=p;
   MakeTableOfCoSin();
}

////////////////////////////////////////////////////////////////////////////////
/// Create PCON points

void TPCON::SetPoints(Double_t *points) const
{
   Int_t i, j;
   Int_t indx = 0;

   Int_t n = GetNumberOfDivisions()+1;

   if (points) {
      if (!fCoTab) MakeTableOfCoSin();
      for (i = 0; i < fNz; i++) {
         for (j = 0; j < n; j++) {
            points[indx++] = fRmin[i] * fCoTab[j];
            points[indx++] = fRmin[i] * fSiTab[j];
            points[indx++] = fDz[i];
         }
         for (j = 0; j < n; j++) {
            points[indx++] = fRmax[i] * fCoTab[j];
            points[indx++] = fRmax[i] * fSiTab[j];
            points[indx++] = fDz[i];
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return total X3D needed by TNode::ls (when called with option "x")

void TPCON::Sizeof3D() const
{
   Int_t n;

   n = GetNumberOfDivisions()+1;

   gSize3D.numPoints += fNz*2*n;
   gSize3D.numSegs   += 4*(fNz*n-1+(fDphi1 == 360));
   gSize3D.numPolys  += 2*(fNz*n-1+(fDphi1 == 360));
}

////////////////////////////////////////////////////////////////////////////////
/// Stream a class object

void TPCON::Streamer(TBuffer &b)
{
   if (b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         b.ReadClassBuffer(TPCON::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TShape::Streamer(b);
      b >> fPhi1;
      b >> fDphi1;
      b >> fNz;
      fRmin  = new Float_t [fNz];
      fRmax  = new Float_t [fNz];
      fDz    = new Float_t [fNz];
      b.ReadArray(fRmin);
      b.ReadArray(fRmax);
      b.ReadArray(fDz);
      b >> fNdiv;
      b.CheckByteCount(R__s, R__c, TPCON::IsA());
      //====end of old versions

   } else {
      b.WriteClassBuffer(TPCON::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get buffer 3d.

const TBuffer3D & TPCON::GetBuffer3D(Int_t reqSections) const
{
   static TBuffer3D buffer(TBuffer3DTypes::kGeneric);

   TShape::FillBuffer3D(buffer, reqSections);

   // No kShapeSpecific or kBoundingBox

   if (reqSections & TBuffer3D::kRawSizes)
   {
      const Int_t n = GetNumberOfDivisions()+1;
      Int_t nbPnts = fNz*2*n;
      Bool_t specialCase = (fDphi1 == 360);
      Int_t nbSegs = 4*(fNz*n-1+(specialCase == kTRUE));
      Int_t nbPols = 2*(fNz*n-1+(specialCase == kTRUE));

      if (buffer.SetRawSizes(nbPnts, 3*nbPnts, nbSegs, 3*nbSegs, nbPols, 6*nbPols)) {
         buffer.SetSectionsValid(TBuffer3D::kRawSizes);
      }
   }
   if ((reqSections & TBuffer3D::kRaw) && buffer.SectionsValid(TBuffer3D::kRawSizes))
   {
      // Points
      SetPoints(buffer.fPnts);
      if (!buffer.fLocalFrame) {
         TransformPoints(buffer.fPnts, buffer.NbPnts());
      }

      // Segments and Polygons
      if (SetSegsAndPols(buffer))
      {
         buffer.SetSectionsValid(TBuffer3D::kRaw);
      }
   }
   return buffer;
}

////////////////////////////////////////////////////////////////////////////////
/// Set segments and polygons.

Bool_t TPCON::SetSegsAndPols(TBuffer3D & buffer) const
{
   if (fNz < 2) return kFALSE;
   const Int_t n = GetNumberOfDivisions()+1;
   Bool_t specialCase = (fDphi1 == 360);

   Int_t c = GetBasicColor();

   Int_t i, j, k;
   Int_t indx = 0;
   Int_t indx2 = 0;

   //inside & outside circles, number of segments: 2*fNz*(n-1)
   //             special case number of segments: 2*fNz*n
   for (i = 0; i < fNz*2; i++) {
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
      indx2 = i*(fNz-1)*2*n;
      for (j = 0; j < n; j++) {
         buffer.fSegs[indx++] = c;
         buffer.fSegs[indx++] = indx2+j;
         buffer.fSegs[indx++] = indx2+n+j;
      }
   }

   //inside & outside cilindres, number of segments: 2*(fNz-1)*n
   for (i = 0; i < (fNz-1); i++) {

      //inside cilinder
      indx2 = i*n*2;
      for (j = 0; j < n; j++) {
         buffer.fSegs[indx++] = c+2;
         buffer.fSegs[indx++] = indx2+j;
         buffer.fSegs[indx++] = indx2+n*2+j;
      }
      //outside cilinder
      indx2 = i*n*2+n;
      for (j = 0; j < n; j++) {
         buffer.fSegs[indx++] = c+3;
         buffer.fSegs[indx++] = indx2+j;
         buffer.fSegs[indx++] = indx2+n*2+j;
      }
   }

   //left & right sections, number of segments: 2*(fNz-2)
   //          special case number of segments: 0
   if (!specialCase) {
      for (i = 1; i < (fNz-1); i++) {
         for (j = 0; j < 2; j++) {
            buffer.fSegs[indx++] = c;
            buffer.fSegs[indx++] =  2*i    * n + j*(n-1);
            buffer.fSegs[indx++] = (2*i+1) * n + j*(n-1);
         }
      }
   }

   Int_t m = n - 1 + (specialCase == kTRUE);
   indx = 0;

   //bottom & top, number of polygons: 2*(n-1)
   // special case number of polygons: 2*n
   for (j = 0; j < n-1; j++) {
      buffer.fPols[indx++] = c+3;
      buffer.fPols[indx++] = 4;
      buffer.fPols[indx++] = 2*fNz*m+j;
      buffer.fPols[indx++] = m+j;
      buffer.fPols[indx++] = 2*fNz*m+j+1;
      buffer.fPols[indx++] = j;
   }
   for (j = 0; j < n-1; j++) {
      buffer.fPols[indx++] = c+3;
      buffer.fPols[indx++] = 4;
      buffer.fPols[indx++] = 2*fNz*m+n+j;
      buffer.fPols[indx++] = (fNz*2-2)*m+j;
      buffer.fPols[indx++] = 2*fNz*m+n+j+1;
      buffer.fPols[indx++] = (fNz*2-2)*m+m+j;
   }
   if (specialCase) {
      buffer.fPols[indx++] = c+3;
      buffer.fPols[indx++] = 4;
      buffer.fPols[indx++] = 2*fNz*m+j;
      buffer.fPols[indx++] = m+j;
      buffer.fPols[indx++] = 2*fNz*m;
      buffer.fPols[indx++] = j;

      buffer.fPols[indx++] = c+3;
      buffer.fPols[indx++] = 4;
      buffer.fPols[indx++] = 2*fNz*m+n+j;
      buffer.fPols[indx++] = (fNz*2-2)*m+j;
      buffer.fPols[indx++] = 2*fNz*m+n;
      buffer.fPols[indx++] = (fNz*2-2)*m+m+j;
   }
   for (k = 0; k < (fNz-1); k++) {
      for (j = 0; j < n-1; j++) {
         buffer.fPols[indx++] = c;
         buffer.fPols[indx++] = 4;
         buffer.fPols[indx++] = 2*k*m+j;
         buffer.fPols[indx++] = fNz*2*m+(2*k+2)*n+j+1;
         buffer.fPols[indx++] = (2*k+2)*m+j;
         buffer.fPols[indx++] = fNz*2*m+(2*k+2)*n+j;
      }
      for (j = 0; j < n-1; j++) {
         buffer.fPols[indx++] = c+1;
         buffer.fPols[indx++] = 4;
         buffer.fPols[indx++] = (2*k+1)*m+j;
         buffer.fPols[indx++] = fNz*2*m+(2*k+3)*n+j;
         buffer.fPols[indx++] = (2*k+3)*m+j;
         buffer.fPols[indx++] = fNz*2*m+(2*k+3)*n+j+1;
      }

      if (specialCase) {
         buffer.fPols[indx++] = c;
         buffer.fPols[indx++] = 4;
         buffer.fPols[indx++] = 2*k*m+j;
         buffer.fPols[indx++] = fNz*2*m+(2*k+2)*n;
         buffer.fPols[indx++] = (2*k+2)*m+j;
         buffer.fPols[indx++] = fNz*2*m+(2*k+2)*n+j;

         buffer.fPols[indx++] = c+1;
         buffer.fPols[indx++] = 4;
         buffer.fPols[indx++] = (2*k+1)*m+j;
         buffer.fPols[indx++] = fNz*2*m+(2*k+3)*n+j;
         buffer.fPols[indx++] = (2*k+3)*m+j;
         buffer.fPols[indx++] = fNz*2*m+(2*k+3)*n;
      }
   }

   if (!specialCase) {
      indx2 = fNz*2*(n-1);
      for (k = 0; k < (fNz-1); k++) {
         buffer.fPols[indx++] = c+2;
         buffer.fPols[indx++] = 4;
         buffer.fPols[indx++] = k==0 ? indx2 : indx2+2*fNz*n+2*(k-1);
         buffer.fPols[indx++] = indx2+2*(k+1)*n;
         buffer.fPols[indx++] = indx2+2*fNz*n+2*k;
         buffer.fPols[indx++] = indx2+(2*k+3)*n;

         buffer.fPols[indx++] = c+2;
         buffer.fPols[indx++] = 4;
         buffer.fPols[indx++] = k==0 ? indx2+n-1 : indx2+2*fNz*n+2*(k-1)+1;
         buffer.fPols[indx++] = indx2+(2*k+3)*n+n-1;
         buffer.fPols[indx++] = indx2+2*fNz*n+2*k+1;
         buffer.fPols[indx++] = indx2+2*(k+1)*n+n-1;
      }
      buffer.fPols[indx-8] = indx2+n;
      buffer.fPols[indx-2] = indx2+2*n-1;
   }

   return kTRUE;
}
