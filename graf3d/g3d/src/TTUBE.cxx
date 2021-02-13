// @(#)root/g3d:$Id$
// Author: Nenad Buncic   18/09/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TTUBE.h"
#include "TNode.h"
#include "TBuffer.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TGeometry.h"
#include "TMath.h"

ClassImp(TTUBE);

/** \class TTUBE
\ingroup g3d
A tube.

\image html g3d_tube.png

It has 6 parameters:

  - name:       name of the shape
  - title:      shape's title
  - material:  (see TMaterial)
  - rmin:       inside radius
  - rmax:       outside radius
  - dz:         half length in z
*/

////////////////////////////////////////////////////////////////////////////////
/// TUBE shape default constructor

TTUBE::TTUBE()
{
   fCoTab       = 0;
   fSiTab       = 0;
   fAspectRatio = 1;
   fDz          = 0.;
   fNdiv        = 0;
   fRmin        = 0.;
   fRmax        = 0.;
}

////////////////////////////////////////////////////////////////////////////////
/// TUBE shape normal constructor

TTUBE::TTUBE(const char *name, const char *title, const char *material, Float_t rmin, Float_t rmax, Float_t dz,Float_t aspect)
      : TShape(name, title,material)
{
   fRmin  = rmin;
   fRmax  = rmax;

   fDz   = dz;
   fNdiv = 0;

   fCoTab = 0;
   fSiTab = 0;

   fAspectRatio = aspect;

   MakeTableOfCoSin();
}

////////////////////////////////////////////////////////////////////////////////
/// TUBE shape "simplified" constructor

TTUBE::TTUBE(const char *name, const char *title, const char *material, Float_t rmax, Float_t dz)
      : TShape(name, title,material)
{
   fRmin  = 0;
   fRmax  = rmax;

   fDz   = dz;
   fNdiv = 0;

   fCoTab = 0;
   fSiTab = 0;

   fAspectRatio = 1;

   MakeTableOfCoSin();
}

////////////////////////////////////////////////////////////////////////////////
///copy constructor

TTUBE::TTUBE(const TTUBE& tu) :
  TShape(tu),
  fRmin(tu.fRmin),
  fRmax(tu.fRmax),
  fDz(tu.fDz),
  fNdiv(tu.fNdiv),
  fAspectRatio(tu.fAspectRatio),
  fSiTab(tu.fSiTab),
  fCoTab(tu.fCoTab)
{
}

////////////////////////////////////////////////////////////////////////////////
/// assignment operator

TTUBE& TTUBE::operator=(const TTUBE& tu)
{
   if(this!=&tu) {
      TShape::operator=(tu);
      fRmin=tu.fRmin;
      fRmax=tu.fRmax;
      fDz=tu.fDz;
      fNdiv=tu.fNdiv;
      fAspectRatio=tu.fAspectRatio;
      fSiTab=tu.fSiTab;
      fCoTab=tu.fCoTab;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Make table of sine and cosine.

void TTUBE::MakeTableOfCoSin() const // Internal cache - const so other const fn can use
{
   const Double_t pi  = TMath::ATan(1) * 4.0;

   Int_t j;
   Int_t n = GetNumberOfDivisions ();
   if (fCoTab) delete [] fCoTab; // Delete the old tab if any
   fCoTab = new Double_t [n];
   if (!fCoTab ) {
      Error("MakeTableOfCoSin()","No cos table done");
      return;
   }

   if (fSiTab) delete [] fSiTab; // Delete the old tab if any
   fSiTab = new Double_t [n];
   if (!fSiTab ) {
      Error("MakeTableOfCoSin()","No sin table done");
      return;
   }

   Double_t range = 2*pi;

   Double_t angstep = range/n;

   Double_t ph = 0;
   for (j = 0; j < n; j++) {
      ph = j*angstep;
      fCoTab[j] = TMath::Cos(ph);
      fSiTab[j] = TMath::Sin(ph);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// TUBE shape default destructor

TTUBE::~TTUBE()
{
   delete [] fCoTab;
   delete [] fSiTab;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from point px,py to a TUBE
///
/// Compute the closest distance of approach from point px,py to each
/// computed outline point of the TUBE.

Int_t TTUBE::DistancetoPrimitive(Int_t px, Int_t py)
{
   Int_t n = GetNumberOfDivisions();
   Int_t numPoints = n*4;
   return ShapeDistancetoPrimitive(numPoints,px,py);
}

////////////////////////////////////////////////////////////////////////////////
/// Set number of divisions used to draw this tube

void TTUBE::SetNumberOfDivisions (Int_t ndiv)
{
   fNdiv = ndiv;
   MakeTableOfCoSin();
}

////////////////////////////////////////////////////////////////////////////////
/// Create TUBE points

void TTUBE::SetPoints(Double_t *points) const
{
   Int_t j, n;
   Int_t indx = 0;

   n = GetNumberOfDivisions();

   if (points) {
      if (!fCoTab)   MakeTableOfCoSin();
      for (j = 0; j < n; j++) {
         points[indx+6*n] = points[indx] = fRmin * fCoTab[j];
         indx++;
         points[indx+6*n] = points[indx] = fAspectRatio*fRmin * fSiTab[j];
         indx++;
         points[indx+6*n] = fDz;
         points[indx]     =-fDz;
         indx++;
      }
      for (j = 0; j < n; j++) {
         points[indx+6*n] = points[indx] = fRmax * fCoTab[j];
         indx++;
         points[indx+6*n] = points[indx] = fAspectRatio*fRmax * fSiTab[j];
         indx++;
         points[indx+6*n]= fDz;
         points[indx]    =-fDz;
         indx++;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set segments and polygons.

void TTUBE::SetSegsAndPols(TBuffer3D & buffer) const
{
   Int_t i, j;
   Int_t n = GetNumberOfDivisions();
   Int_t c = GetBasicColor();

   for (i = 0; i < 4; i++) {
      for (j = 0; j < n; j++) {
         buffer.fSegs[(i*n+j)*3  ] = c;
         buffer.fSegs[(i*n+j)*3+1] = i*n+j;
         buffer.fSegs[(i*n+j)*3+2] = i*n+j+1;
      }
      buffer.fSegs[(i*n+j-1)*3+2] = i*n;
   }
   for (i = 4; i < 6; i++) {
      for (j = 0; j < n; j++) {
         buffer.fSegs[(i*n+j)*3  ] = c+1;
         buffer.fSegs[(i*n+j)*3+1] = (i-4)*n+j;
         buffer.fSegs[(i*n+j)*3+2] = (i-2)*n+j;
      }
   }
   for (i = 6; i < 8; i++) {
      for (j = 0; j < n; j++) {
         buffer.fSegs[(i*n+j)*3  ] = c;
         buffer.fSegs[(i*n+j)*3+1] = 2*(i-6)*n+j;
         buffer.fSegs[(i*n+j)*3+2] = (2*(i-6)+1)*n+j;
      }
   }

   Int_t indx = 0;
   i=0;
   for (j = 0; j < n; j++) {
      indx = 6*(i*n+j);
      buffer.fPols[indx  ] = c;
      buffer.fPols[indx+1] = 4;
      buffer.fPols[indx+5] = i*n+j;
      buffer.fPols[indx+4] = (4+i)*n+j;
      buffer.fPols[indx+3] = (2+i)*n+j;
      buffer.fPols[indx+2] = (4+i)*n+j+1;
   }
   buffer.fPols[indx+2] = (4+i)*n;
   i=1;
   for (j = 0; j < n; j++) {
      indx = 6*(i*n+j);
      buffer.fPols[indx  ] = c;
      buffer.fPols[indx+1] = 4;
      buffer.fPols[indx+2] = i*n+j;
      buffer.fPols[indx+3] = (4+i)*n+j;
      buffer.fPols[indx+4] = (2+i)*n+j;
      buffer.fPols[indx+5] = (4+i)*n+j+1;
   }
   buffer.fPols[indx+5] = (4+i)*n;
   i=2;
   for (j = 0; j < n; j++) {
      indx = 6*(i*n+j);
      buffer.fPols[indx  ] = c+i;
      buffer.fPols[indx+1] = 4;
      buffer.fPols[indx+2] = (i-2)*2*n+j;
      buffer.fPols[indx+3] = (4+i)*n+j;
      buffer.fPols[indx+4] = ((i-2)*2+1)*n+j;
      buffer.fPols[indx+5] = (4+i)*n+j+1;
   }
   buffer.fPols[indx+5] = (4+i)*n;
   i=3;
   for (j = 0; j < n; j++) {
      indx = 6*(i*n+j);
      buffer.fPols[indx  ] = c+i;
      buffer.fPols[indx+1] = 4;
      buffer.fPols[indx+5] = (i-2)*2*n+j;
      buffer.fPols[indx+4] = (4+i)*n+j;
      buffer.fPols[indx+3] = ((i-2)*2+1)*n+j;
      buffer.fPols[indx+2] = (4+i)*n+j+1;
   }
   buffer.fPols[indx+2] = (4+i)*n;
}

////////////////////////////////////////////////////////////////////////////////
/// Return total X3D needed by TNode::ls (when called with option "x")

void TTUBE::Sizeof3D() const
{
   Int_t n = GetNumberOfDivisions();

   gSize3D.numPoints += n*4;
   gSize3D.numSegs   += n*8;
   gSize3D.numPolys  += n*4;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TTUBE.

void TTUBE::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         R__b.ReadClassBuffer(TTUBE::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TShape::Streamer(R__b);
      R__b >> fRmin;
      R__b >> fRmax;
      R__b >> fDz;
      R__b >> fNdiv;
      if (R__v > 1) R__b >> fAspectRatio;
      R__b.CheckByteCount(R__s, R__c, TTUBE::IsA());
      //====end of old versions
   } else {
      R__b.WriteClassBuffer(TTUBE::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get buffer 3d.

const TBuffer3D & TTUBE::GetBuffer3D(Int_t reqSections) const
{
   static TBuffer3D buffer(TBuffer3DTypes::kGeneric);

   TShape::FillBuffer3D(buffer, reqSections);

   // TODO: Although we now have a TBuffer3DTube class for
   // tube shapes, we do not use it for old geometry tube, as
   // OGL viewer needs various rotation matrix info we can't easily
   // pass yet. To be revisited.

   // We also do not fill the bounding box as derived classes can adjust shape
   // leave up to viewer to work out
   if (reqSections & TBuffer3D::kRawSizes) {
      Int_t n = GetNumberOfDivisions();
      Int_t nbPnts = 4*n;
      Int_t nbSegs = 8*n;
      Int_t nbPols = 4*n;
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
