// @(#)root/g3d:$Id$
// Author: Rene Brun, Nenad Buncic, Evgueni Tcherniaev, Olivier Couet   18/08/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "RConfigure.h"

#include "TVirtualPad.h"
#include "TView3D.h"
#include "TAxis3D.h"
#include "TPolyLine3D.h"
#include "TVirtualX.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TClass.h"
#include "TList.h"
#include "TPluginManager.h"
#include "TMath.h"

// Remove when TView3Der3DPad fix in ExecuteRotateView() is removed
#include "TVirtualViewer3D.h"

ClassImp(TView3D);

//const Int_t kPerspective = BIT(14);

const Int_t    kCARTESIAN = 1;
const Int_t    kPOLAR     = 2;
const Double_t kRad       = 3.14159265358979323846/180.0;

/** \class TView3D
\ingroup g3d
The 3D view class.

This package was originally written by Evgueni Tcherniaev from IHEP/Protvino.

The original Fortran implementation was adapted to HIGZ/PAW by Olivier Couet and
Evgueni Tcherniaev.

This View class is a subset of the original system. It has been converted to a
C++ class  by Rene Brun.

TView3D creates a 3-D view in the current pad. In this 3D view Lego and Surface
plots can be drawn and also 3D polyline and markers. Most of the time a TView3D
is created automatically when a 3D object needs to be painted in a pad (for
instance a Lego or a Surface plot).

In some case a TView3D should be explicitly. For instance to paint a 3D simple
scene composed of simple objects like polylines and polymarkers.
The following macro gives an example:

Begin_Macro(source)
{
   cV3D = new TCanvas("cV3D","PolyLine3D & PolyMarker3D Window",200,10,500,500);

   // Creating a view
   TView3D *view = (TView3D*) TView::CreateView(1);
   view->SetRange(5,5,5,25,25,25);

   // Create a first PolyLine3D
   TPolyLine3D *pl3d1 = new TPolyLine3D(6);
   pl3d1->SetPoint(0, 10, 20, 10);
   pl3d1->SetPoint(1, 15, 15, 15);
   pl3d1->SetPoint(2, 20, 20, 20);
   pl3d1->SetPoint(3, 20, 10, 20);
   pl3d1->SetPoint(4, 10, 10, 20);
   pl3d1->SetPoint(5, 10, 10, 10);

   // Create a first PolyMarker3D
   TPolyMarker3D *pm3d1 = new TPolyMarker3D(9);
   pm3d1->SetPoint( 0, 10, 10, 10);
   pm3d1->SetPoint( 1, 20, 20, 20);
   pm3d1->SetPoint( 2, 10, 20, 20);
   pm3d1->SetPoint( 3, 10, 10, 20);
   pm3d1->SetPoint( 4, 20, 20, 10);
   pm3d1->SetPoint( 5, 20, 10, 10);
   pm3d1->SetPoint( 6, 20, 10, 20);
   pm3d1->SetPoint( 7, 10, 20, 10);
   pm3d1->SetPoint( 8, 15, 15, 15);
   pm3d1->SetMarkerSize(2);
   pm3d1->SetMarkerColor(4);
   pm3d1->SetMarkerStyle(2);

   // Draw
   pl3d1->Draw();
   pm3d1->Draw();
}
End_Macro


Several coordinate systems are available:

  - Cartesian
  - Polar
  - Cylindrical
  - Spherical
  - PseudoRapidity/Phi
*/

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TView3D::TView3D() :TView()
{
   fSystem = 0;
   fOutline = 0;
   fDefaultOutline = kFALSE;
   fAutoRange      = kFALSE;
   fChanged        = kFALSE;

   fPsi = 0;
   Int_t i;
   for (i = 0; i < 3; i++) {
      fRmin[i] = 0;
      fRmax[i] = 1;
      fX1[i] = fX2[i] = fY1[i] = fY2[i] = fZ1[i] = fZ2[i] = 0;
   }

   if (gPad) {
      fLongitude = -90 - gPad->GetPhi();
      fLatitude  =  90 - gPad->GetTheta();
   } else {
      fLongitude = 0;
      fLatitude  = 0;
   }
   Int_t irep = 1;
   ResetView(fLongitude, fLatitude, fPsi, irep);
}

////////////////////////////////////////////////////////////////////////////////
/// TView3D constructor
///
///  Creates a 3-D view in the current pad
///  rmin[3], rmax[3] are the limits of the object depending on
///  the selected coordinate system
///
///  Before drawing a 3-D object in a pad, a 3-D view must be created.
///  Note that a view is automatically created when drawing legos or surfaces.
///
///  The coordinate system is selected via system:
///  - system = 1  Cartesian
///  - system = 2  Polar
///  - system = 3  Cylindrical
///  - system = 4  Spherical
///  - system = 5  PseudoRapidity/Phi

TView3D::TView3D(Int_t system, const Double_t *rmin, const Double_t *rmax) : TView()
{
   Int_t irep;

   SetBit(kMustCleanup);

   fSystem = system;
   fOutline = 0;
   fDefaultOutline = kFALSE;
   fAutoRange      = kFALSE;
   fChanged        = kFALSE;

   if (system == kCARTESIAN || system == kPOLAR || system == 11) fPsi = 0;
   else fPsi = 90;

   // By default pad range in 3-D view is (-1,-1,1,1), so ...
   if (gPad) gPad->Range(-1, -1, 1, 1);
   fAutoRange      = kFALSE;

   Int_t i;
   for (i = 0; i < 3; i++) {
      if (rmin) fRmin[i] = rmin[i];
      else      fRmin[i] = 0;
      if (rmax) fRmax[i] = rmax[i];
      else      fRmax[i] = 1;
      fX1[i] = fX2[i] = fY1[i] = fY2[i] = fZ1[i] = fZ2[i] = 0;
   }

   if (gPad) {
      fLongitude = -90 - gPad->GetPhi();
      fLatitude  =  90 - gPad->GetTheta();
   } else {
      fLongitude = 0;
      fLatitude  = 0;
   }
   ResetView(fLongitude, fLatitude, fPsi, irep);

   if (gPad) gPad->SetView(this);
   if (system == 11) SetPerspective();
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TView3D::TView3D(const TView3D& tv)
  :TView(tv),
   fLatitude(tv.fLatitude),
   fLongitude(tv.fLongitude),
   fPsi(tv.fPsi),
   fDview(tv.fDview),
   fDproj(tv.fDproj),
   fUpix(tv.fUpix),
   fVpix(tv.fVpix),
   fSystem(tv.fSystem),
   fOutline(tv.fOutline),
   fDefaultOutline(tv.fDefaultOutline),
   fAutoRange(tv.fAutoRange),
   fChanged(tv.fChanged)
{
   for (Int_t i=0; i<16; i++) {
      fTN[i]=tv.fTN[i];
      fTB[i]=tv.fTB[i];
      fTnorm[i]=tv.fTnorm[i];
      fTback[i]=tv.fTback[i];
   }
   for(Int_t i=0; i<3; i++) {
      fRmax[i]=tv.fRmax[i];
      fRmin[i]=tv.fRmin[i];
      fX1[i]=tv.fX1[i];
      fX2[i]=tv.fX2[i];
      fY1[i]=tv.fY1[i];
      fY2[i]=tv.fY2[i];
      fZ1[i]=tv.fZ1[i];
      fZ2[i]=tv.fZ2[i];
   }
   for(Int_t i=0; i<4; i++)
      fUVcoord[i]=tv.fUVcoord[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TView3D& TView3D::operator=(const TView3D& tv)
{
   if (this!=&tv) {
      TView::operator=(tv);
      fLatitude=tv.fLatitude;
      fLongitude=tv.fLongitude;
      fPsi=tv.fPsi;
      fDview=tv.fDview;
      fDproj=tv.fDproj;
      fUpix=tv.fUpix;
      fVpix=tv.fVpix;
      fSystem=tv.fSystem;
      fOutline=tv.fOutline;
      fDefaultOutline=tv.fDefaultOutline;
      fAutoRange=tv.fAutoRange;
      fChanged=tv.fChanged;
      for(Int_t i=0; i<16; i++) {
         fTN[i]=tv.fTN[i];
         fTB[i]=tv.fTB[i];
         fTnorm[i]=tv.fTnorm[i];
         fTback[i]=tv.fTback[i];
      }
      for(Int_t i=0; i<3; i++) {
         fRmax[i]=tv.fRmax[i];
         fRmin[i]=tv.fRmin[i];
         fX1[i]=tv.fX1[i];
         fX2[i]=tv.fX2[i];
         fY1[i]=tv.fY1[i];
         fY2[i]=tv.fY2[i];
         fZ1[i]=tv.fZ1[i];
         fZ2[i]=tv.fZ2[i];
      }
      for(Int_t i=0; i<4; i++)
         fUVcoord[i]=tv.fUVcoord[i];
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// TView3D default destructor.

TView3D::~TView3D()
{
   if (fOutline) fOutline->Delete();
   delete fOutline;
   fOutline = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Define axis vertices.
///
/// Input:
///  - ANG     - angle between X and Y axis (not used anymore)
///
///    Output:
///  - AV(3,8) - axis vertices
///  - IX1     - 1st point of X-axis (x-min)
///  - IX2     - 2nd point of X-axis (x-max)
///  - IY1     - 1st point of Y-axis (y-min)
///  - IY2     - 2nd point of Y-axis (y-max)
///  - IZ1     - 1st point of Z-axis (z-min)
///  - IZ2     - 2nd point of Z-axis (z-max)

void TView3D::AxisVertex(Double_t, Double_t *av, Int_t &ix1, Int_t &ix2, Int_t &iy1, Int_t &iy2, Int_t &iz1, Int_t &iz2)
{
   Double_t  p[8][3] = {
      { fRmin[0], fRmin[1], fRmin[2] },
      { fRmax[0], fRmin[1], fRmin[2] },
      { fRmax[0], fRmax[1], fRmin[2] },
      { fRmin[0], fRmax[1], fRmin[2] },
      { fRmin[0], fRmin[1], fRmax[2] },
      { fRmax[0], fRmin[1], fRmax[2] },
      { fRmax[0], fRmax[1], fRmax[2] },
      { fRmin[0], fRmax[1], fRmax[2] }
   };
   Int_t inodes[4][8] = {
      { 2,3,4,1, 6,7,8,5 }, // x+, y+
      { 3,4,1,2, 7,8,5,6 }, // x-, y+
      { 1,2,3,4, 5,6,7,8 }, // x+, y-
      { 4,1,2,3, 8,5,6,7 }  // x-, y-
   };
   Int_t ixyminmax[16][4] = { //               8
      { 3,2, 1,2 }, // x+, y+, z+, z-up     5 / \ 7
      { 2,1, 3,2 }, // x-, y+, z+, z-up      |\6/|
      { 1,2, 2,3 }, // x+, y-, z+, z-up      | | |   Top view
      { 2,3, 2,1 }, // x-, y-, z+, z-up     1 \|/ 3
                    //                         2         6
      { 4,1, 4,3 }, // x+, y+, z-, z-up               5 /|\ 7
      { 3,4, 4,1 }, // x-, y+, z-, z-up                | | |
      { 4,3, 1,4 }, // x+, y-, z-, z-up   Bottom view  |/2\|
      { 1,4, 3,4 }, // x-, y-, z-, z-up               1 \ / 3
                    //                         2         4
      { 8,5, 8,7 }, // x+, y+, z+, z-down   3 /|\ 1
      { 7,8, 8,5 }, // x-, y+, z+, z-down    | | |
      { 8,7, 5,8 }, // x+, y-, z+, z-down    |/6\|   Bottom view
      { 5,8, 7,8 }, // x-, y-, z+, z-down   7 \ / 5
                    //                         8         4
      { 7,6, 5,6 }, // x+, y+, z-, z-down             3 / \ 1
      { 6,5, 7,6 }, // x-, y+, z-, z-down              |\2/|
      { 5,6, 6,7 }, // x+, y-, z-, z-down   Top view   | | |
      { 6,7, 6,5 }  // x-, y-, z-, z-down             7 \|/ 5
   };               //                                   6

   //        Set vertices
   Int_t icase = 0;
   if (fTnorm[ 8] <= 0) icase += 1; // z projection of (1,0,0)
   if (fTnorm[ 9] <= 0) icase += 2; // z projection of (0,1,0)
   for (Int_t i=0; i<8; ++i) {
      Int_t k = inodes[icase][i] - 1;
      av[i*3+0] = p[k][0];
      av[i*3+1] = p[k][1];
      av[i*3+2] = p[k][2];
   }

   //        Set indices for min and max
   if (fTnorm[10] < 0) icase += 4; // z projection of (0,0,1)
   if (fTnorm[ 6] < 0) icase += 8; // y projection of (0,0,1)
   ix1 = ixyminmax[icase][0];
   ix2 = ixyminmax[icase][1];
   iy1 = ixyminmax[icase][2];
   iy2 = ixyminmax[icase][3];
   iz1 = (icase < 8) ? 1 : 3;
   iz2 = (icase < 8) ? 5 : 7;
}

////////////////////////////////////////////////////////////////////////////////
/// Define perspective view.
///
/// Compute transformation matrix from world coordinates
/// to normalized coordinates (-1 to +1)
///
/// Input :
///  - theta, phi - spherical angles giving the direction of projection
///  - psi - screen rotation angle
///  - cov[3] - center of view
///  - dview - distance from COV to COP (center of projection)
///  - umin, umax, vmin, vmax - view window in projection plane
///  - dproj - distance from COP to projection plane
///  - bcut, fcut - backward/forward range w.r.t projection plane (fcut<=0)
///
/// Output :
///  - nper[16] - normalizing transformation
///
/// compute tr+rot to get COV in origin, view vector parallel to -Z axis, up
/// vector parallel to Y.
///
/// ~~~ {.cpp}
///                         ^Yv   UP ^  proj. plane
///                        |        |   /|
///                       |        |  /  |
///                      |   dproj  /  x--- center of window (COW)
///                 COV |----------|--x--|------------> Zv
///                              /           | VRP'z
///                       /   --->      |  /
///                /     VPN       |/
///               Xv
/// ~~~
///
///   1. translate COP to origin of MARS : Tper = T(-copx, -copy, -copz)
///   2. rotate VPN : R = Rz(-psi)*Rx(-theta)*Rz(-phi) (inverse Euler)
///   3. left-handed screen reference to right-handed one of MARS : Trl
///
///   T12 = Tper*R*Trl


void TView3D::DefinePerspectiveView()
{
   Double_t t12[16];
   Double_t cov[3];
   Int_t i;
   for (i=0; i<3; i++) cov[i] = 0.5*(fRmax[i]+fRmin[i]);

   Double_t c1 = TMath::Cos(fPsi*kRad);
   Double_t s1 = TMath::Sin(fPsi*kRad);
   Double_t c2 = TMath::Cos(fLatitude*kRad);
   Double_t s2 = TMath::Sin(fLatitude*kRad);
   Double_t s3 = TMath::Cos(fLongitude*kRad);
   Double_t c3 = -TMath::Sin(fLongitude*kRad);

   t12[0] =  c1*c3 - s1*c2*s3;
   t12[4] =  c1*s3 + s1*c2*c3;
   t12[8] =  s1*s2;
   t12[3] =  0;

   t12[1] =  -s1*c3 - c1*c2*s3;
   t12[5] = -s1*s3 + c1*c2*c3;
   t12[9] =  c1*s2;
   t12[7] =  0;

   t12[2] =  s2*s3;
   t12[6] =  -s2*c3;
   t12[10] = c2;      // contains Trl
   t12[11] =  0;

   // translate with -COP (before rotation):
   t12[12] = -(cov[0]*t12[0]+cov[1]*t12[4]+cov[2]*t12[8]);
   t12[13] = -(cov[0]*t12[1]+cov[1]*t12[5]+cov[2]*t12[9]);
   t12[14] = -(cov[0]*t12[2]+cov[1]*t12[6]+cov[2]*t12[10]);
   t12[15] =  1;

   // translate with (0, 0, -dview) after rotation

   t12[14] -= fDview;

   // reflection on Z :
   t12[2]  *= -1;
   t12[6]  *= -1;
   t12[10] *= -1;
   t12[14] *= -1;

   // Now we shear the center of window from (0.5*(umin+umax), 0.5*(vmin+vmax), dproj)
   //                                     to (0, 0, dproj)

   Double_t a2 = -fUVcoord[0]/fDproj;   // shear coef. on x
   Double_t b2 = -fUVcoord[1]/fDproj;   // shear coef. on y

   //               | 1  0  0  0 |
   //  SHz(a2,b2) = | 0  1  0  0 |
   //               | a2 b2 1  0 |
   //               | 0  0  0  1 |

   fTnorm[0] = t12[0] + a2*t12[2];
   fTnorm[1] = t12[1] + b2*t12[2];
   fTnorm[2] = t12[2];
   fTnorm[3] = 0;

   fTnorm[4] = t12[4] + a2*t12[6];
   fTnorm[5] = t12[5] + b2*t12[6];
   fTnorm[6] = t12[6];
   fTnorm[7] = 0;

   fTnorm[8]  = t12[8] + a2*t12[10];
   fTnorm[9]  = t12[9] + b2*t12[10];
   fTnorm[10] = t12[10];
   fTnorm[11] = 0;

   fTnorm[12] = t12[12] + a2*t12[14];
   fTnorm[13] = t12[13] + b2*t12[14];
   fTnorm[14] = t12[14];
   fTnorm[15] = 1;

   // Scale so that the view volume becomes the canonical one
   //
   // Sper = (2/(umax-umin), 2/(vmax-vmin), 1/dproj
   //
   Double_t sz = 1./fDproj;
   Double_t sx = 1./fUVcoord[2];
   Double_t sy = 1./fUVcoord[3];

   fTnorm[0] *= sx;
   fTnorm[4] *= sx;
   fTnorm[8] *= sx;
   fTnorm[1] *= sy;
   fTnorm[5] *= sy;
   fTnorm[9] *= sy;
   fTnorm[2] *= sz;
   fTnorm[6] *= sz;
   fTnorm[10] *= sz;
   fTnorm[12] *= sx;
   fTnorm[13] *= sy;
   fTnorm[14] *= sz;
}

////////////////////////////////////////////////////////////////////////////////
/// Define view direction (in spherical coordinates)
///
/// Compute transformation matrix from world coordinates
/// to normalized coordinates (-1 to +1)
///
/// Input:
///  - S(3)    - scale factors
///  - C(3)    - centre of scope
///  - COSPHI  - longitude COS
///  - SINPHI  - longitude SIN
///  - COSTHE  - latitude COS (angle between +Z and view direction.)
///  - SINTHE  - latitude SIN
///  - COSPSI  - screen plane rotation angle COS
///  - SINPSI  - screen plane rotation angle SIN

void TView3D::DefineViewDirection(const Double_t *s, const Double_t *c,
                                Double_t cosphi, Double_t sinphi,
                                Double_t costhe, Double_t sinthe,
                                Double_t cospsi, Double_t sinpsi,
                                Double_t *tnorm, Double_t *tback)
{
   if (IsPerspective()) {
      DefinePerspectiveView();
      return;
   }
   Int_t i, k;
   Double_t tran[16]   /* was [4][4] */, rota[16]      /* was [4][4] */;
   Double_t c1, c2, c3, s1, s2, s3, scalex, scaley, scalez;

   // Parameter adjustments
   tback -= 5;
   tnorm -= 5;

   scalex = s[0];
   scaley = s[1];
   scalez = s[2];

   //*-*-        S E T   T R A N S L A T I O N   M A T R I X
   tran[0] = 1 / scalex;
   tran[1] = 0;
   tran[2] = 0;
   tran[3] = -c[0] / scalex;

   tran[4] = 0;
   tran[5] = 1 / scaley;
   tran[6] = 0;
   tran[7] = -c[1] / scaley;

   tran[8] = 0;
   tran[9] = 0;
   tran[10] = 1 / scalez;
   tran[11] = -c[2] / scalez;

   tran[12] = 0;
   tran[13] = 0;
   tran[14] = 0;
   tran[15] = 1;

   //*-*-        S E T    R O T A T I O N   M A T R I X
   //    ( C(PSI) S(PSI) 0)   (1      0          0 )   ( C(90+PHI) S(90+PHI) 0)
   //    (-S(PSI) C(PSI) 0) * (0  C(THETA) S(THETA)) * (-S(90+PHI) C(90+PHI) 0)
   //    (   0      0    1)   (0 -S(THETA) C(THETA))   (     0           0   1)
   c1 = cospsi;
   s1 = sinpsi;
   c2 = costhe;
   s2 = sinthe;
   c3 = -sinphi;
   s3 = cosphi;

   rota[0] = c1*c3 - s1*c2*s3;
   rota[1] = c1*s3 + s1*c2*c3;
   rota[2] = s1*s2;
   rota[3] = 0;

   rota[4] = -s1*c3 - c1* c2*s3;
   rota[5] = -s1*s3 + c1* c2*c3;
   rota[6] = c1*s2;
   rota[7] = 0;

   rota[8] = s2*s3;
   rota[9] = -s2*c3;
   rota[10] = c2;
   rota[11] = 0;

   rota[12] = 0;
   rota[13] = 0;
   rota[14] = 0;
   rota[15] = 1;

   //*-*-        F I N D   T R A N S F O R M A T I O N   M A T R I X
   for (i = 1; i <= 3; ++i) {
      for (k = 1; k <= 4; ++k) {
         tnorm[k + (i << 2)] = rota[(i << 2) - 4]*tran[k - 1] + rota[(i
                 << 2) - 3]*tran[k + 3] + rota[(i << 2) - 2]*tran[k +7]
                 + rota[(i << 2) - 1]*tran[k + 11];
      }
   }

   //*-*-        S E T   B A C K   T R A N S L A T I O N   M A T R I X
   tran[0] = scalex;
   tran[3] = c[0];

   tran[5] = scaley;
   tran[7] = c[1];

   tran[10] = scalez;
   tran[11] = c[2];

   //*-*-        F I N D   B A C K   T R A N S F O R M A T I O N
   for (i = 1; i <= 3; ++i) {
      for (k = 1; k <= 4; ++k) {
         tback[k + (i << 2)] = tran[(i << 2) - 4]*rota[(k << 2) - 4] +
            tran[(i << 2) - 3]*rota[(k << 2) - 3] + tran[(i << 2) -2]
            *rota[(k << 2) - 2] + tran[(i << 2) - 1]*rota[(k <<2) - 1];
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the outline of a cube while rotating a 3-d object in the pad.

void TView3D::DrawOutlineCube(TList *outline, Double_t *rmin, Double_t *rmax)
{
   TPolyLine3D::DrawOutlineCube(outline,rmin,rmax);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute action corresponding to one event.

void TView3D::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   ExecuteRotateView(event,px,py);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute action corresponding to one event.
///
/// This member function is called when a object is clicked with the locator
///
/// If Left button clicked in the object area, while the button is kept down
/// the cube representing the surrounding frame for the corresponding
/// new latitude and longitude position is drawn.

void TView3D::ExecuteRotateView(Int_t event, Int_t px, Int_t py)
{
   static Int_t system, framewasdrawn;
   static Double_t xrange, yrange, xmin, ymin, longitude1, latitude1, longitude2, latitude2;
   static Double_t newlatitude, newlongitude, oldlatitude, oldlongitude;
   Double_t dlatitude, dlongitude, x, y;
   Int_t irep = 0;
   Double_t psideg;
   Bool_t opaque  = gPad->OpaqueMoving();

   // All coordinates transformation are from absolute to relative
   if (!gPad->IsEditable()) return;
   gPad->AbsCoordinates(kTRUE);

   switch (event) {

   case kKeyPress :
      fChanged = kTRUE;
      MoveViewCommand(Char_t(px), py);
      break;

   case kMouseMotion:
      gPad->SetCursor(kRotate);
      break;

   case kButton1Down:

      // remember position of the cube
      xmin   = gPad->GetX1();
      ymin   = gPad->GetY1();
      xrange = gPad->GetX2() - xmin;
      yrange = gPad->GetY2() - ymin;
      x      = gPad->PixeltoX(px);
      y      = gPad->PixeltoY(py);
      system = GetSystem();
      framewasdrawn = 0;
      if (system == kCARTESIAN || system == kPOLAR || system == 11 || IsPerspective()) {
         longitude1 = 180*(x-xmin)/xrange;
         latitude1  =  90*(y-ymin)/yrange;
      } else {
         latitude1  =  90*(x-xmin)/xrange;
         longitude1 = 180*(y-ymin)/yrange;
      }
      newlongitude = oldlongitude = -90 - gPad->GetPhi();
      newlatitude  = oldlatitude  =  90 - gPad->GetTheta();
      psideg       = GetPsi();

      // if outline isn't set, make it look like a cube
      if(!fOutline)
         SetOutlineToCube();
      break;

   case kButton1Motion:
   {
      // draw the surrounding frame for the current mouse position
      // first: Erase old frame
      fChanged = kTRUE;
      if (framewasdrawn && !opaque) fOutline->Paint();
      framewasdrawn = 1;
      x = gPad->PixeltoX(px);
      y = gPad->PixeltoY(py);
      if (system == kCARTESIAN || system == kPOLAR || system == 11 || IsPerspective()) {
         longitude2 = 180*(x-xmin)/xrange;
         latitude2  =  90*(y-ymin)/yrange;
      } else {
         latitude2  =  90*(x-xmin)/xrange;
         longitude2 = 180*(y-ymin)/yrange;
      }
      dlongitude   = longitude2   - longitude1;
      dlatitude    = latitude2    - latitude1;
      newlatitude  = oldlatitude  + dlatitude;
      newlongitude = oldlongitude - dlongitude;
      psideg       = GetPsi();
      ResetView(newlongitude, newlatitude, psideg, irep);
      if (!opaque) {
         fOutline->Paint();
      } else {
         psideg = GetPsi();
         SetView(newlongitude, newlatitude, psideg, irep);
         gPad->SetPhi(-90-newlongitude);
         gPad->SetTheta(90-newlatitude);
         gPad->Modified(kTRUE);
      }

      break;
   }
   case kButton1Up:
      if (gROOT->IsEscaped()) {
         gROOT->SetEscape(kFALSE);
         if (opaque) {
            psideg = GetPsi();
            SetView(oldlongitude, oldlatitude, psideg, irep);
            gPad->SetPhi(-90-oldlongitude);
            gPad->SetTheta(90-oldlatitude);
            gPad->Modified(kTRUE);
         }
         break;
      }

      // Temporary fix for 2D drawing problems on pad. fOutline contains
      // a TPolyLine3D object for the rotation box. This will be painted
      // through a newly created TView3Der3DPad instance, which is left
      // behind on pad. This remaining creates 2D drawing problems.
      //
      // This is a TEMPORARY fix - will be removed when proper multiple viewers
      // on pad problems are resolved.
      if (gPad) {
         TVirtualViewer3D *viewer = gPad->GetViewer3D();
         if (viewer && !strcmp(viewer->IsA()->GetName(),"TView3Der3DPad")) {
            gPad->ReleaseViewer3D();
            delete viewer;
         }
      }
      // End fix

      // Recompute new view matrix and redraw
      psideg = GetPsi();
      SetView(newlongitude, newlatitude, psideg, irep);
      gPad->SetPhi(-90-newlongitude);
      gPad->SetTheta(90-newlatitude);
      gPad->Modified(kTRUE);

      // Set line color, style and width
      gVirtualX->SetLineColor(-1);
      gVirtualX->SetLineStyle(-1);
      gVirtualX->SetLineWidth(-1);
      break;
   }

   // set back to default transformation mode
   gPad->AbsCoordinates(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Find Z component of NORMAL in normalized coordinates.
///
/// Input:
///  - X - X-component of NORMAL
///  - Y - Y-component of NORMAL
///  - Z - Z-component of NORMAL
///
/// Output:
///  - ZN - Z-component of NORMAL in normalized coordinates

void TView3D::FindNormal(Double_t x, Double_t  y, Double_t z, Double_t &zn)
{
   zn = x*(fTN[1] * fTN[6] - fTN[2] * fTN[5]) + y*(fTN[2] * fTN[4] -
           fTN[0] * fTN[6]) + z*(fTN[0] * fTN[5] - fTN[1] * fTN[4]);
}

////////////////////////////////////////////////////////////////////////////////
/// Find critical PHI sectors.
///
/// Input:
///  - IOPT    - options:
///      1. from BACK to FRONT 'BF'
///      2. from FRONT to BACK 'FB'
///  - KPHI    - number of phi sectors
///  - APHI(*) - PHI separators (modified internally)
///
/// Output:
///  - IPHI1  - initial sector
///  - IPHI2  - final sector

void TView3D::FindPhiSectors(Int_t iopt, Int_t &kphi, Double_t *aphi, Int_t &iphi1, Int_t &iphi2)
{
   Int_t iphi[2], i, k;
   Double_t dphi;
   Double_t x1, x2, z1, z2, phi1, phi2;

   // Parameter adjustments
   --aphi;

   if (aphi[kphi + 1] == aphi[1]) aphi[kphi + 1] += 360;
   dphi = TMath::Abs(aphi[kphi + 1] - aphi[1]);
   if (dphi != 360) {
      aphi[kphi + 2] = (aphi[1] + aphi[kphi + 1]) / (float)2. + 180;
      aphi[kphi + 3] = aphi[1] + 360;
      kphi += 2;
   }

   //*-*-       F I N D   C R I T I C A L   S E C T O R S
   k = 0;
   for (i = 1; i <= kphi; ++i) {
      phi1 = kRad*aphi[i];
      phi2 = kRad*aphi[i + 1];
      x1 = fTN[0]*TMath::Cos(phi1) + fTN[1]*TMath::Sin(phi1);
      x2 = fTN[0]*TMath::Cos(phi2) + fTN[1]*TMath::Sin(phi2);
      if (x1 >= 0 && x2 > 0) continue;
      if (x1 <= 0 && x2 < 0) continue;
      ++k;
      if (k == 3) break;
      iphi[k - 1] = i;
   }
   if (k != 2) {
      Error("FindPhiSectors", "something strange: num. of critical sector not equal 2");
      iphi1 = 1;
      iphi2 = 2;
      return;
   }

   //*-*-       F I N D   O R D E R   O F   C R I T I C A L   S E C T O R S
   phi1 = kRad*(aphi[iphi[0]] + aphi[iphi[0] + 1]) / (float)2.;
   phi2 = kRad*(aphi[iphi[1]] + aphi[iphi[1] + 1]) / (float)2.;
   z1 = fTN[8]*TMath::Cos(phi1) + fTN[9]*TMath::Sin(phi1);
   z2 = fTN[8]*TMath::Cos(phi2) + fTN[9]*TMath::Sin(phi2);
   if ((z1 <= z2 && iopt == 1) || (z1 > z2 && iopt == 2)) {
      iphi1 = iphi[0];
      iphi2 = iphi[1];
   } else {
      iphi1 = iphi[1];
      iphi2 = iphi[0];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Find critical THETA sectors for given PHI sector.
///
/// Input:
///  - IOPT        - options:
///      1. from BACK to FRONT 'BF'
///      2. from FRONT to BACK 'FB'
///  - PHI         - PHI sector
///  - KTH         - number of THETA sectors
///  - ATH(*)      - THETA separators (modified internally)
///
/// Output:
///  - ITH1  - initial sector
///  - ITH2  - final sector

void TView3D::FindThetaSectors(Int_t iopt, Double_t phi, Int_t &kth, Double_t *ath, Int_t &ith1, Int_t &ith2)
{
   Int_t i, k, ith[2];
   Double_t z1, z2, cosphi, sinphi, tncons, th1, th2, dth;

   // Parameter adjustments
   --ath;

   // Function Body
   dth = TMath::Abs(ath[kth + 1] - ath[1]);
   if (dth != 360) {
      ath[kth + 2] = 0.5*(ath[1] + ath[kth + 1]) + 180;
      ath[kth + 3] = ath[1] + 360;
      kth += 2;
   }

   //*-*-       F I N D   C R I T I C A L   S E C T O R S
   cosphi = TMath::Cos(phi*kRad);
   sinphi = TMath::Sin(phi*kRad);
   k = 0;
   for (i = 1; i <= kth; ++i) {
      th1 = kRad*ath[i];
      th2 = kRad*ath[i + 1];
      FindNormal(TMath::Cos(th1)*cosphi, TMath::Cos(th1)*sinphi, -TMath::Sin(th1), z1);
      FindNormal(TMath::Cos(th2)*cosphi, TMath::Cos(th2)*sinphi, -TMath::Sin(th2), z2);
      if (z1 >= 0 && z2 > 0) continue;
      if (z1 <= 0 && z2 < 0) continue;
      ++k;
      if (k == 3) break;
      ith[k - 1] = i;
   }
   if (k != 2) {
      Error("FindThetaSectors", "Something strange: num. of critical sectors not equal 2");
      ith1 = 1;
      ith2 = 2;
      return;
   }

   //*-*-       F I N D   O R D E R   O F   C R I T I C A L   S E C T O R S
   tncons = fTN[8]*TMath::Cos(phi*kRad) + fTN[9]*TMath::Sin(phi*kRad);
   th1    = kRad*(ath[ith[0]] + ath[ith[0] + 1]) / (float)2.;
   th2    = kRad*(ath[ith[1]] + ath[ith[1] + 1]) / (float)2.;
   z1     = tncons*TMath::Sin(th1) + fTN[10]*TMath::Cos(th1);
   z2     = tncons*TMath::Sin(th2) + fTN[10]*TMath::Cos(th2);
   if ((z1 <= z2 && iopt == 1) || (z1 > z2 && iopt == 2)) {
      ith1 = ith[0];
      ith2 = ith[1];
   } else {
      ith1 = ith[1];
      ith2 = ith[0];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Find centre of a MIN-MAX scope and scale factors
///
/// Output:
///  - SCALE(3)  - scale factors
///  - CENTER(3) - centre
///  - IREP      - reply (-1 if error in min-max)

void TView3D::FindScope(Double_t *scale, Double_t *center, Int_t &irep)
{
   irep = 0;
   Double_t sqrt3 = 0.5*TMath::Sqrt(3.0);

   for (Int_t i = 0; i < 3; i++) {
      if (fRmin[i] >= fRmax[i]) { irep = -1; return;}
      scale[i]  = sqrt3*(fRmax[i] - fRmin[i]);
      center[i] = 0.5*(fRmax[i] + fRmin[i]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return distance to axis from point px,py.
///
/// Algorithm:
///
/// ~~~ {.cpp}
///       A(x1,y1)         P                             B(x2,y2)
///       ------------------------------------------------
///                        I
///                        I
///                        I
///                        I
///                       M(x,y)
///
///  Let us call  a = distance AM     A=a**2
///               b = distance BM     B=b**2
///               c = distance AB     C=c**2
///               d = distance PM     D=d**2
///               u = distance AP     U=u**2
///               v = distance BP     V=v**2     c = u + v
///
///  D = A - U
///  D = B - V  = B -(c-u)**2
///     ==> u = (A -B +C)/2c
/// ~~~

Int_t TView3D::GetDistancetoAxis(Int_t axis, Int_t px, Int_t py, Double_t &ratio)
{
   Double_t x1,y1,x2,y2;
   Double_t x     = px;
   Double_t y     = py;
   ratio = 0;

   if (fSystem != kCARTESIAN) return 9998; // only implemented for Cartesian coordinates
   if (axis == 1) {
      x1 = gPad->XtoAbsPixel(fX1[0]);
      y1 = gPad->YtoAbsPixel(fX1[1]);
      x2 = gPad->XtoAbsPixel(fX2[0]);
      y2 = gPad->YtoAbsPixel(fX2[1]);
   } else if (axis == 2) {
      x1 = gPad->XtoAbsPixel(fY1[0]);
      y1 = gPad->YtoAbsPixel(fY1[1]);
      x2 = gPad->XtoAbsPixel(fY2[0]);
      y2 = gPad->YtoAbsPixel(fY2[1]);
   } else {
      x1 = gPad->XtoAbsPixel(fZ1[0]);
      y1 = gPad->YtoAbsPixel(fZ1[1]);
      x2 = gPad->XtoAbsPixel(fZ2[0]);
      y2 = gPad->YtoAbsPixel(fZ2[1]);
   }
   Double_t xx1   = x  - x1;
   Double_t xx2   = x  - x2;
   Double_t x1x2  = x1 - x2;
   Double_t yy1   = y  - y1;
   Double_t yy2   = y  - y2;
   Double_t y1y2  = y1 - y2;
   Double_t a     = xx1*xx1   + yy1*yy1;
   Double_t b     = xx2*xx2   + yy2*yy2;
   Double_t c     = x1x2*x1x2 + y1y2*y1y2;
   if (c <= 0) return 9999;
   Double_t v     = TMath::Sqrt(c);
   Double_t u     = (a - b + c)/(2*v);
   Double_t d     = TMath::Abs(a - u*u);

   Int_t dist = Int_t(TMath::Sqrt(d) - 0.5);
   ratio = u/v;
   return dist;
}

////////////////////////////////////////////////////////////////////////////////
/// Get maximum view extent.

Double_t TView3D::GetExtent() const
{
   Double_t dx = 0.5*(fRmax[0]-fRmin[0]);
   Double_t dy = 0.5*(fRmax[1]-fRmin[1]);
   Double_t dz = 0.5*(fRmax[2]-fRmin[2]);
   Double_t extent = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
   return extent;
}

////////////////////////////////////////////////////////////////////////////////
/// Get Range function.

void TView3D::GetRange(Float_t *min, Float_t *max)
{
   for (Int_t i = 0; i < 3; max[i] = fRmax[i], min[i] = fRmin[i], i++) { }
}

////////////////////////////////////////////////////////////////////////////////
/// Get Range function.

void TView3D::GetRange(Double_t *min, Double_t *max)
{
   for (Int_t i = 0; i < 3; max[i] = fRmax[i], min[i] = fRmin[i], i++) { }
}

////////////////////////////////////////////////////////////////////////////////
/// Get current window extent.

void TView3D::GetWindow(Double_t &u0, Double_t &v0, Double_t &du, Double_t &dv) const
{
   u0 = fUVcoord[0];
   v0 = fUVcoord[1];
   du = fUVcoord[2];
   dv = fUVcoord[3];
}

////////////////////////////////////////////////////////////////////////////////
/// Check if point is clipped in perspective view.

Bool_t TView3D::IsClippedNDC(Double_t *p) const
{
   if (TMath::Abs(p[0])>p[2]) return kTRUE;
   if (TMath::Abs(p[1])>p[2]) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Transfer point from normalized to world coordinates.
///
/// Input:
///  - PN(3) - point in world coordinate system
///  - PW(3) - point in normalized coordinate system

void TView3D::NDCtoWC(const Float_t* pn, Float_t* pw)
{
   Float_t x = pn[0], y = pn[1], z = pn[2];
   pw[0] = fTback[0]*x + fTback[1]*y + fTback[2]*z  + fTback[3];
   pw[1] = fTback[4]*x + fTback[5]*y + fTback[6]*z  + fTback[7];
   pw[2] = fTback[8]*x + fTback[9]*y + fTback[10]*z + fTback[11];
}

////////////////////////////////////////////////////////////////////////////////
/// Transfer point from normalized to world coordinates.
///
/// Input:
///  - PN(3) - point in world coordinate system
///  - PW(3) - point in normalized coordinate system

void TView3D::NDCtoWC(const Double_t* pn, Double_t* pw)
{
   Double_t x = pn[0], y = pn[1], z = pn[2];
   pw[0] = fTback[0]*x + fTback[1]*y + fTback[2]*z  + fTback[3];
   pw[1] = fTback[4]*x + fTback[5]*y + fTback[6]*z  + fTback[7];
   pw[2] = fTback[8]*x + fTback[9]*y + fTback[10]*z + fTback[11];
}

////////////////////////////////////////////////////////////////////////////////
/// Transfer vector of NORMAL from word to normalized coordinates.
///
/// Input:
///  - PW(3) - vector of NORMAL in word coordinate system
///  - PN(3) - vector of NORMAL in normalized coordinate system

void TView3D::NormalWCtoNDC(const Float_t *pw, Float_t *pn)
{
   Double_t x, y, z, a1, a2, a3, b1, b2, b3, c1, c2, c3;

   x = pw[0];
   y = pw[1];
   z = pw[2];
   a1 = fTnorm[0];
   a2 = fTnorm[1];
   a3 = fTnorm[2];
   b1 = fTnorm[4];
   b2 = fTnorm[5];
   b3 = fTnorm[6];
   c1 = fTnorm[8];
   c2 = fTnorm[9];
   c3 = fTnorm[10];
   pn[0] = x*(b2*c3 - b3*c2) + y*(b3*c1 - b1*c3) + z*(b1*c2 - b2*c1);
   pn[1] = x*(c2*a3 - c3*a2) + y*(c3*a1 - c1*a3) + z*(c1*a2 - c2*a1);
   pn[2] = x*(a2*b3 - a3*b2) + y*(a3*b1 - a1*b3) + z*(a1*b2 - a2*b1);
}

////////////////////////////////////////////////////////////////////////////////
/// Transfer vector of NORMAL from word to normalized coordinates.
///
/// Input:
///  - PW(3) - vector of NORMAL in word coordinate system
///  - PN(3) - vector of NORMAL in normalized coordinate system

void TView3D::NormalWCtoNDC(const Double_t *pw, Double_t *pn)
{
   Double_t x, y, z, a1, a2, a3, b1, b2, b3, c1, c2, c3;

   x = pw[0];
   y = pw[1];
   z = pw[2];
   a1 = fTnorm[0];
   a2 = fTnorm[1];
   a3 = fTnorm[2];
   b1 = fTnorm[4];
   b2 = fTnorm[5];
   b3 = fTnorm[6];
   c1 = fTnorm[8];
   c2 = fTnorm[9];
   c3 = fTnorm[10];
   pn[0] = x*(b2*c3 - b3*c2) + y*(b3*c1 - b1*c3) + z*(b1*c2 - b2*c1);
   pn[1] = x*(c2*a3 - c3*a2) + y*(c3*a1 - c1*a3) + z*(c1*a2 - c2*a1);
   pn[2] = x*(a2*b3 - a3*b2) + y*(a3*b1 - a1*b3) + z*(a1*b2 - a2*b1);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the correct window size for lego and surface plots.
///
///  Set the correct window size for lego and surface plots.
///  And draw the background if necessary.
///
/// Input parameters:
///  - RBACK : Background colour

void TView3D::PadRange(Int_t rback)
{
   Int_t i, k;
   Double_t x, y, z, r1, r2, r3, xx, yy, smax[2];
   Double_t xgraf[6], ygraf[6];

   for (i = 1; i <= 2; ++i) {
      smax[i - 1] = fTnorm[(i << 2) - 1];
      for (k = 1; k <= 3; ++k) {
         if (fTnorm[k + (i << 2) - 5] < 0) {
            smax[i - 1] += fTnorm[k + (i << 2) - 5]*fRmin[k-1];
         } else {
            smax[i - 1] += fTnorm[k + (i << 2) - 5]*fRmax[k-1];
         }
      }
   }

   //*-*- Compute x,y range
   Double_t xmin = -smax[0];
   Double_t xmax = smax[0];
   Double_t ymin = -smax[1];
   Double_t ymax = smax[1];
   Double_t dx   = xmax-xmin;
   Double_t dy   = ymax-ymin;
   Double_t dxr  = dx/(1 - gPad->GetLeftMargin() - gPad->GetRightMargin());
   Double_t dyr  = dy/(1 - gPad->GetBottomMargin() - gPad->GetTopMargin());

   // Range() could change the size of the pad pixmap and therefore should
   // be called before the other paint routines
   gPad->Range(xmin - dxr*gPad->GetLeftMargin(),
      ymin - dyr*gPad->GetBottomMargin(),
      xmax + dxr*gPad->GetRightMargin(),
      ymax + dyr*gPad->GetTopMargin());
   gPad->RangeAxis(xmin, ymin, xmax, ymax);

   //*-*-             Draw the background if necessary
   if (rback > 0) {
      r1 = -1;
      r2 = -1;
      r3 = -1;
      xgraf[0] = -smax[0];
      xgraf[1] = -smax[0];
      xgraf[2] = -smax[0];
      xgraf[3] = -smax[0];
      xgraf[4] =  smax[0];
      xgraf[5] =  smax[0];
      ygraf[0] = -smax[1];
      ygraf[1] =  smax[1];
      ygraf[2] = -smax[1];
      ygraf[3] =  smax[1];
      ygraf[5] =  smax[1];
      ygraf[4] = -smax[1];
      for (i = 1; i <= 8; ++i) {
         x = 0.5*((1 - r1)*fRmin[0] + (r1 + 1)*fRmax[0]);
         y = 0.5*((1 - r2)*fRmin[1] + (r2 + 1)*fRmax[1]);
         z = 0.5*((1 - r3)*fRmin[2] + (r3 + 1)*fRmax[2]);
         xx = fTnorm[0]*x + fTnorm[1]*y + fTnorm[2]*z + fTnorm[3];
         yy = fTnorm[4]*x + fTnorm[5]*y + fTnorm[6]*z + fTnorm[7];
         if (TMath::Abs(xx - xgraf[1]) <= 1e-4) {
            if (ygraf[1] >= yy) ygraf[1] = yy;
            if (ygraf[2] <= yy) ygraf[2] = yy;
         }
         if (TMath::Abs(xx - xgraf[5]) <= 1e-4) {
            if (ygraf[5] >= yy) ygraf[5] = yy;
            if (ygraf[4] <= yy) ygraf[4] = yy;
         }
         if (TMath::Abs(yy - ygraf[0]) <= 1e-4) xgraf[0] = xx;
         if (TMath::Abs(yy - ygraf[3]) <= 1e-4) xgraf[3] = xx;
         r1 = -r1;
         if (i % 2 == 0) r2 = -r2;
         if (i >= 4)     r3 = 1;
      }
      gPad->PaintFillArea(6, xgraf, ygraf);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Store axis coordinates in the NDC system.

void  TView3D::SetAxisNDC(const Double_t *x1, const Double_t *x2, const Double_t *y1, const Double_t *y2, const Double_t *z1, const Double_t *z2)
{
   for (Int_t i=0;i<3;i++) {
      fX1[i] = x1[i];
      fX2[i] = x2[i];
      fY1[i] = y1[i];
      fY2[i] = y2[i];
      fZ1[i] = z1[i];
      fZ2[i] = z2[i];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set default viewing window.

void TView3D::SetDefaultWindow()
{
   if (!gPad) return;
   Double_t screen_factor = 1.;
   Double_t du, dv;
   Double_t extent = GetExtent();
   fDview = 3*extent;
   fDproj = 0.5*extent;

   // width in pixels
   fUpix = gPad->GetWw()*gPad->GetAbsWNDC();

   // height in pixels
   fVpix = gPad->GetWh()*gPad->GetAbsHNDC();
   du = 0.5*screen_factor*fDproj;
   dv = du*fVpix/fUpix;   // keep aspect ratio
   SetWindow(0, 0, du, dv);
}

////////////////////////////////////////////////////////////////////////////////
/// This is a function which creates default outline.
///
/// ~~~ {.cpp}
///      x = fRmin[0]        X = fRmax[0]
///      y = fRmin[1]        Y = fRmax[1]
///      z = fRmin[2]        Z = fRmax[2]
///
///               (x,Y,Z) +---------+ (X,Y,Z)
///                      /         /|
///                     /         / |
///                    /         /  |
///           (x,y,Z) +---------+   |
///                   |         |   + (X,Y,z)
///                   |         |  /
///                   |         | /
///                   |         |/
///                   +---------+
///                (x,y,z)   (X,y,z)
/// ~~~

void TView3D::SetOutlineToCube()
{
   if (!fOutline) {
      fDefaultOutline = kTRUE;
      fOutline = new TList();
   }
   DrawOutlineCube((TList*)fOutline,fRmin,fRmax);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the parallel option (default).

void TView3D::SetParallel()
{
   if (!IsPerspective()) return;
   SetBit(kPerspective, kFALSE);
   Int_t irep;
   ResetView(fLongitude, fLatitude, fPsi, irep);
}

////////////////////////////////////////////////////////////////////////////////
/// Set perspective option.

void TView3D::SetPerspective()
{
   if (IsPerspective()) return;
   SetBit(kPerspective, kTRUE);
   Int_t irep;
   SetDefaultWindow();
   ResetView(fLongitude, fLatitude, fPsi, irep);
}

////////////////////////////////////////////////////////////////////////////////
/// Set Range function.

void TView3D::SetRange(const Double_t *min, const Double_t *max)
{
   Int_t irep;
   for (Int_t i = 0; i < 3; fRmax[i] = max[i], fRmin[i] = min[i], i++) { }
   if (IsPerspective()) SetDefaultWindow();
   ResetView(fLongitude, fLatitude, fPsi, irep);
   if(irep < 0)
      Error("SetRange", "problem setting view");
   if(fDefaultOutline) SetOutlineToCube();
}

////////////////////////////////////////////////////////////////////////////////
/// Set 3-D View range.
///
/// Input:
///  - x0, y0, z0 are minimum coordinates
///  - x1, y1, z1 are maximum coordinates
///
///  - flag values are:
///    - 0 (set always) <- default
///    - 1 (shrink view)
///    - 2 (expand view)

void TView3D::SetRange(Double_t x0, Double_t y0, Double_t z0, Double_t x1, Double_t y1, Double_t z1, Int_t flag)
{
   Double_t rmax[3], rmin[3];

   switch (flag) {
      case 2:                     // expand view
         GetRange(rmin, rmax);
         rmin[0] = x0 < rmin[0] ? x0 : rmin[0];
         rmin[1] = y0 < rmin[1] ? y0 : rmin[1];
         rmin[2] = z0 < rmin[2] ? z0 : rmin[2];
         rmax[0] = x1 > rmax[0] ? x1 : rmax[0];
         rmax[1] = y1 > rmax[1] ? y1 : rmax[1];
         rmax[2] = z1 > rmax[2] ? z1 : rmax[2];
         break;

      case 1:                     // shrink view
         GetRange(rmin, rmax);
         rmin[0] = x0 > rmin[0] ? x0 : rmin[0];
         rmin[1] = y0 > rmin[1] ? y0 : rmin[1];
         rmin[2] = z0 > rmin[2] ? z0 : rmin[2];
         rmax[0] = x1 < rmax[0] ? x1 : rmax[0];
         rmax[1] = y1 < rmax[1] ? y1 : rmax[1];
         rmax[2] = z1 < rmax[2] ? z1 : rmax[2];
         break;

      default:
         rmin[0] = x0; rmax[0] = x1;
         rmin[1] = y0; rmax[1] = y1;
         rmin[2] = z0; rmax[2] = z1;
   }
   SetRange(rmin, rmax);
}

////////////////////////////////////////////////////////////////////////////////
/// Set viewing window.

void TView3D::SetWindow(Double_t u0, Double_t v0, Double_t du, Double_t dv)
{
   fUVcoord[0] = u0;
   fUVcoord[1] = v0;
   fUVcoord[2] = du;
   fUVcoord[3] = dv;
}

////////////////////////////////////////////////////////////////////////////////
/// Set view parameters.

void TView3D::SetView(Double_t longitude, Double_t latitude, Double_t psi, Int_t &irep)
{
   ResetView(longitude, latitude, psi, irep);
}

////////////////////////////////////////////////////////////////////////////////
/// Recompute window for perspective view.

void TView3D::ResizePad()
{
   if (!IsPerspective()) return;
   Double_t upix = fUpix;
   Double_t vpix = fVpix;

   // width in pixels
   fUpix = gPad->GetWw()*gPad->GetAbsWNDC();

   // height in pixels
   fVpix = gPad->GetWh()*gPad->GetAbsHNDC();
   Double_t u0 = fUVcoord[0]*fUpix/upix;
   Double_t v0 = fUVcoord[1]*fVpix/vpix;
   Double_t du = fUVcoord[2]*fUpix/upix;
   Double_t dv = fUVcoord[3]*fVpix/vpix;
   SetWindow(u0, v0, du, dv);
   DefinePerspectiveView();
}

////////////////////////////////////////////////////////////////////////////////
/// Set view direction (in spherical coordinates).
///
/// Input
///  - PHI     - longitude
///  - THETA   - latitude (angle between +Z and view direction)
///  - PSI     - rotation in screen plane
///
/// Output:
///  - IREP   - reply (-1 if error in min-max)
///
/// Errors: error in min-max scope

void TView3D::ResetView(Double_t longitude, Double_t latitude, Double_t psi, Int_t &irep)
{
   Double_t scale[3],  centre[3];
   Double_t c1, c2, c3, s1, s2, s3;

   //*-*-        F I N D   C E N T E R   O F   S C O P E   A N D
   //*-*-        S C A L E   F A C T O R S
   FindScope(scale, centre, irep);
   if (irep < 0) {
      Error("ResetView", "Error in min-max scope");
      return;
   }

   //*-*-        S E T   T R A N S F O R M A T I O N   M A T R I C E S
   fLongitude = longitude;
   fPsi       = psi;
   fLatitude  = latitude;

   if (IsPerspective()) {
      DefinePerspectiveView();
      return;
   }

   c1 = TMath::Cos(longitude*kRad);
   s1 = TMath::Sin(longitude*kRad);
   c2 = TMath::Cos(latitude*kRad);
   s2 = TMath::Sin(latitude*kRad);
   c3 = TMath::Cos(psi*kRad);
   s3 = TMath::Sin(psi*kRad);
   DefineViewDirection(scale, centre, c1, s1, c2, s2, c3, s3, fTnorm, fTback);
   c3 = 1;
   s3 = 0;
   DefineViewDirection(scale, centre, c1, s1, c2, s2, c3, s3, fTN, fTB);
}

////////////////////////////////////////////////////////////////////////////////
/// Transfer point from world to normalized coordinates.
///
/// Input:
///  - PW(3) - point in world coordinate system
///  - PN(3) - point in normalized coordinate system

void TView3D::WCtoNDC(const Float_t *pw, Float_t *pn)
{
   Float_t x = pw[0], y = pw[1], z = pw[2];

   // perspective view
   if (IsPerspective()) {
      for (Int_t i=0; i<3; i++) {
         pn[i] = fTnorm[i]*x + fTnorm[i+4]*y + fTnorm[i+8]*z + fTnorm[i+12];
      }
      if (pn[2]>0) {
         pn[0] /= pn[2];
         pn[1] /= pn[2];
      } else {
         pn[0] *= 1000.;
         pn[1] *= 1000.;
      }
      return;
   }

   // parallel view
   pn[0] = fTnorm[0]*x + fTnorm[1]*y + fTnorm[2]*z  + fTnorm[3];
   pn[1] = fTnorm[4]*x + fTnorm[5]*y + fTnorm[6]*z  + fTnorm[7];
   pn[2] = fTnorm[8]*x + fTnorm[9]*y + fTnorm[10]*z + fTnorm[11];
}

////////////////////////////////////////////////////////////////////////////////
/// Transfer point from world to normalized coordinates.
///
/// Input:
///  - PW(3) - point in world coordinate system
///  - PN(3) - point in normalized coordinate system

void TView3D::WCtoNDC(const Double_t *pw, Double_t *pn)
{
   Double_t x = pw[0], y = pw[1], z = pw[2];

   // perspective view
   if (IsPerspective()) {
      for (Int_t i=0; i<3; i++) {
         pn[i] = fTnorm[i]*x + fTnorm[i+4]*y + fTnorm[i+8]*z + fTnorm[i+12];
      }
      if (pn[2]>0) {
         pn[0] /= pn[2];
         pn[1] /= pn[2];
      } else {
         pn[0] *= 1000.;
         pn[1] *= 1000.;
      }
      return;
   }

   // parallel view
   pn[0] = fTnorm[0]*x + fTnorm[1]*y + fTnorm[2]*z  + fTnorm[3];
   pn[1] = fTnorm[4]*x + fTnorm[5]*y + fTnorm[6]*z  + fTnorm[7];
   pn[2] = fTnorm[8]*x + fTnorm[9]*y + fTnorm[10]*z + fTnorm[11];
}

////////////////////////////////////////////////////////////////////////////////
/// Force the current pad to be updated.

void TView3D::AdjustPad(TVirtualPad *pad)
{
   TVirtualPad *thisPad = pad;
   if (!thisPad) thisPad = gPad;
   if (thisPad) {
#ifdef R__HAS_COCOA
      thisPad->AbsCoordinates(kFALSE);
#endif
      thisPad->Modified();
      thisPad->Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// API to rotate view and adjust the pad provided it the current one.

void TView3D::RotateView(Double_t phi, Double_t theta, TVirtualPad *pad)
{
   Int_t iret;
   Double_t p = phi;
   Double_t t = theta;
   SetView(p, t, 0, iret);

   // Adjust current pad too
   TVirtualPad *thisPad = pad;
   if (!thisPad) thisPad = gPad;
   if (thisPad) {
      thisPad->SetPhi(-90-p);
      thisPad->SetTheta(90-t);
      thisPad->Modified();
      thisPad->Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set to side view.

void TView3D::SideView(TVirtualPad *pad)
{
   RotateView(0,90.0,pad);
}

////////////////////////////////////////////////////////////////////////////////
/// Set to front view.

void TView3D::FrontView(TVirtualPad *pad)
{
   RotateView(270.0,90.0,pad);
}

////////////////////////////////////////////////////////////////////////////////
/// Set to top view.

void TView3D::TopView(TVirtualPad *pad)
{
   RotateView(270.0,0.0,pad);
}

////////////////////////////////////////////////////////////////////////////////
/// Turn on /off 3D axis.

void TView3D::ToggleRulers(TVirtualPad *pad)
{
   TAxis3D::ToggleRulers(pad);
}

////////////////////////////////////////////////////////////////////////////////
/// Turn on /off the interactive option to
///  Zoom / Move / Change attributes of 3D axis correspond this view.

void TView3D::ToggleZoom(TVirtualPad *pad)
{
   TAxis3D::ToggleZoom(pad);
}

////////////////////////////////////////////////////////////////////////////////
/// Adjust all sides of view in respect of the biggest one.

void TView3D::AdjustScales(TVirtualPad *pad)
{
   Double_t min[3],max[3];
   GetRange(min,max);
   int i;
   Double_t maxSide = 0;
   // Find the largest side
   for (i=0;i<3; i++) maxSide = TMath::Max(maxSide,max[i]-min[i]);
   //Adjust scales:
   for (i=0;i<3; i++) max[i] += maxSide - (max[i]-min[i]);
   SetRange(min,max);

   AdjustPad(pad);
}

////////////////////////////////////////////////////////////////////////////////
/// Move view into the center of the scene.

void TView3D::Centered3DImages(TVirtualPad *pad)
{
   Double_t min[3],max[3];
   GetRange(min,max);
   int i;
   for (i=0;i<3; i++) {
      if (max[i] > 0) min[i] = -max[i];
      else            max[i] = -min[i];
   }
   SetRange(min,max);
   AdjustPad(pad);
}

////////////////////////////////////////////////////////////////////////////////
/// unZOOM this view.

void TView3D::UnzoomView(TVirtualPad *pad,Double_t unZoomFactor )
{
   if (TMath::Abs(unZoomFactor) < 0.001) return;
   ZoomView(pad,1./unZoomFactor);
}

////////////////////////////////////////////////////////////////////////////////
/// ZOOM this view.

void TView3D::ZoomView(TVirtualPad *pad,Double_t zoomFactor)
{
   if (TMath::Abs(zoomFactor) < 0.001) return;
   Double_t min[3],max[3];
   GetRange(min,max);
   int i;
   for (i=0;i<3; i++) {
      // Find center
      Double_t c = (max[i]+min[i])/2;
      // Find a new size
      Double_t s = (max[i]-min[i])/(2*zoomFactor);
      // Set a new size
      max[i] = c + s;
      min[i] = c - s;
   }
   SetRange(min,max);
   AdjustPad(pad);
}

////////////////////////////////////////////////////////////////////////////////
/// Move focus to a different box position and extent in nsteps. Perform
/// rotation with dlat,dlong,dpsi at each step.

void TView3D::MoveFocus(Double_t *cov, Double_t dx, Double_t dy, Double_t dz, Int_t nsteps,
                      Double_t dlong, Double_t dlat, Double_t dpsi)
{
   if (!IsPerspective()) return;
   if (nsteps<1) return;
   Double_t fc = 1./Double_t(nsteps);
   Double_t oc[3], od[3], dir[3];
   dir[0] = 0;
   dir[1] = 0;
   dir[2] = 1.;
   Int_t i, j;
   for (i=0; i<3; i++) {
      oc[i] = 0.5*(fRmin[i]+fRmax[i]);
      od[i] = 0.5*(fRmax[i]-fRmin[i]);
   }
   Double_t dox = cov[0]-oc[0];
   Double_t doy = cov[1]-oc[1];
   Double_t doz = cov[2]-oc[2];

   Double_t dd = TMath::Sqrt(dox*dox+doy*doy+doz*doz);
   if (dd!=0) {;
   dir[0] = dox/dd;
   dir[1] = doy/dd;
   dir[2] = doz/dd;
   }
   dd *= fc;
   dox = fc*(dx-od[0]);
   doy = fc*(dy-od[1]);
   doz = fc*(dz-od[2]);
   for (i=0; i<nsteps; i++) {
      oc[0] += dd*dir[0];
      oc[1] += dd*dir[1];
      oc[2] += dd*dir[2];
      od[0]  += dox;
      od[1]  += doy;
      od[2]  += doz;
      for (j=0; j<3; j++) {
         fRmin[j] = oc[j]-od[j];
         fRmax[j] = oc[j]+od[j];
      }
      SetDefaultWindow();
      fLatitude += dlat;
      fLongitude += dlong;
      fPsi += dpsi;
      DefinePerspectiveView();
      if (gPad) {
         gPad->Modified();
         gPad->Update();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
///  - 'a' increase scale factor (clip cube borders)
///  - 's' decrease scale factor (clip cube borders)

void TView3D::MoveViewCommand(Char_t option, Int_t count)
{
   if (count <= 0) count = 1;
   switch (option) {
      case '+':
         ZoomView();
         break;
      case '-':
         UnzoomView();
         break;
      case 's':
      case 'S':
         UnzoomView();
         break;
      case 'a':
      case 'A':
         ZoomView();
         break;
      case 'l':
      case 'L':
      case 'h':
      case 'H':
      case 'u':
      case 'U':
      case 'i':
      case 'I':
         MoveWindow(option);
         break;
      case 'j':
      case 'J':
         ZoomIn();
         break;
      case 'k':
      case 'K':
         ZoomOut();
         break;
      default:
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Move view window :
///  - l,L - left
///  - h,H - right
///  - u,U - down
///  - i,I - up

void TView3D::MoveWindow(Char_t option)
{
   if (!IsPerspective()) return;
   Double_t shiftu = 0.1*fUVcoord[2];
   Double_t shiftv = 0.1*fUVcoord[3];
   switch (option) {
      case 'l':
      case 'L':
         fUVcoord[0] += shiftu;
         break;
      case 'h':
      case 'H':
         fUVcoord[0] -= shiftu;
         break;
      case 'u':
      case 'U':
         fUVcoord[1] += shiftv;
         break;
      case 'i':
      case 'I':
         fUVcoord[1] -= shiftv;
         break;
      default:
         return;
   }
   DefinePerspectiveView();
   if (gPad) {
      gPad->Modified();
      gPad->Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Zoom in.

void TView3D::ZoomIn()
{
   if (!IsPerspective()) return;
   Double_t extent = GetExtent();
   Double_t fc = 0.1;
   if (fDview<extent) {
      fDview -= fc*extent;
   } else {
      fDview /= 1.25;
   }
   DefinePerspectiveView();
   if (gPad) {
      gPad->Modified();
      gPad->Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Zoom out.

void TView3D::ZoomOut()
{
   if (!IsPerspective()) return;
   Double_t extent = GetExtent();
   Double_t fc = 0.1;
   if (fDview<extent) {
      fDview += fc*extent;
   } else {
      fDview *= 1.25;
   }
   DefinePerspectiveView();
   if (gPad) {
      gPad->Modified();
      gPad->Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TView3D.

void TView3D::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         R__b.ReadClassBuffer(TView3D::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      //unfortunately we forgot to increment the TView3D version number
      //when the class was upgraded to double precision in version 2.25.
      //we are forced to use the file version number to recognize old files.
      if (R__b.GetParent() && R__b.GetVersionOwner() < 22500) { //old version in single precision
         TObject::Streamer(R__b);
         TAttLine::Streamer(R__b);
         Float_t single, sa[12];
         Int_t i;
         R__b >> fSystem;
         R__b >> single; fLatitude = single;
         R__b >> single; fLongitude = single;
         R__b >> single; fPsi = single;
         R__b.ReadStaticArray(sa);   for (i=0;i<12;i++) fTN[i] = sa[i];
         R__b.ReadStaticArray(sa);   for (i=0;i<12;i++) fTB[i] = sa[i];
         R__b.ReadStaticArray(sa);   for (i=0;i<3;i++)  fRmax[i] = sa[i];
         R__b.ReadStaticArray(sa);   for (i=0;i<3;i++)  fRmin[i] = sa[i];
         R__b.ReadStaticArray(sa);   for (i=0;i<12;i++) fTnorm[i] = sa[i];
         R__b.ReadStaticArray(sa);   for (i=0;i<12;i++) fTback[i] = sa[i];
         R__b.ReadStaticArray(sa);   for (i=0;i<3;i++)  fX1[i] = sa[i];
         R__b.ReadStaticArray(sa);   for (i=0;i<3;i++)  fX2[i] = sa[i];
         R__b.ReadStaticArray(sa);   for (i=0;i<3;i++)  fY1[i] = sa[i];
         R__b.ReadStaticArray(sa);   for (i=0;i<3;i++)  fY2[i] = sa[i];
         R__b.ReadStaticArray(sa);   for (i=0;i<3;i++)  fZ1[i] = sa[i];
         R__b.ReadStaticArray(sa);   for (i=0;i<3;i++)  fZ2[i] = sa[i];
         R__b >> fOutline;
         R__b >> fDefaultOutline;
         R__b >> fAutoRange;
      } else {
         TObject::Streamer(R__b);
         TAttLine::Streamer(R__b);
         R__b >> fLatitude;
         R__b >> fLongitude;
         R__b >> fPsi;
         R__b.ReadStaticArray(fTN);
         R__b.ReadStaticArray(fTB);
         R__b.ReadStaticArray(fRmax);
         R__b.ReadStaticArray(fRmin);
         R__b.ReadStaticArray(fTnorm);
         R__b.ReadStaticArray(fTback);
         R__b.ReadStaticArray(fX1);
         R__b.ReadStaticArray(fX2);
         R__b.ReadStaticArray(fY1);
         R__b.ReadStaticArray(fY2);
         R__b.ReadStaticArray(fZ1);
         R__b.ReadStaticArray(fZ2);
         R__b >> fSystem;
         R__b >> fOutline;
         R__b >> fDefaultOutline;
         R__b >> fAutoRange;
      }
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TView3D::Class(),this);
   }
}

//      Shortcuts for menus
void TView3D::Centered(){Centered3DImages();}
void TView3D::Front()   {FrontView();}
void TView3D::ShowAxis(){ToggleRulers(); }
void TView3D::Side()    {SideView();}
void TView3D::Top()     {TopView();}
void TView3D::ZoomMove(){ToggleZoom();}
void TView3D::Zoom()    {ZoomView();}
void TView3D::UnZoom()  {UnzoomView();}


