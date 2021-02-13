// @(#)root/gl:$Id$
// Author:  Olivier Couet  17/04/2007

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"

#include "TGLIncludes.h"
#include "TGLUtil.h"
#include "TGLAxis.h"
#include "TGLText.h"
#include "TColor.h"
#include "TString.h"
#include "TMath.h"
#include "THLimitsFinder.h"

/** \class TGLAxis
\ingroup opengl
GL Axis.

To draw a 3D axis in a GL window. The labels are drawn using FTGL.
*/

ClassImp(TGLAxis);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGLAxis::TGLAxis(): TAttLine(1,1,1), TAttText(20,0.,1,42,0.04)
{
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// Default initialization.

void TGLAxis::Init()
{
   fNDiv = fNDiv1 = fNDiv2 = fNDiv3 = 0;
   fNTicks1 = fNTicks2 = 0;
   fTicks1          = 0;
   fTicks2          = 0;
   fLabels          = 0;
   fText            = 0;
   fAngle1          = 90.;
   fAngle2          = 0.;
   fAngle3          = 0.;
   fAxisLength      = 0.;
   fWmin = fWmax    = 0.;
   fTickMarksLength = 0.04; // % of fAxisLength
   fTickMarksOrientation = 2; // can be 0, 1, 2, or 3
   fLabelsOffset    = 0.09; // % of fAxisLength
   fLabelsSize      = 0.06; // % of fAxisLength
   fGridLength      = 0.;   // 0 means no grid
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGLAxis::~TGLAxis()
{
   if (fTicks1) delete [] fTicks1;
   if (fTicks2) delete [] fTicks2;
   if (fLabels) delete [] fLabels;
   if (fText)   delete fText;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint GL Axis.
///
///  - p1, p2     : Axis position in the 3D space.
///  - wmin, wmax : Minimum and maximum values along the axis. wmin < wmax.
///  - ndiv       : Number of axis divisions. It is an integer in the form
///                 "ttsspp" where "tt" is the number of tertiary divisions,
///                 "ss" is the number of secondary divisions and "pp" the
///                 number of primary divisions.
///  - opt        : Options.
///                 "N" - By default the number of divisions is optimized to
///                 get a nice labeling. When option "N" is given, the
///                 number of divisions is not optimized.

void TGLAxis::PaintGLAxis(const Double_t p1[3], const Double_t p2[3],
                          Double_t wmin,  Double_t wmax, Int_t ndiv,
                          Option_t *opt)
{
   fNDiv = ndiv;
   if (wmax<=wmin) {
      fWmax = wmin;
      fWmin = wmax;
   } else {
      fWmax = wmax;
      fWmin = wmin;
   }

   // Compute the axis length.
   Double_t x1 = p1[0], y1 = p1[1], z1 = p1[2];
   Double_t x2 = p2[0], y2 = p2[1], z2 = p2[2];
   fAxisLength = TMath::Sqrt((x2-x1)*(x2-x1)+
                             (y2-y1)*(y2-y1)+
                             (z2-z1)*(z2-z1));

   TicksPositions(opt);

   DoLabels();

   // Paint using GL
   glPushMatrix();

   // Translation
   glTranslatef(x1, y1, z1);

   // Rotation in the Z plane
   Double_t phi=0;
   Double_t normal[3];
   normal[0] = 0;
   normal[1] = 1;
   normal[2] = 0;
   if (z1!=z2) {
      if (y2==y1 && x2==x1) {
         if (z2<z1) phi = 90;
         else       phi = 270;
      } else {
         Double_t p3[3];
         p3[0] = p2[0]; p3[1] = p2[1]; p3[2] = 0;
         TMath::Normal2Plane(p1, p2, p3, normal);
         phi = TMath::ACos(TMath::Abs(z2-z1)/fAxisLength);
         phi = -(90-(180/TMath::Pi())*phi);
      }
      glRotatef(phi, normal[0], normal[1], normal[2]);
   }

   // Rotation in the XY plane
   Double_t theta = 0;
   if (y2!=y1) {
      if ((x2-x1)>0) {
         theta = TMath::ATan((y2-y1)/(x2-x1));
         theta = (180/TMath::Pi())*theta;
      } else if ((x2-x1)<0) {
         theta = TMath::ATan((y2-y1)/(x2-x1));
         theta = 180+(180/TMath::Pi())*theta;
      } else {
         if (y2>y1) theta = 90;
         else       theta = 270;
      }
   } else {
      if (x2<x1) theta = 180;
   }
   glRotatef(theta, 0., 0., 1.);

   PaintGLAxisBody();

   PaintGLAxisTickMarks();

   PaintGLAxisLabels();

   glPopMatrix();
}

////////////////////////////////////////////////////////////////////////////////
/// Paint horizontal axis body at position (0,0,0)

void TGLAxis::PaintGLAxisBody()
{
   TColor *col;
   Float_t red = 1.f, green = 1.f, blue = 1.f;
   if ((col = gROOT->GetColor(GetLineColor())))
      col->GetRGB(red, green, blue);
   glColor3d(red, green, blue);
   TGLUtil::LineWidth(GetLineWidth());
   glBegin(GL_LINES);
   glVertex3d(0., 0., 0.);
   glVertex3d(fAxisLength, 0., 0.);
   glEnd();
}

////////////////////////////////////////////////////////////////////////////////
/// Paint axis tick marks.

void TGLAxis::PaintGLAxisTickMarks()
{
   Int_t i;

   // Ticks marks orientation;
   Double_t yo=0, zo=0;
   switch (fTickMarksOrientation) {
      case 0:
         yo = 0;
         zo = 1;
      break;
      case 1:
         yo = -1;
         zo = 0;
      break;
      case 2:
         yo = 0;
         zo = -1;
      break;
      case 3:
         yo = 1;
         zo = 0;
      break;
   }

   // Paint level 1 tick marks.
   if (fTicks1) {
      // Paint the tick marks, if needed.
      if (fTickMarksLength) {
         Double_t tl = fTickMarksLength*fAxisLength;
         glBegin(GL_LINES);
         for (i=0; i<fNTicks1 ; i++) {
            glVertex3f( fTicks1[i], 0, 0);
            glVertex3f( fTicks1[i], yo*tl, zo*tl);
         }
         glEnd();
      }

      // Paint the grid, if needed, on level 1 tick marks.
      if (fGridLength) {
         const UShort_t stipple = 0x8888;
         glLineStipple(1, stipple);
         glEnable(GL_LINE_STIPPLE);
         glBegin(GL_LINES);
         for (i=0; i<fNTicks1; i++) {
            glVertex3f( fTicks1[i], 0, 0);
            glVertex3f( fTicks1[i], -yo*fGridLength, -zo*fGridLength);
         }
         glEnd();
         glDisable(GL_LINE_STIPPLE);
      }
   }

   // Paint level 2 tick marks.
   if (fTicks2) {
      if (fTickMarksLength) {
         Double_t tl = 0.5*fTickMarksLength*fAxisLength;
         glBegin(GL_LINES);
         for (i=0; i<fNTicks2; i++) {
            glVertex3f( fTicks2[i], 0, 0);
            glVertex3f( fTicks2[i], yo*tl, zo*tl);
         }
         glEnd();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint axis labels on the main tick marks.

void TGLAxis::PaintGLAxisLabels()
{
   if (!fLabelsSize) return;

   Double_t x=0,y=0,z=0;
   Int_t i;

   if (!fText) {
      fText = new TGLText();
      fText->SetTextColor(GetTextColor());
      fText->SetGLTextFont(GetTextFont());
      fText->SetTextSize(fLabelsSize*fAxisLength);
      fText->SetTextAlign(GetTextAlign());
   }
   fText->SetGLTextAngles(fAngle1, fAngle2, fAngle3);

   switch (fTickMarksOrientation) {
      case 0:
         y = 0;
         z = fLabelsOffset*fAxisLength;
      break;
      case 1:
         y = -fLabelsOffset*fAxisLength;
         z = 0;
      break;
      case 2:
         y = 0;
         z = -fLabelsOffset*fAxisLength;
      break;
      case 3:
         y = fLabelsOffset*fAxisLength;
         z = 0;
      break;
   }
   for (i=0; i<fNDiv1+1 ; i++) {
      x = fTicks1[i];
      fText->PaintGLText(x,y,z,fLabels[i]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute ticks positions.

void TGLAxis::TicksPositions(Option_t *opt)
{
   Bool_t optionNoopt /* , optionLog */;

   if (strchr(opt,'N')) optionNoopt = kTRUE;  else optionNoopt = kFALSE;
   // if (strchr(opt,'G')) optionLog   = kTRUE;  else optionLog   = kFALSE;

   // Determine number of tick marks 1, 2 and 3.
   fNDiv3 = fNDiv/10000;
   fNDiv2 = (fNDiv-10000*fNDiv3)/100;
   fNDiv1 = fNDiv%100;

   // Clean the tick marks arrays if they exist.
   if (fTicks1) {
      delete [] fTicks1;
      fTicks1 = 0;
   }
   if (fTicks2) {
      delete [] fTicks2;
      fTicks2 = 0;
   }

   // Compute the tick marks positions according to the options.
   if (optionNoopt) {
      TicksPositionsNoOpt();
   } else {
      TicksPositionsOpt();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute ticks positions. Linear and not optimized.

void TGLAxis::TicksPositionsNoOpt()
{
   Int_t i, j, k;
   Double_t step1 = fAxisLength/(fNDiv1);

   fNTicks1 = fNDiv1+1;
   fTicks1  = new Double_t[fNTicks1];

   // Level 1 tick marks.
   for (i=0; i<fNTicks1; i++) fTicks1[i] = i*step1;

   // Level 2 tick marks.
   if (fNDiv2) {
      Double_t t2;
      Double_t step2 = step1/fNDiv2;
      fNTicks2       = fNDiv1*(fNDiv2-1);
      fTicks2        = new Double_t[fNTicks2];
      k = 0;
      for (i=0; i<fNTicks1-1; i++) {
         t2 = fTicks1[i]+step2;
         for (j=0; j<fNDiv2-1; j++) {
            fTicks2[k] = t2;
            k++;
            t2 = t2+step2;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute ticks positions. Linear and optimized.

void TGLAxis::TicksPositionsOpt()
{
   Int_t i, j, k, nDivOpt;
   Double_t step1=0, step2=0, wmin2=0, wmax2=0;
   Double_t wmin = fWmin;
   Double_t wmax = fWmax;

   // Level 1 tick marks.
   THLimitsFinder::Optimize(wmin,  wmax, fNDiv1,
                            fWmin, fWmax, nDivOpt,
                            step1, "");
   fNDiv1   = nDivOpt;
   fNTicks1 = fNDiv1+1;
   fTicks1  = new Double_t[fNTicks1];

   Double_t r = fAxisLength/(wmax-wmin);
   Double_t w = fWmin;
   i = 0;
   while (w<=fWmax) {
      fTicks1[i] = r*(w-wmin);
      i++;
      w = w+step1;
   }

   // Level 2 tick marks.
   if (fNDiv2) {
      Double_t t2;
      THLimitsFinder::Optimize(fWmin, fWmin+step1, fNDiv2,
                               wmin2, wmax2, nDivOpt,
                               step2, "");
      fNDiv2       = nDivOpt;
      step2        = TMath::Abs((fTicks1[1]-fTicks1[0])/fNDiv2);
      Int_t nTickl = (Int_t)(fTicks1[0]/step2);
      Int_t nTickr = (Int_t)((fAxisLength-fTicks1[fNTicks1-1])/step2);
      fNTicks2     = fNDiv1*(fNDiv2-1)+nTickl+nTickr;
      fTicks2      = new Double_t[fNTicks2];
      k = 0;
      for (i=0; i<fNTicks1-1; i++) {
         t2 = fTicks1[i]+step2;
         for (j=0; j<fNDiv2-1; j++) {
            fTicks2[k] = t2;
            k++;
            t2 = t2+step2;
         }
      }
      if (nTickl) {
         t2 = fTicks1[0]-step2;
         for (i=0; i<nTickl; i++) {
            fTicks2[k] = t2;
            k++;
            t2 = t2-step2;
         }
      }
      if (nTickr) {
         t2 = fTicks1[fNTicks1-1]+step2;
         for (i=0; i<nTickr; i++) {
            fTicks2[k] = t2;
            k++;
            t2 = t2+step2;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Do labels.

void TGLAxis::DoLabels()
{
   if (fLabels) delete [] fLabels;
   fLabels = new TString[fNTicks1];
   Int_t i;

   // Make regular labels between fWmin and fWmax.
   Double_t dw = (fWmax-fWmin)/(fNDiv1);
   for (i=0; i<fNTicks1; i++) {
      fLabels[i] = Form("%g",fWmin+i*dw);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set labels' angles.

void TGLAxis::SetLabelsAngles(Double_t a1, Double_t a2, Double_t a3)
{
   fAngle1 = a1;
   fAngle2 = a2;
   fAngle3 = a3;
}
