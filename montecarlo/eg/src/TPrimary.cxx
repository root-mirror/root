// @(#)root/eg:$Id$
// Author: Ola Nordmann   21/09/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class  TPrimary
    \ingroup eg

Old version of a  dynamic particle class created by event generators.

This class is now obsolete. Use TParticle instead.
*/

#include "TObject.h"
#include "Rtypes.h"
#include "TString.h"
#include "TAttParticle.h"
#include "TPrimary.h"
#include "TView.h"
#include "TMath.h"
#include "TVirtualPad.h"
#include "TPolyLine3D.h"
#include "snprintf.h"

ClassImp(TPrimary);

////////////////////////////////////////////////////////////////////////////////
///
///  Primary vertex particle default constructor
///

TPrimary::TPrimary()
{
   //do nothing
   fPart         = 0;
   fFirstMother  = 0;
   fSecondMother = 0;
   fGeneration   = 0;
   fPx           = 0;
   fPy           = 0;
   fPz           = 0;
   fEtot         = 0;
   fVx           = 0;
   fVy           = 0;
   fVz           = 0;
   fTime         = 0;
   fTimeEnd      = 0;
   fType         = "";

}

////////////////////////////////////////////////////////////////////////////////
///
///  TPrimary vertex particle normal constructor
///

TPrimary::TPrimary(Int_t part, Int_t first, Int_t second, Int_t gener,
                   Double_t px, Double_t py, Double_t pz,
                   Double_t etot, Double_t vx, Double_t vy, Double_t vz,
                   Double_t time, Double_t timend, const char *type)
{
   fPart         = part;
   fFirstMother  = first;
   fSecondMother = second;
   fGeneration   = gener;
   fPx           = px;
   fPy           = py;
   fPz           = pz;
   fEtot         = etot;
   fVx           = vx;
   fVy           = vy;
   fVz           = vz;
   fTime         = time;
   fTimeEnd      = timend;
   fType         = type;
}

////////////////////////////////////////////////////////////////////////////////
///
///   Primaray vertex particle default destructor
///

TPrimary::~TPrimary()
{
   //do nothing
}


////////////////////////////////////////////////////////////////////////////////
/// Compute distance from point px,py to a primary track
///
/// Compute the closest distance of approach from point px,py to each segment
/// of a track.
/// The distance is computed in pixels units.
///

Int_t TPrimary::DistancetoPrimitive(Int_t px, Int_t py)
{
   const Int_t big = 9999;
   Float_t xv[3], xe[3], xndc[3];
   Float_t rmin[3], rmax[3];
   TView *view = gPad->GetView();
   if(!view) return big;

   // compute first and last point in pad coordinates
   Float_t pmom = TMath::Sqrt(fPx*fPx+fPy*fPy+fPz*fPz);
   if (pmom == 0) return big;
   view->GetRange(rmin,rmax);
   Float_t rbox = rmax[2];
   xv[0] = fVx;
   xv[1] = fVy;
   xv[2] = fVz;
   xe[0] = fVx+rbox*fPx/pmom;
   xe[1] = fVy+rbox*fPy/pmom;
   xe[2] = fVz+rbox*fPz/pmom;
   view->WCtoNDC(xv, xndc);
   Float_t x1 = xndc[0];
   Float_t y1 = xndc[1];
   view->WCtoNDC(xe, xndc);
   Float_t x2 = xndc[0];
   Float_t y2 = xndc[1];

   return DistancetoLine(px,py,x1,y1,x2,y2);
}


////////////////////////////////////////////////////////////////////////////////
/// Execute action corresponding to one event
///

void TPrimary::ExecuteEvent(Int_t, Int_t, Int_t)
{
   gPad->SetCursor(kPointer);
}

////////////////////////////////////////////////////////////////////////////////
/// Return name of primary particle

const char *TPrimary::GetName() const
{
   static char def[4] = "XXX";
   const TAttParticle *ap = GetParticle();
   if (ap) return ap->GetName();
   else    return def;
}

////////////////////////////////////////////////////////////////////////////////
///
/// Returning a pointer to the particle attributes
///

const TAttParticle *TPrimary::GetParticle() const
{
   if (!TAttParticle::fgList) TAttParticle::DefinePDG();
   return TAttParticle::GetParticle(fPart);
}

////////////////////////////////////////////////////////////////////////////////
/// Return title of primary particle

const char *TPrimary::GetTitle() const
{
   static char title[128];
   Float_t pmom = TMath::Sqrt(fPx*fPx+fPy*fPy+fPz*fPz);
   snprintf(title,128,"pmom=%f GeV",pmom);
   return title;
}

////////////////////////////////////////////////////////////////////////////////
///
///  Paint a primary track
///

void TPrimary::Paint(Option_t *option)
{
   Float_t rmin[3], rmax[3];
   static TPolyLine3D *pline = 0;
   if (!pline) {
      pline = new TPolyLine3D(2);
   }
   Float_t pmom = TMath::Sqrt(fPx*fPx+fPy*fPy+fPz*fPz);
   if (pmom == 0) return;
   TView *view = gPad->GetView();
   if (!view) return;
   view->GetRange(rmin,rmax);
   Float_t rbox = rmax[2];
   pline->SetPoint(0,fVx, fVy, fVz);
   Float_t xend = fVx+rbox*fPx/pmom;
   Float_t yend = fVy+rbox*fPy/pmom;
   Float_t zend = fVz+rbox*fPz/pmom;
   pline->SetPoint(1, xend, yend, zend);
   pline->SetLineColor(GetLineColor());
   pline->SetLineStyle(GetLineStyle());
   pline->SetLineWidth(GetLineWidth());
   pline->Paint(option);
}

////////////////////////////////////////////////////////////////////////////////
///
///  Print the internals of the primary vertex particle
///

void TPrimary::Print(Option_t *) const
{
   char def[8] = "XXXXXXX";
   const char *name;
   TAttParticle *ap = (TAttParticle*)GetParticle();
   if (ap) name = ap->GetName();
   else    name = def;
   Printf("TPrimary: %-13s  p: %8f %8f %8f Vertex: %8e %8e %8e %5d %5d %s",
   name,fPx,fPy,fPz,fVx,fVy,fVz,
   fFirstMother,fSecondMother,fType.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Return total X3D size of this primary
///

void TPrimary::Sizeof3D() const
{
   Float_t pmom = TMath::Sqrt(fPx*fPx+fPy*fPy+fPz*fPz);
   if (pmom == 0) return;
   Int_t npoints = 2;
   gSize3D.numPoints += npoints;
   gSize3D.numSegs   += (npoints-1);
   gSize3D.numPolys  += 0;

}

