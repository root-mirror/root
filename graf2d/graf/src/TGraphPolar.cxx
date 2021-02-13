// @(#)root/graf:$Id$
// Author: Sebastian Boser, Mathieu Demaret 02/02/06

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGraphPolar
\ingroup BasicGraphics

To draw a polar graph.

TGraphPolar creates a polar graph (including error bars). A TGraphPolar is
a TGraphErrors represented in polar coordinates.
It uses the class TGraphPolargram to draw the polar axis.

Example:

Begin_Macro(source)
{
   TCanvas * CPol = new TCanvas("CPol","TGraphPolar Example",500,500);

   Double_t theta[8];
   Double_t radius[8];
   Double_t etheta[8];
   Double_t eradius[8];

   for (int i=0; i<8; i++) {
      theta[i]   = (i+1)*(TMath::Pi()/4.);
      radius[i]  = (i+1)*0.05;
      etheta[i]  = TMath::Pi()/8.;
      eradius[i] = 0.05;
   }

   TGraphPolar * grP1 = new TGraphPolar(8, theta, radius, etheta, eradius);
   grP1->SetTitle("TGraphPolar Example");

   grP1->SetMarkerStyle(20);
   grP1->SetMarkerSize(2.);
   grP1->SetMarkerColor(4);
   grP1->SetLineColor(2);
   grP1->SetLineWidth(3);
   grP1->Draw("PE");

   // Update, otherwise GetPolargram returns 0
   CPol->Update();
   grP1->GetPolargram()->SetToRadian();

   return CPol;
}
End_Macro
*/

#include "TGraphPolar.h"
#include "TGraphPolargram.h"

ClassImp(TGraphPolar);

////////////////////////////////////////////////////////////////////////////////
/// TGraphPolar default constructor.

TGraphPolar::TGraphPolar() : TGraphErrors(),
             fOptionAxis(kFALSE),fPolargram(0),fXpol(0),fYpol(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphPolar constructor.
///
/// \param[in] n         number of points.
/// \param[in] theta     angular values.
/// \param[in] r         radial values.
/// \param[in] etheta    errors on angular values.
/// \param[in] er        errors on radial values.

TGraphPolar::TGraphPolar(Int_t n, const Double_t* theta, const Double_t* r,
                                  const Double_t *etheta, const Double_t* er)
  : TGraphErrors(n,theta,r,etheta,er),
             fOptionAxis(kFALSE),fPolargram(0),fXpol(0),fYpol(0)
{
   SetEditable(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// TGraphPolar destructor.

TGraphPolar::~TGraphPolar()
{
   delete []fXpol;
   delete []fYpol;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw TGraphPolar.

void TGraphPolar::Draw(Option_t* options)
{
   // Process options
   TString opt = options;
   opt.ToUpper();

   // Ignore same
   opt.ReplaceAll("SAME","");

   // ReDraw polargram if required by options
   if (opt.Contains("A")) fOptionAxis = kTRUE;
   opt.ReplaceAll("A","");

   AppendPad(opt);
}

////////////////////////////////////////////////////////////////////////////////
/// Return points in polar coordinates.

Double_t *TGraphPolar::GetXpol()
{
   if (!fXpol) fXpol = new Double_t[fNpoints];
   return fXpol;
}

////////////////////////////////////////////////////////////////////////////////
/// Return points in polar coordinates.

Double_t *TGraphPolar::GetYpol()
{
   if (!fYpol) fYpol = new Double_t[fNpoints];
   return fYpol;
}

////////////////////////////////////////////////////////////////////////////////
/// Set maximum Polar.

void TGraphPolar::SetMaxPolar(Double_t maximum)
{
   if (fPolargram) fPolargram->ChangeRangePolar(fPolargram->GetTMin(),maximum);
}

////////////////////////////////////////////////////////////////////////////////
/// Set maximum radial at the intersection of the positive X axis part and
/// the circle.

void TGraphPolar::SetMaxRadial(Double_t maximum)
{
   if (fPolargram) fPolargram->SetRangeRadial(fPolargram->GetRMin(),maximum);
}

////////////////////////////////////////////////////////////////////////////////
/// Set minimum Polar.

void TGraphPolar::SetMinPolar(Double_t minimum)
{
   if (fPolargram) fPolargram->ChangeRangePolar(minimum, fPolargram->GetTMax());
}

////////////////////////////////////////////////////////////////////////////////
/// Set minimum radial in the center of the circle.

void TGraphPolar::SetMinRadial(Double_t minimum)
{
   if (fPolargram) fPolargram->SetRangeRadial(minimum, fPolargram->GetRMax());
}
