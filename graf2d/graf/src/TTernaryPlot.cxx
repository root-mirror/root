// @(#)root/graf:$Id$
// Author: Olivier Couet   27/10/20

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"
#include "TVirtualPad.h"
#include "TTernaryPlot.h"
#include "TMath.h"
#include "TString.h"
#include "TGaxis.h"
#include "TPolyMarker.h"

#include <iostream>

ClassImp(TTernaryPlot);

/** \class TTernaryPlot
\ingroup BasicGraphics

Draw a ternary plot.

Example:

~~~ {.cpp}
void ternary_plot()
{
   TCanvas *cnv = new TCanvas("cnv", "Ternary plot", 600, 600);

   TTernaryPlot *tp = new TTernaryPlot(3);

   tp->SetPoint(0.1, 0.8, "AC");
   tp->SetPoint(0.4, 0.1, "AB");
   tp->SetPoint(0.5, 0.1, "BC");

   tp->Draw();
}
~~~

*/


////////////////////////////////////////////////////////////////////////////////
/// Constructor with only the number of points set
/// the arrays x and y will be set later

TTernaryPlot::TTernaryPlot(Int_t n)
{
   fNpoints = 0;
   fMaxSize = n;
   fX       = new Double_t[fMaxSize];
   fY       = new Double_t[fMaxSize];
   for (Int_t i=0; i<fMaxSize; i++) fX[i] = fY[i] = 0.;
   yC = TMath::Sqrt(3.) / 2;
}


////////////////////////////////////////////////////////////////////////////////
/// TTernaryPlot default destructor.

TTernaryPlot::~TTernaryPlot()
{
   delete [] fX;
   delete [] fY;
}


////////////////////////////////////////////////////////////////////////////////
/// Add a point to the plot

void TTernaryPlot::SetPoint(Double_t u, Double_t v, Option_t *option)
{
   if(fNpoints==fMaxSize) {
      Error("SetPoint", "Invalid number of points");
      return;
   }

   Double_t a=0,b=0,c=0;
   TString opt = option;

   if (opt.Contains("AB")) {
      a = u;
      b = v;
   } else if (opt.Contains("AC")) {
      a = u;
      c = v;
      b = 1.-a-c;
   } else if (opt.Contains("BC")) {
      b = u;
      c = v;
      a = 1.-b-c;
   }

   fX[fNpoints] = b+0.5*a;
   fY[fNpoints] = a*yC;

   fNpoints++;
}


////////////////////////////////////////////////////////////////////////////////
/// Draw

void TTernaryPlot::Draw(Option_t * /*option*/)
{
   gPad->Range(-0.2, -0.2, 1.2, 1.2);

   TPolyMarker *ternaryPlot = new TPolyMarker(fNpoints, fX, fY);
   ternaryPlot->SetMarkerColor(kBlue);
   ternaryPlot->SetMarkerStyle(kFullTriangleDown);
   ternaryPlot->Draw();

   TGaxis *a1 = new TGaxis(0., 0., 0.5, yC, 0., 1.,10);
   a1->SetLineColor(kRed);
   a1->SetLabelFont(40);
   a1->SetLabelSize(0.025);
   a1->SetLabelOffset(0.06);
   a1->Draw();

   TGaxis *a2 = new TGaxis(0.5, yC, 1., 0., 0., 1.,10);
   a2->Draw();
   a2->SetLineColor(kBlue);
   a2->SetLabelFont(40);
   a2->SetLabelSize(0.025);
   a2->SetLabelOffset(0.06);


   TGaxis *a3 = new TGaxis(1., 0., 0., 0., 0., 1.,10);
   a3->Draw();
   a3->SetLabelFont(40);
   a3->SetLabelSize(0.025);
   a3->SetLabelOffset(-0.03);
}
