/// \file
/// \ingroup tutorial_v7
///
/// This macro shows how ROOT objects like TH1, TH2, TGraph can be drawn in RCanvas.
///
/// \macro_code
///
/// \date 2017-06-02
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!
/// \authors Axel Naumann <axel@cern.ch>, Sergey Linev <s.linev@gsi.de>

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/TObjectDrawable.hxx>
#include <ROOT/RPad.hxx>
#include <ROOT/RCanvas.hxx>
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TMath.h"
#include "TStyle.h"

#include <iostream>

void draw_v6()
{
   using namespace ROOT::Experimental;

   static constexpr int npoints = 10;
   double x[npoints] = { 0., 1., 2., 3., 4., 5., 6., 7., 8., 9. };
   double y[npoints] = { .1, .2, .3, .4, .3, .2, .1, .2, .3, .4 };
   auto gr = std::make_shared<TGraph>(npoints, x, y);

   static constexpr int nth1points = 100;
   auto th1 = std::make_shared<TH1I>("gaus", "Example of TH1", nth1points, -5, 5);
   th1->SetDirectory(nullptr);
   for (int n=0;n<nth1points;++n) {
      double x = 10.*n/nth1points-5.;
      th1->SetBinContent(n+1, (int) (1000*TMath::Gaus(x)));
   }

   static constexpr int nth2points = 40;
   auto th2 = std::make_shared<TH2I>("gaus2", "Example of TH1", nth2points, -5, 5, nth2points, -5, 5);
   th2->SetDirectory(nullptr);
   for (int n=0;n<nth2points;++n) {
      for (int k=0;k<nth2points;++k) {
         double x = 10.*n/nth2points-5.;
         double y = 10.*k/nth2points-5.;
         th2->SetBinContent(th2->GetBin(n+1, k+1), (int) (1000*TMath::Gaus(x)*TMath::Gaus(y)));
      }
   }

   gStyle->SetPalette(kRainBow);

   auto canvas = RCanvas::Create("RCanvas showing a v6 objects");

   // place copy of gStyle object, will be applied on JSROOT side
   // set on the canvas before any other object is drawn
   canvas->Draw<TObjectDrawable>(TObjectDrawable::kStyle);

   // copy all existing ROOT colors, required when colors was modified
   // or when colors should be possible from client side
   canvas->Draw<TObjectDrawable>(TObjectDrawable::kColors);

   // copy custom palette to canvas, will be used for col drawings
   // style object does not include color settings
   canvas->Draw<TObjectDrawable>(TObjectDrawable::kPalette);

   // Divide canvas on 2x2 sub-pads to show different draw options
   auto subpads = canvas->Divide(2,2);

   subpads[0][0]->Draw<TObjectDrawable>(gr, "AL");

   subpads[0][1]->Draw<TObjectDrawable>(th1, "");

   subpads[1][0]->Draw<TObjectDrawable>(th2, "colz");

   // show same object again, but with other draw options
   subpads[1][1]->Draw<TObjectDrawable>(th2, "lego2");

   canvas->Show(); // new window in default browser should popup and async update will be triggered

   // synchronous, wait until drawing is really finished
   canvas->Update(false, [](bool res) { std::cout << "First sync update done = " << (res ? "true" : "false") << std::endl; });

   // invalidate canvas and force repainting with next Update()
   canvas->Modified();

   // call Update again, should return immediately if canvas was not modified
   canvas->Update(true, [](bool res) { std::cout << "Second async update done = " << (res ? "true" : "false") << std::endl; });

   std::cout << "This message appear normally before second async update" << std::endl;

   // create SVG file
   // canvas->SaveAs("draw_v6.svg");

   // create PNG file
   // canvas->SaveAs("draw_v6.png");
}
