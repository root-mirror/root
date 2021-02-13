/// \file
/// \ingroup tutorial_v7
///
/// \macro_code
///
/// \date 2015-03-22
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!
/// \author Axel Naumann <axel@cern.ch>

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RCanvas.hxx"
#include "ROOT/RColor.hxx"
#include "ROOT/RHistDrawable.hxx"

// macro must be here while cling is not capable to load
// library automatically for outlined function see ROOT-10336
R__LOAD_LIBRARY(libROOTHistDraw)

void draw()
{
   using namespace ROOT::Experimental;

   // Create the histogram.
   RAxisConfig xaxis("x", 10, 0., 1.);
   RAxisConfig yaxis("y", {0., 1., 2., 3., 10.});
   auto pHist = std::make_shared<RH2D>(xaxis, yaxis);

   // Fill a few points.
   pHist->Fill({0.01, 1.02});
   pHist->Fill({0.54, 3.02});
   pHist->Fill({0.98, 1.02});
   pHist->Fill({1.90, 1.02});
   pHist->Fill({0.75, -0.02});

   // Create a canvas to be displayed.
   auto canvas = RCanvas::Create("Canvas Title");
   auto draw1 = canvas->Draw(pHist);
   // draw1->AttrLine().SetColor(RColor::kRed);

   auto other = std::make_shared<RH2D>(*pHist);
   auto draw2 = canvas->Draw(other);
   // draw2->AttrLine().SetColor(RColor::kBlue).SetWidth(12);

   canvas->Show();
}
