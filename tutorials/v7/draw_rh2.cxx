/// \file
/// \ingroup tutorial_v7
///
/// This macro generates RH2D and draw it with different options in RCanvas
///
/// \macro_code
///
/// \date 2020-06-25
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
/// \author Sergey Linev <s.linev@gsi.de>

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RHist.hxx"
#include "ROOT/RHistDrawable.hxx"
#include "ROOT/RFrameTitle.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RPad.hxx"
#include "TRandom.h"

// macro must be here while cling is not capable to load
// library automatically for outlined function see ROOT-10336
R__LOAD_LIBRARY(libROOTHistDraw)

using namespace ROOT::Experimental;

void draw_rh2()
{
   // Create the histogram.
   RAxisConfig xaxis("x", 20, 0., 10.);
   RAxisConfig yaxis("y", 20, 0., 10.);
   auto pHist = std::make_shared<RH2D>(xaxis, yaxis);

   for (int n=0;n<10000;n++)
      pHist->Fill({gRandom->Gaus(5.,2.), gRandom->Gaus(5.,2.)});

   // Create a canvas to be displayed.
   auto canvas = RCanvas::Create("RH2 drawing options");

   // Divide canvas on 2x3 sub-pads to show different draw options
   auto subpads = canvas->Divide(2,3);

   // default draw option
   subpads[0][0]->Draw<RFrameTitle>("Color() draw option (default)");
   subpads[0][0]->Draw(pHist);

   // contour draw options
   subpads[1][0]->Draw<RFrameTitle>("Contour() draw option");
   subpads[1][0]->Draw(pHist)->Contour();

   // text draw options
   subpads[0][1]->Draw<RFrameTitle>("Text() draw option");
   subpads[0][1]->Draw(pHist)->Text(true).AttrText().SetColor(RColor::kBlue);

   // arrow draw options
   subpads[1][1]->Draw<RFrameTitle>("Arrow() draw option");
   subpads[1][1]->Draw(pHist)->Arrow().AttrLine().SetColor(RColor::kRed);

   // lego draw options
   subpads[0][2]->Draw<RFrameTitle>("Lego() draw option");
   subpads[0][2]->Draw(pHist)->Lego(2);

   // surf draw option
   subpads[1][2]->Draw<RFrameTitle>("Surf() draw option");
   subpads[1][2]->Draw(pHist)->Surf(2);

   canvas->SetSize(1000, 700);
   canvas->Show();
}
