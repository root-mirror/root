/// \file
/// \ingroup tutorial_v7
///
/// This macro generates a small RH3D, fills it with random values and
/// draw it in RCanvas, using configured web browser
///
/// \macro_code
///
/// \date 2020-06-18
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
/// \author Sergey Linev <s.linev@gsi.de>

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RHist.hxx"
#include "ROOT/RHistDrawable.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RPad.hxx"
#include "ROOT/RFrameTitle.hxx"
#include "TRandom.h"

// macro must be here while cling is not capable to load
// library automatically for outlined function see ROOT-10336
R__LOAD_LIBRARY(libROOTHistDraw)

using namespace ROOT::Experimental;

void draw_rh3()
{
   // Create the histogram.
   RAxisConfig xaxis("x", 10, -5., 5.);
   RAxisConfig yaxis("y", 10, -5., 5.);
   RAxisConfig zaxis("z", 10, -5., 5.);
   auto pHist = std::make_shared<RH3D>(xaxis, yaxis, zaxis);

   for (int n=0;n<10000;n++)
      pHist->Fill({gRandom->Gaus(0.,2.), gRandom->Gaus(0.,2.), gRandom->Gaus(0.,2.)});

   // Create a canvas to be displayed.
   auto canvas = RCanvas::Create("RH3D drawing options");

   // Divide canvas on 2x2 sub-pads to show different draw options
   auto subpads = canvas->Divide(2,2);

   // default draw option
   subpads[0][0]->Draw<RFrameTitle>("Box(0) default draw option");
   subpads[0][0]->Draw(pHist)->Box(0).AttrFill().SetColor(RColor::kBlue);

   // sphere draw options
   subpads[1][0]->Draw<RFrameTitle>("Sphere(1) draw option");
   subpads[1][0]->Draw(pHist)->Sphere(1);

   // text draw options
   subpads[0][1]->Draw<RFrameTitle>("Color() draw option");
   subpads[0][1]->Draw(pHist)->Color();

   // arrow draw options
   subpads[1][1]->Draw<RFrameTitle>("Scatter() draw option");
   subpads[1][1]->Draw(pHist)->Scatter().AttrFill().SetColor(RColor::kBlack);

   canvas->SetSize(1000, 700);
   canvas->Show();
}
