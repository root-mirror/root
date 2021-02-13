/// \file
/// \ingroup tutorial_v7
///
/// This macro generates two RH1D, fills them and draw with different options in RCanvas.
/// The canvas is display in the web browser
///
/// \macro_code
///
/// \date 2015-03-22
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
/// \authors Axel Naumann <axel@cern.ch>, Sergey Linev <s.linev@gsi.de>

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

void draw_rh1()
{
   // Create the histogram.
   RAxisConfig xaxis(25, 0., 10.);
   auto pHist1 = std::make_shared<RH1D>(xaxis);
   auto pHist2 = std::make_shared<RH1D>(xaxis);

   for (int n=0;n<1000;n++) {
      pHist1->Fill(gRandom->Gaus(3,0.8));
      pHist2->Fill(gRandom->Gaus(7,1.2));
   }

   // Create a canvas to be displayed.
   auto canvas = RCanvas::Create("RH1 drawing options");

   // histograms colors
   auto col1 = RColor::kRed, col2 = RColor::kBlue;

   // Divide canvas on 2x3 sub-pads to show different draw options
   auto subpads = canvas->Divide(2,3);

   // default draw option
   subpads[0][0]->Draw<RFrameTitle>("Default RH1 drawing");
   subpads[0][0]->Draw(pHist1)->AttrLine().SetColor(col1).SetWidth(2);
   subpads[0][0]->Draw(pHist2)->AttrLine().SetColor(col2).SetWidth(4);

   // errors draw options
   subpads[1][0]->Draw<RFrameTitle>("Error() draw options");
   subpads[1][0]->Draw(pHist1)->Error(1).AttrLine().SetColor(col1);
   subpads[1][0]->Draw(pHist2)->Error(4).AttrFill().SetColor(col2).SetStyle(3003);

   // text and marker draw options
   subpads[0][1]->Draw<RFrameTitle>("Text() and Marker() draw options");
   subpads[0][1]->Draw(pHist1)->Text(true).AttrText().SetColor(col1);
   subpads[0][1]->Draw(pHist2)->Marker().AttrMarker().SetColor(col2).SetStyle(30).SetSize(1.5);

   // text and marker draw options
   subpads[1][1]->Draw<RFrameTitle>("Bar() draw options");
   subpads[1][1]->Draw(pHist1)->Bar(0,0.5).AttrFill().SetColor(col1);
   subpads[1][1]->Draw(pHist2)->Bar(0.5,0.5,true).AttrFill().SetColor(col2);

   // line draw option
   subpads[0][2]->Draw<RFrameTitle>("Line() draw option");
   subpads[0][2]->Draw(pHist1)->Line().AttrLine().SetColor(col1);
   subpads[0][2]->Draw(pHist2)->Line().AttrLine().SetColor(col2);

   // lego draw option
   subpads[1][2]->Draw<RFrameTitle>("Lego() draw option");
   subpads[1][2]->Draw(pHist1)->Lego().AttrFill().SetColor(col1);

   canvas->SetSize(1000, 700);
   canvas->Show();
}
