/// \file
/// \ingroup tutorial_v7
///
/// This macro generates really large RH2D histogram, fills it with predefined pattern and
/// draw it in a RCanvas, using Optmize() drawing mode
///
/// \macro_code
///
/// \date 2020-06-26
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
/// \author Sergey Linev <s.linev@gsi.de>

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RHistDrawable.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RFrameTitle.hxx"
#include "ROOT/RHistStatBox.hxx"
#include "ROOT/RFrame.hxx"

// macro must be here while cling is not capable to load
// library automatically for outlined function see ROOT-10336
R__LOAD_LIBRARY(libROOTHistDraw)

using namespace ROOT::Experimental;

void draw_rh2_large()
{
   const int nbins = 100;

   // Create the histogram.
   RAxisConfig xaxis("x", nbins, 0., nbins);
   RAxisConfig yaxis("y", nbins, 0., nbins);
   auto pHist = std::make_shared<RH2D>(xaxis, yaxis);

   for(int i=0;i<nbins;++i)
      for(int j=0;j<nbins;++j)
         pHist->Fill({1.*i,1.*j}, i+j);

   // Create a canvas to be displayed.
   auto canvas = RCanvas::Create("Canvas Title");

   auto frame = canvas->GetOrCreateFrame();

   // should we made special style for frame with palette?
   // frame->Margins().SetRight(0.2_normal);

   frame->SetGridX(false).SetGridY(false);
   frame->AttrX().SetZoomMinMax(nbins*0.2, nbins*0.8);
   frame->AttrY().SetZoomMinMax(nbins*0.2, nbins*0.8);

   canvas->Draw<RFrameTitle>(Form("Large RH2D histogram with %d x %d bins",nbins,nbins));

   auto draw = canvas->Draw(pHist);

   // draw->AttrLine().SetColor(RColor::kLime);
   // draw->Surf(2); // configure surf4 draw option
   // draw->Lego(2); // configure lego2 draw option
   // draw->Contour(); // configure cont draw option
   // draw->Scatter(); // configure color draw option (default)
   // draw->Arrow(); // configure arrow draw option
   // draw->Color(); // configure color draw option (default)
   // draw->Text(true); // configure text drawing (can be enabled with most 2d options)

   draw->Optimize(); // enable draw optimization, reduced data set will be send to clients

   auto stat = canvas->Draw<RHist2StatBox>(pHist, "hist2");
   stat->AttrFill().SetColor(RColor::kRed);

   canvas->SetSize(1000, 700);
   canvas->Show();
}
