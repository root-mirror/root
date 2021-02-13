/// \file
/// \ingroup tutorial_v7
///
/// This macro generates a small V7 TH2D, fills it with random values and
/// draw it in a V7 canvas, using configured web browser
///
/// \macro_code
///
/// \date 2020-03-04
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
#include "ROOT/RPaletteDrawable.hxx"
#include "ROOT/RHistStatBox.hxx"
#include "ROOT/RFrame.hxx"
#include "TRandom.h"

// macro must be here while cling is not capable to load
// library automatically for outlined function see ROOT-10336
R__LOAD_LIBRARY(libROOTHistDraw)

using namespace ROOT::Experimental;

void draw_rh2_colz()
{
   // Create the histogram.
   RAxisConfig xaxis("x", 20, 0., 10.);
   RAxisConfig yaxis("y", 20, 0., 10.);
   auto pHist = std::make_shared<RH2D>(xaxis, yaxis);

   for (int n=0;n<10000;n++)
      pHist->Fill({gRandom->Gaus(5.,2.), gRandom->Gaus(5.,2.)});

   // Create a canvas to be displayed.
   auto canvas = RCanvas::Create("Canvas Title");

   auto frame = canvas->GetOrCreateFrame();

   // should we made special style for frame with palette?
   frame->Margins().SetRight(0.2_normal);

   frame->SetGridX(false).SetGridY(false);

   // draw ticks on both sides
   frame->SetTicksX(2).SetTicksY(2);

   // swap frame side where axes are drawn
   // frame->SetSwapX(true).SetSwapY(true);

   frame->AttrX().SetZoom(2.,8.);

   frame->AttrY().SetZoom(2.,8.);

   canvas->Draw<RFrameTitle>("2D histogram with color palette")->SetMargin(0.01_normal).SetHeight(0.09_normal);

   canvas->Draw<RPaletteDrawable>(RPalette::GetPalette(), true);

   auto draw = canvas->Draw(pHist);
   // draw->AttrLine().SetColor(RColor::kLime);
   // draw->Surf(2); // configure surf4 draw option
   // draw->Lego(2); // configure lego2 draw option
   // draw->Contour(); // configure cont draw option
   // draw->Scatter(); // configure color draw option (default)
   // draw->Arrow(); // configure arrow draw option
   draw->Color(); // configure color draw option (default)
   draw->Text(true); // configure text drawing (can be enabled with most 2d options)

   auto stat = canvas->Draw<RHist2StatBox>(pHist, "hist2");
   stat->AttrFill().SetColor(RColor::kRed);

   canvas->SetSize(1000, 700);
   canvas->Show();

   //canvas->Show("1000x700");

   // canvas->SaveAs("rh2_colz.png");
}
