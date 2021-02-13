/// \file
/// \ingroup tutorial_v7
///
/// This macro generates a small V7 TH1D, fills it and draw it in a V7 canvas.
/// The canvas is display in the web browser
///
/// \macro_code
///
/// \date 2015-03-22
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
/// \author Axel Naumann <axel@cern.ch>

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RCanvas.hxx"
#include "ROOT/RPave.hxx"
#include "ROOT/RPaveText.hxx"

// macro must be here while cling is not capable to load
// library automatically for outlined function see ROOT-10336
R__LOAD_LIBRARY(ROOTGpadv7)

using namespace ROOT::Experimental;

void draw_pave()
{
   // Create a canvas to be displayed.
   auto canvas = RCanvas::Create("Canvas Title");

   // RFrame will be automatically created as well
   auto pave = canvas->Draw<RPave>();
   pave->AttrFill().SetColor(RColor::kBlue);
   pave->AttrBorder().SetColor(RColor::kGreen).SetWidth(3);
   pave->SetCornerY(-0.03_normal);
   pave->SetHeight(0.2_normal);

   auto text = canvas->Draw<RPaveText>();
   text->AddLine("This is RTextPave");
   text->AddLine("It can have several lines");
   text->AddLine("It should be positioned below RPave");
   text->AttrFill().SetColor(RColor::kYellow);
   text->SetCornerY(0.25_normal);
   text->SetHeight(0.3_normal);

   canvas->SetSize(1000, 700);
   canvas->Show();
}
