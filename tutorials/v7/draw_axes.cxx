/// \file
/// \ingroup tutorial_v7
///
/// This ROOT7 example demonstrates how to create a RCanvas and
/// draw several RAxis objects with different options.
///
/// \macro_image (axes.png)
/// \macro_code
///
/// \date 2020-11-03
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!
/// \authors Sergey Linev <S.Linev@gsi.de>

#include <ROOT/RCanvas.hxx>
#include <ROOT/RAxisDrawable.hxx>

#include "TDatime.h"

// macro must be here while cling is not capable to load
// library automatically for outlined function see ROOT-10336
R__LOAD_LIBRARY(libROOTGpadv7)


void draw_axes()
{
   using namespace ROOT::Experimental;

   // Create a canvas to be displayed.
   auto canvas = RCanvas::Create("RAxis examples");

   auto x1 = 0.08_normal, w1 = 0.36_normal, x2 = 0.57_normal, w2 = 0.36_normal;

   auto draw0 = canvas->Draw<RAxisDrawable>(RPadPos(0.03_normal,0.1_normal), true, 0.8_normal);
   draw0->AttrAxis().SetTitle("vertical").SetTicksInvert().SetTicksSize(0.02_normal);

   auto draw1 = canvas->Draw<RAxisDrawable>(RPadPos(x1, 0.9_normal), false, w1);
   draw1->AttrAxis().SetTitle("horizontal").SetTitleCenter();

   auto draw2 = canvas->Draw<RAxisDrawable>(RPadPos(x1, 0.7_normal), false, w1);
   draw2->AttrAxis().SetTicksBoth().SetTicksSize(0.02_normal)
                    .SetTitle("both side ticks").SetTitleCenter();

   auto draw3 = canvas->Draw<RAxisDrawable>(RPadPos(x1, 0.5_normal), false, w1);
   draw3->AttrAxis().SetTitle("center labels").SetLabelsCenter();

   auto draw4 = canvas->Draw<RAxisDrawable>(RPadPos(x1, 0.3_normal), false, w1);
   draw4->SetMinMax(TDatime(2020,11,12,9,0,0).Convert(), TDatime(2020,11,12,12,0,0).Convert())
          .AttrAxis().SetTimeDisplay("%d/%m/%y %H:%M").SetTitle("time display").LabelsAttr().SetSize(0.01).SetColor(RColor::kRed);

   std::vector<std::string> labels = {"first", "second", "third", "forth", "fifth"};
   auto draw5 = canvas->Draw<RAxisLabelsDrawable>(RPadPos(x1, 0.1_normal), false, w1);
   draw5->SetLabels(labels).AttrAxis().SetTitle("labels, swap ticks side").SetTitleLeft().SetTicksInvert();

   auto draw6 = canvas->Draw<RAxisDrawable>(RPadPos(0.5_normal,0.9_normal), true, -0.8_normal);
   draw6->SetMinMax(0, 10).AttrAxis().SetTitle("vertical negative length").SetEndingArrow();

   auto draw7 = canvas->Draw<RAxisDrawable>(RPadPos(x2, 0.9_normal), false, w2);
   draw7->SetMinMax(1, 100).AttrAxis().SetLog(10).SetTitle("log10 scale").SetTitleCenter().TitleAttr().SetFont(12).SetColor(RColor::kGreen);

   auto draw8 = canvas->Draw<RAxisDrawable>(RPadPos(x2, 0.7_normal), false, w2);
   draw8->SetMinMax(0.125, 128).AttrAxis().SetLog(2).SetTitle("log2 scale").SetTitleCenter();

   auto draw9 = canvas->Draw<RAxisDrawable>(RPadPos(x2, 0.5_normal), false, w2);
   draw9->SetMinMax(1, 100).AttrAxis().SetLog(2.7182).SetTitle("ln scale").SetTitleCenter();

   auto draw10 = canvas->Draw<RAxisDrawable>(RPadPos(x2+w2, 0.3_normal), false, -w2);
   draw10->AttrAxis().SetTitle("horizontal negative length").SetTitleCenter().SetEndingArrow();

   auto draw11 = canvas->Draw<RAxisDrawable>(RPadPos(x2, 0.1_normal), false, w2);
   draw11->AttrAxis().SetReverse().SetTitle("horizontal reverse").SetTitleCenter().SetEndingArrow();

   auto draw12 = canvas->Draw<RAxisDrawable>(RPadPos(0.97_normal, 0.1_normal), true, 0.8_normal);
   draw12->AttrAxis().SetTitle("vertical axis with arrow").SetTitleCenter()
                     .SetTicksBoth().SetTicksSize(0.01_normal).SetTicksColor(RColor::kBlue)
                     .SetEndingArrow().SetEndingSize(0.01_normal);

   canvas->SetSize(1000, 800);

   // requires Chrome browser, runs in headless mode
   // canvas->SaveAs("axes.png");

   canvas->Show();
}
