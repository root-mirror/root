/// \file
/// \ingroup tutorial_v7
///
/// \macro_code
///
/// \date 2018-03-18
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!
/// \authors Olivier couet, Iliana Betsou

R__LOAD_LIBRARY(libGpad);

#include "ROOT/TCanvas.hxx"
#include "ROOT/TColor.hxx"
#include "ROOT/TLine.hxx"
#include <ROOT/TPadPos.hxx>

void line()
{
   using namespace ROOT;

   // Create a canvas to be displayed.
   auto canvas = Experimental::TCanvas::Create("Canvas Title");

   for (double i = 0; i < 360; i+=1) {
      double ang = i * TMath::Pi() / 180;
      double x = 0.3*TMath::Cos(ang) + 0.5;
      double y = 0.3*TMath::Sin(ang) + 0.5;

      auto line = std::make_shared<Experimental::TLine>({0.5_normal, 0.5_normal} , {x_normal, y_normal});

      auto col = Experimental::TColor(0.0025*i, 0, 0);
      line->GetOptions().SetLineColor(col);
      line->GetOptions().SetLineWidth(1);

      canvas->Draw(line);
    }

  auto line  = std::make_shared<Experimental::TLine>({0.0_normal, 0.0_normal} , {1.0_normal, 1.0_normal}); canvas->Draw(line);
  auto line1 = std::make_shared<Experimental::TLine>({0.1_normal, 0.1_normal} , {0.9_normal, 0.1_normal}); canvas->Draw(line1);
  auto line2 = std::make_shared<Experimental::TLine>({0.9_normal, 0.1_normal} , {0.9_normal, 0.9_normal}); canvas->Draw(line2);
  auto line3 = std::make_shared<Experimental::TLine>({0.9_normal, 0.9_normal} , {0.1_normal, 0.9_normal}); canvas->Draw(line3);
  auto line4 = std::make_shared<Experimental::TLine>({0.1_normal, 0.1_normal} , {0.1_normal, 0.9_normal}); canvas->Draw(line4);
  auto line0 = std::make_shared<Experimental::TLine>({0.0_normal, 1.0_normal} , {1.0_normal, 0.0_normal}); canvas->Draw(line0);

   canvas->Show();
}
