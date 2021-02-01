/// \file
/// \ingroup tutorial_graphics
/// \notebook
/// Draw a ternary plot
///
/// \macro_image
/// \macro_code
///
/// \author Olivier Couet

void ternaryplot()
{
   TCanvas *cnv = new TCanvas("cnv", "Ternary plot", 600, 600);

   TTernaryPlot *tp = new TTernaryPlot(3);

   tp->SetPoint(0.1, 0.8, "AC");
   tp->SetPoint(0.4, 0.1, "AB");
   tp->SetPoint(0.5, 0.1, "BC");

   tp->Draw();
}
