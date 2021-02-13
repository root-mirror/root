/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Basic functionality: various plotting styles of data, functions in a RooPlot
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \date July 2008
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit;

void rf107_plotstyles()
{

   // S e t u p   m o d e l
   // ---------------------

   // Create observables
   RooRealVar x("x", "x", -10, 10);

   // Create Gaussian
   RooRealVar sigma("sigma", "sigma", 3, 0.1, 10);
   RooRealVar mean("mean", "mean", -3, -10, 10);
   RooGaussian gauss("gauss", "gauss", x, mean, sigma);

   // Generate a sample of 100 events with sigma=3
   RooDataSet *data = gauss.generate(x, 100);

   // Fit pdf to data
   gauss.fitTo(*data);

   // M a k e   p l o t   f r a m e s
   // -------------------------------

   // Make four plot frames to demonstrate various plotting features
   RooPlot *frame1 = x.frame(Name("xframe"), Title("Red Curve / SumW2 Histo errors"), Bins(20));
   RooPlot *frame2 = x.frame(Name("xframe"), Title("Dashed Curve / No XError bars"), Bins(20));
   RooPlot *frame3 = x.frame(Name("xframe"), Title("Filled Curve / Blue Histo"), Bins(20));
   RooPlot *frame4 = x.frame(Name("xframe"), Title("Partial Range / Filled Bar chart"), Bins(20));

   // D a t a   p l o t t i n g   s t y l e s
   // ---------------------------------------

   // Use sqrt(sum(weights^2)) error instead of Poisson errors
   data->plotOn(frame1, DataError(RooAbsData::SumW2));

   // Remove horizontal error bars
   data->plotOn(frame2, XErrorSize(0));

   // Blue markers and error bors
   data->plotOn(frame3, MarkerColor(kBlue), LineColor(kBlue));

   // Filled bar chart
   data->plotOn(frame4, DrawOption("B"), DataError(RooAbsData::None), XErrorSize(0), FillColor(kGray));

   // F u n c t i o n   p l o t t i n g   s t y l e s
   // -----------------------------------------------

   // Change line color to red
   gauss.plotOn(frame1, LineColor(kRed));

   // Change line style to dashed
   gauss.plotOn(frame2, LineStyle(kDashed));

   // Filled shapes in green color
   gauss.plotOn(frame3, DrawOption("F"), FillColor(kOrange), MoveToBack());

   //
   gauss.plotOn(frame4, Range(-8, 3), LineColor(kMagenta));

   TCanvas *c = new TCanvas("rf107_plotstyles", "rf107_plotstyles", 800, 800);
   c->Divide(2, 2);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   frame1->GetYaxis()->SetTitleOffset(1.6);
   frame1->Draw();
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   frame2->GetYaxis()->SetTitleOffset(1.6);
   frame2->Draw();
   c->cd(3);
   gPad->SetLeftMargin(0.15);
   frame3->GetYaxis()->SetTitleOffset(1.6);
   frame3->Draw();
   c->cd(4);
   gPad->SetLeftMargin(0.15);
   frame4->GetYaxis()->SetTitleOffset(1.6);
   frame4->Draw();
}
