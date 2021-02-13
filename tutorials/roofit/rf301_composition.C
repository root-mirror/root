/// \file
/// \ingroup tutorial_roofit
/// \notebook
/// Multidimensional models: multi-dimensional pdfs through composition
/// e.g. substituting a pdf parameter with a function that depends on other observables
///
///  `pdf = gauss(x,f(y),s)` with `f(y) = a0 + a1*y`
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
#include "RooPolyVar.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
using namespace RooFit;

void rf301_composition()
{
   // S e t u p   c o m p o s e d   m o d e l   g a u s s ( x , m ( y ) , s )
   // -----------------------------------------------------------------------

   // Create observables
   RooRealVar x("x", "x", -5, 5);
   RooRealVar y("y", "y", -5, 5);

   // Create function f(y) = a0 + a1*y
   RooRealVar a0("a0", "a0", -0.5, -5, 5);
   RooRealVar a1("a1", "a1", -0.5, -1, 1);
   RooPolyVar fy("fy", "fy", y, RooArgSet(a0, a1));

   // Create gauss(x,f(y),s)
   RooRealVar sigma("sigma", "width of gaussian", 0.5);
   RooGaussian model("model", "Gaussian with shifting mean", x, fy, sigma);

   // S a m p l e   d a t a ,   p l o t   d a t a   a n d   p d f   o n   x   a n d   y
   // ---------------------------------------------------------------------------------

   // Generate 10000 events in x and y from model
   RooDataSet *data = model.generate(RooArgSet(x, y), 10000);

   // Plot x distribution of data and projection of model on x = Int(dy) model(x,y)
   RooPlot *xframe = x.frame();
   data->plotOn(xframe);
   model.plotOn(xframe);

   // Plot x distribution of data and projection of model on y = Int(dx) model(x,y)
   RooPlot *yframe = y.frame();
   data->plotOn(yframe);
   model.plotOn(yframe);

   // Make two-dimensional plot in x vs y
   TH1 *hh_model = model.createHistogram("hh_model", x, Binning(50), YVar(y, Binning(50)));
   hh_model->SetLineColor(kBlue);

   // Make canvas and draw RooPlots
   TCanvas *c = new TCanvas("rf301_composition", "rf301_composition", 1200, 400);
   c->Divide(3);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   xframe->GetYaxis()->SetTitleOffset(1.4);
   xframe->Draw();
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   yframe->GetYaxis()->SetTitleOffset(1.4);
   yframe->Draw();
   c->cd(3);
   gPad->SetLeftMargin(0.20);
   hh_model->GetZaxis()->SetTitleOffset(2.5);
   hh_model->Draw("surf");
}
