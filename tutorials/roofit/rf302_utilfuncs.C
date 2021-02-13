/// \file
/// \ingroup tutorial_roofit
/// \notebook
/// Multidimensional models: utility functions classes available for use in tailoring of
/// composite (multidimensional) pdfs
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
#include "RooFormulaVar.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooPolyVar.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"

using namespace RooFit;

void rf302_utilfuncs()
{
   // C r e a t e   o b s e r v a b l e s ,   p a r a m e t e r s
   // -----------------------------------------------------------

   // Create observables
   RooRealVar x("x", "x", -5, 5);
   RooRealVar y("y", "y", -5, 5);

   // Create parameters
   RooRealVar a0("a0", "a0", -1.5, -5, 5);
   RooRealVar a1("a1", "a1", -0.5, -1, 1);
   RooRealVar sigma("sigma", "width of gaussian", 0.5);

   // U s i n g   R o o F o r m u l a V a r   t o   t a i l o r   p d f
   // -----------------------------------------------------------------------

   // Create interpreted function f(y) = a0 - a1*sqrt(10*abs(y))
   RooFormulaVar fy_1("fy_1", "a0-a1*sqrt(10*abs(y))", RooArgSet(y, a0, a1));

   // Create gauss(x,f(y),s)
   RooGaussian model_1("model_1", "Gaussian with shifting mean", x, fy_1, sigma);

   // U s i n g   R o o P o l y V a r   t o   t a i l o r   p d f
   // -----------------------------------------------------------------------

   // Create polynomial function f(y) = a0 + a1*y
   RooPolyVar fy_2("fy_2", "fy_2", y, RooArgSet(a0, a1));

   // Create gauss(x,f(y),s)
   RooGaussian model_2("model_2", "Gaussian with shifting mean", x, fy_2, sigma);

   // U s i n g   R o o A d d i t i o n   t o   t a i l o r   p d f
   // -----------------------------------------------------------------------

   // Create sum function f(y) = a0 + y
   RooAddition fy_3("fy_3", "a0+y", RooArgSet(a0, y));

   // Create gauss(x,f(y),s)
   RooGaussian model_3("model_3", "Gaussian with shifting mean", x, fy_3, sigma);

   // U s i n g   R o o P r o d u c t   t o   t a i l o r   p d f
   // -----------------------------------------------------------------------

   // Create product function f(y) = a1*y
   RooProduct fy_4("fy_4", "a1*y", RooArgSet(a1, y));

   // Create gauss(x,f(y),s)
   RooGaussian model_4("model_4", "Gaussian with shifting mean", x, fy_4, sigma);

   // P l o t   a l l   p d f s
   // ----------------------------

   // Make two-dimensional plots in x vs y
   TH1 *hh_model_1 = model_1.createHistogram("hh_model_1", x, Binning(50), YVar(y, Binning(50)));
   TH1 *hh_model_2 = model_2.createHistogram("hh_model_2", x, Binning(50), YVar(y, Binning(50)));
   TH1 *hh_model_3 = model_3.createHistogram("hh_model_3", x, Binning(50), YVar(y, Binning(50)));
   TH1 *hh_model_4 = model_4.createHistogram("hh_model_4", x, Binning(50), YVar(y, Binning(50)));
   hh_model_1->SetLineColor(kBlue);
   hh_model_2->SetLineColor(kBlue);
   hh_model_3->SetLineColor(kBlue);
   hh_model_4->SetLineColor(kBlue);

   // Make canvas and draw RooPlots
   TCanvas *c = new TCanvas("rf302_utilfuncs", "rf302_utilfuncs", 800, 800);
   c->Divide(2, 2);
   c->cd(1);
   gPad->SetLeftMargin(0.20);
   hh_model_1->GetZaxis()->SetTitleOffset(2.5);
   hh_model_1->Draw("surf");
   c->cd(2);
   gPad->SetLeftMargin(0.20);
   hh_model_2->GetZaxis()->SetTitleOffset(2.5);
   hh_model_2->Draw("surf");
   c->cd(3);
   gPad->SetLeftMargin(0.20);
   hh_model_3->GetZaxis()->SetTitleOffset(2.5);
   hh_model_3->Draw("surf");
   c->cd(4);
   gPad->SetLeftMargin(0.20);
   hh_model_4->GetZaxis()->SetTitleOffset(2.5);
   hh_model_4->Draw("surf");
}
