/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Likelihood and minimization: working with the profile likelihood estimator
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
#include "RooConstVar.h"
#include "RooAddPdf.h"
#include "RooMinimizer.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit;

void rf605_profilell()
{
   // C r e a t e   m o d e l   a n d   d a t a s e t
   // -----------------------------------------------

   // Observable
   RooRealVar x("x", "x", -20, 20);

   // Model (intentional strong correlations)
   RooRealVar mean("mean", "mean of g1 and g2", 0, -10, 10);
   RooRealVar sigma_g1("sigma_g1", "width of g1", 3);
   RooGaussian g1("g1", "g1", x, mean, sigma_g1);

   RooRealVar sigma_g2("sigma_g2", "width of g2", 4, 3.0, 6.0);
   RooGaussian g2("g2", "g2", x, mean, sigma_g2);

   RooRealVar frac("frac", "frac", 0.5, 0.0, 1.0);
   RooAddPdf model("model", "model", RooArgList(g1, g2), frac);

   // Generate 1000 events
   RooDataSet *data = model.generate(x, 1000);

   // C o n s t r u c t   p l a i n   l i k e l i h o o d
   // ---------------------------------------------------

   // Construct unbinned likelihood
   RooAbsReal *nll = model.createNLL(*data, NumCPU(2));

   // Minimize likelihood w.r.t all parameters before making plots
   RooMinimizer(*nll).migrad();

   // Plot likelihood scan frac
   RooPlot *frame1 = frac.frame(Bins(10), Range(0.01, 0.95), Title("LL and profileLL in frac"));
   nll->plotOn(frame1, ShiftToZero());

   // Plot likelihood scan in sigma_g2
   RooPlot *frame2 = sigma_g2.frame(Bins(10), Range(3.3, 5.0), Title("LL and profileLL in sigma_g2"));
   nll->plotOn(frame2, ShiftToZero());

   // C o n s t r u c t   p r o f i l e   l i k e l i h o o d   i n   f r a c
   // -----------------------------------------------------------------------

   // The profile likelihood estimator on nll for frac will minimize nll w.r.t
   // all floating parameters except frac for each evaluation

   RooAbsReal *pll_frac = nll->createProfile(frac);

   // Plot the profile likelihood in frac
   pll_frac->plotOn(frame1, LineColor(kRed));

   // Adjust frame maximum for visual clarity
   frame1->SetMinimum(0);
   frame1->SetMaximum(3);

   // C o n s t r u c t   p r o f i l e   l i k e l i h o o d   i n   s i g m a _ g 2
   // -------------------------------------------------------------------------------

   // The profile likelihood estimator on nll for sigma_g2 will minimize nll
   // w.r.t all floating parameters except sigma_g2 for each evaluation
   RooAbsReal *pll_sigmag2 = nll->createProfile(sigma_g2);

   // Plot the profile likelihood in sigma_g2
   pll_sigmag2->plotOn(frame2, LineColor(kRed));

   // Adjust frame maximum for visual clarity
   frame2->SetMinimum(0);
   frame2->SetMaximum(3);

   // Make canvas and draw RooPlots
   TCanvas *c = new TCanvas("rf605_profilell", "rf605_profilell", 800, 400);
   c->Divide(2);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   frame1->GetYaxis()->SetTitleOffset(1.4);
   frame1->Draw();
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   frame2->GetYaxis()->SetTitleOffset(1.4);
   frame2->Draw();

   delete pll_frac;
   delete pll_sigmag2;
   delete nll;
}
