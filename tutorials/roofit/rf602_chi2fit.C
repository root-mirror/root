/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
/// Likelihood and minimization: setting up a chi^2 fit to a binned dataset
///
/// \macro_output
/// \macro_code
///
/// \date July 2008
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooChi2Var.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit;

void rf602_chi2fit()
{

   // S e t u p   m o d e l
   // ---------------------

   // Declare observable x
   RooRealVar x("x", "x", 0, 10);

   // Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their parameters
   RooRealVar mean("mean", "mean of gaussians", 5);
   RooRealVar sigma1("sigma1", "width of gaussians", 0.5);
   RooRealVar sigma2("sigma2", "width of gaussians", 1);

   RooGaussian sig1("sig1", "Signal component 1", x, mean, sigma1);
   RooGaussian sig2("sig2", "Signal component 2", x, mean, sigma2);

   // Build Chebychev polynomial pdf
   RooRealVar a0("a0", "a0", 0.5, 0., 1.);
   RooRealVar a1("a1", "a1", 0.2, 0., 1.);
   RooChebychev bkg("bkg", "Background", x, RooArgSet(a0, a1));

   // Sum the signal components into a composite signal pdf
   RooRealVar sig1frac("sig1frac", "fraction of component 1 in signal", 0.8, 0., 1.);
   RooAddPdf sig("sig", "Signal", RooArgList(sig1, sig2), sig1frac);

   // Sum the composite signal and background
   RooRealVar bkgfrac("bkgfrac", "fraction of background", 0.5, 0., 1.);
   RooAddPdf model("model", "g1+g2+a", RooArgList(bkg, sig), bkgfrac);

   // C r e a t e   b i n n e d   d a t a s e t
   // -----------------------------------------

   RooDataSet *d = model.generate(x, 10000);
   RooDataHist *dh = d->binnedClone();

   // Construct a chi^2 of the data and the model.
   // When a pdf is used in a chi^2 fit, the probability density scaled
   // by the number of events in the dataset to obtain the fit function
   // If model is an extended pdf, the expected number events is used
   // instead of the observed number of events.
   model.chi2FitTo(*dh);

   // NB: It is also possible to fit a RooAbsReal function to a RooDataHist
   // using chi2FitTo().

   // Note that entries with zero bins are _not_ allowed
   // for a proper chi^2 calculation and will give error
   // messages
   RooDataSet *dsmall = (RooDataSet *)d->reduce(EventRange(1, 100));
   RooDataHist *dhsmall = dsmall->binnedClone();
   RooChi2Var chi2_lowstat("chi2_lowstat", "chi2", model, *dhsmall);
   cout << chi2_lowstat.getVal() << endl;
}
