/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
/// Data and categories: latex printing of lists and sets of RooArgSets
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
#include "RooExponential.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit;

void rf407_latextables()
{
   // S e t u p   c o m p o s i t e    p d f
   // --------------------------------------

   // Declare observable x
   RooRealVar x("x", "x", 0, 10);

   // Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their parameters
   RooRealVar mean("mean", "mean of gaussians", 5);
   RooRealVar sigma1("sigma1", "width of gaussians", 0.5);
   RooRealVar sigma2("sigma2", "width of gaussians", 1);
   RooGaussian sig1("sig1", "Signal component 1", x, mean, sigma1);
   RooGaussian sig2("sig2", "Signal component 2", x, mean, sigma2);

   // Sum the signal components into a composite signal pdf
   RooRealVar sig1frac("sig1frac", "fraction of component 1 in signal", 0.8, 0., 1.);
   RooAddPdf sig("sig", "Signal", RooArgList(sig1, sig2), sig1frac);

   // Build Chebychev polynomial pdf
   RooRealVar a0("a0", "a0", 0.5, 0., 1.);
   RooRealVar a1("a1", "a1", 0.2, 0., 1.);
   RooChebychev bkg1("bkg1", "Background 1", x, RooArgSet(a0, a1));

   // Build expontential pdf
   RooRealVar alpha("alpha", "alpha", -1);
   RooExponential bkg2("bkg2", "Background 2", x, alpha);

   // Sum the background components into a composite background pdf
   RooRealVar bkg1frac("sig1frac", "fraction of component 1 in background", 0.2, 0., 1.);
   RooAddPdf bkg("bkg", "Signal", RooArgList(bkg1, bkg2), sig1frac);

   // Sum the composite signal and background
   RooRealVar bkgfrac("bkgfrac", "fraction of background", 0.5, 0., 1.);
   RooAddPdf model("model", "g1+g2+a", RooArgList(bkg, sig), bkgfrac);

   // M a k e   l i s t   o f   p a r a m e t e r s   b e f o r e   a n d   a f t e r   f i t
   // ----------------------------------------------------------------------------------------

   // Make list of model parameters
   RooArgSet *params = model.getParameters(x);

   // Save snapshot of prefit parameters
   RooArgSet *initParams = (RooArgSet *)params->snapshot();

   // Do fit to data, to obtain error estimates on parameters
   RooDataSet *data = model.generate(x, 1000);
   model.fitTo(*data);

   // P r i n t   l a t ex   t a b l e   o f   p a r a m e t e r s   o f   p d f
   // --------------------------------------------------------------------------

   // Print parameter list in LaTeX for (one column with names, one column with values)
   params->printLatex();

   // Print parameter list in LaTeX for (names values|names values)
   params->printLatex(Columns(2));

   // Print two parameter lists side by side (name values initvalues)
   params->printLatex(Sibling(*initParams));

   // Print two parameter lists side by side (name values initvalues|name values initvalues)
   params->printLatex(Sibling(*initParams), Columns(2));

   // Write LaTex table to file
   params->printLatex(Sibling(*initParams), OutputFile("rf407_latextables.tex"));
}
