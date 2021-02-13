/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
/// Addition and convolution: tools and utilities for manipulation of composite objects
///
/// \macro_output
/// \macro_code
///
/// \date July 2008
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooCustomizer.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
using namespace RooFit;

void rf207_comptools()
{

   // S e t u p   c o m p o s i t e    p d f,   d a t a s e t
   // --------------------------------------------------------

   // Declare observable x
   RooRealVar x("x", "x", 0, 10);

   // Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their parameters
   RooRealVar mean("mean", "mean of gaussians", 5);
   RooRealVar sigma("sigma", "width of gaussians", 0.5);
   RooGaussian sig("sig", "Signal component 1", x, mean, sigma);

   // Build Chebychev polynomial pdf
   RooRealVar a0("a0", "a0", 0.5, 0., 1.);
   RooRealVar a1("a1", "a1", 0.2, 0., 1.);
   RooChebychev bkg1("bkg1", "Background 1", x, RooArgSet(a0, a1));

   // Build exponential pdf
   RooRealVar alpha("alpha", "alpha", -1);
   RooExponential bkg2("bkg2", "Background 2", x, alpha);

   // Sum the background components into a composite background pdf
   RooRealVar bkg1frac("bkg1frac", "fraction of component 1 in background", 0.2, 0., 1.);
   RooAddPdf bkg("bkg", "Signal", RooArgList(bkg1, bkg2), bkg1frac);

   // Sum the composite signal and background
   RooRealVar bkgfrac("bkgfrac", "fraction of background", 0.5, 0., 1.);
   RooAddPdf model("model", "g1+g2+a", RooArgList(bkg, sig), bkgfrac);

   // Create dummy dataset that has more observables than the above pdf
   RooRealVar y("y", "y", -10, 10);
   RooDataSet data("data", "data", RooArgSet(x, y));

   // ---------------------------------------------------
   // B a s i c   i n f o r m a t i o n   r e q u e s t s
   // ===================================================

   // G e t   l i s t   o f   o b s e r v a b l e s
   // ---------------------------------------------

   // Get list of observables of pdf in context of a dataset
   //
   // Observables are define each context as the variables
   // shared between a model and a dataset. In this case
   // that is the variable 'x'

   RooArgSet *model_obs = model.getObservables(data);
   model_obs->Print("v");

   // G e t   l i s t   o f   p a r a m e t e r s
   // -------------------------------------------

   // Get list of parameters, given list of observables
   RooArgSet *model_params = model.getParameters(x);
   model_params->Print("v");

   // Get list of parameters, given a dataset
   // (Gives identical results to operation above)
   RooArgSet *model_params2 = model.getParameters(data);
   model_params2->Print();

   // G e t   l i s t   o f   c o m p o n e n t s
   // -------------------------------------------

   // Get list of component objects, including top-level node
   RooArgSet *model_comps = model.getComponents();
   model_comps->Print("v");

   // -------------------------------------------------------------------------------
   // M o d i f i c a t i o n s   t o   s t r u c t u r e   o f   c o m p o s i t e s
   // ===============================================================================

   // Create a second Gaussian
   RooRealVar sigma2("sigma2", "width of gaussians", 1);
   RooGaussian sig2("sig2", "Signal component 1", x, mean, sigma2);

   // Create a sum of the original Gaussian plus the new second Gaussian
   RooRealVar sig1frac("sig1frac", "fraction of component 1 in signal", 0.8, 0., 1.);
   RooAddPdf sigsum("sigsum", "sig+sig2", RooArgList(sig, sig2), sig1frac);

   // Construct a customizer utility to customize model
   RooCustomizer cust(model, "cust");

   // Instruct the customizer to replace node 'sig' with node 'sigsum'
   cust.replaceArg(sig, sigsum);

   // Build a clone of the input pdf according to the above customization
   // instructions. Each node that requires modified is clone so that the
   // original pdf remained untouched. The name of each cloned node is that
   // of the original node suffixed by the name of the customizer object
   //
   // The returned head node own all nodes that were cloned as part of
   // the build process so when cust_clone is deleted so will all other
   // nodes that were created in the process.
   RooAbsPdf *cust_clone = (RooAbsPdf *)cust.build(kTRUE);

   // Print structure of clone of model with sig->sigsum replacement.
   cust_clone->Print("t");

   delete cust_clone;
}
