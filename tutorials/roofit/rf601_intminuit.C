/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Likelihood and minimization: interactive minimization with MINUIT
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
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooMinimizer.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
using namespace RooFit;

void rf601_intminuit()
{
   // S e t u p   p d f   a n d   l i k e l i h o o d
   // -----------------------------------------------

   // Observable
   RooRealVar x("x", "x", -20, 20);

   // Model (intentional strong correlations)
   RooRealVar mean("mean", "mean of g1 and g2", 0);
   RooRealVar sigma_g1("sigma_g1", "width of g1", 3);
   RooGaussian g1("g1", "g1", x, mean, sigma_g1);

   RooRealVar sigma_g2("sigma_g2", "width of g2", 4, 3.0, 6.0);
   RooGaussian g2("g2", "g2", x, mean, sigma_g2);

   RooRealVar frac("frac", "frac", 0.5, 0.0, 1.0);
   RooAddPdf model("model", "model", RooArgList(g1, g2), frac);

   // Generate 1000 events
   RooDataSet *data = model.generate(x, 1000);

   // Construct unbinned likelihood of model w.r.t. data
   RooAbsReal *nll = model.createNLL(*data);

   // I n t e r a c t i v e   m i n i m i z a t i o n ,   e r r o r   a n a l y s i s
   // -------------------------------------------------------------------------------

   // Create MINUIT interface object
   RooMinimizer m(*nll);

   // Activate verbose logging of MINUIT parameter space stepping
   m.setVerbose(kTRUE);

   // Call MIGRAD to minimize the likelihood
   m.migrad();

   // Print values of all parameters, that reflect values (and error estimates)
   // that are back propagated from MINUIT
   model.getParameters(x)->Print("s");

   // Disable verbose logging
   m.setVerbose(kFALSE);

   // Run HESSE to calculate errors from d2L/dp2
   m.hesse();

   // Print value (and error) of sigma_g2 parameter, that reflects
   // value and error back propagated from MINUIT
   sigma_g2.Print();

   // Run MINOS on sigma_g2 parameter only
   m.minos(sigma_g2);

   // Print value (and error) of sigma_g2 parameter, that reflects
   // value and error back propagated from MINUIT
   sigma_g2.Print();

   // S a v i n g   r e s u l t s ,   c o n t o u r   p l o t s
   // ---------------------------------------------------------

   // Save a snapshot of the fit result. This object contains the initial
   // fit parameters, the final fit parameters, the complete correlation
   // matrix, the EDM, the minimized FCN , the last MINUIT status code and
   // the number of times the RooFit function object has indicated evaluation
   // problems (e.g. zero probabilities during likelihood evaluation)
   RooFitResult *r = m.save();

   // Make contour plot of mx vs sx at 1,2,3 sigma
   RooPlot *frame = m.contour(frac, sigma_g2, 1, 2, 3);
   frame->SetTitle("Minuit contour plot");

   // Print the fit result snapshot
   r->Print("v");

   // C h a n g e   p a r a m e t e r   v a l u e s ,   f l o a t i n g
   // -----------------------------------------------------------------

   // At any moment you can manually change the value of a (constant)
   // parameter
   mean = 0.3;

   // Rerun MIGRAD,HESSE
   m.migrad();
   m.hesse();
   frac.Print();

   // Now fix sigma_g2
   sigma_g2.setConstant(kTRUE);

   // Rerun MIGRAD,HESSE
   m.migrad();
   m.hesse();
   frac.Print();

   new TCanvas("rf601_intminuit", "rf601_intminuit", 600, 600);
   gPad->SetLeftMargin(0.15);
   frame->GetYaxis()->SetTitleOffset(1.4);
   frame->Draw();
}
