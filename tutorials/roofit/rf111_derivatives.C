/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Basic functionality: numerical 1st,2nd and 3rd order derivatives w.r.t. observables and parameters
///
/// ```
///  pdf = gauss(x,m,s)
/// ```
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

void rf111_derivatives()
{
   // S e t u p   m o d e l
   // ---------------------

   // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
   RooRealVar x("x", "x", -10, 10);
   RooRealVar mean("mean", "mean of gaussian", 1, -10, 10);
   RooRealVar sigma("sigma", "width of gaussian", 1, 0.1, 10);

   // Build gaussian pdf in terms of x,mean and sigma
   RooGaussian gauss("gauss", "gaussian PDF", x, mean, sigma);

   // C r e a t e   a n d   p l o t  d e r i v a t i v e s   w . r . t .   x
   // ----------------------------------------------------------------------

   // Derivative of normalized gauss(x) w.r.t. observable x
   RooAbsReal *dgdx = gauss.derivative(x, 1);

   // Second and third derivative of normalized gauss(x) w.r.t. observable x
   RooAbsReal *d2gdx2 = gauss.derivative(x, 2);
   RooAbsReal *d3gdx3 = gauss.derivative(x, 3);

   // Construct plot frame in 'x'
   RooPlot *xframe = x.frame(Title("d(Gauss)/dx"));

   // Plot gauss in frame (i.e. in x)
   gauss.plotOn(xframe);

   // Plot derivatives in same frame
   dgdx->plotOn(xframe, LineColor(kMagenta));
   d2gdx2->plotOn(xframe, LineColor(kRed));
   d3gdx3->plotOn(xframe, LineColor(kOrange));

   // C r e a t e   a n d   p l o t  d e r i v a t i v e s   w . r . t .   s i g m a
   // ------------------------------------------------------------------------------

   // Derivative of normalized gauss(x) w.r.t. parameter sigma
   RooAbsReal *dgds = gauss.derivative(sigma, 1);

   // Second and third derivative of normalized gauss(x) w.r.t. parameter sigma
   RooAbsReal *d2gds2 = gauss.derivative(sigma, 2);
   RooAbsReal *d3gds3 = gauss.derivative(sigma, 3);

   // Construct plot frame in 'sigma'
   RooPlot *sframe = sigma.frame(Title("d(Gauss)/d(sigma)"), Range(0., 2.));

   // Plot gauss in frame (i.e. in x)
   gauss.plotOn(sframe);

   // Plot derivatives in same frame
   dgds->plotOn(sframe, LineColor(kMagenta));
   d2gds2->plotOn(sframe, LineColor(kRed));
   d3gds3->plotOn(sframe, LineColor(kOrange));

   // Draw all frames on a canvas
   TCanvas *c = new TCanvas("rf111_derivatives", "rf111_derivatives", 800, 400);
   c->Divide(2);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   xframe->GetYaxis()->SetTitleOffset(1.6);
   xframe->Draw();
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   sframe->GetYaxis()->SetTitleOffset(1.6);
   sframe->Draw();
}
