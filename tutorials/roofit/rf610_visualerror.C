/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Likelihood and minimization: visualization of errors from a covariance matrix
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \date April 2009
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TAxis.h"
using namespace RooFit;

void rf610_visualerror()
{
   // S e t u p   e x a m p l e   f i t
   // ---------------------------------------

   // Create sum of two Gaussians pdf with factory
   RooRealVar x("x", "x", -10, 10);

   RooRealVar m("m", "m", 0, -10, 10);
   RooRealVar s("s", "s", 2, 1, 50);
   RooGaussian sig("sig", "sig", x, m, s);

   RooRealVar m2("m2", "m2", -1, -10, 10);
   RooRealVar s2("s2", "s2", 6, 1, 50);
   RooGaussian bkg("bkg", "bkg", x, m2, s2);

   RooRealVar fsig("fsig", "fsig", 0.33, 0, 1);
   RooAddPdf model("model", "model", RooArgList(sig, bkg), fsig);

   // Create binned dataset
   x.setBins(25);
   RooAbsData *d = model.generateBinned(x, 1000);

   // Perform fit and save fit result
   RooFitResult *r = model.fitTo(*d, Save());

   // V i s u a l i z e   f i t   e r r o r
   // -------------------------------------

   // Make plot frame
   RooPlot *frame = x.frame(Bins(40), Title("P.d.f with visualized 1-sigma error band"));
   d->plotOn(frame);

   // Visualize 1-sigma error encoded in fit result 'r' as orange band using linear error propagation
   // This results in an error band that is by construction symmetric
   //
   // The linear error is calculated as
   // error(x) = Z* F_a(x) * Corr(a,a') F_a'(x)
   //
   // where     F_a(x) = [ f(x,a+da) - f(x,a-da) ] / 2,
   //
   //         with f(x) = the plotted curve
   //              'da' = error taken from the fit result
   //        Corr(a,a') = the correlation matrix from the fit result
   //                Z = requested significance 'Z sigma band'
   //
   // The linear method is fast (required 2*N evaluations of the curve, where N is the number of parameters),
   // but may not be accurate in the presence of strong correlations (~>0.9) and at Z>2 due to linear and
   // Gaussian approximations made
   //
   model.plotOn(frame, VisualizeError(*r, 1), FillColor(kOrange));

   // Calculate error using sampling method and visualize as dashed red line.
   //
   // In this method a number of curves is calculated with variations of the parameter values, as sampled
   // from a multi-variate Gaussian pdf that is constructed from the fit results covariance matrix.
   // The error(x) is determined by calculating a central interval that capture N% of the variations
   // for each value of x, where N% is controlled by Z (i.e. Z=1 gives N=68%). The number of sampling curves
   // is chosen to be such that at least 100 curves are expected to be outside the N% interval, and is minimally
   // 100 (e.g. Z=1->Ncurve=356, Z=2->Ncurve=2156)) Intervals from the sampling method can be asymmetric,
   // and may perform better in the presence of strong correlations, but may take (much) longer to calculate
   model.plotOn(frame, VisualizeError(*r, 1, kFALSE), DrawOption("L"), LineWidth(2), LineColor(kRed));

   // Perform the same type of error visualization on the background component only.
   // The VisualizeError() option can generally applied to _any_ kind of plot (components, asymmetries, efficiencies
   // etc..)
   model.plotOn(frame, VisualizeError(*r, 1), FillColor(kOrange), Components("bkg"));
   model.plotOn(frame, VisualizeError(*r, 1, kFALSE), DrawOption("L"), LineWidth(2), LineColor(kRed), Components("bkg"),
                LineStyle(kDashed));

   // Overlay central value
   model.plotOn(frame);
   model.plotOn(frame, Components("bkg"), LineStyle(kDashed));
   d->plotOn(frame);
   frame->SetMinimum(0);

   // V i s u a l i z e   p a r t i a l   f i t   e r r o r
   // ------------------------------------------------------

   // Make plot frame
   RooPlot *frame2 = x.frame(Bins(40), Title("Visualization of 2-sigma partial error from (m,m2)"));

   // Visualize partial error. For partial error visualization the covariance matrix is first reduced as follows
   //        ___                   -1
   // Vred = V22  = V11 - V12 * V22   * V21
   //
   // Where V11,V12,V21,V22 represent a block decomposition of the covariance matrix into observables that
   // are propagated (labeled by index '1') and that are not propagated (labeled by index '2'), and V22bar
   // is the Shur complement of V22, calculated as shown above
   //
   // (Note that Vred is _not_ a simple sub-matrix of V)

   // Propagate partial error due to shape parameters (m,m2) using linear and sampling method
   model.plotOn(frame2, VisualizeError(*r, RooArgSet(m, m2), 2), FillColor(kCyan));
   model.plotOn(frame2, Components("bkg"), VisualizeError(*r, RooArgSet(m, m2), 2), FillColor(kCyan));

   model.plotOn(frame2);
   model.plotOn(frame2, Components("bkg"), LineStyle(kDashed));
   frame2->SetMinimum(0);

   // Make plot frame
   RooPlot *frame3 = x.frame(Bins(40), Title("Visualization of 2-sigma partial error from (s,s2)"));

   // Propagate partial error due to yield parameter using linear and sampling method
   model.plotOn(frame3, VisualizeError(*r, RooArgSet(s, s2), 2), FillColor(kGreen));
   model.plotOn(frame3, Components("bkg"), VisualizeError(*r, RooArgSet(s, s2), 2), FillColor(kGreen));

   model.plotOn(frame3);
   model.plotOn(frame3, Components("bkg"), LineStyle(kDashed));
   frame3->SetMinimum(0);

   // Make plot frame
   RooPlot *frame4 = x.frame(Bins(40), Title("Visualization of 2-sigma partial error from fsig"));

   // Propagate partial error due to yield parameter using linear and sampling method
   model.plotOn(frame4, VisualizeError(*r, RooArgSet(fsig), 2), FillColor(kMagenta));
   model.plotOn(frame4, Components("bkg"), VisualizeError(*r, RooArgSet(fsig), 2), FillColor(kMagenta));

   model.plotOn(frame4);
   model.plotOn(frame4, Components("bkg"), LineStyle(kDashed));
   frame4->SetMinimum(0);

   TCanvas *c = new TCanvas("rf610_visualerror", "rf610_visualerror", 800, 800);
   c->Divide(2, 2);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   frame->GetYaxis()->SetTitleOffset(1.4);
   frame->Draw();
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
