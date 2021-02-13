/// \file
/// \ingroup tutorial_roofit
/// \notebook
/// Addition and convolution: working with a pdf with a convolution operator in terms of a parameter
///
/// This tutorial requires FFT3 to be enabled.
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
#include "RooGenericPdf.h"
#include "RooFormulaVar.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH2.h"
using namespace RooFit;

void rf211_paramconv()
{
   // S e t u p   c o m p o n e n t   p d f s
   // ---------------------------------------

   // Gaussian g(x ; mean,sigma)
   RooRealVar x("x", "x", -10, 10);
   RooRealVar mean("mean", "mean", -3, 3);
   RooRealVar sigma("sigma", "sigma", 0.5, 0.1, 10);
   RooGaussian modelx("gx", "gx", x, mean, sigma);

   // Block function in mean
   RooRealVar a("a", "a", 2, 1, 10);
   RooGenericPdf model_mean("model_mean", "abs(mean)<a", RooArgList(mean, a));

   // Convolution in mean parameter model = g(x,mean,sigma) (x) block(mean)
   x.setBins(1000, "cache");
   mean.setBins(50, "cache");
   RooFFTConvPdf model("model", "model", mean, modelx, model_mean);

   // Configure convolution to construct a 2-D cache in (x,mean)
   // rather than a 1-d cache in mean that needs to be recalculated
   // for each value of x
   model.setCacheObservables(x);
   model.setBufferFraction(1.0);

   // Integrate model over mean projModel = Int model dmean
   RooAbsPdf *projModel = model.createProjection(mean);

   // Generate 1000 toy events
   RooDataHist *d = projModel->generateBinned(x, 1000);

   // Fit pdf to toy data
   projModel->fitTo(*d, Verbose());

   // Plot data and fitted pdf
   RooPlot *frame = x.frame(Bins(25));
   d->plotOn(frame);
   projModel->plotOn(frame);

   // Make 2d histogram of model(x;mean)
   TH1 *hh = model.createHistogram("hh", x, Binning(50), YVar(mean, Binning(50)), ConditionalObservables(mean));
   hh->SetTitle("histogram of model(x|mean)");
   hh->SetLineColor(kBlue);

   // Draw frame on canvas
   TCanvas *c = new TCanvas("rf211_paramconv", "rf211_paramconv", 800, 400);
   c->Divide(2);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   frame->GetYaxis()->SetTitleOffset(1.4);
   frame->Draw();
   c->cd(2);
   gPad->SetLeftMargin(0.20);
   hh->GetZaxis()->SetTitleOffset(2.5);
   hh->Draw("surf");
}
