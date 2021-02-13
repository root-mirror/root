/// \file
/// \ingroup tutorial_roofit
/// \notebook
/// Special pdf's: using a pdf defined by a sum of real-valued amplitude components
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
#include "RooTruthModel.h"
#include "RooFormulaVar.h"
#include "RooRealSumPdf.h"
#include "RooPolyVar.h"
#include "RooProduct.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit;

void rf704_amplitudefit()
{
   // S e t u p   2 D   a m p l i t u d e   f u n c t i o n s
   // -------------------------------------------------------

   // Observables
   RooRealVar t("t", "time", -1., 15.);
   RooRealVar cosa("cosa", "cos(alpha)", -1., 1.);

   // Use RooTruthModel to obtain compiled implementation of sinh/cosh modulated decay functions
   RooRealVar tau("tau", "#tau", 1.5);
   RooRealVar deltaGamma("deltaGamma", "deltaGamma", 0.3);
   RooTruthModel truthModel("tm", "tm", t);
   RooFormulaVar coshGBasis("coshGBasis", "exp(-@0/ @1)*cosh(@0*@2/2)", RooArgList(t, tau, deltaGamma));
   RooFormulaVar sinhGBasis("sinhGBasis", "exp(-@0/ @1)*sinh(@0*@2/2)", RooArgList(t, tau, deltaGamma));
   RooAbsReal *coshGConv = truthModel.convolution(&coshGBasis, &t);
   RooAbsReal *sinhGConv = truthModel.convolution(&sinhGBasis, &t);

   // Construct polynomial amplitudes in cos(a)
   RooPolyVar poly1("poly1", "poly1", cosa, RooArgList(RooConst(0.5), RooConst(0.2), RooConst(0.2)), 0);
   RooPolyVar poly2("poly2", "poly2", cosa, RooArgList(RooConst(1), RooConst(-0.2), RooConst(3)), 0);

   // Construct 2D amplitude as uncorrelated product of amp(t)*amp(cosa)
   RooProduct ampl1("ampl1", "amplitude 1", RooArgSet(poly1, *coshGConv));
   RooProduct ampl2("ampl2", "amplitude 2", RooArgSet(poly2, *sinhGConv));

   // C o n s t r u c t   a m p l i t u d e   s u m   p d f
   // -----------------------------------------------------

   // Amplitude strengths
   RooRealVar f1("f1", "f1", 1, 0, 2);
   RooRealVar f2("f2", "f2", 0.5, 0, 2);

   // Construct pdf
   RooRealSumPdf pdf("pdf", "pdf", RooArgList(ampl1, ampl2), RooArgList(f1, f2));

   // Generate some toy data from pdf
   RooDataSet *data = pdf.generate(RooArgSet(t, cosa), 10000);

   // Fit pdf to toy data with only amplitude strength floating
   pdf.fitTo(*data);

   // P l o t   a m p l i t u d e   s u m   p d f
   // -------------------------------------------

   // Make 2D plots of amplitudes
   TH1 *hh_cos = ampl1.createHistogram("hh_cos", t, Binning(50), YVar(cosa, Binning(50)));
   TH1 *hh_sin = ampl2.createHistogram("hh_sin", t, Binning(50), YVar(cosa, Binning(50)));
   hh_cos->SetLineColor(kBlue);
   hh_sin->SetLineColor(kRed);

   // Make projection on t, plot data, pdf and its components
   // Note component projections may be larger than sum because amplitudes can be negative
   RooPlot *frame1 = t.frame();
   data->plotOn(frame1);
   pdf.plotOn(frame1);
   pdf.plotOn(frame1, Components(ampl1), LineStyle(kDashed));
   pdf.plotOn(frame1, Components(ampl2), LineStyle(kDashed), LineColor(kRed));

   // Make projection on cosa, plot data, pdf and its components
   // Note that components projection may be larger than sum because amplitudes can be negative
   RooPlot *frame2 = cosa.frame();
   data->plotOn(frame2);
   pdf.plotOn(frame2);
   pdf.plotOn(frame2, Components(ampl1), LineStyle(kDashed));
   pdf.plotOn(frame2, Components(ampl2), LineStyle(kDashed), LineColor(kRed));

   TCanvas *c = new TCanvas("rf704_amplitudefit", "rf704_amplitudefit", 800, 800);
   c->Divide(2, 2);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   frame1->GetYaxis()->SetTitleOffset(1.4);
   frame1->Draw();
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   frame2->GetYaxis()->SetTitleOffset(1.4);
   frame2->Draw();
   c->cd(3);
   gPad->SetLeftMargin(0.20);
   hh_cos->GetZaxis()->SetTitleOffset(2.3);
   hh_cos->Draw("surf");
   c->cd(4);
   gPad->SetLeftMargin(0.20);
   hh_sin->GetZaxis()->SetTitleOffset(2.3);
   hh_sin->Draw("surf");
}
