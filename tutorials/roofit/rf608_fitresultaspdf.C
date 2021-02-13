/// \file
/// \ingroup tutorial_roofit
/// \notebook
/// Likelihood and minimization: representing the parabolic approximation of the fit as a multi-variate Gaussian on the
/// parameters of the fitted pdf
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
#include "RooChebychev.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH2.h"
#include "TH3.h"

using namespace RooFit;

void rf608_fitresultaspdf()
{
   // C r e a t e   m o d e l   a n d   d a t a s e t
   // -----------------------------------------------

   // Observable
   RooRealVar x("x", "x", -20, 20);

   // Model (intentional strong correlations)
   RooRealVar mean("mean", "mean of g1 and g2", 0, -1, 1);
   RooRealVar sigma_g1("sigma_g1", "width of g1", 2);
   RooGaussian g1("g1", "g1", x, mean, sigma_g1);

   RooRealVar sigma_g2("sigma_g2", "width of g2", 4, 3.0, 5.0);
   RooGaussian g2("g2", "g2", x, mean, sigma_g2);

   RooRealVar frac("frac", "frac", 0.5, 0.0, 1.0);
   RooAddPdf model("model", "model", RooArgList(g1, g2), frac);

   // Generate 1000 events
   RooDataSet *data = model.generate(x, 1000);

   // F i t   m o d e l   t o   d a t a
   // ----------------------------------

   RooFitResult *r = model.fitTo(*data, Save());

   // C r e a t e M V   G a u s s i a n   p d f   o f   f i t t e d    p a r a m e t e r s
   // ------------------------------------------------------------------------------------

   RooAbsPdf *parabPdf = r->createHessePdf(RooArgSet(frac, mean, sigma_g2));

   // S o m e   e x e c e r c i s e s   w i t h   t h e   p a r a m e t e r   p d f
   // -----------------------------------------------------------------------------

   // Generate 100K points in the parameter space, sampled from the MVGaussian pdf
   RooDataSet *d = parabPdf->generate(RooArgSet(mean, sigma_g2, frac), 100000);

   // Sample a 3-D histogram of the pdf to be visualized as an error ellipsoid using the GLISO draw option
   TH3 *hh_3d = (TH3 *)parabPdf->createHistogram("mean,sigma_g2,frac", 25, 25, 25);
   hh_3d->SetFillColor(kBlue);

   // Project 3D parameter pdf down to 3 permutations of two-dimensional pdfs
   // The integrations corresponding to these projections are performed analytically
   // by the MV Gaussian pdf
   RooAbsPdf *pdf_sigmag2_frac = parabPdf->createProjection(mean);
   RooAbsPdf *pdf_mean_frac = parabPdf->createProjection(sigma_g2);
   RooAbsPdf *pdf_mean_sigmag2 = parabPdf->createProjection(frac);

   // Make 2D plots of the 3 two-dimensional pdf projections
   TH2 *hh_sigmag2_frac = (TH2 *)pdf_sigmag2_frac->createHistogram("sigma_g2,frac", 50, 50);
   TH2 *hh_mean_frac = (TH2 *)pdf_mean_frac->createHistogram("mean,frac", 50, 50);
   TH2 *hh_mean_sigmag2 = (TH2 *)pdf_mean_sigmag2->createHistogram("mean,sigma_g2", 50, 50);
   hh_mean_frac->SetLineColor(kBlue);
   hh_sigmag2_frac->SetLineColor(kBlue);
   hh_mean_sigmag2->SetLineColor(kBlue);

   // Draw the 'sigar'
   new TCanvas("rf608_fitresultaspdf_1", "rf608_fitresultaspdf_1", 600, 600);
   hh_3d->Draw("iso");

   // Draw the 2D projections of the 3D pdf
   TCanvas *c2 = new TCanvas("rf608_fitresultaspdf_2", "rf608_fitresultaspdf_2", 900, 600);
   c2->Divide(3, 2);
   c2->cd(1);
   gPad->SetLeftMargin(0.15);
   hh_mean_sigmag2->GetZaxis()->SetTitleOffset(1.4);
   hh_mean_sigmag2->Draw("surf3");
   c2->cd(2);
   gPad->SetLeftMargin(0.15);
   hh_sigmag2_frac->GetZaxis()->SetTitleOffset(1.4);
   hh_sigmag2_frac->Draw("surf3");
   c2->cd(3);
   gPad->SetLeftMargin(0.15);
   hh_mean_frac->GetZaxis()->SetTitleOffset(1.4);
   hh_mean_frac->Draw("surf3");

   // Draw the distributions of parameter points sampled from the pdf
   TH1 *tmp1 = d->createHistogram("mean,sigma_g2", 50, 50);
   TH1 *tmp2 = d->createHistogram("sigma_g2,frac", 50, 50);
   TH1 *tmp3 = d->createHistogram("mean,frac", 50, 50);

   c2->cd(4);
   gPad->SetLeftMargin(0.15);
   tmp1->GetZaxis()->SetTitleOffset(1.4);
   tmp1->Draw("lego3");
   c2->cd(5);
   gPad->SetLeftMargin(0.15);
   tmp2->GetZaxis()->SetTitleOffset(1.4);
   tmp2->Draw("lego3");
   c2->cd(6);
   gPad->SetLeftMargin(0.15);
   tmp3->GetZaxis()->SetTitleOffset(1.4);
   tmp3->Draw("lego3");
}
