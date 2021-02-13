/// \file
/// \ingroup tutorial_roofit
/// \notebook
/// Likelihood and minimization: demonstration of options of the RooFitResult class
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
#include "TMatrixDSym.h"

using namespace RooFit;

void rf607_fitresult()
{
   // C r e a t e   p d f ,   d a t a
   // --------------------------------

   // Declare observable x
   RooRealVar x("x", "x", 0, 10);

   // Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their parameters
   RooRealVar mean("mean", "mean of gaussians", 5, -10, 10);
   RooRealVar sigma1("sigma1", "width of gaussians", 0.5, 0.1, 10);
   RooRealVar sigma2("sigma2", "width of gaussians", 1, 0.1, 10);

   RooGaussian sig1("sig1", "Signal component 1", x, mean, sigma1);
   RooGaussian sig2("sig2", "Signal component 2", x, mean, sigma2);

   // Build Chebychev polynomial pdf
   RooRealVar a0("a0", "a0", 0.5, 0., 1.);
   RooRealVar a1("a1", "a1", -0.2);
   RooChebychev bkg("bkg", "Background", x, RooArgSet(a0, a1));

   // Sum the signal components into a composite signal pdf
   RooRealVar sig1frac("sig1frac", "fraction of component 1 in signal", 0.8, 0., 1.);
   RooAddPdf sig("sig", "Signal", RooArgList(sig1, sig2), sig1frac);

   // Sum the composite signal and background
   RooRealVar bkgfrac("bkgfrac", "fraction of background", 0.5, 0., 1.);
   RooAddPdf model("model", "g1+g2+a", RooArgList(bkg, sig), bkgfrac);

   // Generate 1000 events
   RooDataSet *data = model.generate(x, 1000);

   // F i t   p d f   t o   d a t a ,   s a v e   f i t r e s u l t
   // -------------------------------------------------------------

   // Perform fit and save result
   RooFitResult *r = model.fitTo(*data, Save());

   // P r i n t   f i t   r e s u l t s
   // ---------------------------------

   // Summary printing: Basic info plus final values of floating fit parameters
   r->Print();

   // Verbose printing: Basic info, values of constant parameters, initial and
   // final values of floating parameters, global correlations
   r->Print("v");

   // V i s u a l i z e   c o r r e l a t i o n   m a t r i x
   // -------------------------------------------------------

   // Construct 2D color plot of correlation matrix
   gStyle->SetOptStat(0);
   TH2 *hcorr = r->correlationHist();

   // Visualize ellipse corresponding to single correlation matrix element
   RooPlot *frame = new RooPlot(sigma1, sig1frac, 0.45, 0.60, 0.65, 0.90);
   frame->SetTitle("Covariance between sigma1 and sig1frac");
   r->plotOn(frame, sigma1, sig1frac, "ME12ABHV");

   // A c c e s s   f i t   r e s u l t   i n f o r m a t i o n
   // ---------------------------------------------------------

   // Access basic information
   cout << "EDM = " << r->edm() << endl;
   cout << "-log(L) at minimum = " << r->minNll() << endl;

   // Access list of final fit parameter values
   cout << "final value of floating parameters" << endl;
   r->floatParsFinal().Print("s");

   // Access correlation matrix elements
   cout << "correlation between sig1frac and a0 is  " << r->correlation(sig1frac, a0) << endl;
   cout << "correlation between bkgfrac and mean is " << r->correlation("bkgfrac", "mean") << endl;

   // Extract covariance and correlation matrix as TMatrixDSym
   const TMatrixDSym &cor = r->correlationMatrix();
   const TMatrixDSym &cov = r->covarianceMatrix();

   // Print correlation, covariance matrix
   cout << "correlation matrix" << endl;
   cor.Print();
   cout << "covariance matrix" << endl;
   cov.Print();

   // P e r s i s t   f i t   r e s u l t   i n   r o o t   f i l e
   // -------------------------------------------------------------

   // Open new ROOT file save save result
   TFile f("rf607_fitresult.root", "RECREATE");
   r->Write("rf607");
   f.Close();

   // In a clean ROOT session retrieve the persisted fit result as follows:
   // RooFitResult* r = gDirectory->Get("rf607") ;

   TCanvas *c = new TCanvas("rf607_fitresult", "rf607_fitresult", 800, 400);
   c->Divide(2);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   hcorr->GetYaxis()->SetTitleOffset(1.4);
   hcorr->Draw("colz");
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   frame->GetYaxis()->SetTitleOffset(1.6);
   frame->Draw();
}
