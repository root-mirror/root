/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Validation and MC studies: RooMCStudy - Using the randomizer and profile likelihood add-on models
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
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooMCStudy.h"
#include "RooRandomizeParamMCSModule.h"
#include "RooDLLSignificanceMCSModule.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TDirectory.h"

using namespace RooFit;

void rf803_mcstudy_addons2()
{
   // C r e a t e   m o d e l
   // -----------------------

   // Simulation of signal and background of top quark decaying into
   // 3 jets with background

   // Observable
   RooRealVar mjjj("mjjj", "m(3jet) (GeV)", 100, 85., 350.);

   // Signal component (Gaussian)
   RooRealVar mtop("mtop", "m(top)", 162);
   RooRealVar wtop("wtop", "m(top) resolution", 15.2);
   RooGaussian sig("sig", "top signal", mjjj, mtop, wtop);

   // Background component (Chebychev)
   RooRealVar c0("c0", "Chebychev coefficient 0", -0.846, -1., 1.);
   RooRealVar c1("c1", "Chebychev coefficient 1", 0.112, -1., 1.);
   RooRealVar c2("c2", "Chebychev coefficient 2", 0.076, -1., 1.);
   RooChebychev bkg("bkg", "combinatorial background", mjjj, RooArgList(c0, c1, c2));

   // Composite model
   RooRealVar nsig("nsig", "number of signal events", 53, 0, 1e3);
   RooRealVar nbkg("nbkg", "number of background events", 103, 0, 5e3);
   RooAddPdf model("model", "model", RooArgList(sig, bkg), RooArgList(nsig, nbkg));

   // C r e a t e   m a n a g e r
   // ---------------------------

   // Configure manager to perform binned extended likelihood fits (Binned(),Extended()) on data generated
   // with a Poisson fluctuation on Nobs (Extended())
   RooMCStudy *mcs = new RooMCStudy(model, mjjj, Binned(), Silence(), Extended(kTRUE),
                                    FitOptions(Extended(kTRUE), PrintEvalErrors(-1)));

   // C u s t o m i z e   m a n a g e r
   // ---------------------------------

   // Add module that randomizes the summed value of nsig+nbkg
   // sampling from a uniform distribution between 0 and 1000
   //
   // In general one can randomize a single parameter, or a
   // sum of N parameters, using either a uniform or a Gaussian
   // distribution. Multiple randomization can be executed
   // by a single randomizer module

   RooRandomizeParamMCSModule randModule;
   randModule.sampleSumUniform(RooArgSet(nsig, nbkg), 50, 500);
   mcs->addModule(randModule);

   // Add profile likelihood calculation of significance. Redo each
   // fit while keeping parameter nsig fixed to zero. For each toy,
   // the difference in -log(L) of both fits is stored, as well
   // a simple significance interpretation of the delta(-logL)
   // using Dnll = 0.5 sigma^2

   RooDLLSignificanceMCSModule sigModule(nsig, 0);
   mcs->addModule(sigModule);

   // R u n   m a n a g e r ,   m a k e   p l o t s
   // ---------------------------------------------

   // Run 1000 experiments. This configuration will generate a fair number
   // of (harmless) MINUIT warnings due to the instability of the Chebychev polynomial fit
   // at low statistics.
   mcs->generateAndFit(500);

   // Make some plots
   TH1 *dll_vs_ngen = mcs->fitParDataSet().createHistogram("ngen,dll_nullhypo_nsig", -40, -40);
   TH1 *z_vs_ngen = mcs->fitParDataSet().createHistogram("ngen,significance_nullhypo_nsig", -40, -40);
   TH1 *errnsig_vs_ngen = mcs->fitParDataSet().createHistogram("ngen,nsigerr", -40, -40);
   TH1 *errnsig_vs_nsig = mcs->fitParDataSet().createHistogram("nsig,nsigerr", -40, -40);

   // Draw plots on canvas
   TCanvas *c = new TCanvas("rf803_mcstudy_addons2", "rf802_mcstudy_addons2", 800, 800);
   c->Divide(2, 2);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   dll_vs_ngen->GetYaxis()->SetTitleOffset(1.6);
   dll_vs_ngen->Draw("box");
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   z_vs_ngen->GetYaxis()->SetTitleOffset(1.6);
   z_vs_ngen->Draw("box");
   c->cd(3);
   gPad->SetLeftMargin(0.15);
   errnsig_vs_ngen->GetYaxis()->SetTitleOffset(1.6);
   errnsig_vs_ngen->Draw("box");
   c->cd(4);
   gPad->SetLeftMargin(0.15);
   errnsig_vs_nsig->GetYaxis()->SetTitleOffset(1.6);
   errnsig_vs_nsig->Draw("box");

   // Make RooMCStudy object available on command line after
   // macro finishes
   gDirectory->Add(mcs);
}
