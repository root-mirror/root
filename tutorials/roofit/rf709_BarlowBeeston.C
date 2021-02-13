/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Implementing the Barlow-Beeston method for taking into account the statistical
/// uncertainty of a Monte-Carlo fit template.
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// Based on a demo by Wouter Verkerke
//
/// \date June 2019
/// \author Stephan Hageboeck, CERN

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooUniform.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistFunc.h"
#include "RooRealSumPdf.h"
#include "RooParamHistFunc.h"
#include "RooHistConstraint.h"
#include "RooProdPdf.h"
#include "RooPlot.h"

#include "TCanvas.h"
#include "TPaveText.h"

#include <iostream>
#include <memory>

using namespace RooFit;

void rf709_BarlowBeeston()
{
  // First, construct a likelihood model with a Gaussian signal on top of a uniform background
  RooRealVar x("x", "x", -20, 20);
  x.setBins(25);

  RooRealVar meanG("meanG", "meanG", 1, -10, 10);
  RooRealVar sigG("sigG", "sigG", 1.5, -10, 10);
  RooGaussian g("g", "Gauss", x, meanG, sigG);
  RooUniform u("u", "Uniform", x);


  // Generate the data to be fitted
  std::unique_ptr<RooDataSet> sigData(g.generate(x, 50));
  std::unique_ptr<RooDataSet> bkgData(u.generate(x, 1000));

  RooDataSet sumData("sumData", "Gauss + Uniform", x);
  sumData.append(*sigData);
  sumData.append(*bkgData);


  // Make histogram templates for signal and background.
  // Let's take a signal distribution with low statistics and a more accurate
  // background distribution.
  // Normally, these come from Monte Carlo simulations, but we will just generate them.
  std::unique_ptr<RooDataHist> dh_sig( g.generateBinned(x, 50) );
  std::unique_ptr<RooDataHist> dh_bkg( u.generateBinned(x, 10000) );


  // ***** Case 0 - 'Rigid templates' *****

  // Construct histogram shapes for signal and background
  RooHistFunc p_h_sig("p_h_sig","p_h_sig",x,*dh_sig);
  RooHistFunc p_h_bkg("p_h_bkg","p_h_bkg",x,*dh_bkg);

  // Construct scale factors for adding the two distributions
  RooRealVar Asig0("Asig","Asig",1,0.01,5000);
  RooRealVar Abkg0("Abkg","Abkg",1,0.01,5000);

  // Construct the sum model
  RooRealSumPdf model0("model0","model0",
      RooArgList(p_h_sig,p_h_bkg),
      RooArgList(Asig0,Abkg0),
      true);



  // ***** Case 1 - 'Barlow Beeston' *****

  // Construct parameterized histogram shapes for signal and background
  RooParamHistFunc p_ph_sig1("p_ph_sig","p_ph_sig",*dh_sig);
  RooParamHistFunc p_ph_bkg1("p_ph_bkg","p_ph_bkg",*dh_bkg);

  RooRealVar Asig1("Asig","Asig",1,0.01,5000);
  RooRealVar Abkg1("Abkg","Abkg",1,0.01,5000);

  // Construct the sum of these
  RooRealSumPdf model_tmp("sp_ph", "sp_ph",
      RooArgList(p_ph_sig1,p_ph_bkg1),
      RooArgList(Asig1,Abkg1),
      true);

  // Construct the subsidiary poisson measurements constraining the histogram parameters
  // These ensure that the bin contents of the histograms are only allowed to vary within
  // the statistical uncertainty of the Monte Carlo.
  RooHistConstraint hc_sig("hc_sig","hc_sig",p_ph_sig1);
  RooHistConstraint hc_bkg("hc_bkg","hc_bkg",p_ph_bkg1);

  // Construct the joint model with template PDFs and constraints
  RooProdPdf model1("model1","model1",RooArgSet(hc_sig,hc_bkg),Conditional(model_tmp,x));



  // ***** Case 2 - 'Barlow Beeston' light (one parameter per bin for all samples) *****

  // Construct the histogram shapes, using the same parameters for signal and background
  // This requires passing the first histogram to the second, so that their common parameters
  // can be re-used.
  // The first ParamHistFunc will create one parameter per bin, such as `p_ph_sig2_gamma_bin_0`.
  // This allows bin 0 to fluctuate up and down.
  // Then, the SAME parameters are connected to the background histogram, so the bins flucutate
  // synchronously. This reduces the number of parameters.
  RooParamHistFunc p_ph_sig2("p_ph_sig2", "p_ph_sig2", *dh_sig);
  RooParamHistFunc p_ph_bkg2("p_ph_bkg2", "p_ph_bkg2", *dh_bkg, p_ph_sig2, true);

  RooRealVar Asig2("Asig","Asig",1,0.01,5000);
  RooRealVar Abkg2("Abkg","Abkg",1,0.01,5000);

  // As before, construct the sum of signal2 and background2
  RooRealSumPdf model2_tmp("sp_ph","sp_ph",
      RooArgList(p_ph_sig2,p_ph_bkg2),
      RooArgList(Asig2,Abkg2),
      true);

  // Construct the subsidiary poisson measurements constraining the statistical fluctuations
  RooHistConstraint hc_sigbkg("hc_sigbkg","hc_sigbkg",RooArgSet(p_ph_sig2,p_ph_bkg2));

  // Construct the joint model
  RooProdPdf model2("model2","model2",hc_sigbkg, RooFit::Conditional(model2_tmp,x));



  // ************ Fit all models to data and plot *********************

  auto result0 = model0.fitTo(sumData, PrintLevel(0), Save());
  auto result1 = model1.fitTo(sumData, PrintLevel(0), Save());
  auto result2 = model2.fitTo(sumData, PrintLevel(0), Save());


  TCanvas* can = new TCanvas("can", "", 1500, 600);
  can->Divide(3,1);

  TPaveText pt(-19.5, 1, -2, 25);
  pt.SetFillStyle(0);
  pt.SetBorderSize(0);


  can->cd(1);
  auto frame = x.frame(Title("No template uncertainties"));
  // Plot data to enable automatic determination of model0 normalisation:
  sumData.plotOn(frame);
  model0.plotOn(frame, LineColor(kBlue), VisualizeError(*result0));
  // Plot data again to show it on top of model0 error bands:
  sumData.plotOn(frame);
  // Plot model components
  model0.plotOn(frame, LineColor(kBlue));
  model0.plotOn(frame, Components(p_h_sig), LineColor(kAzure));
  model0.plotOn(frame, Components(p_h_bkg), LineColor(kRed));
  model0.paramOn(frame);

  sigData->plotOn(frame, MarkerColor(kBlue));
  frame->Draw();

  for (auto text : {
    "No template uncertainties",
    "are taken into account.",
    "This leads to low errors",
    "for the parameters A, since",
    "the only source of errors",
    "are the data statistics."}) {
    pt.AddText(text);
  }
  pt.DrawClone();


  can->cd(2);
  frame = x.frame(Title("Barlow Beeston for Sig & Bkg separately"));
  sumData.plotOn(frame);
  model1.plotOn(frame, LineColor(kBlue), VisualizeError(*result1));
  // Plot data again to show it on top of error bands:
  sumData.plotOn(frame);
  model1.plotOn(frame, LineColor(kBlue));
  model1.plotOn(frame, Components(p_ph_sig1), LineColor(kAzure));
  model1.plotOn(frame, Components(p_ph_bkg1), LineColor(kRed));
  model1.paramOn(frame, Parameters(RooArgSet(Asig1, Abkg1)));

  sigData->plotOn(frame, MarkerColor(kBlue));
  frame->Draw();

  pt.Clear();
  for (auto text : {
    "With gamma parameters, the",
    "signal & background templates",
    "can adapt to the data.",
    "Note how the blue signal",
    "template changes its shape.",
    "This leads to higher errors",
    "of the scale parameters A."}) {
    pt.AddText(text);
  }
  pt.DrawClone();

  can->cd(3);
  frame = x.frame(Title("Barlow Beeston light for (Sig+Bkg)"));
  sumData.plotOn(frame);
  model2.plotOn(frame, LineColor(kBlue), VisualizeError(*result2));
  // Plot data again to show it on top of model0 error bands:
  sumData.plotOn(frame);
  model2.plotOn(frame, LineColor(kBlue));
  model2.plotOn(frame, Components(p_ph_sig2), LineColor(kAzure));
  model2.plotOn(frame, Components(p_ph_bkg2), LineColor(kRed));
  model2.paramOn(frame, Parameters(RooArgSet(Asig2, Abkg2)));

  sigData->plotOn(frame, MarkerColor(kBlue));
  frame->Draw();

  pt.Clear();
  for (auto text : {
    "When signal and background",
    "template share one gamma para-",
    "meter per bin, they adapt less.",
    "The errors of the A parameters",
    "also shrink slightly."}) {
    pt.AddText(text);
  }
  pt.DrawClone();


  std::cout << "Asig [normal ] = " << Asig0.getVal() << " +/- " << Asig0.getError() << std::endl;
  std::cout << "Asig [BB     ] = " << Asig1.getVal() << " +/- " << Asig1.getError() << std::endl;
  std::cout << "Asig [BBlight] = " << Asig2.getVal() << " +/- " << Asig2.getError() << std::endl;

}
