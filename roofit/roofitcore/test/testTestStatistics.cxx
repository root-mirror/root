/*
 * testNLLVar.cxx
 *
 *  Created on: 8 Oct 2020
 *      Author: Stephan Hageboeck
 */

#include <RooRealVar.h>
#include <RooGenericPdf.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooBinning.h>
#include <RooPlot.h>
#include <RooRandom.h>

#include <gtest/gtest.h>

#include <memory>

TEST(RooNLLVar, IntegrateBins) {
  RooRandom::randomGenerator()->SetSeed(1337ul);

  RooRealVar x("x", "x", 0.1, 5.1);
  x.setBins(10);

  RooRealVar a("a", "a", -0.3, -5., 5.);
  RooArgSet targetValues;
  RooArgSet(a).snapshot(targetValues);

  RooGenericPdf pdf("pow", "std::pow(x, a)", RooArgSet(x, a));
  std::unique_ptr<RooDataHist> dataH(pdf.generateBinned(x,  10000));
  RooRealVar w("w", "weight", 0., 0., 10000.);
  RooDataSet data("data", "data", RooArgSet(x, w), RooFit::WeightVar(w));
  for (int i=0; i < dataH->numEntries(); ++i) {
    auto coords = dataH->get(i);
    data.add(*coords, dataH->weight());
  }

  std::unique_ptr<RooPlot> frame( x.frame() );
  dataH->plotOn(frame.get(), RooFit::MarkerColor(kRed));
  data.plotOn(frame.get(), RooFit::Name("data"));


  a.setVal(3.);
  std::unique_ptr<RooFitResult> fit1( pdf.fitTo(data, RooFit::Save(), RooFit::PrintLevel(-1)) );
  pdf.plotOn(frame.get(), RooFit::LineColor(kRed), RooFit::Name("standard"));

  a.setVal(3.);
  std::unique_ptr<RooFitResult> fit2( pdf.fitTo(data, RooFit::Save(), RooFit::PrintLevel(-1),
      RooFit::BatchMode(true),
      RooFit::IntegrateBins(1.E-3)) );
  pdf.plotOn(frame.get(), RooFit::LineColor(kBlue), RooFit::Name("highRes"));


  auto getVal = [](const char* name, const RooArgSet& set) {
    return dynamic_cast<const RooRealVar&>(set[name]).getVal();
  };
  auto getErr = [](const char* name, const RooArgSet& set) {
    return dynamic_cast<const RooRealVar&>(set[name]).getError();
  };

  EXPECT_GT(fabs(getVal("a", targetValues) - getVal("a", fit1->floatParsFinal())), 1. * getErr("a", fit1->floatParsFinal()))
      << "Expecting a bias when sampling PDF in bin centre.";

  EXPECT_NEAR(getVal("a", targetValues), getVal("a", fit2->floatParsFinal()), 1. * getErr("a", fit2->floatParsFinal()))
      << "Expect reduced bias with high-resolution sampling.";

  EXPECT_GT(frame->chiSquare("standard", "data", 1) * 0.9, frame->chiSquare("highRes",  "data", 1))
      << "Expect chi2/ndf at least 10% better.";

//  fit1->Print();
//  fit2->Print();
}


/// Prepare a RooDataSet that looks like the one that HistFactory uses:
/// It pretends to be an unbinned dataset, but instead of single events,
/// events are aggregated in the bin centres using weights.
TEST(RooNLLVar, IntegrateBins_SubRange) {
  RooRandom::randomGenerator()->SetSeed(1337ul);

  RooRealVar x("x", "x", 0.1, 5.1);
  x.setBins(10);
  x.setRange("range", 0.1, 4.1);
  x.setBins(8, "range"); // consistent binning

  RooRealVar a("a", "a", -0.7, -5., 5.);
  RooArgSet targetValues;
  RooArgSet(a).snapshot(targetValues);

  RooGenericPdf pdf("pow", "std::pow(x, a)", RooArgSet(x, a));
  std::unique_ptr<RooDataHist> dataH(pdf.generateBinned(x,  10000));
  RooRealVar w("w", "weight", 0., 0., 10000.);
  RooDataSet data("data", "data", RooArgSet(x, w), RooFit::WeightVar(w));
  for (int i=0; i < dataH->numEntries(); ++i) {
    auto coords = dataH->get(i);
    data.add(*coords, dataH->weight());
  }

  std::unique_ptr<RooPlot> frame( x.frame() );
  dataH->plotOn(frame.get(), RooFit::MarkerColor(kRed));
  data.plotOn(frame.get(), RooFit::Name("data"));


  a.setVal(3.);
  std::unique_ptr<RooFitResult> fit1( pdf.fitTo(data, RooFit::Save(), RooFit::PrintLevel(-1),
      RooFit::Optimize(0),
      RooFit::Range("range"),
      RooFit::BatchMode(true))  );
  pdf.plotOn(frame.get(), RooFit::LineColor(kRed), RooFit::Name("standard"));

  a.setVal(3.);
  std::unique_ptr<RooFitResult> fit2( pdf.fitTo(data, RooFit::Save(), RooFit::PrintLevel(-1),
      RooFit::Optimize(0),
      RooFit::Range("range"),
      RooFit::BatchMode(true),
      RooFit::IntegrateBins(1.E-3)) );
  pdf.plotOn(frame.get(), RooFit::LineColor(kBlue), RooFit::Name("highRes"));


  auto getVal = [](const char* name, const RooArgSet& set) {
    return dynamic_cast<const RooRealVar&>(set[name]).getVal();
  };
  auto getErr = [](const char* name, const RooArgSet& set) {
    return dynamic_cast<const RooRealVar&>(set[name]).getError();
  };

  EXPECT_GT(fabs(getVal("a", targetValues) - getVal("a", fit1->floatParsFinal())), 1. * getErr("a", fit1->floatParsFinal()))
      << "Expecting a bias when sampling PDF in bin centre.";

  EXPECT_NEAR(getVal("a", targetValues), getVal("a", fit2->floatParsFinal()), 1. * getErr("a", fit2->floatParsFinal()))
      << "Expect reduced bias with high-resolution sampling.";

  EXPECT_GT(frame->chiSquare("standard", "data", 1) * 0.9, frame->chiSquare("highRes",  "data", 1))
      << "Expect chi2/ndf at least 10% better.";

//  fit1->Print();
//  fit2->Print();
}

/// Prepare a RooDataSet that looks like the one that HistFactory uses:
/// It pretends to be an unbinned dataset, but instead of single events,
/// events are aggregated in the bin centres using weights.
TEST(RooNLLVar, IntegrateBins_CustomBinning) {
  RooRandom::randomGenerator()->SetSeed(1337ul);

  RooRealVar x("x", "x", 1., 5.);
  RooBinning binning(1., 5.);
  binning.addBoundary(1.5);
  binning.addBoundary(2.0);
  binning.addBoundary(3.);
  binning.addBoundary(4.);
  x.setBinning(binning);

  RooRealVar a("a", "a", -0.3, -5., 5.);
  RooArgSet targetValues;
  RooArgSet(a).snapshot(targetValues);

  RooGenericPdf pdf("pow", "std::pow(x, a)", RooArgSet(x, a));
  std::unique_ptr<RooDataHist> dataH(pdf.generateBinned(x,  50000));
  RooRealVar w("w", "weight", 0., 0., 1000000.);
  RooDataSet data("data", "data", RooArgSet(x, w), RooFit::WeightVar(w));
  for (int i=0; i < dataH->numEntries(); ++i) {
    auto coords = dataH->get(i);
    data.add(*coords, dataH->weight());
  }

  std::unique_ptr<RooPlot> frame( x.frame() );
  dataH->plotOn(frame.get(), RooFit::Name("dataHist"), RooFit::MarkerColor(kRed));
  data.plotOn(frame.get(), RooFit::Name("data"));


  a.setVal(3.);
  std::unique_ptr<RooFitResult> fit1( pdf.fitTo(data, RooFit::Save(), RooFit::PrintLevel(-1),
      RooFit::Optimize(0)) );
  pdf.plotOn(frame.get(), RooFit::LineColor(kRed), RooFit::Name("standard"));

  a.setVal(3.);
  std::unique_ptr<RooFitResult> fit2( pdf.fitTo(data, RooFit::Save(), RooFit::PrintLevel(-1),
      RooFit::Optimize(0),
      RooFit::BatchMode(true),
      RooFit::IntegrateBins(1.E-3)) );
  pdf.plotOn(frame.get(), RooFit::LineColor(kBlue), RooFit::Name("highRes"));


  auto getVal = [](const char* name, const RooArgSet& set) {
    return dynamic_cast<const RooRealVar&>(set[name]).getVal();
  };
  auto getErr = [](const char* name, const RooArgSet& set) {
    return dynamic_cast<const RooRealVar&>(set[name]).getError();
  };

  EXPECT_GT(fabs(getVal("a", targetValues) - getVal("a", fit1->floatParsFinal())), 1. * getErr("a", fit1->floatParsFinal()))
      << "Expecting a bias when sampling PDF in bin centre.";

  EXPECT_NEAR(getVal("a", targetValues), getVal("a", fit2->floatParsFinal()), 1. * getErr("a", fit2->floatParsFinal()))
      << "Expect reduced bias with high-resolution sampling.";

  // Note: We cannot compare with the unbinned dataset here, because when it's plotted, it's filled into a
  // histogram with uniform binning. It therefore creates a jumpy distribution. When comparing with the original
  // data hist, we don't get those jumps.
  EXPECT_GT(frame->chiSquare("standard", "dataHist", 1) * 0.9, frame->chiSquare("highRes",  "dataHist", 1))
      << "Expect chi2/ndf at least 10% better.";

//  fit1->Print();
//  fit2->Print();
}


/// Test the same, but now with RooDataHist. Here, the feature should switch on automatically.
TEST(RooNLLVar, IntegrateBins_RooDataHist) {
  RooRealVar x("x", "x", 0.1, 5.);
  x.setBins(10);

  RooRealVar a("a", "a", -0.3, -5., 5.);
  RooArgSet targetValues;
  RooArgSet(a).snapshot(targetValues);

  RooGenericPdf pdf("pow", "std::pow(x, a)", RooArgSet(x, a));
  std::unique_ptr<RooDataHist> data(pdf.generateBinned(x,  10000));

  std::unique_ptr<RooPlot> frame( x.frame() );
  data->plotOn(frame.get(), RooFit::Name("data"));


  a.setVal(3.);
  std::unique_ptr<RooFitResult> fit1( pdf.fitTo(*data, RooFit::Save(), RooFit::PrintLevel(-1),
      RooFit::BatchMode(true),
      RooFit::IntegrateBins(-1.) // Disable forcefully
      ) );
  pdf.plotOn(frame.get(), RooFit::LineColor(kRed), RooFit::Name("standard"));

  a.setVal(3.);
  std::unique_ptr<RooFitResult> fit2( pdf.fitTo(*data, RooFit::Save(), RooFit::PrintLevel(-1),
      RooFit::BatchMode(true),
      RooFit::IntegrateBins(0.) // Auto-enable for all RooDataHists.
      ) );
  pdf.plotOn(frame.get(), RooFit::LineColor(kBlue), RooFit::Name("highRes"));


  auto getVal = [](const char* name, const RooArgSet& set) {
    return dynamic_cast<const RooRealVar&>(set[name]).getVal();
  };
  auto getErr = [](const char* name, const RooArgSet& set) {
    return dynamic_cast<const RooRealVar&>(set[name]).getError();
  };

  EXPECT_GT(fabs(getVal("a", targetValues) - getVal("a", fit1->floatParsFinal())), 1. * getErr("a", fit1->floatParsFinal()))
      << "Expecting a bias when sampling PDF in bin centre.";

  EXPECT_NEAR(getVal("a", targetValues), getVal("a", fit2->floatParsFinal()), 1. * getErr("a", fit2->floatParsFinal()))
      << "Expect reduced bias with high-resolution sampling.";

  EXPECT_GT(frame->chiSquare("standard", "data", 1) * 0.9, frame->chiSquare("highRes",  "data", 1))
      << "Expect chi2/ndf at least 10% better.";

//  fit1->Print();
//  fit2->Print();
}


TEST(RooChi2Var, IntegrateBins) {
  RooRandom::randomGenerator()->SetSeed(1337ul);

  RooRealVar x("x", "x", 0.1, 5.1);
  x.setBins(10);

  RooRealVar a("a", "a", -0.3, -5., 5.);
  RooArgSet targetValues;
  RooArgSet(a).snapshot(targetValues);

  RooGenericPdf pdf("pow", "std::pow(x, a)", RooArgSet(x, a));
  std::unique_ptr<RooDataHist> dataH(pdf.generateBinned(x,  10000));

  std::unique_ptr<RooPlot> frame( x.frame() );
  dataH->plotOn(frame.get(), RooFit::MarkerColor(kRed));


  a.setVal(3.);
  std::unique_ptr<RooFitResult> fit1( pdf.chi2FitTo(*dataH, RooFit::Save(), RooFit::PrintLevel(-1)) );
  pdf.plotOn(frame.get(), RooFit::LineColor(kRed), RooFit::Name("standard"));

  a.setVal(3.);
  std::unique_ptr<RooFitResult> fit2( pdf.chi2FitTo(*dataH, RooFit::Save(), RooFit::PrintLevel(-1),
      RooFit::BatchMode(true),
      RooFit::IntegrateBins(1.E-3)) );
  pdf.plotOn(frame.get(), RooFit::LineColor(kBlue), RooFit::Name("highRes"));


  auto getVal = [](const char* name, const RooArgSet& set) {
    return dynamic_cast<const RooRealVar&>(set[name]).getVal();
  };
  auto getErr = [](const char* name, const RooArgSet& set) {
    return dynamic_cast<const RooRealVar&>(set[name]).getError();
  };

  EXPECT_GT(fabs(getVal("a", targetValues) - getVal("a", fit1->floatParsFinal())), 1. * getErr("a", fit1->floatParsFinal()))
      << "Expecting a bias when sampling PDF in bin centre.";

  EXPECT_NEAR(getVal("a", targetValues), getVal("a", fit2->floatParsFinal()), 1. * getErr("a", fit2->floatParsFinal()))
      << "Expect reduced bias with high-resolution sampling.";

  EXPECT_GT(frame->chiSquare("standard", nullptr, 1) * 0.9, frame->chiSquare("highRes", nullptr, 1))
      << "Expect chi2/ndf at least 10% better.";

//  fit1->Print();
//  fit2->Print();
}
