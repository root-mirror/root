/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   PB, Patrick Bos,     NL eScience Center, p.bos@esciencecenter.nl        *
 *   IP, Inti Pelupessy,  NL eScience Center, i.pelupessy@esciencecenter.nl  *
 *   VC, Vince Croft,     DIANA / NYU,        vincent.croft@cern.ch          *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#include <MultiProcess/NLLVar.h>

#include <RooRandom.h>
#include <RooWorkspace.h>
#include <RooAbsPdf.h>
#include <RooDataSet.h>
#include <RooMinimizer.h>
#include <RooFitResult.h>
#include <RooProdPdf.h>
#include <RooChebychev.h>
#include <RooStats/RooStatsUtils.h>
#include <RooAddition.h>
#include <RooConstraintSum.h>
#include <RooPolynomial.h>
#include <RooDataHist.h>
#include <RooRealSumPdf.h>

#include "gtest/gtest.h"
#include "../test_lib.h" // Hex


class MultiProcessVectorNLL : public ::testing::TestWithParam<std::tuple<std::size_t, RooFit::MultiProcess::NLLVarTask, std::size_t>> {};


TEST_P(MultiProcessVectorNLL, getVal) {
  // Real-life test: calculate a NLL using event-based parallelization. This
  // should replicate RooRealMPFE results.
  RooRandom::randomGenerator()->SetSeed(std::get<2>(GetParam()));
  RooWorkspace w;
  w.factory("Gaussian::g(x[-5,5],mu[0,-3,3],sigma[1])");
  auto x = w.var("x");
  RooAbsPdf *pdf = w.pdf("g");
  RooDataSet *data = pdf->generate(RooArgSet(*x), 10000);
  auto nll = pdf->createNLL(*data);

  auto nominal_result = nll->getVal();

  std::size_t NumCPU = std::get<0>(GetParam());
  RooFit::MultiProcess::NLLVarTask mp_task_mode = std::get<1>(GetParam());

  RooFit::MultiProcess::NLLVar nll_mp(NumCPU, mp_task_mode, *dynamic_cast<RooNLLVar*>(nll));

  auto mp_result = nll_mp.getVal();

  EXPECT_DOUBLE_EQ(Hex(nominal_result), Hex(mp_result));
  if (HasFailure()) {
    std::cout << "failed test had parameters NumCPU = " << NumCPU << ", task_mode = " << mp_task_mode << ", seed = " << std::get<2>(GetParam()) << std::endl;
  }
}

void check_NLL_type(RooAbsReal *nll) {
  if (dynamic_cast<RooAddition*>(nll) != nullptr) {
    std::cout << "the NLL object is a RooAddition*..." << std::endl;
    bool has_rooconstraintsum = false;
    RooFIter nll_component_iter = nll->getComponents()->fwdIterator();
    RooAbsArg *nll_component;
    while ((nll_component = nll_component_iter.next())) {
      if (nll_component->IsA() == RooConstraintSum::Class()) {
        has_rooconstraintsum = true;
        break;
      } else if (nll_component->IsA() != RooNLLVar::Class() && nll_component->IsA() != RooAddition::Class()) {
        std::cerr << "... containing an unexpected component class: " << nll_component->ClassName() << std::endl;
        throw std::runtime_error("RooAddition* type NLL object contains unexpected component class!");
      }
    }
    if (has_rooconstraintsum) {
      std::cout << "...containing a RooConstraintSum component: " << nll_component->GetName() << std::endl;
    } else {
      std::cout << "...containing only RooNLLVar components." << std::endl;
    }
  } else if (dynamic_cast<RooNLLVar*>(nll) != nullptr) {
    std::cout << "the NLL object is a RooNLLVar*" << std::endl;
  }
}


void count_NLL_components(RooAbsReal *nll) {
  if (dynamic_cast<RooAddition*>(nll) != nullptr) {
    std::cout << "the NLL object is a RooAddition*..." << std::endl;
    unsigned nll_component_count = 0;
    RooFIter nll_component_iter = nll->getComponents()->fwdIterator();
    RooAbsArg *nll_component;
    while ((nll_component = nll_component_iter.next())) {
      if (nll_component->IsA() != RooNLLVar::Class()) {
        ++nll_component_count;
      }
    }
    std::cout << "...containing " << nll_component_count << " RooNLLVar components." << std::endl;
  } else if (dynamic_cast<RooNLLVar*>(nll) != nullptr) {
    std::cout << "the NLL object is a RooNLLVar*" << std::endl;
  }
}


TEST_P(MultiProcessVectorNLL, getValRooAddition) {
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  std::size_t NumCPU = std::get<0>(GetParam());

  RooRandom::randomGenerator()->SetSeed(std::get<2>(GetParam()));

  RooWorkspace w;
  w.factory("Gaussian::g(x[-10,10],mu[0,-3,3],sigma[1])");

  RooRealVar *x = w.var("x");
  x->setRange("x_range",-3,0);
  x->setRange("another_range",1,7);

  RooAbsPdf *pdf = w.pdf("g");
  RooDataSet *data = pdf->generate(*x, 10000);

  RooAbsReal *nll = pdf->createNLL(*data, RooFit::NumCPU(NumCPU),
                                   RooFit::Range("x_range"), RooFit::Range("another_range"));

  check_NLL_type(nll);
  count_NLL_components(nll);

  delete nll;
  delete data;
}


TEST_P(MultiProcessVectorNLL, getValRooConstraintSumAddition) {
  // modified from https://github.com/roofit-dev/rootbench/blob/43d12f33e8dac7af7d587b53a2804ddf6717e92f/root/roofit/roofit/RooFitASUM.cxx#L417

  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  int cpu = 1;
  int bins = 10000;

  RooRealVar x("x","x",0,bins);
  x.setBins(bins);
// Parameters
  RooRealVar a0("a0","a0",0);
  RooRealVar a1("a1","a1",1,0,2);
  RooRealVar a2("a2","a2",0);

  RooPolynomial p0("p0","p0",x);
  RooPolynomial p1("p1","p1",x,RooArgList(a0,a1,a2),0);

  RooDataHist *dh_bkg = p0.generateBinned(x, 1000000000);
  RooDataHist *dh_sig = p1.generateBinned(x, 100000000);
  dh_bkg->SetName("dh_bkg");
  dh_sig->SetName("dh_sig");

  a1.setVal(2);
  RooDataHist *dh_sig_up = p1.generateBinned(x, 1100000000);
  dh_sig_up->SetName("dh_sig_up");
  a1.setVal(.5);
  RooDataHist *dh_sig_down = p1.generateBinned(x, 900000000);
  dh_sig_down->SetName("dh_sig_down");

  RooWorkspace w = RooWorkspace("w");
  w.import(x);
  w.import(*dh_sig);
  w.import(*dh_bkg);
  w.import(*dh_sig_up);
  w.import(*dh_sig_down);
  w.factory("HistFunc::hf_sig(x,dh_sig)");
  w.factory("HistFunc::hf_bkg(x,dh_bkg)");
  w.factory("HistFunc::hf_sig_up(x,dh_sig_up)");
  w.factory("HistFunc::hf_sig_down(x,dh_sig_down)");
  w.factory("PiecewiseInterpolation::pi_sig(hf_sig,hf_sig_down,hf_sig_up,alpha[-5,5])");

  w.factory("ASUM::model(mu[1,0,5]*pi_sig,nu[1]*hf_bkg)");
  w.factory("Gaussian::constraint(alpha,0,1)");
  w.factory("PROD::model2(model,constraint)");

  RooAbsPdf *pdf = w.pdf("model2");

  RooDataHist *data = pdf->generateBinned(x, 1100000);
  RooAbsReal *nll = pdf->createNLL(*data, RooFit::NumCPU(cpu, 0));

  check_NLL_type(nll);
  count_NLL_components(nll);

  delete nll;
}

TEST_P(MultiProcessVectorNLL, setVal) {
  // calculate the NLL twice with different parameters

  RooRandom::randomGenerator()->SetSeed(std::get<2>(GetParam()));
  RooWorkspace w;
  w.factory("Gaussian::g(x[-5,5],mu[0,-3,3],sigma[1])");
  auto x = w.var("x");
  RooAbsPdf *pdf = w.pdf("g");
  RooDataSet *data = pdf->generate(RooArgSet(*x), 10000);
  auto nll = pdf->createNLL(*data);

  std::size_t NumCPU = std::get<0>(GetParam());
  RooFit::MultiProcess::NLLVarTask mp_task_mode = std::get<1>(GetParam());

  RooFit::MultiProcess::NLLVar nll_mp(NumCPU, mp_task_mode, *dynamic_cast<RooNLLVar*>(nll));

  // calculate first results
  nll->getVal();
  nll_mp.getVal();

  w.var("mu")->setVal(2);

  // calculate second results after parameter change
  auto nominal_result2 = nll->getVal();
  auto mp_result2 = nll_mp.getVal();

  EXPECT_DOUBLE_EQ(Hex(nominal_result2), Hex(mp_result2));
  if (HasFailure()) {
    std::cout << "failed test had parameters NumCPU = " << NumCPU << ", task_mode = " << mp_task_mode << ", seed = " << std::get<2>(GetParam()) << std::endl;
  }
}


INSTANTIATE_TEST_SUITE_P(NworkersModeSeed,
                        MultiProcessVectorNLL,
                        ::testing::Combine(::testing::Values(1,2,3),  // number of workers
                                           ::testing::Values(RooFit::MultiProcess::NLLVarTask::all_events,
                                                             RooFit::MultiProcess::NLLVarTask::single_event,
                                                             RooFit::MultiProcess::NLLVarTask::bulk_partition,
                                                             RooFit::MultiProcess::NLLVarTask::interleave),
                                           ::testing::Values(2,3)));  // random seed



class NLLMultiProcessVsMPFE : public ::testing::TestWithParam<std::tuple<std::size_t, RooFit::MultiProcess::NLLVarTask, std::size_t>> {};

TEST_P(NLLMultiProcessVsMPFE, getVal) {
  // Compare our MP NLL to actual RooRealMPFE results using the same strategies.

  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  // parameters
  std::size_t NumCPU = std::get<0>(GetParam());
  RooFit::MultiProcess::NLLVarTask mp_task_mode = std::get<1>(GetParam());
  std::size_t seed = std::get<2>(GetParam());

  RooRandom::randomGenerator()->SetSeed(seed);

  RooWorkspace w;
  w.factory("Gaussian::g(x[-5,5],mu[0,-3,3],sigma[1])");
  auto x = w.var("x");
  RooAbsPdf *pdf = w.pdf("g");
  RooDataSet *data = pdf->generate(RooArgSet(*x), 10000);

  int mpfe_task_mode = 0;
  if (mp_task_mode == RooFit::MultiProcess::NLLVarTask::interleave) {
    mpfe_task_mode = 1;
  }

  auto nll_mpfe = pdf->createNLL(*data, RooFit::NumCPU(NumCPU, mpfe_task_mode));

  auto mpfe_result = nll_mpfe->getVal();

  // create new nll without MPFE for creating nll_mp (an MPFE-enabled RooNLLVar interferes with MP::Vector's bipe use)
  auto nll = pdf->createNLL(*data);
  RooFit::MultiProcess::NLLVar nll_mp(NumCPU, mp_task_mode, *dynamic_cast<RooNLLVar*>(nll));

  auto mp_result = nll_mp.getVal();

  EXPECT_EQ(Hex(mpfe_result), Hex(mp_result));
  if (HasFailure()) {
    std::cout << "failed test had parameters NumCPU = " << NumCPU << ", task_mode = " << mp_task_mode << ", seed = " << seed << std::endl;
  }
}


TEST_P(NLLMultiProcessVsMPFE, minimize) {
  // do a minimization (e.g. like in GradMinimizer_Gaussian1D test)

  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  // TODO: see whether it performs adequately

  // parameters
  std::size_t NumCPU = std::get<0>(GetParam());
  RooFit::MultiProcess::NLLVarTask mp_task_mode = std::get<1>(GetParam());
  std::size_t seed = std::get<2>(GetParam());

  RooRandom::randomGenerator()->SetSeed(seed);

  RooWorkspace w = RooWorkspace();

  w.factory("Gaussian::g(x[-5,5],mu[0,-3,3],sigma[1])");
  auto x = w.var("x");
  RooAbsPdf *pdf = w.pdf("g");
  RooRealVar *mu = w.var("mu");

  RooDataSet *data = pdf->generate(RooArgSet(*x), 10000);
  mu->setVal(-2.9);

  int mpfe_task_mode;
  switch (mp_task_mode) {
    case RooFit::MultiProcess::NLLVarTask::bulk_partition: {
      mpfe_task_mode = 0;
      break;
    }
    case RooFit::MultiProcess::NLLVarTask::interleave: {
      mpfe_task_mode = 1;
      break;
    }
    default: {
      throw std::logic_error("can only compare bulk_partition and interleave strategies to MPFE NLL");
    }
  }

  auto nll_mpfe = pdf->createNLL(*data, RooFit::NumCPU(NumCPU, mpfe_task_mode));
  auto nll_nominal = pdf->createNLL(*data);
  RooFit::MultiProcess::NLLVar nll_mp(NumCPU, mp_task_mode, *dynamic_cast<RooNLLVar*>(nll_nominal));

  // save initial values for the start of all minimizations
  RooArgSet values = RooArgSet(*mu, *pdf);

  RooArgSet *savedValues = dynamic_cast<RooArgSet *>(values.snapshot());
  if (savedValues == nullptr) {
    throw std::runtime_error("params->snapshot() cannot be casted to RooArgSet!");
  }

  // --------

  RooMinimizer m0(*nll_mpfe);
  m0.setMinimizerType("Minuit2");

  m0.setStrategy(0);
  m0.setPrintLevel(-1);

  m0.migrad();

  RooFitResult *m0result = m0.lastMinuitFit();
  double minNll0 = m0result->minNll();
  double edm0 = m0result->edm();
  double mu0 = mu->getVal();
  double muerr0 = mu->getError();

  values = *savedValues;

  RooMinimizer m1(nll_mp);
  m1.setMinimizerType("Minuit2");

  m1.setStrategy(0);
  m1.setPrintLevel(-1);

  m1.migrad();

  RooFitResult *m1result = m1.lastMinuitFit();
  double minNll1 = m1result->minNll();
  double edm1 = m1result->edm();
  double mu1 = mu->getVal();
  double muerr1 = mu->getError();

  EXPECT_EQ(minNll0, minNll1);
  EXPECT_EQ(mu0, mu1);
  EXPECT_EQ(muerr0, muerr1);
  EXPECT_EQ(edm0, edm1);

  m1.cleanup();  // necessary in tests to clean up global _theFitter
}


INSTANTIATE_TEST_SUITE_P(NworkersModeSeed,
                        NLLMultiProcessVsMPFE,
                        ::testing::Combine(::testing::Values(2,3),  // number of workers
                                           ::testing::Values(RooFit::MultiProcess::NLLVarTask::bulk_partition,
                                                             RooFit::MultiProcess::NLLVarTask::interleave),
                                           ::testing::Values(2,3)));  // random seed


TEST(NLLMultiProcessVsMPFE, throwOnCreatingMPwithMPFE) {
  // Using an MPFE-enabled NLL should throw when creating an MP NLL.
  RooWorkspace w;
  w.factory("Gaussian::g(x[-5,5],mu[0,-3,3],sigma[1])");
  auto x = w.var("x");
  RooAbsPdf *pdf = w.pdf("g");
  RooDataSet *data = pdf->generate(RooArgSet(*x), 10);

  RooAbsReal* nll_mpfe = pdf->createNLL(*data, RooFit::NumCPU(2));

  EXPECT_THROW({
                 RooFit::MultiProcess::NLLVar nll_mp(2, RooFit::MultiProcess::NLLVarTask::bulk_partition, *dynamic_cast<RooNLLVar*>(nll_mpfe));
               }, std::logic_error);

  delete nll_mpfe;
}
