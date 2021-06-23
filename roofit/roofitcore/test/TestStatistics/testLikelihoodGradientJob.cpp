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

#include <stdexcept> // runtime_error

#include <TFile.h>

#include <RooRandom.h>
#include <RooWorkspace.h>
#include <RooTimer.h>

#include "RooDataHist.h" // complete type in Binned test
#include "RooCategory.h" // complete type in MultiBinnedConstraint test

#include <RooMinimizer.h>
#include <RooGradMinimizerFcn.h>
#include <RooFitResult.h>

#include <RooStats/ModelConfig.h>

#include <TestStatistics/LikelihoodGradientJob.h>
#include <TestStatistics/LikelihoodSerial.h>
#include <TestStatistics/RooUnbinnedL.h>
#include <TestStatistics/likelihood_builders.h>
#include <RooFit/MultiProcess/JobManager.h>
#include <RooFit/MultiProcess/ProcessManager.h> // need to complete type for debugging

#include "gtest/gtest.h"
#include "../test_lib.h" // generate_1D_gaussian_pdf_nll


class Environment : public testing::Environment {
public:
   void SetUp() override {
//      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
   }
};

testing::Environment* const test_env =
   testing::AddGlobalTestEnvironment(new Environment);


class LikelihoodGradientJob : public ::testing::TestWithParam<std::tuple<std::size_t, std::size_t>> {
};

TEST_P(LikelihoodGradientJob, Gaussian1D)
{
   // do a minimization, but now using GradMinimizer and its MP version

   // parameters
   std::size_t NWorkers = std::get<0>(GetParam());
   //  RooFit::MultiProcess::NLLVarTask mp_task_mode = std::get<1>(GetParam());
   std::size_t seed = std::get<1>(GetParam());

   RooRandom::randomGenerator()->SetSeed(seed);

   RooWorkspace w = RooWorkspace();

   std::unique_ptr<RooAbsReal> nll;
   std::unique_ptr<RooArgSet> values;
   RooAbsPdf *pdf;
   RooDataSet *data;
   std::tie(nll, pdf, data, values) = generate_1D_gaussian_pdf_nll(w, 10000);
   // when c++17 support arrives, change to this:
   //  auto [nll, pdf, data, values] = generate_1D_gaussian_pdf_nll(w, 10000);
   RooRealVar *mu = w.var("mu");

   RooArgSet *savedValues = dynamic_cast<RooArgSet *>(values->snapshot());
   if (savedValues == nullptr) {
      throw std::runtime_error("params->snapshot() cannot be casted to RooArgSet!");
   }

   // --------

   std::unique_ptr<RooMinimizer> m0 = RooMinimizer::create<RooGradMinimizerFcn>(*nll);
   m0->setMinimizerType("Minuit2");

   m0->setStrategy(0);
   m0->setVerbose(true);
   m0->setPrintLevel(1);

   m0->migrad();

   RooFitResult *m0result = m0->lastMinuitFit();
   double minNll0 = m0result->minNll();
   double edm0 = m0result->edm();
   double mu0 = mu->getVal();
   double muerr0 = mu->getError();

   *values = *savedValues;

   RooFit::MultiProcess::JobManager::default_N_workers = NWorkers;
   auto likelihood = std::make_shared<RooFit::TestStatistics::RooUnbinnedL>(pdf, data);
   std::unique_ptr<RooMinimizer> m1 =
      RooMinimizer::create<RooFit::TestStatistics::LikelihoodSerial, RooFit::TestStatistics::LikelihoodGradientJob>(
         likelihood);
   m1->setMinimizerType("Minuit2");

   m1->setStrategy(0);
   m1->setVerbose(true);
   m1->setPrintLevel(1);

   m1->migrad();

   RooFitResult *m1result = m1->lastMinuitFit();
   double minNll1 = m1result->minNll();
   double edm1 = m1result->edm();
   double mu1 = mu->getVal();
   double muerr1 = mu->getError();

   EXPECT_EQ(minNll0, minNll1);
   EXPECT_EQ(mu0, mu1);
   EXPECT_EQ(muerr0, muerr1);
   EXPECT_EQ(edm0, edm1);

   m1->cleanup(); // necessary in tests to clean up global _theFitter
}

TEST(LikelihoodGradientJobDEBUGGING, DISABLED_Gaussian1DNominal)
{
   //  std::size_t NWorkers = 1;
   std::size_t seed = 1;

   RooRandom::randomGenerator()->SetSeed(seed);

   RooWorkspace w = RooWorkspace();

   std::unique_ptr<RooAbsReal> nll;
   std::unique_ptr<RooArgSet> _;
   RooAbsPdf *pdf;
   RooDataSet *data;
   std::tie(nll, pdf, data, _) = generate_1D_gaussian_pdf_nll(w, 10000);

   std::unique_ptr<RooMinimizer> m0 = RooMinimizer::create<RooGradMinimizerFcn>(*nll);
   m0->setMinimizerType("Minuit2");

   m0->setStrategy(0);
   m0->setPrintLevel(2);

   m0->migrad();
   m0->cleanup(); // necessary in tests to clean up global _theFitter
}

TEST(LikelihoodGradientJobDEBUGGING, DISABLED_Gaussian1DMultiProcess)
{
   std::size_t NWorkers = 1;
   std::size_t seed = 1;

   RooRandom::randomGenerator()->SetSeed(seed);

   RooWorkspace w = RooWorkspace();

   std::unique_ptr<RooAbsReal> nll;
   std::unique_ptr<RooArgSet> values;
   RooAbsPdf *pdf;
   RooDataSet *data;
   std::tie(nll, pdf, data, values) = generate_1D_gaussian_pdf_nll(w, 10000);

   RooFit::MultiProcess::JobManager::default_N_workers = NWorkers;
   auto likelihood = std::make_shared<RooFit::TestStatistics::RooUnbinnedL>(pdf, data);
   std::unique_ptr<RooMinimizer> m1 =
      RooMinimizer::create<RooFit::TestStatistics::LikelihoodSerial, RooFit::TestStatistics::LikelihoodGradientJob>(
         likelihood);
   m1->setMinimizerType("Minuit2");

   m1->setStrategy(0);
   m1->setPrintLevel(2);

   m1->migrad();
   m1->cleanup(); // necessary in tests to clean up global _theFitter
}

TEST(LikelihoodGradientJob, RepeatMigrad)
{
   // do multiple minimizations using MP::GradMinimizer, testing breakdown and rebuild

   // parameters
   std::size_t NWorkers = 2;
   std::size_t seed = 5;

   RooRandom::randomGenerator()->SetSeed(seed);

   RooWorkspace w = RooWorkspace();

   std::unique_ptr<RooAbsReal> nll;
   std::unique_ptr<RooArgSet> values;
   RooAbsPdf *pdf;
   RooDataSet *data;
   std::tie(nll, pdf, data, values) = generate_1D_gaussian_pdf_nll(w, 10000);
   // when c++17 support arrives, change to this:
   //  auto [nll, values] = generate_1D_gaussian_pdf_nll(w, 10000);

   RooArgSet *savedValues = dynamic_cast<RooArgSet *>(values->snapshot());
   if (savedValues == nullptr) {
      throw std::runtime_error("params->snapshot() cannot be casted to RooArgSet!");
   }

   // --------

   RooFit::MultiProcess::JobManager::default_N_workers = NWorkers;
   auto likelihood = std::make_shared<RooFit::TestStatistics::RooUnbinnedL>(pdf, data);
   std::unique_ptr<RooMinimizer> m1 =
      RooMinimizer::create<RooFit::TestStatistics::LikelihoodSerial, RooFit::TestStatistics::LikelihoodGradientJob>(
         likelihood);

   m1->setMinimizerType("Minuit2");

   m1->setStrategy(0);
   m1->setPrintLevel(-1);

   std::cout << "... running migrad first time ..." << std::endl;
   m1->migrad();

//   std::cout << "... terminating JobManager instance ..." << std::endl;
//   RooFit::MultiProcess::JobManager::instance()->terminate();

   *values = *savedValues;

   std::cout << "... running migrad second time ..." << std::endl;
   m1->migrad();

   std::cout << "... cleaning up minimizer ..." << std::endl;
   m1->cleanup(); // necessary in tests to clean up global _theFitter
}

TEST_P(LikelihoodGradientJob, GaussianND)
{
   // do a minimization, but now using GradMinimizer and its MP version

   // parameters
   std::size_t NWorkers = std::get<0>(GetParam());
   //  RooFit::MultiProcess::NLLVarTask mp_task_mode = std::get<1>(GetParam());
   std::size_t seed = std::get<1>(GetParam());

   unsigned int N = 4;

   RooRandom::randomGenerator()->SetSeed(seed);

   RooWorkspace w = RooWorkspace();

   std::unique_ptr<RooAbsReal> nll;
   std::unique_ptr<RooArgSet> values;
   RooAbsPdf *pdf;
   RooDataSet *data;
   std::tie(nll, pdf, data, values) = generate_ND_gaussian_pdf_nll(w, N, 1000);
   // when c++17 support arrives, change to this:
   //  auto [nll, all_values] = generate_ND_gaussian_pdf_nll(w, N, 1000);

   RooArgSet *savedValues = dynamic_cast<RooArgSet *>(values->snapshot());
   if (savedValues == nullptr) {
      throw std::runtime_error("params->snapshot() cannot be casted to RooArgSet!");
   }

   RooWallTimer wtimer;

   // --------

   std::unique_ptr<RooMinimizer> m0 = RooMinimizer::create<RooGradMinimizerFcn>(*nll);
   m0->setMinimizerType("Minuit2");

   m0->setStrategy(0);
   m0->setPrintLevel(-1);

   wtimer.start();
   m0->migrad();
   wtimer.stop();
   std::cout << "\nwall clock time RooGradMinimizer.migrad (NWorkers = " << NWorkers << ", seed = " << seed
             << "): " << wtimer.timing_s() << " s" << std::endl;

   RooFitResult *m0result = m0->lastMinuitFit();
   double minNll0 = m0result->minNll();
   double edm0 = m0result->edm();
   double mean0[N];
   double std0[N];
   for (unsigned ix = 0; ix < N; ++ix) {
      {
         std::ostringstream os;
         os << "m" << ix;
         mean0[ix] = dynamic_cast<RooRealVar *>(w.arg(os.str().c_str()))->getVal();
      }
      {
         std::ostringstream os;
         os << "s" << ix;
         std0[ix] = dynamic_cast<RooRealVar *>(w.arg(os.str().c_str()))->getVal();
      }
   }

   // --------

   *values = *savedValues;

   // --------

   RooFit::MultiProcess::JobManager::default_N_workers = NWorkers;
   auto likelihood = std::make_shared<RooFit::TestStatistics::RooUnbinnedL>(pdf, data);
   std::unique_ptr<RooMinimizer> m1 =
      RooMinimizer::create<RooFit::TestStatistics::LikelihoodSerial, RooFit::TestStatistics::LikelihoodGradientJob>(
         likelihood);
   m1->setMinimizerType("Minuit2");

   m1->setStrategy(0);
   m1->setPrintLevel(-1);

   wtimer.start();
   m1->migrad();
   wtimer.stop();
   std::cout << "wall clock time MP::GradMinimizer.migrad (NWorkers = " << NWorkers << ", seed = " << seed
             << "): " << wtimer.timing_s() << " s\n"
             << std::endl;

   RooFitResult *m1result = m1->lastMinuitFit();
   double minNll1 = m1result->minNll();
   double edm1 = m1result->edm();
   double mean1[N];
   double std1[N];
   for (unsigned ix = 0; ix < N; ++ix) {
      {
         std::ostringstream os;
         os << "m" << ix;
         mean1[ix] = dynamic_cast<RooRealVar *>(w.arg(os.str().c_str()))->getVal();
      }
      {
         std::ostringstream os;
         os << "s" << ix;
         std1[ix] = dynamic_cast<RooRealVar *>(w.arg(os.str().c_str()))->getVal();
      }
   }

   EXPECT_EQ(minNll0, minNll1);
   EXPECT_EQ(edm0, edm1);

   for (unsigned ix = 0; ix < N; ++ix) {
      EXPECT_EQ(mean0[ix], mean1[ix]);
      EXPECT_EQ(std0[ix], std1[ix]);
   }

   m1->cleanup(); // necessary in tests to clean up global _theFitter
}

INSTANTIATE_TEST_SUITE_P(NworkersSeed, LikelihoodGradientJob,
                        ::testing::Combine(::testing::Values(1, 2, 3), // number of workers
                                           ::testing::Values(2, 3)));  // random seed


class BasicTest: public ::testing::Test {
protected:
   void SetUp() override {
      RooRandom::randomGenerator()->SetSeed(seed);
      clean_flags = std::make_shared<RooFit::TestStatistics::WrapperCalculationCleanFlags>();
   }

   std::size_t seed = 23;
   RooWorkspace w;
   std::unique_ptr<RooAbsReal> nll;
   RooAbsPdf *pdf;
   RooAbsData *data;
   std::shared_ptr<RooFit::TestStatistics::RooAbsL> likelihood;
   std::shared_ptr<RooFit::TestStatistics::WrapperCalculationCleanFlags> clean_flags;
};


class LikelihoodSimBinnedConstrainedTest : public BasicTest {
protected:
   void SetUp() override {
      BasicTest::SetUp();
      // Unbinned pdfs that define template histograms

      w.factory("Gaussian::gA(x[-10,10],-2,3)") ;
      w.factory("Gaussian::gB(x[-10,10],2,1)") ;
      w.factory("Uniform::u(x)");

      // Generate template histograms

      RooDataHist* h_sigA = w.pdf("gA")->generateBinned(*w.var("x"),1000) ;
      RooDataHist* h_sigB = w.pdf("gB")->generateBinned(*w.var("x"),1000) ;
      RooDataHist *h_bkg = w.pdf("u")->generateBinned(*w.var("x"), 1000);

      w.import(*h_sigA, RooFit::Rename("h_sigA"));
      w.import(*h_sigB, RooFit::Rename("h_sigB"));
      w.import(*h_bkg, RooFit::Rename("h_bkg"));

      // Construct binned pdf as sum of amplitudes
      w.factory("HistFunc::hf_sigA(x,h_sigA)") ;
      w.factory("HistFunc::hf_sigB(x,h_sigB)") ;
      w.factory("HistFunc::hf_bkg(x,h_bkg)") ;

      w.factory("ASUM::model_phys_A(mu_sig[1,-1,10]*hf_sigA,expr::mu_bkg_A('1+0.02*alpha_bkg_A',alpha_bkg_A[-5,5])*hf_bkg)") ;
      w.factory("ASUM::model_phys_B(mu_sig*hf_sigB,expr::mu_bkg_B('1+0.05*alpha_bkg_B',alpha_bkg_B[-5,5])*hf_bkg)") ;

      // Construct L_subs: Gaussian subsidiary measurement that constrains alpha_bkg
      w.factory("Gaussian:model_subs_A(alpha_bkg_obs_A[0],alpha_bkg_A,1)") ;
      w.factory("Gaussian:model_subs_B(alpha_bkg_obs_B[0],alpha_bkg_B,1)") ;

      // Construct full pdfs for each component (A,B)
      w.factory("PROD::model_A(model_phys_A,model_subs_A)") ;
      w.factory("PROD::model_B(model_phys_B,model_subs_B)") ;

      // Construct simulatenous pdf
      w.factory("SIMUL::model(index[A,B],A=model_A,B=model_B)") ;

      pdf = w.pdf("model");
      // Construct dataset from physics pdf
      data = pdf->generate(RooArgSet(*w.var("x"), *w.cat("index")), RooFit::AllBinned());
   }
};

TEST_F(LikelihoodSimBinnedConstrainedTest, BasicParameters)
{
   // original test:
   nll.reset(pdf->createNLL(*data, RooFit::GlobalObservables(RooArgSet(*w.var("alpha_bkg_obs_A"), *w.var("alpha_bkg_obs_B")))));

   // --------

   auto nll0 = nll->getVal();

   likelihood = RooFit::TestStatistics::build_simultaneous_likelihood(
      pdf, data, RooFit::TestStatistics::GlobalObservables({*w.var("alpha_bkg_obs_A"), *w.var("alpha_bkg_obs_B")}));
   RooFit::TestStatistics::LikelihoodSerial nll_ts(likelihood, clean_flags/*, nullptr*/);

   nll_ts.evaluate();
   auto nll1 = nll_ts.return_result();

   EXPECT_DOUBLE_EQ(nll0, nll1);
}

TEST_F(LikelihoodSimBinnedConstrainedTest, Minimize)
{
   // do a minimization, but now using GradMinimizer and its MP version
   nll.reset(pdf->createNLL(*data, RooFit::Constrain(RooArgSet(*w.var("alpha_bkg_obs_A"))),
                            RooFit::GlobalObservables(RooArgSet(*w.var("alpha_bkg_obs_B"))), RooFit::Offset(kTRUE)));

   // parameters
   std::size_t NWorkers = 2; //std::get<0>(GetParam());

   RooArgSet *values = pdf->getParameters(data);

   values->add(*pdf);
   values->add(*nll);

   RooArgSet *savedValues = dynamic_cast<RooArgSet *>(values->snapshot());
   if (savedValues == nullptr) {
      throw std::runtime_error("params->snapshot() cannot be casted to RooArgSet!");
   }

   // --------

   std::unique_ptr<RooMinimizer> m0 = RooMinimizer::create<RooGradMinimizerFcn>(*nll);

   m0->setMinimizerType("Minuit2");
   m0->setStrategy(0);
//   m0->setVerbose(true);
   m0->setPrintLevel(1);

   m0->migrad();

   RooFitResult *m0result = m0->lastMinuitFit();
   double minNll_nominal = m0result->minNll();
   double edm_nominal = m0result->edm();
   double alpha_bkg_A_nominal = w.var("alpha_bkg_A")->getVal();
   double alpha_bkg_A_error_nominal = w.var("alpha_bkg_A")->getError();
   double alpha_bkg_B_nominal = w.var("alpha_bkg_B")->getVal();
   double alpha_bkg_B_error_nominal = w.var("alpha_bkg_B")->getError();
   double mu_sig_nominal = w.var("mu_sig")->getVal();
   double mu_sig_error_nominal = w.var("mu_sig")->getError();

   *values = *savedValues;

   RooFit::MultiProcess::JobManager::default_N_workers = NWorkers;

   auto likelihood = RooFit::TestStatistics::build_simultaneous_likelihood(
      pdf, data, RooFit::TestStatistics::ConstrainedParameters({*w.var("alpha_bkg_obs_A")}),
      RooFit::TestStatistics::GlobalObservables({*w.var("alpha_bkg_obs_B")}));

   std::unique_ptr<RooMinimizer> m1 = RooMinimizer::create<RooFit::TestStatistics::LikelihoodSerial, RooFit::TestStatistics::LikelihoodGradientJob>(likelihood);
   m1->enable_likelihood_offsetting(true);

   m1->setMinimizerType("Minuit2");
   m1->setStrategy(0);
//   m1->setVerbose(true);
   m1->setPrintLevel(1);
   m1->optimizeConst(2);

   m1->migrad();

   RooFitResult *m1result = m1->lastMinuitFit();
   double minNll_GradientJob = m1result->minNll();
   double edm_GradientJob = m1result->edm();
   double alpha_bkg_A_GradientJob = w.var("alpha_bkg_A")->getVal();
   double alpha_bkg_A_error_GradientJob = w.var("alpha_bkg_A")->getError();
   double alpha_bkg_B_GradientJob = w.var("alpha_bkg_B")->getVal();
   double alpha_bkg_B_error_GradientJob = w.var("alpha_bkg_B")->getError();
   double mu_sig_GradientJob = w.var("mu_sig")->getVal();
   double mu_sig_error_GradientJob = w.var("mu_sig")->getError();

   EXPECT_EQ(minNll_nominal, minNll_GradientJob);
   EXPECT_EQ(edm_nominal, edm_GradientJob);
   EXPECT_EQ(alpha_bkg_A_nominal, alpha_bkg_A_GradientJob);
   EXPECT_EQ(alpha_bkg_A_error_nominal, alpha_bkg_A_error_GradientJob);
   EXPECT_EQ(alpha_bkg_B_nominal, alpha_bkg_B_GradientJob);
   EXPECT_EQ(alpha_bkg_B_error_nominal, alpha_bkg_B_error_GradientJob);
   EXPECT_EQ(mu_sig_nominal, mu_sig_GradientJob);
   EXPECT_EQ(mu_sig_error_nominal, mu_sig_error_GradientJob);

   m1->cleanup(); // necessary in tests to clean up global _theFitter
}

class CarstenGGFWorkspaceTest: public ::testing::Test {
protected:
   void SetUp() override {
      RooRandom::randomGenerator()->SetSeed(seed);

      TFile *_file0 = TFile::Open("/Users/pbos/projects/apcocsm/carsten/lxplus/ggF/ggF-stxs1-v1.root");

      w = static_cast<RooWorkspace*>(gDirectory->Get("HWWRun2GGF"));

      data = w->data("obsData");
      auto mc = dynamic_cast<RooStats::ModelConfig *>(w->genobj("ModelConfig"));
      global_observables = mc->GetGlobalObservables();
      nuisance_parameters = mc->GetNuisanceParameters();
      pdf = w->pdf(mc->GetPdf()->GetName());
   }

   std::size_t seed = 23;
   RooWorkspace* w;
   RooAbsPdf *pdf;
   RooAbsData *data;
   const RooArgSet *global_observables;
   const RooArgSet *nuisance_parameters;
   std::unique_ptr<RooMinimizer> m;
};

TEST_F(CarstenGGFWorkspaceTest, DISABLED_NoMultiProcess)
{
   RooAbsReal *nll = pdf->createNLL(*data,
                                    RooFit::GlobalObservables(*global_observables),
                                    RooFit::Constrain(*nuisance_parameters),
                                    RooFit::Offset(kTRUE));

   m = RooMinimizer::create(*nll);

   m->setPrintLevel(1);
   m->setStrategy(0);
   m->setProfile(false);
   m->optimizeConst(2);
   m->setMinimizerType("Minuit2");
//    m->setVerbose(kTRUE);
   m->setEps(1);

   m->migrad();

   m->cleanup(); // necessary in tests to clean up global _theFitter
}

TEST_F(CarstenGGFWorkspaceTest, DISABLED_MultiProcess)
{
   RooFit::MultiProcess::JobManager::default_N_workers = 4;
   auto likelihood = RooFit::TestStatistics::build_simultaneous_likelihood(pdf, data, RooFit::TestStatistics::ConstrainedParameters(*nuisance_parameters), RooFit::TestStatistics::GlobalObservables(*global_observables));
   m = RooMinimizer::create<RooFit::TestStatistics::LikelihoodSerial, RooFit::TestStatistics::LikelihoodGradientJob>(likelihood);
   m->enable_likelihood_offsetting(true);

   m->setPrintLevel(1);
   m->setStrategy(0);
   m->setProfile(false);
   m->optimizeConst(2);
   m->setMinimizerType("Minuit2");
//    m->setVerbose(kTRUE);
   m->setEps(1);

   m->migrad();

   m->cleanup(); // necessary in tests to clean up global _theFitter
}