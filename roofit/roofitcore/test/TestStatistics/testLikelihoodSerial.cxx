/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *
 * Copyright (c) 2016-2020, Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */

#include <stdexcept> // runtime_error

#include <RooRandom.h>
#include <RooWorkspace.h>
#include <RooTimer.h>

#include <RooMinimizer.h>
#include <RooGradMinimizerFcn.h>
#include <RooFitResult.h>

#include "RooDataHist.h" // complete type in Binned test
#include "RooCategory.h" // complete type in MultiBinnedConstraint test

#include <TestStatistics/LikelihoodGradientJob.h>
#include <TestStatistics/LikelihoodSerial.h>
#include <TestStatistics/RooUnbinnedL.h>
#include <TestStatistics/RooBinnedL.h>
#include <TestStatistics/RooSimultaneousL.h>
#include <RooFit/MultiProcess/JobManager.h>
#include <RooFit/MultiProcess/ProcessManager.h> // need to complete type for debugging
#include <RooNLLVar.h>

#include "gtest/gtest.h"
#include "../test_lib.h" // generate_1D_gaussian_pdf_nll

class Environment : public testing::Environment {
public:
   void SetUp() override {
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
   }
};

testing::Environment* const test_env =
   testing::AddGlobalTestEnvironment(new Environment);

class LikelihoodSerialTest: public ::testing::Test {
protected:
   void SetUp() override {
      RooRandom::randomGenerator()->SetSeed(seed);
      clean_flags = std::make_shared<RooFit::TestStatistics::WrapperCalculationCleanFlags>();
   }

   std::size_t seed = 23;
   RooWorkspace w;
   std::unique_ptr<RooAbsReal> nll;
   std::unique_ptr<RooArgSet> values;
   RooAbsPdf *pdf;
   RooAbsData *data;
   std::shared_ptr<RooFit::TestStatistics::RooAbsL> likelihood;
   std::shared_ptr<RooFit::TestStatistics::WrapperCalculationCleanFlags> clean_flags;
};

class LikelihoodSerialBinnedDatasetTest : public LikelihoodSerialTest {
protected:
   void SetUp() override {
      LikelihoodSerialTest::SetUp();

      // Unbinned pdfs that define template histograms
      w.factory("Gaussian::g(x[-10,10],0,2)");
      w.factory("Uniform::u(x)");

      // Generate template histograms
      RooDataHist *h_sig = w.pdf("g")->generateBinned(*w.var("x"), 1000);
      RooDataHist *h_bkg = w.pdf("u")->generateBinned(*w.var("x"), 1000);

      w.import(*h_sig, RooFit::Rename("h_sig"));
      w.import(*h_bkg, RooFit::Rename("h_bkg"));

      // Construct binned pdf as sum of amplitudes
      w.factory("HistFunc::hf_sig(x,h_sig)");
      w.factory("HistFunc::hf_bkg(x,h_bkg)");
      w.factory("ASUM::model(mu_sig[1,-1,10]*hf_sig,mu_bkg[1,-1,10]*hf_bkg)");

      pdf = w.pdf("model");
   }
};

TEST_F(LikelihoodSerialTest, UnbinnedGaussian1D)
{
   std::tie(nll, pdf, data, values) = generate_1D_gaussian_pdf_nll(w, 10000);
   likelihood = std::make_shared<RooFit::TestStatistics::RooUnbinnedL>(pdf, data);
   RooFit::TestStatistics::LikelihoodSerial nll_ts(likelihood, clean_flags, nullptr);

   auto nll0 = nll->getVal();

   nll_ts.evaluate();
   auto nll1 = nll_ts.return_result();

   EXPECT_EQ(nll0, nll1);
}

TEST_F(LikelihoodSerialTest, UnbinnedGaussianND)
{
   unsigned int N = 1;

   std::tie(nll, pdf, data, values) = generate_ND_gaussian_pdf_nll(w, N, 1000);
   likelihood = std::make_shared<RooFit::TestStatistics::RooUnbinnedL>(pdf, data);
   RooFit::TestStatistics::LikelihoodSerial nll_ts(likelihood, clean_flags, nullptr);

   auto nll0 = nll->getVal();

   nll_ts.evaluate();
   auto nll1 = nll_ts.return_result();

   EXPECT_EQ(nll0, nll1);
}

TEST_F(LikelihoodSerialBinnedDatasetTest, UnbinnedPdf)
{
   data = pdf->generateBinned(*w.var("x"));

   nll.reset(pdf->createNLL(*data));

   likelihood = std::make_shared<RooFit::TestStatistics::RooUnbinnedL>(pdf, data);
   RooFit::TestStatistics::LikelihoodSerial nll_ts(likelihood, clean_flags, nullptr);

   auto nll0 = nll->getVal();

   nll_ts.evaluate();
   auto nll1 = nll_ts.return_result();

   EXPECT_EQ(nll0, nll1);
}


TEST_F(LikelihoodSerialBinnedDatasetTest, UnbinnedPdfWithBinnedLikelihoodAttribute)
{
   pdf->setAttribute("BinnedLikelihood");
   data = pdf->generateBinned(*w.var("x"));

   nll.reset(pdf->createNLL(*data));

   likelihood = std::make_shared<RooFit::TestStatistics::RooUnbinnedL>(pdf, data);
   RooFit::TestStatistics::LikelihoodSerial nll_ts(likelihood, clean_flags, nullptr);

   auto nll0 = nll->getVal();

   nll_ts.evaluate();
   auto nll1 = nll_ts.return_result();

   EXPECT_EQ(nll0, nll1);
}



TEST_F(LikelihoodSerialBinnedDatasetTest, BinnedManualNLL)
{
   data = pdf->generateBinned(*w.var("x"));

   // manually create NLL, ripping all relevant parts from RooAbsPdf::createNLL, except here we also set binnedL = true
   RooArgSet projDeps;
   const char* rangeName = 0;
   const char* addCoefRangeName = 0;
   Int_t numcpu = 1;
   Int_t ext = 2;
   RooFit::MPSplit interl = (RooFit::MPSplit) 0;
   Bool_t cpuAffinity = static_cast<Bool_t>(1);
   Bool_t verbose = 0;
   Int_t splitr = 0;
   Int_t cloneData = 0;
   RooNLLVar nll_manual("nlletje", "-log(likelihood)", *pdf, *data, projDeps, ext, rangeName, addCoefRangeName,
                        numcpu, interl, cpuAffinity, verbose, splitr, cloneData, /*binnedL=*/true);

   likelihood = std::make_shared<RooFit::TestStatistics::RooBinnedL>(pdf, data);
   RooFit::TestStatistics::LikelihoodSerial nll_ts(likelihood, clean_flags, nullptr);

   auto nll0 = nll_manual.getVal();

   nll_ts.evaluate();
   auto nll1 = nll_ts.return_result();

   EXPECT_EQ(nll0, nll1);
}


TEST_F(LikelihoodSerialTest, SimBinned)
{
   // Unbinned pdfs that define template histograms
   w.factory("Gaussian::gA(x[-10,10],-2,3)");
   w.factory("Gaussian::gB(x[-10,10],2,1)");
   w.factory("Uniform::u(x)");

   // Generate template histograms
   RooDataHist *h_sigA = w.pdf("gA")->generateBinned(*w.var("x"), 1000);
   RooDataHist *h_sigB = w.pdf("gB")->generateBinned(*w.var("x"), 1000);
   RooDataHist *h_bkg = w.pdf("u")->generateBinned(*w.var("x"), 1000);

   w.import(*h_sigA, RooFit::Rename("h_sigA"));
   w.import(*h_sigB, RooFit::Rename("h_sigB"));
   w.import(*h_bkg, RooFit::Rename("h_bkg"));

   // Construct L_phys: binned pdf as sum of amplitudes
   w.factory("HistFunc::hf_sigA(x,h_sigA)");
   w.factory("HistFunc::hf_sigB(x,h_sigB)");
   w.factory("HistFunc::hf_bkg(x,h_bkg)");

   w.factory("ASUM::model_A(mu_sig[1,-1,10]*hf_sigA,mu_bkg_A[1,-1,10]*hf_bkg)");
   w.factory("ASUM::model_B(mu_sig*hf_sigB,mu_bkg_B[1,-1,10]*hf_bkg)");

   w.pdf("model_A")->setAttribute("BinnedLikelihood");
   w.pdf("model_B")->setAttribute("BinnedLikelihood");

   // Construct simulatenous pdf
   w.factory("SIMUL::model(index[A,B],A=model_A,B=model_B)");

   // Construct dataset
   pdf = w.pdf("model");
   data = pdf->generate(RooArgSet(*w.var("x"), *w.cat("index")), RooFit::AllBinned());

   nll.reset(pdf->createNLL(*data));

   likelihood = std::make_shared<RooFit::TestStatistics::RooSimultaneousL>(pdf, data);
   RooFit::TestStatistics::LikelihoodSerial nll_ts(likelihood, clean_flags, nullptr);

   auto nll0 = nll->getVal();

   nll_ts.evaluate();
   auto nll1 = nll_ts.return_result();

   EXPECT_EQ(nll0, nll1);
}

TEST_F(LikelihoodSerialTest, BinnedConstrained)
{
   // Unbinned pdfs that define template histograms

   w.factory("Gaussian::g(x[-10,10],0,2)");
   w.factory("Uniform::u(x)");

   // Generate template histograms

   RooDataHist *h_sig = w.pdf("g")->generateBinned(*w.var("x"), 1000);
   RooDataHist *h_bkg = w.pdf("u")->generateBinned(*w.var("x"), 1000);

   w.import(*h_sig, RooFit::Rename("h_sig"));
   w.import(*h_bkg, RooFit::Rename("h_bkg"));

   // Construct binned pdf as sum of amplitudes
   w.factory("HistFunc::hf_sig(x,h_sig)");
   w.factory("HistFunc::hf_bkg(x,h_bkg)");
   w.factory("ASUM::model_phys(mu_sig[1,-1,10]*hf_sig,expr::mu_bkg('1+0.02*alpha_bkg',alpha_bkg[-5,5])*hf_bkg)") ;

   // Construct L_subs: Gaussian subsidiary measurement that constrains alpha_bkg
   w.factory("Gaussian:model_subs(alpha_bkg_obs[0],alpha_bkg,1)") ;

   // Construct full pdf
   w.factory("PROD::model(model_phys,model_subs)") ;

   pdf = w.pdf("model");
   // Construct dataset from physics pdf
   data = w.pdf("model_phys")->generateBinned(*w.var("x"));

   nll.reset(pdf->createNLL(*data, RooFit::GlobalObservables(*w.var("alpha_bkg_obs"))));

   // --------

   auto nll0 = nll->getVal();

   likelihood =
//      std::make_shared<RooFit::TestStatistics::RooBinnedL>(pdf, data, RooFit::GlobalObservables(*w.var("alpha_bkg_obs")));
      std::make_shared<RooFit::TestStatistics::RooUnbinnedL>(pdf, data);
   RooFit::TestStatistics::LikelihoodSerial nll_ts(likelihood, clean_flags, nullptr);

   nll_ts.evaluate();
   auto nll1 = nll_ts.return_result();

   EXPECT_EQ(nll0, nll1);
}

TEST_F(LikelihoodSerialTest, SimUnbinned)
{
   // SIMULTANEOUS FIT OF 2 UNBINNED DATASETS

   w.factory("ExtendPdf::egA(Gaussian::gA(x[-10,10],mA[2,-10,10],s[3,0.1,10]),nA[1000])") ;
   w.factory("ExtendPdf::egB(Gaussian::gB(x,mB[-2,-10,10],s),nB[100])") ;
   w.factory("SIMUL::model(index[A,B],A=egA,B=egB)") ;

   pdf = w.pdf("model");
   // Construct dataset from physics pdf
   data = pdf->generate(RooArgSet(*w.var("x"),*w.cat("index")));

   nll.reset(pdf->createNLL(*data));

   // --------

   auto nll0 = nll->getVal();

   likelihood =
      std::make_shared<RooFit::TestStatistics::RooSimultaneousL>(pdf, data);
   RooFit::TestStatistics::LikelihoodSerial nll_ts(likelihood, clean_flags, nullptr);

   nll_ts.evaluate();
   auto nll1 = nll_ts.return_result();

   EXPECT_EQ(nll0, nll1);
}

TEST_F(LikelihoodSerialTest, SimUnbinnedNonExtended)
{
   // SIMULTANEOUS FIT OF 2 UNBINNED DATASETS

   w.factory("ExtendPdf::egA(Gaussian::gA(x[-10,10],mA[2,-10,10],s[3,0.1,10]),nA[1000])") ;
   w.factory("ExtendPdf::egB(Gaussian::gB(x,mB[-2,-10,10],s),nB[100])") ;
   w.factory("SIMUL::model(index[A,B],A=gA,B=gB)") ;

   RooDataSet* dA = w.pdf("gA")->generate(*w.var("x"),1) ;
   RooDataSet* dB = w.pdf("gB")->generate(*w.var("x"),1) ;
   w.cat("index")->setLabel("A") ;
   dA->addColumn(*w.cat("index")) ;
   w.cat("index")->setLabel("B") ;
   dB->addColumn(*w.cat("index")) ;

   data = (RooDataSet*) dA->Clone() ;
   static_cast<RooDataSet*>(data)->append(*dB) ;

   pdf = w.pdf("model");

   nll.reset(pdf->createNLL(*data));

   likelihood =
      std::make_shared<RooFit::TestStatistics::RooSimultaneousL>(pdf, data);
   RooFit::TestStatistics::LikelihoodSerial nll_ts(likelihood, clean_flags, nullptr);

   auto nll0 = nll->getVal();

   nll_ts.evaluate();
   auto nll1 = nll_ts.return_result();

   EXPECT_EQ(nll0, nll1);
}


TEST_F(LikelihoodSerialTest, SimBinnedConstrained)
{
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

   nll.reset(pdf->createNLL(*data, RooFit::GlobalObservables(RooArgSet(*w.var("alpha_bkg_obs_A"),*w.var("alpha_bkg_obs_B")))));

   // --------

   auto nll0 = nll->getVal();

   likelihood =
//      std::make_shared<RooFit::TestStatistics::RooSimultaneousL>(pdf, data, RooFit::GlobalObservables(RooArgSet(*w.var("alpha_bkg_obs_A"),*w.var("alpha_bkg_obs_B"))));
      std::make_shared<RooFit::TestStatistics::RooSimultaneousL>(pdf, data);
   RooFit::TestStatistics::LikelihoodSerial nll_ts(likelihood, clean_flags, nullptr);

   nll_ts.evaluate();
   auto nll1 = nll_ts.return_result();

   EXPECT_EQ(nll0, nll1);
}
