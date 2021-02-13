#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TRandom.h"
#include "TH1.h"
#include "Math/QuantFuncMathCore.h"

#include <iostream>

#include "gtest/gtest.h"

bool testTEfficiency_vs_TGA(int nexp = 1000, TEfficiency::EStatOption statOpt = TEfficiency::kBUniform,
                            bool mode = true, bool central = false)
{
   gRandom->SetSeed(111);

   bool ok = true;
   for (int i = 0; i < nexp; ++i) {

      // if (i>0 && i%500==0) std::cout << i << std::endl;

      // loop on the experiment
      double n = int(std::abs(gRandom->BreitWigner(0, 5))) + 1;
      double cut = ROOT::Math::beta_quantile(gRandom->Rndm(), 0.5, 0.5);
      double k = int(cut * n);
      TH1D *h1 = new TH1D("h1", "h1", 1, 0, 1);
      TH1D *h2 = new TH1D("h2", "h2", 1, 0, 1);
      h1->SetDirectory(0);
      h2->SetDirectory(0);
      h1->SetBinContent(1, k);
      h2->SetBinContent(1, n);

      // test the various option : case no mode (average) and shortes (no central)
      // cannot be done with TGraphAsymmErrors. ROOT-10324 is missing mode central
      // that is now fixed
      TGraphAsymmErrors *g = new TGraphAsymmErrors();
      if (statOpt == TEfficiency::kBUniform && mode && !central)
         g->BayesDivide(h1, h2); // this is mode not central for uniform prior
      else if (statOpt == TEfficiency::kBUniform && !mode && central)
         g->Divide(h1, h2, "cl=0.683 b(1,1)"); // default is no mode and central
      else if (statOpt == TEfficiency::kBUniform && mode && central)
         g->Divide(h1, h2, "cl=0.683 b(1,1) mode central");
      else if (statOpt == TEfficiency::kBUniform && !mode && !central)
         g->Divide(h1, h2, "cl=0.683 b(1,1) shortest");

      else if (statOpt == TEfficiency::kBJeffrey && mode && !central)
         g->Divide(h1, h2, "cl=0.683 b(0.5,0.5) mode");
      else if (statOpt == TEfficiency::kBJeffrey && !mode && central)
         g->Divide(h1, h2, "cl=0.683 b(0.5,0.5) central"); // adding central is actually useless here but so we test this
      else if (statOpt == TEfficiency::kBJeffrey && mode && central)
         g->Divide(h1, h2, "cl=0.683 b(0.5,0.5) mode central");
      else if (statOpt == TEfficiency::kBJeffrey && !mode && !central)
         g->Divide(h1, h2, "cl=0.683 b(0.5,0.5) shortest");

      else if (statOpt == TEfficiency::kFCP)
         g->Divide(h1, h2, "cl=0.683 cp");
      else if (statOpt == TEfficiency::kFAC)
         g->Divide(h1, h2, "cl=0.683 ac");
      else if (statOpt == TEfficiency::kFFC)
         g->Divide(h1, h2, "cl=0.683 fc");
      else if (statOpt == TEfficiency::kFWilson)
         g->Divide(h1, h2, "cl=0.683 w");
      else if (statOpt == TEfficiency::kFNormal)
         g->Divide(h1, h2, "cl=0.683 n");
      else {
         std::cout << "\n\nERROR:  invalid statistic options - exit \n" << std::endl;
         return false;
      }
      double eff = g->GetY()[0];
      double eu = g->GetEYhigh()[0];
      double el = g->GetEYlow()[0];

      TEfficiency *e = new TEfficiency(*h1, *h2);
      // eff->SetPosteriorMode(false);
      e->SetStatisticOption(statOpt);
      e->SetPosteriorMode(mode);
      if (central)
         e->SetCentralInterval();
      else
         e->SetShortestInterval();

      e->SetConfidenceLevel(0.683);

      double eff2 = e->GetEfficiency(1);
      double el2 = e->GetEfficiencyErrorLow(1);
      double eu2 = e->GetEfficiencyErrorUp(1);

      if (eff2 != eff) {
         std::cerr << "Different efficiency " << eff2 << "  vs  " << eff << std::endl;
         ok = false;
      }
      if (el2 != el) {
         std::cerr << "Different low error " << el2 << "  vs  " << el << std::endl;
         ok = false;
      }
      if (eu2 != eu) {
         std::cerr << "Different up  error " << eu2 << "  vs " << eu << std::endl;
         ok = false;
      }
      EXPECT_EQ(ok, true);
      if (!ok) {
         std::cerr << "Iteration " << i << ":\t Error for (k,n) " << int(k) << " , " << int(n) << std::endl;
         break;
      }
      delete e;

      delete h1;
      delete h2;
      delete g;
   }

   if (ok)
      std::cout << "Comparison TEfficiency-TGraphAsymError :  OK for nevt = " << nexp << std::endl;
   else
      std::cout << "Comparison TEfficiency-TGraphAsymError :  FAILED ! " << std::endl;

   return ok;
}

bool testNormal()
{
   float tol = 1e-3;
   bool ok = true;

   // test the 95% confidence intervals
   // taken from: http://www.measuringusability.com/wald.htm
   //
   // format: (k,n) -> lower bound, upper bound
   // (0,0) -> 0, 1
   // (3,7) -> 0.062, 0.795
   // (0,8) -> 0, 0
   // (3,12) -> 0.005, 0.495
   // (2,14) -> 0, 0.326
   // (5,18) -> 0.071, 0.485
   // (15,30) -> 0.321, 0.679
   // (10,10) -> 1, 1

   const int max = 8;
   Double_t k[max] = {0, 3, 0, 3, 2, 5, 15, 10};
   Double_t n[max] = {0, 7, 8, 12, 14, 18, 30, 10};
   Double_t low[max] = {0, 0.062, 0, 0.005, 0, 0.071, 0.321, 1};
   Double_t up[max] = {1, 0.795, 0, 0.495, 0.326, 0.485, 0.679, 1};

   for (int i = 0; i < max; ++i) {
      if (fabs(TEfficiency::Normal(n[i], k[i], 0.95, true) - up[i]) > tol) {
         std::cerr << "different upper bound for Normal interval (" << k[i] << "," << n[i] << ")" << std::endl;
         std::cerr << "got: " << TEfficiency::Normal(n[i], k[i], 0.95, true) << " expecting: " << up[i] << std::endl;
         ok = false;
      }
      if (fabs(TEfficiency::Normal(n[i], k[i], 0.95, false) - low[i]) > tol) {
         std::cerr << "different lower bound for Normal interval (" << k[i] << "," << n[i] << ")" << std::endl;
         std::cerr << "got: " << TEfficiency::Normal(n[i], k[i], 0.95, false) << " expecting: " << low[i] << std::endl;
         ok = false;
      }
   }

   std::cout << "confidence interval for Normal ";
   (ok) ? std::cout << "OK" : std::cout << "FAILED";
   std::cout << std::endl;
   EXPECT_EQ(ok, true);
   return ok;
}

bool testWilson()
{
   float tol = 1e-3;
   bool ok = true;

   // test the 95% confidence intervals
   // taken from: http://www.measuringusability.com/wald.htm
   //
   // format: (k,n) -> lower bound, upper bound
   // (0,0) -> 0, 1
   // (3,7) -> 0.158, 0.750
   // (0,8) -> 0, 0.324
   // (3,12) -> 0.089, 0.532
   // (2,14) -> 0.040, 0.399
   // (5,18) -> 0.125, 0.509
   // (15,30) -> 0.332, 0.669
   // (10,10) -> 0.722, 1.000

   const int max = 8;
   Double_t k[max] = {0, 3, 0, 3, 2, 5, 15, 10};
   Double_t n[max] = {0, 7, 8, 12, 14, 18, 30, 10};
   Double_t low[max] = {0, 0.158, 0, 0.089, 0.040, 0.125, 0.332, 0.722};
   Double_t up[max] = {1, 0.750, 0.324, 0.532, 0.399, 0.509, 0.669, 1};

   for (int i = 0; i < max; ++i) {
      if (fabs(TEfficiency::Wilson(n[i], k[i], 0.95, true) - up[i]) > tol) {
         std::cerr << "different upper bound for Wilson interval (" << k[i] << "," << n[i] << ")" << std::endl;
         std::cerr << "got: " << TEfficiency::Wilson(n[i], k[i], 0.95, true) << " expecting: " << up[i] << std::endl;
         ok = false;
      }
      if (fabs(TEfficiency::Wilson(n[i], k[i], 0.95, false) - low[i]) > tol) {
         std::cerr << "different lower bound for Wilson interval (" << k[i] << "," << n[i] << ")" << std::endl;
         std::cerr << "got: " << TEfficiency::Wilson(n[i], k[i], 0.95, false) << " expecting: " << low[i] << std::endl;
         ok = false;
      }
   }

   std::cout << "confidence interval for Wilson ";
   (ok) ? std::cout << "OK" : std::cout << "FAILED";
   std::cout << std::endl;
   EXPECT_EQ(ok, true);
   return ok;
}

bool testFeldmanCousins()
{
   float tol = 1e-3;
   bool ok = true;

   // test the 95% confidence intervals
   // taken from: http://people.na.infn.it/~lista/cgi/binomial/binomial.pl
   //
   // format: (k,n) -> lower bound, upper bound
   // (0,0) -> 0, 1
   // (3,7) -> 0.129, 0.775
   // (0,8) -> 0, 0.321
   // (3,12) -> 0.072, 0.548
   // (2,14) -> 0.026, 0.418
   // (5,18) -> 0.106, 0.531
   // (15,30) -> 0.324, 0.676
   // (10,10) -> 0.733, 1.000

   const int max = 8;
   Double_t k[max] = {0, 3, 0, 3, 2, 5, 15, 10};
   Double_t n[max] = {0, 7, 8, 12, 14, 18, 30, 10};
   Double_t low[max] = {0, 0.129, 0, 0.072, 0.026, 0.106, 0.324, 0.733};
   Double_t up[max] = {1, 0.775, 0.321, 0.548, 0.418, 0.531, 0.676, 1};

   for (int i = 0; i < max; ++i) {
      if (fabs(TEfficiency::FeldmanCousins(n[i], k[i], 0.95, true) - up[i]) > tol) {
         std::cerr << "different upper bound for Feldman-Cousins interval (" << k[i] << "," << n[i] << ")" << std::endl;
         std::cerr << "got: " << TEfficiency::FeldmanCousins(n[i], k[i], 0.95, true) << " expecting: " << up[i]
                   << std::endl;
         ok = false;
      }
      if (fabs(TEfficiency::FeldmanCousins(n[i], k[i], 0.95, false) - low[i]) > tol) {
         std::cerr << "different lower bound for Feldman-Cousins interval (" << k[i] << "," << n[i] << ")" << std::endl;
         std::cerr << "got: " << TEfficiency::FeldmanCousins(n[i], k[i], 0.95, false) << " expecting: " << low[i]
                   << std::endl;
         ok = false;
      }
   }

   std::cout << "confidence interval for Feldman-Cousins ";
   (ok) ? std::cout << "OK" : std::cout << "FAILED";
   std::cout << std::endl;
   EXPECT_EQ(ok, true);
   return ok;
}

bool testClopperPearson()
{
   float tol = 1e-3;
   bool ok = true;

   // test the 95% confidence intervals
   // taken from: http://people.na.infn.it/~lista/cgi/binomial/binomial.pl
   //
   // format: (k,n) -> lower bound, upper bound
   // (0,0) -> 0, 1
   // (3,7) -> 0.099, 0.816
   // (0,8) -> 0, 0.369
   // (3,12) -> 0.055, 0.572
   // (2,14) -> 0.018, 0.428
   // (5,18) -> 0.097, 0.535
   // (15,30) -> 0.313, 0.687
   // (10,10) -> 0.692, 1.000

   const int max = 8;
   Double_t k[max] = {0, 3, 0, 3, 2, 5, 15, 10};
   Double_t n[max] = {0, 7, 8, 12, 14, 18, 30, 10};
   Double_t low[max] = {0, 0.099, 0, 0.055, 0.018, 0.097, 0.313, 0.692};
   Double_t up[max] = {1, 0.816, 0.369, 0.572, 0.428, 0.535, 0.687, 1};

   for (int i = 0; i < max; ++i) {
      if (fabs(TEfficiency::ClopperPearson(n[i], k[i], 0.95, true) - up[i]) > tol) {
         std::cerr << "different upper bound for Clopper-Pearson interval (" << k[i] << "," << n[i] << ")" << std::endl;
         std::cerr << "got: " << TEfficiency::ClopperPearson(n[i], k[i], 0.95, true) << " expecting: " << up[i]
                   << std::endl;
         ok = false;
      }
      if (fabs(TEfficiency::ClopperPearson(n[i], k[i], 0.95, false) - low[i]) > tol) {
         std::cerr << "different lower bound for Clopper=Pearson interval (" << k[i] << "," << n[i] << ")" << std::endl;
         std::cerr << "got: " << TEfficiency::ClopperPearson(n[i], k[i], 0.95, false) << " expecting: " << low[i]
                   << std::endl;
         ok = false;
      }
   }

   std::cout << "confidence interval for Clopper-Pearson ";
   (ok) ? std::cout << "OK" : std::cout << "FAILED";
   std::cout << std::endl;
   EXPECT_EQ(ok, true);
   return ok;
}

bool testJeffreyPrior()
{
   float tol = 1e-3;
   bool ok = true;

   // test the 95% confidence intervals
   // taken from:
   // "Interval Estimation for a Binomial Proportion" Brown, Cai, DasGupta
   // Table 5
   //
   // format: (k,n) -> lower bound, upper bound
   // (0,0) -> 0.002, 0.998
   // (3,7) -> 0.139, 0.766
   // (0,8) -> 0, 0.262
   // (3,12) -> 0.076, 0.529
   // (2,14) -> 0.031, 0.385
   // (5,18) -> 0.115, 0.506
   // (15,30) -> 0.328, 0.672
   // (10,10) -> 0.783, 1.000
   //
   // alpha = k + 0.5
   // beta = n - k + 0.5

   const int max = 8;
   Double_t k[max] = {0, 3, 0, 3, 2, 5, 15, 10};
   Double_t n[max] = {0, 7, 8, 12, 14, 18, 30, 10};
   Double_t low[max] = {0.002, 0.139, 0, 0.076, 0.031, 0.115, 0.328, 0.783};
   Double_t up[max] = {0.998, 0.766, 0.262, 0.529, 0.385, 0.506, 0.672, 1};

   Double_t alpha, beta;
   for (int i = 0; i < max; ++i) {
      alpha = k[i] + 0.5;
      beta = n[i] - k[i] + 0.5;

      if (fabs(TEfficiency::BetaCentralInterval(0.95, alpha, beta, true) - up[i]) > tol) {
         std::cerr << "different upper bound for Jeffrey interval (" << k[i] << "," << n[i] << ")" << std::endl;
         std::cerr << "got: " << TEfficiency::BetaCentralInterval(0.95, alpha, beta, true) << " expecting: " << up[i]
                   << std::endl;
         ok = false;
      }
      if (fabs(TEfficiency::BetaCentralInterval(0.95, alpha, beta, false) - low[i]) > tol) {
         std::cerr << "different lower bound for Jeffrey interval (" << k[i] << "," << n[i] << ")" << std::endl;
         std::cerr << "got: " << TEfficiency::BetaCentralInterval(0.95, alpha, beta, false) << " expecting: " << low[i]
                   << std::endl;
         ok = false;
      }
   }

   std::cout << "confidence interval for Jeffrey prior ";
   (ok) ? std::cout << "OK" : std::cout << "FAILED";
   std::cout << std::endl;
   EXPECT_EQ(ok, true);
   return ok;
}

void testConsistencyWithTGraph()
{
   // check consistency between TEfficiency and TGraphAsymmErrors
   std::cout << "uniform prior with mode: ";
   testTEfficiency_vs_TGA(1000, TEfficiency::kBUniform, true, false);
   std::cout << "uniform prior with mean and central: ";
   testTEfficiency_vs_TGA(1000, TEfficiency::kBUniform, false, true);
   std::cout << "uniform prior with mode and central interval: ";
   testTEfficiency_vs_TGA(1000, TEfficiency::kBUniform, true, true); // case of ROOT-10324
   std::cout << "uniform prior with mean and shortest ";
   testTEfficiency_vs_TGA(1000, TEfficiency::kBUniform, false, false);
   std::cout << "Jeffrey prior with mode: ";
   testTEfficiency_vs_TGA(1000, TEfficiency::kBJeffrey, true, false);
   std::cout << "Jeffrey prior with mean and central: ";
   testTEfficiency_vs_TGA(1000, TEfficiency::kBJeffrey, false, true);
   std::cout << "Jeffrey prior with mode and central: ";
   testTEfficiency_vs_TGA(1000, TEfficiency::kBJeffrey, true, true);
   std::cout << "Jeffrey prior with mean and shortest: ";
   testTEfficiency_vs_TGA(1000, TEfficiency::kBJeffrey, false, false);
   std::cout << "Clopper-Pearson: ";
   testTEfficiency_vs_TGA(1000, TEfficiency::kFCP);
   std::cout << "Agresti-Coull: ";
   testTEfficiency_vs_TGA(1000, TEfficiency::kFAC);
   std::cout << "Feldman-Cousin: ";
   testTEfficiency_vs_TGA(1000, TEfficiency::kFFC);
   std::cout << "Wilson: ";
   testTEfficiency_vs_TGA(1000, TEfficiency::kFWilson);
   std::cout << "Normal: ";
   testTEfficiency_vs_TGA(1000, TEfficiency::kFNormal);
}

void testConfIntervals()
{
   // check confidence intervals for a few points
   testClopperPearson();
   testNormal();
   testWilson();
   testFeldmanCousins();
   testJeffreyPrior();
}

TEST(TFEfficiency, ConfIntervals)
{
   testConfIntervals();
}
TEST(TFEfficiency, ConsistencyWithTGraph)
{
   testConsistencyWithTGraph();
}