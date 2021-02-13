#include "gtest/gtest.h"


#include "TKDE.h"
#include "TString.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TVirtualPad.h"
#include "TF1.h"
#include "TH1.h"
#include "Math/DistFuncMathCore.h"


struct TestKDE {

   // using n < 10000 use a default nbins=100 for KDE
   int n = 10000;
   bool makePlot = true;
   bool useSetters = false;
   bool binned = false;
   bool adaptive = false;
   bool mirroring = false;
   bool debug = false;
   bool dataInRange = false;
   bool independentData = true;
   int nbins = 100;
   int seed = 1111;
   const char *mirror = "mirrorBoth";

   // sigma of the two gaussian used for generating the data
   double sigma1 = 1;
   double sigma2 = 7;

   double pcut = 0.01; // cut on minimum of Chi2 and KS p values
   double pcut1 = 1.E-4; // cut on min p value
   double pcut2 = 0.05;  // cut on second p value

   std::vector<double> data;  // kde data
   std::vector<double> data2;  // second data set for histogram

   std::vector<double> xtest;

   int nTestPoints = 21;
   double values1[21];
   double values2[21];

   double pval;
   double pval2;

   TKDE * Create() {
      TRandom3 r(seed);
      data.resize(n);
      data2.resize(n);
      if (!dataInRange) {
         for (int i = 0; i < n; ++i) {
            data[i] =  (r.Rndm() < 0.2) ? r.Gaus(10, sigma1) : r.Gaus(10, sigma2);
            data2[i] = (r.Rndm() < 0.2 ) ? r.Gaus(10, sigma1) : r.Gaus(10, sigma2);
         }
      } else {
         int i = 0;
         while (i < 2*n) {
            double x = (r.Rndm() < 0.2) ? r.Gaus(10, sigma1) : r.Gaus(10, sigma2);
            if (x >= 0 && x < 20.) {
               if (i < n)
                  data[i] = x;
               else
                  data2[i-n] = x;
               i++;
            }
         }
      }
      const char *iterationType = (adaptive) ? "Adaptive" : "Fixed";
      const char *binningType = (binned) ? "ForcedBinning" : "Unbinned";
      const char *mirrorType = (mirroring) ? mirror : "noMirror";
      TString opt;
      if (useSetters)
         opt = "KernelType:Gaussian;Iteration:Fixed;Mirror:noMirror;Binning:Unbinned";
      else
         opt = TString::Format("KernelType:Gaussian;Iteration:%s;Mirror:%s;Binning:%s",iterationType, mirrorType, binningType);

      assert((int)data.size() == n);
      TKDE *kde = new TKDE(n, data.data(), 0., 20., opt, 1);
      if (useSetters) {
         if (adaptive) {
            kde->SetIteration(TKDE::kAdaptive);
         }
         if (binned) {
            // use relaxed binned mode
            kde->SetBinning(TKDE::kRelaxedBinning);
            kde->SetUseBinsNEvents(10); // this is like using forced binning
            kde->SetNBins(nbins);  // when nbins=100 (using setters or not should give same results)->
         }
         if (mirroring) {
            kde->SetMirror(TKDE::kMirrorBoth);
         }
      }
      xtest.resize(20);
      TAxis a(nTestPoints,-0.5,20.5);
      for (size_t i = 0; i < xtest.size(); ++i) {
         xtest[i] = a.GetBinCenter(i+1);
         values1[i] = (*kde)(xtest[i] );
      }
      if (debug) {
         std::cout << "Before writing,  kde values :  ";
         for (size_t i = 0; i < xtest.size(); ++i) {
            std::cout << values1[i];
            if (i <  xtest.size()-1) std::cout << " , ";
            else std::cout << std::endl;
         }
      }

      if (makePlot) {
         //auto c1 = new TCanvas();
         kde->Draw();
         kde->GetDrawnFunction()->SetLineColor(kBlue);
      }
      return kde;
   }

   void CompareWithHist(TString name = "tkde") {
      auto kde = Create();
      int nhistbins = 100;
      auto h1 = new TH1D("h1", "h1", nhistbins, 0., 20.);
      // neglect correction for sigma=1 gaussian
      double corrNorm = 2.* ROOT::Math::normal_cdf(0, sigma2, 10);
      double a1 = 0.2;
      double a2 = 0.8;
      if (dataInRange) {
         // correction will increase a1 since events of a2 will be outside range
         a1 = 0.2 / (0.2 + 0.8 * (1.-corrNorm));
         a2 = 1. - a1;
      }
      auto pdf = new TF1("pdf", "[a1]*ROOT::Math::normal_pdf(x,[s1],10)+[a2]*ROOT::Math::normal_pdf(x,[s2],10)", 0, 20);
      pdf->SetParameters(a1, a2, sigma1,sigma2);

      if (independentData) {
         // should use half of data for histogram since
         // we do not count statistical error from TKDE?
         //for (int i = 0; i < 0.5*n; ++i)
         //   h1->Fill(data2[i]);
         for (int i = 1; i < h1->GetNbinsX(); ++i) {
            h1->SetBinContent(i, pdf->Eval(h1->GetBinCenter(i))*h1->GetBinWidth(i)*n);
         }
      } else {
         for (int i = 0; i < n; ++i)
            h1->Fill(data[i]);
      }

      h1->Sumw2();
      // scale to n - if data not in range scale will be also correct
      if (debug)
         std::cout << dataInRange << " ..." << n << " true integral - " << h1->Integral() << std::endl;
      double nhist = h1->Integral();
       //(dataInRange || binned) ? h1->Integral() : n;
      h1->Scale(1. / nhist, "width");

      auto fkde = kde->GetFunction(nhistbins);
      if (debug) {
         auto fkde1000 = kde->GetFunction(1000);
         std::cout << "KDE integral is " << fkde->GetHistogram()->Integral() * fkde->GetHistogram()->GetBinWidth(1)
                   << " with more bins " << fkde1000->GetHistogram()->Integral() * fkde1000->GetHistogram()->GetBinWidth(1) << std::endl;
      }

      if (makePlot) {
         //auto c1 = new TCanvas();
         kde->Draw();
         h1->SetLineColor(kBlue);
         h1->Draw("SAME");
         TString fname = name + "_hist.pdf";
         if (gPad) gPad->SaveAs(fname);
      }


      // make chi2 test
      //double pvalChi2  = h1->Chi2Test(h2,"WW P");
      // note Chisquare called SetParameters and delete histogram returned from the function
      double chi2 = h1->Chisquare(fkde);
      double pvalChi2 = TMath::Prob(chi2, nhistbins);

      auto h2 = fkde->GetHistogram();
      h2->Sumw2();
      // set correct bin error
      for (int i = 1; i <= h2->GetNbinsX() + 1; ++i)
         h2->SetBinError(i, sqrt(h2->GetBinContent(i) / (n * h2->GetBinWidth(i))));

      // compute KS test
      double pvalKS = h1->KolmogorovTest(h2);
      std::cout << "CompareWithHist: " << name << " - Chi2 test = " << chi2 << "/" << nhistbins << " pvalue = " << pvalChi2 <<
                   "  KS test p value = " << pvalKS << std::endl;
      pval = std::min(pvalKS, pvalChi2);
      pval2 = std::max(pvalKS, pvalChi2);

      delete kde;
      delete h1;
   }

   bool IsPValid() {
      bool valid = pval > pcut;
      if (!valid)
         std::cerr << "Error : pvalue " << pval << " is smaller than " << pcut << std::endl;
      return valid;
   }
   bool ArePsValid(){
      bool valid = pval > pcut1 && pval2 > pcut2;
      if (!valid) {
         if (pval < pcut1)
            std::cerr << "Error : pvalue " << pval << " is smaller than " << pcut << std::endl;
         if (pval2 < pcut2)
            std::cerr << "Error : pvalue " << pval2 << " is smaller than " << pcut2 << std::endl;
      }
      return valid;
   }

   void Write(TString name = "tkde") {
      TFile *f = new TFile("test_tkde.root", "RECREATE");

      auto kde = Create();

      kde->Write(name);
      f->Close();
   }

   void Read(TString name = "tkde") {
      TFile *f = new TFile("test_tkde.root");
      auto kde = (TKDE*) f->Get(name);
      //kde->Dump();
      for (size_t i = 0; i < xtest.size(); ++i) values2[i] = (*kde)(xtest[i] );

      if (debug) {
         std::cout << "After reading,  kde values : ";
         for (size_t i = 0; i < xtest.size(); ++i) {
            std::cout << values2[i];
            if (i <  xtest.size()-1) std::cout << " , ";
            else std::cout << std::endl;
         }
      }

      if (makePlot) {
         kde->Draw("SAME");
         TString gfname = name + ".pdf";
         gPad->SaveAs(gfname);
      }
      //f->Close();
      delete f;
   }
};

/// Hstogram comparison tests
/// In this test we compare the TKDE with an histogram
/// filled with the same type of data. A GoF test (chi2 and KS) is applied to
/// check for compatibility
TEST(TKDE, tkde_hist)
{
   TestKDE t;
   t.n = 20000;
   t.CompareWithHist("tkde");
   EXPECT_TRUE(t.IsPValid());
   t.dataInRange = true;
   t.CompareWithHist("tkde_in");
   // increase p -value because without mirroring this is expected to fail
   t.pcut = 1.E-5;
   EXPECT_TRUE(t.IsPValid());
}
TEST(TKDE, tkde_hist_adaptive)
{
   TestKDE t;
   // reduce n here since unbinned adaptive is slow
   t.n = 2000;
   t.adaptive = true;
   t.CompareWithHist("tkde_adaptive");
   EXPECT_TRUE(t.IsPValid());
   // test also default constructors + usage of setters functions
   t.useSetters = true;
   t.dataInRange = true;
   t.CompareWithHist("tkde_adaptive_in");
   // this is
   EXPECT_TRUE(t.IsPValid());
}
TEST(TKDE, tkde_hist_binned)
{
   TestKDE t;
   t.dataInRange = false;
   t.adaptive = false;
   t.binned = true;
   // this test will fail miserably if sigma2 is large since
   // binning does not use events outside range
   // data.InRange=false is not doing since the bin data will
   // consider only data in the range !!!
   t.sigma2 = 3;
   t.n = 100000;
   t.nbins = 5000;
   // with binned not adaptive it is good for high statistics
   t.CompareWithHist("tkde_binned");
   EXPECT_TRUE(t.ArePsValid());
   t.useSetters = true;
   t.dataInRange = true;
   t.CompareWithHist("tkde_binned_in");
   EXPECT_TRUE(t.ArePsValid());
}
TEST(TKDE, tkde_hist_adaptive_binned)
{
   TestKDE t;
   t.dataInRange = true;
   t.nbins = t.n/10;
   //t.nbins = 3000;
   t.adaptive = true;
   t.binned = true;
   t.CompareWithHist("tkde_adaptive_binned");
   EXPECT_TRUE(t.IsPValid());
   t.useSetters = true;
   t.dataInRange = false;
   t.CompareWithHist("tkde_adaptive_binned_all");
   EXPECT_TRUE(t.IsPValid());
}

TEST(TKDE, tkde_hist_mirror)
{
   TestKDE t;
   t.mirroring = true;
   t.n = 10000;
   // you need data in range for mirroring !!!
   t.dataInRange = true;
   t.CompareWithHist("tkde_mirrorLR");
   // for mirror case KS test is not so good but Chi2 is better
   // test both p-values
   EXPECT_TRUE(t.ArePsValid());
   // t.useSetters = true;
   // t.CompareWithHist("tkde_mirrorLR_2");
   // EXPECT_TRUE(t.IsPValid());
   // check other mirrors types:
   t.useSetters = false;
   t.mirror = "mirrorLeft";
   t.CompareWithHist("tkde_mirrorL");
   EXPECT_TRUE(t.ArePsValid());
   //EXPECT_TRUE(t.IsPValid());
   t.mirror = "mirrorRight";
   t.CompareWithHist("tkde_mirrorR");
   EXPECT_TRUE(t.ArePsValid());
   //EXPECT_TRUE(t.IsPValid());
   // for asym case also data should be in range
   // asym does not work very well
   // increase p value cut
   t.pcut = 1.E-4;
   t.mirror = "mirrorAsymBoth";
   t.CompareWithHist("tkde_mirrorALR");
   EXPECT_TRUE(t.IsPValid());
   t.mirror = "mirrorAsymLeft";
   t.CompareWithHist("tkde_mirrorAL");
   EXPECT_TRUE(t.IsPValid());
   t.mirror = "mirrorAsymRight";
   t.CompareWithHist("tkde_mirrorAR");
   EXPECT_TRUE(t.IsPValid());
   t.mirror = "mirrorRightAsymLeft";
   t.CompareWithHist("tkde_mirrorRAL");
   EXPECT_TRUE(t.IsPValid());
   t.mirror = "mirrorLeftAsymRight";
   t.CompareWithHist("tkde_mirrorLAR");
   EXPECT_TRUE(t.IsPValid());
}

TEST(TKDE, tkde_hist_adaptive_mirror)
{
   // when using adaptive mirror a much
   // smaller bandwidth is used. Not sure this is correct
   // seee line TKDE.cxx:770
   TestKDE t;
   t.n = 2000;
   t.adaptive = true;
   t.mirroring = true;
   // you need data in range for mirroring !!!
   t.dataInRange = true;
   t.CompareWithHist("tkde_adaptive_mirrorLR");
   EXPECT_TRUE(t.IsPValid());
   // t.useSetters = true;
   // t.CompareWithHist("tkde_adaptive_mirrorLR_2");
   // EXPECT_TRUE(t.IsPValid());
   // check other mirrors types:
   t.useSetters = false;
   t.mirror = "mirrorLeft";
   t.CompareWithHist("tkde_adaptive_mirrorL");
   EXPECT_TRUE(t.IsPValid());
   t.mirror = "mirrorRight";
   t.CompareWithHist("tkde_adaptive_mirrorR");
   EXPECT_TRUE(t.IsPValid());
   // for asym case data can be outside range
   t.mirror = "mirrorAsymBoth";
   t.CompareWithHist("tkde_adaptive_mirrorALR");
   EXPECT_TRUE(t.IsPValid());
   t.mirror = "mirrorAsymLeft";
   t.CompareWithHist("tkde_adaptive_mirrorAL");
   EXPECT_TRUE(t.IsPValid());
   t.mirror = "mirrorAsymRight";
   t.CompareWithHist("tkde_adaptive_mirrorAR");
   EXPECT_TRUE(t.IsPValid());
   t.mirror = "mirrorRightAsymLeft";
   t.CompareWithHist("tkde_adaptive_mirrorRAL");
   EXPECT_TRUE(t.IsPValid());
   t.mirror = "mirrorLeftAsymRight";
   t.CompareWithHist("tkde_adaptive_mirrorLAR");
   EXPECT_TRUE(t.IsPValid());
}

TEST(TKDE, tkde_hist_binned_mirror)
{
   TestKDE t;
   t.n = 9000;
   t.nbins = 100;
   t.binned = true;
   t.mirroring = true;
   t.adaptive = true;
   // you need data in range for mirroring !!!
   t.dataInRange = true;
   t.CompareWithHist("tkde_binned_mirrorLR");
   EXPECT_TRUE(t.IsPValid());
   t.useSetters = true;
   //note using setters use nbins bins instead of 1000
   t.CompareWithHist("tkde_binned_mirrorLR_2");
   EXPECT_TRUE(t.IsPValid());
   // check other mirrors types:
   t.useSetters = false;
   t.mirror = "mirrorLeft";
   t.CompareWithHist("tkde_binned_mirrorL");
   EXPECT_TRUE(t.IsPValid());
   t.mirror = "mirrorRight";
   t.CompareWithHist("tkde_binned_mirrorR");
   EXPECT_TRUE(t.IsPValid());
   // for asym case data can be outside range
   t.mirror = "mirrorAsymBoth";
   t.CompareWithHist("tkde_binned_mirrorALR");
   EXPECT_TRUE(t.IsPValid());
   t.mirror = "mirrorAsymLeft";
   t.CompareWithHist("tkde_binned_mirrorAL");
   EXPECT_TRUE(t.IsPValid());
   t.mirror = "mirrorAsymRight";
   t.CompareWithHist("tkde_binned_mirrorAR");
   EXPECT_TRUE(t.IsPValid());
   t.mirror = "mirrorRightAsymLeft";
   t.CompareWithHist("tkde_binned_mirrorRAL");
   EXPECT_TRUE(t.IsPValid());
   t.mirror = "mirrorLeftAsymRight";
   t.CompareWithHist("tkde_binned_mirrorLAR");
   EXPECT_TRUE(t.IsPValid());
}

/// IO tests
/// In this test we compare the value before writing and after reading of the TKDE
TEST(TKDE, tkde_io)
{
   TestKDE t;
   t.Write();
   t.Read();

   for (size_t i = 0; i < t.xtest.size(); ++i) {
      EXPECT_DOUBLE_EQ(t.values1[i], t.values2[i]);
      // double delta = 1.E-15;
      // EXPECT_NEAR(t.values1[0], t.values2[0], delta);
   }
}

TEST(TKDE, tkde_io_adaptive)
{
   TString name = "tkde_adaptive";
   TestKDE t;
   t.adaptive = true;
   t.n = 2000;
   t.Write(name);
   t.Read(name);
   double delta = 1.E-3;
   for (size_t i = 0; i < t.xtest.size(); ++i) {
      EXPECT_NEAR(t.values1[i], t.values2[i], delta);
      //EXPECT_DOUBLE_EQ(t.values1[i], t.values2[i]);
   }
}

TEST(TKDE, tkde_io_binned)
{
   TestKDE t;
   t.binned = true;
   TString name = "tkde_binned";
   t.Write(name);
   t.Read(name);
   double delta = 1.E-3;
   for (size_t i = 0; i < t.xtest.size(); ++i) {
      EXPECT_NEAR(t.values1[i], t.values2[i], delta);
   }
}

TEST(TKDE, tkde_io_mirror)
{
   TestKDE t;
   t.mirroring = true;
   TString name = "tkde_mirror";
   t.Write(name);
   t.Read(name);
   double delta = 1.E-3;
   for (size_t i = 0; i < t.xtest.size(); ++i) {
      EXPECT_NEAR(t.values1[i], t.values2[i], delta);
   }
}