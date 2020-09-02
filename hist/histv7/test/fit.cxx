#include "gtest/gtest.h"

#include "ROOT/RHist.hxx"

#include "TRandom3.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TVirtualFitter.h"
#include "Fit/Fitter.h"
#include "Math/WrappedMultiTF1.h"

using namespace ROOT::Experimental;

TEST(FittingTest, RH1D) {
   // Standard Chi2
   auto test = [](ROOT::Experimental::RH1D & hist) {
      // Set gauss function
      std::string fname("gaus");
      TF1 * f1 = (TF1*)gROOT->GetFunction(fname.c_str());
      f1->SetParameter(0,10);
      f1->SetParameter(1,0);
      f1->SetParameter(2,3.0);

      TRandom3 rndm;
      double chi2ref = 0;

      // Fill a histogram
      for (int i = 0; i <1000; ++i)
         hist.Fill(rndm.Gaus(0,1));
      
      int ndim = hist.GetNDim();
      
      // Fill a BinData structure with content of histogram
      ROOT::Fit::BinData fitData;
      ROOT::Fit::DataOptions fitOption;
      ROOT::Fit::FitConfig fitConfig;
      ROOT::Fit::DataRange range(ndim);
      if (fitOption.fUseRange) {
         RFit::GetFunctionRange(*f1,range);
      }

      RFit::BinContentToBinData(hist, fitData, f1, fitOption, range);

      // Create gauss function for fitting
      ROOT::Math::WrappedMultiTF1 wF1(*f1);
      // Create IParamMultiFunction from TF1 to avoid gradient calculation
      ROOT::Math::IParamMultiFunction & iParamFunc = wF1;

      // Set parameters of IParamMultiFunction
      double p[3] = {100,0,3.};
      iParamFunc.SetParameters(p);

      // Create the fitter
      ROOT::Fit::Fitter fitter;

      bool ret = fitter.Fit(fitData, iParamFunc);
      if (ret)
         fitter.Result().Print(std::cout);
      else {
         std::cout << "Chi2 Fit Failed " << std::endl;
         return -1;
      }
      chi2ref = fitter.Result().Chi2();

      // Compare with TH1::Fit
      TVirtualFitter::SetDefaultFitter("Minuit2");
      std::cout << "\n******************************\n\t RHist::Fit Result \n" << std::endl;
      f1->SetParameters(p);
      auto res = ROOT::Experimental::Fit(hist, f1, fitOption, fitConfig);

      // iret |= compareResult( res->Chi2(), chi2ref,"1D histogram chi2 fit");

      // RAxisEquidistant axis = Internal::AxisConfigToType<RAxisConfig::kEquidistant>()(cfg);
      // EXPECT_EQ(axis.GetTitle(), title);
      // EXPECT_EQ(axis.GetNBinsNoOver(), 10);
      // EXPECT_EQ(axis.GetMinimum(), 1.2);
      // EXPECT_DOUBLE_EQ(axis.GetMaximum(), 3.4);
   };

   {
      SCOPED_TRACE("RH1D with standard Chi2 fit");
      ROOT::Experimental::RH1D hist1d({30, -5., 5.});
      test(hist1d);
   }
}
