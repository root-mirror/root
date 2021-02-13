// @(#)root/test:$Id$
// Author: Lorenzo Moneta       06/2005
//         Martin Storø Nyfløtt 05/2017
///////////////////////////////////////////////////////////////////////////////////
//
//  MathMore test suite
//  ==============================
//
//  This program performs tests :
//     - numerical integration, derivation and root finders
//     - it compares for various values of the gamma and beta distribution)
//          - the numerical calculated integral of pdf  with cdf function,
//          - the calculated derivative of cdf with pdf
//          - the inverse (using root finder) of cdf with quantile
//
//     to run the program outside ROOT do:
//        > make stressMathMoreUnit
//        > ctest -R gtest-math-mathmore-test-stressMathMoreUnit
//

#include "StatFunction.h"

#include "Math/DistFunc.h"

#include "TBenchmark.h"

#include "gtest/gtest.h"

using ::testing::TestWithParam;
using ::testing::Values;

using namespace ROOT::Math;

template <typename T>
class MathMoreStress : public testing::Test {
private:
   T fStatFunctionFactory;

protected:
   template <typename J>
   void RunRangedTest(J ref)
   {
      for (int i = this->fStatFunctionFactory.GetRangeStart(); i < this->fStatFunctionFactory.GetRangeEnd(); i++) {
         StatFunction dist = fStatFunctionFactory.Produce(i);
         ref(dist);
      }
   }
};

class GammaTestFactory {
public:
   StatFunction Produce(const int step)
   {
      auto pdf = [](double x, double a, double b) { return gamma_pdf(x, a, b, 0.); };
      auto cdf = [](double x, double a, double b) { return gamma_cdf(x, a, b, 0.); };
      auto quantile = [](double x, double a, double b) { return gamma_quantile(x, a, b); };

      StatFunction dist(pdf, cdf, quantile, 0.);

      dist.SetNTest(10000);

      dist.SetTestRange(0., 10.);
      dist.SetScaleDer(10); // few tests fail here
      dist.SetScaleInv(10000);
      double k = std::pow(2., double(step - 1));
      double theta = 2. / double(step);
      dist.SetParameters(k, theta);

      if (k <= 1)
         dist.SetStartRoot(0.1);
      else
         dist.SetStartRoot(k * theta - 1.);

      std::string name = "Gamma(" + Util::ToString(int(k)) + "," + Util::ToString(theta) + ") ";
      std::cout << "\nTest " << name << " distribution\n";

      return dist;
   }

   int GetRangeStart() { return 1; }
   int GetRangeEnd() { return 5; }
};

class BetaTestFactory {
public:
   StatFunction Produce(const int step)
   {
      auto pdf = [](double x, double a, double b) { return beta_pdf(x, a, b); };
      auto cdf = [](double x, double a, double b) { return beta_cdf(x, a, b); };
      auto quantile = [](double x, double a, double b) { return beta_quantile(x, a, b); };

      StatFunction dist(pdf, cdf, quantile, 0., 1.);

      dist.SetNTest(10000);
      dist.SetTestRange(0., 1.);

      // avoid the case where alpha or beta = 1
      double alpha = step + 2;
      double beta = 6 - step;
      dist.SetParameters(alpha, beta);
      dist.SetStartRoot(alpha / (alpha + beta)); // use mean value

      std::string name = "Beta(" + Util::ToString(int(alpha)) + "," + Util::ToString(beta) + ") ";
      std::cout << "\nTest " << name << " distribution\n";

      return dist;
   }

   int GetRangeStart() { return 0; }
   int GetRangeEnd() { return 4; }
};

TYPED_TEST_SUITE_P(MathMoreStress);

TYPED_TEST_P(MathMoreStress, kADAPTIVESINGULAR)
{
   this->RunRangedTest([](StatFunction dist) { dist.TestIntegral(IntegrationOneDim::kADAPTIVESINGULAR); });
}

TYPED_TEST_P(MathMoreStress, kGAUSS)
{
   this->RunRangedTest([](StatFunction dist) {
      dist.SetScaleIg(100); // relax for Gauss integral
      dist.TestIntegral(IntegrationOneDim::kGAUSS);
   });
}

TYPED_TEST_P(MathMoreStress, TestDerivative)
{
   this->RunRangedTest([](StatFunction dist) { dist.TestDerivative(); });
}

TYPED_TEST_P(MathMoreStress, kBRENT)
{
   this->RunRangedTest([](StatFunction dist) { dist.TestInverse1(RootFinder::kBRENT); });
}

TYPED_TEST_P(MathMoreStress, kGSLBRENT)
{
   this->RunRangedTest([](StatFunction dist) { dist.TestInverse1(RootFinder::kGSL_BRENT); });
}

TYPED_TEST_P(MathMoreStress, kGSLSTEFFENSON)
{
   this->RunRangedTest([](StatFunction dist) { dist.TestInverse2(RootFinder::kGSL_STEFFENSON); });
}

REGISTER_TYPED_TEST_SUITE_P(MathMoreStress, kADAPTIVESINGULAR, kGAUSS, TestDerivative, kBRENT, kGSLBRENT,
                           kGSLSTEFFENSON);

typedef testing::Types<BetaTestFactory, GammaTestFactory> Factories_t;

INSTANTIATE_TYPED_TEST_SUITE_P(StressMathMore, MathMoreStress, Factories_t);
