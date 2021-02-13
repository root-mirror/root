// @(#)root/test:$Id$
// Author: Lorenzo Moneta       06/2005
//         Martin Storø Nyfløtt 05/2017

#include "StatFunction.h"

#include "Math/Integrator.h"
#include "Math/Derivator.h"
#include "TBenchmark.h"

#include "gtest/gtest.h"

using namespace ROOT::Math;

void StatFunction::TestIntegral(IntegrationOneDim::Type algoType)
{
   // scan all values from fXMin to fXMax
   double dx = (fXMax - fXMin) / fNFuncTest;

   // create Integrator
   Integrator ig(algoType, 1.E-12, 1.E-12, 100000);
   ig.SetFunction(*this);

   for (int i = 0; i < fNFuncTest; ++i) {
      double v1 = fXMin + dx * i; // value used for testing
      double q1 = Cdf(v1);
      // calculate integral of pdf
      double q2 = 0;

      // lower integral (cdf)
      if (!fHasLowRange)
         q2 = ig.IntegralLow(v1);
      else
         q2 = ig.Integral(fXLow, v1);

      EXPECT_EQ(ig.Status(), 0);

      // use a larger scale (integral error is 10-9)
      double err = ig.Error();

      // Gauss integral sometimes returns an error of 0
      err = std::max(err, std::numeric_limits<double>::epsilon());
      double scale = std::max(fScaleIg * err / std::numeric_limits<double>::epsilon(), 1.);

      // test integral
      EXPECT_NEAR(q1, q2, scale * std::numeric_limits<double>::epsilon());
   }
}

void StatFunction::TestDerivative()
{
   // scan all values from fXMin to fXMax
   double dx = (fXMax - fXMin) / fNFuncTest;
   // create CDF function
   Functor1D func(this, &StatFunction::Cdf);
   Derivator d(func);
   for (int i = 0; i < fNFuncTest; ++i) {
      double v1 = fXMin + dx * i; // value used  for testing
      double q1 = Pdf(v1);
      // calculate derivative of cdf
      double q2 = 0;
      if (fHasLowRange && v1 == fXLow)
         q2 = d.EvalForward(v1);
      else if (fHasUpRange && v1 == fXUp)
         q2 = d.EvalBackward(v1);
      else
         q2 = d.Eval(v1);

      EXPECT_EQ(d.Status(), 0);

      double err = d.Error();
      double scale = std::max(1., fScaleDer * err / std::numeric_limits<double>::epsilon());

      EXPECT_NEAR(q1, q2, scale * std::numeric_limits<double>::epsilon());
   }
}

void StatFunction::TestInverse1(RootFinder::EType algoType)
{
   int maxitr = 2000;
   double abstol = 1.E-15;
   double reltol = 1.E-15;

   // scan all values from 0.05 to 0.95  to avoid problem at the border of definitions
   double x1 = 0.05;
   double x2 = 0.95;
   double dx = (x2 - x1) / fNFuncTest;
   double vmin = Quantile(dx / 2);
   double vmax = Quantile(1. - dx / 2);

   // test ROOT finder algorithm function without derivative
   RootFinder rf1(algoType);
   for (int i = 1; i < fNFuncTest; ++i) {
      double v1 = x1 + dx * i; // value used  for testing
      auto fInv = [this, v1](double x) { return this->Cdf(x) - v1; };

      Functor1D func(fInv);
      rf1.SetFunction(func, vmin, vmax);
      EXPECT_TRUE(rf1.Solve(maxitr, abstol, reltol));

      // test that quantile value correspond:
      double q1 = rf1.Root();
      double q2 = Quantile(v1);

      EXPECT_NEAR(q1, q2, fScaleInv * std::numeric_limits<double>::epsilon());
   }
}

void StatFunction::TestInverse2(RootFinder::EType algoType)
{
   int maxitr = 2000;
   // put lower tolerance
   double abstol = 1.E-12;
   double reltol = 1.E-12;

   // scan all values from 0.05 to 0.95  to avoid problem at the border of definitions
   double x1 = 0.05;
   double x2 = 0.95;
   double dx = (x2 - x1) / fNFuncTest;
   // starting root is always on the left to avoid to go negative
   // it is very sensible at the starting point
   double vstart = fStartRoot; // depends on function shape
   // test ROOT finder algorithm function with derivative
   RootFinder rf1(algoType);

   for (int i = 1; i < fNFuncTest; ++i) {
      double v1 = x1 + dx * i; // value used  for testing

      auto fInv = [this, v1](double x) { return this->Cdf(x) - v1; };
      // make a gradient function using inv function and derivative (which is pdf)
      GradFunctor1D func(fInv, *this);
      // use as estimate the quantile at 0.5
      rf1.SetFunction(func, vstart);
      EXPECT_TRUE(rf1.Solve(maxitr, abstol, reltol));

      // test that quantile value correspond:
      double q1 = rf1.Root();
      double q2 = Quantile(v1);

      EXPECT_NEAR(q1, q2, fScaleInv * std::numeric_limits<double>::epsilon());
   }
}
