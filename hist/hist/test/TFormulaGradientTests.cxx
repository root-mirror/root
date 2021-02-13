/// \file TFormulaGradientTests.cxx
///
/// \brief The file contain unit tests which test the clad-based gradient
///        computations.
///
/// \author Vassil Vassilev <vvasilev@cern.ch>
///
/// \date Oct, 2018
///
/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOTUnitTestSupport.h"

#include <Math/MinimizerOptions.h>
#include <TFormula.h>
#include <TF1.h>
#include <TFitResult.h>

TEST(TFormulaGradientPar, Sanity)
{
   TFormula f("f", "x*std::sin([0]) - y*std::cos([1])");
   double p[] = {30, 60};
   f.SetParameters(p);
   double x[] = {1, 2};
   TFormula::GradientStorage result(2);
   f.GradientPar(x, result);

   ASSERT_FLOAT_EQ(x[0] * std::cos(30), result[0]);
   ASSERT_FLOAT_EQ(-x[1] * -std::sin(60), result[1]);
}

TEST(TFormulaGradientPar, ResultUpsize)
{
   TFormula f("f", "std::sin([1]) - std::cos([0])");
   double p[] = {60, 30};
   f.SetParameters(p);
   TFormula::GradientStorage result;
   double x[] = {2, 1};

   ASSERT_TRUE(0 == result.size());
   ROOT_EXPECT_WARNING(f.GradientPar(x, result),
                       "TFormula::GradientPar",
                       "The size of gradient result is 0 but 2 is required. Resizing."
                       );

   ASSERT_FLOAT_EQ(std::cos(30), result[1]);
   ASSERT_FLOAT_EQ(std::sin(60), result[0]);
   ASSERT_TRUE(2 == result.size());
}

TEST(TFormulaGradientPar, ResultDownsize)
{
   TFormula f("f", "std::sin([0])");
   double p[] = {60};
   f.SetParameters(p);
   TFormula::GradientStorage result(2);
   double x[] = {1};

   ASSERT_TRUE(2 == result.size());

   ROOT_EXPECT_NODIAG(f.GradientPar(x, result));

   ASSERT_FLOAT_EQ(std::cos(60), result[0]);
   ASSERT_TRUE(2 == result.size());
}

TEST(TFormulaGradientPar, GausCrossCheck)
{
   auto h = new TF1("f1", "gaus");
   // auto h = new TF1("f1", "landau"); -- inheritently does not work. See DIFLAN
   //crystalball, breitwigner, cheb3, bigaus,
   //auto h = new TF1("f1", "");
   double p[] = {3, 1, 2};
   h->SetParameters(p);
   double x[] = {0};
   TFormula::GradientStorage result_clad(3);
   h->GetFormula()->GradientPar(x, result_clad);

   TFormula::GradientStorage result_num(3);
   h->GradientPar(x, result_num.data());

   ASSERT_FLOAT_EQ(result_num[0], result_clad[0]);
   ASSERT_FLOAT_EQ(result_num[1], result_clad[1]);
   ASSERT_FLOAT_EQ(result_num[2], result_clad[2]);
}

TEST(TFormulaGradientPar, BreitWignerCrossCheck)
{
   auto h = new TF1("f1", "breitwigner");
   double p[] = {3, 1, 2.1};
   h->SetParameters(p);
   double x[] = {0};
   TFormula::GradientStorage result_clad(3);
   TFormula* formula = h->GetFormula();
   formula->GradientPar(x, result_clad);
   TFormula::GradientStorage result_num(3);
   h->GradientPar(x, result_num.data());

   ASSERT_FLOAT_EQ(result_num[0], result_clad[0]);
   ASSERT_FLOAT_EQ(result_num[1], result_clad[1]);
   ASSERT_FLOAT_EQ(result_num[2], result_clad[2]);
}

TEST(TFormulaGradientPar, BreitWignerCrossCheckAccuracyDemo)
{
   auto h = new TF1("f1", "breitwigner");
   double p[] = {3, 1, 2};
   h->SetParameters(p);
   double x[] = {0};
   TFormula::GradientStorage result_clad(3);
   TFormula* formula = h->GetFormula();
   formula->GradientPar(x, result_clad);
   TFormula::GradientStorage result_num(3);
   h->GradientPar(x, result_num.data());

   // This is a classical example why clad is better.
   // The gradient with respect to gamma leads to a cancellation when gamma is
   // set to 2. This is not a problem for clad yielding the correct result of 0
   ASSERT_FLOAT_EQ(0, result_clad[2]);

   // However, that is not the case for the numerical approach where we give
   // a small but non-zero result.
   EXPECT_NEAR(0, result_num[2], /*abs_error*/1e-13);
}

// FIXME: Add more: crystalball, cheb3, bigaus?

// FIXME: Disable because of a known failure in -Druntime_cxxmodules=On.
// This is identical to getFuncBody test failure in roottest.
TEST(DISABLED_TFormulaGradientPar, GetGradFormula)
{
   TFormula f("f", "gaus");
   double p[] = {3, 1, 2};
   f.SetParameters(p);
   ASSERT_TRUE(f.GenerateGradientPar());
   std::string s = f.GetGradientFormula().Data();
   ASSERT_THAT(s, testing::ContainsRegex("void TFormula____id[0-9]*_grad"));
}

