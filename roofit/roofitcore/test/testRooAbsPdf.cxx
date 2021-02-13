// Tests for RooAbsPdf
// Author: Stephan Hageboeck, CERN 04/2020

#include "RooRealVar.h"
#include "RooGenericPdf.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"

#include "gtest/gtest.h"

#include <memory>

// ROOT-10668: Asympt. correct errors don't work when title and name differ
TEST(RooAbsPdf, AsymptoticallyCorrectErrors)
{
  RooRealVar x("x", "xxx", 0, 0, 10);
  RooRealVar a("a", "aaa", 2, 0, 10);
  // Cannot play with RooAbsPdf, since abstract.
  RooGenericPdf pdf("pdf", "std::pow(x,a)", RooArgSet(x, a));
  RooFormulaVar formula("w", "(x-5)*(x-5)*1.2", RooArgSet(x));

  std::unique_ptr<RooDataSet> data(pdf.generate(x, 5000));
  data->addColumn(formula);
  RooRealVar w("w", "weight", 1, 0, 20);
  RooDataSet weightedData("weightedData", "weightedData",
      RooArgSet(x, w), RooFit::Import(*data), RooFit::WeightVar(w));

  ASSERT_TRUE(weightedData.isWeighted());
  weightedData.get(0);
  ASSERT_NE(weightedData.weight(), 1);

  a = 1.2;
  auto result = pdf.fitTo(weightedData, RooFit::Save(), RooFit::AsymptoticError(true));
  const double aError = a.getError();
  a = 1.2;
  auto result2 = pdf.fitTo(weightedData, RooFit::Save());

  EXPECT_TRUE(result->isIdentical(*result2)) << "Fit results should be very similar.";
  EXPECT_GT(aError, a.getError()*2.) << "Asymptotically correct errors should be significantly larger.";
}

