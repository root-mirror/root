// Tests for the RooProdPdf
// Authors: Stephan Hageboeck, CERN  02/2019

#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooGenericPdf.h"
#include "RooDataSet.h"

#include "gtest/gtest.h"

class TestProdPdf : public ::testing::Test {
protected:
  TestProdPdf() :
  Test()
  {
    datap.reset(prod.generate(RooArgSet(x), 1000));
    a.setConstant(true);
  }

  ~TestProdPdf() {

  }

  constexpr static double bTruth = -0.5;

  RooRealVar x{"x", "x", 2., 0., 5.};
  RooRealVar a{"a", "a", -0.2, -5., 0.};
  RooRealVar b{"b", "b", bTruth, -5., 0.};

  RooGenericPdf c1{"c1", "exp(x*a)", RooArgSet(x,a)};
  RooGenericPdf c2{"c2", "exp(x*b)", RooArgSet(x,b)};
  RooProdPdf prod{"mypdf", "mypdf", RooArgList(c1,c2)};
  std::unique_ptr<RooDataSet> datap{nullptr};
};

TEST_F(TestProdPdf, CachingOpt2) {
  prod.fitTo(*datap, RooFit::Optimize(2), RooFit::PrintLevel(-1));
  EXPECT_LT(fabs(b.getVal() - bTruth), b.getError()*1.1) << "b=" << b.getVal()
          << " +- " << b.getError() << " doesn't match truth value with O2.";
}

TEST_F(TestProdPdf, CachingOpt1) {
  prod.fitTo(*datap, RooFit::Optimize(1), RooFit::PrintLevel(-1));
  EXPECT_LT(fabs(b.getVal() - bTruth), b.getError()*1.1) << "b=" << b.getVal()
        << " +- " << b.getError() << " doesn't match truth value with O1.";
}

TEST_F(TestProdPdf, CachingOpt0) {
  prod.fitTo(*datap, RooFit::Optimize(0), RooFit::PrintLevel(-1));
  EXPECT_LT(fabs(b.getVal() - bTruth), b.getError()*1.1) << "b=" << b.getVal()
        << " +- " << b.getError() << " doesn't match truth value with O0.";
}
