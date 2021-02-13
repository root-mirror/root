// Test for RooBernstein
// Authors: Rahul Balasubramanian, CERN  05/2020

#include "RooRealVar.h"
#include "RooNumIntConfig.h"
#include "RooBernstein.h"
#include "TMath.h"
#include "RooPlot.h"
#include "TAxis.h"

#include "gtest/gtest.h"

using namespace RooFit;

void IntegrationChecker(double a0, double a1, double a2, double a3)
{
  RooRealVar x("x", "x", 0., 100.);

  // Set ranges for the variable
  x.setRange("range1", 0., 30.);
  x.setRange("range2", 30., 70.);
  x.setRange("range3", 70., 100.);
  x.setRange("FULL", 0., 100.);

  // Define coefficients for a bernstein polynomial of order 3
  RooRealVar c0("c0", "c1 coeff", a0, 0., 10.);
  RooRealVar c1("c1", "c1 coeff", a1, 0., 10.);
  RooRealVar c2("c2", "c2 coeff", a2, 0., 10.);
  RooRealVar c3("c3", "c3 coeff", a3, 0., 10.);

  // Build bernstein p.d.f in terms of coefficients
  RooBernstein bern("bern", "bernstein PDF", x, RooArgList(c0, c1, c2, c3));
  RooBernstein bernNumInt(bern);
  // Set normalization range
  bern.selectNormalizationRange("FULL",kTRUE);
  bernNumInt.selectNormalizationRange("FULL",kTRUE);

  RooNumIntConfig intConfig(*RooAbsReal::defaultIntegratorConfig());
  intConfig.setEpsAbs(1.E-15);
  intConfig.setEpsRel(1.E-12);

  intConfig.getConfigSection("RooIntegrator1D").setRealValue("maxSteps", 100);
  bernNumInt.setIntegratorConfig(intConfig);
  bernNumInt.forceNumInt(true);

  // Test ranged integration
  auto int_range1 = bern.createIntegral(x, "range1");
  auto int_range2 = bern.createIntegral(x, "range2");
  auto int_range3 = bern.createIntegral(x, "range3");
  auto int_full = bern.createIntegral(x, "FULL");

  auto numInt_range1 = bern.createIntegral(x, "range1");
  auto numInt_range2 = bern.createIntegral(x, "range2");
  auto numInt_range3 = bern.createIntegral(x, "range3");
  auto numInt_full = bern.createIntegral(x, "FULL");

  // closure
  EXPECT_LT(fabs(int_full->getVal() - int_range1->getVal() - int_range2->getVal() - int_range3->getVal()), 1e-10);
  EXPECT_LT(fabs(numInt_full->getVal() - numInt_range1->getVal() - numInt_range2->getVal() - numInt_range3->getVal()), 1e-10);

  // comparision with polynomial
  double accAnaVsNum = 1.;
  EXPECT_NEAR(int_range1->getVal(), numInt_range1->getVal(),
      accAnaVsNum/100.*numInt_range1->getVal())
  << "Analytical vs numerical integral"
  << " within " << accAnaVsNum << "%. ";
  
  EXPECT_NEAR(int_range1->getVal(), numInt_range1->getVal(),
      accAnaVsNum/100.*numInt_range2->getVal())
  << "Analytical vs numerical integral"
  << " within " << accAnaVsNum << "%. ";

  EXPECT_NEAR(int_range1->getVal(), numInt_range1->getVal(),
      accAnaVsNum/100.*numInt_range3->getVal())
  << "Analytical vs numerical integral"
  << " within " << accAnaVsNum << "%. ";

  EXPECT_NEAR(int_full->getVal(), numInt_full->getVal(),
      accAnaVsNum/100.*numInt_full->getVal())
  << "Analytical vs numerical integral"
  << " within " << accAnaVsNum << "%. ";
}

TEST(RooBernstein, NumInt)
{
  IntegrationChecker(0.3, 0.7, 0.2, 0.5);
}
