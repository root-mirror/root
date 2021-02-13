// Tests for the RooNaNPacker
// Authors: Stephan Hageboeck, CERN  04/2020
#include "RooNaNPacker.h"
#include "RooRealVar.h"
#include "RooGenericPdf.h"
#include "RooMinimizer.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooRealSumPdf.h"
#include "RooRandom.h"

#include "gtest/gtest.h"

#include <cmath>
#include <bitset>
#include <cstdint>

void dumpFloats(double val) {
  float tmp[2];
  std::memcpy(tmp, &val, sizeof(double));

  for (int i = 1; i >= 0; --i) {
    unsigned long long ull = 0ull;
    std::memcpy(&ull, &tmp[i], 4);
    std::bitset<32> bits(ull);
    std::cout << bits << " " << std::flush;
  }
  std::cout << std::endl;
}

// Test if we can pack floats into NaNs, and recover them.
TEST(RooNaNPacker, CanPackStuffIntoNaNs)
{
  static_assert((RooNaNPacker::magicTag & RooNaNPacker::magicTagMask) == RooNaNPacker::magicTag, "Bit mask wrong.");
  constexpr bool dump = false;

  // Create a NaN that has 1.337f as payload
  RooNaNPacker rnp;
  rnp.setPayload(1.337f);
  EXPECT_TRUE(std::isnan(rnp._payload));
  EXPECT_TRUE(rnp.isNaNWithPayload());
  EXPECT_FLOAT_EQ(rnp.getPayload(), 1.337f);
  if (dump) dumpFloats(rnp._payload);

  // Create a normal double
  rnp._payload = 1.337;
  EXPECT_FALSE(std::isnan(rnp._payload));
  EXPECT_FALSE(rnp.isNaNWithPayload());
  EXPECT_DOUBLE_EQ(rnp._payload, 1.337);
  EXPECT_FLOAT_EQ(rnp.getPayload(), 0.f);
  if (dump) dumpFloats(rnp._payload);

  // Create a normal double
  rnp._payload = 4.;
  EXPECT_FALSE(std::isnan(rnp._payload));
  EXPECT_FALSE(rnp.isNaNWithPayload());
  EXPECT_DOUBLE_EQ(rnp._payload, 4.);
  EXPECT_FLOAT_EQ(rnp.getPayload(), 0.f);
  if (dump) dumpFloats(rnp._payload);

  // Create a simple payload
  rnp.setPayload(0.);
  EXPECT_TRUE(std::isnan(rnp._payload));
  EXPECT_TRUE(rnp.isNaNWithPayload());
  EXPECT_FLOAT_EQ(rnp.getPayload(), 0.f);
  if (dump) dumpFloats(rnp._payload);

  // Create a simple payload
  rnp.setPayload(2.);
  EXPECT_TRUE(std::isnan(rnp._payload));
  EXPECT_TRUE(rnp.isNaNWithPayload());
  EXPECT_FLOAT_EQ(rnp.getPayload(), 2.f);
  if (dump) dumpFloats(rnp._payload);

  // Create a simple payload
  rnp.setPayload(4.);
  EXPECT_TRUE(std::isnan(rnp._payload));
  EXPECT_TRUE(rnp.isNaNWithPayload());
  EXPECT_FLOAT_EQ(rnp.getPayload(), 4.f);
  if (dump) dumpFloats(rnp._payload);

  // Create a NaN that doesn't have the magic tag,
  // so no information encoded
  rnp._payload = std::numeric_limits<double>::quiet_NaN();
  const float tmp = 1234.5f;
  std::memcpy(&rnp._payload, &tmp, sizeof(float));
  EXPECT_TRUE(std::isnan(rnp._payload));
  EXPECT_FALSE(rnp.isNaNWithPayload());
  EXPECT_FLOAT_EQ(rnp.getPayload(), 0.f);
  if (dump) dumpFloats(rnp._payload);

}

/// Fit a simple linear function, that starts in the negative.
TEST(RooNaNPacker, FitSimpleLinear) {
  RooRealVar x("x", "x", -10, 10);
  RooRealVar a1("a1", "a1", 12., -5., 15.);
  RooGenericPdf pdf("pdf", "a1 + x", RooArgSet(x, a1));
  std::unique_ptr<RooDataSet> data(pdf.generate(x, 1000));
  std::unique_ptr<RooAbsReal> nll(pdf.createNLL(*data));

  ASSERT_FALSE(std::isnan(pdf.getVal(RooArgSet(x))));
  a1.setVal(-9.);
  ASSERT_TRUE(std::isnan(pdf.getVal(RooArgSet(x))));

  RooMinimizer minim(*nll);
  minim.setPrintLevel(-1);
  minim.setPrintEvalErrors(-1);
  minim.migrad();
  minim.hesse();
  auto fitResult = minim.save();

  EXPECT_EQ(fitResult->status(), 0);
  EXPECT_NEAR(a1.getVal(), 12., a1.getError());
}


/// Fit a parabola, where parameters are set up such that negative function values are obtained.
/// The minimiser needs to recover from that.
/// Test also that when recovery with NaN packing is switched off, the minimiser fails to recover.
TEST(RooNaNPacker, FitParabola) {
  constexpr bool verbose = false;

  RooRealVar x("x", "x", -10, 10);
  RooRealVar a1("a1", "a1", 12., -10., 20.);
  RooRealVar a2("a2", "a2", 1.1, -10., 20.);
  RooGenericPdf pdf("pdf", "a1 + x + a2 *x*x", RooArgSet(x, a1, a2));
  std::unique_ptr<RooDataSet> data(pdf.generate(x, 10000));
  auto nll = pdf.createNLL(*data);

  RooArgSet params(a1, a2);
  RooArgSet paramsInit;
  params.snapshot(paramsInit);
  RooArgSet evilValues;
  a1.setVal(-9.);
  a2.setVal(-1.);
  params.snapshot(evilValues);

  RooFitResult *fitResult1 = nullptr, *fitResult2 = nullptr;
  for (auto tryRecover : std::initializer_list<double>{0., 10.}) {
    params = evilValues;

    RooMinimizer::cleanup();
    RooMinimizer minim(*nll);
    minim.setRecoverFromNaNStrength(tryRecover);
    minim.setPrintLevel(-1);
    minim.setPrintEvalErrors(-1);
    minim.migrad();
    minim.hesse();
    minim.minos();
    auto fitResult = minim.save();
    (tryRecover != 0. ? fitResult1 : fitResult2) = fitResult;

    std::string config = (tryRecover != 0. ? "With recovery" : "Without recovery");
    const auto& a1Final = static_cast<RooRealVar&>(fitResult->floatParsFinal()[0]);
    const auto& a2Final = static_cast<RooRealVar&>(fitResult->floatParsFinal()[1]);

    if (tryRecover != 0.) {
      EXPECT_EQ(fitResult->status(), 0) << config;
      EXPECT_NEAR(a1Final.getVal(), static_cast<RooAbsReal&>(paramsInit["a1"]).getVal(), a1Final.getError()) << config;
      EXPECT_NEAR(a2Final.getVal(), static_cast<RooAbsReal&>(paramsInit["a2"]).getVal(), a2Final.getError()) << config;
    } else {
      EXPECT_LT(a1Final.getVal(), 0.);
      EXPECT_LT(a2Final.getVal(), 0.);
    }

    if (verbose) {
      std::cout << "Recovery strength:" << tryRecover << std::endl;
      fitResult->Print();
    }
  }

  // This makes clang-tidy happy:
  ASSERT_NE(fitResult1, nullptr);
  ASSERT_NE(fitResult2, nullptr);
  if (fitResult1 && fitResult2) { // makes clang-tidy happy
    EXPECT_LT(fitResult1->numInvalidNLL(), fitResult2->numInvalidNLL());
  }
}

/// Make coefficients of RooAddPdf sum to more than 1. Fitter should recover from this.
TEST(RooNaNPacker, FitAddPdf_DegenerateCoeff) {
  constexpr bool verbose = false;
  RooRandom::randomGenerator()->SetSeed(100);

  RooRealVar x("x", "x", 0., 10);
  RooRealVar a1("a1", "a1", 0.4, -10., 10.);
  RooRealVar a2("a2", "a2", 0.4, -10., 10.);
  RooGenericPdf pdf1("gen1", "exp(-2.*x)", RooArgSet(x));
  RooGenericPdf pdf2("gen2", "TMath::Gaus(x, 3, 2)", RooArgSet(x));
  RooGenericPdf pdf3("gen3", "x*x*x+1", RooArgSet(x));
  RooAddPdf pdf("sum", "a1*gen1 + a2*gen2 + (1-a1-a2)*gen3", RooArgList(pdf1, pdf2, pdf3), RooArgList(a1, a2));
  std::unique_ptr<RooDataSet> data(pdf.generate(x, 2000));
  auto nll = pdf.createNLL(*data);

  RooArgSet params(a1, a2);
  RooArgSet paramsInit;
  params.snapshot(paramsInit);

  RooArgSet evilValues;
  a1.setVal(0.6);
  a2.setVal(0.7);
  params.snapshot(evilValues);

  params = evilValues;

  RooFitResult *fitResult1 = nullptr, *fitResult2 = nullptr;
  for (auto tryRecover : std::initializer_list<double>{0., 10.}) {
    params = evilValues;

    RooMinimizer::cleanup();
    RooMinimizer minim(*nll);
    minim.setRecoverFromNaNStrength(tryRecover);
    minim.setPrintLevel(-1);
    minim.setPrintEvalErrors(-1);
    minim.migrad();
    minim.hesse();
    auto fitResult = minim.save();
    (tryRecover != 0. ? fitResult1 : fitResult2) = fitResult;

    const auto& a1Final = static_cast<RooRealVar&>(fitResult->floatParsFinal()[0]);
    const auto& a2Final = static_cast<RooRealVar&>(fitResult->floatParsFinal()[1]);

    if (tryRecover != 0.) {
      EXPECT_EQ(fitResult->status(), 0) << "Recovery strength=" << tryRecover;
      EXPECT_NEAR(a1Final.getVal(), static_cast<RooAbsReal&>(paramsInit["a1"]).getVal(), a1Final.getError()) << "Recovery strength=" << tryRecover;
      EXPECT_NEAR(a2Final.getVal(), static_cast<RooAbsReal&>(paramsInit["a2"]).getVal(), a2Final.getError()) << "Recovery strength=" << tryRecover;
      EXPECT_NEAR(a1Final.getVal() + a2Final.getVal(), 0.8, 0.02) << "Check that coefficients sum to 1. " << "Recovery strength=" << tryRecover;
    } else {
      EXPECT_TRUE(a1Final.getVal() < 0. || a1Final.getVal() > 1. || a2Final.getVal() < 0. || a2Final.getVal() > 1.) << "Recovery strength=" << tryRecover;
    }

    if (verbose) {
      std::cout << "Recovery strength:" << tryRecover << "\n";
      fitResult->Print();
    }
  }

  // This makes clang-tidy happy:
  ASSERT_NE(fitResult1, nullptr);
  ASSERT_NE(fitResult2, nullptr);
  if (fitResult1 && fitResult2) { // makes clang-tidy happy
    EXPECT_LT(fitResult1->numInvalidNLL(), fitResult2->numInvalidNLL());
  }
}

/// Make coefficients of RooRealSumPdf sum to more than 1. Fitter should recover from this.
TEST(RooNaNPacker, Interface_RooAbsPdf_fitTo_RooRealSumPdf_DegenerateCoeff) {
  constexpr bool verbose = false;
  RooRandom::randomGenerator()->SetSeed(100);

  RooRealVar x("x", "x", 0., 10);
  RooRealVar a1("a1", "a1", 0.3, -10., 10.);
  RooRealVar a2("a2", "a2", 0.4, -10., 10.);
  RooGenericPdf pdf1("gen1", "exp(-0.5*x)", RooArgSet(x));
  RooGenericPdf pdf2("gen2", "TMath::Gaus(x, 5, 0.7)", RooArgSet(x));
  RooGenericPdf pdf3("gen3", "TMath::Gaus(x, 8, 0.8)", RooArgSet(x));
  RooRealSumPdf pdf("sum", "a1*gen1 + a2*gen2 + (1-a1-a2)*gen3", RooArgList(pdf1, pdf2, pdf3), RooArgList(a1, a2));
  std::unique_ptr<RooDataSet> data(pdf.generate(x, 5000));

  RooArgSet params(a1, a2);
  RooArgSet paramsInit;
  params.snapshot(paramsInit);

  RooArgSet evilValues;
  a1.setVal(0.6);
  a2.setVal(0.7);
  params.snapshot(evilValues);

  params = evilValues;

  RooFitResult *fitResult1 = nullptr, *fitResult2 = nullptr;
  for (auto tryRecover : std::initializer_list<double>{0., 10.}) {
    params = evilValues;

    auto fitResult = pdf.fitTo(*data, RooFit::PrintLevel(-1), RooFit::PrintEvalErrors(-1), RooFit::Save(), RooFit::RecoverFromUndefinedRegions(tryRecover));
    (tryRecover != 0. ? fitResult1 : fitResult2) = fitResult;

    const auto& a1Final = static_cast<RooRealVar&>(fitResult->floatParsFinal()[0]);
    const auto& a2Final = static_cast<RooRealVar&>(fitResult->floatParsFinal()[1]);

    if (tryRecover != 0.) {
      EXPECT_EQ(fitResult->status(), 0) << "Recovery strength=" << tryRecover;
      EXPECT_NEAR(a1Final.getVal(), static_cast<RooAbsReal&>(paramsInit["a1"]).getVal(), 3.*a1Final.getError()) << "Recovery strength=" << tryRecover;
      EXPECT_NEAR(a2Final.getVal(), static_cast<RooAbsReal&>(paramsInit["a2"]).getVal(), 3.*a2Final.getError()) << "Recovery strength=" << tryRecover;
      EXPECT_GE(a1Final.getVal() + a2Final.getVal(), 0.) << "Check that coefficients are in [0, 1]. " << "Recovery strength=" << tryRecover;
      EXPECT_LE(a1Final.getVal() + a2Final.getVal(), 1.) << "Check that coefficients sum to [0, 1]. " << "Recovery strength=" << tryRecover;
    } else {
      EXPECT_TRUE(a1Final.getVal() < 0. || a1Final.getVal() > 1. || a2Final.getVal() < 0. || a2Final.getVal() > 1.) << "Recovery strength=" << tryRecover;
    }

    if (verbose) {
      std::cout << "Recovery strength:" << tryRecover << "\n";
      fitResult->Print();
    }
  }

  // This makes clang-tidy happy:
  ASSERT_NE(fitResult1, nullptr);
  ASSERT_NE(fitResult2, nullptr);

  if (fitResult1 && fitResult2) { // makes clang-tidy happy
    EXPECT_LT(fitResult1->numInvalidNLL(), fitResult2->numInvalidNLL());
  }
}

