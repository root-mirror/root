// Tests for the RooRealVar
// Authors: Stephan Hageboeck, CERN  07/2020

#include "RooRealVar.h"
#include "RooUniformBinning.h"
#include "RooBinning.h"

#include <TFile.h>

#include "gtest/gtest.h"

/// ROOT-10781
/// Searching binning in linked lists is slow, so these were replaced by unordered maps.
/// Here, we test that sharing alternative binning still works.
TEST(RooRealVar, AlternativeBinnings)
{ RooArgSet survivingX;
  {
    RooRealVar x("x", "x", -10, 10);
    x.setBinning(RooUniformBinning(-9, 10, 20, "uniform"), "uniform");

    RooArgSet xSet(x);
    auto clones = xSet.snapshot(true);

    auto& uniform = dynamic_cast<RooRealVar&>((*clones)["x"]).getBinning("uniform");
    EXPECT_TRUE(uniform.isUniform());
    EXPECT_EQ(uniform.lowBound(), -9.);
    EXPECT_EQ(uniform.highBound(), 10.);
    uniform.setRange(-10., 10.);

    double boundaries[] = {-5., 5., 10.};
    x.setBinning(RooBinning(2, boundaries, "custom"), "uniform");

    auto& overwrittenBinning = dynamic_cast<RooRealVar&>((*clones)["x"]).getBinning("uniform");
    EXPECT_EQ(overwrittenBinning.lowBound(), -5.);
    EXPECT_EQ(overwrittenBinning.highBound(), 10.);
    EXPECT_EQ(overwrittenBinning.binWidth(0), 10.);

    delete clones;

    auto& uniformFromX = x.getBinning("uniform");
    EXPECT_EQ(&uniformFromX, &overwrittenBinning);
    EXPECT_EQ(uniformFromX.lowBound(), -5.);

    RooArgSet(x).snapshot(survivingX);
  }

  auto& transferredOwnership = dynamic_cast<RooRealVar&>(survivingX["x"]).getBinning("uniform");
  EXPECT_EQ(transferredOwnership.numBins(), 2);
}

/// ROOT-10781
/// Searching binning in linked lists is slow, so these were replaced by unordered maps.
/// Here, we test that sharing works also after writing to files.
TEST(RooRealVar, AlternativeBinningsPersistency) {
  RooArgSet localCopy;
  {
    TFile infile("testRooRealVar_data1.root");
    ASSERT_TRUE(infile.IsOpen());

    RooArgSet* fromFile = nullptr;
    infile.GetObject("SetWithX", fromFile);
    ASSERT_NE(fromFile, nullptr);

    fromFile->snapshot(localCopy);
  }

  RooRealVar& x = dynamic_cast<RooRealVar&>(localCopy["x"]);
  ASSERT_TRUE(x.hasBinning("sharedBinning"));

  auto& binningFromFile = x.getBinning("sharedBinning");
  EXPECT_EQ(binningFromFile.lowBound(), -5.);
  EXPECT_EQ(binningFromFile.highBound(), 10.);
  EXPECT_EQ(binningFromFile.binWidth(0), 10.);

}
