#include "gtest/gtest.h"
#include "ROOT/RHist.hxx"
#include "ROOT/RHistBinIter.hxx"
#include <cmath>

using namespace ROOT::Experimental;

// Tests the number of bins
TEST(BinIterNBins, NumBins) {
  RH2F h({2, -10., 0.}, {3, -1., 10.});
  EXPECT_EQ(4*5, h.GetImpl()->GetNBins());

  int nBins = 0;
  for (auto bin: h) {
    (void)bin;
    ++nBins;
  }
  EXPECT_EQ(h.GetImpl()->GetNBinsNoOver(), nBins);
}

// Tests the bin ref: fill, retrieve.
TEST(BinIterNBins, BinRef) {
  RH2F h({10, -1., 1.}, {10, -1., 1.});

  double x = 0.11;
  double y = 0.33;
  double w = 3.141;
  h.Fill({x, y}, w);
  h.Fill({x, y}, 2. * w);
  EXPECT_FLOAT_EQ(3. * w, h.GetBinContent({x,y}));

  int xbin = h.GetImpl()->GetAxis(0).FindBin(x);
  EXPECT_EQ(6, xbin);
  int ybin = h.GetImpl()->GetAxis(1).FindBin(y);
  EXPECT_EQ(7, ybin);

  int linbin =  h.GetImpl()->GetBinIndex({x, y});
  const int expectedLinBin = xbin + (ybin - 1) * 10;
  EXPECT_EQ(expectedLinBin, linbin);

  EXPECT_FLOAT_EQ(3. * w, h.GetImpl()->GetBinContent(linbin));

  int foundbin = -1;
  double foundweight = -1.;
  double founduncert = -1.;
  RH2F::CoordArray_t foundcoord{};

  int nBins = 1;
  for (auto bin: h) {
    auto binCenter = bin.GetCenter();
    if (std::fabs(binCenter[0] - x) < 0.1 && std::fabs(binCenter[1] - y) < 0.1) {
      foundbin = nBins;
      foundcoord = binCenter;
      foundweight = bin.GetContent();
      founduncert = bin.GetUncertainty();
    }
    ++nBins;
  }

  EXPECT_FLOAT_EQ(0.1, foundcoord[0]);
  EXPECT_FLOAT_EQ(0.3, foundcoord[1]);
  EXPECT_EQ(expectedLinBin, foundbin);

  EXPECT_FLOAT_EQ(3. * w, foundweight);
  EXPECT_FLOAT_EQ(std::sqrt(w * w + (2 * w) * (2 * w)), founduncert);
}

