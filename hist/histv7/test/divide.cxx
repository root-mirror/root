#include "gtest/gtest.h"

#include "ROOT/RHist.hxx"

// Test "0 / x = 0"
TEST(HistDivideTest, DivideEmptyHist) {
  ROOT::Experimental::RH1F hFrom({100,0.,1});
  ROOT::Experimental::RH1F hTo({100,0.,1});
  hFrom.Fill({0.1111}, 0.12f);
  hFrom.Fill({0.1111}, 0.34f);
  ROOT::Experimental::Divide(hTo, hFrom);
  EXPECT_EQ(2, hTo.GetEntries());
  EXPECT_FLOAT_EQ(.00f, hTo.GetBinContent({0.1111}));
}

// Test "x / +0 = +∞"
TEST(HistDivideTest, DivideEmptySelf) {
   ROOT::Experimental::RH1F hFrom({100,0.,1});
   ROOT::Experimental::RH1F hTo({100,0.,1});
   hTo.Fill({0.1111}, 0.12f);
   hTo.Fill({0.1111}, 0.34f);
   ROOT::Experimental::Divide(hTo, hFrom);
   EXPECT_EQ(2, hTo.GetEntries());
   EXPECT_EQ(INFINITY, hTo.GetBinContent({0.1111}));
}

// Test "x / x = 1"
TEST(HistDivideTest, DivideSelfHist) {
   ROOT::Experimental::RH1F hFrom({100,0.,1});
   ROOT::Experimental::RH1F hTo({100,0.,1});
   hFrom.Fill({0.1111}, .42f);
   hTo.Fill({0.1111}, .42f);
   ROOT::Experimental::Divide(hTo, hFrom);
   EXPECT_EQ(2, hTo.GetEntries());
   EXPECT_FLOAT_EQ(1.0f, hTo.GetBinContent({0.1111}));
}

// Test "x / y" with less STAT
TEST(HistDivideTest, DivideHist) {
   ROOT::Experimental::RH1F hFrom({100,0.,1});
   // Implicitly tests ROOT-8485.
   ROOT::Experimental::RHist<1, float> hTo({100,0.,1});
   hFrom.Fill({0.1111}, .42f);
   hTo.Fill({0.1111}, .17f);
   ROOT::Experimental::Divide(hTo, hFrom);
   EXPECT_EQ(2, hTo.GetEntries());
   EXPECT_FLOAT_EQ(.17f / .42f, hTo.GetBinContent({0.1111}));
}

// Test "2 * x / x = 2"
TEST(HistDivideTest, DivideMultipliedHistTo) {
   ROOT::Experimental::RH1F hFrom({100,0.,1});
   ROOT::Experimental::RH1F hTo({100,0.,1});
   hFrom.Fill({0.1111}, .42f);
   hTo.Fill({0.1111}, .42f);
   hTo *= 2;
   ROOT::Experimental::Divide(hTo, hFrom);
   EXPECT_EQ(2, hTo.GetEntries());
   EXPECT_FLOAT_EQ(2.0f, hTo.GetBinContent({0.1111}));
}

// Test "x / (2 * x) = 1/2"
TEST(HistDivideTest, DivideMultipliedHistFrom) {
   ROOT::Experimental::RH1F hFrom({100,0.,1});
   ROOT::Experimental::RH1F hTo({100,0.,1});
   hFrom.Fill({0.1111}, .42f);
   hTo.Fill({0.1111}, .42f);
   hFrom *= 2;
   ROOT::Experimental::Divide(hTo, hFrom);
   EXPECT_EQ(2, hTo.GetEntries());
   EXPECT_FLOAT_EQ(0.5f, hTo.GetBinContent({0.1111}));
}

// Test "(3 * x) / (2 * x) = 3/2"
TEST(HistDivideTest, DivideMultipliedHistToFrom) {
   ROOT::Experimental::RH1F hFrom({100,0.,1});
   ROOT::Experimental::RH1F hTo({100,0.,1});
   hFrom.Fill({0.1111}, .42f);
   hTo.Fill({0.1111}, .42f);
   hTo *= 3;
   hFrom *= 2;
   ROOT::Experimental::Divide(hTo, hFrom);
   EXPECT_EQ(2, hTo.GetEntries());
   EXPECT_FLOAT_EQ(3/2.0f, hTo.GetBinContent({0.1111}));
}
