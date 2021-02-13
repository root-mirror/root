#include "gtest/gtest.h"

#include "ROOT/RHist.hxx"

// Test "0 + x = x"
TEST(HistAddTest, AddEmptyHist) {
  ROOT::Experimental::RH1F hFrom({100,0.,1});
  ROOT::Experimental::RH1F hTo({100,0.,1});
  hFrom.Fill({0.1111}, 0.12f);
  hFrom.Fill({0.1111}, 0.34f);
  ROOT::Experimental::Add(hTo, hFrom);
  EXPECT_EQ(2, hTo.GetEntries());
  EXPECT_FLOAT_EQ(0.46f, hTo.GetBinContent({0.1111}));
}

// Test "x + 0 = x"
TEST(HistAddTest, AddEmptySelf) {
   ROOT::Experimental::RH1F hFrom({100,0.,1});
   ROOT::Experimental::RH1F hTo({100,0.,1});
   hTo.Fill({0.1111}, 0.12f);
   hTo.Fill({0.1111}, 0.34f);
   ROOT::Experimental::Add(hTo, hFrom);
   EXPECT_EQ(2, hTo.GetEntries());
   EXPECT_FLOAT_EQ(0.46f, hTo.GetBinContent({0.1111}));
}

// Test "x + x = 2*x"
TEST(HistAddTest, AddSelfHist) {
   ROOT::Experimental::RH1F hist({100,0.,1});
   hist.Fill({0.1111}, .12f);
   hist.Fill({0.1111}, .34f);
   ROOT::Experimental::Add(hist, hist);
   EXPECT_EQ(4, hist.GetEntries());
   EXPECT_FLOAT_EQ(0.92f, hist.GetBinContent({0.1111}));
}

// Test "x - x = 0"
TEST(HistAddTest, SubstractSelfHist) {
   ROOT::Experimental::RH1F hFrom({100,0.,1});
   ROOT::Experimental::RH1F hTo({100,0.,1});
   hFrom.Fill({0.1111}, -.42f);
   hTo.Fill({0.1111}, .42f);
   ROOT::Experimental::Add(hTo, hFrom);
   EXPECT_EQ(2, hTo.GetEntries());
   EXPECT_FLOAT_EQ(.00f, hTo.GetBinContent({0.1111}));
}

// Test "x + y" with less STAT
TEST(HistAddTest, AddHist) {
   ROOT::Experimental::RH1F hFrom({100,0.,1});
   // Implicitly tests ROOT-8485.
   ROOT::Experimental::RHist<1, float> hTo({100,0.,1});
   hFrom.Fill({0.1111}, .42f);
   hTo.Fill({0.1111}, .17f);
   ROOT::Experimental::Add(hTo, hFrom);
   EXPECT_EQ(2, hTo.GetEntries());
   EXPECT_FLOAT_EQ(0.59f, hTo.GetBinContent({0.1111}));
}
