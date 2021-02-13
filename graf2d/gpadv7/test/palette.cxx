#include "gtest/gtest.h"

#include "ROOT/RPalette.hxx"

// Test palette interpolation
TEST(Palette, DiscreteTrivial)
{
   using namespace ROOT::Experimental;
   RPalette p1(RPalette::kDiscrete, {RColor::kWhite, RColor::kBlack});
   EXPECT_EQ(p1.IsDiscrete(), true);
   EXPECT_EQ(p1.GetColor(0.), RColor::kWhite);
   EXPECT_EQ(p1.GetColor(1.), RColor::kBlack);
}

TEST(Palette, DiscretePlaced)
{
   using namespace ROOT::Experimental;
   RPalette p1(RPalette::kDiscrete, {{0., RColor::kWhite}, {.3, RColor::kRed}, {.7, RColor::kBlue}, {1., RColor::kBlack}});
   EXPECT_EQ(p1.IsDiscrete(), true);
   EXPECT_EQ(p1.GetColor(0.), RColor::kWhite);
   EXPECT_EQ(p1.GetColor(.3), RColor::kRed);
   EXPECT_EQ(p1.GetColor(.7), RColor::kBlue);
   EXPECT_EQ(p1.GetColor(1.), RColor::kBlack);
}

TEST(Palette, InterpolateSimple)
{
   using namespace ROOT::Experimental;
   RPalette p1({RColor::kRed, RColor::kBlue});

   EXPECT_EQ(p1.IsGradient(), true);
   EXPECT_EQ(p1.GetColor(0.), RColor::kRed);
   EXPECT_EQ(p1.GetColor(0.2), RColor(0xcc,0x00,0x33));
   EXPECT_EQ(p1.GetColor(0.4), RColor(0x99,0x00,0x66));
   EXPECT_EQ(p1.GetColor(0.5), RColor(0x80,0x00,0x80));
   EXPECT_EQ(p1.GetColor(0.6), RColor(0x66,0x00,0x99));
   EXPECT_EQ(p1.GetColor(0.8), RColor(0x33,0x00,0xcc));
   EXPECT_EQ(p1.GetColor(1.), RColor::kBlue);
}

