#include "gtest/gtest.h"
#include "ROOT/RHistDrawable.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RColor.hxx"

using namespace ROOT::Experimental;

// Test drawing of histograms.
TEST(DrawTest, OneD)
{
   RAxisConfig xaxis{10, 0., 1.};
   auto h = std::make_shared<RH1D>(xaxis);
   RCanvas canv;
   canv.Draw(h);
   EXPECT_EQ(canv.NumPrimitives(), 2u);
}

TEST(DrawTest, TwoD)
{
   RAxisConfig xaxis{10, 0., 1.};
   RAxisConfig yaxis{{0., 1., 10., 100.}};
   auto h1 = std::make_shared<RH2I>(xaxis, yaxis);
   auto h2 = std::make_shared<RH2C>(xaxis, yaxis);
   RCanvas canv;
   canv.Draw(h1);
   canv.Draw(h2);
   // No RHist copt c'tor:
   // canv.Draw(RH2F(xaxis, yaxis));
   // canv.Draw(std::make_unique<RH2C>(xaxis, yaxis));
   EXPECT_EQ(canv.NumPrimitives(), 3u);
}

TEST(DrawTest, ThreeD)
{
   RAxisConfig xaxis{{0., 1.}};
   RAxisConfig yaxis{10, 0., 1.};
   RAxisConfig zaxis{{0., 1., 10., 100.}};
   auto h = std::make_shared<RH3F>(xaxis, yaxis, zaxis);
   RCanvas canv;
   canv.Draw(h);
   EXPECT_EQ(canv.NumPrimitives(), 2u);
}

// Drawing options:
TEST(DrawOptTest, OneD)
{
   RAxisConfig xaxis{10, 0., 1.};
   auto h = std::make_shared<RH1D>(xaxis);
   RCanvas canv;
   auto drawable = canv.Draw(h);
   drawable->AttrLine().SetColor(RColor::kRed);
   RColor shouldBeRed = drawable->AttrLine().GetColor();
   EXPECT_EQ(shouldBeRed, RColor::kRed);
}
