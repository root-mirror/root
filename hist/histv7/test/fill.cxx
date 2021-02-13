#include "gtest/gtest.h"

#include "ROOT/RHist.hxx"

// Test Fill(), FillN(), GetEntries(), GetBinContent(), GetBinUncertainty()




// RHist::Fill for 1D

// Test Fill() without weight and GetEntries()
TEST(HistFillTest, FillCoordEntries)
{
   ROOT::Experimental::RH1F hist({100, 0., 1});
   EXPECT_EQ(0, hist.GetEntries());
   hist.Fill({0.1111});
   EXPECT_EQ(1, hist.GetEntries());
   hist.Fill({0.1111});
   hist.Fill({0.2222});
   EXPECT_EQ(3, hist.GetEntries());
}

// Test Fill() without weight and GetBinContent()
TEST(HistFillTest, FillCoordContent)
{
   ROOT::Experimental::RH1F hist({100, 0., 1});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.1111}));
   hist.Fill({0.1111});
   EXPECT_FLOAT_EQ(1.f, hist.GetBinContent({0.1111}));
   hist.Fill({0.1111});
   hist.Fill({0.2222});
   EXPECT_FLOAT_EQ(2.f, hist.GetBinContent({0.1111}));
}

// Test Fill() without weight and GetBinUncertainty()
TEST(HistFillTest, FillCoordUncertainty)
{
   ROOT::Experimental::RH1F hist({100, 0., 1});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.1111}));
   hist.Fill({0.1111});
   EXPECT_FLOAT_EQ(1.f, hist.GetBinUncertainty({0.1111}));
   hist.Fill({0.1111});
   hist.Fill({0.2222});
   EXPECT_FLOAT_EQ(std::sqrt(2.f), hist.GetBinUncertainty({0.1111}));
}

// Test Fill() with weight and GetEntries()
TEST(HistFillTest, FillCoordWeightEntries)
{
   ROOT::Experimental::RH1F hist({100, 0., 1});
   EXPECT_EQ(0, hist.GetEntries());
   hist.Fill({0.1111}, .42f);
   EXPECT_EQ(1, hist.GetEntries());
   hist.Fill({0.1111}, .32f);
   hist.Fill({0.2222}, .52f);
   EXPECT_EQ(3, hist.GetEntries());
}

// Test Fill() with weight and GetBinContent()
TEST(HistFillTest, FillCoordWeightContent)
{
   ROOT::Experimental::RH1F hist({100, 0., 1});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.1111}));
   hist.Fill({0.1111}, .42f);
   EXPECT_FLOAT_EQ(.42f, hist.GetBinContent({0.1111}));
   hist.Fill({0.1111}, .32f);
   hist.Fill({0.2222}, .52f);
   EXPECT_FLOAT_EQ(.42f + .32f, hist.GetBinContent({0.1111}));
}

// Test Fill() with weight and GetBinUncertainty()
TEST(HistFillTest, FillCoordWeightUncertainty)
{
   ROOT::Experimental::RH1F hist({100, 0., 1});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.1111}));

   float weight1 = .42f;
   hist.Fill({0.1111}, weight1);
   EXPECT_FLOAT_EQ(std::sqrt(weight1 * weight1), hist.GetBinUncertainty({0.1111}));

   float weight2 = .32f;
   hist.Fill({0.1111}, weight2);
   hist.Fill({0.2222}, .52f);
   EXPECT_FLOAT_EQ(std::sqrt((weight1 * weight1) + (weight2 * weight2)), hist.GetBinUncertainty({0.1111}));
}

// Test Fill() with NaN and GetEntries()
TEST(HistFillTest, FillCoordNaNEntries)
{
   ROOT::Experimental::RH1F hist({100, 0., 1});
   EXPECT_EQ(0, hist.GetEntries());
   hist.Fill({0.1111}, NAN);
   EXPECT_EQ(1, hist.GetEntries());
   hist.Fill({0.1111}, .32f);
   hist.Fill({0.2222}, NAN);
   EXPECT_EQ(3, hist.GetEntries());
}

// Test Fill() with NaN and GetBinContent()
TEST(HistFillTest, FillCoordNaNContent)
{
   ROOT::Experimental::RH1F hist({100, 0., 1});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.1111}));
   hist.Fill({0.1111}, NAN);
   EXPECT_EQ(true, std::isnan(hist.GetBinContent({0.1111})));
   hist.Fill({0.1111}, .32f);
   hist.Fill({0.2222}, NAN);
   EXPECT_EQ(true, std::isnan(hist.GetBinContent({0.1111})));
}

// Test Fill() with NaN and GetBinUncertainty()
TEST(HistFillTest, FillCoordNaNUncertainty)
{
   ROOT::Experimental::RH1F hist({100, 0., 1});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.1111}));
   hist.Fill({0.1111}, NAN);
   EXPECT_EQ(true, std::isnan(hist.GetBinUncertainty({0.1111})));

   float weight1 = .42f;
   hist.Fill({0.1111}, weight1);
   EXPECT_EQ(true, std::isnan(hist.GetBinUncertainty({0.1111})));

   float weight2 = .32f;
   hist.Fill({0.1111}, weight2);
   hist.Fill({0.2222}, NAN);
   EXPECT_EQ(true, std::isnan(hist.GetBinUncertainty({0.1111})));
}

// Test Fill() with NaN and GetEntries()
TEST(HistFillTest, FillCoordInfEntries)
{
   ROOT::Experimental::RH1F hist({100, 0., 1});
   EXPECT_EQ(0, hist.GetEntries());
   hist.Fill({0.1111}, INFINITY);
   EXPECT_EQ(1, hist.GetEntries());
   hist.Fill({0.1111}, .32f);
   hist.Fill({0.2222}, INFINITY);
   EXPECT_EQ(3, hist.GetEntries());
}

// Test Fill() with NaN and GetBinContent()
TEST(HistFillTest, FillCoordInfContent)
{
   ROOT::Experimental::RH1F hist({100, 0., 1});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.1111}));
   hist.Fill({0.1111}, INFINITY);
   EXPECT_EQ(INFINITY, hist.GetBinContent({0.1111}));
   hist.Fill({0.1111}, .32f);
   hist.Fill({0.2222}, INFINITY);
   EXPECT_EQ(INFINITY, hist.GetBinContent({0.1111}));
   hist.Fill({0.2222}, -INFINITY);
   EXPECT_EQ(true, std::isnan(hist.GetBinContent({0.2222})));
}

// Test Fill() with NaN and GetBinUncertainty()
TEST(HistFillTest, FillCoordInfUncertainty)
{
   ROOT::Experimental::RH1F hist({100, 0., 1});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.1111}));
   hist.Fill({0.1111}, INFINITY);
   EXPECT_EQ(INFINITY, hist.GetBinUncertainty({0.1111}));

   float weight1 = .42f;
   hist.Fill({0.1111}, weight1);
   EXPECT_EQ(INFINITY, hist.GetBinUncertainty({0.1111}));

   float weight2 = .32f;
   hist.Fill({0.1111}, weight2);
   hist.Fill({0.2222}, INFINITY);
   EXPECT_EQ(INFINITY, hist.GetBinUncertainty({0.1111}));
   
   hist.Fill({0.2222}, -INFINITY);
   EXPECT_EQ(INFINITY, hist.GetBinUncertainty({0.2222}));
}

// Test FillN() without weights and GetEntries()
TEST(HistFillTest, FillNCoordsEntries)
{
   ROOT::Experimental::RH1F hist({100, 0., 1});
   EXPECT_EQ(0, hist.GetEntries());
   hist.FillN({{0.1111}, {0.2222}, {0.3333}});
   EXPECT_EQ(3, hist.GetEntries());
   hist.FillN({{0.1111}, {0.3333}});
   EXPECT_EQ(5, hist.GetEntries());
}

// Test FillN() without weights and GetBinContent()
TEST(HistFillTest, FillNCoordsContent)
{
   ROOT::Experimental::RH1F hist({100, 0., 1});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.1111}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.2222}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.3333}));

   hist.FillN({{0.1111}, {0.2222}, {0.3333}});
   EXPECT_FLOAT_EQ(1.f, hist.GetBinContent({0.1111}));
   EXPECT_FLOAT_EQ(1.f, hist.GetBinContent({0.2222}));
   EXPECT_FLOAT_EQ(1.f, hist.GetBinContent({0.3333}));

   hist.FillN({{0.1111}, {0.3333}});
   EXPECT_FLOAT_EQ(2.f, hist.GetBinContent({0.1111}));
   EXPECT_FLOAT_EQ(1.f, hist.GetBinContent({0.2222}));
   EXPECT_FLOAT_EQ(2.f, hist.GetBinContent({0.3333}));
}

// Test FillN() without weights and GetBinUncertainty()
TEST(HistFillTest, FillNCoordsUncertainty)
{
   ROOT::Experimental::RH1F hist({100, 0., 1});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.1111}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.2222}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.3333}));

   hist.FillN({{0.1111}, {0.2222}, {0.3333}});
   EXPECT_FLOAT_EQ(1.f, hist.GetBinUncertainty({0.1111}));
   EXPECT_FLOAT_EQ(1.f, hist.GetBinUncertainty({0.2222}));
   EXPECT_FLOAT_EQ(1.f, hist.GetBinUncertainty({0.3333}));

   hist.FillN({{0.1111}, {0.3333}});
   EXPECT_FLOAT_EQ(std::sqrt(2.f), hist.GetBinUncertainty({0.1111}));
   EXPECT_FLOAT_EQ(1.f, hist.GetBinUncertainty({0.2222}));
   EXPECT_FLOAT_EQ(std::sqrt(2.f), hist.GetBinUncertainty({0.3333}));
}

// Test FillN() with weights and GetEntries()
TEST(HistFillTest, FillNCoordsWeightsEntries)
{
   ROOT::Experimental::RH1F hist({100, 0., 1});
   EXPECT_EQ(0, hist.GetEntries());
   hist.FillN({{0.1111}, {0.2222}, {0.3333}}, {.42f, .32f, .52f});
   EXPECT_EQ(3, hist.GetEntries());
   hist.FillN({{0.1111}, {0.3333}}, {.42f, .32f});
   EXPECT_EQ(5, hist.GetEntries());
}

// Test FillN() with weights and GetBinContent()
TEST(HistFillTest, FillNCoordsWeightsContent)
{
   ROOT::Experimental::RH1F hist({100, 0., 1});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.1111}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.2222}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.3333}));

   hist.FillN({{0.1111}, {0.2222}, {0.3333}}, {.42f, .32f, .52f});
   EXPECT_FLOAT_EQ(.42f, hist.GetBinContent({0.1111}));
   EXPECT_FLOAT_EQ(.32f, hist.GetBinContent({0.2222}));
   EXPECT_FLOAT_EQ(.52f, hist.GetBinContent({0.3333}));

   hist.FillN({{0.1111}, {0.3333}}, {.42f, .32f});
   EXPECT_FLOAT_EQ(.42f + .42f, hist.GetBinContent({0.1111}));
   EXPECT_FLOAT_EQ(.32f, hist.GetBinContent({0.2222}));
   EXPECT_FLOAT_EQ(.52f + .32f, hist.GetBinContent({0.3333}));
}

// Test FillN() and GetBinUncertainty()
TEST(HistFillTest, FillNCoordsWeightsUncertainty)
{
   ROOT::Experimental::RH1F hist({100, 0., 1});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.1111}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.2222}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.3333}));

   float weight1 = .42f;
   float weight2 = .32f;
   float weight3 = .52f;

   hist.FillN({{0.1111}, {0.2222}, {0.3333}}, {weight1, weight2, weight3});
   EXPECT_FLOAT_EQ(std::sqrt(weight1 * weight1), hist.GetBinUncertainty({0.1111}));
   EXPECT_FLOAT_EQ(std::sqrt(weight2 * weight2), hist.GetBinUncertainty({0.2222}));
   EXPECT_FLOAT_EQ(std::sqrt(weight3 * weight3), hist.GetBinUncertainty({0.3333}));

   hist.FillN({{0.1111}, {0.3333}}, {weight1, weight2});
   EXPECT_FLOAT_EQ(std::sqrt((weight1 * weight1) + (weight1 * weight1)), hist.GetBinUncertainty({0.1111}));
   EXPECT_FLOAT_EQ(std::sqrt(weight2 * weight2), hist.GetBinUncertainty({0.2222}));
   EXPECT_FLOAT_EQ(std::sqrt((weight3 * weight3) + (weight2 * weight2)), hist.GetBinUncertainty({0.3333}));
}




// RHist::Fill for 2D

// Test Fill() without weight and GetEntries()
TEST(HistFillTest, Fill2DCoordEntries)
{
   ROOT::Experimental::RH2F hist({100, 0., 1}, {10, 3., 5.});
   EXPECT_EQ(0, hist.GetEntries());
   hist.Fill({0.1111, 4.22});
   EXPECT_EQ(1, hist.GetEntries());
   hist.Fill({0.1111, 4.22});
   EXPECT_EQ(2, hist.GetEntries());
}

// Test Fill() without weight and GetBinContent()
TEST(HistFillTest, Fill2DCoordContent)
{
   ROOT::Experimental::RH2F hist({100, 0., 1}, {10, 3., 5.});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.1111, 4.22}));
   hist.Fill({0.1111, 4.22});
   EXPECT_FLOAT_EQ(1.f, hist.GetBinContent({0.1111, 4.22}));
   hist.Fill({0.1111, 4.22});
   EXPECT_FLOAT_EQ(2.f, hist.GetBinContent({0.1111, 4.22}));
}

// Test Fill() without weight and GetBinUncertainty()
TEST(HistFillTest, Fill2DCoordUncertainty)
{
   ROOT::Experimental::RH2F hist({100, 0., 1}, {10, 3., 5.});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.1111, 4.22}));
   hist.Fill({0.1111, 4.22});
   EXPECT_FLOAT_EQ(1.f, hist.GetBinUncertainty({0.1111, 4.22}));
   hist.Fill({0.1111, 4.22});
   EXPECT_FLOAT_EQ(std::sqrt(2.f), hist.GetBinUncertainty({0.1111, 4.22}));
}

// Test Fill() with weight and GetEntries()
TEST(HistFillTest, Fill2DCoordWeightEntries)
{
   ROOT::Experimental::RH2F hist({100, 0., 1}, {10, 3., 5.});
   EXPECT_EQ(0, hist.GetEntries());
   hist.Fill({0.1111, 4.22}, .42f);
   EXPECT_EQ(1, hist.GetEntries());
   hist.Fill({0.1111, 4.22}, .32f);
   EXPECT_EQ(2, hist.GetEntries());
}

// Test Fill() with weight and GetBinContent()
TEST(HistFillTest, Fill2DCoordWeightContent)
{
   ROOT::Experimental::RH2F hist({100, 0., 1}, {10, 3., 5.});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.1111, 4.22}));
   hist.Fill({0.1111, 4.22}, .42f);
   EXPECT_FLOAT_EQ(.42f, hist.GetBinContent({0.1111, 4.22}));
   hist.Fill({0.1111, 4.22}, .32f);
   EXPECT_FLOAT_EQ(.42f + .32f, hist.GetBinContent({0.1111, 4.22}));
}

// Test Fill() with weight and GetBinUncertainty()
TEST(HistFillTest, Fill2DCoordWeightUncertainty)
{
   ROOT::Experimental::RH2F hist({100, 0., 1}, {10, 3., 5.});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.1111, 4.22}));

   float weight1 = .42f;
   hist.Fill({0.1111, 4.22}, weight1);
   EXPECT_FLOAT_EQ(std::sqrt(weight1 * weight1), hist.GetBinUncertainty({0.1111, 4.22}));

   float weight2 = .32f;
   hist.Fill({0.1111, 4.22}, weight2);
   EXPECT_FLOAT_EQ(std::sqrt((weight1 * weight1) + (weight2 * weight2)), hist.GetBinUncertainty({0.1111, 4.22}));
}

// Test FillN() without weights and GetEntries()
TEST(HistFillTest, FillN2DCoordsEntries)
{
   ROOT::Experimental::RH2F hist({100, 0., 1}, {10, 3., 5.});
   EXPECT_EQ(0, hist.GetEntries());
   hist.FillN({{0.1111, 4.22}, {0.2222, 4.33}, {0.3333, 4.11}});
   EXPECT_EQ(3, hist.GetEntries());
   hist.FillN({{0.1111, 4.22}, {0.3333, 4.11}});
   EXPECT_EQ(5, hist.GetEntries());
}

// Test FillN() without weights and GetBinContent()
TEST(HistFillTest, FillN2DCoordsContent)
{
   ROOT::Experimental::RH2F hist({100, 0., 1}, {10, 3., 5.});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.1111, 4.22}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.2222, 4.33}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.3333, 4.11}));

   hist.FillN({{0.1111, 4.22}, {0.2222, 4.33}, {0.3333, 4.11}});
   EXPECT_FLOAT_EQ(1.f, hist.GetBinContent({0.1111, 4.22}));
   EXPECT_FLOAT_EQ(1.f, hist.GetBinContent({0.2222, 4.33}));
   EXPECT_FLOAT_EQ(1.f, hist.GetBinContent({0.3333, 4.11}));

   hist.FillN({{0.1111, 4.22}, {0.3333, 4.11}});
   EXPECT_FLOAT_EQ(2.f, hist.GetBinContent({0.1111, 4.22}));
   EXPECT_FLOAT_EQ(1.f, hist.GetBinContent({0.2222, 4.33}));
   EXPECT_FLOAT_EQ(2.f, hist.GetBinContent({0.3333, 4.11}));
}

// Test FillN() without weights and GetBinUncertainty()
TEST(HistFillTest, FillN2DCoordsUncertainty)
{
   ROOT::Experimental::RH2F hist({100, 0., 1}, {10, 3., 5.});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.1111, 4.22}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.2222, 4.33}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.3333, 4.11}));

   hist.FillN({{0.1111, 4.22}, {0.2222, 4.33}, {0.3333, 4.11}});
   EXPECT_FLOAT_EQ(1.f, hist.GetBinUncertainty({0.1111, 4.22}));
   EXPECT_FLOAT_EQ(1.f, hist.GetBinUncertainty({0.2222, 4.33}));
   EXPECT_FLOAT_EQ(1.f, hist.GetBinUncertainty({0.3333, 4.11}));

   hist.FillN({{0.1111, 4.22}, {0.3333, 4.11}});
   EXPECT_FLOAT_EQ(std::sqrt(2.f), hist.GetBinUncertainty({0.1111, 4.22}));
   EXPECT_FLOAT_EQ(1.f, hist.GetBinUncertainty({0.2222, 4.33}));
   EXPECT_FLOAT_EQ(std::sqrt(2.f), hist.GetBinUncertainty({0.3333, 4.11}));
}

// Test FillN() with weights and GetEntries()
TEST(HistFillTest, FillN2DCoordsWeightsEntries)
{
   ROOT::Experimental::RH2F hist({100, 0., 1}, {10, 3., 5.});
   EXPECT_EQ(0, hist.GetEntries());
   hist.FillN({{0.1111, 4.22}, {0.2222, 4.33}, {0.3333, 4.11}}, {.42f, .32f, .52f});
   EXPECT_EQ(3, hist.GetEntries());
   hist.FillN({{0.1111, 4.22}, {0.3333, 4.11}}, {.42f, .32f});
   EXPECT_EQ(5, hist.GetEntries());
}

// Test FillN() with weights and GetBinContent()
TEST(HistFillTest, FillN2DCoordsWeightsContent)
{
   ROOT::Experimental::RH2F hist({100, 0., 1}, {10, 3., 5.});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.1111, 4.22}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.2222, 4.33}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.3333, 4.11}));

   hist.FillN({{0.1111, 4.22}, {0.2222, 4.33}, {0.3333, 4.11}}, {.42f, .32f, .52f});
   EXPECT_FLOAT_EQ(.42f, hist.GetBinContent({0.1111, 4.22}));
   EXPECT_FLOAT_EQ(.32f, hist.GetBinContent({0.2222, 4.33}));
   EXPECT_FLOAT_EQ(.52f, hist.GetBinContent({0.3333, 4.11}));

   hist.FillN({{0.1111, 4.22}, {0.3333, 4.11}}, {.42f, .32f});
   EXPECT_FLOAT_EQ(.42f + .42f, hist.GetBinContent({0.1111, 4.22}));
   EXPECT_FLOAT_EQ(.32f, hist.GetBinContent({0.2222, 4.33}));
   EXPECT_FLOAT_EQ(.52f + .32f, hist.GetBinContent({0.3333, 4.11}));
}

// Test FillN() and GetBinUncertainty()
TEST(HistFillTest, FillN2DCoordsWeightsUncertainty)
{
   ROOT::Experimental::RH2F hist({100, 0., 1}, {10, 3., 5.});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.1111, 4.22}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.2222, 4.33}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.3333, 4.11}));

   float weight1 = .42f;
   float weight2 = .32f;
   float weight3 = .52f;

   hist.FillN({{0.1111, 4.22}, {0.2222, 4.33}, {0.3333, 4.11}}, {weight1, weight2, weight3});
   EXPECT_FLOAT_EQ(std::sqrt(weight1 * weight1), hist.GetBinUncertainty({0.1111, 4.22}));
   EXPECT_FLOAT_EQ(std::sqrt(weight2 * weight2), hist.GetBinUncertainty({0.2222, 4.33}));
   EXPECT_FLOAT_EQ(std::sqrt(weight3 * weight3), hist.GetBinUncertainty({0.3333, 4.11}));

   hist.FillN({{0.1111, 4.22}, {0.3333, 4.11}}, {weight1, weight2});
   EXPECT_FLOAT_EQ(std::sqrt((weight1 * weight1) + (weight1 * weight1)), hist.GetBinUncertainty({0.1111, 4.22}));
   EXPECT_FLOAT_EQ(std::sqrt(weight2 * weight2), hist.GetBinUncertainty({0.2222, 4.33}));
   EXPECT_FLOAT_EQ(std::sqrt((weight3 * weight3) + (weight2 * weight2)), hist.GetBinUncertainty({0.3333, 4.11}));
}




// RHist::Fill for 3D

// Test Fill() without weight and GetEntries()
TEST(HistFillTest, Fill3DCoordEntries)
{
   ROOT::Experimental::RH3F hist({100, 0., 1}, {10, 3., 5.}, {5, 7., 9.});
   EXPECT_EQ(0, hist.GetEntries());
   hist.Fill({0.1111, 4.22, 7.33});
   EXPECT_EQ(1, hist.GetEntries());
   hist.Fill({0.1111, 4.22, 7.33});
   EXPECT_EQ(2, hist.GetEntries());
}

// Test Fill() without weight and GetBinContent()
TEST(HistFillTest, Fill3DCoordContent)
{
   ROOT::Experimental::RH3F hist({100, 0., 1}, {10, 3., 5.}, {5, 7., 9.});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.1111, 4.22, 7.33}));
   hist.Fill({0.1111, 4.22, 7.33});
   EXPECT_FLOAT_EQ(1.f, hist.GetBinContent({0.1111, 4.22, 7.33}));
   hist.Fill({0.1111, 4.22, 7.33});
   EXPECT_FLOAT_EQ(2.f, hist.GetBinContent({0.1111, 4.22, 7.33}));
}

// Test Fill() without weight and GetBinUncertainty()
TEST(HistFillTest, Fill3DCoordUncertainty)
{
   ROOT::Experimental::RH3F hist({100, 0., 1}, {10, 3., 5.}, {5, 7., 9.});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.1111, 4.22, 7.33}));
   hist.Fill({0.1111, 4.22, 7.33});
   EXPECT_FLOAT_EQ(1.f, hist.GetBinUncertainty({0.1111, 4.22, 7.33}));
   hist.Fill({0.1111, 4.22, 7.33});
   EXPECT_FLOAT_EQ(std::sqrt(2.f), hist.GetBinUncertainty({0.1111, 4.22, 7.33}));
}

// Test Fill() with weight and GetEntries()
TEST(HistFillTest, Fill3DCoordWeightEntries)
{
   ROOT::Experimental::RH3F hist({100, 0., 1}, {10, 3., 5.}, {5, 7., 9.});
   EXPECT_EQ(0, hist.GetEntries());
   hist.Fill({0.1111, 4.22, 7.33}, .42f);
   EXPECT_EQ(1, hist.GetEntries());
   hist.Fill({0.1111, 4.22, 7.33}, .32f);
   EXPECT_EQ(2, hist.GetEntries());
}

// Test Fill() with weight and GetBinContent()
TEST(HistFillTest, Fill3DCoordWeightContent)
{
   ROOT::Experimental::RH3F hist({100, 0., 1}, {10, 3., 5.}, {5, 7., 9.});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.1111, 4.22, 7.33}));
   hist.Fill({0.1111, 4.22, 7.33}, .42f);
   EXPECT_FLOAT_EQ(.42f, hist.GetBinContent({0.1111, 4.22, 7.33}));
   hist.Fill({0.1111, 4.22, 7.33}, .32f);
   EXPECT_FLOAT_EQ(.42f + .32f, hist.GetBinContent({0.1111, 4.22, 7.33}));
}

// Test Fill() with weight and GetBinUncertainty()
TEST(HistFillTest, Fill3DCoordWeightUncertainty)
{
   ROOT::Experimental::RH3F hist({100, 0., 1}, {10, 3., 5.}, {5, 7., 9.});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.1111, 4.22, 7.33}));

   float weight1 = .42f;
   hist.Fill({0.1111, 4.22, 7.33}, weight1);
   EXPECT_FLOAT_EQ(std::sqrt(weight1 * weight1), hist.GetBinUncertainty({0.1111, 4.22, 7.33}));

   float weight2 = .32f;
   hist.Fill({0.1111, 4.22, 7.33}, weight2);
   EXPECT_FLOAT_EQ(std::sqrt((weight1 * weight1) + (weight2 * weight2)), hist.GetBinUncertainty({0.1111, 4.22, 7.33}));
}

// Test FillN() without weights and GetEntries()
TEST(HistFillTest, FillN3DCoordsEntries)
{
   ROOT::Experimental::RH3F hist({100, 0., 1}, {10, 3., 5.}, {5, 7., 9.});
   EXPECT_EQ(0, hist.GetEntries());
   hist.FillN({{0.1111, 4.22, 7.33}, {0.2222, 4.33, 7.11}, {0.3333, 4.11, 7.22}});
   EXPECT_EQ(3, hist.GetEntries());
   hist.FillN({{0.1111, 4.22, 7.33}, {0.3333, 4.11, 7.22}});
   EXPECT_EQ(5, hist.GetEntries());
}

// Test FillN() without weights and GetBinContent()
TEST(HistFillTest, FillN3DCoordsContent)
{
   ROOT::Experimental::RH3F hist({100, 0., 1}, {10, 3., 5.}, {5, 7., 9.});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.1111, 4.22, 7.33}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.2222, 4.33, 7.11}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.3333, 4.11, 7.22}));

   hist.FillN({{0.1111, 4.22, 7.33}, {0.2222, 4.33, 7.11}, {0.3333, 4.11, 7.22}});
   EXPECT_FLOAT_EQ(1.f, hist.GetBinContent({0.1111, 4.22, 7.33}));
   EXPECT_FLOAT_EQ(1.f, hist.GetBinContent({0.2222, 4.33, 7.11}));
   EXPECT_FLOAT_EQ(1.f, hist.GetBinContent({0.3333, 4.11, 7.22}));

   hist.FillN({{0.1111, 4.22, 7.33}, {0.3333, 4.11, 7.22}});
   EXPECT_FLOAT_EQ(2.f, hist.GetBinContent({0.1111, 4.22, 7.33}));
   EXPECT_FLOAT_EQ(1.f, hist.GetBinContent({0.2222, 4.33, 7.11}));
   EXPECT_FLOAT_EQ(2.f, hist.GetBinContent({0.3333, 4.11, 7.22}));
}

// Test FillN() without weights and GetBinUncertainty()
TEST(HistFillTest, FillN3DCoordsUncertainty)
{
   ROOT::Experimental::RH3F hist({100, 0., 1}, {10, 3., 5.}, {5, 7., 9.});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.1111, 4.22, 7.33}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.2222, 4.33, 7.11}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.3333, 4.11, 7.22}));
   
   hist.FillN({{0.1111, 4.22, 7.33}, {0.2222, 4.33, 7.11}, {0.3333, 4.11, 7.22}});
   EXPECT_FLOAT_EQ(1.f, hist.GetBinUncertainty({0.1111, 4.22, 7.33}));
   EXPECT_FLOAT_EQ(1.f, hist.GetBinUncertainty({0.2222, 4.33, 7.11}));
   EXPECT_FLOAT_EQ(1.f, hist.GetBinUncertainty({0.3333, 4.11, 7.22}));

   hist.FillN({{0.1111, 4.22, 7.33}, {0.3333, 4.11, 7.22}});
   EXPECT_FLOAT_EQ(std::sqrt(2.f), hist.GetBinUncertainty({0.1111, 4.22, 7.33}));
   EXPECT_FLOAT_EQ(1.f, hist.GetBinUncertainty({0.2222, 4.33, 7.11}));
   EXPECT_FLOAT_EQ(std::sqrt(2.f), hist.GetBinUncertainty({0.3333, 4.11, 7.22}));
}

// Test FillN() with weights and GetEntries()
TEST(HistFillTest, FillN3DCoordsWeightsEntries)
{
   ROOT::Experimental::RH3F hist({100, 0., 1}, {10, 3., 5.}, {5, 7., 9.});
   EXPECT_EQ(0, hist.GetEntries());
   hist.FillN({{0.1111, 4.22, 7.33}, {0.2222, 4.33, 7.11}, {0.3333, 4.11, 7.22}}, {.42f, .32f, .52f});
   EXPECT_EQ(3, hist.GetEntries());
   hist.FillN({{0.1111, 4.22, 7.33}, {0.3333, 4.11, 7.22}}, {.42f, .32f});
   EXPECT_EQ(5, hist.GetEntries());
}

// Test FillN() with weights and GetBinContent()
TEST(HistFillTest, FillN3DCoordsWeightsContent)
{
   ROOT::Experimental::RH3F hist({100, 0., 1}, {10, 3., 5.}, {5, 7., 9.});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.1111, 4.22, 7.33}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.2222, 4.33, 7.11}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinContent({0.3333, 4.11, 7.22}));

   hist.FillN({{0.1111, 4.22, 7.33}, {0.2222, 4.33, 7.11}, {0.3333, 4.11, 7.22}}, {.42f, .32f, .52f});
   EXPECT_FLOAT_EQ(.42f, hist.GetBinContent({0.1111, 4.22, 7.33}));
   EXPECT_FLOAT_EQ(.32f, hist.GetBinContent({0.2222, 4.33, 7.11}));
   EXPECT_FLOAT_EQ(.52f, hist.GetBinContent({0.3333, 4.11, 7.22}));

   hist.FillN({{0.1111, 4.22, 7.33}, {0.3333, 4.11, 7.22}}, {.42f, .32f});
   EXPECT_FLOAT_EQ(.42f + .42f, hist.GetBinContent({0.1111, 4.22, 7.33}));
   EXPECT_FLOAT_EQ(.32f, hist.GetBinContent({0.2222, 4.33, 7.11}));
   EXPECT_FLOAT_EQ(.52f + .32f, hist.GetBinContent({0.3333, 4.11, 7.22}));
}

// Test FillN() and GetBinUncertainty()
TEST(HistFillTest, FillN3DCoordsWeightsUncertainty)
{
   ROOT::Experimental::RH3F hist({100, 0., 1}, {10, 3., 5.}, {5, 7., 9.});
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.1111, 4.22, 7.33}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.2222, 4.33, 7.11}));
   EXPECT_FLOAT_EQ(.0f, hist.GetBinUncertainty({0.3333, 4.11, 7.22}));

   float weight1 = .42f;
   float weight2 = .32f;
   float weight3 = .52f;

   hist.FillN({{0.1111, 4.22, 7.33}, {0.2222, 4.33, 7.11}, {0.3333, 4.11, 7.22}}, {weight1, weight2, weight3});
   EXPECT_FLOAT_EQ(std::sqrt(weight1 * weight1), hist.GetBinUncertainty({0.1111, 4.22, 7.33}));
   EXPECT_FLOAT_EQ(std::sqrt(weight2 * weight2), hist.GetBinUncertainty({0.2222, 4.33, 7.11}));
   EXPECT_FLOAT_EQ(std::sqrt(weight3 * weight3), hist.GetBinUncertainty({0.3333, 4.11, 7.22}));

   hist.FillN({{0.1111, 4.22, 7.33}, {0.3333, 4.11, 7.22}}, {weight1, weight2});
   EXPECT_FLOAT_EQ(std::sqrt((weight1 * weight1) + (weight1 * weight1)), hist.GetBinUncertainty({0.1111, 4.22, 7.33}));
   EXPECT_FLOAT_EQ(std::sqrt(weight2 * weight2), hist.GetBinUncertainty({0.2222, 4.33, 7.11}));
   EXPECT_FLOAT_EQ(std::sqrt((weight3 * weight3) + (weight2 * weight2)), hist.GetBinUncertainty({0.3333, 4.11, 7.22}));
}
