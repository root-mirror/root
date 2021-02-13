#include "gtest/gtest.h"

#include "ROOT/RHist.hxx"

// Test RH.. constructors and especially GetNBins(), GetNDim() and GetBinContent() from RHistImpl.hxx

// Test RH1F constructor
TEST(HistConstructorTest, ConstructorRH1F)
{
   ROOT::Experimental::RH1F hist({100, 0., 1});
   EXPECT_EQ(100 + 2, hist.GetImpl()->GetNBins());
   EXPECT_EQ(100, hist.GetImpl()->GetNBinsNoOver());
   EXPECT_EQ(2, hist.GetImpl()->GetNOverflowBins());
   EXPECT_EQ(0, hist.GetImpl()->GetBinContent(42));
   EXPECT_EQ(1, hist.GetImpl()->GetNDim());
}

// Test RH1D constructor
TEST(HistConstructorTest, ConstructorRH1D)
{
   ROOT::Experimental::RH1D hist({100, 0, 1});
   EXPECT_EQ(100 + 2, hist.GetImpl()->GetNBins());
   EXPECT_EQ(100, hist.GetImpl()->GetNBinsNoOver());
   EXPECT_EQ(2, hist.GetImpl()->GetNOverflowBins());
   EXPECT_EQ(0, hist.GetImpl()->GetBinContent(42));
   EXPECT_EQ(1, hist.GetImpl()->GetNDim());
}

// Test RH1I constructor
TEST(HistConstructorTest, ConstructorRH1I)
{
   ROOT::Experimental::RH1I hist({100, 0, 1});
   EXPECT_EQ(100 + 2, hist.GetImpl()->GetNBins());
   EXPECT_EQ(100, hist.GetImpl()->GetNBinsNoOver());
   EXPECT_EQ(2, hist.GetImpl()->GetNOverflowBins());
   EXPECT_EQ(0, hist.GetImpl()->GetBinContent(42));
   EXPECT_EQ(1, hist.GetImpl()->GetNDim());
}

// Test RH1C constructor
TEST(HistConstructorTest, ConstructorRH1C)
{
   ROOT::Experimental::RH1C hist({100, 0, 1});
   EXPECT_EQ(100 + 2, hist.GetImpl()->GetNBins());
   EXPECT_EQ(100, hist.GetImpl()->GetNBinsNoOver());
   EXPECT_EQ(2, hist.GetImpl()->GetNOverflowBins());
   EXPECT_EQ(0, hist.GetImpl()->GetBinContent(42));
   EXPECT_EQ(1, hist.GetImpl()->GetNDim());
}

// Test RH1LL constructor
TEST(HistConstructorTest, ConstructorRH1LL)
{
   ROOT::Experimental::RH1LL hist({100, 0, 1});
   EXPECT_EQ(100 + 2, hist.GetImpl()->GetNBins());
   EXPECT_EQ(100, hist.GetImpl()->GetNBinsNoOver());
   EXPECT_EQ(2, hist.GetImpl()->GetNOverflowBins());
   EXPECT_EQ(0, hist.GetImpl()->GetBinContent(42));
   EXPECT_EQ(1, hist.GetImpl()->GetNDim());
}



// Test RH2F constructor
TEST(HistConstructorTest, ConstructorRH2F)
{
   ROOT::Experimental::RH2F hist({100, 0., 1}, {{0., 1., 2., 3., 10.}});
   EXPECT_EQ((100 + 2) * (4 + 2), hist.GetImpl()->GetNBins());
   EXPECT_EQ(100 * 4, hist.GetImpl()->GetNBinsNoOver());
   EXPECT_EQ(((100 + 2) * (4 + 2)) - (100 * 4), hist.GetImpl()->GetNOverflowBins());
   EXPECT_EQ(0, hist.GetImpl()->GetBinContent(42));
   EXPECT_EQ(2, hist.GetImpl()->GetNDim());
}

// Test RH2D constructor
TEST(HistConstructorTest, ConstructorRH2D)
{
   ROOT::Experimental::RH2D hist({100, 0, 1}, {{0, 1, 2, 3, 10}});
   EXPECT_EQ((100 + 2) * (4 + 2), hist.GetImpl()->GetNBins());
   EXPECT_EQ(100 * 4, hist.GetImpl()->GetNBinsNoOver());
   EXPECT_EQ(((100 + 2) * (4 + 2)) - (100 * 4), hist.GetImpl()->GetNOverflowBins());
   EXPECT_EQ(0, hist.GetImpl()->GetBinContent(42));
   EXPECT_EQ(2, hist.GetImpl()->GetNDim());
}

// Test RH2I constructor
TEST(HistConstructorTest, ConstructorRH2I)
{
   ROOT::Experimental::RH2I hist({100, 0, 1}, {{0, 1, 2, 3, 10}});
   EXPECT_EQ((100 + 2) * (4 + 2), hist.GetImpl()->GetNBins());
   EXPECT_EQ(100 * 4, hist.GetImpl()->GetNBinsNoOver());
   EXPECT_EQ(((100 + 2) * (4 + 2)) - (100 * 4), hist.GetImpl()->GetNOverflowBins());
   EXPECT_EQ(0, hist.GetImpl()->GetBinContent(42));
   EXPECT_EQ(2, hist.GetImpl()->GetNDim());
}

// Test RH2C constructor
TEST(HistConstructorTest, ConstructorRH2C)
{
   ROOT::Experimental::RH2C hist({100, 0, 1}, {{0, 1, 2, 3, 10}});
   EXPECT_EQ((100 + 2) * (4 + 2), hist.GetImpl()->GetNBins());
   EXPECT_EQ(100 * 4, hist.GetImpl()->GetNBinsNoOver());
   EXPECT_EQ(((100 + 2) * (4 + 2)) - (100 * 4), hist.GetImpl()->GetNOverflowBins());
   EXPECT_EQ(0, hist.GetImpl()->GetBinContent(42));
   EXPECT_EQ(2, hist.GetImpl()->GetNDim());
}

// Test RH2LL constructor
TEST(HistConstructorTest, ConstructorRH2LL)
{
   ROOT::Experimental::RH2LL hist({100, 0, 1}, {{0, 1, 2, 3, 10}});
   EXPECT_EQ((100 + 2) * (4 + 2), hist.GetImpl()->GetNBins());
   EXPECT_EQ(100 * 4, hist.GetImpl()->GetNBinsNoOver());
   EXPECT_EQ(((100 + 2) * (4 + 2)) - (100 * 4), hist.GetImpl()->GetNOverflowBins());
   EXPECT_EQ(0, hist.GetImpl()->GetBinContent(42));
   EXPECT_EQ(2, hist.GetImpl()->GetNDim());
}



// Test RH3F constructor
TEST(HistConstructorTest, ConstructorRH3F)
{
   ROOT::Experimental::RH3F hist({100, 0., 1}, {{0., 1., 2., 3., 10.}}, {{0., 1., 2., 3., 4., 10.}});
   EXPECT_EQ((100 + 2) * (4 + 2) * (5 + 2), hist.GetImpl()->GetNBins());
   EXPECT_EQ(100 * 4 * 5, hist.GetImpl()->GetNBinsNoOver());
   EXPECT_EQ(((100 + 2) * (4 + 2) * (5 + 2)) - (100 * 4 * 5), hist.GetImpl()->GetNOverflowBins());
   EXPECT_EQ(0, hist.GetImpl()->GetBinContent(42));
   EXPECT_EQ(3, hist.GetImpl()->GetNDim());
}

// Test RH3D constructor
TEST(HistConstructorTest, ConstructorRH3D)
{
   ROOT::Experimental::RH3D hist({100, 0, 1}, {{0, 1, 2, 3, 10}}, {{0, 1, 2, 3, 4, 10}});
   EXPECT_EQ((100 + 2) * (4 + 2) * (5 + 2), hist.GetImpl()->GetNBins());
   EXPECT_EQ(100 * 4 * 5, hist.GetImpl()->GetNBinsNoOver());
   EXPECT_EQ(((100 + 2) * (4 + 2) * (5 + 2)) - (100 * 4 * 5), hist.GetImpl()->GetNOverflowBins());
   EXPECT_EQ(0, hist.GetImpl()->GetBinContent(42));
   EXPECT_EQ(3, hist.GetImpl()->GetNDim());
}

// Test RH3I constructor
TEST(HistConstructorTest, ConstructorRH3I)
{
   ROOT::Experimental::RH3I hist({100, 0, 1}, {{0, 1, 2, 3, 10}}, {{0, 1, 2, 3, 4, 10}});
   EXPECT_EQ((100 + 2) * (4 + 2) * (5 + 2), hist.GetImpl()->GetNBins());
   EXPECT_EQ(100 * 4 * 5, hist.GetImpl()->GetNBinsNoOver());
   EXPECT_EQ(((100 + 2) * (4 + 2) * (5 + 2)) - (100 * 4 * 5), hist.GetImpl()->GetNOverflowBins());
   EXPECT_EQ(0, hist.GetImpl()->GetBinContent(42));
   EXPECT_EQ(3, hist.GetImpl()->GetNDim());
}

// Test RH3C constructor
TEST(HistConstructorTest, ConstructorRH3C)
{
   ROOT::Experimental::RH3C hist({100, 0, 1}, {{0, 1, 2, 3, 10}}, {{0, 1, 2, 3, 4, 10}});
   EXPECT_EQ((100 + 2) * (4 + 2) * (5 + 2), hist.GetImpl()->GetNBins());
   EXPECT_EQ(100 * 4 * 5, hist.GetImpl()->GetNBinsNoOver());
   EXPECT_EQ(((100 + 2) * (4 + 2) * (5 + 2)) - (100 * 4 * 5), hist.GetImpl()->GetNOverflowBins());
   EXPECT_EQ(0, hist.GetImpl()->GetBinContent(42));
   EXPECT_EQ(3, hist.GetImpl()->GetNDim());
}

// Test RH3LL constructor
TEST(HistConstructorTest, ConstructorRH3LL)
{
   ROOT::Experimental::RH3LL hist({100, 0, 1}, {{0, 1, 2, 3, 10}}, {{0, 1, 2, 3, 4, 10}});
   EXPECT_EQ((100 + 2) * (4 + 2) * (5 + 2), hist.GetImpl()->GetNBins());
   EXPECT_EQ(100 * 4 * 5, hist.GetImpl()->GetNBinsNoOver());
   EXPECT_EQ(((100 + 2) * (4 + 2) * (5 + 2)) - (100 * 4 * 5), hist.GetImpl()->GetNOverflowBins());
   EXPECT_EQ(0, hist.GetImpl()->GetBinContent(42));
   EXPECT_EQ(3, hist.GetImpl()->GetNDim());
}



// Test RHist constructor with 4 dimensions
TEST(HistConstructorTest, ConstructorRH4)
{
   ROOT::Experimental::RHist<4, float> hist({{ {2, 0., 1}, {{0., 1., 2., 3., 10.}}, {{0., 1., 2., 3., 4., 10.}}, {3, 0., 5.} }});
   EXPECT_EQ((2 + 2) * (4 + 2) * (5 + 2) * (3 + 2), hist.GetImpl()->GetNBins());
   EXPECT_EQ(2 * 4 * 5 * 3, hist.GetImpl()->GetNBinsNoOver());
   EXPECT_EQ(((2 + 2) * (4 + 2) * (5 + 2) * (3 + 2)) - (2 * 4 * 5 * 3), hist.GetImpl()->GetNOverflowBins());
   EXPECT_EQ(0, hist.GetImpl()->GetBinContent(42));
   EXPECT_EQ(4, hist.GetImpl()->GetNDim());
}



// Test RHist constructor with title
TEST(HistConstructorTest, ConstructorTitle)
{
   ROOT::Experimental::RH1F hist1d("1D hist", {{ {100, 0., 1} }});
   ROOT::Experimental::RH2F hist2d("2D hist", {{ {100, 0., 1}, {{0., 1., 2., 3., 10.}} }});
   ROOT::Experimental::RH3F hist3d("3D hist", {{ {100, 0., 1}, {{0., 1., 2., 3., 10.}}, {{0., 1., 2., 3., 4., 10.}} }});

   EXPECT_EQ("1D hist", hist1d.GetImpl()->GetTitle());
   EXPECT_EQ("2D hist", hist2d.GetImpl()->GetTitle());
   EXPECT_EQ("3D hist", hist3d.GetImpl()->GetTitle());

   // Checking that adding a title doesn't create inconsistancies 
   
   EXPECT_EQ((100 + 2), hist1d.GetImpl()->GetNBins());
   EXPECT_EQ(100, hist1d.GetImpl()->GetNBinsNoOver());
   EXPECT_EQ(2, hist1d.GetImpl()->GetNOverflowBins());
   EXPECT_EQ(0, hist1d.GetImpl()->GetBinContent(42));
   EXPECT_EQ(1, hist1d.GetImpl()->GetNDim());

   EXPECT_EQ((100 + 2) * (4 + 2), hist2d.GetImpl()->GetNBins());
   EXPECT_EQ(100 * 4, hist2d.GetImpl()->GetNBinsNoOver());
   EXPECT_EQ(((100 + 2) * (4 + 2)) - (100 * 4), hist2d.GetImpl()->GetNOverflowBins());
   EXPECT_EQ(0, hist2d.GetImpl()->GetBinContent(42));
   EXPECT_EQ(2, hist2d.GetImpl()->GetNDim());

   EXPECT_EQ((100 + 2) * (4 + 2) * (5 + 2), hist3d.GetImpl()->GetNBins());
   EXPECT_EQ(100 * 4 * 5, hist3d.GetImpl()->GetNBinsNoOver());
   EXPECT_EQ(((100 + 2) * (4 + 2) * (5 + 2)) - (100 * 4 * 5), hist3d.GetImpl()->GetNOverflowBins());
   EXPECT_EQ(0, hist3d.GetImpl()->GetBinContent(42));
   EXPECT_EQ(3, hist3d.GetImpl()->GetNDim());
}
