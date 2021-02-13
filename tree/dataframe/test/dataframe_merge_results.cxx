#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx> // GetMergeableValue
#include <stdexcept>
#include <gtest/gtest.h>

using ROOT::Detail::RDF::GetMergeableValue;
using ROOT::Detail::RDF::MergeValues;

TEST(RDataFrameMergeResults, MergeCount)
{
   ROOT::RDataFrame df{100};

   auto count = df.Count();

   auto mc1 = GetMergeableValue(count);
   auto mc2 = GetMergeableValue(count);

   auto mergedptr = MergeValues(std::move(mc1), std::move(mc2));
   const auto &mc = mergedptr->GetValue();

   EXPECT_EQ(mc, 200);
}

TEST(RDataFrameMergeResults, MergeGraph)
{
   ROOT::RDataFrame df{100};

   auto col1 = df.Define("x", [](ULong64_t e) { return double(e); }, {"rdfentry_"});
   auto col2 = col1.Define("y", [](ULong64_t e) { return double(e); }, {"rdfentry_"});

   auto graph1 = col2.Graph<double, double>("x", "y");
   auto graph2 = col2.Graph<double, double>("x", "y");

   auto mg1 = GetMergeableValue(graph1);
   auto mg2 = GetMergeableValue(graph2);

   auto mergedptr = MergeValues(std::move(mg1), std::move(mg2));

   const auto &mg = mergedptr->GetValue();

   EXPECT_EQ(mg.GetN(), 200);
   EXPECT_DOUBLE_EQ(mg.GetMean(), 49.5);
}

TEST(RDataFrameMergeResults, MergeMean)
{
   ROOT::RDataFrame df1{10};
   ROOT::RDataFrame df2{20};

   auto col1 = df1.Define("x", [](ULong64_t e) { return double(e); }, {"rdfentry_"});
   auto col2 = df2.Define("x", [](ULong64_t e) { return double(e); }, {"rdfentry_"});

   auto mean1 = col1.Mean<double>("x");
   auto mean2 = col2.Mean<double>("x");

   auto mm1 = GetMergeableValue(mean1);
   auto mm2 = GetMergeableValue(mean2);

   auto mergedptr = MergeValues(std::move(mm1), std::move(mm2));
   const auto &mm = mergedptr->GetValue();

   const double truemean = (4.5 * 10 + 9.5 * 20) / (30.);

   EXPECT_DOUBLE_EQ(mm, truemean);
}

TEST(RDataFrameMergeResults, MergeHisto1DSimple)
{
   ROOT::RDataFrame df1{100};
   ROOT::RDataFrame df2{100};

   auto col1 = df1.Define("x", [](ULong64_t e) { return double(e); }, {"rdfentry_"});
   auto col2 = df2.Define("x", [](ULong64_t e) { return double(e); }, {"rdfentry_"});

   auto hist1 = col1.Histo1D<double>("x");
   auto hist2 = col2.Histo1D<double>("x");

   auto mh1 = GetMergeableValue(hist1);
   auto mh2 = GetMergeableValue(hist2);

   auto mergedptr = MergeValues(std::move(mh1), std::move(mh2));

   const auto &mh = mergedptr->GetValue();

   EXPECT_EQ(mh.GetEntries(), 200);
   EXPECT_DOUBLE_EQ(mh.GetMean(), 49.5);
}

TEST(RDataFrameMergeResults, MergeHisto1DModel)
{
   ROOT::RDataFrame df1{100};
   ROOT::RDataFrame df2{100};

   auto col1 = df1.Define("x", [](ULong64_t e) { return double(e); }, {"rdfentry_"});
   auto col2 = df2.Define("x", [](ULong64_t e) { return double(e); }, {"rdfentry_"});

   auto hist1 = col1.Histo1D<double>({"name", "title", 100, 0, 100}, "x");
   auto hist2 = col2.Histo1D<double>({"name", "title", 100, 0, 100}, "x");

   auto mh1 = GetMergeableValue(hist1);
   auto mh2 = GetMergeableValue(hist2);

   auto mergedptr = MergeValues(std::move(mh1), std::move(mh2));

   const auto &mh = mergedptr->GetValue();

   EXPECT_EQ(mh.GetEntries(), 200);
   EXPECT_DOUBLE_EQ(mh.GetMean(), 49.5);
}

TEST(RDataFrameMergeResults, MergeHisto2D)
{
   ROOT::RDataFrame df{100};

   auto col1 = df.Define("x", [](ULong64_t e) { return double(e); }, {"rdfentry_"});
   auto col2 = col1.Define("y", [](ULong64_t e) { return double(e); }, {"rdfentry_"});

   auto hist1 = col2.Histo2D<double, double>({"name", "title", 100, 0, 100, 100, 0, 100}, "x", "y");
   auto hist2 = col2.Histo2D<double, double>({"name", "title", 100, 0, 100, 100, 0, 100}, "x", "y");

   auto mh1 = GetMergeableValue(hist1);
   auto mh2 = GetMergeableValue(hist2);

   auto mergedptr = MergeValues(std::move(mh1), std::move(mh2));

   const auto &mh = mergedptr->GetValue();

   EXPECT_EQ(mh.GetEntries(), 200);
   EXPECT_DOUBLE_EQ(mh.GetMean(), 49.5);
}

TEST(RDataFrameMergeResults, MergeProfile1D)
{
   ROOT::RDataFrame df{100};

   auto col1 = df.Define("x", [](ULong64_t e) { return double(e); }, {"rdfentry_"});
   auto col2 = col1.Define("y", [](ULong64_t e) { return double(e); }, {"rdfentry_"});

   auto prof1 = col2.Profile1D<double, double>({"name", "title", 100, 0, 100, 0, 100}, "x", "y");
   auto prof2 = col2.Profile1D<double, double>({"name", "title", 100, 0, 100, 0, 100}, "x", "y");

   auto mp1 = GetMergeableValue(prof1);
   auto mp2 = GetMergeableValue(prof2);

   auto mergedptr = MergeValues(std::move(mp1), std::move(mp2));

   const auto &mp = mergedptr->GetValue();

   EXPECT_EQ(mp.GetEntries(), 200);
   EXPECT_DOUBLE_EQ(mp.GetMean(), 49.5);
   EXPECT_EQ(mp.GetMaximum(), 99);
   EXPECT_EQ(mp.GetMinimum(), 0);
}

TEST(RDataFrameMergeResults, MergeProfile2D)
{
   ROOT::RDataFrame df{100};

   auto col1 = df.Define("x", [](ULong64_t e) { return double(e); }, {"rdfentry_"});
   auto col2 = col1.Define("y", [](ULong64_t e) { return double(e); }, {"rdfentry_"});
   auto col3 = col2.Define("z", [](ULong64_t e) { return double(e); }, {"rdfentry_"});

   auto prof1 = col3.Profile2D<double, double, double>({"name", "title", 100, 0, 100, 100, 0, 100}, "x", "y", "z");
   auto prof2 = col3.Profile2D<double, double, double>({"name", "title", 100, 0, 100, 100, 0, 100}, "x", "y", "z");

   auto mp1 = GetMergeableValue(prof1);
   auto mp2 = GetMergeableValue(prof2);

   auto mergedptr = MergeValues(std::move(mp1), std::move(mp2));

   const auto &mp = mergedptr->GetValue();

   EXPECT_EQ(mp.GetEntries(), 200);
   EXPECT_DOUBLE_EQ(mp.GetMean(), 49.5);
   EXPECT_EQ(mp.GetMaximum(), 99);
   EXPECT_EQ(mp.GetMinimum(), 0);
}

TEST(RDataFrameMergeResults, MergeSum)
{
   ROOT::RDataFrame df1{100};
   ROOT::RDataFrame df2{100};

   auto col1 = df1.Define("x", [](ULong64_t e) { return double(e); }, {"rdfentry_"});
   auto col2 = df2.Define("x", [](ULong64_t e) { return double(e); }, {"rdfentry_"});

   auto sum1 = col1.Sum<double>("x");
   auto sum2 = col2.Sum<double>("x");

   auto ms1 = GetMergeableValue(sum1);
   auto ms2 = GetMergeableValue(sum2);

   auto mergedptr = MergeValues(std::move(ms1), std::move(ms2));

   const auto &ms = mergedptr->GetValue();
   EXPECT_DOUBLE_EQ(ms, 9900.);
}

TEST(RDataFrameMergeResults, MergeMin)
{
   ROOT::RDataFrame df1{100};
   ROOT::RDataFrame df2{100};

   auto filter1 = df1.Filter([](ULong64_t e) { return (e < 50); }, {"rdfentry_"});
   auto filter2 = df2.Filter([](ULong64_t e) { return (e > 50); }, {"rdfentry_"});

   auto min1 = filter1.Min<ULong64_t>("rdfentry_");
   auto min2 = filter2.Min<ULong64_t>("rdfentry_");

   auto mm1 = GetMergeableValue(min1);
   auto mm2 = GetMergeableValue(min2);

   auto mergedptr = MergeValues(std::move(mm1), std::move(mm2));

   const auto &mm = mergedptr->GetValue();
   EXPECT_EQ(mm, 0ull);
}

TEST(RDataFrameMergeResults, MergeMax)
{
   ROOT::RDataFrame df1{100};
   ROOT::RDataFrame df2{100};

   auto filter1 = df1.Filter([](ULong64_t e) { return (e < 50); }, {"rdfentry_"});
   auto filter2 = df2.Filter([](ULong64_t e) { return (e > 50); }, {"rdfentry_"});

   auto max1 = filter1.Max<ULong64_t>("rdfentry_");
   auto max2 = filter2.Max<ULong64_t>("rdfentry_");

   auto mm1 = GetMergeableValue(max1);
   auto mm2 = GetMergeableValue(max2);

   auto mergedptr = MergeValues(std::move(mm1), std::move(mm2));

   const auto &mm = mergedptr->GetValue();

   EXPECT_EQ(mm, 99ull);
}

TEST(RDataFrameMergeResults, MergeStdDev)
{
   ROOT::RDataFrame df1{100};
   ROOT::RDataFrame df2{100};

   auto col1 = df1.Define("x", [](ULong64_t e) { return double(e); }, {"rdfentry_"});
   auto col2 = df2.Define("x", [](ULong64_t e) { return double(e); }, {"rdfentry_"});

   auto dev1 = col1.StdDev<double>("x");
   auto dev2 = col2.StdDev<double>("x");

   auto md1 = GetMergeableValue(dev1);
   auto md2 = GetMergeableValue(dev2);

   auto mergedptr = MergeValues(std::move(md1), std::move(md2));

   const auto &md = mergedptr->GetValue();

   const double truestddev = 28.938506974784449; // True std dev computed separately
   EXPECT_DOUBLE_EQ(md, truestddev);
}

TEST(RDataFrameMergeResults, MergeStats)
{
   ROOT::RDataFrame df1{100};
   ROOT::RDataFrame df2{100};

   auto col1 = df1.Define("x", [](ULong64_t e) { return double(e); }, {"rdfentry_"});
   auto col2 = df2.Define("x", [](ULong64_t e) { return double(e); }, {"rdfentry_"});

   auto stats1 = col1.Stats<double>("x");
   auto stats2 = col2.Stats<double>("x");

   auto ms1 = GetMergeableValue(stats1);
   auto ms2 = GetMergeableValue(stats2);

   auto mergedptr = MergeValues(std::move(ms1), std::move(ms2));

   const auto &ms = mergedptr->GetValue();

   EXPECT_EQ(ms.GetN(), 200);
   EXPECT_DOUBLE_EQ(ms.GetMean(), 49.5);
   EXPECT_EQ(ms.GetMin(), 0);
   EXPECT_EQ(ms.GetMax(), 99);
}

TEST(RDataFrameMergeResults, Merge5Hists)
{
   ROOT::RDataFrame df{100};

   auto col1 = df.Define("x", [](ULong64_t e) { return double(e); }, {"rdfentry_"});

   auto hist1 = col1.Histo1D<double>("x");
   auto hist2 = col1.Histo1D<double>("x");
   auto hist3 = col1.Histo1D<double>("x");
   auto hist4 = col1.Histo1D<double>("x");
   auto hist5 = col1.Histo1D<double>("x");

   auto mh1 = GetMergeableValue(hist1);
   auto mh2 = GetMergeableValue(hist2);
   auto mh3 = GetMergeableValue(hist3);
   auto mh4 = GetMergeableValue(hist4);
   auto mh5 = GetMergeableValue(hist5);

   auto mergedptr = MergeValues(std::move(mh1), std::move(mh2), std::move(mh3), std::move(mh4), std::move(mh5));

   EXPECT_FALSE(mh1);
   EXPECT_FALSE(mh2);
   EXPECT_FALSE(mh3);
   EXPECT_FALSE(mh4);
   EXPECT_FALSE(mh5);
   EXPECT_TRUE(mergedptr);

   const auto &mh = mergedptr->GetValue();

   EXPECT_EQ(mh.GetEntries(), 500);
   EXPECT_DOUBLE_EQ(mh.GetMean(), 49.5);
}

TEST(RDataFrameMergeResults, WrongMergeMinMax)
{
   // Tricky case: two results of the same type with the same RMergeableValue subclass, different action helper.
   ROOT::RDataFrame df{100};

   auto min = df.Min<ULong64_t>("rdfentry_");
   auto max = df.Max<ULong64_t>("rdfentry_");

   auto mergeablemin = GetMergeableValue(min); // std::unique_ptr<RMergeableGenericAction<ULong64_t>>
   auto mergeablemax = GetMergeableValue(max); // std::unique_ptr<RMergeableGenericAction<ULong64_t>>

   try {
      auto mergedptr = MergeValues(std::move(mergeablemin), std::move(mergeablemax));
      FAIL() << "Expected std::invalid_argument error.";
   } catch (const std::invalid_argument &e) {
      EXPECT_STREQ(e.what(), "Results from different actions cannot be merged together.");
   } catch (...) {
      FAIL() << "Expected std::invalid_argument error.";
   }
}
