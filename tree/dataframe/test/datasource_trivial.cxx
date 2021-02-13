#include <ROOT/RTrivialDS.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDataSource.hxx>
#include <ROOT/TSeq.hxx>
#include <TROOT.h>
#include <TSystem.h>

#include "gtest/gtest.h"

using namespace ROOT;
using namespace ROOT::RDF;

TEST(RTrivialDS, ColTypeNames)
{
   RTrivialDS tds(32);
   tds.SetNSlots(1);

   auto colName = tds.GetColumnNames()[0]; // We know it's one.
   EXPECT_STREQ("col0", colName.c_str());
   EXPECT_STREQ("ULong64_t", tds.GetTypeName("col0").c_str());

   EXPECT_TRUE(tds.HasColumn("col0"));
   EXPECT_FALSE(tds.HasColumn("col1"));
}

TEST(RTrivialDS, EntryRanges)
{
   RTrivialDS tds(32);
   const auto nSlots = 4U;
   tds.SetNSlots(nSlots);
   tds.Initialise();
   auto ranges = tds.GetEntryRanges();

   EXPECT_EQ(4U, ranges.size());
   EXPECT_EQ(0U, ranges[0].first);
   EXPECT_EQ(8U, ranges[0].second);
   EXPECT_EQ(8U, ranges[1].first);
   EXPECT_EQ(16U, ranges[1].second);
   EXPECT_EQ(16U, ranges[2].first);
   EXPECT_EQ(24U, ranges[2].second);
   EXPECT_EQ(24U, ranges[3].first);
   EXPECT_EQ(32U, ranges[3].second);
}

TEST(RTrivialDS, ColumnReaders)
{
   RTrivialDS tds(32);
   const auto nSlots = 4U;
   tds.SetNSlots(nSlots);
   auto vals = tds.GetColumnReaders<ULong64_t>("col0");
   tds.Initialise();
   auto ranges = tds.GetEntryRanges();
   auto slot = 0U;
   for (auto &&range : ranges) {
      for (auto i : ROOT::TSeq<ULong64_t>(range.first, range.second)) {
         tds.SetEntry(slot, i);
         auto val = **vals[slot];
         EXPECT_EQ(i, val);
      }
      slot++;
   }
}

TEST(RTrivialDS, ColumnReadersWrongType)
{
   RTrivialDS tds(32);
   const auto nSlots = 4U;
   tds.SetNSlots(nSlots);
   int res = 1;
   try {
      auto vals = tds.GetColumnReaders<float>("col0");
   } catch (const std::runtime_error &e) {
      EXPECT_STREQ("The type specified for the column \"col0\" is not ULong64_t.", e.what());
      res = 0;
   }
   EXPECT_EQ(0, res);
}

#ifndef NDEBUG
TEST(RTrivialDS, SetNSlotsTwice)
{
   auto theTest = []() {
      RTrivialDS tds(1);
      tds.SetNSlots(1);
      tds.SetNSlots(1);
   };
   ASSERT_DEATH(theTest(), "Setting the number of slots even if the number of slots is different from zero.");
}
#endif

TEST(RTrivialDS, FromARDF)
{
   std::unique_ptr<RDataSource> tds(new RTrivialDS(32));
   RDataFrame tdf(std::move(tds));
   auto max = tdf.Max<ULong64_t>("col0");
   auto min = tdf.Min<ULong64_t>("col0");
   auto c = tdf.Count();
   auto max2 = tdf.Filter([](ULong64_t col0) { return col0 < 10; }, {"col0"}).Max<ULong64_t>("col0");
   auto min2 = tdf.Filter([](ULong64_t col0) { return col0 > 10; }, {"col0"})
                  .Define("j", [](ULong64_t col0) { return col0 * 2; }, {"col0"})
                  .Min<ULong64_t>("j");

   EXPECT_EQ(32U, *c);
   EXPECT_DOUBLE_EQ(31., *max);
   EXPECT_DOUBLE_EQ(0., *min);
   EXPECT_DOUBLE_EQ(9., *max2);
   EXPECT_DOUBLE_EQ(22., *min2);
}

TEST(RTrivialDS, SkipEntries)
{
   auto nevts = 8;
   RTrivialDS tdsOdd(nevts, true);
   tdsOdd.SetNSlots(1);
   auto retVal = false;
   for (auto ievt : ROOT::TSeqI(nevts)) {
      EXPECT_EQ(retVal, tdsOdd.SetEntry(0, ievt));
      retVal = !retVal;
   }

   RTrivialDS tdsAll(nevts);
   tdsAll.SetNSlots(1);
   retVal = true;
   for (auto ievt : ROOT::TSeqI(nevts)) {
      EXPECT_TRUE(retVal == tdsAll.SetEntry(0, ievt));
   }


   auto tdfOdd = ROOT::RDF::MakeTrivialDataFrame(20ULL, true);
   EXPECT_EQ(*tdfOdd.Count(), 10ULL);
   auto tdfAll = ROOT::RDF::MakeTrivialDataFrame(20ULL);
   EXPECT_EQ(*tdfAll.Count(), 20ULL);
}

// Test for issue #6455, "RDS does not early-quit event loops when all Ranges are exhausted"
TEST(RTrivialDS, EarlyQuitWithRange)
{
   // this is a data-source that returns an infinite amount of entries
   auto df = ROOT::RDF::MakeTrivialDataFrame();
   // here we check that the event loop early-quits when the Range is exhausted
   EXPECT_EQ(df.Range(10).Count().GetValue(), 10);
}

#ifdef R__B64

TEST(RTrivialDS, FromARDFWithJitting)
{
   std::unique_ptr<RDataSource> tds(new RTrivialDS(32));
   RDataFrame tdf(std::move(tds));
   auto max = tdf.Filter("col0 < 10").Max("col0");
   auto min = tdf.Filter("col0 > 10").Define("j", "col0*2").Min("j");

   EXPECT_DOUBLE_EQ(9., *max);
   EXPECT_DOUBLE_EQ(22., *min);
}

// NOW MT!-------------

#ifdef R__USE_IMT

TEST(RTrivialDS, DefineSlotCheckMT)
{
   const auto nSlots = 4U;
   ROOT::EnableImplicitMT(nSlots);

   std::vector<unsigned int> ids(nSlots, 0u);
   std::unique_ptr<RDataSource> tds(new RTrivialDS(nSlots));
   RDataFrame d(std::move(tds));
   auto m = d.DefineSlot("x", [&](unsigned int slot) {
                ids[slot] = 1u;
                return 1;
             }).Max("x");
   EXPECT_EQ(1, *m); // just in case

   const auto nUsedSlots = std::accumulate(ids.begin(), ids.end(), 0u);
   EXPECT_GT(nUsedSlots, 0u);
   EXPECT_LE(nUsedSlots, nSlots);
   ROOT::DisableImplicitMT();
}

TEST(RTrivialDS, FromARDFMT)
{
   std::unique_ptr<RDataSource> tds(new RTrivialDS(320));
   RDataFrame tdf(std::move(tds));
   auto max = tdf.Max<ULong64_t>("col0");
   auto min = tdf.Min<ULong64_t>("col0");
   auto c = tdf.Count();
   auto max2 = tdf.Filter([](ULong64_t col0) { return col0 < 10; }, {"col0"}).Max<ULong64_t>("col0");
   auto min2 = tdf.Filter([](ULong64_t col0) { return col0 > 10; }, {"col0"})
                  .Define("j", [](ULong64_t col0) { return col0 * 2; }, {"col0"})
                  .Min<ULong64_t>("j");

   EXPECT_EQ(320U, *c);
   EXPECT_DOUBLE_EQ(319., *max);
   EXPECT_DOUBLE_EQ(0., *min);
   EXPECT_DOUBLE_EQ(9., *max2);
   EXPECT_DOUBLE_EQ(22., *min2);
}

TEST(RTrivialDS, FromARDFWithJittingMT)
{
   std::unique_ptr<RDataSource> tds(new RTrivialDS(320));
   RDataFrame tdf(std::move(tds));
   auto max = tdf.Filter("col0 < 10").Max("col0");
   auto min = tdf.Filter("col0 > 10").Define("j", "col0*2").Min("j");

   EXPECT_DOUBLE_EQ(9., *max);
   EXPECT_DOUBLE_EQ(22., *min);
}

TEST(RTrivialDS, Snapshot)
{
   std::unique_ptr<RDataSource> tds(new RTrivialDS(10));
   const auto fname = "datasource_trivial_snapshot.root";
   RDataFrame tdf(std::move(tds));
   auto tdf2 = tdf.Snapshot("t", fname, "col0");
   auto c = tdf2->Take<ULong64_t>("col0");
   auto i = 0u;
   for (auto e : c) {
      EXPECT_EQ(e, i);
      ++i;
   }
   gSystem->Unlink(fname);
}

TEST(RTrivialDS, Cache)
{
   std::unique_ptr<RDataSource> tds(new RTrivialDS(10));
   RDataFrame tdf(std::move(tds));
   auto tdfCached = tdf.Cache("col0");
   auto c = tdfCached.Take<ULong64_t>("col0");
   auto i = 0u;
   for (auto e : c) {
      EXPECT_EQ(e, i);
      ++i;
   }
}

TEST(RTrivialDS, SkipEntriesMT)
{
   auto tdfOdd = ROOT::RDF::MakeTrivialDataFrame(80ULL, true);
   EXPECT_EQ(*tdfOdd.Count(), 40ULL);
   auto tdfAll = ROOT::RDF::MakeTrivialDataFrame(80ULL);
   EXPECT_EQ(*tdfAll.Count(), 80ULL);
}

#endif // R__USE_IMT

#endif // R__B64
