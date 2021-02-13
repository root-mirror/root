#include <ROOT/RDataFrame.hxx>
#include <ROOT/RArrowDS.hxx>
#include <ROOT/TSeq.hxx>
#include <TROOT.h>

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#endif
#include <arrow/builder.h>
#include <arrow/memory_pool.h>
#include <arrow/record_batch.h>
#include <arrow/table.h>
#include <arrow/testing/gtest_util.h>
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#include <gtest/gtest.h>

#include <iostream>

using namespace ROOT;
using namespace ROOT::RDF;
using namespace arrow;

std::shared_ptr<Schema> exampleSchema()
{
   return schema({field("Name", arrow::utf8()), field("Age", arrow::int64()), field("Height", arrow::float64()),
                  field("Married", arrow::boolean()), field("Babies", arrow::uint32())});
}

template <typename T>
std::shared_ptr<T> makeColumn(std::shared_ptr<Field>, std::shared_ptr<arrow::Array> array) {
  return std::make_shared<T>(field, array);
}

template <>
std::shared_ptr<arrow::ChunkedArray> makeColumn<arrow::ChunkedArray>(std::shared_ptr<Field>, std::shared_ptr<arrow::Array> array) {
  return std::make_shared<arrow::ChunkedArray>(array);
}

std::shared_ptr<Table> createTestTable()
{
   auto schema_ = exampleSchema();

   std::vector<bool> is_valid(6, true);
   std::vector<std::string> names = {"Harry", "Bob,Bob", "\"Joe\"", "Tom", " John  ", " Mary Ann "};
   std::vector<int64_t> ages = {64, 50, 40, 30, 2, 0};
   std::vector<double> heights = {180.0, 200.5, 1.7, 1.9, 1.0, 0.8};
   std::vector<bool> marriageStatus = {true, true, false, true, false, false};
   std::vector<unsigned int> babies = {1, 0, 2, 3, 4, 21};

   std::shared_ptr<Array> arrays_[5];

   arrow::ArrayFromVector<StringType, std::string>(names, &arrays_[0]);
   arrow::ArrayFromVector<Int64Type, int64_t>(ages, &arrays_[1]);
   arrow::ArrayFromVector<DoubleType, double>(heights, &arrays_[2]);
   arrow::ArrayFromVector<BooleanType, bool>(marriageStatus, &arrays_[3]);
   arrow::ArrayFromVector<UInt32Type, unsigned int>(babies, &arrays_[4]);

   using ColumnType = typename decltype(std::declval<arrow::Table>().column(0))::element_type;

   std::vector<std::shared_ptr<ColumnType>> columns_ = {
      makeColumn<ColumnType>(schema_->field(0), arrays_[0]),
      makeColumn<ColumnType>(schema_->field(1), arrays_[1]),
      makeColumn<ColumnType>(schema_->field(2), arrays_[2]),
      makeColumn<ColumnType>(schema_->field(3), arrays_[3]),
      makeColumn<ColumnType>(schema_->field(4), arrays_[4])};

   auto table_ = Table::Make(schema_, columns_);
   return table_;
}

TEST(RArrowDS, ColTypeNames)
{
   RArrowDS tds(createTestTable(), {"Name", "Age", "Height", "Married", "Babies"});
   tds.SetNSlots(1);

   auto colNames = tds.GetColumnNames();

   EXPECT_TRUE(tds.HasColumn("Name"));
   EXPECT_TRUE(tds.HasColumn("Age"));
   EXPECT_FALSE(tds.HasColumn("Address"));

   ASSERT_EQ(colNames.size(), 5U);
   EXPECT_STREQ("Height", colNames[2].c_str());
   EXPECT_STREQ("Married", colNames[3].c_str());

   EXPECT_STREQ("string", tds.GetTypeName("Name").c_str());
   EXPECT_STREQ("Long64_t", tds.GetTypeName("Age").c_str());
   EXPECT_STREQ("double", tds.GetTypeName("Height").c_str());
   EXPECT_STREQ("bool", tds.GetTypeName("Married").c_str());
   EXPECT_STREQ("UInt_t", tds.GetTypeName("Babies").c_str());
}

TEST(RArrowDS, EntryRanges)
{
   RArrowDS tds(createTestTable(), {});
   tds.SetNSlots(3U);
   tds.Initialise();

   // Still dividing in equal parts...
   auto ranges = tds.GetEntryRanges();

   ASSERT_EQ(3U, ranges.size());
   EXPECT_EQ(0U, ranges[0].first);
   EXPECT_EQ(2U, ranges[0].second);
   EXPECT_EQ(2U, ranges[1].first);
   EXPECT_EQ(4U, ranges[1].second);
   EXPECT_EQ(4U, ranges[2].first);
   EXPECT_EQ(6U, ranges[2].second);
}

TEST(RArrowDS, ColumnReaders)
{
   RArrowDS tds(createTestTable(), {});

   const auto nSlots = 3U;
   tds.SetNSlots(nSlots);
   auto valsAge = tds.GetColumnReaders<Long64_t>("Age");
   auto valsBabies = tds.GetColumnReaders<unsigned int>("Babies");

   tds.Initialise();
   auto ranges = tds.GetEntryRanges();
   auto slot = 0U;
   std::vector<Long64_t> RefsAge = {64, 50, 40, 30, 2, 0};
   std::vector<unsigned int> RefsBabies = {1, 0, 2, 3, 4, 21};
   for (auto &&range : ranges) {
      tds.InitSlot(slot, range.first);
      ASSERT_LT(slot, valsAge.size());
      for (auto i : ROOT::TSeq<int>(range.first, range.second)) {
         tds.SetEntry(slot, i);
         auto valAge = **valsAge[slot];
         EXPECT_EQ(RefsAge[i], valAge);
         auto valBabies = **valsBabies[slot];
         EXPECT_EQ(RefsBabies[i], valBabies);
      }
      slot++;
   }
}

TEST(RArrowDS, ColumnReadersString)
{
   RArrowDS tds(createTestTable(), {});

   const auto nSlots = 3U;
   tds.SetNSlots(nSlots);
   auto vals = tds.GetColumnReaders<std::string>("Name");
   tds.Initialise();
   auto ranges = tds.GetEntryRanges();
   auto slot = 0U;
   std::vector<std::string> names = {"Harry", "Bob,Bob", "\"Joe\"", "Tom", " John  ", " Mary Ann "};
   for (auto &&range : ranges) {
      tds.InitSlot(slot, range.first);
      ASSERT_LT(slot, vals.size());
      for (auto i : ROOT::TSeqU(range.first, range.second)) {
         tds.SetEntry(slot, i);
         auto val = *((std::string *)*vals[slot]);
         ASSERT_LT(i, names.size());
         EXPECT_EQ(names[i], val);
      }
      slot++;
   }
}

#ifndef NDEBUG

TEST(RArrowDS, SetNSlotsTwice)
{
   auto theTest = []() {
      RArrowDS tds(createTestTable(), {});
      tds.SetNSlots(1);
      tds.SetNSlots(1);
   };
   ASSERT_DEATH(theTest(), "Setting the number of slots even if the number of slots is different from zero.");
}
#endif

#ifdef R__B64

TEST(RArrowDS, FromARDF)
{
   std::unique_ptr<RDataSource> tds(new RArrowDS(createTestTable(), {}));
   ROOT::RDataFrame rdf(std::move(tds));
   auto max = rdf.Max<double>("Height");
   auto min = rdf.Min<double>("Height");
   auto c = rdf.Count();

   EXPECT_EQ(6U, *c);
   EXPECT_DOUBLE_EQ(200.5, *max);
   EXPECT_DOUBLE_EQ(0.8, *min);
}

TEST(RArrowDS, FromARDFWithJitting)
{
   std::unique_ptr<RDataSource> tds(new RArrowDS(createTestTable(), {}));
   ROOT::RDataFrame rdf(std::move(tds));
   auto max = rdf.Filter("Age<40").Max("Age");
   auto min = rdf.Define("Age2", "Age").Filter("Age2>30").Min("Age2");

   EXPECT_EQ(30, *max);
   EXPECT_EQ(40, *min);
}

// NOW MT!-------------
#ifdef R__USE_IMT

TEST(RArrowDS, DefineSlotCheckMT)
{
   const auto nSlots = 4U;
   ROOT::EnableImplicitMT(nSlots);

   std::vector<unsigned int> ids(nSlots, 0u);
   std::unique_ptr<RDataSource> tds(new RArrowDS(createTestTable(), {}));
   ROOT::RDataFrame d(std::move(tds));
   auto m = d.DefineSlot("x", [&](unsigned int slot) {
                ids[slot] = 1u;
                return 1;
             }).Max("x");
   EXPECT_EQ(1, *m); // just in case

   const auto nUsedSlots = std::accumulate(ids.begin(), ids.end(), 0u);
   EXPECT_GT(nUsedSlots, 0u);
   EXPECT_LE(nUsedSlots, nSlots);
}

TEST(RArrowDS, FromARDFMT)
{
   std::unique_ptr<RDataSource> tds(new RArrowDS(createTestTable(), {}));
   ROOT::RDataFrame tdf(std::move(tds));
   auto max = tdf.Max<double>("Height");
   auto min = tdf.Min<double>("Height");
   auto c = tdf.Count();

   EXPECT_EQ(6U, *c);
   EXPECT_DOUBLE_EQ(200.5, *max);
   EXPECT_DOUBLE_EQ(.8, *min);
}

TEST(RArrowDS, FromARDFWithJittingMT)
{
   std::unique_ptr<RDataSource> tds(new RArrowDS(createTestTable(), {}));
   ROOT::RDataFrame tdf(std::move(tds));
   auto max = tdf.Filter("Age<40").Max("Age");
   auto min = tdf.Define("Age2", "Age").Filter("Age2>30").Min("Age2");

   EXPECT_EQ(30, *max);
   EXPECT_EQ(40, *min);
}

#endif // R__USE_IMT

#endif // R__B64
