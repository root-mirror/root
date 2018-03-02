#include "TRandom.h"
#include "ROOT/TDataFrame.hxx"
#include "ROOT/TSeq.hxx"
#include "gtest/gtest.h"

using namespace ROOT::Experimental;

TEST(TDataFrameReport, AnalyseCuts)
{
   // Full coverage :) ?
   TDataFrame d(128);
   TRandom r(1);
   auto gen = [&r]() { return r.Gaus(0, 1); };
   auto cut0 = [](double x) { return x > 0; };
   auto cut1 = [](double x) { return x > 0.1; };
   auto cut2 = [](double x) { return x > 0.2; };
   auto colName = "col0";
   auto dd = d.Define(colName, gen)
                .Filter(cut0, {colName}, "cut0")
                .Filter(cut1, {colName}, "cut1")
                .Filter(cut2, {colName}, "cut2");

   testing::internal::CaptureStdout();
   auto rep = dd.Report();
   std::string output = testing::internal::GetCapturedStdout();
   auto expOut = "cut0      : pass=67         all=128        --   52.344 %\n"
                 "cut1      : pass=59         all=67         --   88.060 %\n"
                 "cut2      : pass=50         all=59         --   84.746 %\n";

   EXPECT_STREQ(output.c_str(), expOut);

   std::vector<const char *> cutNames{"cut0", "cut1", "cut2"};
   std::vector<ULong64_t> allEvts{128, 67, 59};
   std::vector<ULong64_t> passEvts{67, 59, 50};
   std::vector<float> effs{52.34375f, 88.0597f, 84.745766f};
   unsigned int i = 0;

   for (auto &&cut : rep) {
      EXPECT_STREQ(cut.GetName().c_str(), cutNames[i]);
      EXPECT_EQ(cut.GetAll(), allEvts[i]);
      EXPECT_EQ(cut.GetPass(), passEvts[i]);
      ASSERT_FLOAT_EQ(cut.GetEff(), effs[i]);
      i++;
   }

   std::vector<TDF::TCutInfo> cutis{rep["cut0"], rep["cut1"], rep["cut2"]};

   for (auto j : ROOT::TSeqI(3)) {
      EXPECT_STREQ(cutis[j].GetName().c_str(), cutNames[j]);
      EXPECT_EQ(cutis[j].GetAll(), allEvts[j]);
      EXPECT_EQ(cutis[j].GetPass(), passEvts[j]);
      ASSERT_FLOAT_EQ(cutis[j].GetEff(), effs[j]);
   }

   int ret = 1;
   try {
      rep["NonExisting"];
   } catch (...) {
      ret = 0;
   }
   ASSERT_EQ(0, ret) << "No exception thrown when trying to get a non-existing cut.\n";

   ret = 1;
   try {
      rep[""];
   } catch (...) {
      ret = 0;
   }
   ASSERT_EQ(0, ret) << "No exception thrown when trying to get an unnamed cut.\n";
}

TEST(TDataFrameReport, Printing)
{
   // Full coverage :) ?
   TDataFrame d(8);
   TRandom r(1);
   auto gen = [&r]() { return r.Gaus(0, 1); };
   auto cut0 = [](double x) { return x > 0; };
   auto colName = "col0";
   auto dd = d.Define(colName, gen).Filter(cut0, {colName}, "cut0");

   testing::internal::CaptureStdout();
   dd.Report();
   std::string output0 = testing::internal::GetCapturedStdout();
   EXPECT_FALSE(output0.empty());

   testing::internal::CaptureStdout();
   auto rep = dd.Report(false);
   std::string output1 = testing::internal::GetCapturedStdout();
   EXPECT_TRUE(output1.empty());

   testing::internal::CaptureStdout();
   rep.Print();
   output1 = testing::internal::GetCapturedStdout();
   EXPECT_STREQ(output1.c_str(), output0.c_str());
}
