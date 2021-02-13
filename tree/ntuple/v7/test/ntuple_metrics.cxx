#include "ntuple_test.hxx"

#include <cmath>

TEST(Metrics, Counters)
{
   RNTupleMetrics metrics("test");
   EXPECT_FALSE(metrics.IsEnabled());

   RNTuplePlainCounter *ctrOne = nullptr;
   RNTupleAtomicCounter *ctrTwo = nullptr;
   ctrOne = metrics.MakeCounter<RNTuplePlainCounter *>("plain", "s", "example 1");
   ctrTwo = metrics.MakeCounter<RNTupleAtomicCounter *>("atomic", "s", "example 2");
   ASSERT_NE(nullptr, ctrOne);
   ASSERT_NE(nullptr, ctrTwo);
   EXPECT_FALSE(ctrOne->IsEnabled());
   EXPECT_FALSE(ctrTwo->IsEnabled());

   EXPECT_EQ(0, ctrOne->GetValue());
   ctrOne->Inc();
   ctrTwo->Inc();
   EXPECT_EQ(1, ctrOne->GetValue());
   EXPECT_EQ(0, ctrTwo->GetValue());
   metrics.Enable();
   EXPECT_TRUE(metrics.IsEnabled());
   ctrTwo->Inc();
   EXPECT_EQ(1, ctrTwo->XAdd(5));
   EXPECT_EQ(1, ctrOne->GetValue());
   EXPECT_EQ(6, ctrTwo->GetValue());

   RNTupleCalcPerf *ctrCalc = metrics.MakeCounter<RNTupleCalcPerf *>("calc", "s/s", "example 1/example2",
      metrics, [](const RNTupleMetrics &met) -> std::pair<bool, double> {
         auto ctr1 = met.GetCounter("test.plain");
         EXPECT_NE(ctr1, nullptr);
         auto ctr2 = met.GetCounter("test.atomic");
         EXPECT_NE(ctr2, nullptr);
         EXPECT_NE(ctr2->GetValueAsInt(), 0);
         return {true, (1.*ctr1->GetValueAsInt()) / ctr2->GetValueAsInt()};
      }
   );
   EXPECT_NE(ctrCalc, nullptr);
   EXPECT_DOUBLE_EQ(ctrCalc->GetValue(), 1./6.);
   EXPECT_NE(ctrCalc->ToString().find("calc"), std::string::npos);

   RNTupleCalcPerf *ctrCalcBad = metrics.MakeCounter<RNTupleCalcPerf *>("calcBad", "apples or oranges", "just bad",
      metrics, [](const RNTupleMetrics &) -> std::pair<bool, double> {
         return {false, 42.};
      }
   );
   EXPECT_NE(ctrCalcBad, nullptr);
   EXPECT_TRUE(std::isnan(ctrCalcBad->GetValue()));
   EXPECT_NE(ctrCalcBad->ToString(), ""); // whatever it is, it should not be empty or crash.
}

TEST(Metrics, Nested)
{
   RNTupleMetrics inner("inner");
   auto ctr = inner.MakeCounter<RNTuplePlainCounter *>("plain", "s", "example 1");

   RNTupleMetrics outer("outer");
   outer.ObserveMetrics(inner);

   outer.Enable();
   EXPECT_TRUE(ctr->IsEnabled());
   ctr->SetValue(42);

   EXPECT_EQ(nullptr, outer.GetCounter("a.b.c.d"));
   EXPECT_EQ(nullptr, outer.GetCounter("outer.xyz"));
   auto ctest = outer.GetCounter("outer.inner.plain");
   ASSERT_EQ(ctr, ctest);
   EXPECT_EQ(std::string("42"), ctest->GetValueAsString());
}

TEST(Metrics, Timer)
{
   RNTupleAtomicCounter ctrWallTime("wall time", "ns", "");
   ROOT::Experimental::Detail::RNTupleTickCounter<RNTupleAtomicCounter> ctrCpuTicks("cpu time", "ns", "");
   {
      RNTupleAtomicTimer timer(ctrWallTime, ctrCpuTicks);
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
   }
   EXPECT_EQ(0U, ctrWallTime.GetValue());
   EXPECT_EQ(0U, ctrCpuTicks.GetValue());
   ctrWallTime.Enable();
   ctrCpuTicks.Enable();
   {
      RNTupleAtomicTimer timer(ctrWallTime, ctrCpuTicks);
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
   }
   EXPECT_GT(ctrWallTime.GetValue(), 0U);
}
