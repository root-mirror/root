#include "TROOT.h"
#include "ROOT/RTaskArena.hxx"
#include "ROOT/TThreadExecutor.hxx"
#include "ROOTUnitTestSupport.h"
#include "gtest/gtest.h"
#define TBB_PREVIEW_GLOBAL_CONTROL 1 // required for TBB versions preceding 2019_U4
#include "tbb/global_control.h"

#ifdef R__USE_IMT

const unsigned maxConcurrency = ROOT::Internal::LogicalCPUBandwithControl();

TEST(TBBGlobalControl, RTaskArena)
{
   if (maxConcurrency <= 1)
      return; // skip this test on systems with only 1 core
   auto gTAInstance = ROOT::Internal::GetGlobalTaskArena(maxConcurrency);
   tbb::global_control c(tbb::global_control::max_allowed_parallelism, 1u);
   gTAInstance.reset();
   ROOT_EXPECT_WARNING(gTAInstance = ROOT::Internal::GetGlobalTaskArena(maxConcurrency);,
                       "RTaskArenaWrapper",
                       "tbb::global_control is active, limiting the number of parallel workers"
                       "from this task arena available for execution.");
}

TEST(TBBGlobalControl, TThreadExecutor)
{
   // ***See them pass***
   ROOT::TThreadExecutor executor{maxConcurrency};
   executor.Map([]() { return 1; }, 10); // ParallelFor
   std::vector<double> vd{0., 1., 2.};
   executor.Reduce(vd, std::plus<double>{}); // ParallelReduce double
   std::vector<float> vf{0., 1., 2.};
   executor.Reduce(vf, std::plus<float>{}); // ParallelReduce float

   tbb::global_control c(tbb::global_control::max_allowed_parallelism, 1u);

   if (maxConcurrency <= 1)
      return; // skip this test on systems with only 1 core
   // ***See them warn***
   ROOT_EXPECT_WARNING(executor.Map([]() { return 1; }, 10), "TThreadExecutor::ParallelFor",
                       "tbb::global_control is limiting the number of parallel workers."
                       " Proceeding with 1 threads this time");

   // ParallelReduce double
   ROOT_EXPECT_WARNING(executor.Reduce(vd, std::plus<double>{}), "TThreadExecutor::ParallelReduce",
                       "tbb::global_control is limiting the number of parallel workers."
                       " Proceeding with 1 threads this time");

   // ParallelReduce float
   ROOT_EXPECT_WARNING(executor.Reduce(vf, std::plus<float>{}), "TThreadExecutor::ParallelReduce",
                       "tbb::global_control is limiting the number of parallel workers."
                       " Proceeding with 1 threads this time");
}

#endif
