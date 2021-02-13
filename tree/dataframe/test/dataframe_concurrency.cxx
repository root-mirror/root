#include <ROOT/RConfig.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/TThreadExecutor.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <atomic>
#include <chrono>
#include <thread>
#include <vector>
#include <string>

#include "gtest/gtest.h"

#ifndef R__B64
static const auto NUM_THREADS = 8u;
#else
static const auto NUM_THREADS = 0u;
#endif

#ifdef R__USE_IMT
TEST(RDFConcurrency, NestedParallelismBetweenDefineCalls)
{
   ROOT::EnableImplicitMT(NUM_THREADS);

   // this Define will return unique values from 0 to nEntries - 1 (over all threads)
   const auto nEntries = 100000u;
   std::atomic_int i(0);
   auto df = ROOT::RDataFrame(nEntries).Define("i", [&] { return i++; });

   // this lambda will be used to introduce nested parallelism via a dummy Filter
   auto manysleeps = [&] {
      ROOT::TThreadExecutor().Foreach(
         [] { std::this_thread::sleep_for(std::chrono::milliseconds(std::rand() / RAND_MAX * 200)); }, 8);
      return true;
   };

   // first Take triggers the computation of i, then the Filter executes, then Take accesses the cached value of i
   // hopefully it won't have changed in the meanwhile!
   auto res1 = df.Take<int>("i");
   auto res2 = df.Filter(manysleeps).Take<int>("i");

   // check both Takes return vectors with all numbers from 0 to nEntries - 1
   auto checkAllNumbersAreThere = [](std::vector<int> &vec) {
      std::sort(vec.begin(), vec.end());
      const int s = vec.size();
      for (auto j = 0; j < s; ++j)
         EXPECT_EQ(vec[j], j);
   };

   checkAllNumbersAreThere(*res1);
   checkAllNumbersAreThere(*res2);

   ROOT::DisableImplicitMT();
}

// ROOT-10346:
// [DF] Warn on mismatch between slot pool size and effective number of slots
TEST(RDFSimpleTests, ThrowOnPoolSizeMismatch)
{
   // pool created after RDF
   {
      ROOT::RDataFrame df(1);
      ROOT::EnableImplicitMT(3);
      try {
         df.Count().GetValue();
      } catch (const std::runtime_error &e) {
         const auto expected_msg = "RLoopManager::Run: when the RDataFrame was constructed the number of slots required "
                                   "was 1, but when starting the event loop it was 3. Maybe EnableImplicitMT() was "
                                   "called after the RDataFrame was constructed?";
         EXPECT_STREQ(e.what(), expected_msg);
      }
      ROOT::DisableImplicitMT();
   }

   // pool deleted after RDF creation
   {
      ROOT::EnableImplicitMT(3);
      ROOT::RDataFrame df(1);
      ROOT::DisableImplicitMT();
      try {
         df.Count().GetValue();
      } catch (const std::runtime_error &e) {
         const auto expected_msg = "RLoopManager::Run: when the RDataFrame was constructed the number of slots required "
                                   "was 3, but when starting the event loop it was 1. Maybe DisableImplicitMT() was "
                                   "called after the RDataFrame was constructed?";
         EXPECT_STREQ(e.what(), expected_msg);
      }
   }
}


void SimpleParallelRDFs()
{
   // Run the RDF construction and the event loop in parallel
   auto func = [] {
         std::this_thread::sleep_for(std::chrono::milliseconds(std::rand() / RAND_MAX * 200));
         ROOT::RDataFrame df(10);
         return df.Define("x", "rdfentry_").Mean("x").GetValue();
      };

   ROOT::TThreadExecutor pool;
   auto res = pool.Map(func, 64);

   const auto ref = func();

   for(auto& r : res) {
      EXPECT_EQ(r, ref);
   }
}

TEST(RDFConcurrency, SimpleParallelRDFsEnableThreadSafety)
{
   ROOT::EnableThreadSafety();
   SimpleParallelRDFs();
}

TEST(RDFConcurrency, SimpleParallelRDFsEnableImplicitMT)
{
   ROOT::EnableImplicitMT(NUM_THREADS);
   SimpleParallelRDFs();
   ROOT::DisableImplicitMT();
}

void SimpleParallelRDFLoops()
{
   // Run only the event loops in parallel
   auto create_df = [] {
         std::this_thread::sleep_for(std::chrono::milliseconds(std::rand() / RAND_MAX * 200));
         ROOT::RDataFrame df(10);
         return df.Define("x", "rdfentry_").Mean("x");
      };

   std::vector<ROOT::RDF::RResultPtr<double>> vals(64);
   for (auto i = 0u; i < vals.size(); i++)
       vals[i] = create_df();

   ROOT::TThreadExecutor pool;
   auto func = [](ROOT::RDF::RResultPtr<double> rptr){ return rptr.GetValue(); };
   auto res = pool.Map(func, vals);

   auto ref = func(create_df());

   for(auto& r : res) {
      EXPECT_EQ(r, ref);
   }
}

TEST(RDFConcurrency, SimpleParallelRDFLoopsEnableThreadSafety)
{
   ROOT::EnableThreadSafety();
   SimpleParallelRDFLoops();
}

TEST(RDFConcurrency, SimpleParallelRDFLoopsEnableImplicitMT)
{
   ROOT::EnableImplicitMT(NUM_THREADS);
   SimpleParallelRDFLoops();
   ROOT::DisableImplicitMT();
}

void ParallelRDFSnapshots()
{
   // Run the RDF construction and the event loop in parallel
   const auto nevts = 100u;
   auto func = [&] (int i) {
         std::this_thread::sleep_for(std::chrono::milliseconds(std::rand() / RAND_MAX * 200));
         ROOT::RDataFrame df(nevts);
         df.Define("x", "(int)rdfentry_").Snapshot("tree", "dataframe_parallel_snapshots_" + std::to_string(i) + ".root");
      };

   ROOT::TThreadExecutor pool;
   std::vector<int> vals = {0, 1, 2, 3, 4, 5, 6, 7};
   pool.Foreach(func, vals);

   // Because a parallel snapshot does not preserve ordering of the events, we can just test
   // an aggregate such as the sum of all values.
   auto sum_ref = 0u;
   for (auto i = 0u; i < nevts; i++) sum_ref += i;
   for(auto& v : vals) {
      ROOT::RDataFrame df("tree", "dataframe_parallel_snapshots_" + std::to_string(v) + ".root");
      auto col = df.Take<int>("x");
      int sum_col = 0;
      for (auto i = 0u; i < nevts; i++)
         sum_col += col->at(i);
      EXPECT_EQ(sum_col, sum_ref);
   }
}

TEST(RDFConcurrency, ParallelRDFSnapshotsEnableThreadSafety)
{
   ROOT::EnableThreadSafety();
   ParallelRDFSnapshots();
}

TEST(RDFConcurrency, ParallelRDFSnapshotsEnableImplicitMT)
{
   ROOT::EnableImplicitMT(NUM_THREADS);
   ParallelRDFSnapshots();
   ROOT::DisableImplicitMT();
}

void ParallelRDFCaches()
{
   // Run the RDF construction and the event loop in parallel
   // and use an intermediate cache in the computation graph
   const auto nevts = 100u;
   auto func = [&] {
         std::this_thread::sleep_for(std::chrono::milliseconds(std::rand() / RAND_MAX * 200));
         ROOT::RDataFrame df(nevts);
         auto cache = df.Define("x", "(int)rdfentry_").Cache("x");
         return cache.Sum("x").GetValue();
      };

   ROOT::TThreadExecutor pool;
   auto res = pool.Map(func, 8);

   auto sum_ref = 0u;
   for (auto i = 0u; i < nevts; i++) sum_ref += i;
   for (auto i = 0u; i < res.size(); i++) {
      EXPECT_EQ(res[i], sum_ref);
   }
}

TEST(RDFConcurrency, ParallelRDFCachesEnableThreadSafety)
{
   ROOT::EnableThreadSafety();
   ParallelRDFCaches();
}

TEST(RDFConcurrency, ParallelRDFCachesEnableImplicitMT)
{
   ROOT::EnableImplicitMT(NUM_THREADS);
   ParallelRDFCaches();
   ROOT::DisableImplicitMT();
}

#endif
