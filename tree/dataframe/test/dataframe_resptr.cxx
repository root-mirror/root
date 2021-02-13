#include "ROOTUnitTestSupport.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultHandle.hxx"

#include "gtest/gtest.h"

#include <vector>
#include <cstddef>

using namespace ROOT::RDF;

class Dummy {
};

TEST(RResultPtr, DefCtor)
{
   RResultPtr<Dummy> p1, p2;
   EXPECT_TRUE(p1 == nullptr);
   EXPECT_TRUE(nullptr == p1);
   EXPECT_TRUE(p1 == p2);
   EXPECT_TRUE(p2 == p1);
}

TEST(RResultPtr, CopyCtor)
{
   ROOT::RDataFrame d(1);
   auto hasRun = false;
   auto m = d.Define("i", [&hasRun]() {
                hasRun = true;
                return (int)1;
             }).Mean<int>("i");
   auto mc = m;
   EXPECT_TRUE(mc == m);
   EXPECT_TRUE(m == mc);

   EXPECT_FALSE(hasRun);

   EXPECT_EQ(*mc, *m);

   EXPECT_TRUE(hasRun);
}

TEST(RResultPtr, MoveCtor)
{
   ROOT::RDataFrame df(1);
   ROOT::RDF::RResultPtr<ULong64_t> res(df.Count());

   // also test move-assignment
   res = df.Count();

   EXPECT_EQ(*res, 1u);
}

TEST(RResultPtr, ImplConv)
{
   RResultPtr<Dummy> p1;
   EXPECT_FALSE(p1);

   ROOT::RDataFrame d(1);
   auto hasRun = false;
   auto m = d.Define("i", [&hasRun]() {
                hasRun = true;
                return (int)1;
             }).Histo1D<int>("i");

   EXPECT_TRUE(m != nullptr);
   EXPECT_FALSE(hasRun);

   *m;

   EXPECT_TRUE(m != nullptr);
   EXPECT_TRUE(hasRun);
}

TEST(RResultPtr, IsReady)
{
   ROOT::RDataFrame df(1);
   auto p = df.Sum<ULong64_t>("rdfentry_");
   EXPECT_FALSE(p.IsReady());

   p.GetValue();
   EXPECT_TRUE(p.IsReady());
}

TEST(RResultHandle, Ctor)
{
   ROOT::RDataFrame df(3);
   auto df2 = df.Define("x", [] { return 1.f; }, {});
   auto r1 = df2.Histo1D<float>("x");
   auto r2 = df2.Sum<float>("x");
   auto r3 = df2.Count();

   // Test constructors
   RResultHandle h1(r1);
   RResultHandle h2(r2);
   RResultHandle h3(r3);
}

TEST(RResultHandle, StdVector)
{
   ROOT::RDataFrame df(3);
   auto df2 = df.Define("x", [] { return 1.f; }, {});
   auto r1 = df2.Histo1D<float>("x");
   auto r2 = df2.Sum<float>("x");
   auto r3 = df2.Count();

   // Test push_back and emplace_back
   std::vector<RResultHandle> v1;
   v1.emplace_back(r1);
   v1.push_back(r1);

   // Test initializer list
   std::vector<RResultHandle> v2 = {r1, r2, r3};
}

TEST(RResultHandle, TriggerRun)
{
   ROOT::RDataFrame df(3);
   auto r1 = df.Sum<ULong64_t>("rdfentry_");
   auto r2 = df.Count();

   // Check that no event loops has run
   EXPECT_EQ(df.GetNRuns(), 0u);

   // Create handle and check that no event loop is triggered
   RResultHandle h(r1);
   EXPECT_EQ(df.GetNRuns(), 0u);

   // Trigger event loop
   h.GetValue<ULong64_t>();
   EXPECT_EQ(df.GetNRuns(), 1u);

   // Trigger event loop again via the RResultPtr and check that loop has not run
   r1.GetPtr();
   r1.GetValue();
   EXPECT_EQ(df.GetNRuns(), 1u);

   // Trigger event loop again via the second RResultPtr and check that loop has not run
   r2.GetValue();
   EXPECT_EQ(df.GetNRuns(), 1u);
}

TEST(RResultHandle, IsReady)
{
   ROOT::RDataFrame df(3);
   auto r1 = df.Sum<ULong64_t>("rdfentry_");
   auto r2 = df.Count();

   RResultHandle h1(r1);
   RResultHandle h2(r2);

   // Event loop has not yet run
   EXPECT_FALSE(h1.IsReady());
   EXPECT_FALSE(h2.IsReady());

   // Trigger loop and check that event loop is marked as run in both actions
   h1.GetValue<ULong64_t>();
   EXPECT_TRUE(h1.IsReady());
   EXPECT_TRUE(h2.IsReady());
}

TEST(RResultHandle, GetPtrOrValue)
{
   ROOT::RDataFrame df(3);
   auto df2 = df.Define("x", [] { return 1.f; }, {});
   auto r1 = df2.Sum<float>("x");
   auto r2 = df2.Count();

   RResultHandle h1(r1);
   RResultHandle h2(r2);

   EXPECT_EQ(h1.GetPtr<float>(), r1.GetPtr());
   EXPECT_EQ(h1.GetValue<float>(), r1.GetValue());

   EXPECT_EQ(h2.GetPtr<ULong64_t>(), r2.GetPtr());
   EXPECT_EQ(h2.GetValue<ULong64_t>(), r2.GetValue());
}

TEST(RResultHandle, Comparison)
{
   ROOT::RDataFrame df(3);
   auto r1 = df.Sum<ULong64_t>("rdfentry_");
   auto r2 = df.Count();

   RResultHandle h1(r1);
   RResultHandle h2(r2);
   RResultHandle h3(r1);

   EXPECT_FALSE(h1 == h2);
   EXPECT_TRUE(h1 != h2);

   EXPECT_TRUE(h1 == h1);
   EXPECT_FALSE(h1 != h1);

   EXPECT_TRUE(h1 == h3);
   EXPECT_FALSE(h1 != h3);
}

TEST(RResultHandle, CheckType)
{
   ROOT::RDataFrame df(3);
   auto r = df.Sum<ULong64_t>("rdfentry_");
   RResultHandle h(r);
   EXPECT_THROW(h.GetValue<float>(), std::runtime_error);
   EXPECT_THROW(h.GetPtr<float>(), std::runtime_error);
}
