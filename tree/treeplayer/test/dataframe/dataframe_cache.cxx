#include "ROOT/TDataFrame.hxx"
#include "TH1F.h"
#include "TRandom.h"

#include "gtest/gtest.h"

#include <algorithm>

using namespace ROOT::Experimental;

TEST(Cache, FundType)
{
   TDataFrame tdf(1);
   int i = 1;
   auto cached = tdf.Define("c0", [&i]() { return i++; })
                    .Define("c1", []() { return 1.; })
                    .Cache<int, double>({"c0","c1"});

   auto c = cached.Count();
   auto m = cached.Mean<int>("c0");

   EXPECT_EQ(1, *m);
   EXPECT_EQ(1UL, *c);
}


TEST(Cache, Contiguity)
{
   TDataFrame tdf(2);
   auto f = 0.f;
   auto cached = tdf.Define("float", [&f]() { return f++; })
                    .Cache<float>({"float"});
   int counter = 0;
   float *fPrec = nullptr;
   auto count = [&counter, &fPrec](float &ff){
      if ( 1 == counter++){
         EXPECT_EQ(1U, std::distance(fPrec, &ff));
      }
      fPrec = &ff;
   };
   cached.Foreach(count,{"float"});
}


TEST(Cache, Class)
{
   TH1F h("","h",64,0,1);
   gRandom->SetSeed(1);
   h.FillRandom("gaus",10);
   TDataFrame tdf(1);
   auto cached = tdf.Define("c0", [&h]() { return h; })
                    .Cache<TH1F>({"c0"});

   auto c = cached.Count();
   auto d = cached.Define("Mean", [](TH1F& hh){return hh.GetMean();}, {"c0"})
                  .Define("StdDev", [](TH1F& hh){return hh.GetStdDev();}, {"c0"});
   auto m = d.Max<double>("Mean");
   auto s = d.Max<double>("StdDev");

   EXPECT_EQ(h.GetMean(), *m);
   EXPECT_EQ(h.GetStdDev(), *s);
   EXPECT_EQ(1UL, *c);
}
