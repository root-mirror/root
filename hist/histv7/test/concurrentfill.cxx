#include "gtest/gtest.h"
#include "ROOT/RHist.hxx"
#include "ROOT/RHistConcurrentFill.hxx"

#include <iostream>
#include <future>

using namespace ROOT;

using Filler_t = Experimental::RHistConcurrentFiller<Experimental::RH2D, 1024>;

// Functions for testing concurrentfill
void fillWithoutWeight(Filler_t filler)
{
   for (int i = 0; i < 3000; ++i) {
      filler.Fill({(double)i / 100, (double)i / 10});
   }
}

void fillWithWeights(Filler_t filler)
{
   for (int i = 0; i < 3000; ++i) {
      filler.Fill({(double)i / 100, (double)i / 10}, (float)i);
   }
}

void concurrentHistFillWithoutWeight(Experimental::RH2D &hist)
{
   Experimental::RHistConcurrentFillManager<Experimental::RH2D> fillMgr(hist);

   std::array<std::thread, 2> threads;

   for (auto &thr : threads) {
      thr = std::thread(fillWithoutWeight, fillMgr.MakeFiller());
   }

   for (auto &thr : threads)
      thr.join();
}

void concurrentHistFillWithWeights(Experimental::RH2D &hist)
{
   Experimental::RHistConcurrentFillManager<Experimental::RH2D> fillMgr(hist);

   std::array<std::thread, 2> threads;

   for (auto &thr : threads) {
      thr = std::thread(fillWithWeights, fillMgr.MakeFiller());
   }

   for (auto &thr : threads)
      thr.join();
}

// Test consistancy of the hist after concurrentfill
TEST(ConcurrentFillTest, HistConsistancy)
{
   Experimental::RH2D hist{{100, 0., 1.}, {{0., 1., 2., 3., 10.}}};

   concurrentHistFillWithoutWeight(hist);
   EXPECT_EQ(2 * 3000, hist.GetEntries());
   EXPECT_FLOAT_EQ(2.f, hist.GetBinContent({(double)42 / 100, (double)42 / 10}));

   concurrentHistFillWithWeights(hist);
   EXPECT_EQ(2 * (2 * 3000), hist.GetEntries());
   EXPECT_FLOAT_EQ(2 * (1.f + 42.f), hist.GetBinContent({(double)42 / 100, (double)42 / 10}));
}


// Test flush and fill manually
TEST(ConcurrentFillTest, FillFlush)
{
   Experimental::RH2D hist{{100, 0., 1.}, {{0., 1., 2., 3., 10.}}};

   Experimental::RHistConcurrentFillManager<Experimental::RH2D> fillMgr(hist);

   Filler_t Filler_1 = fillMgr.MakeFiller();
   Filler_t Filler_2 = fillMgr.MakeFiller();

   // Filling both buffers with data
   Filler_1.Fill({0.1111, 4.22});
   Filler_1.Fill({0.3333, 4.44}, .42f);

   Filler_2.Fill({0.2222, 4.11});
   Filler_2.Fill({0.4444, 4.33}, .32f);

   // Checking the consistancy of the first buffer
   EXPECT_EQ(2, (int)Filler_1.GetCoords().size());
   EXPECT_FLOAT_EQ(0.1111, Filler_1.GetCoords().front()[0]);
   EXPECT_FLOAT_EQ(4.22, Filler_1.GetCoords()[0][1]);

   EXPECT_FLOAT_EQ(0.3333, Filler_1.GetCoords().back()[0]);
   EXPECT_FLOAT_EQ(4.44, Filler_1.GetCoords()[1][1]);

   EXPECT_EQ(2, (int)Filler_1.GetWeights().size());
   EXPECT_FLOAT_EQ(1.f, Filler_1.GetWeights()[0]);
   EXPECT_FLOAT_EQ(.42f, Filler_1.GetWeights()[1]);

   // Checking the consistancy of the second buffer
   EXPECT_EQ(2, (int)Filler_2.GetCoords().size());
   EXPECT_FLOAT_EQ(0.2222, Filler_2.GetCoords().front()[0]);
   EXPECT_FLOAT_EQ(4.11, Filler_2.GetCoords()[0][1]);

   EXPECT_FLOAT_EQ(0.4444, Filler_2.GetCoords().back()[0]);
   EXPECT_FLOAT_EQ(4.33, Filler_2.GetCoords()[1][1]);

   EXPECT_EQ(2, (int)Filler_2.GetWeights().size());
   EXPECT_FLOAT_EQ(1.f, Filler_2.GetWeights()[0]);
   EXPECT_FLOAT_EQ(.32f, Filler_2.GetWeights()[1]);

   // Flushing one buffer, and checking consistancy of both buffers and of the hist
   Filler_1.Flush();
   EXPECT_EQ(2, hist.GetEntries());
   EXPECT_FLOAT_EQ(1.f, hist.GetBinContent({0.1111, 4.22}));
   EXPECT_FLOAT_EQ(.42f, hist.GetBinContent({0.3333, 4.44}));
   EXPECT_FLOAT_EQ(0.f, hist.GetBinContent({0.2222, 4.11}));
   EXPECT_FLOAT_EQ(0.f, hist.GetBinContent({0.4444, 4.33}));

   EXPECT_EQ(0, (int)Filler_1.GetWeights().size());
   EXPECT_EQ(2, (int)Filler_2.GetWeights().size());

   EXPECT_EQ(0, (int)Filler_1.GetCoords().size());
   EXPECT_EQ(2, (int)Filler_2.GetCoords().size());

   // Filling the first buffer with new data
   Filler_1.Fill({0.1111, 4.22}, .52f);
   Filler_1.Fill({0.3333, 4.44}, .32f);

   EXPECT_EQ(2, (int)Filler_1.GetWeights().size());

   // Flushing both buffers, and checking consistancy of both buffers and of the hist
   Filler_1.Flush();
   Filler_2.Flush();
   EXPECT_EQ(6, hist.GetEntries());
   EXPECT_FLOAT_EQ(1.f + .52f, hist.GetBinContent({0.1111, 4.22}));
   EXPECT_FLOAT_EQ(.42f + .32f, hist.GetBinContent({0.3333, 4.44}));
   EXPECT_FLOAT_EQ(1.f, hist.GetBinContent({0.2222, 4.11}));
   EXPECT_FLOAT_EQ(.32f, hist.GetBinContent({0.4444, 4.33}));

   EXPECT_EQ(0, (int)Filler_1.GetWeights().size());
   EXPECT_EQ(0, (int)Filler_2.GetWeights().size());

   EXPECT_EQ(0, (int)Filler_1.GetCoords().size());
   EXPECT_EQ(0, (int)Filler_2.GetCoords().size());
}
