#include "ROOT/TThreadedObject.hxx"
#include "TH1F.h"
#include "TRandom.h"

#include "gtest/gtest.h"

#include <condition_variable>
#include <mutex>
#include <thread>

void IsHistEqual(const TH1F &a, const TH1F &b)
{
   EXPECT_STREQ(a.GetName(), b.GetName()) << "The names of the histograms differ: " << a.GetName() << " " << b.GetName()
                                          << std::endl;

   EXPECT_STREQ(a.GetTitle(), b.GetTitle()) << "The title of the histograms differ: " << a.GetTitle() << " "
                                            << b.GetTitle() << std::endl;

   auto nbinsa = a.GetNbinsX();
   auto nbinsb = b.GetNbinsX();
   EXPECT_DOUBLE_EQ(nbinsa, nbinsb) << "The # of bins of the histograms differ: " << nbinsa << " " << nbinsb
                                    << std::endl;

   for (int i = 0; i < a.GetNbinsX(); ++i) {
      auto binca = a.GetBinContent(i);
      auto bincb = b.GetBinContent(i);
      EXPECT_DOUBLE_EQ(binca, bincb) << "The content of bin " << i << "  of the histograms differ: " << binca << " "
                                     << bincb << std::endl;

      auto binea = a.GetBinError(i);
      auto bineb = b.GetBinError(i);
      EXPECT_DOUBLE_EQ(binea, bineb) << "The error of bin " << i << "  of the histograms differ: " << binea << " "
                                     << bineb << std::endl;
   }
}

TEST(TThreadedObject, CreateAndDestroy)
{
   ROOT::TThreadedObject<TH1F> tto("h", "h", 64, -4, 4);
}

TEST(TThreadedObject, Get)
{
   TH1F model("h", "h", 64, -4, 4);
   ROOT::TThreadedObject<TH1F> tto("h", "h", 64, -4, 4);
   auto h = tto.Get();
   IsHistEqual(model, *h);
}

TEST(TThreadedObject, GetAtSlot)
{
   TH1F model("h", "h", 64, -4, 4);
   ROOT::TThreadedObject<TH1F> tto("h", "h", 64, -4, 4);
   auto h = tto.GetAtSlot(0);
   IsHistEqual(model, *h);
}

TEST(TThreadedObject, GetAtSlotUnchecked)
{
   TH1F model("h", "h", 64, -4, 4);
   ROOT::TThreadedObject<TH1F> tto("h", "h", 64, -4, 4);
   tto->SetName("h");
   auto h = tto.GetAtSlot(0);
   IsHistEqual(model, *h);
}

TEST(TThreadedObject, GetAtSlotRaw)
{
   TH1F model("h", "h", 64, -4, 4);
   ROOT::TThreadedObject<TH1F> tto("h", "h", 64, -4, 4);
   tto->SetName("h");
   auto h = tto.GetAtSlotRaw(0);
   IsHistEqual(model, *h);
}

TEST(TThreadedObject, SetAtSlot)
{
   ROOT::TThreadedObject<TH1F> tto("h", "h", 64, -4, 4);
   tto.SetAtSlot(1, std::make_shared<TH1F>("h", "h", 64, -4, 4));
   auto h0 = tto.GetAtSlot(0);
   auto h1 = tto.GetAtSlot(1);
   IsHistEqual(*h0, *h1);
}

TEST(TThreadedObject, Merge)
{
   TH1::AddDirectory(false);

   TH1F m0("h", "h", 64, -4, 4);
   TH1F m1("h", "h", 64, -4, 4);
   gRandom->SetSeed(1);
   m0.FillRandom("gaus");
   m1.FillRandom("gaus");
   m0.Add(&m1);

   ROOT::TThreadedObject<TH1F> tto("h", "h", 64, -4, 4);
   tto->SetName("h");
   tto.SetAtSlot(1, std::make_shared<TH1F>("h", "h", 64, -4, 4));
   gRandom->SetSeed(1);
   tto->FillRandom("gaus");
   tto.GetAtSlot(1)->FillRandom("gaus");
   auto hsum = tto.Merge();
   IsHistEqual(*hsum, m0);
}

TEST(TThreadedObject, SnapshotMerge)
{
   TH1::AddDirectory(false);

   TH1F m0("h", "h", 64, -4, 4);
   TH1F m1("h", "h", 64, -4, 4);
   gRandom->SetSeed(1);
   m0.FillRandom("gaus", 100);
   m1.FillRandom("gaus", 100);
   m0.Add(&m1);

   ROOT::TThreadedObject<TH1F> tto("h", "h", 64, -4, 4);
   tto->SetName("h");
   tto.SetAtSlot(1, std::make_shared<TH1F>("h", "h", 64, -4, 4));
   gRandom->SetSeed(1);
   tto->FillRandom("gaus", 100);
   tto.GetAtSlot(1)->FillRandom("gaus", 100);
   auto hsum0 = tto.SnapshotMerge();
   IsHistEqual(*hsum0, m0);
   auto hsum1 = tto.SnapshotMerge();
   IsHistEqual(*hsum1, m0);
   IsHistEqual(*hsum1, *hsum0);
   EXPECT_TRUE(hsum1 != hsum0);
}

TEST(TThreadedObject, GrowSlots)
{
   // create a TThreadedObject with 3 slots...
   ROOT::TThreadedObject<int> tto(ROOT::TNumSlots{3}, 42);

   // and then use it from 4 threads
   auto task = [&tto] { *tto.Get() = 1; };
   std::vector<std::thread> threads;
   for (int i = 0; i < 4; ++i)
      threads.emplace_back(task);
   for (auto &t : threads)
      t.join();

   auto sum_ints = [](std::shared_ptr<int> first, std::vector<std::shared_ptr<int>> &all) {
      for (auto &e : all)
         if (e != first)
            *first += *e;
   };
   EXPECT_EQ(*tto.Merge(sum_ints), 4);
}

TEST(TThreadedObject, GetNSlots)
{
   // default ctor produces fgMaxSlots slots
   EXPECT_EQ(ROOT::TThreadedObject<int>(42).GetNSlots(), ROOT::TThreadedObject<int>::fgMaxSlots.fVal);

   ROOT::TThreadedObject<int> tto(ROOT::TNumSlots{0});
   EXPECT_EQ(tto.GetNSlots(), 0u);

   // we need to make sure all threads will be in flight at the same time so the OS cannot re-use thread IDs,
   // otherwise TThreadedObject could count less than 4 threads.
   std::mutex m;
   std::condition_variable cv;
   bool ready = false;
   auto task = [&] {
      std::unique_lock<std::mutex> lk(m);
      cv.wait(lk, [&] { return ready; });
      tto.Get();
   };
   std::vector<std::thread> threads;
   for (int i = 0; i < 4; ++i)
      threads.emplace_back(task);
   {
      std::lock_guard<std::mutex> lk(m);
      ready = true;
   }
   cv.notify_all();
   for (auto &t : threads)
      t.join();

   EXPECT_EQ(tto.GetNSlots(), 4u);
}
