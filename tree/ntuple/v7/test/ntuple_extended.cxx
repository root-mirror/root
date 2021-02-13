#include <ROOT/RConfig.hxx>

#include "ntuple_test.hxx"

#include "TROOT.h"

TEST(RNTuple, RealWorld1)
{
   ROOT::EnableImplicitMT();
   FileRaii fileGuard("test_ntuple_realworld1.root");

   // See https://github.com/olifre/root-io-bench/blob/master/benchmark.cpp
   auto modelWrite = RNTupleModel::Create();
   auto wrEvent   = modelWrite->MakeField<std::uint32_t>("event");
   auto wrSignal  = modelWrite->MakeField<bool>("signal");
   auto wrEnergy  = modelWrite->MakeField<double>("energy");
   auto wrTimes   = modelWrite->MakeField<std::vector<double>>("times");
   auto wrIndices = modelWrite->MakeField<std::vector<std::uint32_t>>("indices");

   TRandom3 rnd(42);
   double chksumWrite = 0.0;
   {
      auto ntuple = RNTupleWriter::Recreate(std::move(modelWrite), "myNTuple", fileGuard.GetPath());
      constexpr unsigned int nEvents = 60000;
      for (unsigned int i = 0; i < nEvents; ++i) {
         *wrEvent = i;
         *wrEnergy = rnd.Rndm() * 1000.;
         *wrSignal = i % 2;

         chksumWrite += double(*wrEvent);
         chksumWrite += double(*wrSignal);
         chksumWrite += *wrEnergy;

         auto nTimes = 1 + floor(rnd.Rndm() * 1000.);
         wrTimes->resize(nTimes);
         for (unsigned int n = 0; n < nTimes; ++n) {
            wrTimes->at(n) = 1 + rnd.Rndm()*1000. - 500.;
            chksumWrite += wrTimes->at(n);
         }

         auto nIndices = 1 + floor(rnd.Rndm() * 1000.);
         wrIndices->resize(nIndices);
         for (unsigned int n = 0; n < nIndices; ++n) {
            wrIndices->at(n) = 1 + floor(rnd.Rndm() * 1000.);
            chksumWrite += double(wrIndices->at(n));
         }

         ntuple->Fill();
      }
   }

   auto modelRead  = RNTupleModel::Create();
   auto rdEvent   = modelRead->MakeField<std::uint32_t>("event");
   auto rdSignal  = modelRead->MakeField<bool>("signal");
   auto rdEnergy  = modelRead->MakeField<double>("energy");
   auto rdTimes   = modelRead->MakeField<std::vector<double>>("times");
   auto rdIndices = modelRead->MakeField<std::vector<std::uint32_t>>("indices");

   double chksumRead = 0.0;
   auto ntuple = RNTupleReader::Open(std::move(modelRead), "myNTuple", fileGuard.GetPath());
   for (auto entryId : *ntuple) {
      ntuple->LoadEntry(entryId);
      chksumRead += double(*rdEvent) + double(*rdSignal) + *rdEnergy;
      for (auto t : *rdTimes) chksumRead += t;
      for (auto ind : *rdIndices) chksumRead += double(ind);
   }

   // The floating point arithmetic should have been executed in the same order for reading and writing,
   // thus we expect the checksums to be bitwise identical
   EXPECT_EQ(chksumRead, chksumWrite);
}


// Stress test the asynchronous cluster pool by a deliberately unfavourable read pattern
TEST(RNTuple, RandomAccess)
{
   ROOT::EnableImplicitMT();
   FileRaii fileGuard("test_ntuple_random_access.root");

   auto modelWrite = RNTupleModel::Create();
   auto wrValue   = modelWrite->MakeField<std::int32_t>("value", 42);

   constexpr unsigned int nEvents = 1000000;
   {
      RNTupleWriteOptions options;
      options.SetCompression(0);
      auto ntuple = RNTupleWriter::Recreate(std::move(modelWrite), "myNTuple", fileGuard.GetPath(), options);
      for (unsigned int i = 0; i < nEvents; ++i)
         ntuple->Fill();
   }

   RNTupleReadOptions options;
   options.SetClusterCache(RNTupleReadOptions::EClusterCache::kOn);
   auto ntuple = RNTupleReader::Open("myNTuple", fileGuard.GetPath(), options);
   EXPECT_GT(ntuple->GetDescriptor().GetNClusters(), 10);

   auto viewValue = ntuple->GetView<std::int32_t>("value");

   std::int32_t sum = 0;
   constexpr unsigned int nSamples = 1000;
   TRandom3 rnd(42);
   for (unsigned int i = 0; i < 1000; ++i) {
      auto entryId = floor(rnd.Rndm() * (nEvents - 1));
      sum += viewValue(entryId);
   }
   EXPECT_EQ(42 * nSamples, sum);
}


#if !defined(_MSC_VER) || defined(R__ENABLE_BROKEN_WIN_TESTS)
TEST(RNTuple, LargeFile)
{
   ROOT::EnableImplicitMT();
   FileRaii fileGuard("test_large_file.root");

   auto modelWrite = RNTupleModel::Create();
   auto& wrEnergy  = *modelWrite->MakeField<double>("energy");

   TRandom3 rnd(42);
   double chksumWrite = 0.0;
   {
      RNTupleWriteOptions options;
      options.SetCompression(0);
      auto ntuple = RNTupleWriter::Recreate(std::move(modelWrite), "myNTuple", fileGuard.GetPath(), options);
      constexpr unsigned long nEvents = 1024 * 1024 * 256; // Exceed 2GB file size
      for (unsigned int i = 0; i < nEvents; ++i) {
         wrEnergy = rnd.Rndm();
         chksumWrite += wrEnergy;
         ntuple->Fill();
      }
   }
#ifdef R__SEEK64
   FILE *file = fopen64(fileGuard.GetPath().c_str(), "rb");
   ASSERT_TRUE(file != nullptr);
   EXPECT_EQ(0, fseeko64(file, 0, SEEK_END));
   EXPECT_GT(ftello64(file), 2048LL * 1024LL * 1024LL);
#else
   FILE *file = fopen(fileGuard.GetPath().c_str(), "rb");
   ASSERT_TRUE(file != nullptr);
   EXPECT_EQ(0, fseek(file, 0, SEEK_END));
   EXPECT_GT(ftell(file), 2048LL * 1024LL * 1024LL);
#endif
   fclose(file);

   auto ntuple = RNTupleReader::Open("myNTuple", fileGuard.GetPath());
   auto rdEnergy  = ntuple->GetView<double>("energy");
   double chksumRead = 0.0;

   for (auto i : ntuple->GetEntryRange()) {
      chksumRead += rdEnergy(i);
   }

   EXPECT_EQ(chksumRead, chksumWrite);
   auto f = TFile::Open(fileGuard.GetPath().c_str(), "READ");
   EXPECT_TRUE(f != nullptr);
   delete f;
}
#endif
