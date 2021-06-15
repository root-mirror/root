#include "ntuple_test.hxx"

TEST(RNTuple, Basics)
{
   FileRaii fileGuard("test_ntuple_barefile.ntuple");

   auto model = RNTupleModel::Create();
   auto wrPt = model->MakeField<float>("pt", 42.0);

   {
      RNTupleWriteOptions options;
      options.SetContainerFormat(ENTupleContainerFormat::kBare);
      auto ntuple = RNTupleWriter::Recreate(std::move(model), "f", fileGuard.GetPath(), options);
      ntuple->Fill();
      ntuple->CommitCluster();
      *wrPt = 24.0;
      ntuple->Fill();
      *wrPt = 12.0;
      ntuple->Fill();
   }

   auto ntuple = RNTupleReader::Open("f", fileGuard.GetPath());
   EXPECT_EQ(3U, ntuple->GetNEntries());
   auto rdPt = ntuple->GetModel()->GetDefaultEntry()->Get<float>("pt");

   ntuple->LoadEntry(0);
   EXPECT_EQ(42.0, *rdPt);
   ntuple->LoadEntry(1);
   EXPECT_EQ(24.0, *rdPt);
   ntuple->LoadEntry(2);
   EXPECT_EQ(12.0, *rdPt);
}

TEST(RNTuple, Extended)
{
   FileRaii fileGuard("test_ntuple_barefile_ext.ntuple");

   auto model = RNTupleModel::Create();
   auto wrVector = model->MakeField<std::vector<double>>("vector");

   TRandom3 rnd(42);
   double chksumWrite = 0.0;
   {
      RNTupleWriteOptions options;
      options.SetContainerFormat(ENTupleContainerFormat::kBare);
      auto ntuple = RNTupleWriter::Recreate(std::move(model), "f", fileGuard.GetPath(), options);
      constexpr unsigned int nEvents = 32000;
      for (unsigned int i = 0; i < nEvents; ++i) {
         auto nVec = 1 + floor(rnd.Rndm() * 1000.);
         wrVector->resize(nVec);
         for (unsigned int n = 0; n < nVec; ++n) {
            auto val = 1 + rnd.Rndm()*1000. - 500.;
            (*wrVector)[n] = val;
            chksumWrite += val;
         }
         ntuple->Fill();
         if (i % 1000 == 0)
            ntuple->CommitCluster();
      }
   }

   auto ntuple = RNTupleReader::Open("f", fileGuard.GetPath());
   auto rdVector = ntuple->GetModel()->GetDefaultEntry()->Get<std::vector<double>>("vector");

   double chksumRead = 0.0;
   for (auto entryId : *ntuple) {
      ntuple->LoadEntry(entryId);
      for (auto v : *rdVector)
         chksumRead += v;
   }
   EXPECT_EQ(chksumRead, chksumWrite);
}

TEST(RPageSinkBuf, Basics)
{
   struct TestModel {
      std::unique_ptr<RNTupleModel> fModel;
      std::shared_ptr<float> fFloatField;
      std::shared_ptr<std::vector<CustomStruct>> fFieldKlassVec;
      TestModel() {
         fModel = RNTupleModel::Create();
         fFloatField = fModel->MakeField<float>("pt");
         fFieldKlassVec = fModel->MakeField<std::vector<CustomStruct>>("klassVec");
      }
   };

   FileRaii fileGuardBuf("test_ntuple_sinkbuf_basics_buf.root");
   FileRaii fileGuard("test_ntuple_sinkbuf_basics.root");
   {
      TestModel bufModel;
      // PageSinkBuf wraps a concrete page source
      auto ntupleBuf = std::make_unique<RNTupleWriter>(std::move(bufModel.fModel),
         std::make_unique<RPageSinkBuf>(std::make_unique<RPageSinkFile>(
            "buf", fileGuardBuf.GetPath(), RNTupleWriteOptions()
      )));
      ntupleBuf->EnableMetrics();

      TestModel unbufModel;
      auto ntuple = std::make_unique<RNTupleWriter>(std::move(unbufModel.fModel),
         std::make_unique<RPageSinkFile>("unbuf", fileGuard.GetPath(), RNTupleWriteOptions()
      ));

      for (int i = 0; i < 20000; i++) {
         *bufModel.fFloatField = static_cast<float>(i);
         *unbufModel.fFloatField = static_cast<float>(i);
         CustomStruct klass;
         klass.a = 42.0;
         klass.v1.emplace_back(static_cast<float>(i));
         klass.v2.emplace_back(std::vector<float>(3, static_cast<float>(i)));
         klass.s = "hi" + std::to_string(i);
         *bufModel.fFieldKlassVec = std::vector<CustomStruct>{klass};
         *unbufModel.fFieldKlassVec = std::vector<CustomStruct>{klass};

         ntupleBuf->Fill();
         ntuple->Fill();

         if (i && i % 15000 == 0) {
            ntupleBuf->CommitCluster();
            ntuple->CommitCluster();
            auto *parallel_zip = ntupleBuf->GetMetrics().GetCounter(
               "RNTupleWriter.RPageSinkBuf.ParallelZip");
            ASSERT_FALSE(parallel_zip == nullptr);
            EXPECT_EQ(0, parallel_zip->GetValueAsInt());
         }
      }
   }

   auto ntupleBuf = RNTupleReader::Open("buf", fileGuardBuf.GetPath());
   auto ntuple = RNTupleReader::Open("unbuf", fileGuard.GetPath());
   EXPECT_EQ(ntuple->GetNEntries(), ntupleBuf->GetNEntries());

   auto viewPtBuf = ntupleBuf->GetView<float>("pt");
   auto viewKlassVecBuf = ntupleBuf->GetView<std::vector<CustomStruct>>("klassVec");
   auto viewPt = ntuple->GetView<float>("pt");
   auto viewKlassVec = ntuple->GetView<std::vector<CustomStruct>>("klassVec");
   for (auto i : ntupleBuf->GetEntryRange()) {
      EXPECT_EQ(static_cast<float>(i), viewPtBuf(i));
      EXPECT_EQ(viewPt(i), viewPtBuf(i));
      EXPECT_EQ(viewKlassVec(i).at(0).v1, viewKlassVecBuf(i).at(0).v1);
      EXPECT_EQ(viewKlassVec(i).at(0).v2, viewKlassVecBuf(i).at(0).v2);
      EXPECT_EQ(viewKlassVec(i).at(0).s, viewKlassVecBuf(i).at(0).s);
   }

   std::vector<std::pair<DescriptorId_t, std::int64_t>> pagePositions;
   std::size_t num_columns = 10;
   const auto &cluster0 = ntupleBuf->GetDescriptor().GetClusterDescriptor(0);
   for (std::size_t i = 0; i < num_columns; i++) {
      const auto &columnPages = cluster0.GetPageRange(i);
      for (const auto &page: columnPages.fPageInfos) {
         pagePositions.push_back(std::make_pair(i, page.fLocator.fPosition));
      }
   }

   auto sortedPages = pagePositions;
   std::sort(begin(sortedPages), end(sortedPages),
      [](const auto &a, const auto &b) { return a.second < b.second; });

   // For this test, ensure at least some columns have multiple pages
   ASSERT_TRUE(sortedPages.size() > num_columns);
   // Buffered sink cluster column pages are written out together
   for (std::size_t i = 0; i < pagePositions.size() - 1; i++) {
      // if the next page belongs to another column, skip the check
      if (pagePositions.at(i+1).first != pagePositions.at(i).first) {
         continue;
      }
      auto page = std::find(begin(sortedPages), end(sortedPages), pagePositions[i]);
      ASSERT_TRUE(page != sortedPages.end());
      auto next_page = page + 1;
      auto column = pagePositions[i].first;
      ASSERT_EQ(column, next_page->first);
   }
}

TEST(RPageSinkBuf, ParallelZip) {
   ROOT::EnableImplicitMT();

   FileRaii fileGuard("test_ntuple_sinkbuf_pzip.root");
   {
      auto model = RNTupleModel::Create();
      auto floatField = model->MakeField<float>("pt");
      auto fieldKlassVec = model->MakeField<std::vector<CustomStruct>>("klassVec");
      auto ntuple = std::make_unique<RNTupleWriter>(std::move(model),
         std::make_unique<RPageSinkBuf>(std::make_unique<RPageSinkFile>(
            "buf_pzip", fileGuard.GetPath(), RNTupleWriteOptions()
      )));
      ntuple->EnableMetrics();
      for (int i = 0; i < 20000; i++) {
         *floatField = static_cast<float>(i);
         CustomStruct klass;
         klass.a = 42.0;
         klass.v1.emplace_back(static_cast<float>(i));
         klass.v2.emplace_back(std::vector<float>(3, static_cast<float>(i)));
         klass.s = "hi" + std::to_string(i);
         *fieldKlassVec = std::vector<CustomStruct>{klass};
         ntuple->Fill();
         if (i && i % 15000 == 0) {
            ntuple->CommitCluster();
            auto *parallel_zip = ntuple->GetMetrics().GetCounter(
               "RNTupleWriter.RPageSinkBuf.ParallelZip");
            ASSERT_FALSE(parallel_zip == nullptr);
            EXPECT_EQ(1, parallel_zip->GetValueAsInt());
         }
      }
   }

   auto ntuple = RNTupleReader::Open("buf_pzip", fileGuard.GetPath());
   EXPECT_EQ(20000, ntuple->GetNEntries());

   auto viewPt = ntuple->GetView<float>("pt");
   auto viewKlassVec = ntuple->GetView<std::vector<CustomStruct>>("klassVec");
   for (auto i : ntuple->GetEntryRange()) {
      float fi = static_cast<float>(i);
      EXPECT_EQ(fi, viewPt(i));
      EXPECT_EQ(std::vector<float>{fi}, viewKlassVec(i).at(0).v1);
      EXPECT_EQ((std::vector<float>(3, fi)), viewKlassVec(i).at(0).v2.at(0));
      EXPECT_EQ("hi" + std::to_string(i), viewKlassVec(i).at(0).s);
   }
   ROOT::DisableImplicitMT();
}

TEST(RPageSinkBuf, SelectiveMT) {
   FileRaii fileGuard("test_ntuple_selective_mt.root");
   {
      auto model = RNTupleModel::Create();
      auto floatField = model->MakeField<float>("pt");
      auto fieldKlassVec = model->MakeField<std::vector<CustomStruct>>("klassVec");
      auto ntuple = std::make_unique<RNTupleWriter>(std::move(model),
         std::make_unique<RPageSinkBuf>(std::make_unique<RPageSinkFile>(
            "buf_pzip", fileGuard.GetPath(), RNTupleWriteOptions()
      )));
      ASSERT_FALSE(ROOT::IsImplicitMTEnabled());
      ntuple->EnableMT(); // only enable MT for the RNTupleWriter
      ASSERT_FALSE(ROOT::IsImplicitMTEnabled());
      ntuple->EnableMetrics();
      for (int i = 0; i < 20000; i++) {
         *floatField = static_cast<float>(i);
         CustomStruct klass;
         klass.a = 42.0;
         klass.v1.emplace_back(static_cast<float>(i));
         klass.v2.emplace_back(std::vector<float>(3, static_cast<float>(i)));
         klass.s = "hi" + std::to_string(i);
         *fieldKlassVec = std::vector<CustomStruct>{klass};
         ntuple->Fill();
         if (i && i % 15000 == 0) {
            ntuple->CommitCluster();
            auto *parallel_zip = ntuple->GetMetrics().GetCounter(
               "RNTupleWriter.RPageSinkBuf.ParallelZip");
            ASSERT_FALSE(parallel_zip == nullptr);
            EXPECT_EQ(1, parallel_zip->GetValueAsInt());
         }
      }
   }

   auto ntuple = RNTupleReader::Open("buf_pzip", fileGuard.GetPath());
   EXPECT_EQ(20000, ntuple->GetNEntries());

   auto viewPt = ntuple->GetView<float>("pt");
   auto viewKlassVec = ntuple->GetView<std::vector<CustomStruct>>("klassVec");
   for (auto i : ntuple->GetEntryRange()) {
      float fi = static_cast<float>(i);
      EXPECT_EQ(fi, viewPt(i));
      EXPECT_EQ(std::vector<float>{fi}, viewKlassVec(i).at(0).v1);
      EXPECT_EQ((std::vector<float>(3, fi)), viewKlassVec(i).at(0).v2.at(0));
      EXPECT_EQ("hi" + std::to_string(i), viewKlassVec(i).at(0).s);
   }
}
