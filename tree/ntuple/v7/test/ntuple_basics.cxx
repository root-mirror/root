#include "ntuple_test.hxx"

#if __cplusplus >= 201703L
TEST(RNTuple, ReconstructModel)
{
   FileRaii fileGuard("test_ntuple_reconstruct.root");
   auto model = RNTupleModel::Create();
   auto fieldPt = model->MakeField<float>("pt", 42.0);
   auto fieldNnlo = model->MakeField<std::vector<std::vector<float>>>("nnlo");
   auto fieldKlass = model->MakeField<CustomStruct>("klass");
   auto fieldArray = model->MakeField<std::array<double, 2>>("array");
   auto fieldVariant = model->MakeField<std::variant<double, std::variant<std::string, double>>>("variant");
   {
      RPageSinkFile sink("myNTuple", fileGuard.GetPath(), RNTupleWriteOptions());
      sink.Create(*model.get());
      sink.CommitDataset();
      model = nullptr;
   }

   RPageSourceFile source("myNTuple", fileGuard.GetPath(), RNTupleReadOptions());
   source.Attach();

   auto modelReconstructed = source.GetDescriptor().GenerateModel();
   EXPECT_EQ(nullptr, modelReconstructed->GetDefaultEntry()->Get<float>("xyz"));
   auto vecPtr = modelReconstructed->GetDefaultEntry()->Get<std::vector<std::vector<float>>>("nnlo");
   EXPECT_TRUE(vecPtr != nullptr);
   // Don't crash
   vecPtr->push_back(std::vector<float>{1.0});
   auto array = modelReconstructed->GetDefaultEntry()->Get<std::array<double, 2>>("array");
   EXPECT_TRUE(array != nullptr);
   auto variant = modelReconstructed->GetDefaultEntry()->Get<
      std::variant<double, std::variant<std::string, double>>>("variant");
   EXPECT_TRUE(variant != nullptr);
}
#endif // __cplusplus >= 201703L

TEST(RNTuple, MultipleInFile)
{
   FileRaii fileGuard("test_ntuple_multi.root");
   auto file = TFile::Open(fileGuard.GetPath().c_str(), "RECREATE");
   {
      auto model = RNTupleModel::Create();
      auto fieldPt = model->MakeField<float>("pt", 42.0);
      auto ntuple = RNTupleWriter::Append(std::move(model), "first", *file);
      ntuple->Fill();
   }
   {
      auto model = RNTupleModel::Create();
      auto fieldPt = model->MakeField<float>("E", 1.0);
      auto ntuple = RNTupleWriter::Append(std::move(model), "second", *file);
      ntuple->Fill();
   }
   file->Close();
   delete file;

   auto ntupleFirst = RNTupleReader::Open("first", fileGuard.GetPath());
   auto viewPt = ntupleFirst->GetView<float>("pt");
   int n = 0;
   for (auto i : ntupleFirst->GetEntryRange()) {
      EXPECT_EQ(42.0, viewPt(i));
      n++;
   }
   EXPECT_EQ(1, n);

   auto ntupleSecond = RNTupleReader::Open("second", fileGuard.GetPath());
   auto viewE = ntupleSecond->GetView<float>("E");
   n = 0;
   for (auto i : ntupleSecond->GetEntryRange()) {
      EXPECT_EQ(1.0, viewE(i));
      n++;
   }
   EXPECT_EQ(1, n);
}


TEST(RNTuple, WriteRead)
{
   FileRaii fileGuard("test_ntuple_writeread.root");

   auto modelWrite = RNTupleModel::Create();
   auto wrSignal = modelWrite->MakeField<bool>("signal", true);
   auto wrPt = modelWrite->MakeField<float>("pt", 42.0);
   auto wrEnergy = modelWrite->MakeField<float>("energy", 7.0);
   auto wrTag = modelWrite->MakeField<std::string>("tag", "xyz");
   auto wrJets = modelWrite->MakeField<std::vector<float>>("jets");
   wrJets->push_back(1.0);
   wrJets->push_back(2.0);
   auto wrNnlo = modelWrite->MakeField<std::vector<std::vector<float>>>("nnlo");
   wrNnlo->push_back(std::vector<float>());
   wrNnlo->push_back(std::vector<float>{1.0});
   wrNnlo->push_back(std::vector<float>{1.0, 2.0, 4.0, 8.0});
   auto wrKlass = modelWrite->MakeField<CustomStruct>("klass");
   wrKlass->s = "abc";

   auto modelRead = modelWrite->Clone();

   {
      RNTupleWriter ntuple(std::move(modelWrite),
         std::make_unique<RPageSinkFile>("myNTuple", fileGuard.GetPath(), RNTupleWriteOptions()));
      ntuple.Fill();
   }

   auto rdSignal = modelRead->Get<bool>("signal");
   auto rdPt = modelRead->Get<float>("pt");
   auto rdEnergy = modelRead->Get<float>("energy");
   auto rdTag = modelRead->Get<std::string>("tag");
   auto rdJets = modelRead->Get<std::vector<float>>("jets");
   auto rdNnlo = modelRead->Get<std::vector<std::vector<float>>>("nnlo");
   auto rdKlass = modelRead->Get<CustomStruct>("klass");

   RNTupleReader ntuple(std::move(modelRead),
      std::make_unique<RPageSourceFile>("myNTuple", fileGuard.GetPath(), RNTupleReadOptions()));
   EXPECT_EQ(1U, ntuple.GetNEntries());
   ntuple.LoadEntry(0);

   EXPECT_TRUE(*rdSignal);
   EXPECT_EQ(42.0, *rdPt);
   EXPECT_EQ(7.0, *rdEnergy);
   EXPECT_STREQ("xyz", rdTag->c_str());

   EXPECT_EQ(2U, rdJets->size());
   EXPECT_EQ(1.0, (*rdJets)[0]);
   EXPECT_EQ(2.0, (*rdJets)[1]);

   EXPECT_EQ(3U, rdNnlo->size());
   EXPECT_EQ(0U, (*rdNnlo)[0].size());
   EXPECT_EQ(1U, (*rdNnlo)[1].size());
   EXPECT_EQ(4U, (*rdNnlo)[2].size());
   EXPECT_EQ(1.0, (*rdNnlo)[1][0]);
   EXPECT_EQ(1.0, (*rdNnlo)[2][0]);
   EXPECT_EQ(2.0, (*rdNnlo)[2][1]);
   EXPECT_EQ(4.0, (*rdNnlo)[2][2]);
   EXPECT_EQ(8.0, (*rdNnlo)[2][3]);

   EXPECT_STREQ("abc", rdKlass->s.c_str());
}

TEST(RNTuple, Clusters)
{
   FileRaii fileGuard("test_ntuple_clusters.root");

   auto modelWrite = RNTupleModel::Create();
   auto wrPt = modelWrite->MakeField<float>("pt", 42.0);
   auto wrTag = modelWrite->MakeField<std::string>("tag", "xyz");
   auto wrNnlo = modelWrite->MakeField<std::vector<std::vector<float>>>("nnlo");
   auto wrFourVec = modelWrite->MakeField<std::array<float, 4>>("fourVec");
   wrNnlo->push_back(std::vector<float>());
   wrNnlo->push_back(std::vector<float>{1.0});
   wrNnlo->push_back(std::vector<float>{1.0, 2.0, 4.0, 8.0});
   wrFourVec->at(0) = 0.0;
   wrFourVec->at(1) = 1.0;
   wrFourVec->at(2) = 2.0;
   wrFourVec->at(3) = 3.0;

   auto modelRead = modelWrite->Clone();

   {
      RNTupleWriter ntuple(std::move(modelWrite),
         std::make_unique<RPageSinkFile>("myNTuple", fileGuard.GetPath(), RNTupleWriteOptions()));
      ntuple.Fill();
      ntuple.CommitCluster();
      *wrPt = 24.0;
      wrNnlo->clear();
      *wrTag = "";
      wrFourVec->at(2) = 42.0;
      ntuple.Fill();
      *wrPt = 12.0;
      wrNnlo->push_back(std::vector<float>{42.0});
      *wrTag = "12345";
      wrFourVec->at(1) = 24.0;
      ntuple.Fill();
   }

   auto rdPt = modelRead->Get<float>("pt");
   auto rdTag = modelRead->Get<std::string>("tag");
   auto rdNnlo = modelRead->Get<std::vector<std::vector<float>>>("nnlo");
   auto rdFourVec = modelRead->Get<std::array<float, 4>>("fourVec");

   RNTupleReader ntuple(std::move(modelRead),
      std::make_unique<RPageSourceFile>("myNTuple", fileGuard.GetPath(), RNTupleReadOptions()));
   EXPECT_EQ(3U, ntuple.GetNEntries());

   ntuple.LoadEntry(0);
   EXPECT_EQ(42.0, *rdPt);
   EXPECT_STREQ("xyz", rdTag->c_str());
   EXPECT_EQ(3U, rdNnlo->size());
   EXPECT_EQ(0U, (*rdNnlo)[0].size());
   EXPECT_EQ(1U, (*rdNnlo)[1].size());
   EXPECT_EQ(4U, (*rdNnlo)[2].size());
   EXPECT_EQ(1.0, (*rdNnlo)[1][0]);
   EXPECT_EQ(1.0, (*rdNnlo)[2][0]);
   EXPECT_EQ(2.0, (*rdNnlo)[2][1]);
   EXPECT_EQ(4.0, (*rdNnlo)[2][2]);
   EXPECT_EQ(8.0, (*rdNnlo)[2][3]);
   EXPECT_EQ(0.0, (*rdFourVec)[0]);
   EXPECT_EQ(1.0, (*rdFourVec)[1]);
   EXPECT_EQ(2.0, (*rdFourVec)[2]);
   EXPECT_EQ(3.0, (*rdFourVec)[3]);

   ntuple.LoadEntry(1);
   EXPECT_EQ(24.0, *rdPt);
   EXPECT_STREQ("", rdTag->c_str());
   EXPECT_TRUE(rdNnlo->empty());
   EXPECT_EQ(42.0, (*rdFourVec)[2]);

   ntuple.LoadEntry(2);
   EXPECT_EQ(12.0, *rdPt);
   EXPECT_STREQ("12345", rdTag->c_str());
   EXPECT_EQ(1U, rdNnlo->size());
   EXPECT_EQ(1U, (*rdNnlo)[0].size());
   EXPECT_EQ(42.0, (*rdNnlo)[0][0]);
   EXPECT_EQ(24.0, (*rdFourVec)[1]);
}

TEST(RNTuple, ClusterEntries)
{
   FileRaii fileGuard("test_ntuple_cluster_entries.root");
   auto model = RNTupleModel::Create();
   auto field = model->MakeField<float>({"pt", "transverse momentum"}, 42.0);

   {
      RNTupleWriteOptions opt;
      opt.SetNEntriesPerCluster(5);
      auto ntuple = RNTupleWriter::Recreate(
         std::move(model), "ntuple", fileGuard.GetPath(), opt
      );
      for (int i = 0; i < 100; i++) {
         ntuple->Fill();
      }
   }

   auto ntuple = RNTupleReader::Open("ntuple", fileGuard.GetPath());
   // 100 entries / 5 entries per cluster
   EXPECT_EQ(20, ntuple->GetDescriptor().GetNClusters());
}

TEST(RNTuple, ElementsPerPage)
{
   FileRaii fileGuard("test_ntuple_elements_per_page.root");
   auto model = RNTupleModel::Create();
   auto field = model->MakeField<float>({"pt", "transverse momentum"}, 42.0);

   {
      RNTupleWriteOptions opt;
      opt.SetNElementsPerPage(5);
      auto ntuple = RNTupleWriter::Recreate(
         std::move(model), "ntuple", fileGuard.GetPath(), opt
      );
      for (int i = 0; i < 100; i++) {
         ntuple->Fill();
      }
   }

   auto ntuple = RNTupleReader::Open("ntuple", fileGuard.GetPath());
   const auto &col0_pages = ntuple->GetDescriptor().GetClusterDescriptor(0).GetPageRange(0);
   // 100 column elements / 5 elements per page
   EXPECT_EQ(20, col0_pages.fPageInfos.size());
}

TEST(RNTupleModel, EnforceValidFieldNames)
{
   auto model = RNTupleModel::Create();

   auto field = model->MakeField<float>("pt", 42.0);

   // MakeField
   try {
      auto field2 = model->MakeField<float>("pt", 42.0);
      FAIL() << "repeated field names should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("field name 'pt' already exists"));
   }
   try {
      auto field3 = model->MakeField<float>("", 42.0);
      FAIL() << "empty string as field name should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("name cannot be empty string"));
   }
   try {
      auto field3 = model->MakeField<float>("pt.pt", 42.0);
      FAIL() << "field name with periods should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("name 'pt.pt' cannot contain dot characters '.'"));
   }

   // AddField
   try {
      model->AddField(std::make_unique<RField<float>>(RField<float>("pt")));
      FAIL() << "repeated field names should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("field name 'pt' already exists"));
   }
   try {
      float num = 10.0;
      model->AddField("pt", &num);
      FAIL() << "repeated field names should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("field name 'pt' already exists"));
   }

   // MakeCollection
   try {
      auto otherModel = RNTupleModel::Create();
      auto collection = model->MakeCollection("pt", std::move(otherModel));
      FAIL() << "repeated field names should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("field name 'pt' already exists"));
   }

   // subfield names don't throw because full name differs (otherModel.pt)
   auto otherModel = RNTupleModel::Create();
   auto otherField = otherModel->MakeField<float>("pt", 42.0);
   auto collection = model->MakeCollection("otherModel", std::move(otherModel));
}

TEST(RNTupleModel, FieldDescriptions)
{
   FileRaii fileGuard("test_ntuple_field_descriptions.root");
   auto model = RNTupleModel::Create();

   auto pt = model->MakeField<float>({"pt", "transverse momentum"}, 42.0);

   float num = 10.0;
   model->AddField({"mass", "mass"}, &num);

   auto charge = std::make_unique<RField<float>>(RField<float>("charge"));
   charge->SetDescription("electric charge");
   model->AddField(std::move(charge));

   {
      RNTupleWriter::Recreate(std::move(model), "ntuple", fileGuard.GetPath());
   }

   auto ntuple = RNTupleReader::Open("ntuple", fileGuard.GetPath());
   std::vector<std::string> fieldDescriptions;
   for (auto& f: ntuple->GetDescriptor().GetTopLevelFields()) {
      fieldDescriptions.push_back(f.GetFieldDescription());
   }
   ASSERT_EQ(3, fieldDescriptions.size());
   EXPECT_EQ(std::string("transverse momentum"), fieldDescriptions[0]);
   EXPECT_EQ(std::string("mass"), fieldDescriptions[1]);
   EXPECT_EQ(std::string("electric charge"), fieldDescriptions[2]);
}

TEST(RNTupleModel, CollectionFieldDescriptions)
{
   FileRaii fileGuard("test_ntuple_collection_field_descriptions.root");
   {
      auto muon = RNTupleModel::Create();
      muon->SetDescription("muons after basic selection");

      auto model = RNTupleModel::Create();
      model->MakeCollection("Muon", std::move(muon));
      RNTupleWriter::Recreate(std::move(model), "ntuple", fileGuard.GetPath());
   }

   auto ntuple = RNTupleReader::Open("ntuple", fileGuard.GetPath());
   const auto& muon_desc = *ntuple->GetDescriptor().GetTopLevelFields().begin();
   EXPECT_EQ(std::string("muons after basic selection"), muon_desc.GetFieldDescription());
}

TEST(RNTupleModel, GetField)
{
   auto m = RNTupleModel::Create();
   m->MakeField<int>("x");
   m->MakeField<CustomStruct>("cs");
   EXPECT_EQ(m->GetField("x")->GetName(), "x");
   EXPECT_EQ(m->GetField("x")->GetType(), "std::int32_t");
   EXPECT_EQ(m->GetField("cs.v1")->GetName(), "v1");
   EXPECT_EQ(m->GetField("cs.v1")->GetType(), "std::vector<float>");
   EXPECT_EQ(m->GetField("nonexistent"), nullptr);
}

TEST(RNTuple, EmptyString)
{
   // empty storage string
   try {
      auto model = RNTupleModel::Create();
      auto ntuple = RNTupleWriter::Recreate(std::move(model), "myNTuple", "");
      FAIL() << "empty writer storage location should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("empty storage location"));
   }
   try {
      auto ntuple = RNTupleReader::Open("myNTuple", "");
      FAIL() << "empty reader storage location should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("empty storage location"));
   }

   // empty RNTuple name
   try {
      auto model = RNTupleModel::Create();
      auto ntuple = RNTupleWriter::Recreate(std::move(model), "", "file.root");
      FAIL() << "empty RNTuple name should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("empty RNTuple name"));
   }
   try {
      auto ntuple = RNTupleReader::Open("", "file.root");
      FAIL() << "empty RNTuple name should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("empty RNTuple name"));
   }
}

TEST(RNTuple, NullSafety)
{
   // RNTupleModel
   auto model = RNTupleModel::Create();
   try {
      auto bad = model->MakeCollection("bad", nullptr);
      FAIL() << "null submodels should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("null collectionModel"));
   }
   try {
      model->AddField(nullptr);
      FAIL() << "null fields should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("null field"));
   }
   try {
      model->AddField<float>("pt", nullptr);
      FAIL() << "null fields should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("null field fromWhere"));
   }

   // RNTupleReader and RNTupleWriter
   FileRaii fileGuard("test_ntuple_null_safety.root");
   try {
      RNTupleWriter ntuple(nullptr,
         std::make_unique<RPageSinkFile>("myNTuple", fileGuard.GetPath(), RNTupleWriteOptions()));
      FAIL() << "null models should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("null model"));
   }
   try {
      RNTupleWriter ntuple(std::move(model), nullptr);
      FAIL() << "null sinks should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("null sink"));
   }
   try {
      auto ntuple = RNTupleWriter::Recreate(nullptr, "myNtuple", fileGuard.GetPath());
      FAIL() << "null models should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("null model"));
   }
   try {
      auto file = std::make_unique<TFile>(fileGuard.GetPath().c_str(), "RECREATE", "", 101);
      auto ntuple = RNTupleWriter::Append(nullptr, "myNtuple", *file);
      FAIL() << "null models should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("null model"));
   }

   try {
      auto ntuple = RNTupleReader::Open(nullptr, "myNTuple", fileGuard.GetPath());
      FAIL() << "null models should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("null model"));
   }
   try {
      RNTupleReader ntuple(nullptr,
         std::make_unique<RPageSourceFile>("myNTuple", fileGuard.GetPath(), RNTupleReadOptions()));
      FAIL() << "null models should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("null model"));
   }
   try {
      RNTupleReader ntuple(RNTupleModel::Create(), nullptr);
      FAIL() << "null sources should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("null source"));
   }
   try {
      RNTupleReader ntuple(nullptr);
      FAIL() << "null sources should throw";
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("null source"));
   }
}
