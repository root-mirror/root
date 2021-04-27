#include "ntuple_test.hxx"

TEST(RPageStorageFriends, Null)
{
   FileRaii fileGuard1("test_ntuple_friends_null1.root");
   FileRaii fileGuard2("test_ntuple_friends_null2.root");

   auto model1 = RNTupleModel::Create();
   auto model2 = RNTupleModel::Create();
   {
      auto ntuple1 = RNTupleWriter::Recreate(std::move(model1), "ntpl1", fileGuard1.GetPath());
      auto ntuple2 = RNTupleWriter::Recreate(std::move(model2), "ntpl2", fileGuard2.GetPath());
   }

   std::vector<std::unique_ptr<RPageSource>> realSources;
   realSources.emplace_back(std::make_unique<RPageSourceFile>("ntpl1", fileGuard1.GetPath(), RNTupleReadOptions()));
   realSources.emplace_back(std::make_unique<RPageSourceFile>("ntpl2", fileGuard2.GetPath(), RNTupleReadOptions()));
   RPageSourceFriends friendSource("myNTuple", realSources);
   friendSource.Attach();
   EXPECT_EQ(0u, friendSource.GetNEntries());
}


TEST(RPageStorageFriends, Basic)
{
   FileRaii fileGuard1("test_ntuple_friends_basic1.root");
   FileRaii fileGuard2("test_ntuple_friends_basic2.root");

   auto model1 = RNTupleModel::Create();
   auto fieldPt = model1->MakeField<float>("pt", 42.0);

   auto model2 = RNTupleModel::Create();
   auto fieldEta = model2->MakeField<float>("eta", 24.0);

   {
      auto ntuple = RNTupleWriter::Recreate(std::move(model1), "ntpl1", fileGuard1.GetPath());
      *fieldPt = 1.0;
      ntuple->Fill();
      ntuple->CommitCluster();
      *fieldPt = 2.0;
      ntuple->Fill();
      *fieldPt = 3.0;
      ntuple->Fill();
   }
   {
      auto ntuple = RNTupleWriter::Recreate(std::move(model2), "ntpl2", fileGuard2.GetPath());
      *fieldEta = 4.0;
      ntuple->Fill();
      *fieldEta = 5.0;
      ntuple->Fill();
      ntuple->CommitCluster();
      *fieldEta = 6.0;
      ntuple->Fill();
   }

   std::vector<std::unique_ptr<RPageSource>> realSources;
   realSources.emplace_back(std::make_unique<RPageSourceFile>("ntpl1", fileGuard1.GetPath(), RNTupleReadOptions()));
   realSources.emplace_back(std::make_unique<RPageSourceFile>("ntpl2", fileGuard2.GetPath(), RNTupleReadOptions()));
   RNTupleReader ntuple(std::make_unique<RPageSourceFriends>("myNTuple", realSources));
   EXPECT_EQ(3u, ntuple.GetNEntries());

   auto viewPt = ntuple.GetView<float>("ntpl1.pt");
//   auto viewEta = ntuple.GetView<float>("ntpl2.eta");

//   EXPECT_DOUBLE_EQ(1.0, viewPt(0));
//   EXPECT_DOUBLE_EQ(2.0, viewPt(1));
//   EXPECT_DOUBLE_EQ(3.0, viewPt(2));
//
//   EXPECT_DOUBLE_EQ(4.0, viewEta(0));
//   EXPECT_DOUBLE_EQ(5.0, viewEta(1));
//   EXPECT_DOUBLE_EQ(6.0, viewEta(2));
}


TEST(RPageStorageFriends, FailOnNtupleNameClash)
{
   FileRaii fileGuard1("test_ntuple_friends_name1.root");
   FileRaii fileGuard2("test_ntuple_friends_name2.root");

   auto model1 = RNTupleModel::Create();
   auto model2 = RNTupleModel::Create();
   {
      auto ntuple1 = RNTupleWriter::Recreate(std::move(model1), "ntpl", fileGuard1.GetPath());
      auto ntuple2 = RNTupleWriter::Recreate(std::move(model2), "ntpl", fileGuard2.GetPath());
   }

   std::vector<std::unique_ptr<RPageSource>> realSources;
   realSources.emplace_back(std::make_unique<RPageSourceFile>("ntpl", fileGuard1.GetPath(), RNTupleReadOptions()));
   realSources.emplace_back(std::make_unique<RPageSourceFile>("ntpl", fileGuard2.GetPath(), RNTupleReadOptions()));
   RPageSourceFriends friendSource("myNTuple", realSources);
   EXPECT_THROW(friendSource.Attach(), ROOT::Experimental::RException);
}

TEST(RPageStorageFriends, FailOnEntryCountMismatch)
{
   FileRaii fileGuard1("test_ntuple_friends_count1.root");
   FileRaii fileGuard2("test_ntuple_friends_count2.root");

   auto model1 = RNTupleModel::Create();
   auto fieldPt = model1->MakeField<float>("pt", 42.0);
   auto model2 = RNTupleModel::Create();

   {
      auto ntuple1 = RNTupleWriter::Recreate(std::move(model1), "ntpl1", fileGuard1.GetPath());
      *fieldPt = 1.0;
      ntuple1->Fill();
      auto ntuple2 = RNTupleWriter::Recreate(std::move(model2), "ntpl2", fileGuard2.GetPath());
   }

   std::vector<std::unique_ptr<RPageSource>> realSources;
   realSources.emplace_back(std::make_unique<RPageSourceFile>("ntpl1", fileGuard1.GetPath(), RNTupleReadOptions()));
   realSources.emplace_back(std::make_unique<RPageSourceFile>("ntpl2", fileGuard2.GetPath(), RNTupleReadOptions()));
   RPageSourceFriends friendSource("myNTuple", realSources);
   EXPECT_THROW(friendSource.Attach(), ROOT::Experimental::RException);
}
