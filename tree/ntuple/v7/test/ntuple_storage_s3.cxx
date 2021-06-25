#include "ntuple_test.hxx"

TEST(RNTuple, S3Basics)
{
   std::string s3Uri("s3://ntpl0.s3.us-east-2.amazonaws.com");
   {
      auto model = RNTupleModel::Create();
      auto pt = model->MakeField<float>("pt");
      auto vec = model->MakeField<std::vector<int>>("vec");
      auto writer = RNTupleWriter::Recreate(std::move(model), "my_ntuple", s3Uri);
      for (int i = 0; i < 100; i++) {
         *pt = 42.0;
         *vec = {1, 2, 3};
         writer->Fill();
      }
   }

   auto ntuple = RNTupleReader::Open("my_ntuple", s3Uri);
   ntuple->PrintInfo();
   ntuple->PrintInfo(ROOT::Experimental::ENTupleInfo::kStorageDetails);
   EXPECT_EQ(100U, ntuple->GetNEntries());

   auto pt = ntuple->GetView<float>("pt");
   auto vec = ntuple->GetView<std::vector<int>>("vec");
   for (auto i : ntuple->GetEntryRange()) {
      EXPECT_EQ(42.0, pt(i));
      EXPECT_EQ((std::vector<int>{1, 2, 3}), vec(i));
   }
}
