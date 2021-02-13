#include "ntuple_test.hxx"

// Tests ReadV() in RColumn.hxx (the case where a std::string overflows to the next page)
TEST(RNTuple, ReadString)
{
   const std::string_view ntupleName = "rs";
   constexpr int numEntries = 2500;
   const std::string contentString = "foooooo";

   FileRaii fileGuard("test_ntuple_readstring.root");
   {
      auto model = RNTupleModel::Create();
      auto st = model->MakeField<std::string>("st");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleName, fileGuard.GetPath());

      for (int i = 0; i < numEntries; ++i) {
         *st = contentString;
         ntuple->Fill();
      }
   }

   auto ntuple = RNTupleReader::Open(ntupleName, fileGuard.GetPath());
   auto viewSt = ntuple->GetView<std::string>("st");
   if (ntuple->GetDescriptor().GetClusterDescriptor(0).GetPageRange(1).fPageInfos.size() < 2) {
      FAIL(); // This means all entries are inside the same page and numEntries should be increased.
   }
   int nElementsPerPage = ntuple->GetDescriptor().GetClusterDescriptor(0).GetPageRange(1).fPageInfos.at(1).fNElements;
   EXPECT_EQ(contentString, viewSt(nElementsPerPage/7));
}
