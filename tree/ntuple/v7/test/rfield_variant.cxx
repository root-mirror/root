#include "ntuple_test.hxx"

#if __cplusplus >= 201703L
TEST(RNTuple, Variant)
{
   FileRaii fileGuard("test_ntuple_variant.root");

   auto modelWrite = RNTupleModel::Create();
   auto wrVariant = modelWrite->MakeField<std::variant<double, int>>("variant");
   *wrVariant = 2.0;

   auto modelRead = std::unique_ptr<RNTupleModel>(modelWrite->Clone());

   {
      RNTupleWriter ntuple(std::move(modelWrite),
         std::make_unique<RPageSinkFile>("myNTuple", fileGuard.GetPath(), RNTupleWriteOptions()));
      ntuple.Fill();
      ntuple.CommitCluster();
      *wrVariant = 4;
      ntuple.Fill();
      *wrVariant = 8.0;
      ntuple.Fill();
   }
   auto rdVariant = modelRead->Get<std::variant<double, int>>("variant");

   RNTupleReader ntuple(std::move(modelRead),
      std::make_unique<RPageSourceFile>("myNTuple", fileGuard.GetPath(), RNTupleReadOptions()));
   EXPECT_EQ(3U, ntuple.GetNEntries());

   ntuple.LoadEntry(0);
   EXPECT_EQ(2.0, *std::get_if<double>(rdVariant));
   ntuple.LoadEntry(1);
   EXPECT_EQ(4, *std::get_if<int>(rdVariant));
   ntuple.LoadEntry(2);
   EXPECT_EQ(8.0, *std::get_if<double>(rdVariant));
}

#else
TEST(RNTuple, VariantNotSupported)
{
}
#endif // __cplusplus >= 201703L
