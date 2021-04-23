#include "ntuple_test.hxx"

namespace {

// Reads an integer from a little-endian 4 byte buffer
std::int32_t ReadRawInt(void *ptr)
{
   std::int32_t val = *reinterpret_cast<std::int32_t *>(ptr);
#ifndef R__BYTESWAP
   // on big endian system
   auto x = (val & 0x0000FFFF) << 16 | (val & 0xFFFF0000) >> 16;
   return (x & 0x00FF00FF) << 8 | (x & 0xFF00FF00) >> 8;
#else
   return val;
#endif
}

} // anonymous namespace

TEST(RPageStorage, ReadSealedPages)
{
   FileRaii fileGuard("test_ntuple_sealed_pages.root");

   // Create an ntuple at least 2 clusters, one with 1 entry and one with 100000 entries.
   // Hence the second cluster should have more than a single page per column.  We write uncompressed
   // pages so that we can meaningfully peek into the content of read sealed pages later on.
   auto model = RNTupleModel::Create();
   auto wrPt = model->MakeField<std::int32_t>("pt", 42);
   {
      RNTupleWriteOptions options;
      options.SetCompression(0);
      RNTupleWriter ntuple(std::move(model),
         std::make_unique<RPageSinkFile>("myNTuple", fileGuard.GetPath(), RNTupleWriteOptions()));
      ntuple.Fill();
      ntuple.CommitCluster();
      for (unsigned i = 0; i < 100000; ++i) {
         *wrPt = i;
         ntuple.Fill();
      }
   }

   RPageSourceFile source("myNTuple", fileGuard.GetPath(), RNTupleReadOptions());
   source.Attach();
   const auto &descriptor = source.GetDescriptor();
   auto columnId = descriptor.FindColumnId(descriptor.FindFieldId("pt"), 0);

   // Check first cluster consisting of a single entry
   RClusterIndex index(descriptor.FindClusterId(columnId, 0), 0);
   auto sealedPage = source.ReadSealedPage(columnId, index);
   ASSERT_EQ(1U, sealedPage.fNElements);
   ASSERT_EQ(4U, sealedPage.fSize);
   EXPECT_EQ(42, ReadRawInt(sealedPage.fBuffer.get()));

   // Check second, big cluster
   auto clusterId = descriptor.FindClusterId(columnId, 1);
   ASSERT_NE(clusterId, index.GetClusterId());
   const auto &clusterDesc = descriptor.GetClusterDescriptor(clusterId);
   const auto &pageRange = clusterDesc.GetPageRange(columnId);
   EXPECT_GT(pageRange.fPageInfos.size(), 1U);
   std::uint32_t firstElementInPage = 0;
   for (const auto &pi : pageRange.fPageInfos) {
      sealedPage = source.ReadSealedPage(columnId, RClusterIndex(clusterId, firstElementInPage));
      ASSERT_GE(sealedPage.fSize, 4U);
      EXPECT_EQ(firstElementInPage, ReadRawInt(sealedPage.fBuffer.get()));
      firstElementInPage += pi.fNElements;
   }
}


TEST(RFieldMerger, FieldMergeEmpty)
{
   auto model = RNTupleModel::Create();
   RPageSinkNull sinkNull;
   sinkNull.Create(*model);

   const auto &desc = sinkNull.GetDescriptor();
   RFieldMerger merger({desc, desc.GetFieldZeroId()});
   EXPECT_EQ(desc.GetFieldZeroId(), merger.GetReferenceFieldId(desc.GetFieldZeroId(), 0));
   EXPECT_EQ(desc.GetFieldZeroId(), merger.GetInputFieldId(desc.GetFieldZeroId(), 0));

   auto result = merger.Merge({desc, desc.GetFieldZeroId()});
   EXPECT_EQ(1, result.Inspect());
   EXPECT_EQ(desc.GetFieldZeroId(), merger.GetReferenceFieldId(desc.GetFieldZeroId(), 1));
   EXPECT_EQ(desc.GetFieldZeroId(), merger.GetInputFieldId(desc.GetFieldZeroId(), 1));
}


TEST(RFieldMerger, FieldMergeSame)
{
   auto model = RNTupleModel::Create();
   model->MakeField<float>("pt");
   model->MakeField<std::string>("tag");
   model->MakeField<std::vector<std::string>>("vec");
   RPageSinkNull sinkNull;
   sinkNull.Create(*model);
   // Needs to be destructed before the sink
   model = nullptr;

   const auto &desc = sinkNull.GetDescriptor();
   RFieldMerger merger({desc, desc.GetFieldZeroId()});
   auto result = merger.Merge({desc, desc.GetFieldZeroId()});
   EXPECT_EQ(1, result.Inspect());

   // Check that both string column ids are present
   auto tagId = desc.FindFieldId("tag");
   auto tagColIdx0 = desc.FindColumnId(tagId, 0);
   auto tagColIdx1 = desc.FindColumnId(tagId, 1);
   EXPECT_EQ(tagColIdx0, merger.GetReferenceColumnId(tagColIdx0, 1));
   EXPECT_EQ(tagColIdx1, merger.GetReferenceColumnId(tagColIdx1, 1));

   // Check that sub fields are present
   auto vecId = desc.FindFieldId("vec");
   DescriptorId_t vecSubId = ROOT::Experimental::kInvalidDescriptorId;
   for (const auto &f : desc.GetFieldRange(vecId))
      vecSubId = f.GetId();
   ASSERT_NE(vecSubId, ROOT::Experimental::kInvalidDescriptorId);
   ASSERT_EQ(vecId, desc.GetFieldDescriptor(vecSubId).GetParentId());
   EXPECT_EQ(vecSubId, merger.GetInputFieldId(vecSubId, 0));
   EXPECT_EQ(vecSubId, merger.GetInputFieldId(vecSubId, 1));
}


TEST(RFieldMerger, FieldMergeReorder)
{
   auto model = RNTupleModel::Create();
   model->MakeField<float>("pt");
   model->MakeField<std::string>("tag");
   model->MakeField<std::vector<std::string>>("vec");
   RPageSinkNull sinkRef;
   sinkRef.Create(*model);

   // Same model but different order of fields
   model = RNTupleModel::Create();
   model->MakeField<std::vector<std::string>>("vec");
   model->MakeField<std::string>("tag");
   model->MakeField<float>("pt");
   RPageSinkNull sinkInput1;
   sinkInput1.Create(*model);
   model = nullptr;

   const auto &descRef = sinkRef.GetDescriptor();
   RFieldMerger merger({descRef, descRef.GetFieldZeroId()});

   const auto &descInput1 = sinkInput1.GetDescriptor();
   auto result = merger.Merge({descInput1, descInput1.GetFieldZeroId()});
   ASSERT_EQ(1, result.Inspect());

   auto ptRefId = descRef.FindFieldId("pt");
   auto ptInput1Id = descInput1.FindFieldId("pt");
   EXPECT_NE(ptRefId, ptInput1Id);
   EXPECT_EQ(ptRefId, merger.GetReferenceFieldId(ptInput1Id, 1));
   EXPECT_EQ(ptInput1Id, merger.GetInputFieldId(ptRefId, 1));
}


TEST(RFieldMerger, FieldMergeSlice)
{
   auto model = RNTupleModel::Create();
   model->MakeField<float>("pt");
   RPageSinkNull sinkRef;
   sinkRef.Create(*model);

   // 2nd model has one more field which gets ignored during the merge
   model = RNTupleModel::Create();
   model->MakeField<std::string>("tag");
   model->MakeField<float>("pt");
   RPageSinkNull sinkInput1;
   sinkInput1.Create(*model);
   model = nullptr;

   const auto &descRef = sinkRef.GetDescriptor();
   RFieldMerger merger({descRef, descRef.GetFieldZeroId()});

   const auto &descInput1 = sinkInput1.GetDescriptor();
   auto result = merger.Merge({descInput1, descInput1.GetFieldZeroId()});
   ASSERT_EQ(1, result.Inspect());

   auto ptRefId = descRef.FindFieldId("pt");
   auto ptInput1Id = descInput1.FindFieldId("pt");
   EXPECT_EQ(ptRefId, merger.GetReferenceFieldId(ptInput1Id, 1));
   EXPECT_EQ(ptInput1Id, merger.GetInputFieldId(ptRefId, 1));
}


TEST(RFieldMerger, FieldMergeErrors)
{
   auto model = RNTupleModel::Create();
   model->MakeField<float>("pt");
   RPageSinkNull sinkRef;
   sinkRef.Create(*model);

   const auto &descRef = sinkRef.GetDescriptor();
   RFieldMerger merger({descRef, descRef.GetFieldZeroId()});

   // Failure 1: missing field
   model = RNTupleModel::Create();
   RPageSinkNull sinkInput1;
   sinkInput1.Create(*model);
   const auto &descInput1 = sinkInput1.GetDescriptor();
   auto result1 = merger.Merge({descInput1, descInput1.GetFieldZeroId()});
   EXPECT_FALSE(result1);

   // Failure 2: type mismatch
   model = RNTupleModel::Create();
   model->MakeField<int>("pt");
   RPageSinkNull sinkInput2;
   sinkInput2.Create(*model);
   model = nullptr;
   const auto &descInput2 = sinkInput2.GetDescriptor();
   auto result2 = merger.Merge({descInput2, descInput2.GetFieldZeroId()});
   EXPECT_FALSE(result2);

   // Eventually: successfully merge after errors
   model = RNTupleModel::Create();
   model->MakeField<float>("pt");
   RPageSinkNull sinkInput3;
   sinkInput3.Create(*model);
   model = nullptr;
   const auto &descInput3 = sinkInput3.GetDescriptor();
   auto result3 = merger.Merge({descInput3, descInput3.GetFieldZeroId()});
   EXPECT_EQ(1, result3.Inspect());
}
