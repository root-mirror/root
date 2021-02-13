#include "ntuple_test.hxx"

TEST(RNTuple, Descriptor)
{
   RNTupleDescriptorBuilder descBuilder;
   descBuilder.SetNTuple("MyTuple", "Description", "Me", RNTupleVersion(1, 2, 3), ROOT::Experimental::RNTupleUuid());
   descBuilder.AddField(RDanglingFieldDescriptor()
      .FieldId(0)
      .FieldName("")
      .Structure(ENTupleStructure::kRecord)
      .MakeDescriptor()
      .Unwrap());
   descBuilder.AddField(RDanglingFieldDescriptor()
      .FieldId(1)
      .FieldName("list")
      .TypeName("std::vector<std::int32_t>")
      .Structure(ENTupleStructure::kCollection)
      .MakeDescriptor()
      .Unwrap());
   descBuilder.AddFieldLink(0, 1);

   descBuilder.AddField(RDanglingFieldDescriptor()
      .FieldId(2)
      .FieldName("list") // at different levels, duplicate names are fine
      .TypeName("std::int32_t")
      .Structure(ENTupleStructure::kLeaf)
      .MakeDescriptor()
      .Unwrap());
   descBuilder.AddFieldLink(1, 2);

   descBuilder.AddField(RDanglingFieldDescriptor()
      .FieldId(42)
      .FieldName("x")
      .TypeName("std::string")
      .Structure(ENTupleStructure::kLeaf)
      .MakeDescriptor()
      .Unwrap());
   descBuilder.AddFieldLink(0, 42);

   descBuilder.AddColumn(3, 42, RNTupleVersion(), RColumnModel(EColumnType::kIndex, true), 0);
   descBuilder.AddColumn(4, 42, RNTupleVersion(), RColumnModel(EColumnType::kByte, true), 1);

   ROOT::Experimental::RClusterDescriptor::RColumnRange columnRange;
   ROOT::Experimental::RClusterDescriptor::RPageRange::RPageInfo pageInfo;
   // Description of cluster #0
   descBuilder.AddCluster(0, RNTupleVersion(), 0, ROOT::Experimental::ClusterSize_t(100));
   columnRange.fColumnId = 3;
   columnRange.fFirstElementIndex = 0;
   columnRange.fNElements = 100;
   descBuilder.AddClusterColumnRange(0, columnRange);
   ROOT::Experimental::RClusterDescriptor::RPageRange pageRange0;
   pageRange0.fPageInfos.clear();
   pageRange0.fColumnId = 3;
   pageInfo.fNElements = 40;
   pageInfo.fLocator.fPosition = 0;
   pageRange0.fPageInfos.emplace_back(pageInfo);
   pageInfo.fNElements = 60;
   pageInfo.fLocator.fPosition = 1024;
   pageRange0.fPageInfos.emplace_back(pageInfo);
   descBuilder.AddClusterPageRange(0, std::move(pageRange0));

   columnRange.fColumnId = 4;
   columnRange.fFirstElementIndex = 0;
   columnRange.fNElements = 300;
   descBuilder.AddClusterColumnRange(0, columnRange);
   ROOT::Experimental::RClusterDescriptor::RPageRange pageRange1;
   pageRange1.fColumnId = 4;
   pageInfo.fNElements = 200;
   pageInfo.fLocator.fPosition = 2048;
   pageRange1.fPageInfos.emplace_back(pageInfo);
   pageInfo.fNElements = 100;
   pageInfo.fLocator.fPosition = 4096;
   pageRange1.fPageInfos.emplace_back(pageInfo);
   descBuilder.AddClusterPageRange(0, std::move(pageRange1));

   // Description of cluster #1
   descBuilder.AddCluster(1, RNTupleVersion(), 100, ROOT::Experimental::ClusterSize_t(1000));
   columnRange.fColumnId = 3;
   columnRange.fFirstElementIndex = 100;
   columnRange.fNElements = 1000;
   descBuilder.AddClusterColumnRange(1, columnRange);
   ROOT::Experimental::RClusterDescriptor::RPageRange pageRange2;
   pageRange2.fColumnId = 3;
   pageInfo.fNElements = 1000;
   pageInfo.fLocator.fPosition = 8192;
   pageRange2.fPageInfos.emplace_back(pageInfo);
   descBuilder.AddClusterPageRange(1, std::move(pageRange2));

   columnRange.fColumnId = 4;
   columnRange.fFirstElementIndex = 300;
   columnRange.fNElements = 3000;
   descBuilder.AddClusterColumnRange(1, columnRange);
   ROOT::Experimental::RClusterDescriptor::RPageRange pageRange3;
   pageRange3.fColumnId = 4;
   pageInfo.fNElements = 3000;
   pageInfo.fLocator.fPosition = 16384;
   pageRange3.fPageInfos.emplace_back(pageInfo);
   descBuilder.AddClusterPageRange(1, std::move(pageRange3));

   const auto &reference = descBuilder.GetDescriptor();
   EXPECT_EQ("MyTuple", reference.GetName());
   EXPECT_EQ(1U, reference.GetVersion().GetVersionUse());
   EXPECT_EQ(2U, reference.GetVersion().GetVersionMin());
   EXPECT_EQ(3U, reference.GetVersion().GetFlags());

   EXPECT_EQ(1U, reference.FindNextClusterId(0));
   EXPECT_EQ(ROOT::Experimental::kInvalidDescriptorId, reference.FindNextClusterId(1));
   EXPECT_EQ(0U, reference.FindPrevClusterId(1));
   EXPECT_EQ(ROOT::Experimental::kInvalidDescriptorId, reference.FindPrevClusterId(0));

   auto szHeader = reference.GetHeaderSize();
   auto headerBuffer = new unsigned char[szHeader];
   reference.SerializeHeader(headerBuffer);
   auto szFooter = reference.GetFooterSize();
   auto footerBuffer = new unsigned char[szFooter];
   reference.SerializeFooter(footerBuffer);

   const auto nbytesPostscript = RNTupleDescriptor::kNBytesPostscript;
   ASSERT_GE(szFooter, nbytesPostscript);
   std::uint32_t szPsHeader;
   std::uint32_t szPsFooter;
   RNTupleDescriptor::LocateMetadata(footerBuffer + szFooter - nbytesPostscript, szPsHeader, szPsFooter);
   EXPECT_EQ(szHeader, szPsHeader);
   EXPECT_EQ(szFooter, szPsFooter);

   RNTupleDescriptorBuilder reco;
   reco.SetFromHeader(headerBuffer);
   reco.AddClustersFromFooter(footerBuffer);
   EXPECT_EQ(reference, reco.GetDescriptor());

   EXPECT_EQ(NTupleSize_t(1100), reference.GetNEntries());
   EXPECT_EQ(NTupleSize_t(1100), reference.GetNElements(3));
   EXPECT_EQ(NTupleSize_t(3300), reference.GetNElements(4));

   EXPECT_EQ(DescriptorId_t(0), reference.GetFieldZeroId());
   EXPECT_EQ(DescriptorId_t(1), reference.FindFieldId("list", 0));
   EXPECT_EQ(DescriptorId_t(1), reference.FindFieldId("list"));
   EXPECT_EQ(DescriptorId_t(2), reference.FindFieldId("list", 1));
   EXPECT_EQ(DescriptorId_t(42), reference.FindFieldId("x", 0));
   EXPECT_EQ(DescriptorId_t(42), reference.FindFieldId("x"));
   EXPECT_EQ(ROOT::Experimental::kInvalidDescriptorId, reference.FindFieldId("listX", 1));
   EXPECT_EQ(ROOT::Experimental::kInvalidDescriptorId, reference.FindFieldId("list", 1024));

   EXPECT_EQ(DescriptorId_t(3), reference.FindColumnId(42, 0));
   EXPECT_EQ(DescriptorId_t(4), reference.FindColumnId(42, 1));
   EXPECT_EQ(ROOT::Experimental::kInvalidDescriptorId, reference.FindColumnId(42, 2));
   EXPECT_EQ(ROOT::Experimental::kInvalidDescriptorId, reference.FindColumnId(43, 0));

   EXPECT_EQ(DescriptorId_t(0), reference.FindClusterId(3, 0));
   EXPECT_EQ(DescriptorId_t(1), reference.FindClusterId(3, 100));
   EXPECT_EQ(ROOT::Experimental::kInvalidDescriptorId, reference.FindClusterId(3, 40000));

   delete[] footerBuffer;
   delete[] headerBuffer;
}

TEST(RDanglingFieldDescriptor, MakeDescriptorErrors)
{
   // minimum requirements for making a field descriptor from scratch
   RFieldDescriptor fieldDesc = RDanglingFieldDescriptor()
      .FieldId(1)
      .Structure(ENTupleStructure::kCollection)
      .FieldName("someField")
      .MakeDescriptor()
      .Unwrap();

   // MakeDescriptor() returns an RResult<RFieldDescriptor>
   // -- here we check the error cases

   // must set field id
   RResult<RFieldDescriptor> fieldDescRes = RDanglingFieldDescriptor().MakeDescriptor();
   ASSERT_FALSE(fieldDescRes) << "default constructed dangling descriptors should throw";
   EXPECT_THAT(fieldDescRes.GetError()->GetReport(), testing::HasSubstr("invalid field id"));

   // must set field structure
   fieldDescRes = RDanglingFieldDescriptor()
      .FieldId(1)
      .MakeDescriptor();
   ASSERT_FALSE(fieldDescRes) << "field descriptors without structure should throw";
   EXPECT_THAT(fieldDescRes.GetError()->GetReport(), testing::HasSubstr("invalid field structure"));

   // must set field name
   fieldDescRes = RDanglingFieldDescriptor()
      .FieldId(1)
      .Structure(ENTupleStructure::kCollection)
      .MakeDescriptor();
   ASSERT_FALSE(fieldDescRes) << "unnamed field descriptors should throw";
   EXPECT_THAT(fieldDescRes.GetError()->GetReport(), testing::HasSubstr("name cannot be empty string"));
}

TEST(RNTupleDescriptorBuilder, CatchBadLinks)
{
   RNTupleDescriptorBuilder descBuilder;
   descBuilder.AddField(RDanglingFieldDescriptor()
      .FieldId(0)
      .Structure(ENTupleStructure::kRecord)
      .MakeDescriptor()
      .Unwrap());
   descBuilder.AddField(RDanglingFieldDescriptor()
      .FieldId(1)
      .FieldName("field")
      .TypeName("int32_t")
      .Structure(ENTupleStructure::kLeaf)
      .MakeDescriptor()
      .Unwrap());
   try {
      descBuilder.AddFieldLink(1, 0);
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("cannot make FieldZero a child field"));
   }
   try {
      descBuilder.AddFieldLink(1, 1);
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("cannot make field '1' a child of itself"));
   }
   try {
      descBuilder.AddFieldLink(1, 10);
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("child field with id '10' doesn't exist in NTuple"));
   }
}

TEST(RNTupleDescriptorBuilder, CatchInvalidDescriptors)
{
   RNTupleDescriptorBuilder descBuilder;

   // empty string is not a valid NTuple name
   descBuilder.SetNTuple("", "", "", RNTupleVersion(1, 2, 3), ROOT::Experimental::RNTupleUuid());
   try {
      descBuilder.EnsureValidDescriptor();
   } catch (const RException& err) {
      EXPECT_THAT(err.what(), testing::HasSubstr("name cannot be empty string"));
   }
   descBuilder.SetNTuple("something", "", "", RNTupleVersion(1, 2, 3), ROOT::Experimental::RNTupleUuid());
   descBuilder.EnsureValidDescriptor();
}

TEST(RNTupleDescriptor, QualifiedFieldName)
{
   auto model = RNTupleModel::Create();
   model->MakeField<std::int32_t>("ints");
   model->MakeField<std::vector<float>>("jets");
   FileRaii fileGuard("test_field_qualified.root");
   {
      RNTupleWriter ntuple(std::move(model),
         std::make_unique<RPageSinkFile>("ntuple", fileGuard.GetPath(), RNTupleWriteOptions()));
   }

   auto ntuple = RNTupleReader::Open("ntuple", fileGuard.GetPath());
   const auto &desc = ntuple->GetDescriptor();
   EXPECT_TRUE(desc.GetQualifiedFieldName(desc.GetFieldZeroId()).empty());
   auto fldIdInts = desc.FindFieldId("ints");
   EXPECT_STREQ("ints", desc.GetQualifiedFieldName(fldIdInts).c_str());
   auto fldIdJets = desc.FindFieldId("jets");
   auto fldIdInner = desc.FindFieldId("float", fldIdJets);
   EXPECT_STREQ("jets", desc.GetQualifiedFieldName(fldIdJets).c_str());
   EXPECT_STREQ("jets.float", desc.GetQualifiedFieldName(fldIdInner).c_str());
}

TEST(RFieldDescriptorRange, IterateOverFieldNames)
{
   auto model = RNTupleModel::Create();
   auto floats = model->MakeField<std::vector<float>>("jets");
   auto bools = model->MakeField<std::vector<bool>>("bools");
   auto bool_vec_vec = model->MakeField<std::vector<std::vector<bool>>>("bool_vec_vec");
   auto ints = model->MakeField<std::int32_t>("ints");

   FileRaii fileGuard("test_field_iterator.root");
   {
      RNTupleWriter ntuple(std::move(model),
         std::make_unique<RPageSinkFile>("ntuple", fileGuard.GetPath(), RNTupleWriteOptions()));
      ntuple.Fill();
   }

   auto ntuple = RNTupleReader::Open("ntuple", fileGuard.GetPath());
   // iterate over top-level fields
   std::vector<std::string> names{};
   for (auto& f: ntuple->GetDescriptor().GetTopLevelFields()) {
      names.push_back(f.GetFieldName());
   }
   ASSERT_EQ(names.size(), 4);
   EXPECT_EQ(names[0], std::string("jets"));
   EXPECT_EQ(names[1], std::string("bools"));
   EXPECT_EQ(names[2], std::string("bool_vec_vec"));
   EXPECT_EQ(names[3], std::string("ints"));

   const auto& ntuple_desc = ntuple->GetDescriptor();
   auto top_level_fields = ntuple_desc.GetTopLevelFields();

   // iterate over child field ranges
   const auto& float_vec_desc = *top_level_fields.begin();
   EXPECT_EQ(float_vec_desc.GetFieldName(), std::string("jets"));
   auto float_vec_child_range = ntuple_desc.GetFieldRange(float_vec_desc);
   std::vector<std::string> child_names{};
   for (auto& child_field: float_vec_child_range) {
      child_names.push_back(child_field.GetFieldName());
      // check the empty range
      auto float_child_range = ntuple_desc.GetFieldRange(child_field);
      EXPECT_EQ(float_child_range.begin(), float_child_range.end());
   }
   ASSERT_EQ(child_names.size(), 1);
   EXPECT_EQ(child_names[0], std::string("float"));

   // check if canonical iterator methods work
   auto iter = top_level_fields.begin();
   std::advance(iter, 2);
   const auto& bool_vec_vec_desc = *iter;
   EXPECT_EQ(bool_vec_vec_desc.GetFieldName(), std::string("bool_vec_vec"));

   child_names.clear();
   for (auto& child_field: ntuple_desc.GetFieldRange(bool_vec_vec_desc)) {
      child_names.push_back(child_field.GetFieldName());
   }
   ASSERT_EQ(child_names.size(), 1);
   EXPECT_EQ(child_names[0], std::string("std::vector<bool>"));
}

TEST(RFieldDescriptorRange, SortByLambda)
{
   auto model = RNTupleModel::Create();
   auto floats = model->MakeField<std::vector<float>>("jets");
   auto bools = model->MakeField<std::vector<bool>>("bools");
   auto bool_vec_vec = model->MakeField<std::vector<std::vector<bool>>>("bool_vec_vec");
   auto ints = model->MakeField<std::int32_t>("ints");

   FileRaii fileGuard("test_field_iterator.root");
   {
      RNTupleWriter ntuple(std::move(model),
         std::make_unique<RPageSinkFile>("ntuple", fileGuard.GetPath(), RNTupleWriteOptions()));
      ntuple.Fill();
   }

   auto ntuple = RNTupleReader::Open("ntuple", fileGuard.GetPath());
   const auto& ntuple_desc = ntuple->GetDescriptor();
   auto alpha_order = [&](auto lhs, auto rhs) {
      return ntuple_desc.GetFieldDescriptor(lhs).GetFieldName()
         < ntuple_desc.GetFieldDescriptor(rhs).GetFieldName();
   };

   std::vector<std::string> sorted_names = {};
   for (auto& f: ntuple_desc.GetTopLevelFields(alpha_order)) {
      sorted_names.push_back(f.GetFieldName());
   }
   ASSERT_EQ(sorted_names.size(), 4);
   EXPECT_EQ(sorted_names[0], std::string("bool_vec_vec"));
   EXPECT_EQ(sorted_names[1], std::string("bools"));
   EXPECT_EQ(sorted_names[2], std::string("ints"));
   EXPECT_EQ(sorted_names[3], std::string("jets"));

   // reverse alphabetical
   sorted_names.clear();
   for (auto& f: ntuple_desc.GetTopLevelFields(
      [&](auto lhs, auto rhs) { return !alpha_order(lhs, rhs); }))
   {
      sorted_names.push_back(f.GetFieldName());
   }
   ASSERT_EQ(sorted_names.size(), 4);
   EXPECT_EQ(sorted_names[0], std::string("jets"));
   EXPECT_EQ(sorted_names[1], std::string("ints"));
   EXPECT_EQ(sorted_names[2], std::string("bools"));
   EXPECT_EQ(sorted_names[3], std::string("bool_vec_vec"));

   // alphabetical by type name
   std::vector<std::string> sorted_by_typename = {};
   for (auto& f: ntuple_desc.GetTopLevelFields(
      [&](auto lhs, auto rhs) {
         return ntuple_desc.GetFieldDescriptor(lhs).GetTypeName()
            < ntuple_desc.GetFieldDescriptor(rhs).GetTypeName(); }))
   {
      sorted_by_typename.push_back(f.GetFieldName());
   }
   // int32_t, vector<bool>, vector<float>, vector<vector<bool>
   ASSERT_EQ(sorted_by_typename.size(), 4);
   EXPECT_EQ(sorted_by_typename[0], std::string("ints"));
   EXPECT_EQ(sorted_by_typename[1], std::string("bools"));
   EXPECT_EQ(sorted_by_typename[2], std::string("jets"));
   EXPECT_EQ(sorted_by_typename[3], std::string("bool_vec_vec"));
}


TEST(RColumnDescriptorRange, IterateOverColumns)
{
   auto model = RNTupleModel::Create();
   auto floats = model->MakeField<std::vector<float>>("jets");
   auto bools = model->MakeField<std::string>("tag");

   FileRaii fileGuard("test_column_iterator.root");
   {
      RNTupleWriter ntuple(std::move(model),
         std::make_unique<RPageSinkFile>("ntuple", fileGuard.GetPath(), RNTupleWriteOptions()));
      ntuple.Fill();
   }

   auto ntuple = RNTupleReader::Open("ntuple", fileGuard.GetPath());
   const auto &desc = ntuple->GetDescriptor();

   // No column attached to the zero field
   unsigned int counter = 0;
   for (const auto &c : desc.GetColumnRange(desc.GetFieldZeroId())) {
      (void)c;
      counter++;
   }
   EXPECT_EQ(0u, counter);

   const auto tagId = desc.FindFieldId("tag");
   for (const auto &c : desc.GetColumnRange(tagId)) {
      EXPECT_EQ(tagId, c.GetFieldId());
      EXPECT_EQ(counter, c.GetIndex());
      counter++;
   }
   EXPECT_EQ(2, counter);

   const auto jetsId = desc.FindFieldId("jets");
   for (const auto &c : desc.GetColumnRange(desc.FindFieldId("jets"))) {
      EXPECT_EQ(jetsId, c.GetFieldId());
      counter++;
   }
   EXPECT_EQ(3, counter);
}
