#include "ntuple_test.hxx"

TEST(RNTupleZip, Basics)
{
   RNTupleCompressor compressor;
   RNTupleDecompressor decompressor;

   std::string data = "xxxxxxxxxxxxxxxxxxxxxxxx";
   auto szZipped = compressor.Zip(data.data(), data.length(), 101);
   EXPECT_LT(szZipped, data.length());
   auto unzipBuffer = std::unique_ptr<char[]>(new char[data.length()]);
   decompressor.Unzip(compressor.GetZipBuffer(), szZipped, data.length(), unzipBuffer.get());
   EXPECT_EQ(data, std::string(unzipBuffer.get(), data.length()));

   // inplace decompression
   auto zipBuffer = std::unique_ptr<unsigned char[]>(new unsigned char [data.length()]);
   memcpy(zipBuffer.get(), compressor.GetZipBuffer(), szZipped);
   decompressor.Unzip(zipBuffer.get(), szZipped, data.length());
   EXPECT_EQ(data, std::string(reinterpret_cast<char *>(zipBuffer.get()), data.length()));
}


TEST(RNTupleZip, Empty)
{
   RNTupleCompressor compressor;

   char x;
   EXPECT_EQ(0U, compressor.Zip(&x, 0, 0));
   EXPECT_EQ(0U, compressor.Zip(&x, 0, 101));

   // Don't crash
   RNTupleDecompressor().Unzip(&x, 0, 0, &x);
}


TEST(RNTupleZip, Uncompressed)
{
   RNTupleCompressor compressor;
   char X = 'x';
   EXPECT_EQ(1U, compressor.Zip(&X, 1, 0));
   RNTupleDecompressor().Unzip(compressor.GetZipBuffer(), 1, 1, &X);
   EXPECT_EQ('x', X);
}


TEST(RNTupleZip, Small)
{
   RNTupleCompressor compressor;
   char X = 'x';
   EXPECT_EQ(1U, compressor.Zip(&X, 1, 101));
   RNTupleDecompressor().Unzip(compressor.GetZipBuffer(), 1, 1, &X);
   EXPECT_EQ('x', X);
}


TEST(RNTupleZip, Large)
{
   constexpr unsigned int N = kMAXZIPBUF + 32;
   auto zipBuffer = std::make_unique<unsigned char[]>(N);
   auto unzipBuffer = std::make_unique<char[]>(N);
   std::string data(N, 'x');

   RNTupleCompressor compressor;
   RNTupleDecompressor decompressor;

   /// Trailing byte cannot be compressed, entire buffer returns uncompressed
   int nWrites = 0;
   auto szZip = compressor.Zip(data.data(), kMAXZIPBUF + 1, 101,
      [&zipBuffer, &nWrites](const void *buffer, size_t nbytes, size_t offset) {
         memcpy(zipBuffer.get() + offset, buffer, nbytes);
         nWrites++;
      });
   EXPECT_EQ(2, nWrites);
   EXPECT_EQ(static_cast<unsigned int>(kMAXZIPBUF) + 1, szZip);

   nWrites = 0;
   szZip = compressor.Zip(data.data(), data.length(), 101,
      [&zipBuffer, &nWrites](const void *buffer, size_t nbytes, size_t offset) {
         memcpy(zipBuffer.get() + offset, buffer, nbytes);
         nWrites++;
      });
   EXPECT_LT(szZip, N);
   EXPECT_EQ(2, nWrites);
   decompressor.Unzip(zipBuffer.get(), szZip, N, unzipBuffer.get());
   EXPECT_EQ(data, std::string(unzipBuffer.get(), N));
}

TEST(RNTupleZip, CompressionOverride)
{
   RNTupleWriteOptions options;
   EXPECT_EQ(404, options.GetCompression());
   EXPECT_FALSE(options.IsCompressionOverride());
   options.SetCompression(404);
   EXPECT_TRUE(options.IsCompressionOverride());
}

TEST(RNTupleZip, TFilePtrCompressionSettings)
{
   FileRaii fileGuard("test_ntuple_zip_tfileptr_comp.root");
   {
      std::unique_ptr<TFile> file = nullptr;
      auto model = RNTupleModel::Create();
      auto field = model->MakeField<float>("field");
      auto klassVec = model->MakeField<std::vector<CustomStruct>>("klassVec");
      RNTupleWriteOptions options;
      options.SetCompression(407);
      auto ntuple = std::make_unique<RNTupleWriter>(std::move(model),
         std::make_unique<RPageSinkFile>("ntuple", fileGuard.GetPath(), options, file
      ));
      for (int i = 0; i < 2; i++) {
         *field = static_cast<float>(i);
         CustomStruct klass;
         klass.s = std::to_string(i);
         *klassVec = {klass};
         ntuple->Fill();
         ntuple->CommitCluster();
      }
   }
   // ... ntuple uses the explicitly set compression level
   auto ntuple = RNTupleReader::Open("ntuple", fileGuard.GetPath());
#if !defined(_MSC_VER) || defined(R__ENABLE_BROKEN_WIN_TESTS)
   std::ostringstream oss;
   ntuple->PrintInfo(ROOT::Experimental::ENTupleInfo::kStorageDetails, oss);
   EXPECT_THAT(oss.str(), testing::HasSubstr("Compression: 407"));
#endif
   auto rdField = ntuple->GetView<float>("field");
   auto klassVecField = ntuple->GetView<std::vector<CustomStruct>>("klassVec");

   EXPECT_EQ(2, ntuple->GetNEntries());

   for (auto i : ntuple->GetEntryRange()) {
      ASSERT_EQ(static_cast<float>(i), rdField(i));
      ASSERT_EQ(std::to_string(i), klassVecField(i).at(0).s);
   }
}

TEST(RNTupleZip, TFileCompressionSettings)
{
   // TFile compression will be set to 101 (zlib), but some RNTuples will override this setting
   FileRaii fileGuard("test_ntuple_zip_tfile_comp.root");
   RNTupleWriteOptions overrideCompression;
   overrideCompression.SetCompression(505);
   auto file = std::make_unique<TFile>(fileGuard.GetPath().c_str(), "RECREATE", "", 101);
   // test using RPageSinkFile constructor taking TFile&
   {
      auto model = RNTupleModel::Create();
      auto field = model->MakeField<float>("field");
      auto ntuple1 = std::make_unique<RNTupleWriter>(std::move(model),
         std::make_unique<RPageSinkFile>("ntuple1", *file, overrideCompression));
      ntuple1->Fill();
   }
   // test using RNTupleWriter::Append (which calls the RPageSinkFile TFile& constructor)
   {
      auto make_ntuple = [&](std::string name, const RNTupleWriteOptions& opt) {
         auto model = RNTupleModel::Create();
         auto field = model->MakeField<float>("field");
         return RNTupleWriter::Append(std::move(model), name, *file, opt);
      };
      // ntuple2 has the default compression setting and so will use the file's setting
      auto ntuple2 = make_ntuple("ntuple2", RNTupleWriteOptions());
      // ntuple3 has a specific compression setting and will use it
      auto ntuple3 = make_ntuple("ntuple3", overrideCompression);
      // ntuple4 has the default explicitly set to avoid using the file's setting
      RNTupleWriteOptions defaultCompression;
      defaultCompression.SetCompression(404);
      auto ntuple4 = make_ntuple("ntuple4", defaultCompression);
      ntuple2->Fill();
      ntuple3->Fill();
      ntuple4->Fill();
   }
   file.reset();

#if !defined(_MSC_VER) || defined(R__ENABLE_BROKEN_WIN_TESTS)
   std::ostringstream oss;
   // ... ntuple1 uses a specific compression level
   auto ntuple1 = RNTupleReader::Open("ntuple1", fileGuard.GetPath());
   ntuple1->PrintInfo(ROOT::Experimental::ENTupleInfo::kStorageDetails, oss);
   EXPECT_THAT(oss.str(), testing::HasSubstr("Compression: 505"));
   oss.str("");
   // ... ntuple2 uses the TFile's compression level
   auto ntuple2 = RNTupleReader::Open("ntuple2", fileGuard.GetPath());
   ntuple2->PrintInfo(ROOT::Experimental::ENTupleInfo::kStorageDetails, oss);
   EXPECT_THAT(oss.str(), testing::HasSubstr("Compression: 101"));
   oss.str("");
   // ... ntuple3 uses a specific compression level
   auto ntuple3 = RNTupleReader::Open("ntuple3", fileGuard.GetPath());
   ntuple3->PrintInfo(ROOT::Experimental::ENTupleInfo::kStorageDetails, oss);
   EXPECT_THAT(oss.str(), testing::HasSubstr("Compression: 505"));
   oss.str("");
   // ... ntuple4 uses a specific compression level
   auto ntuple4 = RNTupleReader::Open("ntuple4", fileGuard.GetPath());
   ntuple4->PrintInfo(ROOT::Experimental::ENTupleInfo::kStorageDetails, oss);
   EXPECT_THAT(oss.str(), testing::HasSubstr("Compression: 404"));
   oss.str("");
#endif
}

TEST(RNTupleZip, TFileCompressionUpdated)
{
   FileRaii fileGuard("test_ntuple_zip_tfile_comp_updated.root");
   auto file = std::make_unique<TFile>(fileGuard.GetPath().c_str(), "RECREATE", "", 101);
   {
      auto model = RNTupleModel::Create();
      auto field = model->MakeField<float>("field");
      // ntuple is created when TFile has compression setting 101
      auto ntuple = std::make_unique<RNTupleWriter>(std::move(model),
         std::make_unique<RPageSinkFile>("ntuple", *file, RNTupleWriteOptions()));
      // if the TFile compression is later adjusted, this will be picked up by the ntuple
      file->SetCompressionSettings(404);
      ntuple->Fill();
   }
#if !defined(_MSC_VER) || defined(R__ENABLE_BROKEN_WIN_TESTS)
   std::ostringstream oss;
   auto ntuple = RNTupleReader::Open("ntuple", fileGuard.GetPath());
   ntuple->PrintInfo(ROOT::Experimental::ENTupleInfo::kStorageDetails, oss);
   EXPECT_THAT(oss.str(), testing::HasSubstr("Compression: 404"));
#endif
}
