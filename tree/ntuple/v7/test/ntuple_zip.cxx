#include "ntuple_test.hxx"

TEST(RNTupleZip, Basics)
{
   RNTupleCompressor compressor;
   RNTupleDecompressor decompressor;

   std::string data = "xxxxxxxxxxxxxxxxxxxxxxxx";
   auto szZipped = compressor(data.data(), data.length(), 101);
   EXPECT_LT(szZipped, data.length());
   auto unzipBuffer = std::unique_ptr<char[]>(new char[data.length()]);
   decompressor(compressor.GetZipBuffer(), szZipped, data.length(), unzipBuffer.get());
   EXPECT_EQ(data, std::string(unzipBuffer.get(), data.length()));

   // inplace decompression
   auto zipBuffer = std::unique_ptr<unsigned char[]>(new unsigned char [data.length()]);
   memcpy(zipBuffer.get(), compressor.GetZipBuffer(), szZipped);
   decompressor(zipBuffer.get(), szZipped, data.length());
   EXPECT_EQ(data, std::string(reinterpret_cast<char *>(zipBuffer.get()), data.length()));
}


TEST(RNTupleZip, Empty)
{
   RNTupleCompressor compressor;

   char x;
   EXPECT_EQ(0U, compressor(&x, 0, 0));
   EXPECT_EQ(0U, compressor(&x, 0, 101));

   // Don't crash
   RNTupleDecompressor()(&x, 0, 0, &x);
}


TEST(RNTupleZip, Uncompressed)
{
   RNTupleCompressor compressor;
   char X = 'x';
   EXPECT_EQ(1U, compressor(&X, 1, 0));
   RNTupleDecompressor()(compressor.GetZipBuffer(), 1, 1, &X);
   EXPECT_EQ('x', X);
}


TEST(RNTupleZip, Small)
{
   RNTupleCompressor compressor;
   char X = 'x';
   EXPECT_EQ(1U, compressor(&X, 1, 101));
   RNTupleDecompressor()(compressor.GetZipBuffer(), 1, 1, &X);
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
   auto szZip = compressor(data.data(), kMAXZIPBUF + 1, 101,
      [&zipBuffer, &nWrites](const void *buffer, size_t nbytes, size_t offset) {
         memcpy(zipBuffer.get() + offset, buffer, nbytes);
         nWrites++;
      });
   EXPECT_EQ(2, nWrites);
   EXPECT_EQ(static_cast<unsigned int>(kMAXZIPBUF) + 1, szZip);

   nWrites = 0;
   szZip = compressor(data.data(), data.length(), 101,
      [&zipBuffer, &nWrites](const void *buffer, size_t nbytes, size_t offset) {
         memcpy(zipBuffer.get() + offset, buffer, nbytes);
         nWrites++;
      });
   EXPECT_LT(szZip, N);
   EXPECT_EQ(2, nWrites);
   decompressor(zipBuffer.get(), szZip, N, unzipBuffer.get());
   EXPECT_EQ(data, std::string(unzipBuffer.get(), N));
}
