#include <gtest/gtest.h>
#include <TMVA/RTensor.hxx>

using namespace TMVA::Experimental;

TEST(RTensor, GetElement)
{
   float data[6] = {0, 1, 2, 3, 4, 5};
   RTensor<float> x(data, {2, 3});
   auto shape = x.GetShape();
   float count = 0.0;
   for (size_t i = 0; i < shape[0]; i++) {
      for (size_t j = 0; j < shape[1]; j++) {
         EXPECT_EQ(count, x({i, j}));
         EXPECT_EQ(count, x(i, j));
         count++;
      }
   }
}

TEST(RTensor, SetElement)
{
   RTensor<float> x({2, 3});
   auto shape = x.GetShape();
   float count = 0.0;
   for (size_t i = 0; i < shape[0]; i++) {
      for (size_t j = 0; j < shape[1]; j++) {
         x({i, j}) = count;
         x(i, j) = count;
         count++;
      }
   }
   auto data = x.GetData();
   for (size_t i = 0; i < shape.size(); i++)
      EXPECT_EQ((float)i, *(data + i));
}

TEST(RTensor, AdoptMemory)
{
   float data[4] = {0, 0, 0, 0};
   RTensor<float> x(data, {4});
   for (size_t i = 0; i < 4; i++)
      x(i) = (float)i;
   for (size_t i = 0; i < 4; i++)
      EXPECT_EQ((float)i, data[i]);
}

TEST(RTensor, GetShape)
{
   RTensor<int> x({2, 3});
   const auto s = x.GetShape();
   EXPECT_EQ(s.size(), 2u);
   EXPECT_EQ(s[0], 2u);
   EXPECT_EQ(s[1], 3u);
}

TEST(RTensor, Reshape)
{
   RTensor<int> x({2, 3});
   const auto s = x.GetShape();
   EXPECT_EQ(s.size(), 2u);
   EXPECT_EQ(s[0], 2u);
   EXPECT_EQ(s[1], 3u);

   auto x2 = x.Reshape({1, 6});
   const auto s2 = x2.GetShape();
   EXPECT_EQ(s2.size(), 2u);
   EXPECT_EQ(s2[0], 1u);
   EXPECT_EQ(s2[1], 6u);

   auto x3 = x.Reshape({6});
   const auto s3 = x3.GetShape();
   EXPECT_EQ(s3.size(), 1u);
   EXPECT_EQ(s3[0], 6u);
}

TEST(RTensor, ExpandDims)
{
   RTensor<int> x({2, 3});
   auto xa = x.ExpandDims(0);
   const auto s = xa.GetShape();
   EXPECT_EQ(s.size(), 3u);
   EXPECT_EQ(s[0], 1u);
   EXPECT_EQ(s[1], 2u);
   EXPECT_EQ(s[2], 3u);

   RTensor<int> x1({2, 3});
   auto xb = x1.ExpandDims(1);
   const auto s1 = xb.GetShape();
   EXPECT_EQ(s1.size(), 3u);
   EXPECT_EQ(s1[0], 2u);
   EXPECT_EQ(s1[1], 1u);
   EXPECT_EQ(s1[2], 3u);

   RTensor<int> x2({2, 3});
   auto xc = x2.ExpandDims(-1);
   const auto s2 = xc.GetShape();
   EXPECT_EQ(s2.size(), 3u);
   EXPECT_EQ(s2[0], 2u);
   EXPECT_EQ(s2[1], 3u);
   EXPECT_EQ(s2[2], 1u);
}

TEST(RTensor, Squeeze)
{
   RTensor<int> x({1, 2, 3});
   auto xa = x.Squeeze();
   const auto s = xa.GetShape();
   EXPECT_EQ(s.size(), 2u);
   EXPECT_EQ(s[0], 2u);
   EXPECT_EQ(s[1], 3u);

   RTensor<int> x1({1, 2, 1, 3, 1});
   auto xb = x1.Squeeze();
   const auto s1 = xb.GetShape();
   EXPECT_EQ(s1.size(), 2u);
   EXPECT_EQ(s1[0], 2u);
   EXPECT_EQ(s1[1], 3u);

   RTensor<int> x2({1, 1, 1});
   auto xc = x2.Squeeze();
   const auto s2 = xc.GetShape();
   EXPECT_EQ(s2.size(), 1u);
   EXPECT_EQ(s2[0], 1u);

   RTensor<int> x3({});
   auto xd = x3.Squeeze();
   const auto s3 = xd.GetShape();
   EXPECT_EQ(s3.size(), 0u);
}

TEST(RTensor, Transpose)
{
   // Layout (physical):
   // 0, 1, 2, 3, 4, 5
   float data[6] = {0, 1, 2, 3, 4, 5};
   RTensor<float> x(data, {2, 3});

   // Layout (logical):
   // 0, 1, 2
   // 3, 4, 5
   float count = 0;
   for (size_t i = 0; i < 2; i++) {
      for (size_t j = 0; j < 3; j++) {
         x(i, j) = count;
         count++;
      }
   }

   // Layout (logical):
   // 0, 3
   // 1, 4
   // 2, 5
   auto x2 = x.Transpose();
   EXPECT_EQ(x2.GetShape().size(), 2u);
   EXPECT_EQ(x2.GetShape()[0], 3u);
   EXPECT_EQ(x2.GetShape()[1], 2u);
   count = 0;
   for (size_t i = 0; i < 2; i++) {
      for (size_t j = 0; j < 3; j++) {
         x2(j, i) = count;
         count++;
      }
   }
}

TEST(RTensor, InitWithZeros)
{
   RTensor<float> x({2, 3});
   for (size_t i = 0; i < 2; i++) {
      for (size_t j = 0; j < 3; j++) {
         EXPECT_EQ(0.0, x(i, j));
      }
   }
}

TEST(RTensor, SetElementRowMajor)
{
   // Layout (logical):
   // 0, 1, 2
   // 3, 4, 5
   RTensor<float> x({2, 3}, MemoryLayout::RowMajor);
   float count = 0;
   for (size_t i = 0; i < 2; i++) {
      for (size_t j = 0; j < 3; j++) {
         x(i, j) = count;
         count++;
      }
   }

   // Layout (physical):
   // 0, 1, 2, 3, 4, 5
   float ref[6] = {0, 1, 2, 3, 4, 5};
   float *data = x.GetData();
   for (size_t i = 0; i < 6; i++) {
      EXPECT_EQ(data[i], ref[i]);
   }
}

TEST(RTensor, SetElementColumnMajor)
{
   // Layout (logical):
   // 0, 1, 2
   // 3, 4, 5
   RTensor<float> x({2, 3}, MemoryLayout::ColumnMajor);
   float count = 0;
   for (size_t i = 0; i < 2; i++) {
      for (size_t j = 0; j < 3; j++) {
         x(i, j) = count;
         count++;
      }
   }

   // Layout (physical):
   // 0, 3, 1, 4, 2, 5
   float ref[6] = {0, 3, 1, 4, 2, 5};
   float *data = x.GetData();
   for (size_t i = 0; i < 6; i++) {
      EXPECT_EQ(data[i], ref[i]);
   }
}

TEST(RTensor, GetElementRowMajor)
{
   // Layout (physical):
   // 0, 1, 2, 3, 4, 5
   float data[6] = {0, 1, 2, 3, 4, 5};
   RTensor<float> x(data, {2, 3}, MemoryLayout::RowMajor);

   // Layout (logical):
   // 0, 1, 2
   // 3, 4, 5
   float count = 0;
   for (size_t i = 0; i < 2; i++) {
      for (size_t j = 0; j < 3; j++) {
         EXPECT_EQ(x(i, j), count);
         count++;
      }
   }
}

TEST(RTensor, GetElementColumnMajor)
{
   // Layout (physical):
   // 0, 3, 1, 4, 2, 5
   float data[6] = {0, 3, 1, 4, 2, 5};
   RTensor<float> x(data, {2, 3}, MemoryLayout::ColumnMajor);

   // Layout (logical):
   // 0, 1, 2
   // 3, 4, 5
   float count = 0;
   for (size_t i = 0; i < 2; i++) {
      for (size_t j = 0; j < 3; j++) {
         EXPECT_EQ(x(i, j), count);
         count++;
      }
   }
}

TEST(RTensor, SliceRowMajor)
{
   // Data layout:
   // [ 0, 1, 2, 3, 4, 5 ] ]
   RTensor<float> x({2, 3}, MemoryLayout::RowMajor);
   float c = 0.f;
   for (auto i = 0; i < 2; i++) {
      for (auto j = 0; j < 3; j++) {
         x(i, j) = c;
         c++;
      }
   }

   // Slice:
   // [ [ 0, 1 ], [ 3, 4 ] ]
   auto s1 = x.Slice({{0, 2}, {0, 2}});
   EXPECT_EQ(s1.GetSize(), 4u);
   EXPECT_EQ(s1.GetShape().size(), 2u);
   EXPECT_EQ(s1.GetShape()[0], 2u);
   EXPECT_EQ(s1.GetShape()[1], 2u);
   EXPECT_EQ(s1(0, 0), 0.f);
   EXPECT_EQ(s1(0, 1), 1.f);
   EXPECT_EQ(s1(1, 0), 3.f);
   EXPECT_EQ(s1(1, 1), 4.f);

   // Slice:
   // [ 5 ]
   auto s2 = x.Slice({{1, 2}, {2, 3}});
   EXPECT_EQ(s2.GetSize(), 1u);
   EXPECT_EQ(s2.GetShape().size(), 1u);
   EXPECT_EQ(s2.GetShape()[0], 1u);
   EXPECT_EQ(s2(0), 5.f);

   // Slice:
   // [ [ 0, 1, 2 ], [ 3, 4, 5 ] ]
   auto s3 = x.Slice({{0, 2}, {0, 3}});
   EXPECT_EQ(s3.GetSize(), 6u);
   EXPECT_EQ(s3.GetShape().size(), 2u);
   EXPECT_EQ(s3.GetShape()[0], 2u);
   EXPECT_EQ(s3.GetShape()[1], 3u);
   for(auto i = 0; i < 2; i ++) {
      for(auto j = 0; j < 3; j ++) {
         EXPECT_EQ(x(i, j), s3(i, j));
      }
   }
}

TEST(RTensor, SliceColumnMajor)
{
   // Data layout:
   // [ 0, 3, 1, 4, 2, 5  ]
   RTensor<float> x({2, 3}, MemoryLayout::ColumnMajor);
   float c = 0.f;
   for (auto i = 0; i < 2; i++) {
      for (auto j = 0; j < 3; j++) {
         x(i, j) = c;
         c++;
      }
   }

   // Slice:
   // [ [ 0, 1 ], [ 3, 4 ] ]
   auto s1 = x.Slice({{0, 2}, {0, 2}});
   EXPECT_EQ(s1.GetSize(), 4u);
   EXPECT_EQ(s1.GetShape().size(), 2u);
   EXPECT_EQ(s1.GetShape()[0], 2u);
   EXPECT_EQ(s1.GetShape()[1], 2u);
   EXPECT_EQ(s1(0, 0), 0.f);
   EXPECT_EQ(s1(0, 1), 1.f);
   EXPECT_EQ(s1(1, 0), 3.f);
   EXPECT_EQ(s1(1, 1), 4.f);

   // Slice:
   // [ 5 ]
   auto s2 = x.Slice({{1, 2}, {2, 3}});
   EXPECT_EQ(s2.GetSize(), 1u);
   EXPECT_EQ(s2.GetShape().size(), 1u);
   EXPECT_EQ(s2.GetShape()[0], 1u);
   EXPECT_EQ(s2(0), 5.f);

   // Slice:
   // [ [ 0, 1, 2 ], [ 3, 4, 5 ] ]
   auto s3 = x.Slice({{0, 2}, {0, 3}});
   EXPECT_EQ(s3.GetSize(), 6u);
   EXPECT_EQ(s3.GetShape().size(), 2u);
   EXPECT_EQ(s3.GetShape()[0], 2u);
   EXPECT_EQ(s3.GetShape()[1], 3u);
   for(auto i = 0; i < 3; i ++) {
      for(auto j = 0; j < 2; j ++) {
         EXPECT_EQ(x(i, j), s3(i, j));
      }
   }
}

TEST(RTensor, Concatenate1)
{
   float data1[4] = {0, 1, 2, 3};
   float data2[4] = {4, 5, 6, 7};
   RTensor<float> x(data1, {2, 2});
   RTensor<float> y(data2, {2, 2});
   RTensor<float> z = x.Concatenate(y);
   float * new_data = z.GetData();
   auto shape = z.GetShape();

   EXPECT_EQ(new_data[0], 0);
   EXPECT_EQ(new_data[1], 1);
   EXPECT_EQ(new_data[2], 2);
   EXPECT_EQ(new_data[3], 3);
   EXPECT_EQ(new_data[4], 4);
   EXPECT_EQ(new_data[5], 5);
   EXPECT_EQ(new_data[6], 6);
   EXPECT_EQ(new_data[7], 7);
   EXPECT_EQ(shape[0], 4);
   EXPECT_EQ(shape[1], 2);
}

TEST(RTensor, Concatenate2)
{
   float data1[4] = {0, 1, 2, 3};
   float data2[4] = {4, 5, 6, 7};
   RTensor<float> x(data1, {2, 2});
   RTensor<float> y(data2, {2, 2});
   RTensor<float> z = x.Concatenate(y, 1);
   float * new_data = z.GetData();
   auto shape = z.GetShape();

   EXPECT_EQ(new_data[0], 0);
   EXPECT_EQ(new_data[1], 1);
   EXPECT_EQ(new_data[2], 4);
   EXPECT_EQ(new_data[3], 5);
   EXPECT_EQ(new_data[4], 2);
   EXPECT_EQ(new_data[5], 3);
   EXPECT_EQ(new_data[6], 6);
   EXPECT_EQ(new_data[7], 7);
   EXPECT_EQ(shape[0], 2);
   EXPECT_EQ(shape[1], 4);
}

TEST(RTensor, Flatten)
{
   float data1[4] = {0, 1, 2, 3};
   RTensor<float> x(data1, {2, 2});
   RTensor<float> z = x.Flatten();
   float * new_data = z.GetData();
   auto shape = z.GetShape();

   EXPECT_EQ(new_data[0], 0);
   EXPECT_EQ(new_data[1], 1);
   EXPECT_EQ(new_data[2], 2);
   EXPECT_EQ(new_data[3], 3);
   EXPECT_EQ(shape[0], 4);
}
