#include <gtest/gtest.h>

#include "BDTHelpers.hxx"
#include "TMVA/RBDT.hxx"

#include "ROOT/RVec.hxx"

#include <cmath>

using namespace TMVA::Experimental;

TEST(RBDT, ClassificationSingleEvent)
{
   const auto maxDepth = 1;
   const auto numInputs = 1;
   const auto numTrees = 1;
   WriteModel("myModel", "TestRBDT0.root", "identity", {0}, {0}, {0.0, 1.0, -1.0}, {maxDepth}, {numTrees}, {numInputs},
              {1});

   RBDT<> bdt("myModel", "TestRBDT0.root");
   auto y = bdt.Compute({-999.0});
   EXPECT_EQ(y.size(), 1u);
   EXPECT_FLOAT_EQ(y[0], 1.0);
   EXPECT_FLOAT_EQ(bdt.Compute({999.0})[0], -1.0);
}

TEST(RBDT, ClassificationSingleEventRVec)
{
   const auto maxDepth = 1;
   const auto numInputs = 1;
   const auto numTrees = 1;
   WriteModel("myModel", "TestRBDT1.root", "identity", {0}, {0}, {0.0, 1.0, -1.0}, {maxDepth}, {numTrees}, {numInputs},
              {1});

   RBDT<> bdt("myModel", "TestRBDT1.root");
   ROOT::RVec<float> x = {-999.0};
   auto y = bdt.Compute(x);
   EXPECT_EQ(y.size(), 1u);
   EXPECT_FLOAT_EQ(y[0], 1.0);
}

TEST(RBDT, ClassificationBatch)
{
   const auto maxDepth = 1;
   const auto numInputs = 1;
   const auto numTrees = 1;
   WriteModel("myModel", "TestRBDT2.root", "identity", {0}, {0}, {0.0, 1.0, -1.0}, {maxDepth}, {numTrees}, {numInputs},
              {1});

   RBDT<> bdt("myModel", "TestRBDT2.root");
   RTensor<float> x({2, 1});
   x(0, 0) = -999.0;
   x(0, 1) = 999.0;
   auto y = bdt.Compute(x);
   const auto shape = y.GetShape();
   EXPECT_EQ(shape[0], 2u);
   EXPECT_EQ(shape[1], 1u);
   EXPECT_FLOAT_EQ(y(0, 0), 1.0);
   EXPECT_FLOAT_EQ(y(1, 0), -1.0);
}

TEST(RBDT, MulticlassSingleEvent)
{
   const auto maxDepth = 1;
   const auto numInputs = 1;
   const auto numOutputs = 3;
   const auto numTrees = 3;
   WriteModel("myModel", "TestRBDT3.root", "softmax", {0, 0, 0}, {0, 1, 2},
              {0.0, 1.0, -1.0, 0.0, -1.0, 1.0, 0.0, 2.0, -2.0}, {maxDepth}, {numTrees}, {numInputs}, {numOutputs});

   RBDT<> bdt("myModel", "TestRBDT3.root");

   auto y = bdt.Compute({-999.0});
   EXPECT_EQ(y.size(), 3u);
   const auto s = std::exp(1.0) + std::exp(-1.0) + std::exp(2.0);
   EXPECT_FLOAT_EQ(y[0], std::exp(1.0) / s);
   EXPECT_FLOAT_EQ(y[1], std::exp(-1.0) / s);
   EXPECT_FLOAT_EQ(y[2], std::exp(2.0) / s);

   auto y2 = bdt.Compute({999.0});
   EXPECT_EQ(y2.size(), 3u);
   const auto s2 = std::exp(-1.0) + std::exp(1.0) + std::exp(-2.0);
   EXPECT_FLOAT_EQ(y2[0], std::exp(-1.0) / s2);
   EXPECT_FLOAT_EQ(y2[1], std::exp(1.0) / s2);
   EXPECT_FLOAT_EQ(y2[2], std::exp(-2.0) / s2);
}

TEST(RBDT, MulticlassBatch)
{
   const auto maxDepth = 1;
   const auto numInputs = 1;
   const auto numOutputs = 3;
   const auto numTrees = 3;
   WriteModel("myModel", "TestRBDT4.root", "identity", {0, 0, 0}, {0, 1, 2},
              {0.0, 1.0, -1.0, 0.0, -1.0, 1.0, 0.0, 2.0, -2.0}, {maxDepth}, {numTrees}, {numInputs}, {numOutputs});

   RBDT<> bdt("myModel", "TestRBDT4.root");
   RTensor<float> x({2, 1});
   x(0, 0) = -999.0;
   x(0, 1) = 999.0;
   auto y = bdt.Compute(x);
   const auto shape = y.GetShape();
   EXPECT_EQ(shape[0], 2u);
   EXPECT_EQ(shape[1], 3u);
   EXPECT_FLOAT_EQ(y(0, 0), 1.0);
   EXPECT_FLOAT_EQ(y(0, 1), -1.0);
   EXPECT_FLOAT_EQ(y(0, 2), 2.0);
   EXPECT_FLOAT_EQ(y(1, 0), -1.0);
   EXPECT_FLOAT_EQ(y(1, 1), 1.0);
   EXPECT_FLOAT_EQ(y(1, 2), -2.0);
}

TEST(RBDT, ColumnMajorInput)
{
   const auto maxDepth = 1;
   const auto numInputs = 2;
   const auto numTrees = 1;
   WriteModel("myModel", "TestRBDT5.root", "identity", {0}, {0}, {0.0, 1.0, -1.0}, {maxDepth}, {numTrees}, {numInputs},
              {1});

   RBDT<> bdt("myModel", "TestRBDT5.root");
   float data[4] = {-999.0, -999.0, 999.0, 999.0};
   RTensor<float> x(data, {2, 2}, MemoryLayout::ColumnMajor);
   auto y = bdt.Compute(x);
   EXPECT_FLOAT_EQ(y(0, 0), 1.0);
   EXPECT_FLOAT_EQ(y(1, 0), 1.0);
}
