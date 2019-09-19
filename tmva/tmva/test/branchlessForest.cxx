#include <gtest/gtest.h>

#include "BDTHelpers.hxx"

#include "TMVA/TreeInference/Forest.hxx"
#include "TMVA/TreeInference/BranchlessTree.hxx"
#include "TMVA/TreeInference/Objectives.hxx"

#include <iostream>
#include <vector>

using namespace TMVA::Experimental;

TEST(BranchlessTree, InferenceFullTreeDepth0)
{
   BranchlessTree<float> tree;
   tree.fTreeDepth = 0;
   tree.fThresholds = {-1.0};
   tree.fInputs = {};
   EXPECT_FLOAT_EQ(tree.Inference(nullptr, 1), -1.0);
}

TEST(BranchlessTree, InferenceFullTreeDepth1)
{
   BranchlessTree<float> tree;
   tree.fTreeDepth = 1;
   tree.fThresholds = {0.0, 1.0, -1.0};
   tree.fInputs = {0};
   float input[1] = {999.0};
   EXPECT_FLOAT_EQ(tree.Inference(input, 1), -1.0);
}

TEST(BranchlessTree, InferenceSparseTreeDepth1)
{
   BranchlessTree<float> tree;
   tree.fTreeDepth = 1;
   tree.fThresholds = {1.0, 0.0, 0.0};
   tree.fInputs = {-1};
   tree.FillSparse();
   EXPECT_FLOAT_EQ(tree.fInputs[0], 0.0);
   EXPECT_FLOAT_EQ(tree.fThresholds[1], 1.0);
   EXPECT_FLOAT_EQ(tree.fThresholds[2], 1.0);
   float input0[1] = {-999.0};
   EXPECT_FLOAT_EQ(tree.Inference(input0, 1), 1.0);
   float input1[1] = {999.0};
   EXPECT_FLOAT_EQ(tree.Inference(input1, 1), 1.0);
}

TEST(BranchlessTree, InferenceFullTreeDepth2)
{
   BranchlessTree<float> tree;
   tree.fTreeDepth = 2;
   tree.fThresholds = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
   tree.fInputs = {0, 1, 2};
   float input0[3] = {-1.0, 0.0, -999.0};
   EXPECT_FLOAT_EQ(tree.Inference(input0, 1), 3.0);
   float input1[3] = {-1.0, 2.0, -999.0};
   EXPECT_FLOAT_EQ(tree.Inference(input1, 1), 4.0);
   float input2[3] = {1.0, -999.0, 1.0};
   EXPECT_FLOAT_EQ(tree.Inference(input2, 1), 5.0);
   float input3[3] = {1.0, -999.0, 3.0};
   EXPECT_FLOAT_EQ(tree.Inference(input3, 1), 6.0);
}

TEST(BranchlessTree, InferenceSparseTreeDepth2)
{
   BranchlessTree<float> tree;
   tree.fTreeDepth = 2;
   tree.fThresholds = {0.0, 1.0, 2.0, 0.0, 0.0, 5.0, 6.0};
   tree.fInputs = {0, -1, 2};
   tree.FillSparse();
   EXPECT_FLOAT_EQ(tree.fInputs[1], 0.0);
   EXPECT_FLOAT_EQ(tree.fThresholds[3], 1.0);
   EXPECT_FLOAT_EQ(tree.fThresholds[4], 1.0);
   float input0[3] = {-1.0, 0.0, -999.0};
   EXPECT_FLOAT_EQ(tree.Inference(input0, 1), 1.0);
   float input1[3] = {-1.0, 2.0, -999.0};
   EXPECT_FLOAT_EQ(tree.Inference(input1, 1), 1.0);
   float input2[3] = {1.0, -999.0, 1.0};
   EXPECT_FLOAT_EQ(tree.Inference(input2, 1), 5.0);
   float input3[3] = {1.0, -999.0, 3.0};
   EXPECT_FLOAT_EQ(tree.Inference(input3, 1), 6.0);
}

TEST(BranchlessForest, InferenceSingleTree)
{
   const auto maxDepth = 1;
   const auto numInputs = 1;
   const auto numTrees = 1;
   WriteModel("myModel", "TestBranchlessForest0.root", "identity", {0}, {0}, {0.0, 1.0, -1.0}, {maxDepth}, {numTrees},
              {numInputs}, {1});

   BranchlessForest<float> forest;
   forest.Load("myModel", "TestBranchlessForest0.root", 0);

   const auto rows = 2;
   float inputs[numInputs * rows] = {-999.0, 999.0};
   float predictions[rows];
   forest.Inference(inputs, rows, true, predictions);
   EXPECT_FLOAT_EQ(predictions[0], 1.0);
   EXPECT_FLOAT_EQ(predictions[1], -1.0);
}

TEST(BranchlessForest, InferenceSingleTreeObjectiveLogistic)
{
   const auto maxDepth = 1;
   const auto numInputs = 1;
   const auto numTrees = 1;
   WriteModel("myModel", "TestBranchlessForest1.root", "logistic", {0}, {0}, {0.0, 1.0, -1.0}, {maxDepth}, {numTrees},
              {numInputs}, {1});

   BranchlessForest<float> forest;
   forest.Load("myModel", "TestBranchlessForest1.root", 0);

   const auto rows = 2;
   float inputs[numInputs * rows] = {-999.0, 999.0};
   float predictions[rows];
   forest.Inference(inputs, rows, true, predictions);
   EXPECT_FLOAT_EQ(predictions[0], Objectives::Logistic<float>(1.0));
   EXPECT_FLOAT_EQ(predictions[1], Objectives::Logistic<float>(-1.0));
}

TEST(BranchlessForest, InferenceTwoTrees)
{
   const auto maxDepth = 1;
   const auto numInputs = 2;
   const auto numTrees = 2;
   const auto outputNode = 1;
   WriteModel("myModel", "TestBranchlessForest2.root", "identity", {0, 1}, {outputNode, outputNode},
              {0.0, 1.0, -1.0, 0.0, 2.0, -2.0}, {maxDepth}, {numTrees}, {numInputs}, {1});

   BranchlessForest<float> forest;
   forest.Load("myModel", "TestBranchlessForest2.root", outputNode);

   const auto rows = 2;
   float inputs[numInputs * rows] = {-999.0, 999.0, 999.0, -999.0};
   float predictions[rows];
   forest.Inference(inputs, rows, true, predictions);
   EXPECT_FLOAT_EQ(predictions[0], 1.0 + -2.0);
   EXPECT_FLOAT_EQ(predictions[1], -1.0 + 2.0);
}

TEST(BranchlessForest, SortTrees)
{
   const auto maxDepth = 1;
   const auto numInputs = 2;
   const auto numTrees = 3;
   WriteModel("myModel", "TestBranchlessForest3.root", "identity", {1, 0, 0}, {0, 0, 0},
              {0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 1.0, 1.0, 1.0}, {maxDepth}, {numTrees}, {numInputs}, {1});

   BranchlessForest<float> forest;
   forest.Load("myModel", "TestBranchlessForest3.root", 0, false);
   EXPECT_EQ(forest.fTrees[0].fInputs[0], 1);
   EXPECT_EQ(forest.fTrees[1].fInputs[0], 0);
   EXPECT_EQ(forest.fTrees[2].fInputs[0], 0);
   EXPECT_EQ(forest.fTrees[0].fThresholds[0], 0.0);
   EXPECT_EQ(forest.fTrees[1].fThresholds[0], 2.0);
   EXPECT_EQ(forest.fTrees[2].fThresholds[0], 1.0);

   BranchlessForest<float> forest2;
   forest2.Load("myModel", "TestBranchlessForest3.root", 0, true);
   EXPECT_EQ(forest2.fTrees[0].fInputs[0], 0);
   EXPECT_EQ(forest2.fTrees[1].fInputs[0], 0);
   EXPECT_EQ(forest2.fTrees[2].fInputs[0], 1);
   EXPECT_EQ(forest2.fTrees[0].fThresholds[0], 1.0);
   EXPECT_EQ(forest2.fTrees[1].fThresholds[0], 2.0);
   EXPECT_EQ(forest2.fTrees[2].fThresholds[0], 0.0);
}
