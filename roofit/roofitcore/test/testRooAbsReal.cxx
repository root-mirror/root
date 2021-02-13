// Tests for RooAbsReal
// Author: Stephan Hageboeck, CERN 05/2020

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooHelpers.h"
#include "RooGlobalFunc.h"

#include "TTree.h"
#include "TFile.h"
#include "gtest/gtest.h"

#include <memory>

// ROOT-6882: Cannot read from ULong64_t branches.
TEST(RooAbsReal, ReadFromTree)
{
  RooRealVar x("ULong64Branch", "xx", 0, 0, 10);
  RooRealVar y("FloatBranch", "yy", 2, 0, 10);
  RooRealVar z("IntBranch", "zz", 2, 0, 10);
  RooRealVar a("DoubleBranch", "aa", 2, 0, 10);

  TFile theFile("testRooAbsReal_1.root", "READ");
  ASSERT_TRUE(theFile.IsOpen());

  TTree* tree = nullptr;
  theFile.GetObject("tree", tree);
  ASSERT_NE(tree, nullptr);

  tree->AddFriend("tree2", "testRooAbsReal_2.root");


  RooHelpers::HijackMessageStream hijack(RooFit::INFO, RooFit::DataHandling);

  RooDataSet data("data", "data", RooArgSet(x,y,z,a), RooFit::Import(*tree));
  ASSERT_EQ(data.numEntries(), 4);

  EXPECT_DOUBLE_EQ(data.moment(x, 1), (1.+3.+3.+7.)/4.)       << "ULong64Branch should contain {1, 3, 3, 7}";
  EXPECT_FLOAT_EQ (data.moment(y, 1), (1.f+0.3f+0.03f+0.007f)/4.f) << "FloatBranch should contain {1., 0.3, 0.03, 0.007}";
  EXPECT_DOUBLE_EQ(data.moment(z, 1), (1.+3.+3.+7.)/4.)       << "IntBranch should contain {1, 3, 3, 7}";
  EXPECT_DOUBLE_EQ(data.moment(a, 1), (1.+0.3+0.03+0.007)/4.) << "DoubleBranch should contain {1, 3, 3, 7}";

  const std::string& msgs = hijack.str();
  const char* targetMsg = " will be converted to double precision";
  EXPECT_NE(msgs.find(std::string(x.GetName()) + targetMsg), std::string::npos) << "Expect to see INFO messages for conversion of integer branch to double.";
  EXPECT_NE(msgs.find(std::string(y.GetName()) + targetMsg), std::string::npos) << "Expect to see INFO messages for conversion of float branch to double.";
  EXPECT_NE(msgs.find(std::string(z.GetName()) + targetMsg), std::string::npos) << "Expect to see INFO messages for conversion of integer branch to double.";
  EXPECT_EQ(msgs.find(std::string(a.GetName()) + targetMsg), std::string::npos) << "Expect not to see INFO messages for conversion of double branch to double.";
}

