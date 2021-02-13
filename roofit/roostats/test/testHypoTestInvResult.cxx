// Author: Stephan Hageboeck, CERN  04/2020
#include "RooStats/HypoTestInverterResult.h"

#include "TFile.h"

#include "gtest/gtest.h"

using namespace RooStats;

/// Test that we can correctly read a HypoTestInverterResult
TEST(HypoTestInvResult, ReadFromFile)
{
  const char* filename = "testHypoTestInvResult_1.root";

  TFile file(filename, "READ");
  ASSERT_TRUE(file.IsOpen());

  HypoTestInverterResult* result;
  file.GetObject("result", result);
  ASSERT_NE(result, nullptr);

  // This just reads members
  EXPECT_NEAR(result->UpperLimit(), 2.4613465, 1.E-7);
  EXPECT_NEAR(result->UpperLimitEstimatedError(), 0.059684301, 2.E-7);

  // This accesses the sampling distribution
  EXPECT_DOUBLE_EQ(result->GetExpectedUpperLimit(0), 1.60988427028569);
  EXPECT_NEAR(result->GetExpectedUpperLimit(1), 2.0901937952514,  1.E-13);
  EXPECT_NEAR(result->GetExpectedUpperLimit(2), 2.79444561874078, 1.E-14);

  // Also test that HypoTestResults are coming back correctly
  HypoTestResult* htr = result->GetResult(1);
  ASSERT_NE(htr, nullptr);
  EXPECT_NEAR(htr->CLs(), 0.819079, 1.E-6);
  EXPECT_NEAR(htr->CLsError(), 0.0188863, 1.E-6);
}
