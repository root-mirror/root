/// \file FoundationUtilsTests.h
///
/// \brief The file contain unit tests which test the ROOT::FoundationUtils
///
/// \author Vassil Vassilev <vvasilev@cern.ch>
///
/// \date Jun, 2020
///
/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/FoundationUtils.hxx>

#include "gtest/gtest.h"

using namespace ROOT::FoundationUtils;

TEST(FoundationUtilsTests, CanConvertEnvValueToBool)
{
  ASSERT_FALSE(CanConvertEnvValueToBool(""));
  ASSERT_TRUE(CanConvertEnvValueToBool("0"));
  ASSERT_TRUE(CanConvertEnvValueToBool("false"));
  ASSERT_TRUE(CanConvertEnvValueToBool("False"));
  ASSERT_TRUE(CanConvertEnvValueToBool("FALSE"));
  ASSERT_TRUE(CanConvertEnvValueToBool("Off"));
  ASSERT_TRUE(CanConvertEnvValueToBool("off"));
  ASSERT_TRUE(CanConvertEnvValueToBool("OFF"));

  ASSERT_TRUE(CanConvertEnvValueToBool("1"));
  ASSERT_TRUE(CanConvertEnvValueToBool("true"));
  ASSERT_TRUE(CanConvertEnvValueToBool("True"));
  ASSERT_TRUE(CanConvertEnvValueToBool("TRUE"));
  ASSERT_TRUE(CanConvertEnvValueToBool("On"));
  ASSERT_TRUE(CanConvertEnvValueToBool("on"));
  ASSERT_TRUE(CanConvertEnvValueToBool("ON"));
}

TEST(FoundationUtilsTests, ConvertEnvValueToBool)
{
  ASSERT_TRUE(ConvertEnvValueToBool("1"));
  ASSERT_TRUE(ConvertEnvValueToBool("TruE"));
  ASSERT_TRUE(ConvertEnvValueToBool("oN"));

  ASSERT_FALSE(ConvertEnvValueToBool("0"));
  ASSERT_FALSE(ConvertEnvValueToBool("FalSe"));
  ASSERT_FALSE(ConvertEnvValueToBool("oFf"));
}
