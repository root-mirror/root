#include "gtest/gtest.h"

#include "ROOTUnitTestSupport.h"

#include "TString.h"

TEST(TString, Basics)
{
   TString n("Test", -5);
   EXPECT_STREQ("", n);
   ROOT_EXPECT_ERROR(s = new TString("Test", -5), "TString::TString", "Negative length!");
   TString p("Test", 1);
   EXPECT_STREQ("T", p);
   TString a = "test";
   a.Append("s", -5);
   EXPECT_STREQ("test", a);
   ROOT_EXPECT_ERROR(a.Append("s", -5, "TString::Append", "Negative length!");
}
