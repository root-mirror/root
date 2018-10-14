#include "TClass.h"
#include "TInterpreter.h"
#include "TSystem.h"


#include "gtest/gtest.h"
#include "gmock/gmock.h"

// Copied from TFileMergerTests.cxx.
// FIXME: Factor out in a new testing library in ROOT.
namespace {
using testing::internal::GetCapturedStderr;
using testing::internal::CaptureStderr;
using testing::internal::RE;
class ExpectedErrorRAII {
   std::string ExpectedRegex;
   void pop()
   {
      std::string Seen = GetCapturedStderr();
      bool match = RE::FullMatch(Seen, RE(ExpectedRegex));
      EXPECT_TRUE(match);
      if (!match) {
         std::string msg = "Match failed!\nSeen: '" + Seen + "'\nRegex: '" + ExpectedRegex + "'\n";
         GTEST_NONFATAL_FAILURE_(msg.c_str());
      }
   }

public:
   ExpectedErrorRAII(std::string E) : ExpectedRegex(E) { CaptureStderr(); }
   ~ExpectedErrorRAII() { pop(); }
};
}

#define EXPECT_ROOT_ERROR(expression, expected_error) \
   {                                                  \
      ExpectedErrorRAII EE(expected_error);           \
      expression;                                     \
   }

// FIXME: We should probably have such a facility in TCling.
static void cleanup()
{
   // Remove AutoDict
   void* dir = gSystem->OpenDirectory(gSystem->pwd());
   const char* name = 0;
   while ((name = gSystem->GetDirEntry(dir)))
      if (!strncmp(name, "AutoDict_", 9))
         gSystem->Unlink(name);

   gSystem->FreeDirectory(dir);
}

class TClingTests : public ::testing::Test {
protected:
   // virtual void SetUp() { }

   // FIXME: We cannot rely on TearDown because it is executed at the end of
   // every test. This triggers another bug in the dictionary generation phase,
   // possibly due to concurrent file system operations.
   //virtual void TearDown() {
      // If there are failures we want to keep the created files.
      //if (!::testing::Test::HasFatalFailure())
      //   cleanup();
   //}

};

// FIXME: Merge with TearDown.
struct CleanupRAII {
   CleanupRAII() {
      if (!::testing::Test::HasFatalFailure())
         cleanup();
   }
} Cleanup;

TEST_F(TClingTests, GenerateDictionaryErrorHandling)
{
   // Check error reporting and handling.
   EXPECT_ROOT_ERROR(ASSERT_FALSE(gInterpreter->GenerateDictionary("", "")),
                     "Error in .* Cannot generate dictionary without passing classes.\n");
   EXPECT_ROOT_ERROR(ASSERT_FALSE(gInterpreter->GenerateDictionary(nullptr, nullptr)),
                     "Error in .* Cannot generate dictionary without passing classes.\n");
}

TEST_F(TClingTests, GenerateDictionaryRegression)
{
   // Make sure we do not crash or go in an infinite loop.
   ASSERT_TRUE(gInterpreter->GenerateDictionary("std::set<int>"));
   ASSERT_TRUE(gInterpreter->GenerateDictionary("std::set<int>", ""));
   ASSERT_TRUE(gInterpreter->GenerateDictionary("std::set<int>", "set"));

   // FIXME: This makes the linkdef parser go in an infinite loop.
   //ASSERT_TRUE(gInterpreter->GenerateDictionary("std::vector<std::array<int, 5>>", ""));
}

TEST_F(TClingTests, GenerateDictionary)
{
   auto cl = TClass::GetClass("vector<TNamed*>");
   ASSERT_FALSE(cl && cl->IsLoaded());

   ASSERT_TRUE(gInterpreter->GenerateDictionary("std::vector<TNamed*>"));
   cl = TClass::GetClass("vector<TNamed*>");
   ASSERT_TRUE(cl != nullptr);
   ASSERT_TRUE(cl->IsLoaded());
}

// Test ROOT-6967
TEST_F(TClingTests, GetEnumWithSameVariableName)
{
   gInterpreter->ProcessLine("int en;enum en{kNone};");
   auto en = gInterpreter->GetEnum(nullptr, "en");
   ASSERT_TRUE(en != nullptr);
}

// Check if we can get the source code of function definitions.
TEST_F(TClingTests, MakeInterpreterValue)
{
   gInterpreter->Declare("void my_func_to_print() {}");
   std::unique_ptr<TInterpreterValue> v = gInterpreter->MakeInterpreterValue();
   gInterpreter->Evaluate("my_func_to_print", *v);
   ASSERT_THAT(v->ToString(), testing::HasSubstr("void my_func_to_print"));
}
