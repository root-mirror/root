#include "gtest/gtest.h"

#include "TError.h"
#include "TROOT.h"

#include <cerrno>
#include <string>

int gTestLastLevel = -1;
bool gTestLastAbort = false;
std::string gTestLastLocation;
std::string gTestLastMsg;

void TestErrorHandler(int level, Bool_t abort, const char *location, const char *msg)
{
   gTestLastLevel = level;
   gTestLastAbort = abort;
   gTestLastLocation = location;
   gTestLastMsg = msg;
}

TEST(TError, Basics) {
   ASSERT_TRUE(gROOT);
   SetErrorHandler(TestErrorHandler);

   Info("location", "message");
   EXPECT_STREQ("location", gTestLastLocation.c_str());
   EXPECT_STREQ("message", gTestLastMsg.c_str());
   EXPECT_FALSE(gTestLastAbort);

   errno = 42;
   SysError("location", "message");
   // We expect the explanation for errno 42
   EXPECT_STRNE("message", gTestLastMsg.c_str());

   ROOT::Internal::gROOTLocal->~TROOT();
   // The TROOT destructor re-installed the minimal error handler
   SetErrorHandler(TestErrorHandler);

   errno = 42;
   SysError("location", "message");
   EXPECT_STREQ("message (errno: 42)", gTestLastMsg.c_str());
}


TEST(TError, LongMessage) {
   SetErrorHandler(TestErrorHandler);
   std::string longMessage(10000, 'X');
   EXPECT_EQ(10000U, longMessage.length());

   Info("location", "%s", longMessage.c_str());
   EXPECT_EQ(longMessage, gTestLastMsg);

   // Again, should reuse the now enlarged, thread-local heap-allocated buffer
   Info("location", "%s", longMessage.c_str());
   EXPECT_EQ(longMessage, gTestLastMsg);
}
