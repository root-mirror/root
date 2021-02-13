#include "ROOT/RLogger.hxx"

#include <ROOT/RMakeUnique.hxx>
#include <TError.h>

#include "gtest/gtest.h"

#include <vector>

using namespace ROOT::Experimental;

struct TestLogger {
   class Handler : public RLogHandler {
      TestLogger *fLogger;

   public:
      Handler(TestLogger &logger) : fLogger(&logger) {}

      bool Emit(const RLogEntry &entry) override
      {
         fLogger->fLogEntries.emplace_back(entry);
         // Stop emission.
         return false;
      }
   };

   std::vector<RLogEntry> fLogEntries;
   Handler *fHandler;

   TestLogger()
   {
      auto uptr = std::make_unique<Handler>(*this);
      fHandler = uptr.get();
      RLogManager::Get().PushFront(std::move(uptr));
   }

   ~TestLogger() { RLogManager::Get().Remove(fHandler); }

   size_t size() const { return fLogEntries.size(); }
   bool empty() const { return fLogEntries.empty(); }
};

TEST(Logger, EmittedEntry)
{
   TestLogger testLogger;
   auto prevErrors = RLogManager::Get().GetNumErrors();
   auto prevWarnings = RLogManager::Get().GetNumWarnings();

   RLogChannel channel("TheChannel");

   // clang-format off
   R__LOG_ERROR(channel) << "The text"; auto logLine = __LINE__;
   // clang-format on

   // Check manager's emission
   EXPECT_EQ(testLogger.size(), 1);
   EXPECT_EQ(RLogManager::Get().GetNumErrors(), prevErrors + 1);
   EXPECT_EQ(RLogManager::Get().GetNumWarnings(), prevWarnings);

   // Check emitted RLogEntry
   EXPECT_EQ(testLogger.fLogEntries[0].fChannel->GetName(), "TheChannel");
   EXPECT_NE(testLogger.fLogEntries[0].fLocation.fFile.find("testLogger.cxx"), std::string::npos);
   EXPECT_NE(testLogger.fLogEntries[0].fLocation.fFuncName.find("EmittedEntry"), std::string::npos);
   EXPECT_EQ(testLogger.fLogEntries[0].fLocation.fLine, logLine);
   EXPECT_EQ(testLogger.fLogEntries[0].fLevel, ELogLevel::kError);
   EXPECT_EQ(testLogger.fLogEntries[0].fMessage, "The text");
   EXPECT_TRUE(testLogger.fLogEntries[0].IsError());
   EXPECT_FALSE(testLogger.fLogEntries[0].IsWarning());
}

TEST(Logger, RLogManagerCounts)
{
   TestLogger testLogger;
   // Check diag counts of RLogManager.
   auto initialWarnings = RLogManager::Get().GetNumWarnings();
   auto initialErrors = RLogManager::Get().GetNumErrors();
   R__LOG_ERROR() << "emitted";
   EXPECT_EQ(RLogManager::Get().GetNumWarnings(), initialWarnings);
   EXPECT_EQ(RLogManager::Get().GetNumErrors(), initialErrors + 1);
}

TEST(Logger, RLogDiagCounter)
{
   TestLogger testLogger;
   R__LOG_ERROR() << "emitted"; // before RAII, should not be counted.
   R__LOG_ERROR() << "emitted"; // before RAII, should not be counted.
   // Check counter seeing what was emitted during its lifetime.
   RLogScopedDiagCount counter;
   EXPECT_EQ(counter.GetAccumulatedWarnings(), 0);
   EXPECT_EQ(counter.GetAccumulatedErrors(), 0);
   EXPECT_FALSE(counter.HasWarningOccurred());
   EXPECT_FALSE(counter.HasErrorOccurred());
   EXPECT_FALSE(counter.HasErrorOrWarningOccurred());

   R__LOG_ERROR() << "emitted";

   EXPECT_EQ(counter.GetAccumulatedWarnings(), 0);
   EXPECT_EQ(counter.GetAccumulatedErrors(), 1);
   EXPECT_FALSE(counter.HasWarningOccurred());
   EXPECT_TRUE(counter.HasErrorOccurred());
   EXPECT_TRUE(counter.HasErrorOrWarningOccurred());
}

TEST(Logger, RLogScopedVerbositySuppress)
{
   RLogChannel channel("ABC");

   auto initialLogLevel = RLogManager::Get().GetVerbosity();
   auto initialLogLevelABC = channel.GetVerbosity();
   EXPECT_EQ(initialLogLevel, ELogLevel::kWarning);
   EXPECT_EQ(initialLogLevelABC, ELogLevel::kUnset);

   {
      // Test simple suppression.
      TestLogger testLogger;
      R__LOG_WARNING(channel) << "emitted";
      EXPECT_EQ(testLogger.size(), 1);
      RLogScopedVerbosity suppress(ELogLevel::kError);
      R__LOG_WARNING(channel) << "suppressed";
      EXPECT_EQ(testLogger.size(), 1);
   }
   {
      // Test channel specific suppression given global higher verbosity.
      TestLogger testLogger;
      RLogScopedVerbosity suppressGlobal(ELogLevel::kInfo);
      RLogScopedVerbosity suppress(channel, ELogLevel::kError);
      R__LOG_WARNING(channel) << "suppressed";
      R__LOG_INFO(channel) << "suppressed, too";
      EXPECT_TRUE(testLogger.empty());
      R__LOG_ERROR(channel) << "emitted";
      EXPECT_EQ(testLogger.size(), 1);

      R__LOG_ERROR() << "emitted";
      R__LOG_WARNING() << "emitted, second";
      R__LOG_INFO() << "emitted, third";
      EXPECT_EQ(testLogger.size(), 4);
   }
   {
      // Check unrelated channel.
      TestLogger testLogger;
      RLogChannel channel123("123");
      RLogScopedVerbosity suppress(channel123, ELogLevel::kFatal);
      R__LOG_ERROR(channel) << "emitted";
      EXPECT_EQ(testLogger.size(), 1);
   }
   {
      // Check global vs specific suppression.
      TestLogger testLogger;
      RLogScopedVerbosity suppressGlobal(ELogLevel::kDebug);
      RLogScopedVerbosity suppress(channel, ELogLevel::kError);
      R__LOG_WARNING(channel) << "suppressed";
      EXPECT_TRUE(testLogger.empty());
      R__LOG_WARNING() << "emitted";
      EXPECT_EQ(testLogger.size(), 1);
   }

   // Check that levels have returned to values before RAII.
   EXPECT_EQ(RLogManager::Get().GetVerbosity(), initialLogLevel);
   EXPECT_EQ(channel.GetVerbosity(), initialLogLevelABC);
}

TEST(Logger, RLogScopedVerbosityVerbose)
{
   RLogChannel channel("ABC");

   auto initialLogLevel = RLogManager::Get().GetVerbosity();
   auto initialLogLevelABC = channel.GetVerbosity();
   EXPECT_EQ(initialLogLevel, ELogLevel::kWarning);
   EXPECT_EQ(initialLogLevelABC, ELogLevel::kUnset);

   {
      // Test same diag level as verbosity, in channel and global, before and after RAII.
      TestLogger testLogger;
      R__LOG_DEBUG(0, channel) << "suppressed";
      EXPECT_TRUE(testLogger.empty());
      R__LOG_DEBUG(0) << "suppressed";
      EXPECT_TRUE(testLogger.empty());
      RLogScopedVerbosity verbose(ELogLevel::kDebug);
      R__LOG_DEBUG(0, channel) << "emitted";
      EXPECT_EQ(testLogger.size(), 1);
      R__LOG_DEBUG(0) << "emitted";
      EXPECT_EQ(testLogger.size(), 2);
   }
   {
      // Test different diag levels, in channel and global, before and after RAII.
      TestLogger testLogger;
      R__LOG_DEBUG(0, channel) << "suppressed";
      EXPECT_TRUE(testLogger.empty());
      R__LOG_WARNING(channel) << "emitted";
      EXPECT_EQ(testLogger.size(), 1);
      R__LOG_DEBUG(0) << "suppressed";
      EXPECT_EQ(testLogger.size(), 1);
      R__LOG_WARNING() << "emitted";
      EXPECT_EQ(testLogger.size(), 2);
      RLogScopedVerbosity verbose(channel, ELogLevel::kDebug);
      R__LOG_DEBUG(0, channel) << "emitted";
      EXPECT_EQ(testLogger.size(), 3);
      R__LOG_WARNING(channel) << "emitted, too";
      EXPECT_EQ(testLogger.size(), 4);
      R__LOG_DEBUG(0) << "suppressed";
      EXPECT_EQ(testLogger.size(), 4);
      R__LOG_WARNING() << "emitted";
      EXPECT_EQ(testLogger.size(), 5);
   }
   {
      // Test Info level verbosity.
      TestLogger testLogger;
      R__LOG_INFO(channel) << "suppressed";
      EXPECT_TRUE(testLogger.empty());
      R__LOG_DEBUG(0, channel) << "suppressed, second";
      EXPECT_TRUE(testLogger.empty());
      RLogScopedVerbosity verbose(channel, ELogLevel::kInfo);
      R__LOG_INFO(channel) << "emitted";
      EXPECT_EQ(testLogger.size(), 1);
      R__LOG_DEBUG(0, channel) << "suppressed, third";
      EXPECT_EQ(testLogger.size(), 1);
   }
   {
      // Test verbosity change on other channel not influcing this one.
      TestLogger testLogger;
      RLogChannel otherChannel("123");
      RLogScopedVerbosity verbose(otherChannel, ELogLevel::kDebug);
      R__LOG_DEBUG(0, channel) << "suppressed";
      EXPECT_TRUE(testLogger.empty());
      R__LOG_DEBUG(0, otherChannel) << "emitted";
      EXPECT_EQ(testLogger.size(), 1);
   }

   // Check that levels have returned to values before RAII.
   EXPECT_EQ(RLogManager::Get().GetVerbosity(), initialLogLevel);
   EXPECT_EQ(channel.GetVerbosity(), initialLogLevelABC);
}

TEST(Logger, ExtraVerbosityLevels)
{
   TestLogger testLogger;
   RLogChannel channel("channel");
   RLogScopedVerbosity verbose(channel, ELogLevel::kDebug + 50);

   R__LOG_DEBUG(0, channel) << "emitted";
   EXPECT_EQ(testLogger.size(), 1);
   R__LOG_DEBUG(0, channel) << "emitted";
   EXPECT_EQ(testLogger.size(), 2);
   R__LOG_DEBUG(0, channel) << "emitted";
   EXPECT_EQ(testLogger.size(), 3);
   R__LOG_DEBUG(50, channel) << "emitted";
   EXPECT_EQ(testLogger.size(), 4);
   R__LOG_DEBUG(51, channel) << "suppressed";
   EXPECT_EQ(testLogger.size(), 4);
   R__LOG_DEBUG(50) << "suppressed";
   EXPECT_EQ(testLogger.size(), 4);
}

TEST(Logger, SuppressStreamEval)
{
   TestLogger testLogger;
   RLogChannel channel("channel");
   bool wasEvaluated = false;
   R__LOG_DEBUG(0, channel) << "It's debug, this should not be called!" << [&]() -> int {
      wasEvaluated = true;
      return 0;
   }();
   EXPECT_FALSE(wasEvaluated);
}

namespace {
struct TestErrorHandler_t {
   int fLevel;
   Bool_t fAbort;
   std::string fLocation;
   std::string fMsg;
};

TestErrorHandler_t testErrorHandlerVal;
void testErrorHandler(int level, Bool_t abort, const char *location, const char *msg)
{
   testErrorHandlerVal.fLevel = level;
   testErrorHandlerVal.fAbort = abort;
   testErrorHandlerVal.fLocation = location;
   testErrorHandlerVal.fMsg = msg;
}
} // namespace

TEST(Logger, ROOTErrorHandlerDiagString)
{
   auto prevErrorHandler = GetErrorHandler();
   SetErrorHandler(testErrorHandler);
   RLogChannel channel("ROOT.checkChannelName");

   // clang-format off
   R__LOG_ERROR(channel) << "check message " << 42; auto lineNumber = __LINE__;
   // clang-format on

   EXPECT_EQ(testErrorHandlerVal.fLevel, kError);
   EXPECT_EQ(testErrorHandlerVal.fAbort, false);
   EXPECT_NE(testErrorHandlerVal.fLocation.find("[ROOT.checkChannelName]"), std::string::npos)
      << "testErrorHandlerVal.fLocation is " << testErrorHandlerVal.fLocation;
   EXPECT_NE(testErrorHandlerVal.fLocation.find(" Error "), std::string::npos)
      << "testErrorHandlerVal.fLocation is " << testErrorHandlerVal.fLocation;
   EXPECT_NE(testErrorHandlerVal.fLocation.find("testLogger.cxx:" + std::to_string(lineNumber) + " "),
             std::string::npos)
      << "testErrorHandlerVal.fLocation is " << testErrorHandlerVal.fLocation;
   EXPECT_NE(testErrorHandlerVal.fLocation.find("ROOTErrorHandlerDiagString"), std::string::npos)
      << "testErrorHandlerVal.fLocation is " << testErrorHandlerVal.fLocation;

   EXPECT_EQ(testErrorHandlerVal.fMsg, "check message 42");

   SetErrorHandler(prevErrorHandler);
}
