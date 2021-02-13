/// \file ROOT/RLogger.hxx
/// \ingroup Base ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2015-03-29
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RLogger
#define ROOT7_RLogger

#include <atomic>
#include <list>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <utility>

namespace ROOT {
namespace Experimental {

class RLogEntry;
class RLogManager;

/**
 Kinds of diagnostics.
 */
enum class ELogLevel : unsigned char {
   kUnset,
   kFatal,   ///< An error which causes further processing to be unreliable
   kError,   ///< An error
   kWarning, ///< Warnings about likely unexpected behavior
   kInfo,    ///< Informational messages; used for instance for tracing
   kDebug    ///< Debug information; only useful for developers; can have added verbosity up to 255-kDebug.
};

inline ELogLevel operator+(ELogLevel severity, int offset)
{
   return static_cast<ELogLevel>(static_cast<int>(severity) + offset);
}

/**
 Keep track of emitted errors and warnings.
 */
class RLogDiagCount {
protected:
   std::atomic<long long> fNumWarnings{0ll};    /// Number of warnings.
   std::atomic<long long> fNumErrors{0ll};      /// Number of errors.
   std::atomic<long long> fNumFatalErrors{0ll}; /// Number of fatal errors.

public:
   /// Returns the current number of warnings.
   long long GetNumWarnings() const { return fNumWarnings; }

   /// Returns the current number of errors.
   long long GetNumErrors() const { return fNumErrors; }

   /// Returns the current number of fatal errors.
   long long GetNumFatalErrors() const { return fNumFatalErrors; }

   /// Increase warning or error count.
   void Increment(ELogLevel severity)
   {
      switch (severity) {
      case ELogLevel::kFatal: ++fNumFatalErrors; break;
      case ELogLevel::kError: ++fNumErrors; break;
      case ELogLevel::kWarning: ++fNumWarnings; break;
      default:;
      }
   }
};

/**
 Abstract RLogHandler base class. ROOT logs everything from info to error
 to entities of this class.
 */
class RLogHandler {
public:
   virtual ~RLogHandler();
   /// Emit a log entry.
   /// \param entry - the RLogEntry to be emitted.
   /// \returns false if further emission of this Log should be suppressed.
   ///
   /// \note This function is called concurrently; log emission must be locked
   /// if needed. (The default log handler using ROOT's DefaultErrorHandler is locked.)
   virtual bool Emit(const RLogEntry &entry) = 0;
};

/**
 A log configuration for a channel, e.g. "RHist".
 Each ROOT module has its own log, with potentially distinct verbosity.
 */
class RLogChannel : public RLogDiagCount {
   /// Name as shown in diagnostics
   std::string fName;

   /// Verbosity of this channel. By default, use the global verbosity.
   ELogLevel fVerbosity = ELogLevel::kUnset;

public:
   /// Construct an anonymous channel.
   RLogChannel() = default;

   /// Construct an anonymous channel with a default vebosity.
   explicit RLogChannel(ELogLevel verbosity) : fVerbosity(verbosity) {}

   /// Construct a log channel given its name, which is part of the diagnostics.
   RLogChannel(const std::string &name) : fName(name) {}

   ELogLevel SetVerbosity(ELogLevel verbosity)
   {
      std::swap(fVerbosity, verbosity);
      return verbosity;
   }
   ELogLevel GetVerbosity() const { return fVerbosity; }
   ELogLevel GetEffectiveVerbosity(const RLogManager &mgr) const;

   const std::string &GetName() const { return fName; }
};

/**
 A RLogHandler that multiplexes diagnostics to different client `RLogHandler`s
 and keeps track of the sum of `RLogDiagCount`s for all channels.

 `RLogHandler::Get()` returns the process's (static) log manager.
 */

class RLogManager : public RLogChannel, public RLogHandler {
   std::mutex fMutex;
   std::list<std::unique_ptr<RLogHandler>> fHandlers;

public:
   /// Initialize taking a RLogHandler.
   RLogManager(std::unique_ptr<RLogHandler> lh) : RLogChannel(ELogLevel::kWarning)
   {
      fHandlers.emplace_back(std::move(lh));
   }

   static RLogManager &Get();

   /// Add a RLogHandler in the front - to be called before all others.
   void PushFront(std::unique_ptr<RLogHandler> handler) { fHandlers.emplace_front(std::move(handler)); }

   /// Add a RLogHandler in the back - to be called after all others.
   void PushBack(std::unique_ptr<RLogHandler> handler) { fHandlers.emplace_back(std::move(handler)); }

   /// Remove and return the given log handler. Returns `nullptr` if not found.
   std::unique_ptr<RLogHandler> Remove(RLogHandler *handler);

   // Emit a `RLogEntry` to the RLogHandlers.
   // Returns false if further emission of this Log should be suppressed.
   bool Emit(const RLogEntry &entry) override;
};

/**
 A diagnostic location, part of an RLogEntry.
 */
struct RLogLocation {
   std::string fFile;
   std::string fFuncName;
   int fLine; // C++11 forbids "= 0" for braced-init-list initialization.
};

/**
 A diagnostic that can be emitted by the RLogManager.
 One can construct a RLogEntry through RLogBuilder, including streaming into
 the diagnostic message and automatic emission.
 */

class RLogEntry {
public:
   RLogLocation fLocation;
   std::string fMessage;
   RLogChannel *fChannel = nullptr;
   ELogLevel fLevel = ELogLevel::kFatal;

   RLogEntry(ELogLevel level, RLogChannel &channel) : fChannel(&channel), fLevel(level) {}
   RLogEntry(ELogLevel level, RLogChannel &channel, const RLogLocation &loc)
      : fLocation(loc), fChannel(&channel), fLevel(level)
   {
   }

   bool IsDebug() const { return fLevel >= ELogLevel::kDebug; }
   bool IsInfo() const { return fLevel == ELogLevel::kInfo; }
   bool IsWarning() const { return fLevel == ELogLevel::kWarning; }
   bool IsError() const { return fLevel == ELogLevel::kError; }
   bool IsFatal() const { return fLevel == ELogLevel::kFatal; }
};

namespace Detail {
/**
 Builds a diagnostic entry, emitted by the static RLogManager upon destruction of this builder,
 where - by definition - the RLogEntry has been completely built.

 This builder can be used through the utility preprocessor macros R__LOG_ERROR,
 R__LOG_WARNING etc like this:
 ```c++
     R__LOG_INFO(ROOT::Experimental::HistLog()) << "all we know is " << 42;
     const int decreasedInfoLevel = 5;
     R__LOG_XDEBUG(ROOT::Experimental::WebGUILog(), decreasedInfoLevel) << "nitty-gritty details";
 ```
 This will automatically capture the current class and function name, the file and line number.
 */

class RLogBuilder : public std::ostringstream {
   /// The log entry to be built.
   RLogEntry fEntry;

public:
   RLogBuilder(ELogLevel level, RLogChannel &channel) : fEntry(level, channel) {}
   RLogBuilder(ELogLevel level, RLogChannel &channel, const std::string &filename, int line,
               const std::string &funcname)
      : fEntry(level, channel, {filename, funcname, line})
   {
   }

   /// Emit the log entry through the static log manager.
   ~RLogBuilder()
   {
      fEntry.fMessage = str();
      RLogManager::Get().Emit(fEntry);
   }
};
} // namespace Detail

/**
 Change the verbosity level (global or specific to the RLogChannel passed to the
 constructor) for the lifetime of this object.
 Example:
 ```c++
 RLogScopedVerbosity debugThis(gFooLog, ELogLevel::kDebug);
 Foo::SomethingToDebug();
 ```
 */
class RLogScopedVerbosity {
   RLogChannel *fChannel;
   ELogLevel fPrevLevel;

public:
   RLogScopedVerbosity(RLogChannel &channel, ELogLevel verbosity)
      : fChannel(&channel), fPrevLevel(channel.SetVerbosity(verbosity))
   {
   }
   explicit RLogScopedVerbosity(ELogLevel verbosity) : RLogScopedVerbosity(RLogManager::Get(), verbosity) {}
   ~RLogScopedVerbosity() { fChannel->SetVerbosity(fPrevLevel); }
};

/**
 Object to count the number of warnings and errors emitted by a section of code,
 after construction of this type.
 */
class RLogScopedDiagCount {
   RLogDiagCount *fCounter = nullptr;
   /// The number of the RLogDiagCount's emitted warnings at construction time of *this.
   long long fInitialWarnings = 0;
   /// The number of the RLogDiagCount's emitted errors at construction time.
   long long fInitialErrors = 0;
   /// The number of the RLogDiagCount's emitted fatal errors at construction time.
   long long fInitialFatalErrors = 0;

public:
   /// Construct the scoped count given a counter (e.g. a channel or RLogManager).
   /// The counter's lifetime must exceed the lifetime of this object!
   explicit RLogScopedDiagCount(RLogDiagCount &cnt)
      : fCounter(&cnt), fInitialWarnings(cnt.GetNumWarnings()), fInitialErrors(cnt.GetNumErrors()),
        fInitialFatalErrors(cnt.GetNumFatalErrors())
   {
   }

   /// Construct the scoped count for any diagnostic, whatever its channel.
   RLogScopedDiagCount() : RLogScopedDiagCount(RLogManager::Get()) {}

   /// Get the number of warnings that the RLogDiagCount has emitted since construction of *this.
   long long GetAccumulatedWarnings() const { return fCounter->GetNumWarnings() - fInitialWarnings; }

   /// Get the number of errors that the RLogDiagCount has emitted since construction of *this.
   long long GetAccumulatedErrors() const { return fCounter->GetNumErrors() - fInitialErrors; }

   /// Get the number of errors that the RLogDiagCount has emitted since construction of *this.
   long long GetAccumulatedFatalErrors() const { return fCounter->GetNumFatalErrors() - fInitialFatalErrors; }

   /// Whether the RLogDiagCount has emitted a warnings since construction time of *this.
   bool HasWarningOccurred() const { return GetAccumulatedWarnings(); }

   /// Whether the RLogDiagCount has emitted an error (fatal or not) since construction time of *this.
   bool HasErrorOccurred() const { return GetAccumulatedErrors() + GetAccumulatedFatalErrors(); }

   /// Whether the RLogDiagCount has emitted an error or a warning since construction time of *this.
   bool HasErrorOrWarningOccurred() const { return HasWarningOccurred() || HasErrorOccurred(); }
};

namespace Internal {

inline RLogChannel &GetChannelOrManager()
{
   return RLogManager::Get();
}
inline RLogChannel &GetChannelOrManager(RLogChannel &channel)
{
   return channel;
}

} // namespace Internal

inline ELogLevel RLogChannel::GetEffectiveVerbosity(const RLogManager &mgr) const
{
   if (fVerbosity == ELogLevel::kUnset)
      return mgr.GetVerbosity();
   return fVerbosity;
}

} // namespace Experimental
} // namespace ROOT

#if defined(_MSC_VER)
#define R__LOG_PRETTY_FUNCTION __FUNCSIG__
#else
#define R__LOG_PRETTY_FUNCTION __PRETTY_FUNCTION__
#endif

/*
 Some implementation details:

 - The conditional `RLogBuilder` use prevents stream operators from being called if
 verbosity is too low, i.e.:
 ````
 RLogScopedVerbosity silence(RLogLevel::kFatal);
 R__LOG_DEBUG(7) << WillNotBeCalled();
 ```
 - To update counts of warnings / errors / fatal errors, those RLogEntries must
 always be created, even if in the end their emission will be silenced. This
 should be fine, performance-wise, as they should not happen frequently.
 - Use `(condition) && RLogBuilder(...)` instead of `if (condition) RLogBuilder(...)`
 to prevent "ambiguous else" in invocations such as `if (something) R__LOG_DEBUG()...`.
 */
#define R__LOG_TO_CHANNEL(SEVERITY, CHANNEL)                                                                        \
   ((SEVERITY < ROOT::Experimental::ELogLevel::kInfo + 0) ||                                                        \
     ROOT::Experimental::Internal::GetChannelOrManager(CHANNEL).GetEffectiveVerbosity(                              \
        ROOT::Experimental::RLogManager::Get()) >= SEVERITY) &&                                                     \
      ROOT::Experimental::Detail::RLogBuilder(SEVERITY, ROOT::Experimental::Internal::GetChannelOrManager(CHANNEL), \
                                              __FILE__, __LINE__, R__LOG_PRETTY_FUNCTION)

/// \name LogMacros
/// Macros to log diagnostics.
/// ```c++
///     R__LOG_INFO(ROOT::Experimental::HistLog()) << "all we know is " << 42;
///
///     RLogScopedVerbosity verbose(kDebug + 5);
///     const int decreasedInfoLevel = 5;
///     R__LOG_DEBUG(decreasedInfoLevel, ROOT::Experimental::WebGUILog()) << "nitty-gritty details";
/// ```
///\{
#define R__LOG_FATAL(...) R__LOG_TO_CHANNEL(ROOT::Experimental::ELogLevel::kFatal, __VA_ARGS__)
#define R__LOG_ERROR(...) R__LOG_TO_CHANNEL(ROOT::Experimental::ELogLevel::kError, __VA_ARGS__)
#define R__LOG_WARNING(...) R__LOG_TO_CHANNEL(ROOT::Experimental::ELogLevel::kWarning, __VA_ARGS__)
#define R__LOG_INFO(...) R__LOG_TO_CHANNEL(ROOT::Experimental::ELogLevel::kInfo, __VA_ARGS__)
#define R__LOG_DEBUG(DEBUGLEVEL, ...) R__LOG_TO_CHANNEL(ROOT::Experimental::ELogLevel::kDebug + DEBUGLEVEL, __VA_ARGS__)
///\}

#endif
