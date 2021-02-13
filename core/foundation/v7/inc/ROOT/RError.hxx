/// \file ROOT/RError.hxx
/// \ingroup Base ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2019-12-11
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RError
#define ROOT7_RError

#include <ROOT/RConfig.hxx> // for R__[un]likely
#include <ROOT/RLogger.hxx> // for R__LOG_PRETTY_FUNCTION
#include <ROOT/RMakeUnique.hxx>

#include <cstddef>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// The RResult<T> class and their related classes are used for call chains that can throw exceptions,
// such as I/O code paths.  Throwing of the exception is deferred to allow for `if (result)` style error
// checking where it makes sense.
//
// A function returning an RResult might look like this:
//
//     RResult<int> MyIOFunc()
//     {
//        int rv = syscall(...);
//        if (rv == -1)
//           return R__FAIL("user-facing error message");
//        if (rv == kShortcut)
//           return 42;
//        return R__FORWARD_RESULT(FuncThatReturnsRResultOfInt());
//     }
//
// Code using MyIOFunc might look like this:
//
//     auto result = MyIOOperation();
//     if (!result) {
//        /* custom error handling or result.Throw() */
//     }
//     switch (result.Inspect()) {
//        ...
//     }
//
// Note that RResult<void> can be used for a function without return value, like this
//
//     RResult<void> DoSomething()
//     {
//        if (failure)
//           return R__FAIL("user-facing error messge");
//        return RResult<void>::Success();
//     }


namespace ROOT {
namespace Experimental {

// clang-format off
/**
\class ROOT::Experimental::RError
\ingroup Base
\brief Captures diagnostics related to a ROOT runtime error
*/
// clang-format on
class RError {
public:
   struct RLocation {
      RLocation() = default;
      RLocation(const char *func, const char *file, int line)
         : fFunction(func), fSourceFile(file), fSourceLine(line) {}

      // TODO(jblomer) use std::source_location once available
      const char *fFunction;
      const char *fSourceFile;
      int fSourceLine;
   };

private:
   /// User-facing error message
   std::string fMessage;
   /// The location of the error related to fMessage plus upper frames if the error is forwarded through the call stack
   std::vector<RLocation> fStackTrace;

public:
   /// Used by R__FAIL
   RError(const std::string &message, RLocation &&sourceLocation);
   /// Used by R__FORWARD_RESULT
   void AddFrame(RLocation &&sourceLocation);
   /// Add more information to the diagnostics
   void AppendToMessage(const std::string &info) { fMessage += info; }
   /// Format a dignostics report, e.g. for an exception message
   std::string GetReport() const;
   const std::vector<RLocation> &GetStackTrace() const { return fStackTrace; }
};

// clang-format off
/**
\class ROOT::Experimental::RException
\ingroup Base
\brief Base class for all ROOT issued exceptions
*/
// clang-format on
class RException : public std::runtime_error {
   RError fError;
public:
   explicit RException(const RError &error) : std::runtime_error(error.GetReport()), fError(error) {}
   const RError &GetError() const { return fError; }
};


namespace Internal {
// clang-format off
/**
\class ROOT::Experimental::Internal::RResultBase
\ingroup Base
\brief Common handling of the error case for RResult<T> (T != void) and RResult<void>

RResultBase captures a possible runtime error that might have occured.  If the RResultBase leaves the scope unchecked,
it will throw an exception.  RResultBase should only be allocated on the stack, which is helped by deleting the
new operator.  RResultBase is movable but not copyable to avoid throwing multiple exceptions about the same failure.
*/
// clang-format on
class RResultBase {
protected:
   /// This is the nullptr for an RResult representing success
   std::unique_ptr<RError> fError;
   /// Switches to true once the user of an RResult object checks the object status
   bool fIsChecked{false};

   RResultBase() = default;
   explicit RResultBase(RError &&error) : fError(std::make_unique<RError>(std::move(error))) {}

   /// Used by the RResult<T> bool operator
   bool Check() {
      fIsChecked = true;
      return !fError;
   }

public:
   RResultBase(const RResultBase &other) = delete;
   RResultBase(RResultBase &&other) = default;
   RResultBase &operator =(const RResultBase &other) = delete;
   RResultBase &operator =(RResultBase &&other) = default;

   ~RResultBase() noexcept(false);

   RError *GetError() { return fError.get(); }
   /// Throws an RException with fError
   void Throw();

   /// Used by R__FORWARD_ERROR in order to keep track of the stack trace.
   static RError ForwardError(RResultBase &&result, RError::RLocation &&sourceLocation) {
      if (!result.fError) {
         return RError("internal error: attempt to forward error of successful operation",
                       std::move(sourceLocation));
      }
      result.fError->AddFrame(std::move(sourceLocation));
      return *result.fError;
   }

   // Help to prevent heap construction of RResult objects. Unchecked RResult objects in failure state should throw
   // an exception close to the error location. For stack allocated RResult objects, an exception is thrown
   // the latest when leaving the scope. Heap allocated RResult objects in failure state can live much longer making it
   // difficult to trace back the original error.
   void *operator new(std::size_t size) = delete;
   void *operator new(std::size_t, void *) = delete;
   void *operator new[](std::size_t) = delete;
   void *operator new[](std::size_t, void *) = delete;
}; // class RResultBase
} // namespace Internal


// clang-format off
/**
\class ROOT::Experimental::RResult
\ingroup Base
\brief The class is used as a return type for operations that can fail; wraps a value of type T or an RError

RResult enforces checking whether it contains a valid value or an error state. If the RResult leaves the scope
unchecked, it will throw an exception (due to ~RResultBase).
*/
// clang-format on
template <typename T>
class RResult : public Internal::RResultBase {
private:
   /// The result value in case of successful execution
   T fValue;

   // Ensure accessor methods throw in case of errors
   inline void ThrowOnError() {
      if (R__unlikely(fError)) {
         // Accessors can be wrapped in a try-catch block, so throwing the
         // exception here is akin to checking the error.
         //
         // Setting fIsChecked to true also avoids a spurious warning in the RResult destructor
         fIsChecked = true;

         fError->AppendToMessage(" (unchecked RResult access!)");
         throw RException(*fError);
      }
   }

public:
   RResult(const T &value) : fValue(value) {}
   RResult(T &&value) : fValue(std::move(value)) {}
   RResult(RError &&error) : RResultBase(std::move(error)) {}

   RResult(const RResult &other) = delete;
   RResult(RResult &&other) = default;
   RResult &operator =(const RResult &other) = delete;
   RResult &operator =(RResult &&other) = default;

   ~RResult() = default;

   /// Used by R__FORWARD_RESULT in order to keep track of the stack trace in case of errors
   RResult &Forward(RError::RLocation &&sourceLocation) {
      if (fError)
         fError->AddFrame(std::move(sourceLocation));
      return *this;
   }

   /// If the operation was successful, returns a const reference to the inner type.
   /// If there was an error, Inspect() instead throws an exception.
   const T &Inspect() {
      ThrowOnError();
      return fValue;
   }

   /// If the operation was successful, returns the inner type by value.
   ///
   /// For move-only types, Unwrap can only be called once, as it yields ownership of
   /// the inner value to the caller using std::move, potentially leaving the
   /// RResult in an unspecified state.
   ///
   /// If there was an error, Unwrap() instead throws an exception.
   T Unwrap() {
      ThrowOnError();
      return std::move(fValue);
   }

   explicit operator bool() { return Check(); }
};

/// RResult<void> has no data member and no Inspect() method but instead a Success() factory method
template<>
class RResult<void> : public Internal::RResultBase {
private:
   RResult() = default;

public:
   /// Returns a RResult<void> that captures the successful execution of the function
   static RResult Success() { return RResult(); }
   RResult(RError &&error) : RResultBase(std::move(error)) {}

   RResult(const RResult &other) = delete;
   RResult(RResult &&other) = default;
   RResult &operator =(const RResult &other) = delete;
   RResult &operator =(RResult &&other) = default;

   ~RResult() = default;

   /// Used by R__FORWARD_RESULT in order to keep track of the stack trace in case of errors
   RResult &Forward(RError::RLocation &&sourceLocation) {
      if (fError)
         fError->AddFrame(std::move(sourceLocation));
      return *this;
   }

   /// Short-hand method to throw an exception in the case of errors. Does nothing for
   /// successful RResults.
   void ThrowOnError() {
      if (!Check())
         Throw();
   }

   explicit operator bool() { return Check(); }
};

/// Short-hand to return an RResult<T> in an error state; the RError is implicitly converted into RResult<T>
#define R__FAIL(msg) ROOT::Experimental::RError(msg, {R__LOG_PRETTY_FUNCTION, __FILE__, __LINE__})
/// Short-hand to return an RResult<T> value from a subroutine to the calling stack frame
#define R__FORWARD_RESULT(res) std::move(res.Forward({R__LOG_PRETTY_FUNCTION, __FILE__, __LINE__}))
/// Short-hand to return an RResult<T> in an error state (i.e. after checking)
#define R__FORWARD_ERROR(res) res.ForwardError(std::move(res), {R__LOG_PRETTY_FUNCTION, __FILE__, __LINE__})
} // namespace Experimental
} // namespace ROOT

#endif
