/// \file ROOT/RError.h
/// \ingroup Base ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2019-12-11
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RError
#define ROOT7_RError

#include <ROOT/RConfig.hxx>

#include <cstddef>
#include <new>
#include <stdexcept>
#include <string>

namespace ROOT {
namespace Experimental {

// clang-format off
/**
\class ROOT::Experimental::RException
\ingroup Base
\brief Base class for all ROOT issued exceptions
*/
// clang-format on
class RException : public std::runtime_error {
public:
   explicit RException(const std::string &what) : std::runtime_error(what) {}
};


namespace Detail {

// clang-format off
/**
\class ROOT::Experimental::Detail::RStatusType
\ingroup Base
\brief Values of this type have range of values dedicated to indicate errors, such as negative ints for system calls.

Derived classes need to implement IsError() to select the range of error values.
*/
// clang-format on
template <typename T, class DerivedT>
class RStatusType {
protected:
   T fValue;

public:
   using ValueType_t = T;

   explicit RStatusType(const T &value) : fValue(value) {};
   const T &Get() const { return fValue; }

   bool IsError() const { return static_cast<DerivedT *>(this)->IsError(); }
};

// clang-format off
/**
\class ROOT::Experimental::Detail::RStatusTypeBool
\ingroup Base
\brief For routines that indicate success by returning true
*/
// clang-format on
class RStatusTypeBool : public RStatusType<bool, RStatusTypeBool> {
public:
   explicit RStatusTypeBool(bool value) : RStatusType<bool, RStatusTypeBool>(value) {}
   bool IsError() const { return fValue == false; }
};

// clang-format off
/**
\class ROOT::Experimental::Detail::RStatusTypeSyscall
\ingroup Base
\brief For system calls that return 0 (or a meaningful integer) on success
*/
// clang-format on
class RStatusTypeSyscall : public RStatusType<int, RStatusTypeSyscall> {
public:
   explicit RStatusTypeSyscall(int value) : RStatusType<int, RStatusTypeSyscall>(value) {}
   bool IsError() const { return fValue < 0; }
};

// clang-format off
/**
\class ROOT::Experimental::Detail::RStatusBase
\ingroup Base
\brief Type-independent logic of the RStatus wrapper
*/
// clang-format on
class RStatusBase {
protected:
   static thread_local bool fgThrowInstantExceptions; /* true by default */
public:
   static void SetThrowInstantExceptions(bool value) { fgThrowInstantExceptions = value; }
};

// clang-format off
/**
\class ROOT::Experimental::Detail::RStatus
\ingroup Base
\brief Wraps a return value such that unchecked error states trigger an exception

The Status wrapper itself is movable but not copyable to prevent multiple exceptions from being thrown.
The wrapper is also only stack allocatable because throwing an exception depends on the object going out of scope.

If the class should be used like a C return type, it has to be cleared by a call to IsError() or IsValid() or
ClearError().  Otherwise, an error state will throw an exception nevertheless when the object goes out of scope.
*/
// clang-format on
template <typename T> // Has to be RStatusType
class RStatus : public RStatusBase {
   T fStatus;
   bool fIsChecked = false;

public:
   // Named constructor for error cases
   static RStatus Fail(typename T::ValueType_t value, const std::string &why)
   {
      // Fast path reserved for return code checking
      if (R__unlikely(fgThrowInstantExceptions))
         throw RException(why);
      return RStatus(value);
   }

   RStatus(const RStatus &other) = delete;
   RStatus(RStatus &&other) = default;
   RStatus &operator =(const RStatus &other) = delete;
   RStatus &operator =(RStatus &&other) = default;

   // Construct from return value
   explicit RStatus(const typename T::ValueType_t &value) : fStatus(value) {}
   // Assign a return value
   RStatus &operator =(const typename T::ValueType_t &value) { this->fStatus = value; }

   ~RStatus() noexcept(false)
   {
      if (!fIsChecked && fStatus.IsError()) {
         // Prevent from throwing if the object is deconstructed in the course of stack unwinding for another exception
#if __cplusplus >= 201703L
         if (std::uncaught_exceptions() == 0)
#else
         if (!std::uncaught_exception())
#endif
         {
            throw RException("unchecked error");
         }
      }
   }

   bool IsError()
   {
      fIsChecked = true;
      return fStatus.IsError();
   }
   bool IsValid() { return !IsError(); }
   void ClearError() { fIsChecked = true; }

   operator typename T::ValueType_t() const { return fStatus.Get(); }

   // Prevent heap construction of RStatus objects
   void *operator new(std::size_t size) = delete;
   void *operator new(std::size_t, void *) = delete;
   void *operator new[](std::size_t) = delete;
   void *operator new[](std::size_t, void *) = delete;
};

} // namespace Detail

void SetThrowInstantExceptions(bool value);

using RStatusBool = Detail::RStatus<Detail::RStatusTypeBool>;
using RStatusSyscall = Detail::RStatus<Detail::RStatusTypeSyscall>;

} // namespace Experimental
} // namespace ROOT

#endif
