/// \file RError.cxx
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

#include <ROOT/RConfig.hxx> // for R__[un]likely
#include <ROOT/RError.hxx>
#include <ROOT/RLogger.hxx> // for R__LOG_WARNING

#include <exception>
#include <string>
#include <utility>


std::string ROOT::Experimental::RError::GetReport() const
{
   auto report = fMessage + "\nAt:\n";
   for (const auto &loc : fStackTrace) {
      report += "  " + std::string(loc.fFunction) + " [" + std::string(loc.fSourceFile) + ":" +
         std::to_string(loc.fSourceLine) + "]\n";
   }
   return report;
}

ROOT::Experimental::RError::RError(
   const std::string &message, RLocation &&sourceLocation)
   : fMessage(message)

{
   // Avoid frequent reallocations as we move up the call stack
   fStackTrace.reserve(32);
   AddFrame(std::move(sourceLocation));
}

void ROOT::Experimental::RError::AddFrame(RLocation &&sourceLocation)
{
   fStackTrace.emplace_back(sourceLocation);
}


ROOT::Experimental::Internal::RResultBase::~RResultBase() noexcept(false)
{
   if (R__unlikely(fError && !fIsChecked)) {
      // Prevent from throwing if the object is deconstructed in the course of stack unwinding for another exception
#if __cplusplus >= 201703L
      if (std::uncaught_exceptions() == 0)
#else
      if (!std::uncaught_exception())
#endif
      {
         throw RException(*fError);
      } else {
         R__LOG_WARNING() << "unhandled RResult exception during stack unwinding";
      }
   }
}


void ROOT::Experimental::Internal::RResultBase::Throw()
{
   throw ROOT::Experimental::RException(*fError);
}
