/// \file ROOTUnitTestSupport.h
///
/// \brief The file contains facilities allowing easier writing of in-tree unit tests
///
/// \author Pratyush Das <reikdas@gmail.com>
/// \author Vassil Vassilev <vvasilev@cern.ch>
///
/// \date April, 2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_UNITTESTSUPPORT_H
#define ROOT_UNITTESTSUPPORT_H

#include "gtest/gtest.h"
#include "gmock/gmock.h"

namespace {
using testing::StartsWith;
using testing::StrEq;
using testing::internal::GetCapturedStderr;
using testing::internal::CaptureStderr;
using testing::internal::RE;
class ExpectedDiagRAII {
public:
   enum ExpectedDiagKind {
      EDK_NoDiag = 0,
      EDK_Info,
      EDK_Warning,
      EDK_Error,
      EDK_SysError
   };
private:
   ExpectedDiagKind fDiagKind;
   std::string fExpectedRoutine;
   std::string fExpectedDiag;
   void pop()
   {
      // Diagnostics in ROOT have the format:
      // Error|Warning|Info|...| in <Routine>: free text
      std::string Seen = GetCapturedStderr();

      // Try to reconstruct the precise expected string.
      std::string Expected;
      switch(fDiagKind) {
      default:
         assert (0 && "Unsupported diag kind.");
         break;
      case EDK_NoDiag:
         EXPECT_THAT(Seen, StrEq(""));
         return;
      case EDK_Error:
         Expected = "Error";
         break;
      case EDK_Warning:
         Expected = "Warning";
         break;
      case EDK_Info:
         Expected = "Info";
         break;
      case EDK_SysError:
         Expected = "SysError";
         break;
      }

      // Check if the Diag kind matches what we saw.
      EXPECT_THAT(Seen, StartsWith(Expected));

      Expected += " in ";
      Expected += "<" + fExpectedRoutine + ">: ";

      // Check if the routine matches what we saw.
      EXPECT_THAT(Seen, StartsWith(Expected));

      Expected += fExpectedDiag;

      // The captured stderr also includes new lines.
      Expected += "\n";

      EXPECT_THAT(Seen, StrEq(Expected));
   }

public:
   ExpectedDiagRAII(ExpectedDiagKind DiagKind): fDiagKind(DiagKind) {
      assert(DiagKind == ExpectedDiagRAII::EDK_NoDiag);
      CaptureStderr();
   }

   ExpectedDiagRAII(ExpectedDiagKind DiagKind, std::string InRoutine,
                    std::string E)
      : fDiagKind(DiagKind), fExpectedRoutine(InRoutine), fExpectedDiag(E) {
      CaptureStderr();
   }
   ~ExpectedDiagRAII() { pop(); }
};
}

#define ROOT_EXPECT_ERROR(expression, where, expected_diag )            \
   {                                                                    \
      ExpectedDiagRAII EE(ExpectedDiagRAII::EDK_Error, where,           \
                          expected_diag);                               \
      expression;                                                       \
   }

#define ROOT_EXPECT_WARNING(expression, where, expected_diag)           \
   {                                                                    \
      ExpectedDiagRAII EE(ExpectedDiagRAII::EDK_Warning, where,         \
                          expected_diag);                               \
      expression;                                                       \
   }

#define ROOT_EXPECT_INFO(expression, where, expected_diag)              \
   {                                                                    \
      ExpectedDiagRAII EE(ExpectedDiagRAII::EDK_Info, where,            \
                          expected_diag);                               \
      expression;                                                       \
   }

#define ROOT_EXPECT_NODIAG(expression)                                  \
   {                                                                    \
      ExpectedDiagRAII EE(ExpectedDiagRAII::EDK_NoDiag);                \
      expression;                                                       \
   }

#define ROOT_EXPECT_SYSERROR(expression, where, expected_diag)          \
   {                                                                    \
      ExpectedDiagRAII EE(ExpectedDiagRAII::EDK_SysError, where,        \
                          expected_diag);                               \
      expression;                                                       \
   }

#endif // ROOT_UNITTESTSUPPORT_H
