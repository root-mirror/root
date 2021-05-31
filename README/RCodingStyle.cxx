/// \file RCodingStyle.cxx
/// \author Axel Naumann <axel@cern.ch>
/// \date 2018-07-24
// The above entries are mostly for giving some context.
// The "author" field gives a hint whom to contact in case of questions, also
// from within the team. The date shows whether this is ancient or only part
// of the latest release. It's the date of creation / last massive rewrite.

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

// This file demonstrates the coding style to be used for new ROOT files.
// For much of the rationale see Scott Meyer's "Efficient Modern C++" and the
// [C++ core guidelines](https://github.com/isocpp/CppCoreGuidelines/blob/master/CppCoreGuidelines.md),
// plus a fairly large portion of "coding style personality" inherited from 20 years
// of ROOT: there is nothing wrong with them; they make ROOT classes easily
// recognizable. Our users are used to it. No need to change.

// First, include the header belonging to the source file
#include "RCodingStyle.hxx"

// Now the same order as for the header:
// ROOT headers in alphabetical order; all new ROOT headers are under `ROOT/`
#include "ROOT/RAxis.hxx"

// Old ones are found in include/ without prefix.
// New headers should reduce the use of "old" ROOT headers, and instead use
// stdlib headers and new-style headers wherever possible.
#include "TROOT.h" // for TROOT::GetVersion()

// Then standard library headers, in alphabetical order.
#include <vector>

// Include non-ROOT, non-stdlib headers last.
// Rationale: non-ROOT, non-stdlib headers often #define like mad. Reduce interference
// with ROOT or stdlib headers.
#ifdef R__LINUX
// Indent nested preprocessor defines, keeping '#' on column 1:
#include <dlfcn.h>
#endif

// Do not use `namespace ROOT {` in sources; do not rely on `using namespace ROOT` for
// finding the implemented overload. This causes ambiguities if the header and the source
// disagree on the parameters of free functions, e.g.
// ```
// namespace ROOT {
//   double Add(const RExampleClass& a, const RExampleClass& b, int i) { return 17.; }
// }
// ```
// will just compile (despite the wrong `int` argument), while
// ```
// double ROOT::Add(const RExampleClass& a, const RExampleClass& b, int i) { return 17.; }
// ```
// will not.

// Unused parameters are commented out; they often carry useful semantical info.
// We use the same parameter names in headers and sources.
double ROOT::Add(const RExampleClass & /*a*/, const RExampleClass & /*b*/)
{
   return 17.;
}

// Pins the vtable.
ROOT::RExampleClass::~RExampleClass()
{
}

/// Instead of static data members (whether private or public), use outlined static
/// functions. This solves the static initializion fiasco, and delays the initialization
/// to first use.
const std::string &ROOT::RExampleClass::AccessStaticVar()
{
   // NOTE that this initialization is thread-safe!
   static std::string sString = gROOT->GetVersion();
   return sString;
}

// End the file with a trailing newline.
