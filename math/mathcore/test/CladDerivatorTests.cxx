/// \file CladDerivatorTests.cxx
///
/// \brief The file contain unit tests which test the CladDerivator facility.
///
/// \author Vassil Vassilev <vvasilev@cern.ch>
///
/// \date July, 2018
///
/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <TInterpreter.h>
#include <TInterpreterValue.h>

#include "gtest/gtest.h"

TEST(CladDerivator, Sanity)
{
   gInterpreter->Declare(R"cpp(
                           #include <Math/CladDerivator.h>
                           static double pow2(double x) { return x * x; })cpp");
   auto value = gInterpreter->CreateTemporary();
   std::string code = "clad::differentiate(pow2, 0).execute(1)";
   ASSERT_TRUE(gInterpreter->Evaluate(code.c_str(), *value));
   ASSERT_FLOAT_EQ(2, value->GetAsDouble());
   delete value;
}
