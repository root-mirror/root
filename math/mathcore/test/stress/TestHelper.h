#ifndef ROOT_TESTHELPER
#define ROOT_TESTHELPER

#include "gtest/gtest.h"

#include <string>

bool OutsideBounds(double v1, double v2, double scale);

// Compared to ASSERT_NEAR, this function takes into account also the relative error
::testing::AssertionResult IsNear(std::string name, double v1, double v2, double scale = 2.0);

#endif // ROOT_TESTHELPER
