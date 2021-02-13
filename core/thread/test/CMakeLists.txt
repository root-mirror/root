# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

ROOT_ADD_UNITTEST_DIR(Core Thread Hist)

ROOT_ADD_GTEST(testTThreadedObject testTThreadedObject.cxx LIBRARIES Hist)

ROOT_ADD_GTEST(testInterpreterLock testInterpreterLock.cxx LIBRARIES Imt REPEATS 100)
