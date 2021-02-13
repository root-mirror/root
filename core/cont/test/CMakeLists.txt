# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

ROOT_ADD_UNITTEST_DIR(Core)

ROOT_ADD_GTEST(testTypedIteration testTypedIteration.cxx LIBRARIES Core)
ROOT_ADD_GTEST(TSeqTests TSeqTests.cxx LIBRARIES Core)
