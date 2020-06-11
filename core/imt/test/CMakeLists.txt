# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

ROOT_ADD_UNITTEST_DIR(Imt Thread ${TBB_LIBRARIES})

ROOT_ADD_GTEST(testImt testTFuture.cxx testTTaskGroup.cxx LIBRARIES Imt)
ROOT_ADD_GTEST(testTaskArena testRTaskArena.cxx LIBRARIES Imt ${TBB_LIBRARIES} FAILREGEX "")
