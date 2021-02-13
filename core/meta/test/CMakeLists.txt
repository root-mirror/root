# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

ROOT_ADD_GTEST(testStatusBitsChecker testStatusBitsChecker.cxx LIBRARIES Core)
ROOT_ADD_GTEST(testHashRecursiveRemove testHashRecursiveRemove.cxx LIBRARIES Core)
if(NOT MSVC OR win_broken_tests)
  ROOT_ADD_GTEST(testTClass testTClass.cxx LIBRARIES Core)
endif()
ROOT_ADD_GTEST(testTEnum testTEnum.cxx LIBRARIES Core)
configure_file(stlDictCheck.h . COPYONLY)
configure_file(stlDictCheckAux.h . COPYONLY)
