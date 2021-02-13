# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

include_directories(
  ../res
  ../../clingutils/res
  ../../foundation/res
  ${CLING_INCLUDE_DIRS}
  )

ROOT_ADD_GTEST(
  TClingTest

  TClingTests.cxx
  TClingCallFuncTests.cxx
  TClingDataMemberInfoTests.cxx
  TClingMethodInfoTests.cxx

  LIBRARIES
  Core)
