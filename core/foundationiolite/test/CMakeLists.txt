# Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

ROOT_ADD_GTEST(RRawFile RRawFile.cxx LIBRARIES RIO)
if(uring AND NOT DEFINED ENV{ROOTTEST_IGNORE_URING})
  ROOT_ADD_GTEST(RIoUring RIoUring.cxx LIBRARIES RIO)
endif()
