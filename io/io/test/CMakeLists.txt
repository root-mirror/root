# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

ROOT_ADD_GTEST(RRawFile RRawFile.cxx LIBRARIES RIO)
ROOT_ADD_GTEST(TFile TFileTests.cxx LIBRARIES RIO)
ROOT_ADD_GTEST(TBufferMerger TBufferMerger.cxx LIBRARIES RIO Imt Tree)
ROOT_ADD_GTEST(TFileMerger TFileMergerTests.cxx LIBRARIES RIO Tree)
ROOT_ADD_GTEST(TROMemFile TROMemFileTests.cxx LIBRARIES RIO Tree)
if(uring AND NOT DEFINED ENV{ROOTTEST_IGNORE_URING})
  ROOT_ADD_GTEST(RIoUring RIoUring.cxx LIBRARIES RIO)
endif()
