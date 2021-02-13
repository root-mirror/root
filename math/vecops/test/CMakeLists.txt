# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

ROOT_ADD_GTEST(vecops_rvec vecops_rvec.cxx LIBRARIES Physics ROOTVecOps GenVector RIO Tree)
ROOT_ADD_GTEST(vecops_radoptallocator vecops_radoptallocator.cxx LIBRARIES Core ROOTVecOps)
