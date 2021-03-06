# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

set(Libraries Physics EG)

ROOT_ADD_GTEST(testTDecay tTDecay.cxx LIBRARIES ${Libraries})
ROOT_ADD_GTEST(testTGenDecay tTGenDecay.cxx LIBRARIES ${Libraries})
ROOT_ADD_GTEST(testTGenFoamDecay tTGenFoamDecay.cxx LIBRARIES ${Libraries})
