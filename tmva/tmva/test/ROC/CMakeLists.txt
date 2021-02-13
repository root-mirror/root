# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building TMVA/ROC tests.
# @author Kim Albertsson
############################################################################

ROOT_ADD_GTEST(testRocWeights
                testRocWeights.cxx
                LIBRARIES TMVA)
