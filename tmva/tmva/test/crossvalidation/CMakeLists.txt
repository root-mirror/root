# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT TMVA tests.
# @author Kim Albertsson
############################################################################

set(Libraries TMVA)

# Tests utilising google test
ROOT_ADD_GTEST(testCrossValidationSplitting
               TestCrossValidationSplitting.cxx
               LIBRARIES ${Libraries})
if(NOT MSVC)
  #---Multiproc is not supported on Windows
  ROOT_ADD_GTEST(testCrossValidationMultiProc
                 TestCrossValidationMultiProc.cxx
                 LIBRARIES ${Libraries})
endif()

# Tests
ROOT_EXECUTABLE(testCrossValidationSerialise
                TestCrossValidationSerialise.cxx
                LIBRARIES ${Libraries})
ROOT_ADD_TEST(TMVA-CrossValidation-Serialise
              COMMAND testCrossValidationSerialise)
