# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building TMVA/Method tests.
# @author Kim Albertsson
############################################################################

# As the reference implementation is deprecated and
# can't be instantiated, we can only run this test
# if there if we can run CPU or CUDA versions.
#if (defined DNNCPU || defined DNNCUDA)
if(tmva-cpu OR tmva-gpu)
    ROOT_ADD_GTEST(TestMethodDNNValidationSize
                TestMethodDNNValidationSize.cxx
                LIBRARIES TMVA)
endif()
