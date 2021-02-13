# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT TMVA tests.
# @author Kim Albertsson
############################################################################

# Tests using google test
ROOT_ADD_GTEST(TestRandomGenerator
               TestRandomGenerator.cxx
               LIBRARIES TMVA)
ROOT_ADD_GTEST(TestOptimizeConfigParameters
               TestOptimizeConfigParameters.cxx
               LIBRARIES TMVA)

if(dataframe)
    # RTensor
    ROOT_ADD_GTEST(rtensor rtensor.cxx LIBRARIES ROOTVecOps TMVA)
    ROOT_ADD_GTEST(rtensor-iterator rtensor_iterator.cxx LIBRARIES ROOTVecOps TMVA)
    ROOT_ADD_GTEST(rtensor-utils rtensor_utils.cxx LIBRARIES ROOTVecOps TMVA ROOTDataFrame)
    # RStandardScaler
    ROOT_ADD_GTEST(rstandardscaler rstandardscaler.cxx LIBRARIES ROOTVecOps TMVA ROOTDataFrame)
    # RReader
    ROOT_ADD_GTEST(rreader rreader.cxx LIBRARIES ROOTVecOps TMVA ROOTDataFrame)
    # Tree inference system and user interface
    ROOT_ADD_GTEST(branchlessForest branchlessForest.cxx LIBRARIES TMVA)
    ROOT_ADD_GTEST(rbdt rbdt.cxx LIBRARIES ROOTVecOps TMVA)
endif()

if(dataframe AND NOT pyroot_legacy)
  find_python_module(xgboost QUIET)
  if (PY_XGBOOST_FOUND)
    ROOT_ADD_PYUNITTEST(rbdt_xgboost rbdt_xgboost.py)
  endif()
endif()
