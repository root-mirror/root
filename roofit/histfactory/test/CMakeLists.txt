# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

# @author Stephan Hageboeck CERN, 2019

ROOT_ADD_GTEST(testHistFactory testHistFactory.cxx
  LIBRARIES RooFitCore RooFit RooStats HistFactory
  COPY_TO_BUILDDIR ${CMAKE_CURRENT_SOURCE_DIR}/ref_6.16_example_UsingC_channel1_meas_model.root ${CMAKE_CURRENT_SOURCE_DIR}/ref_6.16_example_UsingC_combined_meas_model.root)
