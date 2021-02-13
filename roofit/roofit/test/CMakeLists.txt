# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

# @author Stephan Hageboeck, CERN, 2019

ROOT_ADD_GTEST(testRooGaussian testRooGaussian.cxx LIBRARIES RooFitCore RooFit)
ROOT_ADD_GTEST(testRooPoisson testRooPoisson.cxx LIBRARIES RooFitCore RooFit)
ROOT_ADD_GTEST(testRooBernstein testRooBernstein.cxx LIBRARIES RooFitCore RooFit)
ROOT_ADD_GTEST(testRooJohnson testRooJohnson.cxx LIBRARIES Gpad RooFitCore RooFit)
ROOT_ADD_GTEST(testRooExponential testRooExponential.cxx
  COPY_TO_BUILDDIR ${CMAKE_CURRENT_SOURCE_DIR}/exponentialPdf.root
  LIBRARIES Core RooFitCore RooFit)
