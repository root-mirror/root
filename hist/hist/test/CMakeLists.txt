# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

ROOT_ADD_GTEST(testTProfile2Poly test_tprofile2poly.cxx LIBRARIES Hist Matrix MathCore RIO)
ROOT_ADD_GTEST(testTH2PolyBinError test_TH2Poly_BinError.cxx LIBRARIES Hist Matrix MathCore RIO)
ROOT_ADD_GTEST(testTH2PolyAdd test_TH2Poly_Add.cxx LIBRARIES Hist Matrix MathCore RIO)
ROOT_ADD_GTEST(testTHn THn.cxx LIBRARIES Hist Matrix MathCore RIO)
ROOT_ADD_GTEST(testTH1 test_TH1.cxx LIBRARIES Hist)
ROOT_ADD_GTEST(testTFormula test_TFormula.cxx LIBRARIES Hist)
ROOT_ADD_GTEST(testTKDE test_tkde.cxx LIBRARIES Hist)
ROOT_ADD_GTEST(testTH1FindFirstBinAbove test_TH1_FindFirstBinAbove.cxx LIBRARIES Hist)
ROOT_ADD_GTEST(test_TEfficiency test_TEfficiency.cxx LIBRARIES Hist)
ROOT_ADD_GTEST(TGraphMultiErrorsTests TGraphMultiErrorsTests.cxx LIBRARIES Hist RIO)
ROOT_ADD_GTEST(test_TF123_Moments test_TF123_Moments.cxx LIBRARIES Hist)
ROOT_ADD_GTEST(test_THBinIterator test_THBinIterator.cxx LIBRARIES Hist)

if(fftw3)
  ROOT_ADD_GTEST(testTF1 test_tf1.cxx LIBRARIES Hist)
endif()

if(clad)
  ROOT_ADD_GTEST(TFormulaGradientTests TFormulaGradientTests.cxx LIBRARIES Core MathCore Hist)
endif()
