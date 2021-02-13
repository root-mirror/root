# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

project(mathmore-tests)

set(Libraries Core RIO Net Hist Graf Graf3d Gpad Tree
              Rint Postscript Matrix Physics Thread MathCore MathMore)

set(TestMathMoreSource
     testChebyshev.cxx
     testPermute.cxx
     testDerivation.cxx
     testGSLIntegration.cxx
     testMCIntegration.cxx
     testGSLRootFinder.cxx
     testMultiRootFinder.cxx
     testMinimization1D.cxx
     testInterpolation.cxx
     testRandom.cxx
     testRandomDist.cxx
     testQuasiRandom.cxx
     testSpecFunc.cxx
     testStatFunc.cxx
     testFunctor.cxx
     testVavilov.cxx
     simanTSP.cxx)

set(testFunctor_LABELS longtest)
set(testPermute_LABELS longtest)
set(testRandom_LABELS longtest)
set(testMCIntegration_LABELS longtest)

if(unuran)
  list(APPEND Libraries Unuran)
  add_definitions(-DHAVE_UNURAN)
endif()

if(minuit2)
  list(APPEND TestSource fit/testMinim.cxx)
  list(APPEND Libraries Minuit2)
endif()

#some tests requires directly gsl
include_directories(${GSL_INCLUDE_DIR})
add_definitions(-DHAVE_ROOTLIBS)

#---Build and add all the defined test in the list---------------
foreach(file ${TestMathMoreSource})
  get_filename_component(testname ${file} NAME_WE)
  ROOT_EXECUTABLE(${testname} ${file} LIBRARIES ${GSL_LIBRARIES} ${Libraries})
  ROOT_ADD_TEST(mathmore-${testname} COMMAND ${testname} LABELS ${${testname}_LABELS})
endforeach()

ROOT_ADD_GTEST(stressMathMoreUnit testStress.cxx StatFunction.cxx LIBRARIES Core MathCore MathMore)
