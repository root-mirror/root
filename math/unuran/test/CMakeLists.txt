# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

set(Libraries Core MathCore Minuit Unuran Gpad Hist)

if(mathmore)
  list(APPEND Libraries MathMore)
endif()

set(TestUnuranSource
       unuranSimple.cxx
       unuranHist.cxx
       unuranDiscrete.cxx
       unuranDistr.cxx
       unuranMulti2D.cxx
       unuranMultiDim.cxx
)

# Configure some tests as longtest
set(unuranMultiDim_LABELS longtest)
set(unuranDiscrete_LABELS longtest)

#---Build and add all the defined test in the list---------------
foreach(file ${TestUnuranSource})
  get_filename_component(testname ${file} NAME_WE)
  ROOT_EXECUTABLE(${testname} ${file} LIBRARIES ${Libraries})
  ROOT_ADD_TEST(unuran-${testname} COMMAND ${testname} LABELS ${${testname}_LABELS})
endforeach()

# google tests
ROOT_ADD_GTEST(testUnuranSampler unuranSampler.cxx LIBRARIES ${Libraries} )
