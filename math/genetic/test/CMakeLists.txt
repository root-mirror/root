# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

if(NOT DEFINED ROOT_SOURCE_DIR)
   cmake_minimum_required(VERSION 3.9)
   project(genetic-tests)
   find_package(ROOT REQUIRED)
   include(${ROOT_USE_FILE})
endif()

set(Libraries Core RIO Net Hist Graf Graf3d Gpad Tree
              Rint Postscript Matrix Physics MathCore Thread Genetic)

set(TestGeneticSource
      testGAMinimizer.cxx
      GAMinTutorial.cxx)


#---Build and add all the defined test in the list---------------
foreach(file ${TestGeneticSource})
  get_filename_component(testname ${file} NAME_WE)
  ROOT_EXECUTABLE(${testname} ${file} LIBRARIES ${Libraries})
  ROOT_ADD_TEST(genetic-${testname} COMMAND ${testname})
endforeach()
