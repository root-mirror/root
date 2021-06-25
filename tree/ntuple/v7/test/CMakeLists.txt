# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

# @author Jakob Blomer CERN

if(NOT dataframe)
  return()
endif()

ROOT_STANDARD_LIBRARY_PACKAGE(CustomStruct
                              NO_INSTALL_HEADERS
                              HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/CustomStruct.hxx
                              SOURCES CustomStruct.cxx
                              LINKDEF CustomStructLinkDef.h
                              DEPENDENCIES RIO)
if(MSVC)
  add_custom_command(TARGET CustomStruct POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/libCustomStruct.dll
                                     ${CMAKE_CURRENT_BINARY_DIR}/libCustomStruct.dll)
endif()

ROOT_ADD_GTEST(ntuple_basics ntuple_basics.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)
ROOT_ADD_GTEST(ntuple_cluster ntuple_cluster.cxx LIBRARIES ROOTNTuple)
ROOT_ADD_GTEST(ntuple_descriptor ntuple_descriptor.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)
ROOT_ADD_GTEST(ntuple_friends ntuple_friends.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)
ROOT_ADD_GTEST(ntuple_merger ntuple_merger.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)
ROOT_ADD_GTEST(ntuple_metrics ntuple_metrics.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)
ROOT_ADD_GTEST(ntuple_packing ntuple_packing.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)
ROOT_ADD_GTEST(ntuple_pages ntuple_pages.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)
ROOT_ADD_GTEST(ntuple_print ntuple_print.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)
ROOT_ADD_GTEST(ntuple_rdf ntuple_rdf.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)
ROOT_ADD_GTEST(ntuple_types ntuple_types.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)
ROOT_ADD_GTEST(ntuple_view ntuple_view.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)
ROOT_ADD_GTEST(ntuple_zip ntuple_zip.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)

ROOT_ADD_GTEST(rfield_class rfield_class.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)
ROOT_ADD_GTEST(rfield_string rfield_string.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)
ROOT_ADD_GTEST(rfield_variant rfield_variant.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)
ROOT_ADD_GTEST(rfield_vector rfield_vector.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)

ROOT_ADD_GTEST(ntuple_minifile ntuple_minifile.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)
ROOT_ADD_GTEST(ntuple_show ntuple_show.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)
ROOT_ADD_GTEST(ntuple_storage ntuple_storage.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)
ROOT_ADD_GTEST(ntuple_extended ntuple_extended.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)

if(daos OR daos_mock)
  # UUID of the DAOS pool used for testing, if not provided (may be any for libdaos_mock).
  if(NOT daos_test_pool)
    set(daos_test_pool 1b245c52-765e-449d-80b1-8dce93bafc4b)
  endif()
  set_property(SOURCE ntuple_storage_daos.cxx
               APPEND PROPERTY COMPILE_DEFINITIONS R__DAOS_TEST_POOL="${daos_test_pool}")

  ROOT_ADD_GTEST(ntuple_storage_daos ntuple_storage_daos.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)
endif()

if (davix OR builtin_davix)
  ROOT_ADD_GTEST(ntuple_storage_s3 ntuple_storage_s3.cxx LIBRARIES ROOTDataFrame ROOTNTuple MathCore CustomStruct)
  target_link_libraries(ntuple_storage_s3 PRIVATE ${DAOS_LIBRARIES})
endif()
