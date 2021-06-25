# Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT ntuple package
# @author Jakob Blomer CERN
############################################################################

if(NOT root7)
  return()
endif()

if (davix OR builtin_davix)
  list(APPEND ROOTNTuple_EXTRA_DEPS RDAVIX)
endif()

ROOT_STANDARD_LIBRARY_PACKAGE(ROOTNTuple
HEADERS
  ROOT/RCluster.hxx
  ROOT/RClusterPool.hxx
  ROOT/RColumn.hxx
  ROOT/RColumnElement.hxx
  ROOT/RColumnModel.hxx
  ROOT/REntry.hxx
  ROOT/RField.hxx
  ROOT/RFieldValue.hxx
  ROOT/RFieldVisitor.hxx
  ROOT/RMiniFile.hxx
  ROOT/RNTuple.hxx
  ROOT/RNTupleDescriptor.hxx
  ROOT/RNTupleMerger.hxx
  ROOT/RNTupleMetrics.hxx
  ROOT/RNTupleModel.hxx
  ROOT/RNTupleOptions.hxx
  ROOT/RNTupleUtil.hxx
  ROOT/RNTupleView.hxx
  ROOT/RNTupleZip.hxx
  ROOT/RPage.hxx
  ROOT/RPageAllocator.hxx
  ROOT/RPagePool.hxx
  ROOT/RPageSinkBuf.hxx
  ROOT/RPageSourceFriends.hxx
  ROOT/RPageStorage.hxx
  ROOT/RPageStorageFile.hxx
SOURCES
  v7/src/RCluster.cxx
  v7/src/RClusterPool.cxx
  v7/src/RColumn.cxx
  v7/src/RColumnElement.cxx
  v7/src/RField.cxx
  v7/src/RFieldVisitor.cxx
  v7/src/REntry.cxx
  v7/src/RMiniFile.cxx
  v7/src/RNTuple.cxx
  v7/src/RNTupleDescriptor.cxx
  v7/src/RNTupleDescriptorFmt.cxx
  v7/src/RNTupleMerger.cxx
  v7/src/RNTupleMetrics.cxx
  v7/src/RNTupleModel.cxx
  v7/src/RNTupleUtil.cxx
  v7/src/RPage.cxx
  v7/src/RPageAllocator.cxx
  v7/src/RPagePool.cxx
  v7/src/RPageSinkBuf.cxx
  v7/src/RPageSourceFriends.cxx
  v7/src/RPageStorage.cxx
  v7/src/RPageStorageFile.cxx
LINKDEF
  LinkDef.h
DEPENDENCIES
  Imt
  RIO
  ROOTVecOps
  ${ROOTNTuple_EXTRA_DEPS}
)

# Enable RNTuple support for Intel DAOS
if(daos OR daos_mock)
  set(ROOTNTuple_EXTRA_HEADERS ROOT/RPageStorageDaos.hxx)
  target_sources(ROOTNTuple PRIVATE v7/src/RDaos.cxx v7/src/RPageStorageDaos.cxx)
  target_compile_definitions(ROOTNTuple PRIVATE R__ENABLE_DAOS)

  target_include_directories(ROOTNTuple PRIVATE ${UUID_INCLUDE_DIR})
  target_link_libraries(ROOTNTuple PRIVATE ${UUID_LIBRARY})
  if(daos_mock)
    target_include_directories(ROOTNTuple PRIVATE v7/inc/ROOT/libdaos_mock/)
    target_sources(ROOTNTuple PRIVATE v7/src/libdaos_mock/libdaos_mock.cxx)
  else()
    target_include_directories(ROOTNTuple PRIVATE ${DAOS_INCLUDE_DIRS})
    target_link_libraries(ROOTNTuple PRIVATE ${DAOS_LIBRARIES})
  endif()
endif()

# Enable RNTuple support for S3 backend through Davix
if(davix OR builtin_davix)
  set(ROOTNTuple_EXTRA_HEADERS ${ROOTNTuple_EXTRA_HEADERS} ROOT/RPageStorageS3.hxx)
  target_sources(ROOTNTuple PRIVATE v7/src/RPageStorageS3.cxx)
  target_compile_definitions(ROOTNTuple PRIVATE R__ENABLE_DAVIX)

  target_include_directories(ROOTNTuple PRIVATE ${DAVIX_INCLUDE_DIR})
  target_link_libraries(ROOTNTuple PRIVATE ${DAVIX_LIBRARY})
endif()

if(MSVC)
  target_compile_definitions(ROOTNTuple PRIVATE _USE_MATH_DEFINES)
endif()

ROOT_ADD_TEST_SUBDIRECTORY(v7/test)
