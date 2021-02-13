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
  v7/src/RPageStorage.cxx
  v7/src/RPageStorageFile.cxx
LINKDEF
  LinkDef.h
DEPENDENCIES
  Imt
  RIO
  ROOTVecOps
)

if(MSVC)
  target_compile_definitions(ROOTNTuple PRIVATE _USE_MATH_DEFINES)
endif()

ROOT_ADD_TEST_SUBDIRECTORY(v7/test)
