# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
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
  ROOT/RColumn.hxx
  ROOT/RColumnElement.hxx
  ROOT/RColumnModel.hxx
  ROOT/REntry.hxx
  ROOT/RField.hxx
  ROOT/RFieldValue.hxx
  ROOT/RFieldVisitor.hxx
  ROOT/RNTuple.hxx
  ROOT/RNTupleDescriptor.hxx
  ROOT/RNTupleMetrics.hxx
  ROOT/RNTupleModel.hxx
  ROOT/RNTupleOptions.hxx
  ROOT/RNTupleUtil.hxx
  ROOT/RNTupleView.hxx
  ROOT/RPage.hxx
  ROOT/RPageAllocator.hxx
  ROOT/RPagePool.hxx
  ROOT/RPageStorage.hxx
  ROOT/RPageStorageChain.hxx
  ROOT/RPageStorageFriend.hxx
  ROOT/RPageStorageRaw.hxx
  ROOT/RPageStorageRoot.hxx
SOURCES
  v7/src/RColumn.cxx
  v7/src/RColumnElement.cxx
  v7/src/RField.cxx
  v7/src/RFieldVisitor.cxx
  v7/src/REntry.cxx
  v7/src/RNTuple.cxx
  v7/src/RNTupleDescriptor.cxx
  v7/src/RNTupleDescriptorFmt.cxx
  v7/src/RNTupleMetrics.cxx
  v7/src/RNTupleModel.cxx
  v7/src/RPage.cxx
  v7/src/RPageAllocator.cxx
  v7/src/RPagePool.cxx
  v7/src/RPageStorage.cxx
  v7/src/RPageStorageChain.cxx
  v7/src/RPageStorageFriend.cxx
  v7/src/RPageStorageRaw.cxx
  v7/src/RPageStorageRoot.cxx
LINKDEF
  LinkDef.h
DEPENDENCIES
  RIO
  ROOTVecOps
)

ROOT_ADD_TEST_SUBDIRECTORY(v7/test)
