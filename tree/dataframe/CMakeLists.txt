# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT tree/dataframe package
# @author Danilo Piparo CERN, Pere Mato CERN
############################################################################

if(NOT dataframe)
    return()
endif()

if(arrow)
  list(APPEND RDATAFRAME_EXTRA_HEADERS ROOT/RArrowDS.hxx)
  list(APPEND RDATAFRAME_EXTRA_INCLUDES -I${ARROW_INCLUDE_DIR})
endif()

if(sqlite)
  list(APPEND RDATAFRAME_EXTRA_HEADERS ROOT/RSqliteDS.hxx)
endif()

if(root7)
  list(APPEND RDATAFRAME_EXTRA_HEADERS ROOT/RNTupleDS.hxx)
  list(APPEND RDATAFRAME_EXTRA_DEPS ROOTNTuple)
endif()

if (imt)
  list(APPEND RDATAFRAME_EXTRA_DEPS Imt)
endif(imt)

ROOT_STANDARD_LIBRARY_PACKAGE(ROOTDataFrame
  HEADERS
    ROOT/RCsvDS.hxx
    ROOT/RDataFrame.hxx
    ROOT/RDataSource.hxx
    ROOT/RDFHelpers.hxx
    ROOT/RLazyDS.hxx
    ROOT/RResultPtr.hxx
    ROOT/RResultHandle.hxx
    ROOT/RRootDS.hxx
    ROOT/RSnapshotOptions.hxx
    ROOT/RTrivialDS.hxx
    ROOT/RDF/ActionHelpers.hxx
    ROOT/RDF/ColumnReaderUtils.hxx
    ROOT/RDF/GraphNode.hxx
    ROOT/RDF/GraphUtils.hxx
    ROOT/RDF/HistoModels.hxx
    ROOT/RDF/InterfaceUtils.hxx
    ROOT/RDF/RActionBase.hxx
    ROOT/RDF/RAction.hxx
    ROOT/RDF/RBookedDefines.hxx
    ROOT/RDF/RDefineBase.hxx
    ROOT/RDF/RDefine.hxx
    ROOT/RDF/RDefineReader.hxx
    ROOT/RDF/RDSColumnReader.hxx
    ROOT/RDF/RColumnReaderBase.hxx
    ROOT/RDF/RCutFlowReport.hxx
    ROOT/RDF/RDisplay.hxx
    ROOT/RDF/RFilterBase.hxx
    ROOT/RDF/RFilter.hxx
    ROOT/RDF/RInterface.hxx
    ROOT/RDF/RJittedAction.hxx
    ROOT/RDF/RJittedDefine.hxx
    ROOT/RDF/RJittedFilter.hxx
    ROOT/RDF/RLazyDSImpl.hxx
    ROOT/RDF/RLoopManager.hxx
    ROOT/RDF/RMergeableValue.hxx
    ROOT/RDF/RNodeBase.hxx
    ROOT/RDF/RRangeBase.hxx
    ROOT/RDF/RRange.hxx
    ROOT/RDF/RSlotStack.hxx
    ROOT/RDF/RTreeColumnReader.hxx
    ROOT/RDF/Utils.hxx
    ROOT/RDF/PyROOTHelpers.hxx
    ${RDATAFRAME_EXTRA_HEADERS}
  SOURCES
    src/RActionBase.cxx
    src/RCsvDS.cxx
    src/RDefineBase.cxx
    src/RCutFlowReport.cxx
    src/RDataFrame.cxx
    src/RDefineReader.cxx
    src/RDFActionHelpers.cxx
    src/RDFBookedDefines.cxx
    src/RDFDisplay.cxx
    src/RDFGraphUtils.cxx
    src/RDFHistoModels.cxx
    src/RDFInterfaceUtils.cxx
    src/RDFUtils.cxx
    src/RDFHelpers.cxx
    src/RFilterBase.cxx
    src/RJittedAction.cxx
    src/RJittedDefine.cxx
    src/RJittedFilter.cxx
    src/RLoopManager.cxx
    src/RRangeBase.cxx
    src/RRootDS.cxx
    src/RSlotStack.cxx
    src/RTrivialDS.cxx
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
    ${RDATAFRAME_EXTRA_INCLUDES}
  DEPENDENCIES
    Tree
    TreePlayer
    Hist
    RIO
    ROOTVecOps
    ${RDATAFRAME_EXTRA_DEPS}
)

if(arrow)
  target_sources(ROOTDataFrame PRIVATE src/RArrowDS.cxx)
  target_include_directories(ROOTDataFrame PRIVATE ${ARROW_INCLUDE_DIR})
  target_link_libraries(ROOTDataFrame PRIVATE ${ARROW_SHARED_LIB})
endif()

if(sqlite)
  target_sources(ROOTDataFrame PRIVATE src/RSqliteDS.cxx)
  target_include_directories(ROOTDataFrame PRIVATE ${SQLITE_INCLUDE_DIR})
  target_link_libraries(ROOTDataFrame PRIVATE ${SQLITE_LIBRARIES})
endif()

if(root7)
  target_sources(ROOTDataFrame PRIVATE src/RNTupleDS.cxx)
endif(root7)

if(MSVC)
  target_compile_definitions(ROOTDataFrame PRIVATE _USE_MATH_DEFINES)
endif()

ROOT_ADD_TEST_SUBDIRECTORY(test)
