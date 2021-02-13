# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT tree/treeplayer package
# @author Pere Mato CERN, Danilo Piparo CERN
############################################################################

if(NOT MSVC)
  list(APPEND TREEPLAYER_EXTRA_HEADERS
    TMPWorkerTree.h
    ROOT/TTreeProcessorMP.hxx
  )
endif()

if(imt)
  list(APPEND TREEPLAYER_EXTRA_HEADERS ROOT/TTreeProcessorMT.hxx)
  list(APPEND TREEPLAYER_EXTRA_DEPENDENCIES Imt)
endif()

ROOT_STANDARD_LIBRARY_PACKAGE(TreePlayer
  HEADERS
    ROOT/TTreeReaderFast.hxx
    ROOT/TTreeReaderValueFast.hxx
    TBranchProxyClassDescriptor.h
    TBranchProxyDescriptor.h
    TBranchProxyDirector.h
    TBranchProxy.h
    TChainIndex.h
    TFileDrawMap.h
    TFormLeafInfo.h
    TFormLeafInfoReference.h
    TFriendProxyDescriptor.h
    TFriendProxy.h
    TRefArrayProxy.h
    TRefProxy.h
    TSelectorDraw.h
    TSelectorEntries.h
    TSimpleAnalysis.h
    TTreeDrawArgsParser.h
    TTreeFormula.h
    TTreeFormulaManager.h
    TTreeGeneratorBase.h
    TTreeIndex.h
    TTreePerfStats.h
    TTreePlayer.h
    TTreeProxyGenerator.h
    TTreeReaderArray.h
    TTreeReaderGenerator.h
    TTreeReader.h
    TTreeReaderUtils.h
    TTreeReaderValue.h
    TTreeTableInterface.h
    ${TREEPLAYER_EXTRA_HEADERS}
  SOURCES
    src/TBranchProxyClassDescriptor.cxx
    src/TBranchProxy.cxx
    src/TBranchProxyDescriptor.cxx
    src/TBranchProxyDirector.cxx
    src/TChainIndex.cxx
    src/TFileDrawMap.cxx
    src/TFormLeafInfo.cxx
    src/TFormLeafInfoReference.cxx
    src/TFriendProxy.cxx
    src/TFriendProxyDescriptor.cxx
    src/TRefArrayProxy.cxx
    src/TRefProxy.cxx
    src/TSelectorDraw.cxx
    src/TSelectorEntries.cxx
    src/TSimpleAnalysis.cxx
    src/TTreeDrawArgsParser.cxx
    src/TTreeFormula.cxx
    src/TTreeFormulaManager.cxx
    src/TTreeGeneratorBase.cxx
    src/TTreeIndex.cxx
    src/TTreePerfStats.cxx
    src/TTreePlayer.cxx
    src/TTreeProxyGenerator.cxx
    src/TTreeReaderArray.cxx
    src/TTreeReader.cxx
    src/TTreeReaderFast.cxx
    src/TTreeReaderGenerator.cxx
    src/TTreeReaderValue.cxx
    src/TTreeReaderValueFast.cxx
    src/TTreeTableInterface.cxx
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  DEPENDENCIES
    Gpad
    Graf
    Graf3d
    Hist
    ${TREEPLAYER_EXTRA_DEPENDENCIES}
    MathCore
    RIO
    Tree
)

if(NOT MSVC)
  target_sources(TreePlayer PRIVATE
    src/TMPWorkerTree.cxx
    src/TTreeProcessorMP.cxx
  )
  target_link_libraries(TreePlayer PUBLIC MultiProc)
endif()

if(imt)
  target_sources(TreePlayer PRIVATE src/TTreeProcessorMT.cxx)
endif()

ROOT_ADD_TEST_SUBDIRECTORY(test)
