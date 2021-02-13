# Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT gui/browsable
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(ROOTBrowsable
  HEADERS
    ROOT/Browsable/RAnyObjectHolder.hxx
    ROOT/Browsable/RElement.hxx
    ROOT/Browsable/RGroup.hxx
    ROOT/Browsable/RHolder.hxx
    ROOT/Browsable/RItem.hxx
    ROOT/Browsable/RLevelIter.hxx
    ROOT/Browsable/RProvider.hxx
    ROOT/Browsable/RShared.hxx
    ROOT/Browsable/RSysFile.hxx
    ROOT/Browsable/RSysFileItem.hxx
    ROOT/Browsable/RUnique.hxx
    ROOT/Browsable/RWrapper.hxx
    ROOT/Browsable/TKeyItem.hxx
    ROOT/Browsable/TObjectElement.hxx
    ROOT/Browsable/TObjectHolder.hxx
    ROOT/Browsable/TObjectItem.hxx
  SOURCES
    src/RElement.cxx
    src/RGroup.cxx
    src/RLevelIter.cxx
    src/RProvider.cxx
    src/RSysFile.cxx
    src/TDirectoryElement.cxx
    src/TObjectElement.cxx
    src/TObjectHolder.cxx
  DEPENDENCIES
    RIO
    Hist
)

# provider for drawing of v7 histograms on RCanvas

ROOT_LINKER_LIBRARY(ROOTHistDrawProvider
     src/RHistDraw7Provider.cxx
  DEPENDENCIES
     ROOTBrowsable
     ROOTGpadv7
     ROOTHistDraw
)

# provider for generic drawing of TObject on TCanvas

ROOT_LINKER_LIBRARY(ROOTObjectDraw6Provider
     src/TObjectDraw6Provider.cxx
  DEPENDENCIES
     ROOTBrowsable
     Gpad
     Hist
)

# provider for generic drawing of TObject on RCanvas

ROOT_LINKER_LIBRARY(ROOTObjectDraw7Provider
     src/TObjectDraw7Provider.cxx
  DEPENDENCIES
     ROOTBrowsable
     ROOTGpadv7
)

# provider for browsing of TBranchElement

ROOT_LINKER_LIBRARY(ROOTBranchBrowseProvider
     src/TBranchBrowseProvider.cxx
  DEPENDENCIES
     ROOTBrowsable
     Tree
)

# provider for drawing of TLeaf on TCanvas

ROOT_LINKER_LIBRARY(ROOTLeafDraw6Provider
     src/TLeafDraw6Provider.cxx
  DEPENDENCIES
     ROOTBrowsable
     Tree
     Hist
     Gpad
)

# provider for drawing of TLeaf on RCanvas

ROOT_LINKER_LIBRARY(ROOTLeafDraw7Provider
     src/TLeafDraw7Provider.cxx
  DEPENDENCIES
     ROOTBrowsable
     Tree
     Hist
     ROOTGpadv7
)

