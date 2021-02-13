# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT hist/histpainter package
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(HistPainter
  HEADERS
    Hoption.h
    Hparam.h
    TGraph2DPainter.h
    TGraphPainter.h
    THistPainter.h
    TPaletteAxis.h
  SOURCES
    src/TGraph2DPainter.cxx
    src/TGraphPainter.cxx
    src/THistPainter.cxx
    src/TPainter3dAlgorithms.cxx
    src/TPaletteAxis.cxx
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  DEPENDENCIES
    Gpad
    Graf
    Hist
    MathCore
    Matrix
)
