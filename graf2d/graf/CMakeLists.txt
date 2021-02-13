# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT graf2d/graf package
# @author Pere Mato, CERN
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(Graf
  HEADERS
    TArc.h
    TArrow.h
    TAttImage.h
    TBox.h
    TCandle.h
    TCrown.h
    TCurlyArc.h
    TCurlyLine.h
    TCutG.h
    TDiamond.h
    TEllipse.h
    TFrame.h
    TGaxis.h
    TGraphPolar.h
    TGraphPolargram.h
    TGraphQQ.h
    TImage.h
    TImagePlugin.h
    TLatex.h
    TLegend.h
    TLegendEntry.h
    TLine.h
    TLink.h
    TMarker.h
    TMathText.h
    TPave.h
    TPaveLabel.h
    TPaveStats.h
    TPaveText.h
    TPavesText.h
    TPie.h
    TPieSlice.h
    TPoints.h
    TPolyLine.h
    TTF.h
    TText.h
    TWbox.h
  SOURCES
    src/TArc.cxx
    src/TArrow.cxx
    src/TAttImage.cxx
    src/TBox.cxx
    src/TCandle.cxx
    src/TCrown.cxx
    src/TCurlyArc.cxx
    src/TCurlyLine.cxx
    src/TCutG.cxx
    src/TDiamond.cxx
    src/TEllipse.cxx
    src/TFrame.cxx
    src/TGaxis.cxx
    src/TGraphPolar.cxx
    src/TGraphPolargram.cxx
    src/TGraphQQ.cxx
    src/TImage.cxx
    src/TLatex.cxx
    src/TLegend.cxx
    src/TLegendEntry.cxx
    src/TLine.cxx
    src/TLink.cxx
    src/TMarker.cxx
    src/TMathText.cxx
    src/TPave.cxx
    src/TPaveLabel.cxx
    src/TPaveStats.cxx
    src/TPaveText.cxx
    src/TPavesText.cxx
    src/TPie.cxx
    src/TPieSlice.cxx
    src/TPoints.cxx
    src/TPolyLine.cxx
    src/TTF.cxx
    src/TText.cxx
    src/TWbox.cxx
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  LIBRARIES
    ${FREETYPE_LIBRARIES}
    ZLIB::ZLIB
    mathtext
  DEPENDENCIES
    Hist
    Matrix
    MathCore
    RIO
  BUILTINS
    FREETYPE
)

target_include_directories(Graf PRIVATE ${FREETYPE_INCLUDE_DIRS})
