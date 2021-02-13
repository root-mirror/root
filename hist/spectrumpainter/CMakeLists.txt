# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT hist/spectrumpainter package
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(SpectrumPainter
  HEADERS
    TSpectrum2Painter.h
  SOURCES
    src/TSpectrum2Painter.cxx
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  DEPENDENCIES
    Graf
    Hist
    Gpad
)
