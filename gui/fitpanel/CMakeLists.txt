# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT gui/fitpanel package
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(FitPanel
  HEADERS
    TAdvancedGraphicsDialog.h
    TFitEditor.h
    TFitParametersDialog.h
    TTreeInput.h
  SOURCES
    src/TAdvancedGraphicsDialog.cxx
    src/TFitEditor.cxx
    src/TFitParametersDialog.cxx
    src/TTreeInput.cxx
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  DEPENDENCIES
    Gpad
    Graf
    Gui
    Hist
    MathCore
    Tree
    TreePlayer
)
