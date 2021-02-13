# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT core/rint package
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(Rint
  STAGE1
  HEADERS
    TRint.h
    TTabCom.h
  SOURCES
    src/TRint.cxx
    src/TTabCom.cxx
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  DEPENDENCIES
    Core
)

ROOT_ADD_TEST_SUBDIRECTORY(test)
