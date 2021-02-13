# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

if(MSVC AND MSVC_VERSION GREATER_EQUAL 1923)
  # FIXME: using /O2 compiler flag prevent the following error when building in debug mode:
  # axis.obj : fatal error LNK1179: invalid or corrupt file: duplicate COMDAT
  # '??$?8DU?$char_traits@D@std@@@__ROOT@experimental@std@@YA_NV?$basic_string_view@DU?$char_traits@D@std@@@012@0@Z'
  # Try to remove those lines when upgrading Visual Studio
  string(REPLACE "-Od -Z7" "-O2" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
endif()

ROOT_ADD_UNITTEST_DIR(ROOTHist)
