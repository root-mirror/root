# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT io/gfal package
############################################################################

add_definitions(-D_FILE_OFFSET_BITS=64)

ROOT_STANDARD_LIBRARY_PACKAGE(GFAL
                              HEADERS TGFALFile.h
                              SOURCES src/TGFALFile.cxx
                              LIBRARIES ${GFAL_LIBRARIES}
                              DEPENDENCIES Core Net RIO)

target_include_directories(GFAL PRIVATE ${GFAL_INCLUDE_DIRS})
