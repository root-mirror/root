# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT graf2d/quartz package
# @author Pere Mato, CERN
############################################################################

add_definitions("-ObjC++")

# Dictionary not needed.
ROOT_LINKER_LIBRARY(GQuartz
    src/QuartzFillArea.mm
    src/QuartzLine.mm
    src/QuartzMarker.mm
    src/QuartzText.mm
    src/QuartzUtils.mm
  LIBRARIES
    "-framework Cocoa"
  DEPENDENCIES
    Core
    MathCore
)
