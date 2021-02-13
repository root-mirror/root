# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT graph3d/rglew package
############################################################################

# This package is only needed for user backward compatibility! (TGLIncludes.h)

add_library(RGlew INTERFACE)
target_include_directories(RGlew INTERFACE inc/)
target_link_libraries(RGlew INTERFACE GLEW::GLEW)
#target_include_directories(RGlew INTERFACE $<BUILD_INTERFACE:${GLEW_INCLUDE_DIR}>)

# We still need to install TGLIncludes.h into include/
ROOT_INSTALL_HEADERS()
