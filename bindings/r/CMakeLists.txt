# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT R package
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(RInterface
  HEADERS
    RExports.h
    TRDataFrame.h
    TRFunctionExport.h
    TRFunctionImport.h
    TRInterface.h
    TRInternalFunction.h
    TRObject.h
  SOURCES
    src/RExports.cxx
    src/TRCompletion.cxx
    src/TRDataFrame.cxx
    src/TRFunctionExport.cxx
    src/TRFunctionImport.cxx
    src/TRInterface.cxx
    src/TRObject.cxx
  DEPENDENCIES
    Core
    Matrix
    RIO
    Thread
  NO_CXXMODULE
)

ROOT_ADD_CXX_FLAG(_R_FLAGS -Wno-cast-function-type)
ROOT_ADD_CXX_FLAG(_R_FLAGS -Wno-overloaded-virtual)
if(x11)
  ROOT_ADD_CXX_FLAG(_R_FLAGS -DHAS_X11)
endif()

separate_arguments(_R_FLAGS)

target_compile_options(RInterface PUBLIC ${_R_FLAGS})
target_include_directories(RInterface BEFORE PUBLIC ${R_INCLUDE_DIRS})
if(x11)
  target_link_libraries(RInterface PUBLIC ${R_LIBRARIES} ${X11_LIBRARIES} readline)
elseif()
  target_link_libraries(RInterface PUBLIC ${R_LIBRARIES} readline)
endif()