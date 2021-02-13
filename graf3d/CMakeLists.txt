# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

if(cxxmodules AND APPLE)
  # FIXME: Any header of graf3d can trigger building of the OSX OpenGL module.
  # It looks like our glew infrastructure combined with the OpenGL modulemap
  # cannot be compiled with -fmodules-local-submodule-visibility yet.
  # Strip out that flag.
  string(REPLACE "-Xclang -fmodules-local-submodule-visibility" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
endif()

add_subdirectory(g3d) # special CMakeLists.txt
if(NOT WIN32 AND x11)
  add_subdirectory(x3d) # special CMakeLists.txt
endif()

add_subdirectory(csg) # special CMakeLists.txt

if (opengl)
  if(builtin_ftgl)
    add_subdirectory(ftgl)
  endif()
  add_subdirectory(rglew)
  add_subdirectory(eve)
  add_subdirectory(gl)
  add_subdirectory(gviz3d)
endif()

if(webgui AND root7)
  add_subdirectory(eve7) # special CMakeLists.txt
endif()
