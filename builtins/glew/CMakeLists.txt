# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

project(GLEW C)

find_package(OpenGL REQUIRED)

include_directories (${CMAKE_CURRENT_SOURCE_DIR}/inc)

set(GLEW_HEADERS
  inc/GL/glew.h
  inc/GL/glxew.h
  inc/GL/wglew.h
)

set(GLEW_SOURCES
  src/glew.c
)

unset(GLEW_FOUND CACHE)
unset(GLEW_FOUND PARENT_SCOPE)
set(GLEW_FOUND TRUE CACHE BOOL "" FORCE)

set(GLEW_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/inc CACHE INTERNAL "" FORCE)
set(GLEW_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/inc CACHE INTERNAL "" FORCE)

if(APPLE)
  if(NOT cocoa)
    set(GLEW_DEFINITIONS -DGLEW_APPLE_GLX)
  endif()
endif()

if(NOT MSVC)
  add_library(GLEW SHARED ${GLEW_HEADERS} ${GLEW_SOURCES})
  target_compile_options(GLEW PRIVATE -fPIC)
else()
  add_library(GLEW STATIC ${GLEW_HEADERS} ${GLEW_SOURCES})
  set_target_properties(GLEW PROPERTIES COMPILE_DEFINITIONS "GLEW_STATIC" PREFIX "lib")
  target_compile_definitions (GLEW PRIVATE "GLEW_STATIC")
endif()
target_compile_options(GLEW PRIVATE ${GLEW_DEFINITIONS})
target_include_directories(GLEW INTERFACE $<BUILD_INTERFACE:${GLEW_INCLUDE_DIR}>)
target_link_libraries(GLEW PRIVATE OpenGL::GL OpenGL::GLU)

add_library(GLEW::GLEW ALIAS GLEW)

set(GLEW_LIBRARY $<TARGET_FILE:GLEW> CACHE INTERNAL "")
set(GLEW_LIBRARIES GLEW::GLEW CACHE INTERNAL "")

set_property(GLOBAL APPEND PROPERTY ROOT_BUILTIN_TARGETS GLEW::GLEW)

install(TARGETS GLEW
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        )

ROOT_INSTALL_HEADERS()
