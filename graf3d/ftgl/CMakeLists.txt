# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT graf3d/ftgl package
############################################################################

if(MACOSX_GLU_DEPRECATED)
  add_definitions(-Wno-deprecated-declarations)
endif()

ROOT_LINKER_LIBRARY(FTGL
    src/FTBitmapGlyph.cxx
    src/FTCharmap.cxx
    src/FTContour.cxx
    src/FTExtrdGlyph.cxx
    src/FTFace.cxx
    src/FTFont.cxx
    src/FTGLBitmapFont.cxx
    src/FTGLExtrdFont.cxx
    src/FTGLOutlineFont.cxx
    src/FTGLPixmapFont.cxx
    src/FTGLPolygonFont.cxx
    src/FTGLTextureFont.cxx
    src/FTGlyphContainer.cxx
    src/FTGlyph.cxx
    src/FTLibrary.cxx
    src/FTOutlineGlyph.cxx
    src/FTPixmapGlyph.cxx
    src/FTPoint.cxx
    src/FTPolyGlyph.cxx
    src/FTSize.cxx
    src/FTTextureGlyph.cxx
    src/FTVectoriser.cxx
  LIBRARIES
    OpenGL::GL
    OpenGL::GLU
    ${FREETYPE_LIBRARIES}
    ZLIB::ZLIB
    RGlew
  BUILTINS
    FREETYPE
)

target_include_directories(FTGL PRIVATE
      ${FREETYPE_INCLUDE_DIRS}
      ${OPENGL_INCLUDE_DIR}
      ${CMAKE_SOURCE_DIR}/graf3d/gl/inc
)

ROOT_INSTALL_HEADERS()
