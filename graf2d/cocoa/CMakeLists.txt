# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT graf2d/cocoa package
# @author Pere Mato, CERN
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(GCocoa
  HEADERS
    TGCocoa.h
    TGOSXGL.h
    TGQuartz.h
  SOURCES
    src/CocoaConstants.mm
    src/CocoaGuiTypes.mm
    src/CocoaPrivate.mm
    src/FontCache.mm
    src/MenuLoader.mm
    src/QuartzPixmap.mm
    src/QuartzWindow.mm
    src/ROOTApplicationDelegate.mm
    src/ROOTOpenGLView.mm
    src/TGCocoa.mm
    src/TGOSXGL.mm
    src/TGQuartz.mm
    src/X11Atoms.mm
    src/X11Buffer.mm
    src/X11Colors.mm
    src/X11Events.mm
    src/XLFDParser.mm
  DEPENDENCIES
    Gui
  LIBRARIES
    GQuartz
    "-framework Cocoa"
    "-framework OpenGL"
    ${FREETYPE_LIBRARIES}
    RGlew
)

target_include_directories(GCocoa PRIVATE ${CMAKE_SOURCE_DIR}/graf3d/gl/inc)

# We know GL is deprecated on macOS. And our answer is WebGUI/Graphics.
target_compile_options(GCocoa PRIVATE -ObjC++ -Wno-deprecated-declarations)
target_include_directories(GCocoa PRIVATE ${FREETYPE_INCLUDE_DIRS})

if(opengl) # special case when cocoa is enabled but not opengl (i.e. gminimal=ON)
  target_include_directories(GCocoa PRIVATE ${OPENGL_INCLUDE_DIR})
else()
  target_include_directories(GCocoa PRIVATE
    ${CMAKE_SOURCE_DIR}/graf3d/gl/inc
  )
endif(opengl)

