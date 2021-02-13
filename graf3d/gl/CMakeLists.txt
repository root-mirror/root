# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT graf3d/gl package
############################################################################

if(WIN32 OR cocoa)
  set(installoptions FILTER "TX11GL")
endif()

if(x11)
  set(RGL_EXTRA_HEADERS TX11GL.h)
  set(RGL_EXTRA_SOURCES TX11GL.cxx)
endif()


if(builtin_gl2ps)
  set(RGL_EXTRA_SOURCES ${RGL_EXTRA_SOURCES} src/gl2ps.cxx)
endif()

set_source_files_properties(src/TGLFontManager.cxx PROPERTIES COMPILE_FLAGS "${FTGL_CFLAGS}")
set_source_files_properties(src/TGLText.cxx PROPERTIES COMPILE_FLAGS "${FTGL_CFLAGS}")

ROOT_STANDARD_LIBRARY_PACKAGE(RGL
  HEADERS
    TArcBall.h
    TF2GL.h
    TGL5DDataSetEditor.h
    TGL5D.h
    TGLAdapter.h
    TGLAnnotation.h
    TGLAutoRotator.h
    TGLAxis.h
    TGLAxisPainter.h
    TGLBoundingBox.h
    TGLBoxPainter.h
    TGLCameraGuide.h
    TGLCamera.h
    TGLCameraOverlay.h
    TGLClip.h
    TGLClipSetEditor.h
    TGLContext.h
    TGLCylinder.h
    TGLEmbeddedViewer.h
    TGLEventHandler.h
    TGLFaceSet.h
    TGLFBO.h
    TGLFontManager.h
    TGLFormat.h
    TGLH2PolyPainter.h
    TGLHistPainter.h
    TGLLegoPainter.h
    TGLLightSetEditor.h
    TGLLightSet.h
    TGLLockable.h
    TGLLogicalShape.h
    TGLManip.h
    TGLManipSet.h
    TGLObject.h
    TGLOrthoCamera.h
    TGLOutput.h
    TGLOverlayButton.h
    TGLOverlay.h
    TGLPadPainter.h
    TGLPadUtils.h
    TGLParametricEquationGL.h
    TGLParametric.h
    TGLPerspectiveCamera.h
    TGLPhysicalShape.h
    TGLPlot3D.h
    TGLPlotBox.h
    TGLPlotCamera.h
    TGLPlotPainter.h
    TGLPolyLine.h
    TGLPolyMarker.h
    TGLPShapeObjEditor.h
    TGLPShapeObj.h
    TGLPShapeRef.h
    TGLQuadric.h
    TGLRnrCtx.h
    TGLRotateManip.h
    TGLSAFrame.h
    TGLSAViewer.h
    TGLScaleManip.h
    TGLSceneBase.h
    TGLScene.h
    TGLSceneInfo.h
    TGLScenePad.h
    TGLSelectBuffer.h
    TGLSelectRecord.h
    TGLSphere.h
    TGLStopwatch.h
    TGLSurfacePainter.h
    TGLText.h
    TGLTF3Painter.h
    TGLTH3Composition.h
    TGLTransManip.h
    TGLUtil.h
    TGLViewerBase.h
    TGLViewerEditor.h
    TGLViewer.h
    TGLVoxelPainter.h
    TGLWidget.h
    TH2GL.h
    TH3GL.h
    TPointSet3DGL.h
    ${RGL_EXTRA_HEADERS}
  SOURCES
    src/TArcBall.cxx
    src/TF2GL.cxx
    src/TGL5D.cxx
    src/TGL5DDataSetEditor.cxx
    src/TGL5DPainter.cxx
    src/TGLAdapter.cxx
    src/TGLAnnotation.cxx
    src/TGLAutoRotator.cxx
    src/TGLAxis.cxx
    src/TGLAxisPainter.cxx
    src/TGLBoundingBox.cxx
    src/TGLBoxPainter.cxx
    src/TGLCamera.cxx
    src/TGLCameraGuide.cxx
    src/TGLCameraOverlay.cxx
    src/TGLClip.cxx
    src/TGLClipSetEditor.cxx
    src/TGLContext.cxx
    src/TGLContextPrivate.cxx
    src/TGLContextPrivate.h
    src/TGLCylinder.cxx
    src/TGLEmbeddedViewer.cxx
    src/TGLEventHandler.cxx
    src/TGLFaceSet.cxx
    src/TGLFBO.cxx
    src/TGLFontManager.cxx
    src/TGLFormat.cxx
    src/TGLH2PolyPainter.cxx
    src/TGLHistPainter.cxx
    src/TGLIsoMesh.cxx
    src/TGLLegoPainter.cxx
    src/TGLLightSet.cxx
    src/TGLLightSetEditor.cxx
    src/TGLLockable.cxx
    src/TGLLogicalShape.cxx
    src/TGLManip.cxx
    src/TGLManipSet.cxx
    src/TGLMarchingCubes.cxx
    src/TGLObject.cxx
    src/TGLOrthoCamera.cxx
    src/TGLOutput.cxx
    src/TGLOverlayButton.cxx
    src/TGLOverlay.cxx
    src/TGLPadPainter.cxx
    src/TGLPadUtils.cxx
    src/TGLParametric.cxx
    src/TGLParametricEquationGL.cxx
    src/TGLPerspectiveCamera.cxx
    src/TGLPhysicalShape.cxx
    src/TGLPlot3D.cxx
    src/TGLPlotBox.cxx
    src/TGLPlotCamera.cxx
    src/TGLPlotPainter.cxx
    src/TGLPolyLine.cxx
    src/TGLPolyMarker.cxx
    src/TGLPShapeObj.cxx
    src/TGLPShapeObjEditor.cxx
    src/TGLPShapeRef.cxx
    src/TGLQuadric.cxx
    src/TGLRnrCtx.cxx
    src/TGLRotateManip.cxx
    src/TGLSAFrame.cxx
    src/TGLSAViewer.cxx
    src/TGLScaleManip.cxx
    src/TGLSceneBase.cxx
    src/TGLScene.cxx
    src/TGLSceneInfo.cxx
    src/TGLScenePad.cxx
    src/TGLSelectBuffer.cxx
    src/TGLSelectRecord.cxx
    src/TGLSphere.cxx
    src/TGLStopwatch.cxx
    src/TGLSurfacePainter.cxx
    src/TGLText.cxx
    src/TGLTF3Painter.cxx
    src/TGLTH3Composition.cxx
    src/TGLTransManip.cxx
    src/TGLUtil.cxx
    src/TGLViewerBase.cxx
    src/TGLViewer.cxx
    src/TGLViewerEditor.cxx
    src/TGLVoxelPainter.cxx
    src/TGLWidget.cxx
    src/TH2GL.cxx
    src/TH3GL.cxx
    src/TKDEAdapter.cxx
    src/TKDEFGT.cxx
    src/TPointSet3DGL.cxx
    ${RGL_EXTRA_SOURCES}
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  LIBRARIES
    OpenGL::GL
    OpenGL::GLU
    ${FTGL_LIBRARIES}
    ${GL2PS_LIBRARIES}
    ${X11_LIBRARIES}
    RGlew
  DEPENDENCIES
    Hist
    Gui
    Ged
    RCsg
  INSTALL_OPTIONS
    ${installoptions}
)

target_include_directories(RGL PRIVATE
     ${OPENGL_INCLUDE_DIR}
     ${FTGL_INCLUDE_DIR}
     ${FREETYPE_INCLUDE_DIRS}
)

if(builtin_gl2ps)
  target_include_directories(RGL PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/gl2ps)
else()
  target_include_directories(RGL PRIVATE ${GL2PS_INCLUDE_DIRS})
endif()


if(MSVC)
  target_compile_definitions (RGL PRIVATE "GLEW_STATIC")
endif()

if(MACOSX_GLU_DEPRECATED)
  target_compile_options(RGL PRIVATE -Wno-deprecated-declarations)
endif()
