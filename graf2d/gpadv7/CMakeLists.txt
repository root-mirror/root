# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT graf2d/gpadv7 package
# @author Pere Mato, CERN
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(ROOTGpadv7
  HEADERS
    ROOT/RCanvas.hxx
    ROOT/RFrame.hxx
    ROOT/RMenuItems.hxx
    ROOT/RColor.hxx
    ROOT/RDisplayItem.hxx
    ROOT/RAttrMap.hxx
    ROOT/RAttrBase.hxx
    ROOT/RAttrAxis.hxx
    ROOT/RAttrBox.hxx
    ROOT/RAttrColor.hxx
    ROOT/RAttrLine.hxx
    ROOT/RAttrFill.hxx
    ROOT/RAttrMarker.hxx
    ROOT/RAttrMargins.hxx
    ROOT/RAttrText.hxx
    ROOT/RAttrValue.hxx
    ROOT/RAxisDrawable.hxx
    ROOT/RPalette.hxx
    ROOT/RPaletteDrawable.hxx
    ROOT/RDrawable.hxx
    ROOT/RDrawableRequest.hxx
    ROOT/RStyle.hxx
    ROOT/RPadDisplayItem.hxx
    ROOT/RPadExtent.hxx
    ROOT/RPadBase.hxx
    ROOT/RPad.hxx
    ROOT/RPadLength.hxx
    ROOT/RPadPos.hxx
    ROOT/RPadUserAxis.hxx
    ROOT/RPave.hxx
    ROOT/RVirtualCanvasPainter.hxx
    ROOT/TObjectDrawable.hxx
  SOURCES
    src/RCanvas.cxx
    src/RFrame.cxx
    src/RMenuItems.cxx
    src/RColor.cxx
    src/RDisplayItem.cxx
    src/RDrawable.cxx
    src/RDrawableRequest.cxx
    src/RAttrMap.cxx
    src/RAttrBase.cxx
    src/RPalette.cxx
    src/RStyle.cxx
    src/RPadBase.cxx
    src/RPad.cxx
    src/RPadExtent.cxx
    src/RPadLength.cxx
    src/RPadPos.cxx
    src/RPadUserAxis.cxx
    src/RVirtualCanvasPainter.cxx
    src/TObjectDrawable.cxx
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
  DEPENDENCIES
    Core
    RIO
)

ROOT_ADD_TEST_SUBDIRECTORY(test)
