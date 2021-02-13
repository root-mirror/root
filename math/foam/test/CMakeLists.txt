# Testing TFoam
# @author Stephan Hageboeck, CERN, 2020

ROOT_ADD_GTEST(testTFoam testTFoam.cxx
  LIBRARIES Tree Foam
  COPY_TO_BUILDDIR ${CMAKE_CURRENT_SOURCE_DIR}/testTFoam_1.root)

