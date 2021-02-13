ROOT_ADD_GTEST(testHypoTestInvResult testHypoTestInvResult.cxx
  LIBRARIES RooStats
  COPY_TO_BUILDDIR ${CMAKE_CURRENT_SOURCE_DIR}/testHypoTestInvResult_1.root)
ROOT_ADD_GTEST(testSPlot testSPlot.cxx LIBRARIES RooStats)
