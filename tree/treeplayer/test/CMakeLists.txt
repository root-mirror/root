# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

if(NOT MSVC OR win_broken_tests)
  ROOT_ADD_GTEST(treetreeplayertestUnit array.cxx leafs.cxx readerarray_iterator.cxx
    LIBRARIES Core TreePlayer
    COPY_TO_BUILDDIR ${CMAKE_CURRENT_SOURCE_DIR}/data.h)
  ROOT_ADD_GTEST(treetreeplayertestUnit2 basic.cxx LIBRARIES TreePlayer FAILREGEX "")
  #  add_custom_command(TARGET treetreeplayertestUnit POST_BUILD
  #                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/data.h data.h)
endif()

if(imt)
   ROOT_ADD_GTEST(treeprocessormt treeprocmt/treeprocessormt.cxx LIBRARIES TreePlayer)
   if(xrootd)
      ROOT_ADD_GTEST(treeprocessormt_remotefiles treeprocmt/treeprocessormt_remotefiles.cxx LIBRARIES TreePlayer)
   endif()
endif()
