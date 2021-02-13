# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

add_subdirectory(proof)
add_subdirectory(proofplayer) # special CMakeListst.txt
if(NOT WIN32)
  add_subdirectory(proofbench) # special CMakeListst.txt
endif()
if(xproofd AND xrootd AND ssl)
  add_subdirectory(proofd)
  add_subdirectory(proofx)
endif()
