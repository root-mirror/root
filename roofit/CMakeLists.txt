# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

add_subdirectory(batchcompute)
add_subdirectory(builtin_zeromq)
add_subdirectory(roofitZMQ)
add_subdirectory(multiprocess)
add_subdirectory(roofitcore)
add_subdirectory(roofit)
if(mathmore)
  add_subdirectory(roofitmore)
endif()
add_subdirectory(roostats)
if(xml)
  add_subdirectory(histfactory)
endif()
add_subdirectory(RDataFrameHelpers)

