# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

add_subdirectory(net)
add_subdirectory(auth)

if(NOT WIN32 AND ssl)
  add_subdirectory(rpdutils)
endif()

if(netx)
  add_subdirectory(netx)
endif()

if(alien)
  add_subdirectory(alien)
endif()

if(monalisa)
  add_subdirectory(monalisa)
endif()

if(davix)
   add_subdirectory(davix)
endif()

if(netxng)
   add_subdirectory(netxng)
endif()

if(http)
   add_subdirectory(http)
   add_subdirectory(httpsniff)
endif()
