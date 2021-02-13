# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

add_subdirectory(tmva)
add_subdirectory(tmvagui)

# Build pyMVA. Allows users to call python mva's from c++ ROOT/TMVA.
if(tmva-pymva AND NUMPY_FOUND)
   add_subdirectory(pymva)
endif()

# Build RMVA. Allows users to call R mva's from c++ ROOT/TMVA.
if(r OR tmva-rmva)
   add_subdirectory(rmva)
endif()
