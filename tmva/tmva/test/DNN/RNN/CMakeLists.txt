# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building TMVA/DNN/RNN tests.
# @author Saurav Shekhar
############################################################################


set(Libraries TMVA)

if (CUDA_FOUND AND  tmva-gpu)
include_directories(${CUDA_INCLUDE_DIRS})
endif()

#---Reference test

if (Test_Reference)

# RNN - BackPropagation Reference
ROOT_EXECUTABLE(testRecurrentBackpropagation TestRecurrentBackpropagation.cxx LIBRARIES ${Libraries})
ROOT_ADD_TEST(TMVA-DNN-RNN-Backpropagation COMMAND testRecurrentBackpropagation)
#
# RNN - Forward Reference
ROOT_EXECUTABLE(testRecurrentForwardPass TestRecurrentForwardPass.cxx LIBRARIES ${Libraries})
ROOT_ADD_TEST(TMVA-DNN-RNN-Forward COMMAND testRecurrentForwardPass)

# RNN - Full Test Reference
ROOT_EXECUTABLE(testFullRNN TestFullRNN.cxx LIBRARIES ${Libraries})
ROOT_ADD_TEST(TMVA-DNN-RNN-FullRNN COMMAND testFullRNN)

endif()

#--- CUDA tests. ---------------------------
if (CUDA_FOUND AND tmva-gpu)

set(DNN_CUDA_LIBRARIES ${CUDA_CUBLAS_LIBRARIES} ${CUDNN_LIBRARIES} )

CUDA_ADD_EXECUTABLE(testRecurrentBackpropagationCuda TestRecurrentBackpropagationCuda.cxx)
target_link_libraries(testRecurrentBackpropagationCuda ${Libraries} ${DNN_CUDA_LIBRARIES})
ROOT_ADD_TEST(TMVA-DNN-RNN-BackpropagationCuda COMMAND testRecurrentBackpropagationCuda)
#mark the native Cuda test failing since RNN are supported only with Cudnn
set_tests_properties(TMVA-DNN-RNN-BackpropagationCuda PROPERTIES WILL_FAIL true)

if (tmva-cudnn)

CUDA_ADD_EXECUTABLE(testRecurrentForwardPassCudnn TestRecurrentForwardPassCudnn.cxx)
target_link_libraries(testRecurrentForwardPassCudnn ${Libraries} ${DNN_CUDA_LIBRARIES})
ROOT_ADD_TEST(TMVA-DNN-RNN-Forwaed-Cudnn COMMAND testRecurrentForwardPassCudnn)

CUDA_ADD_EXECUTABLE(testRecurrentBackpropagationCudnn TestRecurrentBackpropagationCudnn.cxx)
target_link_libraries(testRecurrentBackpropagationCudnn ${Libraries} ${DNN_CUDA_LIBRARIES})
ROOT_ADD_TEST(TMVA-DNN-RNN-BackpropagationCudnn COMMAND testRecurrentBackpropagationCudnn)

# Full Test GPU
CUDA_ADD_EXECUTABLE(testFullRNNCudnn TestFullRNNCudnn.cxx)
target_link_libraries(testFullRNNCudnn ${Libraries} ${DNN_CUDA_LIBRARIES})
ROOT_ADD_TEST(TMVA-DNN-RNN-Full-Cudnn COMMAND testFullRNNCudnn)

endif()

endif()

#--- CPU tests. ----------------------------
if ((BLAS_FOUND OR mathmore) AND imt AND tmva-cpu)
  include_directories(SYSTEM ${TBB_INCLUDE_DIRS})

  # DNN - Forward CPU
  ROOT_EXECUTABLE(testRecurrentForwardPassCpu TestRecurrentForwardPassCpu.cxx LIBRARIES ${Libraries})
  ROOT_ADD_TEST(TMVA-DNN-RNN-Forward-Cpu COMMAND testRecurrentForwardPassCpu)

 ROOT_EXECUTABLE(testRecurrentBackpropagationCpu TestRecurrentBackpropagationCpu.cxx LIBRARIES ${Libraries})
 ROOT_ADD_TEST(TMVA-DNN-RNN-Backpropagation-Cpu COMMAND testRecurrentBackpropagationCpu)

# RNN - Full Test Reference
ROOT_EXECUTABLE(testFullRNNCpu TestFullRNNCpu.cxx LIBRARIES ${Libraries})
ROOT_ADD_TEST(TMVA-DNN-RNN-Full-Cpu COMMAND testFullRNNCpu)


endif ()
