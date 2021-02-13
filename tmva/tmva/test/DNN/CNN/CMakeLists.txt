# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building TMVA/DNN/CNN tests.
# @author Vladimir Ilievski
############################################################################


set(Libraries TMVA)

#--- CUDA tests. ---------------------------
if (CUDA_FOUND AND tmva-gpu)

set(DNN_CUDA_LIBRARIES ${CUDA_CUBLAS_LIBRARIES} ${CUDNN_LIBRARIES} )

CUDA_ADD_EXECUTABLE(testIm2ColCuda TestIm2ColCuda.cxx)
target_link_libraries(testIm2ColCuda ${Libraries} ${DNN_CUDA_LIBRARIES})
ROOT_ADD_TEST(TMVA-DNN-CNN-Im2ColCuda COMMAND testIm2ColCuda)

CUDA_ADD_EXECUTABLE(testPoolingLayerCuda TestPoolingLayerCuda.cxx)
target_link_libraries(testPoolingLayerCuda ${Libraries} ${DNN_CUDA_LIBRARIES})
ROOT_ADD_TEST(TMVA-DNN-CNN-PoolingLayerCuda COMMAND testPoolingLayerCuda)

CUDA_ADD_EXECUTABLE(testReshapeCuda TestReshapeCuda.cxx)
target_link_libraries(testReshapeCuda ${Libraries} ${DNN_CUDA_LIBRARIES})
ROOT_ADD_TEST(TMVA-DNN-CNN-ReshapeCuda COMMAND testReshapeCuda)

CUDA_ADD_EXECUTABLE(testConvLayerCuda TestConvLayerCuda.cxx)
target_link_libraries(testConvLayerCuda ${Libraries} ${DNN_CUDA_LIBRARIES})
ROOT_ADD_TEST(TMVA-DNN-CNN-ConvLayerCuda COMMAND testConvLayerCuda)

CUDA_ADD_EXECUTABLE(testForwardPassCuda TestForwardPassCuda.cxx)
target_link_libraries(testForwardPassCuda ${Libraries} ${DNN_CUDA_LIBRARIES})
ROOT_ADD_TEST(TMVA-DNN-CNN-ForwardCuda COMMAND testForwardPassCuda)

CUDA_ADD_EXECUTABLE(testRotateWeightsCuda TestRotateWeightsCuda.cxx)
target_link_libraries(testRotateWeightsCuda ${Libraries} ${DNN_CUDA_LIBRARIES})
ROOT_ADD_TEST(TMVA-DNN-CNN-RotateWeightsCuda COMMAND testRotateWeightsCuda)

CUDA_ADD_EXECUTABLE(testConvBackpropagationCuda TestConvBackpropagationCuda.cxx)
target_link_libraries(testConvBackpropagationCuda ${Libraries} ${DNN_CUDA_LIBRARIES})
ROOT_ADD_TEST(TMVA-DNN-CNN-ConvBackpropagationCuda COMMAND testConvBackpropagationCuda)

if (tmva-cudnn)

CUDA_ADD_EXECUTABLE(testForwardPassCudnn TestForwardPassCudnn.cxx)
target_link_libraries(testForwardPassCudnn ${Libraries} ${DNN_CUDA_LIBRARIES})
ROOT_ADD_TEST(TMVA-DNN-CNN-ForwardCudnn COMMAND testForwardPassCudnn)

CUDA_ADD_EXECUTABLE(testConvBackpropagationCudnn TestConvBackpropagationCudnn.cxx)
target_link_libraries(testConvBackpropagationCudnn ${Libraries} ${DNN_CUDA_LIBRARIES})
ROOT_ADD_TEST(TMVA-DNN-CNN-ConvBackpropagationCudnn COMMAND testConvBackpropagationCudnn)

CUDA_ADD_EXECUTABLE(testConvLayerCudnn TestConvLayerCudnn.cxx)
target_link_libraries(testConvLayerCudnn ${Libraries} ${DNN_CUDA_LIBRARIES})
ROOT_ADD_TEST(TMVA-DNN-CNN-ConvLayerCudnn COMMAND testConvLayerCudnn)

CUDA_ADD_EXECUTABLE(testPoolingLayerCudnn TestPoolingLayerCudnn.cxx)
target_link_libraries(testPoolingLayerCudnn ${Libraries} ${DNN_CUDA_LIBRARIES})
ROOT_ADD_TEST(TMVA-DNN-CNN-PoolingLayerCudnn COMMAND testPoolingLayerCudnn)

# test mixed architecture Cudnn Cpu
CUDA_ADD_EXECUTABLE(testMixedArchitectures TestMixedArchitectures.cxx)
target_link_libraries(testMixedArchitectures ${Libraries} ${DNN_CUDA_LIBRARIES})
ROOT_ADD_TEST(TMVA-DNN-CNN-MixedArchitectures COMMAND testMixedArchitectures)

endif()

include_directories(${CUDA_INCLUDE_DIRS})

endif()

#-- -test of reference architecture
if (Test_Reference)

ROOT_EXECUTABLE(testIm2Col TestIm2Col.cxx LIBRARIES ${Libraries})
ROOT_ADD_TEST(TMVA-DNN-CNN-Im2Col COMMAND testIm2Col)

ROOT_EXECUTABLE(testPoolingLayer TestPoolingLayer.cxx LIBRARIES ${Libraries})
ROOT_ADD_TEST(TMVA-DNN-CNN-PoolingLayer COMMAND testPoolingLayer)

ROOT_EXECUTABLE(testReshape TestReshape.cxx LIBRARIES ${Libraries})
ROOT_ADD_TEST(TMVA-DNN-CNN-Reshape COMMAND testReshape)

ROOT_EXECUTABLE(testRotWeights TestRotateWeights.cxx LIBRARIES ${Libraries})
ROOT_ADD_TEST(TMVA-DNN-CNN-RotWeights COMMAND testRotWeights)

#-- features not implemeted in ref architectures
#ROOT_EXECUTABLE(testForwardPass TestForwardPass.cxx LIBRARIES ${Libraries})
#ROOT_ADD_TEST(TMVA-DNN-CNN-Forward COMMAND testForwardPass)

#ROOT_EXECUTABLE(testConvNetLoss TestConvNetLoss.cxx LIBRARIES ${Libraries})
#ROOT_ADD_TEST(TMVA-DNN-CNN-Loss COMMAND testConvNetLoss)

#ROOT_EXECUTABLE(testConvNetPred TestConvNetPrediction.cxx LIBRARIES ${Libraries})
#ROOT_ADD_TEST(TMVA-DNN-CNN-Pred COMMAND testConvNetPred)

#ROOT_EXECUTABLE(testDLMinimization TestMinimization.cxx LIBRARIES ${Libraries})
#ROOT_ADD_TEST(TMVA-DNN-CNN-Minimization COMMAND testDLMinimization)

#ROOT_EXECUTABLE(testTensorDataLoader TestTensorDataLoader.cxx LIBRARIES ${Libraries})
#ROOT_ADD_TEST(TMVA-DNN-Tensor-Data-Loader COMMAND testTensorDataLoader)

endif()

#--- CPU tests. ----------------------------
if ((BLAS_FOUND OR mathmore) AND imt AND tmva-cpu)

include_directories(SYSTEM ${TBB_INCLUDE_DIRS})

ROOT_EXECUTABLE(testIm2ColCpu TestIm2ColCpu.cxx LIBRARIES ${Libraries})
ROOT_ADD_TEST(TMVA-DNN-CNN-Im2Col-CPU COMMAND testIm2ColCpu)

ROOT_EXECUTABLE(testPoolingLayerCpu TestPoolingLayerCpu.cxx LIBRARIES ${Libraries})
ROOT_ADD_TEST(TMVA-DNN-CNN-PoolingLayer-CPU COMMAND testPoolingLayerCpu)

ROOT_EXECUTABLE(testConvLayerCpu TestConvLayerCpu.cxx LIBRARIES ${Libraries})
ROOT_ADD_TEST(TMVA-DNN-CNN-ConvLayer-CPU COMMAND testConvLayerCpu)

ROOT_EXECUTABLE(testRotWeightsCpu TestRotateWeightsCpu.cxx LIBRARIES ${Libraries})
ROOT_ADD_TEST(TMVA-DNN-CNN-RotWeights-CPU COMMAND testRotWeightsCpu)

ROOT_EXECUTABLE(testForwardPassCpu TestForwardPassCpu.cxx LIBRARIES ${Libraries})
ROOT_ADD_TEST(TMVA-DNN-CNN-Forward-CPU COMMAND testForwardPassCpu)

ROOT_EXECUTABLE(testConvNetLossCpu TestConvNetLossCpu.cxx LIBRARIES ${Libraries})
ROOT_ADD_TEST(TMVA-DNN-CNN-Loss-CPU COMMAND testConvNetLossCpu)

ROOT_EXECUTABLE(testConvNetPredCpu TestConvNetPredictionCpu.cxx LIBRARIES ${Libraries})
ROOT_ADD_TEST(TMVA-DNN-CNN-Pred-CPU COMMAND testConvNetPredCpu)

ROOT_EXECUTABLE(testReshapeCpu TestReshapeCpu.cxx LIBRARIES ${Libraries})
ROOT_ADD_TEST(TMVA-DNN-CNN-Reshape-CPU COMMAND testReshapeCpu)

#-- need to be fixed
#ROOT_EXECUTABLE(testTensorDataLoaderCpu TestTensorDataLoaderCpu.cxx LIBRARIES ${Libraries})
#ROOT_ADD_TEST(TMVA-DNN-Tensor-Data-Loader-CPU COMMAND testTensorDataLoaderCpu)

#ROOT_EXECUTABLE(testDLMinimizationCpu TestMinimizationCpu.cxx LIBRARIES ${Libraries})
#ROOT_ADD_TEST(TMVA-DNN-CNN-Minimization-CPU COMMAND testDLMinimizationCpu)

ROOT_EXECUTABLE(testConvBackpropagationCpu TestConvBackpropagation.cxx LIBRARIES ${Libraries})
ROOT_ADD_TEST(TMVA-DNN-CNN-Backpropagation-CPU COMMAND testConvBackpropagationCpu)

ROOT_EXECUTABLE(testMethodDLCNN TestMethodDLCNN.cxx LIBRARIES ${Libraries})
ROOT_ADD_TEST(TMVA-DNN-CNN-MethodDL COMMAND testMethodDLCNN)

endif ()
