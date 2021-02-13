############################################################################
# CMakeLists.txt file for building TMVA/DNN/LSTM tests.
# @author Surya S Dwivedi
############################################################################


set(Libraries TMVA)

# LSTM - Forward Reference
#ROOT_EXECUTABLE(testLSTMForwardPass TestLSTMForwardPass.cxx LIBRARIES ${Libraries})
#ROOT_ADD_TEST(TMVA-DNN-LSTM-Forward COMMAND testLSTMForwardPass)

# LSTM - Backpropagation Reference
#ROOT_EXECUTABLE(testLSTMBackpropagation TestLSTMBackpropagation.cxx LIBRARIES ${Libraries})
#ROOT_ADD_TEST(TMVA-DNN-LSTM-Backpropagation COMMAND testLSTMBackpropagation)

# LSTM - Full Test Reference
#ROOT_EXECUTABLE(testFullLSTM TestFullLSTM.cxx LIBRARIES ${Libraries})
#ROOT_ADD_TEST(TMVA-DNN-LSTM-FullLSTM COMMAND testFullLSTM)


#--- CUDA tests. ---------------------------
if (CUDA_FOUND AND tmva-gpu AND tmva-cudnn)

   set(DNN_CUDA_LIBRARIES ${CUDA_CUBLAS_LIBRARIES} ${CUDNN_LIBRARIES} )

   CUDA_ADD_EXECUTABLE(testLSTMForwardPassCudnn TestLSTMForwardPassCudnn.cxx)
   target_link_libraries(testLSTMForwardPassCudnn ${Libraries} ${DNN_CUDA_LIBRARIES})
   ROOT_ADD_TEST(TMVA-DNN-LSTM-Forwaed-Cudnn COMMAND testLSTMForwardPassCudnn)

   CUDA_ADD_EXECUTABLE(testLSTMBackpropagationCudnn TestLSTMBackpropagationCudnn.cxx)
   target_link_libraries(testLSTMBackpropagationCudnn ${Libraries} ${DNN_CUDA_LIBRARIES})
   ROOT_ADD_TEST(TMVA-DNN-LSTM-BackpropagationCudnn COMMAND testLSTMBackpropagationCudnn)

   # LSTM - Full Test GPU
   CUDA_ADD_EXECUTABLE(testFullLSTMCudnn TestFullLSTMCudnn.cxx)
   target_link_libraries(testFullLSTMCudnn ${Libraries} ${DNN_CUDA_LIBRARIES})
   ROOT_ADD_TEST(TMVA-DNN-LSTM-Full-Cudnn COMMAND testFullLSTMCudnn)

endif (CUDA_FOUND AND tmva-gpu AND tmva-cudnn)

#--- CPU tests. ----------------------------
if (BLAS_FOUND AND imt)
    include_directories(SYSTEM ${TBB_INCLUDE_DIRS})

    # LSTM - Forward CPU
    ROOT_EXECUTABLE(testLSTMForwardPassCpu TestLSTMForwardPassCpu.cxx LIBRARIES ${Libraries})
    ROOT_ADD_TEST(TMVA-DNN-LSTM-Forward-Cpu COMMAND testLSTMForwardPassCpu)

    # LSTM - Backward CPU
    ROOT_EXECUTABLE(testLSTMBackpropagationCpu TestLSTMBackpropagationCpu.cxx LIBRARIES ${Libraries})
    ROOT_ADD_TEST(TMVA-DNN-LSTM-Backward-Cpu COMMAND testLSTMBackpropagationCpu)

    # LSTM - Full Test Reference
    ROOT_EXECUTABLE(testFullLSTMCpu TestFullLSTMCpu.cxx LIBRARIES ${Libraries})
    ROOT_ADD_TEST(TMVA-DNN-LSTM-Full-Cpu COMMAND testFullLSTMCpu)

endif (BLAS_FOUND AND imt)
