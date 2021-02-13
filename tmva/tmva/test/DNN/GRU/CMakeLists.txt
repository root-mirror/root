############################################################################
# CMakeLists.txt file for building TMVA/DNN/GRU tests.
# @author Surya S Dwivedi
############################################################################


set(Libraries TMVA)


#--- CUDA tests. ---------------------------
if (CUDA_FOUND AND tmva-gpu AND tmva-cudnn)

   set(DNN_CUDA_LIBRARIES ${CUDA_CUBLAS_LIBRARIES} ${CUDNN_LIBRARIES} )

   CUDA_ADD_EXECUTABLE(testGRUBackpropagationCudnn TestGRUBackpropagationCudnn.cxx)
   target_link_libraries(testGRUBackpropagationCudnn ${Libraries} ${DNN_CUDA_LIBRARIES})
   ROOT_ADD_TEST(TMVA-DNN-GRU-BackpropagationCudnn COMMAND testGRUBackpropagationCudnn)

   CUDA_ADD_EXECUTABLE(testGRUForwardPassCudnn TestGRUForwardPassCudnn.cxx)
   target_link_libraries(testGRUForwardPassCudnn ${Libraries} ${DNN_CUDA_LIBRARIES})
   ROOT_ADD_TEST(TMVA-DNN-GRU-Forwaed-Cudnn COMMAND testGRUForwardPassCudnn)

   # GRU - Full Test GPU
   CUDA_ADD_EXECUTABLE(testFullGRUCudnn TestFullGRUCudnn.cxx)
   target_link_libraries(testFullGRUCudnn ${Libraries} ${DNN_CUDA_LIBRARIES})
   ROOT_ADD_TEST(TMVA-DNN-GRU-Full-Cudnn COMMAND testFullGRUCudnn)

endif (CUDA_FOUND AND tmva-gpu AND tmva-cudnn)

#--- CPU tests. ----------------------------
if (BLAS_FOUND AND imt)
    include_directories(SYSTEM ${TBB_INCLUDE_DIRS})

    # GRU - Forward CPU
    ROOT_EXECUTABLE(testGRUForwardPassCpu TestGRUForwardPassCpu.cxx LIBRARIES ${Libraries})
    ROOT_ADD_TEST(TMVA-DNN-GRU-Forward-Cpu COMMAND testGRUForwardPassCpu)

    # GRU - Backward CPU
    ROOT_EXECUTABLE(testGRUBackpropagationCpu TestGRUBackpropagationCpu.cxx LIBRARIES ${Libraries})
    ROOT_ADD_TEST(TMVA-DNN-GRU-Backward-Cpu COMMAND testGRUBackpropagationCpu)

    # GRU - Full Test CPU
    ROOT_EXECUTABLE(testFullGRUCpu TestFullGRUCpu.cxx LIBRARIES ${Libraries})
    ROOT_ADD_TEST(TMVA-DNN-GRU-Full-Cpu COMMAND testFullGRUCpu)

endif (BLAS_FOUND AND imt)
