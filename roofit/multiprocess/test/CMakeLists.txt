# @author Patrick Bos, NL eScience Center, 2019

#find_package(ZeroMQ REQUIRED)

add_library(RooFit_multiprocess_testing_utils INTERFACE)
target_link_libraries(RooFit_multiprocess_testing_utils INTERFACE RooFitCore RooBatchCompute)

ROOT_ADD_GTEST(test_RooFitMultiProcess_Messenger test_Messenger.cxx LIBRARIES RooFitMultiProcess ${ZeroMQ_LIBRARY})
ROOT_ADD_GTEST(test_RooFitMultiProcess_ProcessManager test_ProcessManager.cxx LIBRARIES RooFitMultiProcess ${ZeroMQ_LIBRARY})
ROOT_ADD_GTEST(test_RooFitMultiProcess_Job test_Job.cxx LIBRARIES RooFitMultiProcess ${ZeroMQ_LIBRARY} RooFit_multiprocess_testing_utils)

#target_include_directories(test_RooFitMultiProcess_Job PUBLIC ${ZeroMQ_INCLUDE_DIR})
#target_include_directories(test_RooFitMultiProcess_Messenger PUBLIC ${ZeroMQ_INCLUDE_DIR})
