#find_package(ZeroMQ REQUIRED)

ROOT_ADD_GTEST(test_RooFitZMQ test_ZMQ.cpp LIBRARIES
        RooFitZMQ ${ZeroMQ_LIBRARY}
        RooFitMultiProcess # wait_for_child in util
)

ROOT_ADD_GTEST(test_RooFitZMQ_polling test_polling.cxx LIBRARIES RooFitZMQ ${ZeroMQ_LIBRARY})
ROOT_ADD_GTEST(test_RooFitZMQ_HWM test_HWM.cxx LIBRARIES RooFitZMQ ${ZeroMQ_LIBRARY})
ROOT_ADD_GTEST(test_RooFitZMQ_load_balancing test_ZMQ_load_balancing.cxx LIBRARIES RooFitZMQ ${ZeroMQ_LIBRARY})
