############################################################################
# CMakeLists.txt file for building ROOT roofitcore/ZMQ package
# @author Patrick Bos, NL eScience Center
############################################################################

#if(NOT ZeroMQ_FOUND)
#    find_package(ZeroMQ REQUIRED)
#endif()

ROOT_STANDARD_LIBRARY_PACKAGE(RooFitZMQ
    HEADERS
        RooFit_ZMQ/zmq.hxx
        RooFit_ZMQ/Utility.h
        RooFit_ZMQ/ZeroMQSvc.h
        RooFit_ZMQ/functions.h
        RooFit_ZMQ/ZeroMQPoller.h
        RooFit_ZMQ/ppoll.h
    SOURCES
        src/ZeroMQSvc.cpp
        src/ZeroMQPoller.cpp
        src/functions.cpp
        src/ppoll.cpp
    DICTIONARY_OPTIONS
#        -I${ZeroMQ_INCLUDE_DIR}
        -writeEmptyRootPCM
    DEPENDENCIES
        Core # TVersionCheck, RegisterModule, these are always called from G__RooFitZMQ.cxx.o
#        RIO  # TBufferFile, TBufferIO
    LIBRARIES
#        ${ZeroMQ_LIBRARY}
        ZeroMQ
    )


target_link_libraries(RooFitZMQ PRIVATE ZeroMQ)
target_include_directories(RooFitZMQ PUBLIC $<BUILD_INTERFACE:${ZeroMQ_INCLUDE_DIR}>)
target_include_directories(RooFitZMQ PUBLIC $<BUILD_INTERFACE:${ZeroMQ_SOURCE_DIR}>)
target_include_directories(RooFitZMQ PUBLIC $<BUILD_INTERFACE:${ZeroMQ_BUILD_DIR}>)

if(testing)
    ROOT_ADD_TEST_SUBDIRECTORY(test)
endif()
