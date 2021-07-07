############################################################################
# CMakeLists.txt file for building ROOT RooFitMultiProcess package
# @author Patrick Bos, Netherlands eScience Center
############################################################################

#FIND_PACKAGE(ZeroMQ REQUIRED)

ROOT_STANDARD_LIBRARY_PACKAGE(RooFitMultiProcess
    HEADERS
        RooFit/MultiProcess/worker.h
        RooFit/MultiProcess/Messenger_decl.h
        RooFit/MultiProcess/Messenger.h
        RooFit/MultiProcess/ProcessManager.h
        RooFit/MultiProcess/JobManager.h
        RooFit/MultiProcess/util.h
        RooFit/MultiProcess/Queue.h
        RooFit/MultiProcess/types.h
        RooFit/MultiProcess/Job.h
    SOURCES
        src/worker.cxx
        src/Messenger.cxx
        src/ProcessManager.cxx
        src/util.cxx
        src/Queue.cxx
        src/JobManager.cxx
        src/Job.cxx
    DICTIONARY_OPTIONS
        "-writeEmptyRootPCM"
    DEPENDENCIES
        RooFitZMQ
    LIBRARIES
        ZeroMQ
)

target_link_libraries(RooFitMultiProcess PRIVATE ZeroMQ)
target_include_directories(RooFitMultiProcess PUBLIC $<BUILD_INTERFACE:${ZeroMQ_INCLUDE_DIR}>)

if(testing)
    ROOT_ADD_TEST_SUBDIRECTORY(test)
endif()
