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
#        GradMinimizer.h
#        NLLVar.h
#        Vector.h
    SOURCES
        src/worker.cxx
        src/Messenger.cxx
        src/ProcessManager.cxx
        src/util.cxx
        src/Queue.cxx
        src/JobManager.cxx
        src/Job.cxx
#        src/GradMinimizerFcn.cxx
#        src/NLLVar.cxx
    DICTIONARY_OPTIONS
        "-writeEmptyRootPCM"
    DEPENDENCIES
        RooFitZMQ
    LIBRARIES
        ${ZeroMQ_LIBRARY}
)

if(testing)
    ROOT_ADD_TEST_SUBDIRECTORY(test)
endif()
