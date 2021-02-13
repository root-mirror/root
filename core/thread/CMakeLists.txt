# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT core/thread package
############################################################################

if(WIN32)
  set(PLATFORM_FILTER FILTER "Posix")
  set(PLATFORM_HEADERS
    TWin32Condition.h
    TWin32Mutex.h
    TWin32Thread.h
    TWin32ThreadFactory.h
  )
else()
  set(PLATFORM_FILTER FILTER "Win32")
  set(PLATFORM_HEADERS
    TPosixCondition.h
    TPosixMutex.h
    TPosixThread.h
    TPosixThreadFactory.h
    PosixThreadInc.h
  )
endif()

ROOT_STANDARD_LIBRARY_PACKAGE(Thread
  HEADERS
    ${PLATFORM_HEADERS}
    TAtomicCount.h
    TCondition.h
    TConditionImp.h
    ThreadLocalStorage.h
    TMutex.h
    TMutexImp.h
    TRWLock.h
    TSemaphore.h
    TThreadFactory.h
    TThread.h
    TThreadImp.h
    TThreadPool.h
    ROOT/RConcurrentHashColl.hxx
    ROOT/TReentrantRWLock.hxx
    ROOT/TRWSpinLock.hxx
    ROOT/TSpinMutex.hxx
    ROOT/TThreadedObject.hxx
  SOURCES
    src/RConcurrentHashColl.cxx
    src/TCondition.cxx
    src/TConditionImp.cxx
    src/TMutex.cxx
    src/TMutexImp.cxx
    src/TReentrantRWLock.cxx
    src/TRWLock.cxx
    src/TRWMutexImp.cxx
    src/TRWSpinLock.cxx
    src/TSemaphore.cxx
    src/TThread.cxx
    src/TThreadFactory.cxx
    src/TThreadImp.cxx
  STAGE1
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM 
  DEPENDENCIES
    Core
  INSTALL_OPTIONS ${installoptions}
)

target_link_libraries(Thread PUBLIC ${CMAKE_THREAD_LIBS_INIT})

target_include_directories(Thread PRIVATE
   ${CMAKE_SOURCE_DIR}/core/base/inc
   ${CMAKE_SOURCE_DIR}/core/clib/inc
   ${CMAKE_BINARY_DIR}/ginclude
)

# keep include directory for ROOT/RSha256.hxx private
set_source_files_properties(src/RConcurrentHashColl.cxx
  PROPERTIES COMPILE_FLAGS -I${CMAKE_SOURCE_DIR}/core/foundation/res)

if(WIN32)
  target_sources(Thread PRIVATE
    src/TWin32Condition.cxx
    src/TWin32Mutex.cxx
    src/TWin32Thread.cxx
    src/TWin32ThreadFactory.cxx
  )
  target_include_directories(Thread PRIVATE ${CMAKE_SOURCE_DIR}/core/winnt/inc)
else()
  target_sources(Thread PRIVATE
    src/TPosixCondition.cxx
    src/TPosixMutex.cxx
    src/TPosixThread.cxx
    src/TPosixThreadFactory.cxx
  )
endif()

ROOT_ADD_TEST_SUBDIRECTORY(test)
