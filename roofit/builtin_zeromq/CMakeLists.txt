include(ExternalProject)

set(ZeroMQ_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/BUILTIN_ZeroMQ-prefix)
set(ZeroMQ_LIBNAME ${CMAKE_STATIC_LIBRARY_PREFIX}zmq${CMAKE_STATIC_LIBRARY_SUFFIX})

set(ZeroMQ_LIBRARY ${ZeroMQ_PREFIX}/lib/${ZeroMQ_LIBNAME} CACHE INTERNAL "" FORCE)
set(ZeroMQ_LIBRARIES ${ZeroMQ_LIBRARY} CACHE INTERNAL "" FORCE)

ExternalProject_Add(BUILTIN_ZeroMQ
    URL https://github.com/zeromq/libzmq/archive/refs/tags/v4.3.2.tar.gz
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        -DWITH_PERF_TOOL=OFF
        -DZMQ_BUILD_TESTS=OFF
        -DPOLLER=select
    BUILD_BYPRODUCTS ${ZeroMQ_LIBRARIES}
    )

set(ZeroMQ_FOUND TRUE CACHE BOOL "" FORCE)

set(ZeroMQ_INCLUDE_DIR  ${ZeroMQ_PREFIX}/include CACHE INTERNAL "" FORCE)
set(ZeroMQ_INCLUDE_DIRS ${ZeroMQ_PREFIX}/include CACHE INTERNAL "" FORCE)

add_library(ZeroMQ INTERFACE)
target_include_directories(ZeroMQ INTERFACE $<BUILD_INTERFACE:${ZeroMQ_INCLUDE_DIR}>)
target_link_libraries(ZeroMQ INTERFACE $<BUILD_INTERFACE:${ZeroMQ_LIBRARIES}>)
add_dependencies(ZeroMQ BUILTIN_ZeroMQ)

add_library(ZeroMQ::ZeroMQ ALIAS ZeroMQ)

set(ZeroMQ_DIR ${ZeroMQ_PREFIX}/share/cmake/ZeroMQ CACHE INTERNAL "" FORCE)

# also need source and build directories for ppoll in RooFitZMQ
set(ZeroMQ_SOURCE_DIR ${ZeroMQ_PREFIX}/src/BUILTIN_ZeroMQ/src CACHE INTERNAL "" FORCE)
set(ZeroMQ_BUILD_DIR ${ZeroMQ_PREFIX}/src/BUILTIN_ZeroMQ-build CACHE INTERNAL "" FORCE)