include(ExternalProject)

set(ZeroMQ_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/BUILTIN_ZeroMQ-prefix)
set(ZeroMQ_LIBNAME ${CMAKE_STATIC_LIBRARY_PREFIX}zmq${CMAKE_STATIC_LIBRARY_SUFFIX})

set(ZeroMQ_LIBRARY ${ZeroMQ_PREFIX}/lib/${ZeroMQ_LIBNAME} CACHE INTERNAL "" FORCE)
set(ZeroMQ_LIBRARIES ${ZeroMQ_LIBRARY} CACHE INTERNAL "" FORCE)

ExternalProject_Add(BUILTIN_ZeroMQ
    URL https://github.com/zeromq/libzmq/archive/refs/tags/v4.3.2.tar.gz
    URL_HASH SHA256=02ecc88466ae38cf2c8d79f09cfd2675ba299a439680b64ade733e26a349edeb
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        -DWITH_PERF_TOOL=OFF
        -DZMQ_BUILD_TESTS=OFF
        -DPOLLER=select
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        # TODO: remove -w flag when updating to v4.3.5 or higher, current warnings should be gone by then
        -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}\ -w
    BUILD_BYPRODUCTS ${ZeroMQ_LIBRARIES}
    )

ExternalProject_Get_Property(BUILTIN_ZeroMQ install_dir)

set(ZeroMQ_FOUND TRUE CACHE BOOL "" FORCE)

set(ZeroMQ_INCLUDE_DIR  ${ZeroMQ_PREFIX}/include CACHE INTERNAL "" FORCE)
set(ZeroMQ_INCLUDE_DIRS ${ZeroMQ_PREFIX}/include CACHE INTERNAL "" FORCE)

add_library(ZeroMQ INTERFACE)
target_include_directories(ZeroMQ INTERFACE $<BUILD_INTERFACE:${ZeroMQ_INCLUDE_DIR}>)
target_link_libraries(ZeroMQ INTERFACE $<BUILD_INTERFACE:${ZeroMQ_LIBRARIES}>)
add_dependencies(ZeroMQ BUILTIN_ZeroMQ)

add_library(ZeroMQ::ZeroMQ ALIAS ZeroMQ)
add_library(libzmq ALIAS ZeroMQ)

set(ZeroMQ_DIR ${ZeroMQ_PREFIX}/share/cmake/ZeroMQ CACHE INTERNAL "" FORCE)

# also need source and build directories for ppoll in RooFitZMQ
set(ZeroMQ_SOURCE_DIR ${ZeroMQ_PREFIX}/src/BUILTIN_ZeroMQ/src CACHE INTERNAL "" FORCE)
set(ZeroMQ_BUILD_DIR ${ZeroMQ_PREFIX}/src/BUILTIN_ZeroMQ-build CACHE INTERNAL "" FORCE)

ExternalProject_Add(BUILTIN_cppzmq
    URL https://github.com/zeromq/cppzmq/archive/refs/tags/v4.3.0.tar.gz
    URL_HASH SHA256=27d1f56406ba94ee779e639203218820975cf68174f92fbeae0f645df0fcada4
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${install_dir}
        -DCPPZMQ_BUILD_TESTS=OFF
    DEPENDS BUILTIN_ZeroMQ
    )
