# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT core/clingutils package
############################################################################

# These files depend on cling/clang/llvm; they need to be linked into libCling.
# They are used by rootcling_stage1, rootcling and libCling.

set(ClingUtils_dict_headers
  root_std_complex.h
  PARENT_SCOPE
)

ROOT_OBJECT_LIBRARY(ClingUtils
  src/RStl.cxx
  src/TClingUtils.cxx
)

add_dependencies(ClingUtils CLING)

target_include_directories(ClingUtils PRIVATE 
   ${CLING_INCLUDE_DIRS}
   ${CMAKE_SOURCE_DIR}/core/foundation/res
   ${CMAKE_SOURCE_DIR}/core/foundation/inc
   ${CMAKE_SOURCE_DIR}/core/base/inc
   ${CMAKE_SOURCE_DIR}/core/clib/inc
   ${CMAKE_SOURCE_DIR}/core/meta/inc
   ${CMAKE_BINARY_DIR}/ginclude)

# Register the llvm include directories after clangs. This instructs the compiler to resolve
# headers from our builtin clang. That's an issue when we are building with bultin_llvm=Off
# and we have installed clang headers, too.
target_include_directories(ClingUtils SYSTEM PRIVATE ${CLANG_INCLUDE_DIRS} ${LLVM_INCLUDE_DIRS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CLING_CXXFLAGS}")

# This is to avoid warnings coming from GCC 7 in llvm/src/include/llvm/ADT/DenseMap.h:1010
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  ROOT_ADD_CXX_FLAG(CMAKE_CXX_FLAGS -Wno-maybe-uninitialized)
endif()

ROOT_INSTALL_HEADERS()

#### STL dictionary (replacement for cintdlls)##############################

set(stldicts
    vector
    list
    forward_list
    deque
    map map2 unordered_map
    multimap multimap2 unordered_multimap
    set unordered_set
    multiset unordered_multiset
    complex)
if(NOT WIN32)
  list(APPEND stldicts valarray)
endif()
foreach(dict ${stldicts})
  string(REPLACE "2" "" header ${dict})
  string(REPLACE "complex" "root_std_complex.h" header ${header})
  string(REPLACE "multi" "" header ${header})
  ROOT_STANDARD_LIBRARY_PACKAGE(${dict}Dict
                                NO_SOURCES NO_INSTALL_HEADERS NO_CXXMODULE
                                STAGE1
                                NODEPHEADERS ${header}
                                LINKDEF src/${dict}Linkdef.h
                                DEPENDENCIES Core)
  target_include_directories(${dict}Dict PRIVATE ${CMAKE_SOURCE_DIR}/interpreter/cling/include/cling/cint)
endforeach()

set(CLANG_RESOURCE_DIR_STEM)
if (builtin_clang)
  set(CLANG_RESOURCE_DIR_STEM ${CMAKE_BINARY_DIR}/interpreter/llvm/src/${CMAKE_CFG_INTDIR}/lib/clang)
  set(CLANG_RESOURCE_DIR_VERSION ${LLVM_VERSION})
else ()
  set(CLANG_RESOURCE_DIR_STEM ${LLVM_LIBRARY_DIR}/clang)
  # A user can define a clang version to use, otherwise find it (but will error if more than one version is present)
  if (NOT DEFINED CLANG_RESOURCE_DIR_VERSION)
    if (NOT EXISTS ${CLANG_RESOURCE_DIR_STEM})
      message(FATAL_ERROR "${CLANG_RESOURCE_DIR_STEM} does not exist. Please install clang.")
    endif()
    # There is no reasonable way to get the version of clang under which is its resource directory.
    # For example, lib/clang/5.0.0/include. Deduce it.
    file(GLOB CHILDREN RELATIVE ${CLANG_RESOURCE_DIR_STEM} ${CLANG_RESOURCE_DIR_STEM}/*)
    list(LENGTH CHILDREN CHILDREN_LENGTH)
    if (${CHILDREN_LENGTH} GREATER 1)
      message(FATAL_ERROR "Found more than one version of clang. CLANG_RESOURCE_DIR_VERSION contains: '${CHILDREN}'." )
    endif()

    list(GET CHILDREN 0 CLANG_RESOURCE_DIR_VERSION)
  endif()
endif()


set(CLANG_RESOURCE_DIR ${CLANG_RESOURCE_DIR_STEM}/${CLANG_RESOURCE_DIR_VERSION}/include)

#---Deal with clang resource here----------------------------------------------
install(DIRECTORY ${CMAKE_BINARY_DIR}/etc/cling/lib/clang/${CLANG_RESOURCE_DIR_VERSION}/include/
        DESTINATION ${CMAKE_INSTALL_SYSCONFDIR}/cling/lib/clang/${CLANG_RESOURCE_DIR_VERSION}/include USE_SOURCE_PERMISSIONS)


#---Install a bunch of files to /etc/cling------------------------------------
set(clinginclude ${CMAKE_SOURCE_DIR}/interpreter/cling/include)

set(custom_modulemaps)
if (runtime_cxxmodules)
  set(custom_modulemaps boost.modulemap tinyxml2.modulemap cuda.modulemap module.modulemap.build)
  if (NOT libcxx)
    if (MSVC)
      set(custom_modulemaps ${custom_modulemaps} vcruntime.modulemap)
      set(custom_modulemaps ${custom_modulemaps} services_msvc.modulemap)
      set(custom_modulemaps ${custom_modulemaps} std_msvc.modulemap)
    else()
      set(custom_modulemaps ${custom_modulemaps} std.modulemap)
    endif()
  endif()
  # Handle libc. Apple's libc is modularized.
  if (MSVC)
    set(custom_modulemaps ${custom_modulemaps} libc_msvc.modulemap)
  elseif (NOT APPLE)
    set(custom_modulemaps ${custom_modulemaps} libc.modulemap)
  endif()
endif(runtime_cxxmodules)

foreach(file ${custom_modulemaps}
        Interpreter/DynamicExprInfo.h
        Interpreter/DynamicLookupRuntimeUniverse.h
        Interpreter/DynamicLookupLifetimeHandler.h
        Interpreter/Exception.h
        Interpreter/RuntimePrintValue.h
        Interpreter/RuntimeUniverse.h
        Interpreter/RuntimeOptions.h
        Interpreter/Value.h)
  get_filename_component(path ${file} PATH)
  set(dest_file ${file})
  if (${file} STREQUAL "module.modulemap.build")
    set(dest_file "module.modulemap")
  elseif(NOT ${file} MATCHES ".*modulemap")
    # We do not want our modulemap to be considered part of the PCH.
    set_property(GLOBAL APPEND PROPERTY CLINGETCPCH etc/cling/${dest_file})
  endif()
  list(APPEND copy_commands COMMAND ${CMAKE_COMMAND} -E copy ${clinginclude}/cling/${file} ${CMAKE_BINARY_DIR}/etc/cling/${dest_file})
  list(APPEND files_to_copy ${clinginclude}/cling/${file})
  install(FILES ${CMAKE_BINARY_DIR}/etc/cling/${dest_file} DESTINATION ${CMAKE_INSTALL_SYSCONFDIR}/cling/${path})
endforeach()

foreach(file  multimap  multiset)
  list(APPEND copy_commands COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/interpreter/cling/include/cling/cint/${file} ${CMAKE_BINARY_DIR}/etc/cling/cint/${file})
  list(APPEND files_to_copy ${CMAKE_SOURCE_DIR}/interpreter/cling/include/cling/cint/${file})
  install(FILES ${CMAKE_SOURCE_DIR}/interpreter/cling/include/cling/cint/${file} DESTINATION ${CMAKE_INSTALL_SYSCONFDIR}/cling/cint)
endforeach()

# These headers do not have complete include guards on all platforms we
# support. This means that the PCH cannot provide their representation at
# runtime and clang will hit disk, triggering a possible incompatibility
# of that file in build-time versus run-time (different length etc).
# Capture their build-time version here, and inject it into runtime.
foreach(file wchar.h bits/stat.h bits/time.h)
  if(EXISTS ${DEFAULT_SYSROOT}/usr/include/${file})
    list(APPEND copy_commands COMMAND ${CMAKE_COMMAND} -E copy ${DEFAULT_SYSROOT}/usr/include/${file}
      ${CMAKE_BINARY_DIR}/etc/cling/lib/clang/${CLANG_RESOURCE_DIR_VERSION}/include/${file})
  endif()
endforeach()

set(stamp_file ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/LLVMRES.stamp)
if(MSVC)
  add_custom_command(OUTPUT ${stamp_file}
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/etc/cling/lib/clang/${CLANG_RESOURCE_DIR_VERSION}/include
        ${copy_commands}
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CLANG_RESOURCE_DIR}
        ${CMAKE_BINARY_DIR}/etc/cling/lib/clang/${CLANG_RESOURCE_DIR_VERSION}/include		
        COMMAND ${CMAKE_COMMAND} -E touch ${stamp_file}
        DEPENDS ${files_to_copy}
        COMMENT "Copying LLVM resource and header files")
else()
  add_custom_command(OUTPUT ${stamp_file}
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/etc/cling/lib/clang/${CLANG_RESOURCE_DIR_VERSION}/include
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CLANG_RESOURCE_DIR}
        ${CMAKE_BINARY_DIR}/etc/cling/lib/clang/${CLANG_RESOURCE_DIR_VERSION}/include
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_SOURCE_DIR}/interpreter/llvm/ROOT/assert.h
        ${CMAKE_BINARY_DIR}/etc/cling/lib/clang/${CLANG_RESOURCE_DIR_VERSION}/include/assert.h
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_SOURCE_DIR}/interpreter/llvm/ROOT/stdlib.h
        ${CMAKE_BINARY_DIR}/etc/cling/lib/clang/${CLANG_RESOURCE_DIR_VERSION}/include/stdlib.h
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_SOURCE_DIR}/interpreter/llvm/ROOT/unistd.h
        ${CMAKE_BINARY_DIR}/etc/cling/lib/clang/${CLANG_RESOURCE_DIR_VERSION}/include/unistd.h
        ${copy_commands}
        COMMAND ${CMAKE_COMMAND} -E touch ${stamp_file}
        DEPENDS ${files_to_copy}
        COMMENT "Copying LLVM resource and header files")
endif()
add_custom_target(LLVMRES DEPENDS ${stamp_file} CLING)
# CLING is a shorthand for CLING_LIBRARIES and some other clang-specific
# dependencies which ensure the correct order of building. Then the cling header
# files (such as RuntimeUniverse.h) are moved to a semi-private place in ROOT
# #ROOTSYS/etc. This is the place where ROOT will use them from and we should
# add an explcit dependency to something cling-related which ROOT knows.
# ClingUtils seems a good candidate because it is very foundational.
add_dependencies(ClingUtils LLVMRES)
ROOT_ADD_TEST_SUBDIRECTORY(test)

