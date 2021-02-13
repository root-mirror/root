# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT (global) core package
############################################################################

file(WRITE ${CMAKE_BINARY_DIR}/etc/gitinfo.txt
  "${GIT_DESCRIBE_ALL}\n${GIT_DESCRIBE_ALWAYS}\n${GIT_TIMESTAMP}\n")

file(WRITE ${CMAKE_BINARY_DIR}/RGitCommit.h.tmp
"#ifndef ROOT_RGITCOMMIT_H
#define ROOT_RGITCOMMIT_H
  #define ROOT_GIT_BRANCH \"${GIT_DESCRIBE_ALL}\"
  #define ROOT_GIT_COMMIT \"${GIT_DESCRIBE_ALWAYS}\"
#endif"
)

set_source_files_properties(${CMAKE_BINARY_DIR}/ginclude/RConfigure.h
  PROPERTIES GENERATED TRUE)

add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/include/RConfigure.h
  COMMAND
    ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_BINARY_DIR}/ginclude/RConfigure.h
    ${CMAKE_BINARY_DIR}/include/RConfigure.h
  DEPENDS
    ${CMAKE_BINARY_DIR}/ginclude/RConfigure.h
)

add_custom_target(rconfigure ALL DEPENDS ${CMAKE_BINARY_DIR}/include/RConfigure.h)

add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/ginclude/RGitCommit.h
  COMMAND
    ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_BINARY_DIR}/RGitCommit.h.tmp
    ${CMAKE_BINARY_DIR}/ginclude/RGitCommit.h
  DEPENDS
    ${CMAKE_BINARY_DIR}/RGitCommit.h.tmp
)

add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/include/RGitCommit.h
  COMMAND
    ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_BINARY_DIR}/RGitCommit.h.tmp
    ${CMAKE_BINARY_DIR}/include/RGitCommit.h
  DEPENDS
    ${CMAKE_BINARY_DIR}/RGitCommit.h.tmp
)

add_custom_target(gitcommit ALL DEPENDS ${CMAKE_BINARY_DIR}/include/RGitCommit.h)

set_source_files_properties(${CMAKE_BINARY_DIR}/ginclude/RGitCommit.h
  PROPERTIES GENERATED TRUE HEADER_FILE_ONLY TRUE)

set_source_files_properties(${CMAKE_BINARY_DIR}/include/RGitCommit.h
  PROPERTIES GENERATED TRUE HEADER_FILE_ONLY TRUE)

install(FILES ${CMAKE_BINARY_DIR}/ginclude/RGitCommit.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_subdirectory(clib)
add_subdirectory(clingutils)
add_subdirectory(cont)
add_subdirectory(dictgen)
add_subdirectory(foundation)
add_subdirectory(gui)
add_subdirectory(meta)
add_subdirectory(metacling)
if(NOT WIN32)
add_subdirectory(multiproc)
endif()
add_subdirectory(rint)
add_subdirectory(textinput)
add_subdirectory(thread)
add_subdirectory(imt)
add_subdirectory(zip)
add_subdirectory(lzma)
add_subdirectory(lz4)
add_subdirectory(zstd)
if(asan)
  add_subdirectory(sanitizer)
endif()

if(NOT WIN32)
  add_subdirectory(newdelete)
endif()

add_subdirectory(base)

if(UNIX)
  add_subdirectory(unix)
  set(unix_objects $<TARGET_OBJECTS:Unix>)
endif()
if(WIN32)
  add_subdirectory(winnt)
  set(winnt_objects $<TARGET_OBJECTS:Winnt>)
endif()
if(cocoa)
  add_subdirectory(macosx)
  set(macosx_objects $<TARGET_OBJECTS:Macosx>)
endif()

if(CMAKE_CXX_STANDARD GREATER 11)
  set(dict_v7dirs base/v7/inc)
endif()

#---G__Core--------------------------------------------------------------------

# Uses includes path as defined by Core module.

set(objectlibs $<TARGET_OBJECTS:Base>
               $<TARGET_OBJECTS:Clib>
               $<TARGET_OBJECTS:Cont>
               $<TARGET_OBJECTS:Foundation>
               $<TARGET_OBJECTS:Lzma>
               $<TARGET_OBJECTS:Lz4>
               $<TARGET_OBJECTS:Zstd>
               $<TARGET_OBJECTS:Zip>
               $<TARGET_OBJECTS:Meta>
               $<TARGET_OBJECTS:GuiCore>
               $<TARGET_OBJECTS:TextInput>
               ${macosx_objects}
               ${unix_objects}
               ${winnt_objects})

ROOT_OBJECT_LIBRARY(BaseTROOT ${CMAKE_SOURCE_DIR}/core/base/src/TROOT.cxx ${CMAKE_BINARY_DIR}/ginclude/RGitCommit.h)
target_include_directories(BaseTROOT PRIVATE
   ${CMAKE_SOURCE_DIR}/core/foundation/res
   ${CMAKE_SOURCE_DIR}/core/base/inc
   ${CMAKE_SOURCE_DIR}/core/clib/inc
   ${CMAKE_SOURCE_DIR}/core/meta/inc
   ${CMAKE_SOURCE_DIR}/core/cont/inc
   ${CMAKE_SOURCE_DIR}/core/foundation/inc
   ${CMAKE_SOURCE_DIR}/core/thread/inc
   ${CMAKE_SOURCE_DIR}/core/gui/inc
   ${CMAKE_BINARY_DIR}/ginclude
)

add_dependencies(BaseTROOT gitcommit rconfigure)
#----------------------------------------------------------------------------------------

if(WIN32)
  set(corelinklibs shell32.lib WSock32.lib Oleaut32.lib Iphlpapi.lib)
elseif(APPLE)
  if(cocoa)
     set(corelinklibs "-framework Cocoa -F/System/Library/PrivateFrameworks -framework CoreSymbolication")
  else()
     set(corelinklibs "-F/System/Library/PrivateFrameworks -framework CoreSymbolication")
  endif()
endif()

add_subdirectory(rootcling_stage1)

#-------------------------------------------------------------------------------
ROOT_LINKER_LIBRARY(Core $<TARGET_OBJECTS:BaseTROOT> ${objectlibs} BUILTINS LZMA)

if (libcxx AND NOT APPLE)
  # In case we use libcxx and glibc together there is a mismatch of the
  # signatures of functions in the header wchar.h. This macro tweaks the
  # header in rootcling resource directory to be compatible with the one from
  # libc++.
  target_compile_definitions(Core PRIVATE __CORRECT_ISO_CPP_WCHAR_H_PROTO)
endif()

if(MSVC)
  # Definitions of MATH Macros (required in MathCore) like M_PI are only  
  # exposed on Windows after defining the _USE_MATH_DEFINES macro. By   
  # specifying it as a property for Core, we ensure that the macros are 
  # properly exposed when using Modules on Windows.
  target_compile_definitions(Core PRIVATE _USE_MATH_DEFINES)
endif()

#while basic libs do not depend on Core, we have to add includes directly

target_include_directories(Core PUBLIC
   $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/ginclude>
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/core/base/inc>
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/core/foundation/inc>
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/core/cont/inc>
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/core/gui/inc>
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/core/meta/inc>
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/core/clib/inc>
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/core/rint/inc>
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/core/zip/inc>
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/core/thread/inc>
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/core/textinput/inc>
)

target_include_directories(Core PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/core/clingutils/inc>)

if(root7 OR CMAKE_CXX_STANDARD GREATER 11)
  target_include_directories(Core PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/core/base/v7/inc>)
  target_include_directories(Core PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/core/foundation/v7/inc>)
endif()


if(ROOT_ARCHITECTURE MATCHES macosx)
  target_include_directories(Core PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/core/unix/inc> $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/core/macosx/inc>)
  target_include_directories(BaseTROOT PRIVATE ${CMAKE_SOURCE_DIR}/core/unix/inc ${CMAKE_SOURCE_DIR}/core/macosx/inc)
elseif(ROOT_ARCHITECTURE MATCHES win32)
  target_include_directories(Core PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/core/winnt/inc>)
  target_include_directories(BaseTROOT PRIVATE ${CMAKE_SOURCE_DIR}/core/winnt/inc)
else()
  target_include_directories(Core PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/core/unix/inc>)
  target_include_directories(BaseTROOT PRIVATE ${CMAKE_SOURCE_DIR}/core/unix/inc)
endif()

if (runtime_cxxmodules)
  if(MSVC)
    list(APPEND core_implicit_modules "-mByproduct" "services")
  endif()
  # Force generation of _Builtin_intrinsics from Core.
  list(APPEND core_implicit_modules "-m" "_Builtin_intrinsics" "-mByproduct" "_Builtin_intrinsics")
  list(APPEND core_implicit_modules "-mByproduct" "ROOT_Foundation_Stage1_NoRTTI")
  list(APPEND core_implicit_modules "-mByproduct" "ROOT_Foundation_C")
  list(APPEND core_implicit_modules "-mByproduct" "ROOT_Rtypes")
endif(runtime_cxxmodules)

ROOT_GENERATE_DICTIONARY(G__Core
  ${Core_dict_headers}
  ${Clib_dict_headers}
  ${Cont_dict_headers}
  ${Foundation_dict_headers}
  ${Macosx_dict_headers}
  ${Unix_dict_headers}
  ${Winnt_dict_headers}
  ${ClingUtils_dict_headers}
  ${GuiCore_dict_headers}
  ${Meta_dict_headers}
  ${TextInput_dict_headers}
  STAGE1
  MODULE
    Core
  OPTIONS
    -writeEmptyRootPCM
    ${core_implicit_modules}
  LINKDEF
    base/inc/LinkDef.h
)

target_include_directories(G__Core PRIVATE ${CMAKE_SOURCE_DIR}/core/clingutils/inc)

target_link_libraries(Core
  PRIVATE
    PCRE::PCRE
    ${LIBLZMA_LIBRARIES}
    xxHash::xxHash
    LZ4::LZ4
    ZLIB::ZLIB
    ${ZSTD_LIBRARIES}
    ${CMAKE_DL_LIBS}
    ${CMAKE_THREAD_LIBS_INIT}
    ${corelinklibs}
)

add_dependencies(Core CLING)
