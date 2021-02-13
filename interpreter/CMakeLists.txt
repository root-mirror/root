#--- Check if we need to build llvm and clang ------------------------------------------------------
if (NOT builtin_clang)
  message(WARNING "Due to ROOT-specific patches you need a special version of clang. You cannot use vanilla clang.")
endif()

#---Define the way we want to build and what of llvm/clang/cling------------------------------------
set(LLVM_ENABLE_WARNINGS OFF CACHE BOOL "")
set(LLVM_INCLUDE_TESTS OFF CACHE BOOL "")
set(CLANG_INCLUDE_TESTS OFF CACHE BOOL "")
set(LLVM_INCLUDE_EXAMPLES OFF CACHE BOOL "")
set(CLANG_BUILD_TOOLS OFF CACHE BOOL "")
# It looks like that turning off CLANG_BUILD_TOOLS is not enough.
set(CLANG_TOOL_ARCMT_TEST_BUILD OFF CACHE BOOL "")
set(CLANG_TOOL_CLANG_CHECK_BUILD OFF CACHE BOOL "")
set(CLANG_TOOL_CLANG_FORMAT_BUILD OFF CACHE BOOL "")
set(CLANG_TOOL_CLANG_FORMAT_VS_BUILD OFF CACHE BOOL "")
set(CLANG_TOOL_CLANG_FUZZER_BUILD OFF CACHE BOOL "")
set(CLANG_TOOL_CLANG_IMPORT_TEST_BUILD OFF CACHE BOOL "")
set(CLANG_TOOL_CLANG_OFFLOAD_BUNDLER_BUILD OFF CACHE BOOL "")
set(CLANG_TOOL_CLANG_RENAME_BUILD OFF CACHE BOOL "")
set(CLANG_TOOL_C_ARCMT_TEST_BUILD OFF CACHE BOOL "")
set(CLANG_TOOL_C_INDEX_TEST_BUILD OFF CACHE BOOL "")
set(CLANG_TOOL_DIAGTOOL_BUILD OFF CACHE BOOL "")
set(CLANG_TOOL_LIBCLANG_BUILD OFF CACHE BOOL "")
set(CLANG_TOOL_SCAN_BUILD_BUILD OFF CACHE BOOL "")
set(CLANG_TOOL_SCAN_VIEW_BUILD OFF CACHE BOOL "")
set(LLVM_BUILD_TOOLS OFF CACHE BOOL "")
set(LLVM_TOOL_LLVM_AR_BUILD OFF CACHE BOOL "")
set(CLANG_TOOL_CLANG_OFFLOAD_BUNDLER_BUILD OFF CACHE BOOL "")
set(LLVM_FORCE_USE_OLD_TOOLCHAIN ON CACHE BOOL "")

# will be set again in case NOT builtin_llvm
set(LLVM_DIR "${CMAKE_BINARY_DIR}/interpreter/llvm/src")
if (clad)
  set(CLING_BUILD_PLUGINS ON)
endif()
set(LLVM_EXTERNAL_CLING_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cling" CACHE STRING "")
set(LLVM_EXTERNAL_PROJECTS "cling" CACHE STRING "")
# We only use llvm/clang through TCling which is (with the help of core/meta) already taking a lock
# to serialize access to llvm.  We can later review how to make this finer grained by using llvm's own locking
# mechanism.
set(LLVM_ENABLE_THREADS OFF CACHE BOOL "")

if(CMAKE_CXX_STANDARD EQUAL 17)
  set(LLVM_ENABLE_CXX1Z ON CACHE BOOL "" FORCE)
elseif(CMAKE_CXX_STANDARD EQUAL 14)
  set(LLVM_ENABLE_CXX1Y ON CACHE BOOL "" FORCE)
endif()

# The llvm::ReverseIterate<bool>::value symbol from llvm's SmallPtrSet.h
# somehow lands in our cling libraries on OS X and doesn't get hidden
# by visibility-inlines-hidden, so we suddenly have a global weak symbol
# from LLVM in cling which our visiblity=hidden compiled LLVM libraries
# reference. This is triggering some build system warnings like this:
#   ld: warning: direct access in function '(anonymous namespace)::NewGVN::runGVN()'
#   from file 'interpreter/llvm/src/lib/libLLVMScalarOpts.a(NewGVN.cpp.o)' to global weak symbol
#   'llvm::ReverseIterate<bool>::value' from file 'interpreter/llvm/src/lib/libclingUtils.a(AST.cpp.o)'
#   means the weak symbol cannot be overridden at runtime. This was likely caused by different
#   translation units being compiled with different visibility settings.
# There is no apparent reason why this is happening and it looks like a compiler bug,
# so let's just disable the part of the code that provides this symbol.
# As it's in the validation part of LLVM and not in something that providing functionality,
# this shouldn't cause any problems.
# TODO: We maybe can remove this code once we upgrade to LLVM>=6.0 as this symbol
# was introduced quite recently into LLVM 5.0 and probably is also causing problems
# for some other projects.
set(LLVM_ENABLE_ABI_BREAKING_CHECKS OFF CACHE BOOL "" FORCE)
set(LLVM_ABI_BREAKING_CHECKS FORCE_OFF CACHE BOOL "" FORCE)

set(CMAKE_REQUIRED_QUIET 1)  # Make the configuration of LLVM quiet

if(ROOT_ARCHITECTURE MATCHES linuxarm64)
  set(ROOT_CLING_TARGET "AArch64")
elseif(ROOT_ARCHITECTURE MATCHES linuxarm)
  set(ROOT_CLING_TARGET "ARM")
elseif(ROOT_ARCHITECTURE MATCHES linuxppc64gcc)
  set(ROOT_CLING_TARGET "PowerPC")
elseif(ROOT_ARCHITECTURE MATCHES linuxs390)
  set(ROOT_CLING_TARGET "SystemZ")
elseif(ROOT_ARCHITECTURE MATCHES linux)
  set(ROOT_CLING_TARGET "X86")
elseif(ROOT_ARCHITECTURE MATCHES macosxarm64)
  set(ROOT_CLING_TARGET "AArch64")
elseif(ROOT_ARCHITECTURE MATCHES macosx)
  set(ROOT_CLING_TARGET "X86")
elseif(ROOT_ARCHITECTURE MATCHES win32)
  set(ROOT_CLING_TARGET "X86")
else()
  set(ROOT_CLING_TARGET "all")
endif()

if(NOT "${ROOT_CLING_TARGET}" STREQUAL "all")
  string(APPEND ROOT_CLING_TARGET ";NVPTX")
endif()

if(MSVC)
  # FIXME: since Visual Studio v16.4.0 the /O2 optimization flag make many (25%) of the tests failing
  # Try to re-enable /O2 after the upgrade of llvm & clang
  string(REPLACE "-O2" "-O1 -Oi" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
  string(REPLACE "-O2" "-O1 -Oi" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
  # replace dashes in the -EH* and -GR* flags with slashes (/EH* /GR*)
  string(REPLACE " -EH" " /EH" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  string(REPLACE " -GR" " /GR" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /ignore:4049,4206,4217,4221")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /ignore:4049,4206,4217,4221")
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} /ignore:4049,4206,4217,4221")
endif()

set(LLVM_TARGETS_TO_BUILD ${ROOT_CLING_TARGET} CACHE STRING "Semicolon-separated list of targets to build, or \"all\".")

if(clingtest)
  message("-- cling test suite enabled: llvm / clang symbols in libCling will be visible!")
  set(CLING_INCLUDE_TESTS ON CACHE BOOL "" FORCE)
  # The path to cling passed through environment variable only relevant when building
  # against external llvm. In that case, cling's testsuite cannot deduce the binary
  # of cling relatively to the llvm tools folder.
  if (NOT builtin_llvm)
    set(CLINGTEST_EXECUTABLE CLING=${CMAKE_CURRENT_BINARY_DIR}/llvm/src/${CMAKE_CFG_INTDIR}/bin/cling)
  endif()
  ROOT_ADD_TEST(clingtest-check-cling COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target check-cling
                                       ENVIRONMENT ${CLINGTEST_EXECUTABLE})
else()
  #---Build LLVM/Clang with symbol visibility=hidden--------------------------------------------------
  set(CMAKE_CXX_VISIBILITY_PRESET hidden)
  set(CMAKE_C_VISIBILITY_PRESET hidden)
endif()

#--- Build LLVM/Clang with modules -----------------------------------------------------------------
if(cxxmodules)
  # LLVM knows how to configure its modules builds. We cannot just add the flags
  # because the cxxmodules builds in llvm have different build dependency order.
  string(REPLACE "${ROOT_CXXMODULES_CXXFLAGS}" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
  string(REPLACE "${ROOT_CXXMODULES_CFLAGS}" "" CMAKE_C_FLAGS ${CMAKE_C_FLAGS})
  if(libcxx)
    # FIXME: We cannot build LLVM/clang with modules on with libstdc++, yet.
    # FIXME: We cannot build LLVM/clang with modules on with libc++, too.
    #set (LLVM_ENABLE_MODULES ON CACHE BOOL "Override the default LLVM_ENABLE_MODULES option value." )
  endif(libcxx)
endif(cxxmodules)

if(gcctoolchain)
  ROOT_ADD_CXX_FLAG(CMAKE_CXX_FLAGS --gcc-toolchain=${gcctoolchain})
endif()

# We will not fix llvm or clang.
string(REPLACE "-Werror " "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")

if(LLVM_SHARED_LINKER_FLAGS)
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${LLVM_SHARED_LINKER_FLAGS}")
endif()
if(LLVM_EXE_LINKER_FLAGS)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LLVM_EXE_LINKER_FLAGS}")
endif()

if(builtin_llvm)
  # Since debug builds of LLVM are quite large, we want to be able
  # to control the build types of ROOT and LLVM independently. The
  # logic below is to make that possible. LLVM is built in Release
  # mode unless a different build type is chosen via LLVM_BUILD_TYPE.

  if(NOT DEFINED LLVM_BUILD_TYPE)
    set(LLVM_BUILD_TYPE Release CACHE STRING "Build type used for LLVM")
  endif()

  message(STATUS "Building LLVM in '${LLVM_BUILD_TYPE}' mode.")

  if(NOT DEFINED LLVM_ENABLE_ASSERTIONS)
    if(CMAKE_BUILD_TYPE MATCHES "(Debug|RelWithDebInfo)"
     OR LLVM_BUILD_TYPE MATCHES "(Debug|RelWithDebInfo)")
      set(LLVM_ENABLE_ASSERTIONS TRUE)
    else()
      set(LLVM_ENABLE_ASSERTIONS FALSE)
    endif()
  endif()

  # Multi-configuration generators ignore CMAKE_BUILD_TYPE, so
  # in that case we set the flags for all configurations to the
  # flags of the build type assigned to LLVM_BUILD_TYPE.

  if(MSVC OR XCODE)
    string(TOUPPER ${LLVM_BUILD_TYPE} LLVM_BUILD_TYPE)
    set(LLVM_C_FLAGS ${CMAKE_C_FLAGS_${LLVM_BUILD_TYPE}})
    set(LLVM_CXX_FLAGS ${CMAKE_CXX_FLAGS_${LLVM_BUILD_TYPE}})
    foreach(CONFIG ${CMAKE_CONFIGURATION_TYPES})
      string(TOUPPER ${CONFIG} CONFIG)
      set(CMAKE_C_FLAGS_${CONFIG} ${LLVM_C_FLAGS})
      set(CMAKE_CXX_FLAGS_${CONFIG} ${LLVM_CXX_FLAGS})
    endforeach()
  elseif(NOT LLVM_BUILD_TYPE STREQUAL CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE ${LLVM_BUILD_TYPE})
  endif()

  set(BUILD_SHARED_LIBS FALSE)

  set(CLANG_ENABLE_STATIC_ANALYZER OFF CACHE BOOL "")
  set(CLANG_ENABLE_ARCMT OFF CACHE BOOL "")
  set(CLANG_ENABLE_FORMAT OFF CACHE BOOL "")

  #---Remove the inherited include_directories()
  set_directory_properties(PROPERTIES INCLUDE_DIRECTORIES "")

  #---Disable VCS writing to prevent that we always have to rebuild when someone touched the git HEAD
  set(LLVM_APPEND_VC_REV OFF CACHE BOOL "Disabled LLVM revision dependency" FORCE)

#---Add the sub-directory excluding all the targets from all-----------------------------------------
  if(CMAKE_GENERATOR MATCHES "Xcode")
    add_subdirectory(llvm/src)
  else()
    add_subdirectory(llvm/src EXCLUDE_FROM_ALL)
  endif()

  #---Make LLVM compilation flags public---------------------------------------------------------------
  # Extract the compilation flags from LLVM and make them public to the
  # rest of ROOT so that we can compile against LLVM with matching flags.

  # LLVM doesn't really give us a API to get this with an in-source build
  # so we just use the normal way of doing this and read the llvm directory
  # compilation properties.
  get_directory_property(LLVM_DEFS DIRECTORY llvm/src COMPILE_DEFINITIONS)
  # Turns DEFINE1;DEFINE2 to -DDEFINE1 -DDEFINE2
  string (REPLACE ";" " -D" LLVM_DEFS ";${LLVM_DEFS}")
  # CACHE variable that contains the additional flags that LLVM needs.
  set(ROOT_LLVM_FLAGS "${LLVM_DEFS}" CACHE STRING "" FORCE)
else()
  # Rely on llvm-config.
  set(CONFIG_OUTPUT)
  find_program(LLVM_CONFIG NAMES "llvm-config" "llvm-config-5.0")
  if(LLVM_CONFIG)
    message(STATUS "Found LLVM_CONFIG as ${LLVM_CONFIG}")
    set(CONFIG_COMMAND ${LLVM_CONFIG}
      "--assertion-mode"
      "--bindir"
      "--libdir"
      "--includedir"
      "--prefix"
      "--src-root"
      "--cmakedir"
      "--build-mode")
    execute_process(
      COMMAND ${CONFIG_COMMAND}
      RESULT_VARIABLE HAD_ERROR
      OUTPUT_VARIABLE CONFIG_OUTPUT
    )
    if(NOT HAD_ERROR)
      string(REGEX REPLACE
        "[ \t]*[\r\n]+[ \t]*" ";"
        CONFIG_OUTPUT ${CONFIG_OUTPUT})
    else()
      string(REPLACE ";" " " CONFIG_COMMAND_STR "${CONFIG_COMMAND}")
      message(STATUS "${CONFIG_COMMAND_STR}")
      message(FATAL_ERROR "llvm-config failed with status ${HAD_ERROR}")
    endif()
  else()
    message(FATAL_ERROR "llvm-config not found -- ${LLVM_CONFIG}")
  endif()

  list(GET CONFIG_OUTPUT 0 ENABLE_ASSERTIONS)
  list(GET CONFIG_OUTPUT 1 TOOLS_BINARY_DIR)
  list(GET CONFIG_OUTPUT 2 LIBRARY_DIR)
  list(GET CONFIG_OUTPUT 3 INCLUDE_DIR)
  list(GET CONFIG_OUTPUT 4 LLVM_OBJ_ROOT)
  list(GET CONFIG_OUTPUT 5 MAIN_SRC_DIR)
  list(GET CONFIG_OUTPUT 6 LLVM_CONFIG_CMAKE_PATH)
  list(GET CONFIG_OUTPUT 7 LLVM_BUILD_MODE)

  message(STATUS "External llvm built in ${LLVM_BUILD_MODE} mode.")

  if(NOT MSVC_IDE)
    set(LLVM_ENABLE_ASSERTIONS ${ENABLE_ASSERTIONS}
      CACHE BOOL "Enable assertions")
    # Assertions should follow llvm-config's.
    mark_as_advanced(LLVM_ENABLE_ASSERTIONS)
  endif()

  set(LLVM_TOOLS_BINARY_DIR ${TOOLS_BINARY_DIR} CACHE PATH "Path to llvm/bin")
  set(LLVM_LIBRARY_DIR ${LIBRARY_DIR} CACHE PATH "Path to llvm/lib")
  set(LLVM_MAIN_INCLUDE_DIR ${INCLUDE_DIR} CACHE PATH "Path to llvm/include")
  set(LLVM_BINARY_DIR ${LLVM_OBJ_ROOT} CACHE PATH "Path to LLVM build tree")
  set(LLVM_MAIN_SRC_DIR ${MAIN_SRC_DIR} CACHE PATH "Path to LLVM source tree")

  # Normalize LLVM_CMAKE_PATH. --cmakedir might contain backslashes.
  # CMake assumes slashes as PATH.
  file(TO_CMAKE_PATH ${LLVM_CONFIG_CMAKE_PATH} LLVM_CMAKE_PATH)

  find_program(LLVM_TABLEGEN_EXE "llvm-tblgen" ${LLVM_TOOLS_BINARY_DIR}
    NO_DEFAULT_PATH)

  set(LLVMCONFIG_FILE "${LLVM_CMAKE_PATH}/LLVMConfig.cmake")
  if(EXISTS ${LLVMCONFIG_FILE})
    list(APPEND CMAKE_MODULE_PATH "${LLVM_CMAKE_PATH}")
    include(${LLVMCONFIG_FILE})
  else()
    message(FATAL_ERROR "Not found: ${LLVMCONFIG_FILE}")
  endif()

  # They are used as destination of target generators.
  # We try to keep these locations consistent with the builtin_llvm. This is important
  # for the LLVMRES target.
  # FIXME: In longer term, we do not really need this and may want to adjust LLVMRES.
  set(LLVM_RUNTIME_OUTPUT_INTDIR ${CMAKE_CURRENT_BINARY_DIR}/llvm/src/${CMAKE_CFG_INTDIR}/bin)
  set(LLVM_LIBRARY_OUTPUT_INTDIR ${CMAKE_CURRENT_BINARY_DIR}/llvm/src/${CMAKE_CFG_INTDIR}/lib)

  if(WIN32 OR CYGWIN)
    # DLL platform -- put DLLs into bin.
    set(LLVM_SHLIB_OUTPUT_INTDIR ${LLVM_RUNTIME_OUTPUT_INTDIR})
  else()
    set(LLVM_SHLIB_OUTPUT_INTDIR ${LLVM_LIBRARY_OUTPUT_INTDIR})
  endif()

  option(LLVM_ENABLE_WARNINGS "Enable compiler warnings." ON)
  option(LLVM_INSTALL_TOOLCHAIN_ONLY
    "Only include toolchain files in the 'install' target." OFF)

  option(LLVM_FORCE_USE_OLD_HOST_TOOLCHAIN
    "Set to ON to force using an old, unsupported host toolchain." OFF)
  option(CLANG_ENABLE_BOOTSTRAP "Generate the clang bootstrap target" OFF)

  include(AddLLVM)
  include(TableGen)
  include(HandleLLVMOptions)
  include(VersionFromVCS)

  set(PACKAGE_VERSION "${LLVM_PACKAGE_VERSION}")
  if (${PACKAGE_VERSION} MATCHES "5\\.0(|\\.[0-9]+)")
    message(STATUS "Using LLVM external library - ${PACKAGE_VERSION}")
  else()
    message(FATAL_ERROR "LLVM version different from ROOT supported, please try 5.0.x")
  endif()

  if (NOT DEFINED LLVM_INCLUDE_TESTS)
    set(LLVM_INCLUDE_TESTS ON)
  endif()

  include_directories("${LLVM_BINARY_DIR}/include" "${LLVM_MAIN_INCLUDE_DIR}")
  link_directories("${LLVM_LIBRARY_DIR}")

#  set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
#  set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib${LLVM_LIBDIR_SUFFIX} )
#  set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib${LLVM_LIBDIR_SUFFIX} )

  if(LLVM_INCLUDE_TESTS OR clingtest)
    set(Python_ADDITIONAL_VERSIONS 2.7)
    include(FindPythonInterp)
    if(NOT PYTHONINTERP_FOUND)
      message(FATAL_ERROR
"Unable to find Python interpreter, required for builds and testing.

Please install Python or specify the PYTHON_EXECUTABLE CMake variable.")
    endif()

    if( ${PYTHON_VERSION_STRING} VERSION_LESS 2.7 )
      message(FATAL_ERROR "Python 2.7 or newer is required")
    endif()

    # Check prebuilt llvm/utils.
    if(EXISTS ${LLVM_TOOLS_BINARY_DIR}/FileCheck${CMAKE_EXECUTABLE_SUFFIX}
        AND EXISTS ${LLVM_TOOLS_BINARY_DIR}/count${CMAKE_EXECUTABLE_SUFFIX}
        AND EXISTS ${LLVM_TOOLS_BINARY_DIR}/not${CMAKE_EXECUTABLE_SUFFIX})
      set(LLVM_UTILS_PROVIDED ON)
    endif()
    set(ROOT_LLVM_MAIN_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/llvm/src/")
    if(EXISTS ${ROOT_LLVM_MAIN_SRC_DIR}/utils/lit/lit.py)
      # Note: path not really used, except for checking if lit was found
      set(LLVM_LIT ${ROOT_LLVM_MAIN_SRC_DIR}/utils/lit/lit.py CACHE PATH "The location of the lit test runner.")
      if(NOT LLVM_UTILS_PROVIDED)
        add_subdirectory(${ROOT_LLVM_MAIN_SRC_DIR}/utils/FileCheck utils/FileCheck)
        add_subdirectory(${ROOT_LLVM_MAIN_SRC_DIR}/utils/count utils/count)
        add_subdirectory(${ROOT_LLVM_MAIN_SRC_DIR}/utils/not utils/not)
        set(LLVM_UTILS_PROVIDED ON)
        set(CLANG_TEST_DEPS FileCheck count not)
      endif()
      # We do not run llvm tests in ROOT.
      #set(UNITTEST_DIR ${ROOT_LLVM_MAIN_SRC_DIR}/utils/unittest)
      #if(EXISTS ${UNITTEST_DIR}/googletest/include/gtest/gtest.h
      #    AND NOT EXISTS ${LLVM_LIBRARY_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX}
      #    AND EXISTS ${UNITTEST_DIR}/CMakeLists.txt)
      #  add_subdirectory(${UNITTEST_DIR} utils/unittest)
      #endif()
    else()
      # Seek installed Lit.
      find_program(LLVM_LIT
                   NAMES llvm-lit lit.py lit
                   PATHS "${ROOT_LLVM_MAIN_SRC_DIR}/utils/lit"
                   DOC "Path to lit.py")
    endif()

    if(LLVM_LIT)
      # Define the default arguments to use with 'lit', and an option for the user
      # to override.
      set(LIT_ARGS_DEFAULT "-sv")
      if (MSVC OR XCODE)
        set(LIT_ARGS_DEFAULT "${LIT_ARGS_DEFAULT} --no-progress-bar")
      endif()
      set(LLVM_LIT_ARGS "${LIT_ARGS_DEFAULT}" CACHE STRING "Default options for lit")
      set(LIT_COMMAND "${LLVM_LIT}" CACHE STRING "Path to the LLVM LIT.")

      # On Win32 hosts, provide an option to specify the path to the GnuWin32 tools.
      if( WIN32 AND NOT CYGWIN )
        set(LLVM_LIT_TOOLS_DIR "" CACHE PATH "Path to GnuWin32 tools")
      endif()
    else()
      set(LLVM_INCLUDE_TESTS OFF)
    endif()
  endif()

  set( CLANG_BUILT_STANDALONE 1 )
  set(BACKEND_PACKAGE_STRING "LLVM ${LLVM_PACKAGE_VERSION}")

  # must be set before add_subdirectory(cling): it can change the value of LLVM_BINARY_DIR
  set(LLVM_DIR "${LLVM_BINARY_DIR}")
  if (builtin_clang)
    # For builtin LLVM this is set in interpreter/llvm/src/CMakeLists.txt
    set(Clang_DIR ${LLVM_BINARY_DIR}/tools/clang/)
    add_subdirectory(llvm/src/tools/clang EXCLUDE_FROM_ALL)
  else()
    set(Clang_DIR ${LLVM_BINARY_DIR}/lib/cmake/clang/)
    add_subdirectory(cling)
  endif()

endif(builtin_llvm)

#---Export the include directories------------------------------------------------------------------
if (builtin_clang)
  set(CLANG_INCLUDE_DIRS
    ${CMAKE_SOURCE_DIR}/interpreter/llvm/src/tools/clang/include
    ${CMAKE_BINARY_DIR}/interpreter/llvm/src/tools/clang/include
    CACHE STRING "Clang include directories.")
endif()

if(builtin_llvm)
  set(LLVM_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/interpreter/llvm/src/include
    ${CMAKE_BINARY_DIR}/interpreter/llvm/src/include
    CACHE STRING "LLVM include directories."
    )
  #---Set into parent scope LLVM_VERSION --------------------------------------------------------------
  set(llvmconfigfile "${CMAKE_CURRENT_SOURCE_DIR}/llvm/src/CMakeLists.txt")
  file(READ ${llvmconfigfile} _filestr)
  foreach(versionpart LLVM_VERSION_MAJOR LLVM_VERSION_MINOR LLVM_VERSION_PATCH)
    if("${_filestr}" MATCHES "set[(](${versionpart}[^)]+)")
      string(REGEX REPLACE " *${versionpart} ([^)]+).*" "\\1" ${versionpart} "${CMAKE_MATCH_1}")
    else()
      MESSAGE(FATAL_ERROR "Cannot extract LLVM version number (${versionpart}) from ${llvmconfigfile}")
    endif()
  endforeach()
  set(LLVM_VERSION "${LLVM_VERSION_MAJOR}.${LLVM_VERSION_MINOR}.${LLVM_VERSION_PATCH}" PARENT_SCOPE)
else()
  set(LLVM_INCLUDE_DIRS ${LLVM_MAIN_INCLUDE_DIR}
    CACHE STRING "System LLVM include directories."
    )
  # We cannot include LLVMConfig.cmake here. We reply on builtin clang to discover the available
  # llvm version and set it in clang version.
  # FIXME: This is a very ugly workaround. Is there a better way to do this?
  set(LLVM_VERSION ${CLANG_VERSION} PARENT_SCOPE)
endif()

#---Avoid to dependencies to system (e.g. atomic) libraries without modifying the LLVM code----------
if(CMAKE_SYSTEM_NAME STREQUAL Darwin)
  get_target_property(__libs LLVMSupport INTERFACE_LINK_LIBRARIES)
  list(REMOVE_ITEM __libs atomic)
  set_target_properties(LLVMSupport PROPERTIES INTERFACE_LINK_LIBRARIES "${__libs}")
endif()

#---Mark the LLVM/CLANG variables as advanced--------------------------------------------------------
get_cmake_property(variables CACHE_VARIABLES)
foreach(var ${variables})
  if(var MATCHES "^(CLANG|LLVM|CLING)_")
    mark_as_advanced(FORCE ${var})
  endif()
endforeach()
mark_as_advanced(FORCE BUG_REPORT_URL BUILD_CLANG_FORMAT_VS_PLUGIN BUILD_SHARED_LIBS BUILD_TESTING
                       C_INCLUDE_DIRS DEFAULT_SYSROOT FFI_INCLUDE_DIR FFI_LIBRARY_DIR
                       GCC_INSTALL_PREFIX LIBCLANG_BUILD_STATIC TOOL_INFO_PLIST)
mark_as_advanced(CLEAR LLVM_ENABLE_ASSERTIONS LLVM_BUILD_TYPE)
