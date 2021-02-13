#.rst:
# FindLibClang
# ------------
#
# Find LibClang
#
# Find LibClang headers and library
#
# ::
#
#   LibClang_FOUND             - True if libclang is found.
#   LibClang_LIBRARY           - Clang library to link against.
#   LibClang_VERSION           - Version number as a string (e.g. "3.9").
#   LibClang_PYTHON_EXECUTABLE - Compatible python version.

#
# Python support for clang might not be available for Python3. We need to
# find what we have.
#
find_library(LibClang_LIBRARY libclang.so PATH_SUFFIXES $ENV{CONDA_PREFIX}/lib x86_64-linux-gnu llvm llvm/6/lib64 llvm-6.0/lib)
function(_find_libclang_python python_executable)
    #
    # Prefer python3 explicitly or implicitly over python2.
    #
    foreach(exe IN ITEMS python3 python python2)
        execute_process(
            COMMAND ${exe} -c "from clang.cindex import Config; Config.set_library_file(\"${LibClang_LIBRARY}\"); Config().lib"
            ERROR_VARIABLE _stderr
            RESULT_VARIABLE _rc
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        if("${_rc}" STREQUAL "0")
            set(pyexe ${exe})
            break()
        endif()
    endforeach()
    set(${python_executable} "${pyexe}" PARENT_SCOPE)
endfunction(_find_libclang_python)

_find_libclang_python(LibClang_PYTHON_EXECUTABLE)
if(LibClang_LIBRARY)
    set(LibClang_LIBRARY ${LibClang_LIBRARY})
    string(REGEX REPLACE ".*clang-\([0-9]+.[0-9]+\).*" "\\1" LibClang_VERSION_TMP "${LibClang_LIBRARY}")
    set(LibClang_VERSION ${LibClang_VERSION_TMP} CACHE STRING "LibClang version" FORCE)
endif()

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LibClang REQUIRED_VARS  LibClang_LIBRARY LibClang_PYTHON_EXECUTABLE
                                           VERSION_VAR    LibClang_VERSION)

find_program(CLANG_EXE clang++)
EXECUTE_PROCESS( COMMAND ${CLANG_EXE} --version OUTPUT_VARIABLE clang_full_version_string )
string (REGEX REPLACE ".*clang version ([0-9]+\\.[0-9]+\\.[0-9]+).*" "\\1" CLANG_VERSION_STRING "${clang_full_version_string}")

set(CLANG_VERSION_STRING ${CLANG_VERSION_STRING} PARENT_SCOPE)

mark_as_advanced(LibClang_VERSION)
unset(_filename)
unset(_find_libclang_filename)
