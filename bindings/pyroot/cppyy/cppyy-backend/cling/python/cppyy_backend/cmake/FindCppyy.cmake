#.rst:
# FindCppyy
# -------
#
# Find Cppyy
#
# This module finds an installed Cppyy.  It sets the following variables:
#
# ::
#
#   Cppyy_FOUND - set to true if Cppyy is found
#   Cppyy_DIR - the directory where Cppyy is installed
#   Cppyy_EXECUTABLE - the path to the cppyy-generator executable
#   Cppyy_INCLUDE_DIRS - Where to find the ROOT header files.
#   Cppyy_VERSION - the version number of the Cppyy backend.
#
#
# The module also defines the following functions:
#
#   cppyy_add_bindings - Generate a set of bindings from a set of header files.
#
# The minimum required version of Cppyy can be specified using the
# standard syntax, e.g.  find_package(Cppyy 4.19)
#
#

get_filename_component(BACKEND_PREFIX "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(BACKEND_PREFIX "${BACKEND_PREFIX}" PATH)
if(BACKEND_PREFIX STREQUAL "/")
  set(BACKEND_PREFIX "")
endif()

find_program(Cppyy_EXECUTABLE NAMES rootcling)

if(CPPYY_MODULE_PATH)
    #
    # Cppyy_DIR: one level above the installed cppyy cmake module path
    #
    set(Cppyy_DIR ${CPPYY_MODULE_PATH}/../)
    #
    # Cppyy_INCLUDE_DIRS: Directory with cppyy H_FILES
    #
    set(Cppyy_INCLUDE_DIRS ${Cppyy_DIR}/include)
    #
    # Cppyy_VERSION.
    #
    find_package(ROOT QUIET REQUIRED PATHS ${CPPYY_MODULE_PATH})
    if(ROOT_FOUND)
        set(Cppyy_VERSION ${ROOT_VERSION})
    endif()
endif()

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Cppyy
                                  REQUIRED_VARS Cppyy_EXECUTABLE Cppyy_DIR Cppyy_INCLUDE_DIRS CPPYY_MODULE_PATH
                                  VERSION_VAR Cppyy_VERSION
)
mark_as_advanced(Cppyy_VERSION)

#
# Generate setup.py from the setup.py.in template.
#
function(cppyy_generate_setup pkg version lib_so_file rootmap_file pcm_file map_file)
    set(SETUP_PY_FILE ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
    set(CPPYY_PKG ${pkg})
    get_filename_component(CPPYY_LIB_SO ${lib_so_file} NAME)
    get_filename_component(CPPYY_ROOTMAP ${rootmap_file} NAME)
    get_filename_component(CPPYY_PCM ${pcm_file} NAME)
    get_filename_component(CPPYY_MAP ${map_file} NAME)
    configure_file(${BACKEND_PREFIX}/pkg_templates/setup.py.in ${SETUP_PY_FILE})

    set(SETUP_PY_FILE ${SETUP_PY_FILE} PARENT_SCOPE)
endfunction(cppyy_generate_setup)

#
# Generate a packages __init__.py using the __init__.py.in template.
#
function(cppyy_generate_init)
    set(simple_args PKG LIB_FILE MAP_FILE)
    set(list_args NAMESPACES)
    cmake_parse_arguments(ARG
                          ""
                          "${simple_args}"
                          "${list_args}"
                          ${ARGN}
    )

    set(INIT_PY_FILE ${CMAKE_CURRENT_BINARY_DIR}/${ARG_PKG}/__init__.py)
    set(CPPYY_PKG ${ARG_PKG})
    get_filename_component(CPPYY_LIB_SO ${ARG_LIB_FILE} NAME)
    get_filename_component(CPPYY_MAP ${ARG_MAP_FILE} NAME)

    string(REPLACE "${ARG_NAMESPACES}" ";" ", " _namespaces)

    if(NOT "${ARG_NAMESPACES}" STREQUAL "")
        string(REPLACE "${ARG_NAMESPACES}" ";" ", " _namespaces)
        set(NAMESPACE_INJECTIONS "from cppyy.gbl import ${_namespaces}")
    else()
        set(NAMESPACE_INJECTIONS "")
    endif()

    configure_file(${BACKEND_PREFIX}/pkg_templates/__init__.py.in ${INIT_PY_FILE})

    set(INIT_PY_FILE ${INIT_PY_FILE} PARENT_SCOPE)
endfunction(cppyy_generate_init)

#
# Generate a set of bindings from a set of header files. Somewhat like CMake's
# add_library(), the output is a compiler target. In addition ancilliary files
# are also generated to allow a complete set of bindings to be compiled,
# packaged and installed.
#
#   cppyy_add_bindings(
#       pkg
#       pkg_version
#       author
#       author_email
#       [URL url]
#       [LICENSE license]
#       [LANGUAGE_STANDARD std]
#       [LINKDEFS linkdef...]
#       [IMPORTS pcm...]
#       [GENERATE_OPTIONS option...]
#       [COMPILE_OPTIONS option...]
#       [INCLUDE_DIRS dir...]
#       [LINK_LIBRARIES library...]
#       [H_DIRS H_DIRSectory]
#       H_FILES h_file...)
#
# The bindings are based on https://cppyy.readthedocs.io/en/latest/, and can be
# used as per the documentation provided via the cppyy.gbl namespace. First add
# the directory of the <pkg>.rootmap file to the LD_LIBRARY_PATH environment
# variable, then "import cppyy; from cppyy.gbl import <some-C++-entity>".
#
# Alternatively, use "import <pkg>". This convenience wrapper supports
# "discovery" of the available C++ entities using, for example Python 3's command
# line completion support.
#
# This function creates setup.py, setup.cfg, and MANIFEST.in appropriate
# for the package in the build directory. It also creates the package directory PKG,
# and within it a tests subdmodule PKG/tests/test_bindings.py to sanity test the bindings.
# Further, it creates PKG/pythonizors/, which can contain files of the form
# pythonize_*.py, with functions of the form pythonize_<NAMESPACE>_*.py, which will
# be consumed by the initialization routine and added as pythonizors for their associated
# namespace on import.
#
# The setup.py and setup.cfg are prepared to create a Wheel. They can be customized
# for the particular package by modifying the templates in pkg_templates/.
#
# The bindings are generated/built/packaged using 3 environments:
#
#   - One compatible with the header files being bound. This is used to
#     generate the generic C++ binding code (and some ancilliary files) using
#     a modified C++ compiler. The needed options must be compatible with the
#     normal build environment of the header files.
#
#   - One to compile the generated, generic C++ binding code using a standard
#     C++ compiler. The resulting library code is "universal" in that it is
#     compatible with both Python2 and Python3.
#
#   - One to package the library and ancilliary files into standard Python2/3
#     wheel format. The packaging is done using native Python tooling.
#
# Arguments and options:
#
#   pkg                 The name of the package to generate. This can be either
#                       of the form "simplename" (e.g. "Akonadi"), or of the
#                       form "namespace.simplename" (e.g. "KF5.Akonadi").
#
#   pkg_version         The version of the package.
#
#   author              The name of the bindings author.
#
#   author_email        The email address of the bindings author.
#
#   URL url             The home page for the library or bindings. Default is
#                       "https://pypi.python.org/pypi/<pkg>".
#
#   LICENSE license     The license, default is "MIT".
#
#   LICENSE_FILE        Path to license file to include in package. Default is LICENSE.
#
#   README              Path to README file to include in package and use as
#                       text for long_description. Default is README.rst.
#
#   IMPORTS pcm         Files which contain previously-generated bindings
#                       which pkg depends on.
#
#   LANGUAGE_STANDARD std
#                       The version of C++ in use, "14" by default.
#
#   GENERATE_OPTIONS option
#                       Options which will be passed to the rootcling invocation
#                       in the cppyy-generate utility. cppyy-generate is used to
#                       create the bindings map.
#
#   LINKDEFS def        Files or lines which contain extra #pragma content
#                       for the linkdef.h file used by rootcling. See
#                       https://root.cern.ch/root/html/guides/users-guide/AddingaClass.html#the-linkdef.h-file.
#
#   EXTRA_CODES code    Files which contain extra code needed by the bindings.
#                       Customisation is by routines named "c13n_<something>";
#                       each such routine is passed the module for <pkg>:
#
#                           def c13n_doit(pkg_module):
#                               print(pkg_module.__dict__)
#
#                       The files and individual routines within files are
#                       processed in alphabetical order.
#
#   EXTRA_HEADERS hdr   Files which contain extra headers needed by the bindings.
#
#   COMPILE_OPTIONS option
#                       Options which are to be passed into the compile/link
#                       command.
#
#   INCLUDE_DIRS dir    Include directories.
#
#   LINK_LIBRARIES library
#                       Libraries to link against.
#
#   NAMESPACES          List of C++ namespaces which should be imported into the
#                       bindings' __init__.py. This avoids having to write imports
#                       of the form `from PKG import NAMESPACE`.
#
#   EXTRA_PKG_FILES     Extra files to copy into the package. Note that non-python
#                       files will need to be added to the MANIFEST.in.in template.
#
#   H_DIRS directory    Base directories for H_FILES.
#
#   H_FILES h_file      Header files for which to generate bindings in pkg.
#                       Absolute filenames, or filenames relative to H_DIRS. All
#                       definitions found directly in these files will contribute
#                       to the bindings. (NOTE: This means that if "forwarding
#                       headers" are present, the real "legacy" headers must be
#                       specified as H_FILES).
#
#                       All header files which contribute to a given C++ namespace
#                       should be grouped into a single pkg to ensure a 1-to-1
#                       mapping with the implementing Python class.
#
# Returns via PARENT_SCOPE variables:
#
#   CPPYY_LIB_TARGET    The target cppyy bindings shared library.
#
#   SETUP_PY_FILE       The generated setup.py.
#
#   INIT_PY_FILE       The generated setup.py.
#
#   PY_WHEEL_FILE       The finished .whl package file
#
# Examples:
#
# cppyy_add_bindings(
#    "${PROJECT_NAME}" "${PROJECT_VERSION}" "user" "user@gmail.com"
#    LANGUAGE_STANDARD "14"
#    GENERATE_OPTIONS "-D__PIC__;-Wno-macro-redefined"
#    INCLUDE_DIRS     ${PCL_INCLUDE_DIRS}
#    LINKDEFS         LinkDef.h
#    LINK_LIBRARIES   ${PCL_LIBRARIES}
#    H_DIRS           ${HEADER_PATH}
#    H_FILES          ${LIB_HEADERS}
#    NAMESPACES       pcl
# )

function(cppyy_add_bindings pkg pkg_version author author_email)
    set(simple_args URL LICENSE LICENSE_FILE LANGUAGE_STANDARD
        README_FILE)
    set(list_args IMPORTS GENERATE_OPTIONS COMPILE_OPTIONS INCLUDE_DIRS LINK_LIBRARIES H_DIRS H_FILES LINKDEFS EXTRA_CODES EXTRA_HEADERS NAMESPACES)
    cmake_parse_arguments(
        ARG
        ""
        "${simple_args}"
        "${list_args}"
        ${ARGN})
    if(NOT "${ARG_UNPARSED_ARGUMENTS}" STREQUAL "")
        message(SEND_ERROR "Unexpected arguments specified '${ARG_UNPARSED_ARGUMENTS}'")
    endif()
    string(REGEX MATCH "[^\.]+$" pkg_simplename ${pkg})
    string(REGEX REPLACE "\.?${pkg_simplename}" "" pkg_namespace ${pkg})
    set(pkg_dir ${CMAKE_CURRENT_BINARY_DIR})
    string(REPLACE "." "/" tmp ${pkg})
    set(pkg_dir "${pkg_dir}/${tmp}")
    set(lib_name "${pkg_namespace}${pkg_simplename}Cppyy")
    set(lib_file ${CMAKE_SHARED_LIBRARY_PREFIX}${lib_name}${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(cpp_file ${CMAKE_CURRENT_BINARY_DIR}/${pkg_simplename}.cpp)
    set(pcm_file ${pkg_dir}/${pkg_simplename}_rdict.pcm)
    set(rootmap_file ${pkg_dir}/${pkg_simplename}.rootmap)
    set(extra_map_file ${pkg_dir}/${pkg_simplename}.map)

    ############################################################
    #
    # Package metadata.
    #
    # License
    if("${ARG_URL}" STREQUAL "")
        string(REPLACE "." "-" tmp ${pkg})
        set(ARG_URL "https://pypi.python.org/pypi/${tmp}")
    endif()
    if("${ARG_LICENSE}" STREQUAL "")
        set(ARG_LICENSE "LGPL2.1")
    endif()
    set(BINDINGS_LICENSE ${ARG_LICENSE})

    # License file
    if("${ARG_LICENSE_FILE}" STREQUAL "")
        set(ARG_LICENSE_FILE ${CMAKE_SOURCE_DIR}/LICENSE)
    endif()
    set(LICENSE_FILE ${ARG_LICENSE_FILE})

    # ReadMe file
    if("${ARG_README_FILE}" STREQUAL "")
        set(ARG_README_FILE ${CMAKE_SOURCE_DIR}/README.rst)
    endif()
    set(README_FILE ${ARG_README_FILE})

    #
    # Language standard.
    #
    if("${ARG_LANGUAGE_STANDARD}" STREQUAL "")
        set(ARG_LANGUAGE_STANDARD "14")
    endif()

    ################################################################
    #
    # Make H_FILES with absolute paths.
    #
    if("${ARG_H_FILES}" STREQUAL "")
    message(SEND_ERROR "No H_FILES specified")
    endif()
    if("${ARG_H_DIRS}" STREQUAL "")
    set(ARG_H_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
    endif()
    set(tmp)
    foreach(h_file IN LISTS ARG_H_FILES)
    if(NOT IS_ABSOLUTE ${h_file})
        foreach(h_dir IN LISTS ARG_H_DIRS)
        if(EXISTS ${h_dir}/${h_file})
            set(h_file ${h_dir}/${h_file})
            break()
        endif()
        endforeach(h_dir)
    endif()
    if(NOT EXISTS ${h_file})
        message(WARNING "H_FILES ${h_file} does not exist")
    endif()
    list(APPEND tmp ${h_file})
    endforeach(h_file)
    set(ARG_H_FILES ${tmp})

    ###################################################
    #
    # Create the linkdef.h file with a line for each h_file.
    #
    set(out_linkdef ${CMAKE_CURRENT_BINARY_DIR}/linkdef.h)
    file(WRITE ${out_linkdef} "/* Per H_FILES entries: */\n")
    foreach(h_file IN LISTS ARG_H_FILES)
        #
        # Doubled-up path separators "//" causes errors in rootcling.
        #
        string(REGEX REPLACE "/+" "/" h_file ${h_file})
        file(APPEND ${out_linkdef} "#pragma link C++ defined_in ${h_file};\n")
    endforeach(h_file)
    foreach(h_file IN LISTS ARG_EXTRA_HEADERS)
        #
        # Doubled-up path separators "//" causes errors in rootcling.
        #
        string(REGEX REPLACE "/+" "/" h_file ${h_file})
        file(APPEND ${out_linkdef} "#pragma extra_include \"${h_file}\";\n")
    endforeach(h_file)

    #
    # Append any manually-provided linkdef.h content.
    #
    set(LINKDEF_EXTRACTS)
    string(REPLACE "\\" "\\\\" ARG_LINKDEFS "${ARG_LINKDEFS}")
    foreach(in_linkdef IN LISTS ARG_LINKDEFS)
        if("${in_linkdef}" STREQUAL "")
        continue()
        endif()
        if(NOT IS_ABSOLUTE "${in_linkdef}" AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${in_linkdef}")
        set(in_linkdef "${CMAKE_CURRENT_SOURCE_DIR}/${in_linkdef}")
        endif()
        if(EXISTS "${in_linkdef}")
        file(APPEND ${out_linkdef} "/* Copied from ${in_linkdef}: */\n")
        file(STRINGS ${in_linkdef} in_linkdef NEWLINE_CONSUME)
        else()
        file(APPEND ${out_linkdef} "/* Inlined: */\n")
        endif()
        string(REPLACE "\n" ";" in_linkdef "${in_linkdef}")
        foreach(line ${in_linkdef})
        file(APPEND ${out_linkdef} "${line}\n")
        endforeach()
        list(GET in_linkdef 0 in_linkdef)
        list(APPEND LINKDEFS_EXTRACTS ${in_linkdef})
    endforeach(in_linkdef)
    #
    # Record diagnostics.
    #
    file(APPEND ${out_linkdef} "//\n// Diagnostics.\n//\n")
    foreach(arg IN LISTS simple_args list_args)
        if(arg STREQUAL "LINKDEFS")
        file(APPEND ${out_linkdef} "// ${arg}=\n")
        foreach(in_linkdef IN LISTS LINKDEFS_EXTRACTS)
            file(APPEND ${out_linkdef} "//    ${in_linkdef}...\n")
        endforeach(in_linkdef)
        else()
        file(APPEND ${out_linkdef} "// ${arg}=${ARG_${arg}}\n")
        endif()
    endforeach(arg)
  ####################################################
    #
    # Set up common args.
    #
    list(APPEND ARG_GENERATE_OPTIONS "-std=c++${ARG_LANGUAGE_STANDARD}")
    foreach(dir ${ARG_H_DIRS} ${ARG_INCLUDE_DIRS})
        list(APPEND ARG_GENERATE_OPTIONS "-I${dir}")
    endforeach(dir)

    #
    # Set up arguments for rootcling.
    #
    set(cling_args)
    list(APPEND cling_args "-f" ${cpp_file})
    list(APPEND cling_args "-s" ${pkg_simplename})
    list(APPEND cling_args "-rmf" ${rootmap_file} "-rml" ${lib_file})
    foreach(in_pcm IN LISTS ARG_IMPORTS)
        #
        # Create -m options for any imported .pcm files.
        #
        list(APPEND cling_args "-m" "${in_pcm}")
    endforeach(in_pcm)
    list(APPEND cling_args "${ARG_GENERATE_OPTIONS}")

    # run rootcling
    add_custom_command(OUTPUT ${cpp_file} ${pcm_file} ${rootmap_file}
    COMMAND ${Cppyy_EXECUTABLE} ${cling_args} ${ARG_H_FILES} ${out_linkdef} WORKING_DIRECTORY ${pkg_dir})


    ############### cppyy-generator #######################
    find_package(LibClang REQUIRED)
    get_filename_component(Cppyygen_EXECUTABLE ${Cppyy_EXECUTABLE} DIRECTORY)
    set(Cppyygen_EXECUTABLE ${Cppyygen_EXECUTABLE}/cppyy-generator)

    #
    # Set up arguments for cppyy-generator.
    #
    if(${CONDA_ACTIVE})
        set(CLANGDEV_INCLUDE $ENV{CONDA_PREFIX}/lib/clang/${CLANG_VERSION_STRING}/include)
        message(STATUS "adding conda clangdev includes to cppyy-generator options (${CLANGDEV_INCLUDE})")
        list(APPEND ARG_GENERATE_OPTIONS "-I${CLANGDEV_INCLUDE}")
    endif()

    #
    # Run cppyy-generator. First check dependencies.
    #
    set(generator_args)
    foreach(arg IN LISTS ARG_GENERATE_OPTIONS)
        string(REGEX REPLACE "^-" "\\\\-" arg ${arg})
        list(APPEND generator_args ${arg})
    endforeach()

    add_custom_command(OUTPUT ${extra_map_file}
                       COMMAND ${LibClang_PYTHON_EXECUTABLE} ${Cppyygen_EXECUTABLE}
                       --libclang ${LibClang_LIBRARY} --flags "\"${generator_args}\""
                               ${extra_map_file} ${ARG_H_FILES} WORKING_DIRECTORY ${pkg_dir}
                       DEPENDS ${ARG_H_FILES}
                       WORKING_DIRECTORY ${pkg_dir}
    )


    #################################################

    #
    # Compile/link.
    #
    add_library(${lib_name} SHARED ${cpp_file} ${pcm_file} ${rootmap_file} ${extra_map_file} ${ARG_EXTRA_CODES})
    set_target_properties(${lib_name} PROPERTIES LINKER_LANGUAGE CXX)
    set_property(TARGET ${lib_name} PROPERTY VERSION ${version})
    set_property(TARGET ${lib_name} PROPERTY CXX_STANDARD ${ARG_LANGUAGE_STANDARD})
    set_property(TARGET ${lib_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY ${pkg_dir})
    set_property(TARGET ${lib_name} PROPERTY LINK_WHAT_YOU_USE TRUE)
    target_include_directories(${lib_name} PRIVATE ${Cppyy_INCLUDE_DIRS} ${ARG_INCLUDE_DIRS})
    target_compile_options(${lib_name} PRIVATE ${ARG_COMPILE_OPTIONS})
    target_link_libraries(${lib_name} PUBLIC ${ARG_LINK_LIBRARIES})

    ####################################################
    #
    # Generate __init__.py
    #
    cppyy_generate_init(PKG        ${pkg}
                        LIB_FILE   ${lib_file}
                        MAP_FILE   ${extra_map_file}
                        NAMESPACES ${ARG_NAMESPACES}
    )
    set(INIT_PY_FILE ${INIT_PY_FILE} PARENT_SCOPE)

    #
    # Generate setup.py
    #
    cppyy_generate_setup(${pkg}
                         ${pkg_version}
                         ${lib_file}
                         ${rootmap_file}
                         ${pcm_file}
                         ${extra_map_file}
    )

    #
    # Generate setup.cfg
    #
    set(setup_cfg ${CMAKE_CURRENT_BINARY_DIR}/setup.cfg)
    configure_file(${BACKEND_PREFIX}/pkg_templates/setup.cfg.in ${setup_cfg})

    #
    # Copy initializor
    #
    set(initializor ${CMAKE_CURRENT_BINARY_DIR}/initializor.py)
    file(COPY ${BACKEND_PREFIX}/pkg_templates/initializor.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${pkg} USE_SOURCE_PERMISSIONS)

    #
    # Copy README and LICENSE
    #
    if (EXISTS ${README_FILE})
        file(COPY ${README_FILE}  DESTINATION . USE_SOURCE_PERMISSIONS)
    endif()
    if (EXISTS ${LICENSE_FILE})
        file(COPY ${LICENSE_FILE} DESTINATION . USE_SOURCE_PERMISSIONS)
    endif()

    #
    # Generate a pytest/nosetest sanity test script.
    #
    set(PKG ${pkg})
    configure_file(${BACKEND_PREFIX}/pkg_templates/test_bindings.py.in ${pkg_dir}/tests/test_bindings.py)

    #
    # Generate MANIFEST.in
    #
    configure_file(${BACKEND_PREFIX}/pkg_templates/MANIFEST.in.in ${CMAKE_CURRENT_BINARY_DIR}/MANIFEST.in)

    #
    # Copy pure python code
    #
    file(COPY ${CMAKE_SOURCE_DIR}/py/ DESTINATION ${pkg_dir}
         USE_SOURCE_PERMISSIONS
         FILES_MATCHING PATTERN "*.py")

    #
    # Copy any extra files into package.
    #
    file(COPY ${ARG_EXTRA_FILES} DESTINATION ${pkg_dir} USE_SOURCE_PERMISSIONS)

    #
    # Kinda ugly: you'e not really supposed to glob like this. Oh well. Using this to set
    # dependencies for the python wheel building command; the file copy above is done on every
    # cmake invocation anyhow.
    #
    # Then, get the system architecture and build the wheel string based on PEP 427.
    #
    file(GLOB_RECURSE PY_PKG_FILES
         LIST_DIRECTORIES FALSE
         CONFIGURE_DEPENDS
         "${CMAKE_SOURCE_DIR}/py/*.py")
    string(TOLOWER ${CMAKE_SYSTEM_NAME} SYSTEM_STR)
    set(pkg_whl "${CMAKE_BINARY_DIR}/dist/${pkg}-${pkg_version}-py3-none-${SYSTEM_STR}_${CMAKE_SYSTEM_PROCESSOR}.whl")
    add_custom_command(OUTPUT  ${pkg_whl}
                       COMMAND ${LibClang_PYTHON_EXECUTABLE} setup.py bdist_wheel
                       DEPENDS ${SETUP_PY_FILE} ${lib_name} ${setup_cfg}
                       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    add_custom_target(wheel ALL
                      DEPENDS ${pkg_whl}
    )
    add_dependencies(wheel ${lib_name})

    #
    # Return results.
    #
    set(CPPYY_LIB_TARGET    ${lib_name} PARENT_SCOPE)
    set(SETUP_PY_FILE       ${SETUP_PY_FILE} PARENT_SCOPE)
    set(INIT_PY_FILE        ${INIT_PY_FILE} PARENT_SCOPE)
    set(PY_WHEEL_FILE       ${pkg_whl}  PARENT_SCOPE)
endfunction(cppyy_add_bindings)


#
# Return a list of available pip programs.
#
function(cppyy_find_pips)
    execute_process(
        COMMAND python -c "from cppyy_backend import bindings_utils; print(\";\".join(bindings_utils.find_pips()))"
        OUTPUT_VARIABLE _stdout
        ERROR_VARIABLE _stderr
        RESULT_VARIABLE _rc
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT "${_rc}" STREQUAL "0")
        message(FATAL_ERROR "Error finding pips: (${_rc}) ${_stderr}")
    endif()
    set(PIP_EXECUTABLES ${_stdout} PARENT_SCOPE)
endfunction(cppyy_find_pips)
