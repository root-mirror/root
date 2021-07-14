# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.


set_package_properties(Doxygen PROPERTIES
    URL "http://www.doxygen.org"
    DESCRIPTION "Documentation generator"
    TYPE OPTIONAL
    PURPOSE "Requires CMake 3.9, adds docs target")

find_package(Doxygen OPTIONAL_COMPONENTS dot)

# Only add docs target if Doxygen present
if(TARGET Doxygen::doxygen AND TARGET Doxygen::dot)
    # Setting DOXYGEN settings
    set(DOXYGEN_GENERATE_LATEX YES)
    set(DOXYGEN_SOURCE_BROWSER YES)
    set(DOXYGEN_EXTRACT_ALL YES)
    set(DOXYGEN EXTRACT_PRIVATE YES)
    set(DOXYGEN_EXTRACT_STATIC YES)

    get_target_property(MATH_SOURCES Minuit2Math SOURCES)
    get_target_property(MINUIT2_SOURCES Minuit2 SOURCES)

    file(READ Minuit2.md MINUIT2_MAINPAGE)
    string(REPLACE "\\page Minuit2Page" "\\mainpage" MINUIT2_MAINPAGE "${MINUIT2_MAINPAGE}")
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/index.md" "${MINUIT2_MAINPAGE}")

    # Actually add the docs target
    doxygen_add_docs(docs
        "${CMAKE_CURRENT_BINARY_DIR}/index.md"
        ${MATH_SOURCES}
        ${MINUIT2_SOURCES}
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/.."
        )
else()
    message(STATUS "Doxygen and/or Dot not found, not adding docs target")
endif()
