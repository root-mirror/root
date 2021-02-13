# Make a shared library that holds the ROOT default config for address sanitizer.
# This is can be used with LD_PRELOAD
add_library(ROOTSanitizerConfig SHARED SanitizerSetup.cxx)

# Make a static library that holds the ROOT default config for address sanitizer.
# We link this into all executables in ROOT, so they automatically start up with good defaults.
set(library ROOTStaticSanitizerConfig)
add_library(${library} STATIC SanitizerSetup.cxx)
set_property(TARGET ${library} PROPERTY POSITION_INDEPENDENT_CODE ON)

# This pulls in the sanitizer link flags when linking against the config library
target_link_libraries(${library} INTERFACE ${ASAN_EXTRA_EXE_LINKER_FLAGS})

# Make it visible to the outside to sanitize e.g. roottest executables
set_property(GLOBAL APPEND PROPERTY ROOT_EXPORTED_TARGETS ${library})
add_library(ROOT::${library} ALIAS ${library})

# Now sanitize executables that are not created with ROOT_EXECUTABLE():
set_property(TARGET llvm-tblgen clang-tblgen
  APPEND PROPERTY LINK_LIBRARIES ${library})
add_dependencies(llvm-tblgen ROOTStaticSanitizerConfig)
