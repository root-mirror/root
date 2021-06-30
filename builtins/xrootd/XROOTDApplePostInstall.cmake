# Post install routine needed to change the LC_RPATH variable of the XRootD
# libraries on macOS after installation. The function changes the value from
# ${build_libdir} to ${install_libdir} with the `install_name_tool` executable.
# On successive reiterations of the `cmake --install` command,
# `install_name_tool` would error out since LC_RPATH would be already set to
# ${install_libdir}. Since that is what we want in the end, we discard these
# errors.

function(xrootd_libs_change_rpath build_libdir install_libdir)

  file(GLOB XROOTD_ALL_LIBRARIES "${install_libdir}/libXrd*")

  find_program(INSTALL_NAME_TOOL install_name_tool)

  if(INSTALL_NAME_TOOL)
    message(STATUS "Found tool ${INSTALL_NAME_TOOL}")
    message(STATUS "Adjusting LC_RPATH variable of XRootD libraries in ${install_libdir}")

    foreach(XRD_LIB_PATH ${XROOTD_ALL_LIBRARIES})
      execute_process(COMMAND ${INSTALL_NAME_TOOL} -rpath 
                              ${build_libdir} ${install_libdir} 
                              ${install_libdir}/${XRD_LIB_PATH}
                      ERROR_QUIET
      )
    endforeach()
  else()
    message(WARNING "install_name_tool was not found. LC_RPATH variable will not be modified.")
  endif()

endfunction()
