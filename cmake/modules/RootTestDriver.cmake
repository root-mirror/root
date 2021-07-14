# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

#-------------------------------------------------------------------------------
#
# RootCTestDriver.cmake
#
# CTest testdriver. Takes arguments via -DARG.
#
# Script arguments:
#
#   CMD   Command to be executed for the test.
#   PRE   Command to be executed before the test command.
#   POST  Command to be executed after the test command.
#   IN    File to be used as input
#   OUT   File to collect stdout and stderr.
#   ENV   Environment VAR1=Value1;VAR2=Value2.
#   CWD   Current working directory.
#   SYS   Value of ROOTSYS
#   DBG   Debug flag.
#   RC    Return code for success.
#
#-------------------------------------------------------------------------------

if(DBG)
  message(STATUS "ENV=${ENV}")
endif()

if(CMD)
  string(REPLACE "^" ";" _cmd ${CMD})
  if(DBG)
    message(STATUS "testdriver:CMD=${_cmd}")
  endif()
endif()

if(COPY)
  string(REPLACE "^" ";" _copy_files ${COPY})
  if(DBG)
    message(STATUS "files to copy: ${_copy_files}")
  endif()
endif()

if(PRE)
  string(REPLACE "^" ";" _pre ${PRE})
  if(DBG)
    message(STATUS "testdriver:PRE=${_pre}")
  endif()
endif()

if(POST)
  string(REPLACE "^" ";" _post ${POST})
  if(DBG)
    message(STATUS "testdriver:POST=${_post}")
  endif()
endif()

if(CWD)
  set(_cwd WORKING_DIRECTORY ${CWD})
  if(DBG)
    message(STATUS "testdriver:CWD=${CWD}")
  endif()
endif()

find_program(diff_cmd diff)

if(DIFFCMD)
  string(REPLACE "^" ";" diff_cmd ${DIFFCMD})
endif()

#---Set environment --------------------------------------------------------------------------------
if(ENV)
  string(REPLACE "#" ";" _env ${ENV})
  foreach(pair ${_env})
    string(REGEX REPLACE "^([^=]+)=(.*)$" "\\1;\\2" pair ${pair})
    list(GET pair 0 var)
    list(GET pair 1 val)
    set(ENV{${var}} ${val})
    if(DBG)
      message(STATUS "testdriver[ENV]:${var}==>${val}")
    endif()
  endforeach()
endif()

if(SYS)
  if(WIN32)
    file(TO_NATIVE_PATH ${SYS}/bin _path)
    set(ENV{PATH} "${_path};$ENV{PATH}")
  elseif(APPLE)
    set(ENV{PATH} ${SYS}/bin:$ENV{PATH})
    set(ENV{DYLD_LIBRARY_PATH} ${SYS}/lib:$ENV{DYLD_LIBRARY_PATH})
  else()
    set(ENV{PATH} ${SYS}/bin:$ENV{PATH})
    set(ENV{LD_LIBRARY_PATH} ${SYS}/lib:$ENV{LD_LIBRARY_PATH})
  endif()
endif()

#---Copy files to current direcotory----------------------------------------------------------------
if(COPY)
  foreach(copyfile ${_copy_files})
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${copyfile} ${CMAKE_CURRENT_BINARY_DIR}
                    RESULT_VARIABLE _rc)
    if(_rc)
      message(FATAL_ERROR "Copying file ${copyfile} to ${CMAKE_CURRENT_BINARY_DIR} failed! Error code : ${_rc}")
    endif()
  endforeach()
endif()

#---Execute pre-command-----------------------------------------------------------------------------
if(PRE)
  execute_process(COMMAND ${_pre} ${_cwd} RESULT_VARIABLE _rc)
  if(_rc)
    message(FATAL_ERROR "pre-command error code : ${_rc}")
  endif()
endif()

if(CMD)
  #---Execute the actual test ------------------------------------------------------------------------
  if(IN)
    set(_input INPUT_FILE ${IN})
  endif()

  if(OUT)
    # log stdout
    if(CHECKOUT)
      set(_chkout OUTPUT_VARIABLE _outvar)
    else()
      set(_chkout OUTPUT_VARIABLE _outvar2)
    endif()

    # log stderr
    if(ERRREF)
      set(_chkerr ERROR_VARIABLE _errvar)     # Check err reference
    else()
      if(CHECKERR AND CHECKOUT)
        set(_chkerr ERROR_VARIABLE _outvar)   # Both err and out together
      elseif(CHECKERR)
        set(_chkerr ERROR_VARIABLE _errvar0)  # Only check (no reference) the err and and ignore out
      else()
        set(_chkerr ERROR_VARIABLE _errvar2)  # Only check out eventually
      endif()
    endif()
 
    execute_process(COMMAND ${_cmd} ${_input} ${_chkout} ${_chkerr} WORKING_DIRECTORY ${CWD} RESULT_VARIABLE _rc)

    string(REGEX REPLACE "([.]*)[;][-][e][;]([^;]+)([.]*)" "\\1;-e '\\2\\3'" res "${_cmd}")
    string(REPLACE ";" " " res "${res}")
    message("-- TEST COMMAND -- ")
    message("cd ${CWD}")
    message("${res}")
    message("-- BEGIN TEST OUTPUT --")
    message("${_outvar}${_outvar2}")
    message("-- END TEST OUTPUT --")
    if(_errvar0 OR _errvar OR _errvar2)
      message("-- BEGIN TEST ERROR --")
      message("${_errvar0}${_errvar}${_errvar2}")
      message("-- END TEST ERROR --")
    endif()

    file(WRITE ${OUT} "${_outvar}")
    if(ERR)
      file(WRITE ${ERR} "${_errvar}")
    endif()

    if(_errvar0)
      # Filter messages in stderr that are expected
      string(STRIP "${_errvar0}" _errvar0)
      string(REPLACE "\n" ";" _lines "${_errvar0}")
      list(FILTER _lines EXCLUDE REGEX "^Info in <.+::ACLiC>: creating shared library.+")
      string(REPLACE ";" "\n" _errvar0 "${_lines}")
      if(_errvar0)
        message(FATAL_ERROR "Unexpected error output")
      endif()
    endif()

    if(DEFINED RC AND (NOT "${_rc}" STREQUAL "${RC}"))
      message(FATAL_ERROR "got exit code ${_rc} but expected ${RC}")
    elseif(NOT DEFINED RC AND _rc)
      message(FATAL_ERROR "got exit code ${_rc} but expected 0")
    endif()

    if(CNVCMD)
      string(REPLACE "^" ";" _outcnvcmd "${CNVCMD}^${OUT}")
      string(REPLACE "@" "=" _outcnvcmd "${_outcnvcmd}")
      execute_process(COMMAND ${_outcnvcmd} ${_chkout} ${_chkerr} RESULT_VARIABLE _rc)
      file(WRITE ${OUT} "${_outvar}")
      if(_rc)
        message(FATAL_ERROR "out conversion error code: ${_rc}")
      endif()
    endif()

    if(CNV)
      string(REPLACE "^" ";" _outcnv "sh;${CNV}")
      execute_process(COMMAND ${_outcnv} INPUT_FILE "${OUT}" OUTPUT_VARIABLE _outvar RESULT_VARIABLE _rc)
      file(WRITE ${OUT} "${_outvar}")
      if(_rc)
        message(FATAL_ERROR "out conversion error code: ${_rc}")
      endif()

      if(ERR)
        execute_process(COMMAND ${_outcnv} INPUT_FILE "${ERR}" OUTPUT_VARIABLE _errvar RESULT_VARIABLE _rc)
        file(WRITE ${ERR} "${_errvar}")
        if(_rc)
          message(FATAL_ERROR "err conversion error code: ${_rc}")
        endif()
      endif()
    endif()
  else()
    execute_process(COMMAND ${_cmd} ${_out} ${_err} ${_cwd} RESULT_VARIABLE _rc)

    if(_rc STREQUAL "Segmentation fault" OR _rc EQUAL 11)
       message(STATUS "Got ${_rc}, retrying again once more... with coredump enabled")
       string (REPLACE ";" " " _cmd_str "${_cmd}")
       file(WRITE run_with_coredump.sh "
pwd
ulimit -c unlimited
${_cmd_str}
")
       execute_process(COMMAND bash run_with_coredump.sh ${_out} ${_err} ${_cwd} RESULT_VARIABLE _rc)
    endif()

    if(DEFINED RC AND (NOT _rc EQUAL RC))
      message(FATAL_ERROR "error code: ${_rc}")
    elseif(NOT DEFINED RC AND _rc)
      message(FATAL_ERROR "error code: ${_rc}")
    endif()
  endif()

endif()

#---Execute post-command-----------------------------------------------------------------------------
if(POST)
  execute_process(COMMAND ${_post} ${_cwd} OUTPUT_VARIABLE _outvar ERROR_VARIABLE _outvar RESULT_VARIABLE _rc)
  if(_outvar)
    message("-- BEGIN POST OUTPUT --")
    message("${_outvar}")
    message("-- END POST OUTPUT --")
  endif()
  if(_rc)
    message(FATAL_ERROR "post-command error code : ${_rc}")
  endif()
endif()

if(OUTREF)
  execute_process(COMMAND ${diff_cmd} ${OUTREF} ${OUT} OUTPUT_VARIABLE _outvar ERROR_VARIABLE _outvar RESULT_VARIABLE _rc)
  if(_outvar)
    message("-- BEGIN OUTDIFF OUTPUT --")
    message("${_outvar}")
    message("-- END OUTDIFF OUTPUT --")
  endif()
  if(_rc)
    message(FATAL_ERROR "compare 'stdout' error: ${_rc}")
  endif()
endif()

if(ERRREF)
  execute_process(COMMAND ${diff_cmd} ${ERRREF} ${ERR} OUTPUT_VARIABLE _outvar ERROR_VARIABLE _outvar RESULT_VARIABLE _rc)
  if(_outvar)
    message("-- BEGIN ERRDIFF OUTPUT --")
    message("${_outvar}")
    message("-- END ERRDIFF OUTPUT --")
  endif()
  if(_rc)
    message(FATAL_ERROR "compare 'stderr' error: ${_rc}")
  endif()
endif()
