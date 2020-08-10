#---------------------------------------------------------------------------------------------------
# ROOT download test data files
#   Script arguments:
#     DST destination directory

if(NOT DST)
  set(DST .)
endif()

set(rootsite http://root.cern/files)
set(timeout 100)

if(NOT EXISTS ${DST}/h1)
  foreach(f dstarmb.root dstarp1a.root dstarp1b.root dstarp2.root)
    file(DOWNLOAD ${rootsite}/h1/${f} ${DST}/h1/${f} TIMEOUT ${timeout} SHOW_PROGRESS)
  endforeach()
else()
  message(STATUS "Already existing files in ${DST}/h1")
endif()

if(NOT EXISTS ${DST}/event)
  foreach(n 1 2 3 4 5 6 7 8 9 10)
    file(DOWNLOAD ${rootsite}/data/event_${n}.root ${DST}/event/event_${n}.root TIMEOUT ${timeout} SHOW_PROGRESS)
  endforeach()
else()
  message(STATUS "Already existing files in ${DST}/event")
endif()

