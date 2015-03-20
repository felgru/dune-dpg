set(SUBDIRS  optimaltest)

set(optimaltestfunctionsdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/optimaltestfunctions)
set(optimaltestfunctions_HEADERS  optimaltest.hh)

install(FILES ${optimaltestfunctions_HEADERS} DESTINATION ${optimaltestfunctionsdir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt "EXCLUDE_FROM_ALL")
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
  unset(opt)
endforeach(i ${SUBDIRS})
