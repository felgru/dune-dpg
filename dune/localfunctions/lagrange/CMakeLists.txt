set(SUBDIRS  qktrace pktrace2d)

set(lagrangedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange)
set(lagrange_HEADERS  qktrace.hh)

install(FILES ${lagrange_HEADERS} DESTINATION ${lagrangedir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt "EXCLUDE_FROM_ALL")
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
  unset(opt)
endforeach(i ${SUBDIRS})
