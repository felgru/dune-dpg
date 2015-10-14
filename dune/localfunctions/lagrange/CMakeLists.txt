set(SUBDIRS  pkdgsubsampled2d pkface2d pktrace2d qkface2d qksubsampled qktrace)

set(lagrangedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange)
set(lagrange_HEADERS
    pkdgsubsampled2d.hh
    pkface2d.hh
    pktrace2d.hh
    pqkdgsubsampledfactory.hh
    pqkfacefactory.hh
    pqksubsampledfactory.hh
    pqktracefactory.hh
    qkface2d.hh
    qksubsampled.hh
    qktrace.hh)

install(FILES ${lagrange_HEADERS} DESTINATION ${lagrangedir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt "EXCLUDE_FROM_ALL")
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
  unset(opt)
endforeach(i ${SUBDIRS})
