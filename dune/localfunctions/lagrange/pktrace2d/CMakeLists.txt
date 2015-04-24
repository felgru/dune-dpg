set(pktrace2ddir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange/pktrace2d)

set(pktrace2d_HEADERS  pktrace2dlocalbasis.hh
                       pktrace2dlocalcoefficients.hh
                       pktrace2dlocalinterpolation.hh)

install(FILES ${pktrace2d_HEADERS} DESTINATION ${pktrace2ddir})
