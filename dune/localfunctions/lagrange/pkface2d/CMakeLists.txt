set(pkface2ddir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange/pkface2d)

set(pkface2d_HEADERS  pkface2dlocalbasis.hh
                      pkface2dlocalcoefficients.hh
                      pkface2dlocalinterpolation.hh)

install(FILES ${pkface2d_HEADERS} DESTINATION ${pkface2ddir})
