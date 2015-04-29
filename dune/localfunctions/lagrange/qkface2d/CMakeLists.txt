set(qkface2d_HEADERS
  qkface2dlocalbasis.hh
  qkface2dlocalcoefficients.hh
  qkface2dlocalinterpolation.hh)

install(FILES ${qkface2d_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange/qkface2d)
