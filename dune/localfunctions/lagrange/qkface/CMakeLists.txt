set(qkface_HEADERS
  qkfacelocalbasis.hh
  qkfacelocalcoefficients.hh
  qkfacelocalinterpolation.hh)

install(FILES ${qkface_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange/qkface)
