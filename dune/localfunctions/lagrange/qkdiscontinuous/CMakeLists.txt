set(qk_HEADERS
  qkdiscontinuouslocalbasis.hh
  qkdiscontinuouslocalcoefficients.hh
  qkdiscontinuouslocalinterpolation.hh)

install(FILES ${qk_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange/qk)
