set(optimaltest_HEADERS
  optimaltestlocalbasis.hh
  optimaltestlocalcoefficients.hh
  optimaltestlocalinterpolation.hh)

install(FILES ${optimaltest_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/optimaltestfunctions/optimaltest)
