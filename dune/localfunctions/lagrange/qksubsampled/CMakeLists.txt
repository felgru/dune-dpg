set(qksubsampled_HEADERS
  qksubsampledlocalbasis.hh
  qksubsampledlocalcoefficients.hh
  qksubsampledlocalinterpolation.hh)

install(FILES ${qksubsampled_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange/qksubsampled)
