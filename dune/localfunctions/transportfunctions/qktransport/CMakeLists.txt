set(qktransport_HEADERS
  qktransportlocalbasis.hh
  qktransportlocalcoefficients.hh
  qktransportlocalinterpolation.hh)

install(FILES ${qktransport_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/transportfunctions/qktransport)
