set(qktrace_HEADERS
  qktracelocalbasis.hh
  qktracelocalcoefficients.hh
  qktracelocalinterpolation.hh)

install(FILES ${qktrace_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange/qktrace)
