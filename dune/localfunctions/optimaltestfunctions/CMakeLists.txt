add_subdirectory(optimaltest refinedoptimaltest)

install(FILES optimaltest.hh
              refinedoptimaltest.hh
        DESTINATION
          ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/optimaltestfunctions)
