add_subdirectory(optimaltest)

install(FILES optimaltest.hh
        DESTINATION
          ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/optimaltestfunctions)
