add_subdirectory(test)

install(FILES approximate_scattering.hh
              kanschat_scattering.hh
              periter.hh
              periter_uniform.hh
              subgridprojection.hh
              svdkernelapproximation.hh
              waveletkernelapproximation.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/dpg/radiative_transfer)
