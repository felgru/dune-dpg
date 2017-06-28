add_subdirectory(test)

install(FILES approximate_scattering.hh
              kanschat_scattering.hh
              periter_uniform.hh
              svdkernelapproximation.hh
              waveletkernelapproximation.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/dpg/radiative_transfer)
