add_subdirectory(test)

install(FILES approximate_scattering.hh
              boundary_extension.hh
              kanschat_scattering.hh
              passkey.hh
              periter.hh
              periter_uniform.hh
              subgridprojection.hh
              subgrids.hh
              svdkernelapproximation.hh
              waveletkernelapproximation.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/dpg/radiative_transfer)
