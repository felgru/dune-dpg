add_subdirectory(test)

install(FILES approximate_scattering.hh
              asti.hh
              asti_common.hh
              asti_uniform.hh
              boundary_extension.hh
              kanschat_scattering.hh
              passkey.hh
              subgridprojection.hh
              subgrids.hh
              waveletkernelapproximation.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/dpg/radiative_transfer)
