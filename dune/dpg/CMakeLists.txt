add_subdirectory(common)

install(FILES assemble_helper.hh
              assemble_types.hh
              bilinearform.hh
              boundarytools.hh
              cholesky.hh
              errortools.hh
              innerproduct.hh
              integralterm.hh
              quadrature.hh
              system_assembler.hh
              type_traits.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/dpg)
