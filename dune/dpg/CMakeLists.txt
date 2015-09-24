#install headers
install(FILES bilinearform.hh
              cholesky.hh
              innerproduct.hh
              integralterm.hh
              assemble_helper.hh
              assemble_types.hh
              system_assembler.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/dpg)
