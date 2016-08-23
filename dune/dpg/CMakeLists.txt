add_subdirectory(functions)

install(FILES applylocalfunctional_ru_impl.hh
              applylocalfunctional_uu_impl.hh
              assemble_helper.hh
              assemble_types.hh
              bilinearform.hh
              boundarytools.hh
              cholesky.hh
              dpg_system_assembler.hh
              errortools.hh
              getvolumeterm_impl.hh
              innerproduct.hh
              integralterm.hh
              integralterm_rr_impl.hh
              integralterm_ru_impl.hh
              integralterm_uu_impl.hh
              linearform.hh
              linearfunctionalterm.hh
              linearintegralterm.hh
              localevaluation.hh
              quadrature.hh
              rhs_assembler.hh
              saddlepoint_system_assembler.hh
              testspace_coefficient_matrix.hh
              type_traits.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/dpg)
