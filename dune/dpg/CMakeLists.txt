add_subdirectory(functions)
add_subdirectory(radiative_transfer)

add_subdirectory(test)

install(FILES applylocalfunctional_ru_impl.hh
              applylocalfunctional_uu_impl.hh
              assemble_helper.hh
              assemble_types.hh
              bilinearform.hh
              bilinearformfactory.hh
              boundarytools.hh
              cholesky.hh
              dpg_system_assembler.hh
              eigen_conversions.hh
              errorplotter.hh
              errortools.hh
              faces.hh
              functionplotter.hh
              innerproduct.hh
              innerproductfactory.hh
              integralterm.hh
              integralterm_rr_impl.hh
              integralterm_ru_impl.hh
              integralterm_uu_impl.hh
              leastsquares.hh
              linearform.hh
              linearformfactory.hh
              linearfunctionalterm.hh
              linearintegralterm.hh
              localcoefficients.hh
              localevaluation.hh
              locallinearterm_impl.hh
              quadrature.hh
              quadratureorder.hh
              refinedfaces.hh
              rhs_assembler.hh
              saddlepoint_system_assembler.hh
              spacetuple.hh
              testspace_coefficient_matrix.hh
              type_traits.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/dpg)
