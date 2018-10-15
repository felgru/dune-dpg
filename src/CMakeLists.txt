if(dune-subgrid_FOUND)
  add_executable("dune_dpg_error" dune_dpg_error.cc)
  add_dune_suitesparse_flags(dune_dpg_error)
  target_link_dune_default_libraries("dune_dpg_error")
  add_dune_ug_flags(dune_dpg_error)
endif()

add_executable("dune_dpg_saddlepoint" dune_dpg_saddlepoint.cc)
add_dune_suitesparse_flags(dune_dpg_saddlepoint)
target_link_dune_default_libraries("dune_dpg_saddlepoint")

#add_executable("visualization" visualization.cc)
#target_link_dune_default_libraries("visualization")

add_executable("dune_dpg_opt_testspace" dune_dpg_opt_testspace.cc)
add_dune_suitesparse_flags(dune_dpg_opt_testspace)
add_dune_ug_flags(dune_dpg_opt_testspace)
target_link_dune_default_libraries("dune_dpg_opt_testspace")

if(dune-subgrid_FOUND)
  add_executable("plot_solution" plot_solution.cc)
  add_dune_suitesparse_flags(plot_solution)
  add_dune_ug_flags(plot_solution)
  target_link_dune_default_libraries("plot_solution")
endif()

add_executable("plot_solution_simple" plot_solution_simple.cc)
add_dune_suitesparse_flags(plot_solution_simple)
add_dune_ug_flags(plot_solution_simple)
target_link_dune_default_libraries("plot_solution_simple")

add_executable("manufactured_transport_uniform"
               manufactured_transport_uniform.cc)
add_dune_suitesparse_flags(manufactured_transport_uniform)
add_dune_ug_flags(manufactured_transport_uniform)
target_link_dune_default_libraries("manufactured_transport_uniform")

if(dune-subgrid_FOUND)
  add_executable("manufactured_transport" manufactured_transport.cc)
  add_dune_suitesparse_flags(manufactured_transport)
  add_dune_ug_flags(manufactured_transport)
  target_link_dune_default_libraries("manufactured_transport")
endif()

add_executable("plot_extension" plot_extension.cc)
add_dune_ug_flags(plot_extension)
target_link_dune_default_libraries("plot_extension")

add_executable("trial-to-test-condition" trial-to-test-condition.cc)
add_dune_ug_flags(trial-to-test-condition)
target_link_dune_default_libraries("trial-to-test-condition")

foreach(level 0 1 2)
  set(degree 3)
  set(profile_executable
      "profile_testspacecoefficientmatrix_ls${level}_ks${degree}")
  add_executable(${profile_executable}
                 profile_testspacecoefficientmatrix.cc)
  add_dune_suitesparse_flags(${profile_executable})
  add_dune_ug_flags(${profile_executable})
  target_link_dune_default_libraries(${profile_executable})
  target_compile_definitions(${profile_executable}
                             PUBLIC LEVEL_SEARCH=${level};
                                    K_SEARCH=${degree})
endforeach()

foreach(level RANGE 0 3)
  set(degree 3)
  set(alevel 0)
  set(adegree 5)
  set(convergence_test
      "convergence_test_ls${level}_ks${degree}_la${alevel}_ka${adegree}")
  add_executable("${convergence_test}" convergence_test.cc)
  add_dune_suitesparse_flags(${convergence_test})
  add_dune_ug_flags(${convergence_test})
  target_link_dune_default_libraries(${convergence_test})
  target_compile_definitions(${convergence_test}
                             PUBLIC LEVEL_SEARCH=${level};
                                    K_SEARCH=${degree};
                                    LEVEL_APOSTERIORI=${alevel};
                                    K_APOSTERIORI=${adegree})
endforeach()

if(dune-subgrid_FOUND)
  add_executable("asti_peaky" dune_dpg_rad_trans_asti.cc)
  add_dune_suitesparse_flags(asti_peaky)
  add_dune_ug_flags(asti_peaky)
  target_link_dune_default_libraries("asti_peaky")
  target_compile_definitions(asti_peaky
                             PUBLIC ASTI_PEAKY_BV=1)

  add_executable("asti_checkerboard" dune_dpg_rad_trans_asti.cc)
  add_dune_suitesparse_flags(asti_checkerboard)
  add_dune_ug_flags(asti_checkerboard)
  target_link_dune_default_libraries("asti_checkerboard")
  target_compile_definitions(asti_checkerboard
                             PUBLIC ASTI_CHECKERBOARD=1)

  add_executable("asti_peaky_normalized" dune_dpg_rad_trans_asti.cc)
  add_dune_suitesparse_flags(asti_peaky_normalized)
  add_dune_ug_flags(asti_peaky_normalized)
  target_link_dune_default_libraries("asti_peaky_normalized")
  target_compile_definitions(asti_peaky_normalized
                             PUBLIC ASTI_PEAKY_BV=1;
                                    ASTI_NORMALIZED_SPACES=1)

  add_executable("asti_checkerboard_normalized" dune_dpg_rad_trans_asti.cc)
  add_dune_suitesparse_flags(asti_checkerboard_normalized)
  add_dune_ug_flags(asti_checkerboard_normalized)
  target_link_dune_default_libraries("asti_checkerboard_normalized")
  target_compile_definitions(asti_checkerboard_normalized
                             PUBLIC ASTI_CHECKERBOARD=1;
                                    ASTI_NORMALIZED_SPACES=1)
endif()

add_executable("asti_peaky_uniform"
               dune_dpg_rad_trans_asti.cc)
add_dune_suitesparse_flags(asti_peaky_uniform)
add_dune_ug_flags(asti_peaky_uniform)
target_link_dune_default_libraries("asti_peaky_uniform")
target_compile_definitions(asti_peaky_uniform
                           PUBLIC ASTI_USE_UNIFORM_GRID=1;
                                  ASTI_PEAKY_BV=1)

add_executable("asti_checkerboard_uniform"
               dune_dpg_rad_trans_asti.cc)
add_dune_suitesparse_flags(asti_checkerboard_uniform)
add_dune_ug_flags(asti_checkerboard_uniform)
target_link_dune_default_libraries("asti_checkerboard_uniform")
target_compile_definitions(asti_checkerboard_uniform
                           PUBLIC ASTI_USE_UNIFORM_GRID=1;
                                  ASTI_CHECKERBOARD=1)

add_executable("kernel_svds" kernel_svds.cc)
target_link_dune_default_libraries("kernel_svds")

add_executable("kernel_truncations" kernel_truncations.cc)
target_link_dune_default_libraries("kernel_truncations")
