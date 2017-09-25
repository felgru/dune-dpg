#add_executable("dune_dpg_error" dune_dpg_error.cc)
#add_dune_suitesparse_flags(dune_dpg_error)
#target_link_dune_default_libraries("dune_dpg_error")
#add_dune_ug_flags(dune_dpg_error)

#add_executable("dune_dpg_saddlepoint" dune_dpg_saddlepoint.cc)
#add_dune_suitesparse_flags(dune_dpg_saddlepoint)
#target_link_dune_default_libraries("dune_dpg_saddlepoint")

#add_executable("visualization" visualization.cc)
#target_link_dune_default_libraries("visualization")

#add_executable("dune_dpg_opt_testspace" dune_dpg_opt_testspace.cc)
#add_dune_suitesparse_flags(dune_dpg_opt_testspace)
#add_dune_ug_flags(dune_dpg_opt_testspace)
#target_link_dune_default_libraries("dune_dpg_opt_testspace")

#add_executable("plot_solution" plot_solution.cc)
#add_dune_suitesparse_flags(plot_solution)
#add_dune_ug_flags(plot_solution)
#target_link_dune_default_libraries("plot_solution")

#add_executable("plot_solution_simple" plot_solution_simple.cc)
#add_dune_suitesparse_flags(plot_solution_simple)
#add_dune_ug_flags(plot_solution_simple)
#target_link_dune_default_libraries("plot_solution_simple")

add_executable("plot_extension" plot_extension.cc)
add_dune_ug_flags(plot_extension)
target_link_dune_default_libraries("plot_extension")

#foreach(level 0 1 2)
#  set(degree 3)
#  set(profile_executable
#      "profile_testspacecoefficientmatrix_ls${level}_ks${degree}")
#  add_executable(${profile_executable}
#                 profile_testspacecoefficientmatrix.cc)
#  add_dune_suitesparse_flags(${profile_executable})
#  add_dune_ug_flags(${profile_executable})
#  target_link_dune_default_libraries(${profile_executable})
#  target_compile_definitions(${profile_executable}
#                             PUBLIC LEVEL_SEARCH=${level};
#                                    K_SEARCH=${degree})
#endforeach()

#foreach(level RANGE 0 3)
#  set(degree 3)
#  set(alevel 0)
#  set(adegree 5)
#  set(convergence_test
#      "convergence_test_ls${level}_ks${degree}_la${alevel}_ka${adegree}")
#  add_executable("${convergence_test}" convergence_test.cc)
#  add_dune_suitesparse_flags(${convergence_test})
#  add_dune_ug_flags(${convergence_test})
#  target_link_dune_default_libraries(${convergence_test})
#  target_compile_definitions(${convergence_test}
#                             PUBLIC LEVEL_SEARCH=${level};
#                                    K_SEARCH=${degree};
#                                    LEVEL_APOSTERIORI=${alevel};
#                                    K_APOSTERIORI=${adegree})
#endforeach()

add_executable("dune_dpg_rad_trans_periter" dune_dpg_rad_trans_periter.cc)
add_dune_suitesparse_flags(dune_dpg_rad_trans_periter)
add_dune_ug_flags(dune_dpg_rad_trans_periter)
target_link_dune_default_libraries("dune_dpg_rad_trans_periter")

add_executable("dune_dpg_rad_trans_periter_uniform"
               dune_dpg_rad_trans_periter.cc)
add_dune_suitesparse_flags(dune_dpg_rad_trans_periter_uniform)
add_dune_ug_flags(dune_dpg_rad_trans_periter_uniform)
target_link_dune_default_libraries("dune_dpg_rad_trans_periter_uniform")
target_compile_definitions(dune_dpg_rad_trans_periter_uniform
                           PUBLIC PERITER_USE_UNIFORM_GRID=1)

add_executable("kernel_svds" kernel_svds.cc)
target_link_dune_default_libraries("kernel_svds")
