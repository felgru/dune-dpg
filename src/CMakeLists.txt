add_executable("dune_dpg_error" dune_dpg_error.cc)
add_dune_suitesparse_flags(dune_dpg_error)
target_link_dune_default_libraries("dune_dpg_error")
add_dune_ug_flags(dune_dpg_error)

#add_executable("dune_dpg_saddlepoint" dune_dpg_saddlepoint.cc)
#add_dune_suitesparse_flags(dune_dpg_saddlepoint)
#target_link_dune_default_libraries("dune_dpg_saddlepoint")

#add_executable("visualization" visualization.cc)
#target_link_dune_default_libraries("visualization")

add_executable("dune_dpg_opt_testspace" dune_dpg_opt_testspace.cc)
add_dune_suitesparse_flags(dune_dpg_opt_testspace)
add_dune_ug_flags(dune_dpg_opt_testspace)
target_link_dune_default_libraries("dune_dpg_opt_testspace")

foreach(level RANGE 0 3)
  set(degree 2)
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

#add_executable("dune_dpg_conv_diff" dune_dpg_conv_diff.cc)
#add_dune_suitesparse_flags(dune_dpg_conv_diff)
#add_dune_ug_flags(dune_dpg_conv_diff)
#target_link_dune_default_libraries("dune_dpg_conv_diff")

#add_executable("dune_dpg_conv_diff_mod" dune_dpg_conv_diff_mod.cc)
#add_dune_suitesparse_flags(dune_dpg_conv_diff_mod)
#add_dune_ug_flags(dune_dpg_conv_diff_mod)
#target_link_dune_default_libraries("dune_dpg_conv_diff_mod")
