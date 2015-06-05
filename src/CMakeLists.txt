add_executable("dune_dpg_error" dune_dpg_error.cc)
add_dune_umfpack_flags(dune_dpg_error)
target_link_dune_default_libraries("dune_dpg_error")
add_dune_ug_flags(dune_dpg_error)

#add_executable("dune_dpg_saddlepoint" dune_dpg_saddlepoint.cc)
#add_dune_umfpack_flags(dune_dpg_saddlepoint)
#target_link_dune_default_libraries("dune_dpg_saddlepoint")

#add_executable("visualization" visualization.cc)
#target_link_dune_default_libraries("visualization")

add_executable("dune_dpg_opt_testspace" dune_dpg_opt_testspace.cc)
add_dune_umfpack_flags(dune_dpg_opt_testspace)
add_dune_ug_flags(dune_dpg_opt_testspace)
target_link_dune_default_libraries("dune_dpg_opt_testspace")

#add_executable("System_Test" System_Test.cc)
#add_dune_ug_flags(System_Test)
#target_link_dune_default_libraries("System_Test")

add_executable("dune_dpg_conv_diff" dune_dpg_conv_diff.cc)
add_dune_umfpack_flags(dune_dpg_conv_diff)
add_dune_ug_flags(dune_dpg_conv_diff)
target_link_dune_default_libraries("dune_dpg_conv_diff")

add_executable("dune_dpg_conv_diff_mod" dune_dpg_conv_diff_mod.cc)
add_dune_umfpack_flags(dune_dpg_conv_diff_mod)
add_dune_ug_flags(dune_dpg_conv_diff_mod)
target_link_dune_default_libraries("dune_dpg_conv_diff_mod")
