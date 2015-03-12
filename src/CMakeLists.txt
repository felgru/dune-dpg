add_executable("dune_dpg" dune_dpg.cc)
add_dune_umfpack_flags(dune_dpg)
target_link_dune_default_libraries("dune_dpg")

#add_executable("visualization" visualization.cc)
#target_link_dune_default_libraries("visualization")

#add_executable("functions-poisson" functions-poisson.cc)
#target_link_dune_default_libraries("functions-poisson")
