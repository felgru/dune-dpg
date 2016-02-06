add_executable("plot_solution" plot_solution.cc)
add_dune_umfpack_flags(plot_solution)
add_dune_ug_flags(plot_solution)
target_link_dune_default_libraries("plot_solution")

foreach(level RANGE 0 3)
  set(degree 3)
  set(alevel 0)
  set(adegree 5)
  set(convergence_test
      "convergence_test_ls${level}_ks${degree}_la${alevel}_ka${adegree}")
  add_executable("${convergence_test}" convergence_test.cc)
  add_dune_umfpack_flags(${convergence_test})
  add_dune_ug_flags(${convergence_test})
  target_link_dune_default_libraries(${convergence_test})
  target_compile_definitions(${convergence_test}
                             PUBLIC LEVEL_SEARCH=${level};
                                    K_SEARCH=${degree};
                                    LEVEL_APOSTERIORI=${alevel};
                                    K_APOSTERIORI=${adegree})
endforeach()
