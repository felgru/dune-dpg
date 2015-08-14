# which tests to run
set(TESTS test-qktrace test-qksubsampled)

# define the programs
add_executable(test-qktrace test-qktrace.cc)
target_link_libraries(test-qktrace ${DUNE_LIBS})
add_executable(test-qksubsampled test-qksubsampled.cc)
target_link_libraries(test-qksubsampled ${DUNE_LIBS})

foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})

# We do not want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
