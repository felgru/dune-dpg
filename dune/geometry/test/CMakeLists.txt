set(TESTS test-subsampledquadrature)

add_executable(test-subsampledquadrature test-subsampledquadrature.cc)
target_link_libraries(test-subsampledquadrature dunegeometry ${DUNE_LIBS})
add_test(test-subsampledquadrature test-subsampledquadrature)

# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} test-subsampledquadrature)
