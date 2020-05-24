add_definitions(-DDUNE_DEPRECATED_INTERPOLATE_CHECK=1)

dune_add_test(SOURCES test-bernstein.cc
              LABELS nompi)

dune_add_test(SOURCES test-pkdgsubsampled.cc
              LABELS nompi)

dune_add_test(SOURCES test-qktrace.cc
              LABELS nompi)

dune_add_test(SOURCES test-qksubsampled.cc
              LABELS nompi)
