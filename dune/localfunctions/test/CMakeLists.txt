dune_add_test(SOURCES test-bernstein.cc
              LINK_LIBRARIES ${DUNE_LIBS}
              LABELS nompi)

dune_add_test(SOURCES test-pkdgsubsampled.cc
              LINK_LIBRARIES ${DUNE_LIBS}
              LABELS nompi)

dune_add_test(SOURCES test-qktrace.cc
              LINK_LIBRARIES ${DUNE_LIBS}
              LABELS nompi)

dune_add_test(SOURCES test-qksubsampled.cc
              LINK_LIBRARIES ${DUNE_LIBS}
              LABELS nompi)
