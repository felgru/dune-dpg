dune_add_test(SOURCES test-alpertwavelets.cc
              LINK_LIBRARIES ${DUNE_LIBS}
              LABELS nompi)

dune_add_test(SOURCES test-haarwavelets.cc
              LINK_LIBRARIES ${DUNE_LIBS}
              LABELS nompi)

dune_add_test(SOURCES test-subgridprojection.cc
              LINK_LIBRARIES ${DUNE_LIBS}
              LABELS nompi)
