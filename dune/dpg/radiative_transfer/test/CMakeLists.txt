dune_add_test(SOURCES test-alpertwavelets.cc
              LABELS nompi)

dune_add_test(SOURCES test-haarwavelets.cc
              LABELS nompi)

if(dune-subgrid_FOUND)
  dune_add_test(SOURCES test-subgridprojection.cc
                LABELS nompi)
endif()
