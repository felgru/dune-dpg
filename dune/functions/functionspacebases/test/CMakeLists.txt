if(dune-subgrid_FOUND)
  dune_add_test(SOURCES hangingnodelagrangebasistest.cc
                LABELS nompi)

  dune_add_test(SOURCES hangingnodebernsteinbasistest.cc
                LABELS nompi)
endif()

dune_add_test(SOURCES normalizedbasistest.cc
              LABELS nompi)

dune_add_test(SOURCES normalizedrefinedbasistest.cc)

dune_add_test(SOURCES lagrangedgrefineddgbasistest.cc)

dune_add_test(SOURCES lagrangetracebasistest.cc
              LABELS nompi)
