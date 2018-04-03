if(dune-subgrid_FOUND)
  dune_add_test(SOURCES hangingnodelagrangebasistest.cc
                LINK_LIBRARIES ${DUNE_LIBS})

  dune_add_test(SOURCES hangingnodebernsteinbasistest.cc
                LINK_LIBRARIES ${DUNE_LIBS})
endif()

dune_add_test(SOURCES normalizedbasistest.cc
              LINK_LIBRARIES ${DUNE_LIBS})

dune_add_test(SOURCES normalizedrefinedbasistest.cc
              LINK_LIBRARIES ${DUNE_LIBS})

dune_add_test(SOURCES pqkdgrefineddgbasistest.cc
              LINK_LIBRARIES ${DUNE_LIBS})

dune_add_test(SOURCES pqktracenodalbasistest.cc
              LINK_LIBRARIES dunecommon dunegeometry)
