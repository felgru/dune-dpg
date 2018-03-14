if(dune-subgrid_FOUND)
  dune_add_test(SOURCES hangingnodebasistest.cc
                LINK_LIBRARIES ${DUNE_LIBS})
endif()

dune_add_test(SOURCES normalizedrefinedbasistest.cc
              LINK_LIBRARIES ${DUNE_LIBS})

dune_add_test(SOURCES pqkdgrefineddgbasistest.cc
              LINK_LIBRARIES ${DUNE_LIBS})

dune_add_test(SOURCES pqktracenodalbasistest.cc
              LINK_LIBRARIES dunecommon dunegeometry)
