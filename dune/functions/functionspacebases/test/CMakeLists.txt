dune_add_test(SOURCES pqktracenodalbasistest.cc
              LINK_LIBRARIES dunecommon dunegeometry)

if(dune-subgrid_FOUND)
  dune_add_test(SOURCES hangingnodebasistest.cc
                LINK_LIBRARIES ${DUNE_LIBS})
endif()
