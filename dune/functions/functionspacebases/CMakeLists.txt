add_subdirectory(test)

install(FILES
        bernsteinbasis.hh
        constrainedglobalbasis.hh
        constrainedlocalindexset.hh
        hangingnodep2nodalbasis.hh
        optimaltestbasis.hh
        pqkfacenodalbasis.hh
        pqksubsampleddgbasis.hh
        pqktracenodalbasis.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases)
