add_subdirectory(test)

install(FILES
        constrainedglobalbasis.hh
        constrainedlocalindexset.hh
        optimaltestbasis.hh
        pqkfacenodalbasis.hh
        pqksubsampleddgbasis.hh
        pqktracenodalbasis.hh
        pqktransportnodalbasis.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases)
