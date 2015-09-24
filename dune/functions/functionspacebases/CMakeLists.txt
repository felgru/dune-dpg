add_subdirectory("test"  EXCLUDE_FROM_ALL)

install(FILES
        optimaltestbasis.hh
        pqkfacenodalbasis.hh
        pqksubsampleddgbasis.hh
        pqksubsamplednodalbasis.hh
        pqktracenodalbasis.hh
        pqktransportnodalbasis.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases)
