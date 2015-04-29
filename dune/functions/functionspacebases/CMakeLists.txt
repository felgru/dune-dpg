add_subdirectory("test"  EXCLUDE_FROM_ALL)

install(FILES
        lagrangedgbasis.hh
        optimaltestbasis.hh
        pqkfacenodalbasis.hh
        pqktracenodalbasis.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases)
