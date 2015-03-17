add_subdirectory("test"  EXCLUDE_FROM_ALL)

install(FILES
        q2tracenodalbasis.hh
        qkdiscontinuousnodalbasis.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases)
