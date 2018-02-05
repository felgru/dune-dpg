add_subdirectory(normalizedbasisadaptor)
add_subdirectory(test)

install(FILES
        bernsteinbasis.hh
        bernsteindgbasis.hh
        bernsteindgrefineddgnodalbasis.hh
        constrainedglobalbasis.hh
        constrainedlocalindexset.hh
        hangingnodep2nodalbasis.hh
        normalizedbasisadaptor.hh
        optimaltestbasis.hh
        pqkdgrefineddgnodalbasis.hh
        pqksubsampleddgbasis.hh
        pqkfacenodalbasis.hh
        pqksubsampleddgbasis.hh
        pqktracenodalbasis.hh
        referencerefinementfactory.hh
        refinednode.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases)
