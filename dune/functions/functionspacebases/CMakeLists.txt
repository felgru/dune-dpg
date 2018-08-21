add_subdirectory(normalizedbasisadaptor)
add_subdirectory(test)

install(FILES
        bernsteinbasis.hh
        bernsteindgbasis.hh
        bernsteindgrefineddgnodalbasis.hh
        constrainedglobalbasis.hh
        constrainedlocalindexset.hh
        constrainedlocalview.hh
        hangingnodebernsteinp2basis.hh
        hangingnodelagrangep2basis.hh
        normalizedbasisadaptor.hh
        normalizedrefinedbasisadaptor.hh
        optimaltestbasis.hh
        pqkdgrefineddgnodalbasis.hh
        pqksubsampleddgbasis.hh
        pqkfacenodalbasis.hh
        pqksubsampleddgbasis.hh
        pqktracenodalbasis.hh
        referencerefinementfactory.hh
        refinedglobalbasis.hh
        refinedlocalview.hh
        refinednode.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases)
