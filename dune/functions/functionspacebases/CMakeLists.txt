add_subdirectory(normalizedbasisadaptor)
add_subdirectory(test)

install(FILES
        bernsteinbasis.hh
        bernsteindgbasis.hh
        bernsteindgrefineddgnodalbasis.hh
        constrainedglobalbasis.hh
        constrainedlocalview.hh
        hangingnodebernsteinp2basis.hh
        hangingnodelagrangep2basis.hh
        lagrangedgrefineddgbasis.hh
        lagrangedgsubsampleddgbasis.hh
        lagrangefacebasis.hh
        lagrangesubsampleddgbasis.hh
        lagrangetracebasis.hh
        normalizedbasisadaptor.hh
        normalizedrefinedbasisadaptor.hh
        optimaltestbasis.hh
        referencerefinementfactory.hh
        refinedglobalbasis.hh
        refinedlocalview.hh
        refinednode.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases)
