add_subdirectory(io)

install(FILES
        analyticgridviewfunction.hh
        concepts.hh
        constraineddiscreteglobalbasisfunction.hh
        discreteglobalbasisfunction.hh
        gridviewfunctions.hh
        interpolate.hh
        localindexsetiteration.hh
        normalizedspaces.hh
        refinementinterpolation.hh
        subgridinterpolation.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/dpg/functions)
