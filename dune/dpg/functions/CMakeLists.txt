add_subdirectory(io)

install(FILES
        concepts.hh
        constraineddiscreteglobalbasisfunction.hh
        discreteglobalbasisfunction.hh
        gridviewfunctions.hh
        interpolate.hh
        localindexsetiteration.hh
        normedspaces.hh
        refinementinterpolation.hh
        subgridinterpolation.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/dpg/functions)
