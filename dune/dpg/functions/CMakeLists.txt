add_subdirectory(io)

install(FILES
        concepts.hh
        constraineddiscreteglobalbasisfunction.hh
        gridviewfunctions.hh
        interpolate.hh
        localindexsetiteration.hh
        normedspaces.hh
        refinementinterpolation.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/dpg/functions)
