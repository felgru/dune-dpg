add_subdirectory(qktransport)

install(FILES  pqktransportfactory.hh
               qktransport.hh
        DESTINATION
        ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/transportfunctionsdir)
