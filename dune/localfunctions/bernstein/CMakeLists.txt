add_subdirectory(pk2d)

install(FILES
          pk.hh
          pk2d.hh
          pqkfactory.hh
        DESTINATION
          ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/bernstein)
