cmake_minimum_required(VERSION 3.1)
project(dune-dpg CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
  ${dune-common_MODULE_PATH})

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

add_subdirectory("src")
add_subdirectory("dune")
add_subdirectory("doc")
add_subdirectory("cmake/modules")

# finalize the dune project, e.g. generating config.h etc.
if(dune-common_VERSION VERSION_GREATER_EQUAL 2.10)
  finalize_dune_project()
else()
  finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
endif()
