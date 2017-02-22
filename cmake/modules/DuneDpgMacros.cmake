# File for module specific CMake tests.
include(DuneBoost)
find_package(BoostHana 1.61 REQUIRED)
find_package(SuiteSparse COMPONENTS UMFPACK)
include(AddSuiteSparseFlags)
