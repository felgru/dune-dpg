# File for module specific CMake tests.
include(DuneBoost)
find_package(BoostFusion 1.56 REQUIRED)
find_package(SuiteSparse COMPONENTS UMFPACK)
include(AddSuiteSparseFlags)
