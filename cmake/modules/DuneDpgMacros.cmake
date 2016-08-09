# File for module specific CMake tests.
include(DuneBoost)
find_package(BoostFusion REQUIRED)
find_package(SuiteSparse COMPONENTS UMFPACK)
include(AddSuiteSparseFlags)
