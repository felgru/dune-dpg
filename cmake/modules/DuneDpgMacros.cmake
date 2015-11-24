# File for module specific CMake tests.
include(DuneBoost)
find_package(BoostFusion)
find_package(SuiteSparse COMPONENTS UMFPACK)
include(AddSuiteSparseFlags)
