# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project does not adhere to [Semantic Versioning](http://semver.org/).

## 0.3 (Unreleased)
### Added
* New `ErrorTools::l2norm` function to compute the L_2 norm of a
  finite element function.
* New GlobalBasis `HangingNodeP2NodalBasis` that implements P2 elements
  with hanging nodes. Since we only handle hanging nodes of order one,
  you have to make sure that neighboring cells in the grid do not vary
  in level by more than one. This can be assured by using
  [dune-subgrid](http://numerik.mi.fu-berlin.de/dune-subgrid/index.php)
  which by default limits the difference in level to maximally one.
  since dune-subgrid does not implement the whole grid interface, we
  defined missing functionality in `dune/dpg/subgrid_workarounds.hh`.
* Since `HangingNodeP2NodalBasis` does not fit the GlobalBasis interface
  as its index set does not implement DefaultIndexSet, we introduce a
  new `ConstrainedGlobalBasis` interface.
  Use the auxiliary functions `iterateOverLocalIndexSet` and
  `addToGlobalMatrix` to conveniently handle index sets from both
  interfaces.
* New functions to save grid functions before adaptive refinement and to
  restore them to the refined grid, namely `attachDataToGrid` and
  `restoreDataToRefinedGrid`.
* New function `BoundaryTools::getBoundaryMask` that marks the boundary
  DoFs of a global basis.
  This works like `BoundaryTools::getInflowBoundaryMask` but marks the
  complete boundary instead of only the inflow boundary and might be
  useful to define the boundary of an elliptical problem.

### Changed
* We now require version 2.5 of the DUNE core modules and a fully C++14
  compatible compiler, e.g. GCC-6.
* `SystemAssembler`, `RhsAssembler`, `BilinearForm`, `LinearForm` and
  `InnerProduct` now try to share spaces. This allows the user to update
  all spaces at once after refining the grid. To make this sharing of
  spaces possible, the interfaces of those classes now take `shared_ptr`s
  to tuples of spaces. To create such `shared_ptr`s to tuples, you can use
  the new function `make_space_tuple`.
* `BoundaryTools::getInflowBoundaryValue` has been renamed to
  `BoundaryTools::getBoundaryValue` as it also works for other boundary
  values. Additionaly, it now takes a single boundary value function
  instead of a tuple of functions. Since `getBoundaryValue` takes only a
  single FE space, it didn’t make any sense to pass multiple boundary
  value functions to it.
* Replace all usage of Boost::Fusion and MPI with C++14 metaprogramming.
  In many places we now use functions from dune/common/hybridutilities.hh
  or dune/common/tupleutility.hh instead. In the remaining places, where
  we needed to do some more complicated type computations, we now use
  Boost Hana. Thus our Boost dependency has been increased to Boost 1.63
  or newer.
* ErrorTools: The element-wise operations are now marked as private.
* ErrorTools: Mark all methods as static.
  In the future, we might make ErrorTools a namespace instead of a class.
* `ErrorTools::DoerflerMarking` can now be called as
  `DoerflerMarking(grid, ratio, errorEstimates)` where `errorEstimates`
  is a `std::vector<std::tuple<EntitySeed, double>>&&` so that you can use
  Dörfler marking with your own error estimators.
* README, INSTALL: Mention dune-uggrid instead of ug.
* Remove some noise from the API documentation.

### Fixed
* `BoundaryTools::getBoundaryValue` was mixing up local and global functions,
  thus giving wrong results. This has now been fixed.

### Deprecated
* The old `ErrorTools::DoerflerMarking` function which could choose between
  two residual has been deprecated. Instead use the new function described
  in the Changed section. There is now a `ErrorTools::residual(...)` function
  that computes the residual from the old `DoerflerMarking` function. Thus
  you can replace your old function call with
  `ErrorTools::DoerflerMarking(grid, ratio, ErrorTools::residual(...))`.

### Removed
* Remove the `make_DPGLinearForm` and `make_DPG_LinearForm` functions,
  use `make_LinearForm` instead.
* Remove some unused auxiliary types.

## 0.2.1 - 2016-11-07
### Fixed
* Lower Boost requirements to 1.56 and document that 1.60 and 1.61
  cause compilation errors.
* Fix compilation with GCC 6.
* Improve the INSTALL instructions:
    * Fix instructions for MacOS.
    * Fix URL for UG 3.12.1.
    * Describe how to use Boost 1.59 from a non-standard include path.

## 0.2 - 2016-09-18
### Added
* Add the possibility to cache the test space coefficient matrices
  with the newly introduced GeometryBuffer class when using a uniform
  grid and constant coefficients.
* Add LinearForm to better describe the rhs functions.
  This allows more easily to modify the rhs than the previous tuple
  of functions.
* Implement the travel distance weighted inner product from
  Broersen, Dahmen and Stevenson, 2015.
* Add a function to do adaptive grid refinements with Dörfler marking.
  This is not used in our examples yet.
  See function DoerflerMarking in dune/dpg/errortools.hh
* Add new aposteriori error indicator
  ||u-lifting(theta)||_L_2^2 + conforming residual(lifting(theta))^2
  used by Broersen, Dahmen and Stevenson.
* New replaceTestSpaces function to easily create BilinearForm,
  InnerProduct and LinearForm for computation of a posteriori errors.
* scripts/run_convergence_test.sh now takes a list of grid resolutions
  as argument.
* Add .editorconfig to define indentations, etc
  (see http://editorconfig.org for more information)

### Changed
* Split SystemAssembler into DPGSystemAssembler and
  SaddlepointSystemAssembler.
* Compute near optimal test spaces like in Roberts 2014.
  This gives a significant speedup and simplifies the interface of
  make_DPGSystemAssembler.
* Add new test program to profile the runtime difference between
  buffered and unbuffered test space coefficient matrices.
* CMake now complains if Boost Fusion is older than 1.59.
* Bump CMake dependency to 2.8.12 to be compatible with latest Dune
  requirements.
* Remove now unused template parameter FormulationType from BilinearForm.
* Add dummy implementation for partial() in local bases for newer
  versions of dune-localfunctions.
* Improve documentation in INSTALL and README.

## 0.1.2 (Unreleased)
### Added
* Completely rewritten build instructions can now be found in INSTALL.
* Describe the example programs more detailedly in README.
* plot_solution gained additional parameters c, betaX and betaY that
  are used to specify the parameters of the transport problem to solve.

### Changed
* Renamed variables in plot_solution to be compatible with the paper.
* Boost Fusion is now marked as a required dependency and will not appear
  under unfound packages in cmake.

### Fixed
* Don't try to install the removed header pqksubsamplednodalbasis.hh.
* Fixed remaining pedantic warnings (exept for the deprecation warnings
  generated by dune-grid which are out of our control).

## 0.1.1 - 2016-07-22
### Changed
* dune-dpg has been made compatible with Dune 2.4.1.
* The definition of BOOST_FUSION_INVOKE_PROCEDURE_MAX_ARITY has been
  moved to config.h.cmake.
* Consistently name tuples of spaces, localViews etc. with plural s
  to distinguish them from single objects.

### Fixed
* A bug in the computation of the right hand side for refined test spaces
  has been fixed. This also fixed the wrong results we were getting from
  the a posteriori estimator when using refined a posteriori test spaces.
* Added some missing #includes and removed unused ones.
* Fixed a lot of pedantic warnings like signed–unsigned comparison
  and unused variables.

## 0.1 - 2016-02-08
* initial release, compatible with Dune 2.4.0
