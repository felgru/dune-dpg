# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project does not adhere to [Semantic Versioning](http://semver.org/).

## 0.5 (Unreleased)
### Changed
* We now require a C++17-compatible compiler, e.g. GCC 7 or later.
* dune-dpg has been made compatible with version 2.8 of the Dune core
  modules.
* The `DoerflerMarking` function now uses `std::stable_sort` instead
  of `std::sort` to remove a source of potential non-determinism
  when there are elements with the same error estimates and only one
  of them gets marked.
* Fix last remaining compilation warnings given by GCC 8.
* Small optimization of `intersectSubGrids`: If the maxLevel of both
  subgrids is different, we only need to iterate until the smaller
  maxLevel.
* ASTI: we now have one canonical definition of the test space inner
  product. Previously we used different definitions depending on
  whether we used normalized or unnormalized spaces.
* Rename all PQk... spaces to Lagrange... .
  This adds the new spaces `LagrangeDGRefinedDGBasis`,
  `LagrangeSubsampledDGBasis`, `LagrangeDGSubsampledDGBasis`,
  `LagrangeFaceBasis` and `LagrangeTraceBasis`.
  The old PQk... names of these spaces remain as deprecated typedefs.
* Large parts of our code have been refactored for readability and
  to better adhere to modern C++ practices.

### Deprecated
* The spaces `PQkDGRefinedDGBasis`, `PQkSubsampledDGNodalBasis`,
  `PQkDGSubsampledDGNodalBasis`, `PQkFaceNodalBasis` and
  `PQkTraceNodalBasis` have been deprecated. All those spaces are
  now available under a name starting with Lagrange... .

### Removed
* `bindLocalIndexSets` and `getLocalIndexSets` have been removed,
  as IndexSets have been deprecated in dune-functions 2.6.

## 0.4.3 - 2020-02-26
### Added
* Added references to Felix’ thesis and our ASTI Math. Comp. paper
  to the README.

### Changed
* dune-dpg has been made compatible with the 2.7.0 release of the
  Dune core modules.
* All plotting scripts have been made compatible with Python 3.
* The dependency for dune-subgrid has been lowered to version 2.6 which
  should contain all the commits required by dune-dpg.
* When setting a custom output directory for the asti_* example
  programs, we don't append a timestamp to this directory anymore.
* The descriptions of example programs and how to run them have all
  been moved to doc/experiments.md.

### Fixed
* FindBoostHana.cmake now correctly uses the `BOOST_INCLUDEDIR`
  variable. It can thus find Boost::Hana now if you have a non-standard
  Boost include path, provided that you specify it with
  `-DBOOST_INCLUDEDIR=…` when invoking cmake.
* The `OptimalTestBasis` class has been fixed. It went unnoticed for
  a long time, that this class didn't compile as none of our test
  programs is using it anymore.

## 0.4.2 - 2018-11-21
### Added
* Added description of how to plot the condition number of the
  discrete trial to test space mapping to
  [experiments.md](doc/experiments.md).

### Fixed
* ASTI: use the same test inner product independent of whether we
  are using unnormalized or normalized spaces.
* ErrorPlotter: actually plot the error and not its square.
* Fixed plots that had the x label cut off at the bottom.

## 0.4 - 2018-10-18
### Added
* A simpler interface to create (bi)linear forms and inner products
  has been added. This tries to be less verbose than the old interface
  and allows to easily use both constant coefficients and ones that
  are `GridViewFunction`s. You can find examples of using the new
  interfaces in our example programs and in the commit notes of
  commit [a7104b565f35c735336fda764e95d68a97258b8d](https://gitlab.dune-project.org/felix.gruber/dune-dpg/commit/a7104b565f35c735336fda764e95d68a97258b8d).
* The radiative transfer code described in Felix’ thesis has been
  merged. You can use the example program asti_checkerboard_normalized
  to recreate the results in his thesis. See doc/experiments.md for
  more details.
* You can now normalize the basis functions in a GlobalBasis with
  respect to a given inner product by using `NormalizedBasis`
  or `NormalizedRefinedBasis`.
* `DPGSystemAssembler` gained a new method
  `applyHomogeneousDirichletBoundary` with a simpler interface than
  `applyDirichletBoundary`.
* A new `ConstrainedLocalView` was added that provides the indexing
  methods of `ConstrainedLocalIndexSet` like `DefaultLocalView` now
  provides indexing methods in dune-functions 2.6.
* FunctionPlotter can now be used to plot functions defined over a
  finite element space with local grid refinement. This is useful to
  plot data defined in the test search space.

### Changed
* dune-dpg now requires Version 2.6 of the Dune core modules.
  This allowed us to remove a lot of fallback code, making our code
  more readable and easier to maintain.
* You can now use `GridViewFunction`s as coefficients in integral terms.
  The quadrature order needed by a grid view function is queried by our
  matrix and vector assembly routines via the type trait
  `requiredQuadratureOrder<LocalFunction>::value`.
  If you write your own GridViewFunction, make sure to specialize
  `requiredQuadratureOrder` for its `LocalFunction`.
* Our own textbook implementation of the Cholesky decomposition has
  been replaced with the one implemented in Eigen which seems to be
  much more robust when the matrix is very badly conditioned.
* The `interpolate` method of all localfunctions can now be called with
  a function object f that implements `f(x)` instead of the
  `f.evaluate(x,y)` interface. This is in accordance with changes in
  dune-localfunctions 2.7.
* `iterateOverLocalIndexSet` has been renamed to `iterateOverLocalIndices`.
* The `BasisBuilder` namespace has been renamed to `BasisFactory` in
  accordance with a rename in dune-functions 2.6.
  (Since we never use the members of this namespace, this rename will
   probably not affect you.)

### Deprecated
* The `localIndexSet()` method of `ConstrainedGlobalBasis` has been
  marked as deprecated.

### Removed
* The workarounds for missing dune-subgrid functionality have been removed
  as everything is now implemented in SubGrid itself.
  Be sure to use at least commit 5706c70d8c4472f39e5bd61d1de1caf22381b8b9
  of dune-subgrid which implements all required methods in the SubGrid class.
* The version of `make_LinearFunctionalTerm` with domainOfIntegration
  template argument, that has been deprecated in 0.3, has been removed.
* The old `ErrorTools::DoerflerMarking` function which could choose between
  two residuals and which was deprecated in 0.3 has been removed.

## 0.3 - 2018-04-04
### Added
* New `ErrorTools::l2norm` function to compute the L_2 norm of a
  finite element function.
* New `ErrorPlotter` class to plot error estimators.
* New global bases `HangingNodeLagrangeP2Basis` and
  `HangingNodeBernsteinP2Basis` implementing P2 elements with hanging
  nodes. Since we only handle hanging nodes of order one, you have to
  make sure that neighboring cells in the grid do not vary in level by
  more than one. This can be assured by using
  [dune-subgrid](http://numerik.mi.fu-berlin.de/dune-subgrid/index.php)
  which by default limits the difference in level to maximally one.
  since dune-subgrid does not implement the whole grid interface, we
  defined missing functionality in `dune/dpg/subgrid_workarounds.hh`.
* Since the new hanging node bases do not fit the GlobalBasis interface
  as their index set does not implement DefaultIndexSet, we introduce a
  new `ConstrainedGlobalBasis` interface.
  Use the auxiliary functions `iterateOverLocalIndexSet` and
  `addToGlobalMatrix` to conveniently handle index sets from both
  interfaces.
* Add a Bernstein polynomial basis under the name
  `BernsteinPkLocalFiniteElement` with corresponding global bases
  `BernsteinBasis`, `BernsteinDGBasis`, `BernsteinDGRefinedDGBasis` and
  `HangingNodeBernsteinP2Basis`.
  They can be used as a more stable replacement for the standard
  Lagrange basis.
* New functions to save grid functions before adaptive refinement and to
  restore them to the refined grid, namely `attachDataToGrid` and
  `restoreDataToRefinedGrid`.
* New function `BoundaryTools::getBoundaryMask` that marks the boundary
  DoFs of a global basis.
  This works like `BoundaryTools::getInflowBoundaryMask` but marks the
  complete boundary instead of only the inflow boundary and might be
  useful to define the boundary of an elliptical problem.
* New `SkeletalLinearFunctionalTerm` that can be used to implement
  the second alternative for handling boundary values as given in
  Broersen, Dahmen, Stevenson Remark 3.6.
* New `SpaceTupleView` class to use part of a tuple of spaces as a
  space tuple pointer. This is used e.g. in our saddle point problem
  example to define test and trial spaces as part of the same tuple
  of spaces.
* New option to use a least squares solver instead of Cholesky by defining
  the preprocessor variable DUNE_DPG_USE_LEAST_SQUARES_INSTEAD_OF_CHOLESKY.
  This might help on strongly refined grids, where the Gramian in the trial
  to test map is so ill-conditioned that we cannot compute the Cholesky
  decomposition anymore. This has been observed when locally refining
  the unit square domain 19 times.
* `plot_solution.cc` now implements adaptive grid refinements. The old
  version on a fixed grid from v0.2.1 has been moved to
  `plot_solution_simple.cc`.

### Changed
* We now require version 2.5 of the DUNE core modules and a fully C++14
  compatible compiler, e.g. GCC-6.
  dune-dpg is fully forward-compatible to the upcoming 2.6 release of
  the DUNE core modules.
* `SystemAssembler`, `RhsAssembler`, `BilinearForm`, `LinearForm` and
  `InnerProduct` now try to share spaces. This allows the user to update
  all spaces at once after refining the grid. To make this sharing of
  spaces possible, the interfaces of those classes now take `shared_ptr`s
  to tuples of spaces. To create such `shared_ptr`s to tuples, you can use
  the new function `make_space_tuple`.
* the example programs plot_solution and dune_dpg_error now require
  dune-subgrid to compile as they use the new hanging node bases for the
  trace variables now.
* `BoundaryTools::getInflowBoundaryValue` has been renamed to
  `BoundaryTools::getBoundaryValue` as it also works for other boundary
  values. Additionaly, it now takes a single boundary value function
  instead of a tuple of functions. Since `getBoundaryValue` takes only a
  single FE space, it didn’t make any sense to pass multiple boundary
  value functions to it.
* Replace all usage of Boost::Fusion and MPL with C++14 metaprogramming.
  In many places we now use functions from dune/common/hybridutilities.hh
  or dune/common/tupleutility.hh instead. In the remaining places, where
  we needed to do some more complicated type computations, we now use
  Boost Hana. Thus our Boost dependency has been increased to Boost 1.63
  or newer.
* ErrorTools: The element-wise operations are now marked as private.
* ErrorTools: Mark all methods as static.
  In the future, we might make ErrorTools a namespace instead of a class.
* `ErrorTools` and `BoundaryTools` are not constructible anymore.
* `ErrorTools::DoerflerMarking` can now be called as
  `DoerflerMarking(grid, ratio, errorEstimates)` where `errorEstimates`
  is a `std::vector<std::tuple<EntitySeed, double>>&&` where each entry
  gives the seed of an element and its squared error estimate.
  You can use Dörfler marking with your own error estimators or use the
  function `ErrorTools::squaredCellwiseResidual` to compute the error
  estimates.
* `ErrorTools::computeL2error` now takes the exact solution directly
  instead of having it wrapped in a `std::tuple`.
* `LinearFunctionalTerm` now also works with a refined solution space.
* We now check in a static_assert that `defineCharacteristicFaces` is
  not called with refined spaces.
* Using DUNE 2.6 will use the more efficient `indices` interface which
  was recently added to dune-functions. This should give faster access
  to global indices, especially if you are using the `OptimalTestBasis`
  class.
* The boundaryValue argument of the `applyDirichletBoundaryToMatrix`
  method of DPGAssembler and SaddlepointAssembler has been removed as
  it was never used in this method.
* README, INSTALL: Mention dune-uggrid instead of ug.
* Remove some noise from the API documentation and add some more
  documentation for previously undocumented things.
* `ReferenceRefinementCache` now gets shared between refined spaces.
* INSTALL.md: Improve documentation of non-standard search paths.

### Fixed
* `BoundaryTools::getBoundaryValue` was mixing up local and global functions,
  thus giving wrong results. This has now been fixed.
* The computation of `unitOuterNormal` on faces of refined finite elements
  was wrong and has been fixed.
* Fix wrong indices when using `LinearFunctionalTerm` with more than one
  test space.
* Fix indices in `IntegralTerm` and `LinearFunctionalTerm` when both spaces
  are refined.
* Fix a bug in the computation of `splitPoint` used in the travel distance
  weighted inner product. This bug might have manifested in your code by
  throwing an error from the Cholesky solver that said that the matrix is
  not symmetrical positive definite.

### Deprecated
* The old `ErrorTools::DoerflerMarking` function which could choose between
  two residual has been deprecated. Instead use the new function described
  in the Changed section. There is now a function
  `ErrorTools::squaredCellwiseResidual(...)` that computes the squared
  residual from the old `DoerflerMarking` function. Thus you can replace
  your old function call with `ErrorTools::DoerflerMarking(grid, ratio,
  ErrorTools::squaredCellwiseResidual(...))`.

### Removed
* Remove the `make_DPGLinearForm`, `make_DPG_LinearForm` and
  `make_Saddlepoint_LinearForm` functions,
  use `make_LinearForm` instead.
* Remove the `refinedReferenceElement` method of refined global basis
  nodes. As a replacement, a new method `refinedReferenceElementGridView`
  has been added that directly returns the leafGridView.
* Remove the `update` method from `BufferedTestspaceCoefficientMatrix`
  and `UnbufferedTestspaceCoefficientMatrix`.
* Remove the `clone` method from all local bases. This method was never
  used and has been removed from all the local bases in
  dune-localfunctions 2.6.
* Make the auxiliary class `computeIndex` a private member.
  This should have never been a part of the public dune-dpg API.
* Remove `PQkTransportBasis` and `TransportQuadratureRule`.
  Our code did not work with this basis for a long time so obviously
  nobody was using it. This allowed us to simplify the internal
  `ChooseQuadrature` interface.
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
* Don’t try to install the removed header pqksubsamplednodalbasis.hh.
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
