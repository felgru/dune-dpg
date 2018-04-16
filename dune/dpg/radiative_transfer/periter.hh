// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_PERITER_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_PERITER_HH

#include <chrono>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <dune/common/std/optional.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/io.hh>
#include <dune/istl/umfpack.hh>

#include <dune/functions/functionspacebases/hangingnodebernsteinp2basis.hh>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/errortools.hh>
#include <dune/dpg/functionplotter.hh>
#include <dune/dpg/functions/interpolate.hh>
#include <dune/dpg/functions/subgridinterpolation.hh>
#include <dune/dpg/functions/refinementinterpolation.hh>
#include <dune/dpg/linearfunctionalterm.hh>
#include <dune/dpg/radiative_transfer/approximate_scattering.hh>
#include <dune/dpg/radiative_transfer/boundary_extension.hh>
#include <dune/dpg/radiative_transfer/subgridprojection.hh>
#include <dune/dpg/rhs_assembler.hh>
#include <dune/dpg/dpg_system_assembler.hh>
#include <dune/dpg/type_traits.hh>

#include <dune/dpg/radiative_transfer/periter_common.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <dune/subgrid/subgrid.hh>
#pragma GCC diagnostic pop

#include <boost/math/constants/constants.hpp>

namespace Dune {

template<class SubGrid>
std::set<typename SubGrid::HostGridType::GlobalIdSet::IdType>
saveSubGridToIdSet(const SubGrid& subGrid)
{
  auto& idSet = subGrid.getHostGrid().globalIdSet();
  auto subGridView = subGrid.leafGridView();
  std::set<typename std::decay_t<decltype(idSet)>::IdType> subGridElements;
  for(const auto& e : elements(subGridView)) {
    subGridElements.insert(idSet.id(subGrid.template getHostEntity<0>(e)));
  }
  return subGridElements;
}

template<class SubGrid>
std::unique_ptr<SubGrid>
restoreSubGridFromIdSet(
    std::set<typename SubGrid::HostGridType::GlobalIdSet::IdType>&& idSet,
    typename SubGrid::HostGridType& hostGrid)
{
  auto subGrid = std::make_unique<SubGrid>(hostGrid);
  subGrid->createBegin();
  subGrid->insertSet(idSet);
  subGrid->createEnd();
  subGrid->setMaxLevelDifference(1);
  return subGrid;
}

template<class SubGrid>
std::unique_ptr<SubGrid>
restoreSubGridFromIdSet(
    const std::set<typename SubGrid::HostGridType::GlobalIdSet::IdType>& idSet,
    typename SubGrid::HostGridType& hostGrid)
{
  std::set<typename SubGrid::HostGridType::GlobalIdSet::IdType>
    idSetCopy(idSet);
  return restoreSubGridFromIdSet<SubGrid>(std::move(idSetCopy), hostGrid);
}

template<class SubGrid>
std::unique_ptr<SubGrid> copySubGrid(const SubGrid& subGrid) {
  std::unique_ptr<SubGrid> gr
      = std::make_unique<SubGrid>(subGrid.getHostGrid());
  gr->createBegin();
  const auto gridView = subGrid.leafGridView();
  for(const auto& e : elements(gridView)) {
    gr->insert(subGrid.template getHostEntity<0>(e));
  }
  gr->createEnd();
  gr->setMaxLevelDifference(1);
  return gr;
}

template<class SubGrid, class HostGrid>
std::unique_ptr<SubGrid> fullSubGrid(HostGrid& hostGrid) {
  std::unique_ptr<SubGrid> gr
      = std::make_unique<SubGrid>(hostGrid);
  gr->createBegin();
  const auto gridView = hostGrid.leafGridView();
  for(const auto& e : elements(gridView)) {
    gr->insert(e);
  }
  gr->createEnd();
  gr->setMaxLevelDifference(1);
  return gr;
}

template<class GridView>
class SubGridSpaces;

/**
 * This class describes the Periter algorithm for radiative transfer problems
 *
 * \tparam ScatteringKernelApproximation
 *         specifies the method used to approximate the scattering kernel
 * \tparam RHSApproximation  if right hand side and lifting of boundary
 *                           values are finite element functions, set this
 *                           to FeRHS, otherwise set this to ApproximateRHS
 */
template<class ScatteringKernelApproximation, class RHSApproximation>
class Periter {
  public:

  using Direction = typename ScatteringKernelApproximation::Direction;

  /**
   * Solve a radiative transfer problem using the Periter algorithm
   *
   * \param hostGrid
   * \param f  right hand side function
   * \param g  lifting of the boundary values
   * \param is_inflow_boundary_homogeneous
   *            checks if g is 0 on the inflow boundary
   * \param sigma   absorbtion coefficient
   * \param kernel  the scattering kernel, e.g. a Henyey–Greenstein kernel
   * \param rho  the contraction parameter ρ
   * \param CT  an upper bound for the norm of the transport solver
   * \param cB  the inf-sup constant of the operator B = T - K
   * \param targetAccuracy  periter solves up to this accuracy
   * \param maxNumberOfIterations  ... or up to the given number of iterations
   *                               (whatever comes first)
   * \param plotSolutions  specifies when to create .vtu files for plotting
   *                       the solution
   * \note By default, the scattering is computed from the inner degrees
   *       of freedom of the SolutionSpace. To use the skeletal degrees of
   *       freedom instead, define the preprocessor variable
   *       PERITER_SKELETAL_SCATTERING before including periter.hh.
   */
  template<class HostGrid, class F, class G, class HB,
           class Sigma, class Kernel>
  void solve(HostGrid& hostGrid,
             const F& f,
             const G& g,
             const HB& is_inflow_boundary_homogeneous,
             const Sigma sigma,
             const Kernel& kernel,
             double rho,
             double CT,
             double cB,
             double targetAccuracy,
             unsigned int maxNumberOfIterations,
             unsigned int maxNumberOfInnerIterations,
             const std::string& outputfolder,
             PlotSolutions plotSolutions = PlotSolutions::doNotPlot);

  private:
  using VectorType = BlockVector<FieldVector<double,1>>;

  /**
   * Compute the solution of a transport problem
   *
   * This solves a transport problem for a given direction s on a fixed
   * grid. Afterwards it computes a posteriori error estimates and does
   * Dörfler marking.
   *
   * This computes the modified problem from Broersen, Dahmen, Stevenson,
   * Remark 3.6, alternative 2 for non-homogeneous problems. In particular
   * that means that the trace part of the solution will not fit the
   * non-homogeneous boundary values while the inner part does.
   *
   * \param[out] x the solution of the transport problem
   * \param grid the grid needed for Döfler marking
   * \param spaces
   * \param s the transport direction
   * \param sigma the absorbtion coefficient
   * \param rhsData the data for the unmodified rhs
   * \param boundary_is_homogeneous
   * \param bvExtension a lifting of the boundary values
   *
   * \return the squared a posteriori error of the solution
   */
  template<class Spaces, class Grid, class Sigma, class RHSData>
  static double compute_transport_solution(
      VectorType& x,
      Grid& grid,
      const Spaces& spaces,
      const FieldVector<double, 2>& s,
      const Sigma sigma,
      RHSData& rhsData,
      bool boundary_is_homogeneous,
      const VectorType& bvExtension);

  /**
   * Apply the scattering integral to a solution x
   *
   * This corresponds to [K, u_n, κ_1 * η] in the Periter algorithm
   * (see Dahmen, Gruber, Mula).
   *
   * \param kernelApproximation an approximation to the scattering kernel
   * \param x  solution to which we want to apply the scattering kernel
   * \param subGridSpaces
   * \param accuracy
   */
  template<class SubGridView, class HostGrid, class Grid, class GridIdSet>
  static std::vector<VectorType> apply_scattering(
      ScatteringKernelApproximation& kernelApproximation,
      const std::vector<VectorType>& x,
      SubGridSpaces<SubGridView>& subGridSpaces,
      const HostGrid& hostGrid,
      std::vector<Direction>& sVector,
      std::vector<std::unique_ptr<Grid>>& grids,
      std::vector<GridIdSet>& gridIdSets,
      double accuracy);

  /**
   * create grids for new directions created in apply_scattering
   */
  template<class Grid>
  static void create_new_grids(
      std::vector<std::unique_ptr<Grid>>& grids,
      size_t numNewGrids);

  /**
   * insert new gridIdSets in apply_scattering after adding new grids
   */
  template<class GridIdSet, class Grid>
  static void save_grids_to_gridIdSets(
      std::vector<GridIdSet>& gridIdSets,
      const std::vector<std::unique_ptr<Grid>>& grids);
};

template<class GridView>
class SubGridSpaces {
  public:
  using Spaces = TransportSpaces<GridView,
                    Functions::HangingNodeBernsteinP2Basis<GridView>>;

  using FEBasisInterior     = typename Spaces::FEBasisInterior;
  using FEBasisTrace        = typename Spaces::FEBasisTrace;

  using FEBasisTest         = typename Spaces::FEBasisTest;
  using FEBasisEnrichedTest = typename Spaces::FEBasisEnrichedTest;

  using SolutionSpacePtr = typename Spaces::SolutionSpacePtr;
  using TestSpacePtr = typename Spaces::TestSpacePtr;
  using EnrichedTestSpacePtr = typename Spaces::EnrichedTestSpacePtr;

  using GridsVector = std::vector<std::unique_ptr<typename GridView::Grid>>;

  SubGridSpaces(const GridsVector& grids)
  {
    spaces_.reserve(grids.size());
    for(const auto& grid : grids)
    {
      spaces_.emplace_back(grid->leafGridView());
    }
  }

  void update(size_t i, const GridView& gridView) {
    spaces_[i].update(gridView);
  }

  /**
   * insert new spaces in apply_scattering after adding new grids
   */
  void create_new_spaces(const GridsVector& grids) {
    if(spaces_.size() != grids.size()) {
      spaces_.clear();
      spaces_.reserve(grids.size());
      for(const auto& grid : grids) {
        spaces_.emplace_back(grid->leafGridView());
      }
    } else {
      for(size_t i = 0, iend = grids.size(); i < iend; i++) {
        update(i, grids[i]->leafGridView());
      }
    }
  }

  const FEBasisInterior& interiorSolutionSpace(size_t i) const {
    return spaces_[i].interiorSolutionSpace();
  }

  const FEBasisTrace& traceSolutionSpace(size_t i) const {
    return spaces_[i].traceSolutionSpace();
  }

  const FEBasisTest& testSpace(size_t i) const {
    return spaces_[i].testSpace();
  }

  const Spaces& spaces(size_t i) const {
    return spaces_[i];
  }

  private:

  std::vector<Spaces> spaces_;
};

#ifndef DOXYGEN
namespace detail {
  constexpr size_t dim = 2;
  using VectorType = BlockVector<FieldVector<double,1>>;
  using Direction = FieldVector<double, dim>;

  template<class FEHostBasis, class Grids, class F>
  inline void approximate_rhs (
      std::vector<VectorType>& rhsFunctional,
      Grids&,
      double,
      const FEHostBasis& hostGridBasis,
      const std::vector<Direction>& sVector,
      F& f,
      FeRHS) {
    static_assert(!is_RefinedFiniteElement<FEHostBasis>::value,
        "Functions::interpolate won't work for refined finite elements");
    const size_t numS = sVector.size();
    for(unsigned int i = 0; i < numS; ++i)
    {
      VectorType gInterpolation(hostGridBasis.size());
      const Direction s = sVector[i];
      Functions::interpolate(hostGridBasis, gInterpolation,
          [s,&f](const Direction& x) { return f(x,s); });

      // Add gInterpolate to first hostGridBasis.size() entries of
      // rhsFunctional[i].
      using Iterator = std::decay_t<decltype(rhsFunctional[i].begin())>;
      for(Iterator rIt=rhsFunctional[i].begin(),
                   rEnd=rhsFunctional[i].begin()+hostGridBasis.size(),
                   gIt=gInterpolation.begin(); rIt!=rEnd; ++rIt, ++gIt) {
        *rIt += *gIt;
      }
    }
  }


  // Refine grid until accuracy kappa2*eta is reached for
  // approximation of rhs.
  template<class FEBases, class Grids, class F>
  inline void approximate_rhs (
      std::vector<VectorType>& rhsFunctional,
      Grids& grids,
      double accuracy,
      const std::vector<std::shared_ptr<FEBases>>& solutionSpaces,
      const std::vector<Direction>& sVector,
      F& f,
      ApproximateRHS) {
    DUNE_THROW(Dune::NotImplemented,
        "Implementation of approximate_rhs for non-FE right hand side "
        "functions currently broken.");
  }
} // end namespace detail
#endif

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class HostGrid, class F, class G, class HB,
         class Sigma, class Kernel>
void Periter<ScatteringKernelApproximation, RHSApproximation>::solve(
           HostGrid& hostGrid,
           const F& f,
           const G& g,
           const HB& is_inflow_boundary_homogeneous,
           const Sigma sigma,
           const Kernel& kernel,
           double rho,
           double CT,
           double cB,
           double targetAccuracy,
           unsigned int maxNumberOfIterations,
           unsigned int maxNumberOfInnerIterations,
           const std::string& outputfolder,
           PlotSolutions plotSolutions) {
  if((plotSolutions & PlotSolutions::plotLastIteration)
      == PlotSolutions::plotLastIteration) {
    std::cerr
        << "Plotting of only the last iteration is not implemented yet!\n";
    std::exit(1);
  }
  static_assert(std::is_same<RHSApproximation, FeRHS>::value
      || std::is_same<RHSApproximation, ApproximateRHS>::value,
      "Unknown type provided for RHSApproximation!\n"
      "Should be either FeRHS or ApproximateRHS.");

  const unsigned int dim = HostGrid::dimension;
  using Grid = SubGrid<dim, HostGrid, false>;
  using LeafGridView = typename Grid::LeafGridView;
  using HostGridView = typename HostGrid::LeafGridView;
  using Direction = FieldVector<double, dim>;
  using GridIdSet = std::set<typename HostGrid::GlobalIdSet::IdType>;

  /////////////////////////////////////////////
  // To print information in dune-dpg/results/
  /////////////////////////////////////////////

  std::ofstream ofs(outputfolder+"/output");

  ///////////////////////////////////
  // Parameters for adaptivity
  ///////////////////////////////////

  // TODO: estimate norm of rhs f in V'
  // Remark: Here, V=H_{0,+}(D\times S)
  const double fnorm = 1;
  const double err0 = fnorm / cB;

  detail::ApproximationParameters approximationParameters(rho, CT, err0,
                                                          RHSApproximation{});

  ////////////////////////////////////////////
  // Handle directions of discrete ordinates
  ////////////////////////////////////////////

  ScatteringKernelApproximation kernelApproximation(kernel,
      approximationParameters.finalScatteringAccuracy(targetAccuracy));

  ofs << "PERITER algorithm\n"
      << "=================\n"
      << "Prescribed final accuracy: "
      << targetAccuracy << '\n'
      << kernel.info()  << '\n'
      << "Wavelet order: "
      << kernelApproximation.getWltOrder() << '\n'
      << kernelApproximation.typeApprox()  << '\n'
      << "Maximum wavelet level: "
      << kernelApproximation.getMaxLevel() << '\n'
      << "Maximum number of directions: "
      << kernelApproximation.maxNumS()     << '\n'
      << "Periter parameters:" << '\n'
      << approximationParameters
      << "CT = "     << CT     << '\n';

  if(kernelApproximation.typeApprox() == "Kernel approximation with: SVD"){
    std::vector<double> singularValues
      = kernelApproximation.getSingularValues();
    ofs << "Singular values of kernel matrix:\n";
    for(size_t i=0; i<singularValues.size(); i++){
      ofs << singularValues[i] << '\n';
    }
  }

  // As the solution u we use for the initial scattering is 0, and the
  // formula for the accuracy contains a 1/\|u\|, we set the initial
  // accuracy to a large enough value.
  std::vector<Direction>
    sVector(kernelApproximation.setAccuracyAndInputSize(1e5, 0));
  unsigned int numS = sVector.size();

  std::vector<std::unique_ptr<Grid>> grids;
  grids.reserve(sVector.size());
  for(size_t i = 0, imax = sVector.size(); i < imax; i++) {
    std::unique_ptr<Grid> gr = std::make_unique<Grid>(hostGrid);
    gr->createBegin();
    gr->insertLevel(hostGrid.maxLevel());
    gr->createEnd();
    gr->setMaxLevelDifference(1);
    grids.emplace_back(std::move(gr));
  }

  std::vector<GridIdSet> gridIdSets;
  save_grids_to_gridIdSets(gridIdSets, grids);

  //////////////////////////////////
  //   right hand side vector type
  //////////////////////////////////
  typedef BlockVector<FieldVector<double,1> > VectorType;

  //////////////////////////////////
  //   Solution vectors
  //////////////////////////////////
  std::vector<VectorType> x(numS);
  using Spaces = SubGridSpaces<LeafGridView>;
  Spaces spaces(grids);

  for(unsigned int i = 0; i < grids.size(); ++i)
  {
    x[i].resize(spaces.interiorSolutionSpace(i).size()
                 + spaces.traceSolutionSpace(i).size());
    x[i] = 0;
  }

  /////////////////////////////////////////////////////////
  //  Fixed-point iterations
  /////////////////////////////////////////////////////////

  double accuracy = err0;

  std::vector<double> aposterioriIter(maxNumberOfIterations,0.);

  for(unsigned int n = 0; accuracy > targetAccuracy
                          && n < maxNumberOfIterations; ++n)
  {
    ofs << "\nIteration n=" << n
        << "\n================\n";
    std::cout << "\nIteration " << n << "\n\n";

    for(size_t i = 0; i < grids.size(); ++i) {
      grids[i] = restoreSubGridFromIdSet<Grid>(gridIdSets[i],
                                               hostGrid);
      spaces.update(i, grids[i]->leafGridView());
    }

    std::chrono::steady_clock::time_point startScatteringApproximation
        = std::chrono::steady_clock::now();

    const std::vector<double> quadWeight
      = kernelApproximation.getQuadWeightSubinterval();

    double kappaNorm = 1.;
    double uNorm = 0.;
    for(size_t i=0; i<grids.size(); ++i) {
      const double uiNorm =
        ErrorTools::l2norm(spaces.interiorSolutionSpace(i), x[i]);
      uNorm += uiNorm * uiNorm
                / (1 << kernelApproximation.getLevel())
                * quadWeight[i % kernelApproximation.numSperInterval];
    }
    uNorm = std::sqrt(uNorm);
    // To prevent division by zero.
    if(uNorm == 0.) uNorm = 1e-5;

#ifdef PERITER_SKELETAL_SCATTERING
    using FEBasisHostTraceDiscontinuous
        = Functions::BernsteinDGBasis<HostGridView, 2>;
#else
    using FEBasisHostInterior
        = changeGridView_t<typename Spaces::FEBasisInterior, HostGridView>;
#endif
    using RHSData = std::decay_t<decltype(attachDataToSubGrid(
              std::declval<typename Spaces::FEBasisTest>(),
#ifdef PERITER_SKELETAL_SCATTERING
              std::declval<FEBasisHostTraceDiscontinuous>(),
#else
              std::declval<FEBasisHostInterior>(),
#endif
              std::declval<VectorType>()))>;
    std::vector<RHSData> rhsData;
    using FEBasisHostTrace
        = changeGridView_t<typename Spaces::FEBasisTrace, HostGridView>;
    using BVData = std::decay_t<decltype(attachDataToSubGrid(
              std::declval<typename Spaces::FEBasisTest>(),
              std::declval<FEBasisHostTrace>(),
              std::declval<VectorType>()))>;
    std::vector<Std::optional<BVData>> bvData;
    const double accuKernel = approximationParameters.scatteringAccuracy()
                            / (kappaNorm * uNorm);
    {
#ifdef PERITER_SKELETAL_SCATTERING
      FEBasisHostTraceDiscontinuous
          hostGridGlobalBasis(hostGrid.leafGridView());
#else
      FEBasisHostInterior hostGridGlobalBasis(hostGrid.leafGridView());
#endif

      std::vector<VectorType> rhsFunctional =
          apply_scattering (
            kernelApproximation, x, spaces,
            hostGridGlobalBasis,
            sVector, grids, gridIdSets,
            accuKernel);

      numS = sVector.size();
      x.resize(numS);
      rhsData.reserve(numS);

      if((plotSolutions & PlotSolutions::plotScattering)
          == PlotSolutions::plotScattering) {
        std::cout << "Plot scattering:\n";

        for(unsigned int i = 0; i < numS; ++i)
        {
          std::cout << "Direction " << i << '\n';

          std::string name = outputfolder
                          + std::string("/scattering_n")
                          + std::to_string(n)
                          + std::string("_s")
                          + std::to_string(i);
          FunctionPlotter scatteringPlotter(name);
          scatteringPlotter.plot("scattering", rhsFunctional[i],
                                 hostGridGlobalBasis, 0);
        }
      }

      // TODO: restore grids from gridIdSets and update spaces
      //       shouldn't be necessary, as long as RHSApproximation == FeRHS
      detail::approximate_rhs (
          rhsFunctional,
          grids,
          approximationParameters.rhsAccuracy(),
          hostGridGlobalBasis,
          sVector,
          f, RHSApproximation{});
      for(size_t i = 0; i < numS; i++) {
        const auto& subGridGlobalBasis = spaces.testSpace(i);
        rhsData.emplace_back(
          attachDataToSubGrid(
            subGridGlobalBasis,
            hostGridGlobalBasis,
            rhsFunctional[i]));
      }
    }
    std::vector<bool> boundary_is_homogeneous(numS, false);
    {
      for(size_t i = 0; i < numS; i++) {
        boundary_is_homogeneous[i]
            = is_inflow_boundary_homogeneous(sVector[i]);
        // TODO: write a generic test for homogeneous inflow boundary
      }
      // get bv contribution to rhs
      FEBasisHostTrace hostGridGlobalBasis(hostGrid.leafGridView());
      const VectorType boundaryExtension
          = harmonic_extension_of_boundary_values(g, hostGridGlobalBasis);
      bvData.reserve(numS);
      for(size_t i = 0; i < numS; i++) {
        if(boundary_is_homogeneous[i]) {
          bvData.emplace_back();
        } else {
          const auto& subGridGlobalBasis = spaces.testSpace(i);
          bvData.emplace_back(
            attachDataToSubGrid(
              subGridGlobalBasis,
              hostGridGlobalBasis,
              boundaryExtension));
        }
      }
    }

    std::chrono::steady_clock::time_point endScatteringApproximation
        = std::chrono::steady_clock::now();

    ofs << "eta_n = rhobar^{-n}: " << approximationParameters.eta() << '\n'
        << "\n--------------------\n"
        << "Info angular approx:\n"
        << "--------------------\n"
        << "Current wavelet level: "
        << kernelApproximation.getLevel() << '\n'
        << "Number of directions: "
        << kernelApproximation.getNumS()  << '\n'
        << "Directions are:\n";
    for(size_t i = 0; i<numS; i++){
      ofs << sVector[i] << '\n';
    }
    ofs << "\n---------------------\n"
        << "Kernel approximation:\n"
        << "---------------------\n"
        << "Accuracy required: "
          << approximationParameters.scatteringAccuracy() << " (kappa1 * eta)\n"
        << "Accuracy introduced in code: "
        << accuKernel << " (kappa1 * eta / (kappaNorm * uNorm))\n"
        << kernelApproximation.info()      << '\n'
        << "Computing time: "
          << std::chrono::duration_cast<std::chrono::microseconds>
          (endScatteringApproximation - startScatteringApproximation).count()
          << "us\n" << std::flush;


    ////////////////////////////////////////////////////
    // Inner loop
    ////////////////////////////////////////////////////
    ofs << "\n-----------------------------------\n"
             "Inner iterations (transport solves)\n"
             "-----------------------------------\n";
    double aposterioriTransportGlobal = 0.;

    // Loop over the spatial grids.
    // For each angular direction there is one spatial grid.
    for(unsigned int i = 0; i < grids.size(); ++i)
    {
      ofs << "\nAngle " << i
          << "\n--------\n";

      grids[i] = restoreSubGridFromIdSet<Grid>(gridIdSets[i],
                                               hostGrid);
      spaces.update(i, grids[i]->leafGridView());

      double aposteriori_s;
      unsigned int nRefinement = 0;
      for( ; ; )
        // At the end of the loop, we will break if
        // aposteriori_s < kapp3*eta (pointwise in s)
        // or ++nRefinement >= maxNumberOfInnerIterations
        // thus the inner loop terminates eventually.
      {
        VectorType bvExtension;
        if(!boundary_is_homogeneous[i]) {
          const auto& feBasisTest = spaces.testSpace(i);
          auto newGridData
              = bvData[i]->restoreDataToRefinedSubGrid(feBasisTest);
          bvExtension.resize(newGridData.size(),
                             false /* don't copy old values */);
          for(size_t k = 0, kmax = newGridData.size(); k < kmax; k++) {
            bvExtension[k] = newGridData[k];
          }
        }

        {
          std::cout << "Direction " << i
                    << ", inner iteration " << nRefinement << '\n';

          aposteriori_s
              = compute_transport_solution(x[i], *grids[i],
                  spaces.spaces(i),
                  sVector[i], sigma, rhsData[i], boundary_is_homogeneous[i],
                  bvExtension);

          aposteriori_s = std::sqrt(aposteriori_s);

          // TODO: Add (interpolation of) g to theta part of x[i]?

          ofs << "Iteration " << n << '.' << nRefinement
              << " for direction " << i << ": \n"
              << "  - A posteriori estimation of || (u,trace u) - (u_fem,theta) || = "
              << aposteriori_s
              << "\n  - Grid level: "     << grids[i]->maxLevel()
              << "\n  - Number of DOFs: " << x[i].size()
              << std::endl;

          std::cout << "\nIteration " << n << '.' << nRefinement
              << " for direction " << i << ": \n"
              << "  - A posteriori estimation of || (u,trace u) - (u_fem,theta) || = "
              << aposteriori_s
              << "\n  - Grid level: " << grids[i]->maxLevel()
              << "\n  - Number of DOFs: " << x[i].size()
              << std::endl;
        }


        if(++nRefinement >= maxNumberOfInnerIterations
            || aposteriori_s <= approximationParameters.transportAccuracy()) {
          gridIdSets[i] = saveSubGridToIdSet(*grids[i]);

          ofs << "\na posteriori error for current direction: "
              << aposteriori_s;
          if (aposteriori_s <= approximationParameters.transportAccuracy())
          {
            ofs << " (enough";
          } else {
            ofs << " (not enough, required "
                << approximationParameters.transportAccuracy();
          }
          ofs << ")\n\n" << std::flush;

          break;
        } else {
          grids[i]->preAdapt();
          grids[i]->adapt();
          grids[i]->postAdapt();
          spaces.update(i, grids[i]->leafGridView());

          ofs << "\na posteriori error for current direction: "
              << aposteriori_s
              << " (required "
              << approximationParameters.transportAccuracy()
              << ")\n\n" << std::flush;
        }
      } // end of spatial refinement for angular direction i
      aposterioriTransportGlobal = std::max(aposterioriTransportGlobal,
                                            aposteriori_s);

      if((plotSolutions & PlotSolutions::plotOuterIterations)
          == PlotSolutions::plotOuterIterations) {
        //////////////////////////////////////////////////////////////////////
        //  Write result to VTK file
        //  We need to subsample, because VTK cannot natively display
        //  real second-order functions
        //////////////////////////////////////////////////////////////////////

        const auto& feBasisInterior = spaces.interiorSolutionSpace(i);
        const auto& feBasisTrace    = spaces.traceSolutionSpace(i);

        std::cout << "Plot solution for direction " << i << '\n';

        std::string name = outputfolder
                          + std::string("/u_rad_trans_n")
                          + std::to_string(n)
                          + std::string("_s")
                          // + std::to_string(i*stride);
                          + std::to_string(i);
        FunctionPlotter uPlotter(name);
        uPlotter.plot("u", x[i], feBasisInterior, 0, 0);
        name = outputfolder
                          + std::string("/theta_rad_trans_n")
                          + std::to_string(n)
                          + std::string("_s")
                          // + std::to_string(i*stride);
                          + std::to_string(i);
        FunctionPlotter thetaPlotter(name);
        thetaPlotter.plot("theta", x[i], feBasisTrace, 2,
                          feBasisInterior.size());
      }
    }

    const size_t accumulatedDoFs = std::accumulate(x.cbegin(), x.cend(),
        static_cast<size_t>(0),
        [](size_t acc, auto vec) { return acc + vec.size(); });

    // A posteriori estimation of error ||bar u_n -T^{-1}K bar u_{n-1}||
    aposterioriIter[n] = aposterioriTransportGlobal
                       + CT * approximationParameters.scatteringAccuracy();

    // Error bound for || u - \bar u_n || based on a priori errors
    accuracy = approximationParameters.combinedAccuracy();
    // Error bound for || u_n - \bar u_n || based on a posteriori errors
    double errorAPosteriori
        = approximationParameters.aPosterioriError(aposterioriIter);

    ofs << "---------------------\n"
           "End inner iterations \n"
           "---------------------\n"
           "Error transport solves (a posteriori estimation): "
          << aposterioriTransportGlobal                  << '\n'
        << "Accuracy kernel: "
          << approximationParameters.scatteringAccuracy() << '\n'
        << "Error bound ||bar u_n -T^{-1}K bar u_{n-1}|| (a posteriori): "
          << aposterioriIter[n]   << '\n'
        << "Error bound ||u_n - bar u_n|| (a posteriori): "
          << errorAPosteriori << '\n'
        << "Bound global accuracy ||u - bar u_n|| (a priori + a posteriori): "
          << accuracy
          << " (rho * err0 + 2) * rho^n\n"
        << "Total number of DoFs: "
          << accumulatedDoFs
        << "\n\n" << std::flush;

    std::cout << "Error at end of Iteration " << n << ": "
              << aposterioriTransportGlobal << ", using "
              << accumulatedDoFs << " DoFs\n";

    approximationParameters.decreaseEta();
  }
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class Spaces, class Grid, class Sigma, class RHSData>
double
Periter<ScatteringKernelApproximation, RHSApproximation>::
compute_transport_solution(
    VectorType& x,
    Grid& grid,
    const Spaces& spaces,
    const FieldVector<double, 2>& s,
    const Sigma sigma,
    RHSData& rhsData,
    bool boundary_is_homogeneous,
    const VectorType& bvExtension)
{
  auto bilinearForm =
    make_BilinearForm(spaces.testSpacePtr(), spaces.solutionSpacePtr(),
        make_tuple(
            make_IntegralTerm<0,0,IntegrationType::valueValue,
                                  DomainOfIntegration::interior>(sigma),
            make_IntegralTerm<0,0,
                              IntegrationType::gradValue,
                              DomainOfIntegration::interior>(-1., s),
            make_IntegralTerm<0,1,IntegrationType::normalVector,
                                  DomainOfIntegration::face>(1., s)));
  auto bilinearFormEnriched =
      replaceTestSpaces(bilinearForm, spaces.enrichedTestSpacePtr());
  auto innerProduct =
    make_InnerProduct(spaces.testSpacePtr(),
        make_tuple(
            make_IntegralTerm<0,0,IntegrationType::gradGrad,
                                  DomainOfIntegration::interior>(1., s),
            make_IntegralTerm<0,0,
                              IntegrationType::travelDistanceWeighted,
                              DomainOfIntegration::face>(1., s)));
  auto innerProductEnriched =
      replaceTestSpaces(innerProduct, spaces.enrichedTestSpacePtr());


  using LeafGridView = typename Grid::LeafGridView;
  static_assert(std::is_same<typename Spaces::GridView, LeafGridView>::value,
                "GridViews don't match!");
  using Geometry = typename LeafGridView::template Codim<0>::Geometry;
  using GeometryBuffer_t = GeometryBuffer<Geometry>;

  GeometryBuffer_t geometryBuffer;
  auto systemAssembler =
      make_DPGSystemAssembler(bilinearForm,
                              innerProduct,
                              geometryBuffer);

  VectorType rhsFunctional;
  {
    auto& feBasisTest = spaces.testSpace();
    auto newGridData
        = rhsData.restoreDataToRefinedSubGrid(feBasisTest);
    rhsFunctional.resize(newGridData.size(),
                         false /* don't copy old values */);
    for(size_t k = 0, kmax = newGridData.size(); k < kmax; k++) {
      rhsFunctional[k] = newGridData[k];
    }
  }

  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

  VectorType rhs;
  MatrixType stiffnessMatrix;

  /////////////////////////////////////////////////////////
  //  Assemble the systems
  /////////////////////////////////////////////////////////
  if(boundary_is_homogeneous) {
    auto rhsFunction = make_LinearForm(
          systemAssembler.getTestSearchSpaces(),
          std::make_tuple(
            make_LinearFunctionalTerm<0>
              (rhsFunctional, spaces.testSpace())));
    systemAssembler.assembleSystem(
        stiffnessMatrix, rhs,
        rhsFunction);
  } else {
    assert(bvExtension.size() > 0);
    auto rhsFunction = make_LinearForm(
          systemAssembler.getTestSearchSpaces(),
          std::make_tuple(
            make_LinearFunctionalTerm<0>
              (rhsFunctional, spaces.testSpace()),
            make_SkeletalLinearFunctionalTerm
              <0, IntegrationType::normalVector>
              (bvExtension, spaces.testSpace(), -1, s)));
    systemAssembler.assembleSystem(
        stiffnessMatrix, rhs,
        rhsFunction);
  }
  {
    // Determine Dirichlet dofs for theta (inflow boundary)
    std::vector<bool> dirichletNodesInflow;
    BoundaryTools::getInflowBoundaryMask(
                spaces.traceSolutionSpace(),
                dirichletNodesInflow,
                s);
    systemAssembler.template applyDirichletBoundary<1>
        (stiffnessMatrix,
        rhs,
        dirichletNodesInflow,
        0.);
  }
#if 0
  systemAssembler.template defineCharacteristicFaces<1>(
      stiffnessMatrix,
      rhs, s);
#endif

  ////////////////////////////////////
  //   Initialize solution vector
  ////////////////////////////////////
  x.resize(spaces.interiorSolutionSpace().size()
           + spaces.traceSolutionSpace().size());
  x = 0;

  ////////////////////////////
  //   Compute solution
  ////////////////////////////

  std::cout << "rhs size = " << rhs.size()
            << " matrix size = " << stiffnessMatrix.N()
                        << " x " << stiffnessMatrix.M()
            << " solution size = " << x.size() << '\n';


  const int verbosity = 0; // 0: not verbose; >0: verbose
  UMFPack<MatrixType> umfPack(stiffnessMatrix, verbosity);
  InverseOperatorResult statistics;
  umfPack.apply(x, rhs, statistics);

  ////////////////////////////////////
  //  A posteriori error
  ////////////////////////////////////
  std::cout << "Compute a posteriori error\n";

  // We compute the a posteriori error
  // - We compute the rhs with the enriched test space ("rhs=f(v)")
  // -- Contribution of the source term f that has an
  //    analytic expression
  // -- Contribution of the scattering term
  if(boundary_is_homogeneous) {
    auto rhsAssemblerEnriched
      = make_RhsAssembler(spaces.enrichedTestSpacePtr());
    auto rhsFunction = make_LinearForm(
          rhsAssemblerEnriched.getTestSpaces(),
          std::make_tuple(
            make_LinearFunctionalTerm<0>
              (rhsFunctional,
               std::get<0>(*systemAssembler.getTestSearchSpaces()))));
    rhsAssemblerEnriched.assembleRhs(rhs, rhsFunction);
  } else {
    auto rhsAssemblerEnriched
      = make_RhsAssembler(spaces.enrichedTestSpacePtr());
    auto rhsFunction = make_LinearForm(
          rhsAssemblerEnriched.getTestSpaces(),
          std::make_tuple(
            make_LinearFunctionalTerm<0>
              (rhsFunctional,
               std::get<0>(*systemAssembler.getTestSearchSpaces())),
            make_SkeletalLinearFunctionalTerm
              <0, IntegrationType::normalVector>
              (bvExtension,
               std::get<0>(*systemAssembler.getTestSearchSpaces()),
               -1, s)));
    rhsAssemblerEnriched.assembleRhs(rhs, rhsFunction);
  }
  // - Computation of the a posteriori error
  // - Doerfler marking
  const double ratio = .6;
  const double aposteriori_s =
      ErrorTools::DoerflerMarking(grid, ratio,
            ErrorTools::squaredCellwiseResidual(bilinearFormEnriched,
                                                innerProductEnriched,
                                                x, rhs));

  return aposteriori_s;
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class SubGridView, class HostGridBasis, class Grid, class GridIdSet>
std::vector<typename Periter<ScatteringKernelApproximation,
                             RHSApproximation>::VectorType>
Periter<ScatteringKernelApproximation, RHSApproximation>::apply_scattering(
      ScatteringKernelApproximation& kernelApproximation,
      const std::vector<VectorType>& x,
      SubGridSpaces<SubGridView>& subGridSpaces,
      const HostGridBasis& hostGridBasis,
      std::vector<Direction>& sVector,
      std::vector<std::unique_ptr<Grid>>& grids,
      std::vector<GridIdSet>& gridIdSets,
      double accuracy) {
  sVector = kernelApproximation.setAccuracyAndInputSize(accuracy, x.size());

  using Spaces = SubGridSpaces<SubGridView>;
#ifdef PERITER_SKELETAL_SCATTERING
  using FEBasisTrace = typename Spaces::FEBasisTrace;
#else
  using FEBasisInterior = typename Spaces::FEBasisInterior;
#endif

  // Interpolate x[i] to hostGridBasis.
  std::vector<VectorType> xHost(x.size());
  for(size_t i = 0, imax = x.size(); i < imax; i++) {
#ifdef PERITER_SKELETAL_SCATTERING
    const FEBasisTrace& feBasisTrace = subGridSpaces.traceSolutionSpace(i);
    const size_t  feBasisTraceOffset = subGridSpaces.interiorSolutionSpace(i)
                                                    .size();
    interpolateFromSubGrid(
        feBasisTrace, x[i].begin() + feBasisTraceOffset,
        hostGridBasis, xHost[i]);
#else
    const FEBasisInterior& feBasisInterior
        = subGridSpaces.interiorSolutionSpace(i);
    interpolateFromSubGrid(
        feBasisInterior, x[i],
        hostGridBasis, xHost[i]);
#endif
  }

  const auto scatteringAssembler =
      make_ApproximateScatteringAssembler(hostGridBasis,
                                          kernelApproximation);
  const size_t numS = sVector.size();
  std::vector<VectorType> rhsFunctional(numS);
  scatteringAssembler.computeScattering(rhsFunctional, xHost);

  create_new_grids(grids, numS);
  save_grids_to_gridIdSets(gridIdSets, grids);
  subGridSpaces.create_new_spaces(grids);

  return rhsFunctional;
}

/**
 * compute intersection of two SubGrids
 */
template<class SubGrid>
std::unique_ptr<SubGrid> intersectSubGrids(const SubGrid& subGrid1,
                                           const SubGrid& subGrid2)
{
  // assert(subGrid1.getHostGrid() == subGrid2.getHostGrid());
  std::unique_ptr<SubGrid> gr
      = std::make_unique<SubGrid>(subGrid1.getHostGrid());
  gr->createBegin();
  for(int level=0; level <= subGrid1.maxLevel(); ++level)
  {
    const auto levelGridView = subGrid1.levelGridView(level);
    for (auto&& e : elements(levelGridView))
    {
      const auto eHost = subGrid1.template getHostEntity<0>(e);
      if(subGrid2.template contains<0>(eHost))
      {
        gr->insertRaw(eHost);
      }
    }
  }
  gr->createEnd();
  gr->setMaxLevelDifference(1);
  return gr;
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class Grid>
void
Periter<ScatteringKernelApproximation, RHSApproximation>::
create_new_grids(
      std::vector<std::unique_ptr<Grid>>& grids,
      size_t numNewGrids)
{
  std::vector<std::unique_ptr<Grid>> newGrids;
  newGrids.reserve(numNewGrids);
  newGrids.push_back(fullSubGrid<Grid>(grids[0]->getHostGrid()));
  const auto& subGrid = *newGrids.back();
  for(size_t i = 1; i < numNewGrids; i++) {
    newGrids.push_back(copySubGrid(subGrid));
  }

  std::swap(grids, newGrids);
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class GridIdSet, class Grid>
void
Periter<ScatteringKernelApproximation, RHSApproximation>::
save_grids_to_gridIdSets(
      std::vector<GridIdSet>& gridIdSets,
      const std::vector<std::unique_ptr<Grid>>& grids)
{
  gridIdSets.clear();
  gridIdSets.reserve(grids.size());
  for(const auto& gridPtr : grids) {
    const auto& grid = *gridPtr;
    gridIdSets.push_back(saveSubGridToIdSet(grid));
  }
}

} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_PERITER_HH
