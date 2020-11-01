// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_ASTI_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_ASTI_HH

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <iomanip>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/io.hh>
#include <dune/istl/umfpack.hh>

#include <dune/functions/functionspacebases/hangingnodebernsteinp2basis.hh>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/dpg/bilinearformfactory.hh>
#include <dune/dpg/innerproductfactory.hh>
#include <dune/dpg/linearformfactory.hh>
#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/errortools.hh>
#include <dune/dpg/functionplotter.hh>
#include <dune/dpg/functions/interpolate.hh>
#if ASTI_NORMALIZED_SPACES
#  include <dune/dpg/functions/normalizedspaces.hh>
#endif
#include <dune/dpg/functions/refinementinterpolation.hh>
#include <dune/dpg/functions/subgridinterpolation.hh>
#include <dune/dpg/linearfunctionalterm.hh>
#include <dune/dpg/radiative_transfer/approximate_scattering.hh>
#include <dune/dpg/radiative_transfer/boundary_extension.hh>
#include <dune/dpg/radiative_transfer/passkey.hh>
#include <dune/dpg/radiative_transfer/asti_common.hh>
#include <dune/dpg/radiative_transfer/subgridprojection.hh>
#include <dune/dpg/rhs_assembler.hh>
#include <dune/dpg/dpg_system_assembler.hh>
#include <dune/dpg/type_traits.hh>

#include "subgrids.hh"

#include <boost/math/constants/constants.hpp>

namespace Dune {

template<class GridView>
class SubGridSpaces;

class ASTILogger;
class TransportLogger;
class TransportPlotter;

/**
 * This class describes the ASTI algorithm for radiative transfer problems
 *
 * \tparam ScatteringKernelApproximation
 *         specifies the method used to approximate the scattering kernel
 * \tparam RHSApproximation  if right hand side and lifting of boundary
 *                           values are finite element functions, set this
 *                           to FeRHS, otherwise set this to ApproximateRHS
 */
template<class ScatteringKernelApproximation, class RHSApproximation>
class ASTI {
  public:

  using Direction = typename ScatteringKernelApproximation::Direction;

  /**
   * Solve a radiative transfer problem using the ASTI algorithm
   *
   * \param hostGrid
   * \param f  right hand side function
   * \param g  lifting of the boundary values
   * \param is_inflow_boundary_homogeneous
   *            checks if g is 0 on the inflow boundary
   * \param sigma   absorption coefficient
   * \param kernel  the scattering kernel, e.g. a Henyey–Greenstein kernel
   * \param approximationParameters
   * \param targetAccuracy  ASTI solves up to this accuracy
   * \param maxNumberOfIterations  ... or up to the given number of iterations
   *                               (whatever comes first)
   * \param plotFlags  specifies when to create .vtu files for plotting
   *                   the solution, scattering or rhs
   * \note By default, the scattering is computed from the inner degrees
   *       of freedom of the SolutionSpace. To use the skeletal degrees of
   *       freedom instead, define the preprocessor variable
   *       ASTI_SKELETAL_SCATTERING before including asti.hh.
   */
  template<class HostGrid, class F, class G, class HB,
           class Sigma, class Kernel>
  void solve(HostGrid& hostGrid,
             const F& f,
             const G& g,
             const HB& is_inflow_boundary_homogeneous,
             const Sigma sigma,
             const Kernel& kernel,
             ASTIApproximationParameters& approximationParameters,
             double targetAccuracy,
             unsigned int maxNumberOfIterations,
             unsigned int maxNumberOfInnerIterations,
             const std::string& outputfolder,
             ASTIPlotFlags plotFlags = ASTIPlotFlags::doNotPlot);

  private:
  using VectorType = BlockVector<FieldVector<double,1>>;

  template<class HostGridView, class SubGridView, class Grid,
           class GridIdSet, class Plotter>
  static auto computeScatteringData(
        HostGridView hostGridView,
        ScatteringKernelApproximation& kernelApproximation,
        std::vector<VectorType>& x,
        SubGridSpaces<SubGridView>& subGridSpaces,
        std::vector<Direction>& sVector,
        std::vector<std::unique_ptr<Grid>>& grids,
        std::vector<GridIdSet>& gridIdSets,
        double accuKernel,
        const Plotter& plotter,
        unsigned int n);

  template<class HostGridView, class SubGridView, class Grid, class F>
  static auto computeRhsData(
        HostGridView hostGridView,
        SubGridSpaces<SubGridView>& subGridSpaces,
        std::vector<std::unique_ptr<Grid>>& grids,
        double accuracy,
        const std::vector<Direction>& sVector,
        const F& f);

  template<class HostGridView, class SubGridView, class G, class HB>
  static auto computeBvData(
        HostGridView hostGridView,
        const SubGridSpaces<SubGridView>& subGridSpaces,
        const std::vector<Direction>& sVector,
        const G& g,
        const HB& is_inflow_boundary_homogeneous);

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
   * \param sigma the absorption coefficient
   * \param rhsVector the data for the unmodified rhs
   * \param bvExtension a lifting of the boundary values
   *
   * \return the squared a posteriori error of the solution
   */
  template<class Spaces, class Grid, class Sigma, class RHSVector>
  static double compute_transport_solution(
      VectorType& x,
      Grid& grid,
      const Spaces& spaces,
      const FieldVector<double, 2>& s,
      const Sigma sigma,
      RHSVector& rhsVector,
      const std::optional<VectorType>& bvExtension);

  /**
   * Adaptively compute transport solution
   *
   * This function calls compute_transport_solution iteratively
   * until the given accuracy or the maximal number of iterations
   * is reached.
   */
  template<class Spaces, class ScatteringHostGridBasis, class BvHostGridBasis,
           class Grid, class GridIdSet, class Sigma,
           class RHSData, class ScatteringData, class BVData>
  static double compute_adaptive_transport_solution(
      VectorType& x,
      Spaces& spaces,
      ScatteringHostGridBasis&& scatteringHostGridGlobalBasis,
      BvHostGridBasis&& bvHostGridGlobalBasis,
      Grid& grid,
      GridIdSet& gridIdSet,
      const FieldVector<double, 2>& s,
      const Sigma sigma,
      RHSData& rhsData,
      ScatteringData& scatteringData,
      std::optional<BVData>& bvData,
      TransportLogger transportLogger,
      TransportPlotter transportPlotter,
      double transportAccuracy,
      unsigned int maxNumberOfInnerIterations);

  template<class SubGridView, class HostGridBasis>
  static std::vector<VectorType>
  interpolate_solutions_to_hostgrid(
        const std::vector<VectorType>& solution,
        SubGridSpaces<SubGridView>& subGridSpaces,
        const HostGridBasis& hostGridBasis);

  /**
   * Apply the scattering integral to a solution x
   *
   * This corresponds to [K, u_n, κ_1 * η] in the ASTI algorithm
   * (see Dahmen, Gruber, Mula).
   *
   * \param kernelApproximation an approximation to the scattering kernel
   * \param x  solution to which we want to apply the scattering kernel
   * \param subGridSpaces
   * \param accuracy
   */
  template<class SubGridView, class HostGridBasis, class Grid, class GridIdSet>
  static std::vector<VectorType> apply_scattering(
      ScatteringKernelApproximation& kernelApproximation,
      const std::vector<VectorType>& x,
      SubGridSpaces<SubGridView>& subGridSpaces,
      const HostGridBasis& hostGridBasis,
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
  using HostGridView = typename GridView::Grid::HostGridType::LeafGridView;

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

#if ASTI_NORMALIZED_SPACES
  SubGridSpaces(const GridsVector& grids,
                const std::vector<FieldVector<double, 2>>& directions)
  {
    spaces_.reserve(grids.size());
    for(size_t i = 0, iend = grids.size(); i < iend; i++) {
      spaces_.emplace_back(grids[i]->leafGridView(), directions[i]);
    }
  }
#else
  SubGridSpaces(const GridsVector& grids)
  {
    spaces_.reserve(grids.size());
    for(const auto& grid : grids)
    {
      spaces_.emplace_back(grid->leafGridView());
    }
  }
#endif

  void update(size_t i, const GridView& gridView) {
    spaces_[i].update(gridView);
  }

  /**
   * insert new spaces in apply_scattering after adding new grids
   */
#if ASTI_NORMALIZED_SPACES
  void create_new_spaces(const GridsVector& grids,
                         const std::vector<FieldVector<double, 2>>& directions)
  {
    if(spaces_.size() != grids.size()) {
      spaces_.clear();
      spaces_.reserve(grids.size());
      for(size_t i = 0, iend = grids.size(); i < iend; i++) {
        spaces_.emplace_back(grids[i]->leafGridView(), directions[i]);
      }
    }
  }
#else
  void create_new_spaces(const GridsVector& grids)
  {
    if(spaces_.size() != grids.size()) {
      spaces_.clear();
      spaces_.reserve(grids.size());
      for(const auto& grid : grids) {
        spaces_.emplace_back(grid->leafGridView());
      }
    }
  }
#endif

  const FEBasisInterior& interiorSolutionSpace(size_t i) const {
    return spaces_[i].interiorSolutionSpace();
  }

  const FEBasisTrace& traceSolutionSpace(size_t i) const {
    return spaces_[i].traceSolutionSpace();
  }

  const FEBasisTest& testSpace(size_t i) const {
    return spaces_[i].testSpace();
  }

#ifdef ASTI_SKELETAL_SCATTERING
  static auto scatteringHostGridBasis(HostGridView hostGridView) {
    // Discontinuous version of the trace space
    using FEBasisHost = Functions::BernsteinDGBasis<HostGridView, 2>;
    return FEBasisHost(hostGridView);
  }

  const FEBasisTrace& scatteringSubGridBasis(size_t i) const {
    return traceSolutionSpace(i);
  }

  size_t scatteringSubGridBasisOffset(size_t i) const {
    return interiorSolutionSpace(i).size();
  }
#else
  static auto scatteringHostGridBasis(HostGridView hostGridView) {
#if ASTI_NORMALIZED_SPACES
    using FEBasisInteriorHost
        = changeGridView_t<typename Spaces::FEBasisInterior, HostGridView>;
    auto interiorSpace = make_space_tuple<FEBasisInteriorHost>(hostGridView);
    auto l2InnerProduct
      = innerProductWithSpace(interiorSpace)
        .template addIntegralTerm<0,0,IntegrationType::valueValue,
                                      DomainOfIntegration::interior>(1.)
        .create();
    return make_normalized_space(l2InnerProduct);
#else
    using FEBasisHost = changeGridView_t<FEBasisInterior, HostGridView>;
    return FEBasisHost(hostGridView);
#endif
  }

  const FEBasisInterior& scatteringSubGridBasis(size_t i) const {
    return interiorSolutionSpace(i);
  }

  size_t scatteringSubGridBasisOffset(size_t i) const {
    return 0;
  }
#endif

  static auto bvHostGridBasis(HostGridView hostGridView) {
    using FEBasisHostTrace = changeGridView_t<FEBasisTrace, HostGridView>;
    return FEBasisHostTrace(hostGridView);
  }

  Spaces& operator[](size_t i) {
    return spaces_[i];
  }

  private:

  std::vector<Spaces> spaces_;
};

class TransportLogger {
  public:
  explicit TransportLogger(std::ofstream& ofs,
                           unsigned int outerIteration,
                           unsigned int direction,
                           PassKey<ASTILogger>)
    : ofs(ofs)
    , n(outerIteration)
    , i(direction) {};

  void printCurrentIteration(const unsigned int nRefinement) const {
    std::cout << "Direction " << i
              << ", inner iteration " << nRefinement << '\n';
  }

  void logSolverStats(
      const unsigned int nRefinement,
      const double aposteriori_s,
      const int maxGridLevel,
      const size_t numDofs) {
    ofs << "Iteration " << n << '.' << nRefinement
        << " for direction " << i << ":\n"
        << "  - A posteriori estimation of || (u,trace u) - (u_fem,theta) || = "
        << aposteriori_s
        << "\n  - Grid level: "     << maxGridLevel
        << "\n  - Number of DOFs: " << numDofs
        << std::endl;

    std::cout << "\nIteration " << n << '.' << nRefinement
        << " for direction " << i << ":\n"
        << "  - A posteriori estimation of || (u,trace u) - (u_fem,theta) || = "
        << aposteriori_s
        << "\n  - Grid level: " << maxGridLevel
        << "\n  - Number of DOFs: " << numDofs
        << std::endl;
  }

  void logFinalAccuracy(const double aposteriori_s,
                        const double transportAccuracy)
  {
    ofs << "\na posteriori error for current direction: "
        << aposteriori_s;
    if (aposteriori_s <= transportAccuracy)
    {
      ofs << " (enough";
    } else {
      ofs << " (not enough, required "
          << transportAccuracy;
    }
    ofs << ")\n\n" << std::flush;
  }

  void logAccuracyAfterRefinement(const double aposteriori_s,
                                  const double transportAccuracy)
  {
    ofs << "\na posteriori error for current direction: "
        << aposteriori_s
        << " (required "
        << transportAccuracy
        << ")\n\n" << std::flush;
  }

  private:
  std::ofstream& ofs;
  const unsigned int n;
  const unsigned int i;
};

class ASTILogger {
  public:
  template<class ScatteringKernelApproximation, class RHSApproximation>
  explicit ASTILogger(std::string filename,
                         PassKey<ASTI<ScatteringKernelApproximation,
                                         RHSApproximation>>) : ofs(filename) {}

  TransportLogger transportLogger(unsigned int outerIteration,
                                  unsigned int direction)
  {
    return TransportLogger(ofs, outerIteration, direction, {});
  }

  template<class Kernel, class KernelApproximation>
  void logASTIOverview(
      const double targetAccuracy,
      const Kernel& kernel,
      const KernelApproximation& kernelApproximation,
      const ASTIApproximationParameters& approximationParameters)
  {
    ofs << "ASTI algorithm\n"
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
        << "ASTI parameters:" << '\n'
        << approximationParameters;
  }

  void logOuterIterationHeader(const unsigned int n)
  {
    ofs << "\nIteration n=" << n
        << "\n================\n";
    std::cout << "\nIteration " << n << "\n\n";
  }

  template<class KernelApproximation>
  void logKernelApproximationInfo(
      const ASTIApproximationParameters& approximationParameters,
      const KernelApproximation& kernelApproximation,
      const double accuKernel,
      const std::vector<FieldVector<double, 2>>& sVector,
      const std::chrono::steady_clock::time_point startScatteringApproximation,
      const std::chrono::steady_clock::time_point endScatteringApproximation)
  {
    ofs << "eta_n = (1+n)^{-alpha} rho^n: "
        << approximationParameters.eta() << '\n'
        << "\n--------------------\n"
        << "Info angular approx:\n"
        << "--------------------\n"
        << "Current wavelet level: "
        << kernelApproximation.getLevel() << '\n'
        << "Number of directions: "
        << kernelApproximation.getNumS()  << '\n'
        << "Directions are:\n";
    for(const auto& s : sVector) {
      ofs << s << '\n';
    }
    ofs << "\n---------------------\n"
        << "Kernel approximation:\n"
        << "---------------------\n"
        << "Accuracy required: "
          << approximationParameters.scatteringAccuracy() << " (kappa1 * eta)\n"
        << "Accuracy introduced in code: "
        << accuKernel << " (kappa1 * eta / (kappaNorm * uNorm))\n"
        << kernelApproximation.info() << '\n';
    if(kernelApproximation.typeApprox() == "Kernel approximation with: SVD")
    {
      ofs << "Singular values of kernel matrix:\n"
          << kernelApproximation.getSingularValues() << '\n';
    }
    ofs << "Computing time: "
          << std::chrono::duration_cast<std::chrono::microseconds>
          (endScatteringApproximation - startScatteringApproximation).count()
          << "us\n" << std::flush;
  }

  void logInnerIterationsHeader()
  {
    ofs << "\n-----------------------------------\n"
             "Inner iterations (transport solves)\n"
             "-----------------------------------\n";
  }

  void logInnerIterationIndex(const unsigned int i)
  {
    ofs << "\nAngle " << i
        << "\n--------\n";
  }

  template<class VectorType>
  void logInnerIterationStats(
      const std::vector<VectorType>& x,
      const double aposterioriTransportGlobal,
      const ASTIApproximationParameters& approximationParameters,
      const double deviationOfInexactIterate,
      const double accuracy,
      const unsigned int n)
  {
    const size_t accumulatedDoFs = std::accumulate(x.cbegin(), x.cend(),
        static_cast<size_t>(0),
        [](size_t acc, auto vec) { return acc + vec.size(); });

    // Error bound for || u - \bar u_n || based on a posteriori errors
    const double aPosterioriError
        = approximationParameters
          .combinedAPosterioriError(deviationOfInexactIterate);

    ofs << "--------------------\n"
           "End inner iterations\n"
           "--------------------\n"
           "Error transport solves (a posteriori estimation): "
          << aposterioriTransportGlobal                  << '\n'
        << "Accuracy kernel: "
          << approximationParameters.scatteringAccuracy() << '\n'
        << "Error bound ||u_n - bar u_n|| (a posteriori): "
          << deviationOfInexactIterate << '\n'
        << "A priori bound global accuracy ||u - bar u_n||: "
          << accuracy
          << " (rho * err0 + 2) * rho^n\n"
        << "A posteriori bound global accuracy ||u - bar u_n||: "
          << aPosterioriError << '\n'
        << "Total number of DoFs: "
          << accumulatedDoFs
        << "\n\n" << std::flush;

    std::cout << "Error at end of Iteration " << n << ": "
              << aposterioriTransportGlobal << ", using "
              << accumulatedDoFs << " DoFs\n";
  }

  private:
  std::ofstream ofs;
};

class TransportPlotter {
  public:
  TransportPlotter(
      ASTIPlotFlags plotFlags,
      std::string outputfolder,
      unsigned int outerIteration,
      unsigned int direction)
    : plotFlags(plotFlags)
    , outputfolder(outputfolder)
    , n(outerIteration)
    , i(direction)
  {};

  template<class RhsGlobalBasis, class VectorType>
  void plotRhsAndScattering(
      const RhsGlobalBasis& feBasisTest,
      const VectorType& rhsValues,
      const VectorType& scatteringValues,
      const unsigned int nRefinement) const
  {
    if(flagIsSet(plotFlags, ASTIPlotFlags::plotRhs)) {
      std::string name = outputfolder
                        + "/rhs_rad_trans_n"
                        + std::to_string(n)
                        + "_s"
                        + std::to_string(i)
                        + '_' + std::to_string(nRefinement);
      FunctionPlotter rhsPlotter(name);
      rhsPlotter.plot("rhs", rhsValues, feBasisTest, 1);

      name = outputfolder
                        + "/scattering_rad_trans_n"
                        + std::to_string(n)
                        + "_s"
                        + std::to_string(i)
                        + '_' + std::to_string(nRefinement);
      FunctionPlotter scatteringPlotter(name);
      scatteringPlotter.plot("rhs", scatteringValues, feBasisTest, 1);
    }
  }

  private:
  const ASTIPlotFlags plotFlags;
  const std::string outputfolder;
  const unsigned int n;
  const unsigned int i;
};

class ASTIPlotter {
  public:
  ASTIPlotter(ASTIPlotFlags plotFlags, std::string outputfolder)
    : plotFlags(plotFlags)
    , outputfolder(outputfolder)
  {
    if(flagIsSet(plotFlags, ASTIPlotFlags::plotLastIteration)) {
      std::cerr
          << "Plotting of only the last iteration is not implemented yet!\n";
      std::exit(1);
    }
  };

  TransportPlotter transportPlotter(unsigned int outerIteration,
                                    unsigned int direction)
  {
    return TransportPlotter(plotFlags, outputfolder,
                            outerIteration, direction);
  }

  template<class Spaces, class VectorType>
  void plotSolutions(
      const Spaces& spaces,
      const VectorType& x,
      const unsigned int n,
      const unsigned int i) const
  {
    if(flagIsSet(plotFlags, ASTIPlotFlags::plotOuterIterations)) {
      //////////////////////////////////////////////////////////////////////
      //  Write result to VTK file
      //  We need to subsample, because VTK cannot natively display
      //  real second-order functions
      //////////////////////////////////////////////////////////////////////

      const auto& feBasisInterior = spaces.interiorSolutionSpace();
      const auto& feBasisTrace    = spaces.traceSolutionSpace();

      std::cout << "Plot solution for direction " << i << '\n';

      std::string name = outputfolder
                        + "/u_rad_trans_n"
                        + std::to_string(n)
                        + "_s"
                        + std::to_string(i);
      FunctionPlotter uPlotter(name);
      uPlotter.plot("u", x, feBasisInterior, 0, 0);
      name = outputfolder
                        + "/theta_rad_trans_n"
                        + std::to_string(n)
                        + "_s"
                        + std::to_string(i);
      FunctionPlotter thetaPlotter(name);
      thetaPlotter.plot("theta", x, feBasisTrace, 2,
                        feBasisInterior.size());
    }
  }

  template<class ScatteringBasis, class VectorType>
  void plotScatteringOnHostGrid(
      const ScatteringBasis& hostGridGlobalBasis,
      const std::vector<VectorType>& scattering,
      const unsigned int n) const
  {
    if(flagIsSet(plotFlags, ASTIPlotFlags::plotScattering)) {
      std::cout << "Plot scattering:\n";

      for(size_t i = 0, numS = scattering.size(); i < numS; ++i)
      {
        std::cout << "Direction " << i << '\n';

        std::string name = outputfolder
                        + "/scattering_n"
                        + std::to_string(n)
                        + "_s"
                        + std::to_string(i);
        FunctionPlotter scatteringPlotter(name);
        scatteringPlotter.plot("scattering", scattering[i],
                               hostGridGlobalBasis, 0);
      }
    }
  }

  template<class ScatteringBasis, class VectorType, class KernelApproximation>
  void plotIntegratedSolution(
      const ScatteringBasis& hostGridGlobalBasis,
      const std::vector<VectorType>& solutionHost,
      const KernelApproximation& kernelApproximation,
      const unsigned int n) const
  {
    if(plotIntegratedSolutionEnabled()) {
      std::cout << "\nPlotting integrated solution of outer iteration "
                << n << ".\n";

      VectorType integratedSolution(hostGridGlobalBasis.size());
      integratedSolution = 0.;

      const auto quadWeightsSubinterval =
          kernelApproximation.quadWeightsOfSubintervalOnCurrentLevel();
      const auto numSperInterval = kernelApproximation.numSperInterval;
      const size_t numS = solutionHost.size();
      for(unsigned int i = 0; i < numS; ++i)
      {
        const auto quadWeight = quadWeightsSubinterval[i % numSperInterval];
        integratedSolution.axpy(quadWeight, solutionHost[i]);
      }

      std::string name = outputfolder
                      + "/integrated_solution_n"
                      + std::to_string(n);
      FunctionPlotter scatteringPlotter(name);
      scatteringPlotter.plot("integrated solution", integratedSolution,
                             hostGridGlobalBasis, 0);
    }
  }

  bool plotIntegratedSolutionEnabled() const {
    return flagIsSet(plotFlags, ASTIPlotFlags::plotIntegratedSolution);
  }

  private:
  const ASTIPlotFlags plotFlags;
  const std::string outputfolder;
};

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class HostGrid, class F, class G, class HB,
         class Sigma, class Kernel>
void ASTI<ScatteringKernelApproximation, RHSApproximation>::solve(
           HostGrid& hostGrid,
           const F& f,
           const G& g,
           const HB& is_inflow_boundary_homogeneous,
           const Sigma sigma,
           const Kernel& kernel,
           ASTIApproximationParameters& approximationParameters,
           double targetAccuracy,
           unsigned int maxNumberOfIterations,
           unsigned int maxNumberOfInnerIterations,
           const std::string& outputfolder,
           ASTIPlotFlags plotFlags) {
  static_assert(std::is_same<RHSApproximation, FeRHS>::value
      || std::is_same<RHSApproximation, ApproximateRHS>::value,
      "Unknown type provided for RHSApproximation!\n"
      "Should be either FeRHS or ApproximateRHS.");

  const unsigned int dim = HostGrid::dimension;
  using Grid = SubGrid<dim, HostGrid, false>;
  using LeafGridView = typename Grid::LeafGridView;
  using Direction = FieldVector<double, dim>;
  using GridIdSet = std::set<typename HostGrid::GlobalIdSet::IdType>;

  /////////////////////////////////////////////
  // To print information in dune-dpg/results/
  /////////////////////////////////////////////

  ASTILogger logger(outputfolder+"/output",
      PassKey<ASTI<ScatteringKernelApproximation, RHSApproximation>>{});
  ASTIPlotter plotter(plotFlags, outputfolder);

  ////////////////////////////////////////////
  // Handle directions of discrete ordinates
  ////////////////////////////////////////////

  ScatteringKernelApproximation kernelApproximation(kernel,
      2, approximationParameters.scatteringTruthLevel());

  logger.logASTIOverview(targetAccuracy, kernel,
      kernelApproximation, approximationParameters);

  // As the solution u we use for the initial scattering is 0, and the
  // formula for the accuracy contains a 1/\|u\|, we set the initial
  // accuracy to a large enough value.
  std::vector<Direction>
    sVector(kernelApproximation.setInitialAccuracy(1e5));
  size_t numS = sVector.size();

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
  //   Solution vectors
  //////////////////////////////////
  std::vector<VectorType> x(numS);
  using Spaces = SubGridSpaces<LeafGridView>;
#if ASTI_NORMALIZED_SPACES
  Spaces spaces(grids, sVector);
#else
  Spaces spaces(grids);
#endif

  for(unsigned int i = 0; i < grids.size(); ++i)
  {
    x[i].resize(spaces.interiorSolutionSpace(i).size()
                 + spaces.traceSolutionSpace(i).size());
    x[i] = 0;
  }

  double uNorm = 0.;

  /////////////////////////////////////////////////////////
  //  Fixed-point iterations
  /////////////////////////////////////////////////////////

  double accuracy = approximationParameters.aPrioriAccuracy();
  double deviationOfInexactIterate = 0.;

  for(unsigned int n = 0; accuracy > targetAccuracy
                          && n < maxNumberOfIterations; ++n)
  {
    logger.logOuterIterationHeader(n);

    for(size_t i = 0; i < grids.size(); ++i) {
      grids[i] = restoreSubGridFromIdSet<Grid>(gridIdSets[i],
                                               hostGrid);
      spaces.update(i, grids[i]->leafGridView());
    }

    auto startScatteringApproximation = std::chrono::steady_clock::now();

    const double kappaNorm = 1.;

    const double accuKernel = approximationParameters.scatteringAccuracy()
                                           // To prevent division by zero.
                            / (kappaNorm * ((uNorm>0.)?uNorm:1e-5));

    auto scatteringData = computeScatteringData
                            (hostGrid.leafGridView(),
                             kernelApproximation, x, spaces,
                             sVector, grids, gridIdSets, accuKernel,
                             plotter, n);
    auto rhsData = computeRhsData
                     (hostGrid.leafGridView(), spaces,
                      grids, approximationParameters.rhsAccuracy(),
                      sVector, f);
    auto bvData = computeBvData
                    (hostGrid.leafGridView(), spaces,
                     sVector, g, is_inflow_boundary_homogeneous);
    numS = sVector.size();

    auto endScatteringApproximation = std::chrono::steady_clock::now();

    logger.logKernelApproximationInfo(approximationParameters,
        kernelApproximation, accuKernel, sVector,
        startScatteringApproximation, endScatteringApproximation);


    ////////////////////////////////////////////////////
    // Inner loop
    ////////////////////////////////////////////////////
    logger.logInnerIterationsHeader();
    double aposterioriTransportGlobal = 0.;

    // Loop over the spatial grids.
    // For each angular direction there is one spatial grid.
    for(unsigned int i = 0; i < grids.size(); ++i)
    {
      logger.logInnerIterationIndex(i);

      grids[i] = restoreSubGridFromIdSet<Grid>(gridIdSets[i],
                                               hostGrid);
      spaces.update(i, grids[i]->leafGridView());

      const double aposteriori_s = compute_adaptive_transport_solution(
          x[i], spaces[i],
          Spaces::scatteringHostGridBasis(hostGrid.leafGridView()),
          Spaces::bvHostGridBasis(hostGrid.leafGridView()),
          *grids[i], gridIdSets[i], sVector[i], sigma,
          rhsData[i], scatteringData[i], bvData[i],
          logger.transportLogger(n, i),
          plotter.transportPlotter(n, i),
          approximationParameters.transportAccuracy(),
          maxNumberOfInnerIterations);
      aposterioriTransportGlobal = std::max(aposterioriTransportGlobal,
                                            aposteriori_s);
      plotter.plotSolutions(spaces[i], x[i], n, i);
    }

    // Error bound for || u_{n+1} - bar u_{n+1} ||
    deviationOfInexactIterate
      = approximationParameters.errorBetweenExactAndInexactIterate(
                                          deviationOfInexactIterate,
                                          aposterioriTransportGlobal);

    // Error bound for || u - \bar u_n || based on a priori errors
    accuracy = approximationParameters.combinedAccuracy();

    logger.logInnerIterationStats(x, aposterioriTransportGlobal,
        approximationParameters, deviationOfInexactIterate, accuracy, n);

    for(unsigned int i = 0; i < grids.size(); ++i)
    {
      grids[i] = restoreSubGridFromIdSet<Grid>(gridIdSets[i],
                                               hostGrid);
      spaces.update(i, grids[i]->leafGridView());
    }

    // compute L2 norm of u
    {
      uNorm = 0.;
      const std::vector<double> quadWeight
        = kernelApproximation.quadWeightsOfSubintervalOnCurrentLevel();
      for(size_t i=0; i<grids.size(); ++i) {
        const double uiNorm =
          ErrorTools::l2norm(spaces.interiorSolutionSpace(i), x[i]);
        uNorm += uiNorm * uiNorm
                  * quadWeight[i % kernelApproximation.numSperInterval];
      }
      uNorm = std::sqrt(uNorm);
    }

    approximationParameters.decreaseEta(uNorm);

    if(plotter.plotIntegratedSolutionEnabled()) {
      auto hostGridGlobalBasis
          = Spaces::scatteringHostGridBasis(hostGrid.leafGridView());

      const auto xHost =
          interpolate_solutions_to_hostgrid(x, spaces, hostGridGlobalBasis);

      // plotting integrated solution for the last iteration.
      plotter.plotIntegratedSolution(hostGridGlobalBasis,
                                     xHost,
                                     kernelApproximation,
                                     n);
    }
  }
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class HostGridView, class SubGridView, class Grid, class GridIdSet,
         class Plotter>
auto
ASTI<ScatteringKernelApproximation, RHSApproximation>::
computeScatteringData(
      HostGridView hostGridView,
      ScatteringKernelApproximation& kernelApproximation,
      std::vector<VectorType>& x,
      SubGridSpaces<SubGridView>& subGridSpaces,
      std::vector<Direction>& sVector,
      std::vector<std::unique_ptr<Grid>>& grids,
      std::vector<GridIdSet>& gridIdSets,
      double accuKernel,
      const Plotter& plotter,
      unsigned int n)
{
  using Spaces = SubGridSpaces<SubGridView>;

  auto hostGridGlobalBasis = Spaces::scatteringHostGridBasis(hostGridView);
  using FEBasisHost = decltype(hostGridGlobalBasis);

  using ScatteringData = std::decay_t<decltype(attachDataToSubGrid(
          std::declval<typename Spaces::FEBasisTest>(),
          std::declval<FEBasisHost>(),
          std::declval<VectorType>()))>;
  std::vector<ScatteringData> scatteringData;

  std::vector<VectorType> scatteringFunctional =
      apply_scattering (
        kernelApproximation, x, subGridSpaces,
        hostGridGlobalBasis,
        sVector, grids, gridIdSets,
        accuKernel);

  const size_t numS = sVector.size();
  x.resize(numS);
  scatteringData.reserve(numS);

  plotter.plotScatteringOnHostGrid(hostGridGlobalBasis,
                                   scatteringFunctional,
                                   n);

  for(size_t i = 0; i < numS; i++) {
    // TODO: subGridGlobalBasis should be normalized in L2!?
    const auto& subGridGlobalBasis = subGridSpaces.testSpace(i);
    scatteringData.emplace_back(
      attachDataToSubGrid(
        subGridGlobalBasis,
        hostGridGlobalBasis,
        scatteringFunctional[i]));
  }
  return scatteringData;
}

#ifndef DOXYGEN
namespace detail {
  constexpr size_t dim = 2;
  using VectorType = BlockVector<FieldVector<double,1>>;
  using Direction = FieldVector<double, dim>;

  template<class FEHostBasis, class Grids, class F>
  inline std::vector<VectorType> approximate_rhs (
      Grids&,
      double,
      const FEHostBasis& hostGridBasis,
      const std::vector<Direction>& sVector,
      const F& f,
      FeRHS) {
    static_assert(!is_RefinedFiniteElement<FEHostBasis>::value,
        "Functions::interpolate won't work for refined finite elements");
    const size_t numS = sVector.size();
    std::vector<VectorType> rhsFunctional(numS);
    auto rhsFunction = f(hostGridBasis.gridView());
    for(auto& rhsFunctionalEntry : rhsFunctional)
    {
      rhsFunctionalEntry.resize(hostGridBasis.size());
      Functions::interpolate(hostGridBasis, rhsFunctionalEntry, rhsFunction);
    }
    return rhsFunctional;
  }


  // Refine grid until accuracy kappa2*eta is reached for
  // approximation of rhs.
  template<class FEBases, class Grids, class F>
  inline std::vector<VectorType> approximate_rhs (
      Grids& grids,
      double accuracy,
      const std::vector<std::shared_ptr<FEBases>>& solutionSpaces,
      const std::vector<Direction>& sVector,
      const F& f,
      ApproximateRHS) {
    DUNE_THROW(Dune::NotImplemented,
        "Implementation of approximate_rhs for non-FE right hand side "
        "functions currently broken.");
  }
} // end namespace detail
#endif

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class HostGridView, class SubGridView, class Grid, class F>
auto
ASTI<ScatteringKernelApproximation, RHSApproximation>::
computeRhsData(
      HostGridView hostGridView,
      SubGridSpaces<SubGridView>& subGridSpaces,
      std::vector<std::unique_ptr<Grid>>& grids,
      double accuracy,
      const std::vector<Direction>& sVector,
      const F& f)
{
  using Spaces = SubGridSpaces<SubGridView>;

  // TODO: is this the right space for the rhs function?
  auto hostGridGlobalBasis = Spaces::scatteringHostGridBasis(hostGridView);
  using FEBasisHostInterior = decltype(hostGridGlobalBasis);
  using RHSData = std::decay_t<decltype(attachDataToSubGrid(
            std::declval<typename Spaces::FEBasisTest>(),
            std::declval<FEBasisHostInterior>(),
            std::declval<VectorType>()))>;
  std::vector<RHSData> rhsData;

  std::vector<VectorType> rhsFunctional =
      detail::approximate_rhs (
          grids,
          accuracy,
          hostGridGlobalBasis,
          sVector,
          f, RHSApproximation{});

  const size_t numS = sVector.size();
  rhsData.reserve(numS);

  // TODO: restore grids from gridIdSets and update spaces
  //       shouldn't be necessary, as long as RHSApproximation == FeRHS
  for(size_t i = 0; i < numS; i++) {
    // TODO: subGridGlobalBasis should be normalized in L2!
    const auto& subGridGlobalBasis = subGridSpaces.testSpace(i);
    rhsData.emplace_back(
      attachDataToSubGrid(
        subGridGlobalBasis,
        hostGridGlobalBasis,
        rhsFunctional[i]));
  }
  return rhsData;
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class HostGridView, class SubGridView, class G, class HB>
auto
ASTI<ScatteringKernelApproximation, RHSApproximation>::
computeBvData(
      HostGridView hostGridView,
      const SubGridSpaces<SubGridView>& subGridSpaces,
      const std::vector<Direction>& sVector,
      const G& g,
      const HB& is_inflow_boundary_homogeneous)
{
  using Spaces = SubGridSpaces<SubGridView>;
  auto hostGridGlobalBasis = Spaces::bvHostGridBasis(hostGridView);
  using FEBasisHostTrace = decltype(hostGridGlobalBasis);
  using BVData = std::decay_t<decltype(attachDataToSubGrid(
            std::declval<typename Spaces::FEBasisTest>(),
            std::declval<FEBasisHostTrace>(),
            std::declval<VectorType>()))>;
  std::vector<std::optional<BVData>> bvData;

  const VectorType boundaryExtension
      = harmonic_extension_of_boundary_values(g, hostGridGlobalBasis);
  const size_t numS = sVector.size();
  bvData.reserve(numS);
  for(size_t i = 0; i < numS; i++) {
    if(is_inflow_boundary_homogeneous(sVector[i])) {
      bvData.emplace_back();
    } else {
      // TODO: subGridGlobalBasis should be normalized in L2!
      const auto& subGridGlobalBasis = subGridSpaces.testSpace(i);
      bvData.emplace_back(
        attachDataToSubGrid(
          subGridGlobalBasis,
          hostGridGlobalBasis,
          boundaryExtension));
    }
  }
  return bvData;
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class Spaces, class Grid, class Sigma, class RHSVector>
double
ASTI<ScatteringKernelApproximation, RHSApproximation>::
compute_transport_solution(
    VectorType& x,
    Grid& grid,
    const Spaces& spaces,
    const FieldVector<double, 2>& s,
    const Sigma sigma,
    RHSVector& rhsFunctional,
    const std::optional<VectorType>& bvExtension)
{
  auto bilinearForm =
    bilinearFormWithSpaces(spaces.testSpacePtr(), spaces.solutionSpacePtr())
    .template addIntegralTerm<0,0, IntegrationType::valueValue,
                                   DomainOfIntegration::interior>
                             (sigma(spaces.gridView()))
    .template addIntegralTerm<0,0, IntegrationType::gradValue,
                                   DomainOfIntegration::interior>(-1., s)
    .template addIntegralTerm<0,1, IntegrationType::normalVector,
                                   DomainOfIntegration::face>(1., s)
    .create();
  auto bilinearFormEnriched =
      replaceTestSpaces(bilinearForm, spaces.enrichedTestSpacePtr());
  auto innerProduct = spaces.testInnerProduct(s);
  auto innerProductEnriched =
      replaceTestSpaces(innerProduct, spaces.enrichedTestSpacePtr());

  auto systemAssembler =
      make_DPGSystemAssembler(bilinearForm,
                              innerProduct);

  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

  VectorType rhs;
  MatrixType stiffnessMatrix;

  /////////////////////////////////////////////////////////
  //  Assemble the systems
  /////////////////////////////////////////////////////////
  if(!bvExtension) {
    auto rhsFunction
      = linearFormWithSpace(systemAssembler.getTestSearchSpaces())
        .template addFunctionalTerm<0>(rhsFunctional, spaces.testSpace())
        .create();
    systemAssembler.assembleSystem(
        stiffnessMatrix, rhs,
        rhsFunction);
  } else {
    assert(bvExtension->size() == spaces.testSpace().size());
    auto rhsFunction
      = linearFormWithSpace(systemAssembler.getTestSearchSpaces())
        .template addFunctionalTerm<0>(rhsFunctional, spaces.testSpace())
        .template addSkeletalFunctionalTerm<0, IntegrationType::normalVector>
              (*bvExtension, spaces.testSpace(), -1., s)
        .create();
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
    systemAssembler.template applyHomogeneousDirichletBoundary<1>
        (stiffnessMatrix,
         rhs,
         dirichletNodesInflow);
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
  if(!bvExtension) {
    auto rhsAssemblerEnriched
      = make_RhsAssembler(spaces.enrichedTestSpacePtr());
    auto rhsFunction
      = linearFormWithSpace(rhsAssemblerEnriched.getTestSpaces())
        .template addFunctionalTerm<0>
              (rhsFunctional,
               std::get<0>(*systemAssembler.getTestSearchSpaces()))
        .create();
    rhsAssemblerEnriched.assembleRhs(rhs, rhsFunction);
  } else {
    auto rhsAssemblerEnriched
      = make_RhsAssembler(spaces.enrichedTestSpacePtr());
    auto rhsFunction
      = linearFormWithSpace(rhsAssemblerEnriched.getTestSpaces())
        .template addFunctionalTerm<0>
              (rhsFunctional,
               std::get<0>(*systemAssembler.getTestSearchSpaces()))
        .template addSkeletalFunctionalTerm<0, IntegrationType::normalVector>
              (*bvExtension,
               std::get<0>(*systemAssembler.getTestSearchSpaces()), -1., s)
        .create();
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
template<class Spaces, class ScatteringHostGridBasis, class BvHostGridBasis,
         class Grid, class GridIdSet, class Sigma,
         class RHSData, class ScatteringData, class BVData>
double
ASTI<ScatteringKernelApproximation, RHSApproximation>::
compute_adaptive_transport_solution(
    VectorType& x,
    Spaces& spaces,
    ScatteringHostGridBasis&& scatteringHostGridGlobalBasis,
    BvHostGridBasis&& bvHostGridGlobalBasis,
    Grid& grid,
    GridIdSet& gridIdSet,
    const FieldVector<double, 2>& s,
    const Sigma sigma,
    RHSData& rhsData,
    ScatteringData& scatteringData,
    std::optional<BVData>& bvData,
    TransportLogger transportLogger,
    TransportPlotter transportPlotter,
    const double transportAccuracy,
    const unsigned int maxNumberOfInnerIterations)
{
  double aposteriori_s = 0.;

  for(unsigned int nRefinement = 0; ; )
    // At the end of the loop, we will break if
    // aposteriori_s < kapp3*eta (pointwise in s)
    // or ++nRefinement >= maxNumberOfInnerIterations
    // thus the inner loop terminates eventually.
  {
    std::optional<VectorType> bvExtension;
    if(bvData) {
      // TODO: subGridGlobalBasis should be normalized in L2!
      const auto& subGridGlobalBasis = spaces.testSpace();
      const auto newGridData
          = bvData->restoreDataToRefinedSubGrid(subGridGlobalBasis,
                                                bvHostGridGlobalBasis);
      bvExtension = VectorType(newGridData.size());
      std::copy(newGridData.cbegin(), newGridData.cend(),
                bvExtension->begin());
    }

    {
      transportLogger.printCurrentIteration(nRefinement);

      VectorType rhsFunctional;
      {
        auto& feBasisTest = spaces.testSpace();
        const auto rhsValues
            = rhsData.restoreDataToRefinedSubGrid(feBasisTest,
                                    scatteringHostGridGlobalBasis);

        const auto scatteringValues
            = scatteringData.restoreDataToRefinedSubGrid(feBasisTest,
                                           scatteringHostGridGlobalBasis);

        transportPlotter.plotRhsAndScattering(feBasisTest,
                                rhsValues, scatteringValues, nRefinement);

        rhsFunctional.resize(rhsValues.size());
        std::transform(rhsValues.cbegin(), rhsValues.cend(),
                       scatteringValues.cbegin(),
                       rhsFunctional.begin(),
                       std::plus<>());
      }

      aposteriori_s
          = compute_transport_solution(x, grid,
              spaces, s, sigma, rhsFunctional, bvExtension);

      aposteriori_s = std::sqrt(aposteriori_s);

      // TODO: Add (interpolation of) g to theta part of x?

      transportLogger.logSolverStats(nRefinement, aposteriori_s,
          grid.maxLevel(), x.size());
    }


    if(++nRefinement >= maxNumberOfInnerIterations
        || aposteriori_s <= transportAccuracy) {
      gridIdSet = saveSubGridToIdSet(grid);

      transportLogger.logFinalAccuracy(aposteriori_s, transportAccuracy);

      break;
    } else {
      grid.preAdapt();
      grid.adapt();
      grid.postAdapt();
      spaces.update(grid.leafGridView());
      scatteringHostGridGlobalBasis.update(grid.getHostGrid().leafGridView());
      bvHostGridGlobalBasis.update(grid.getHostGrid().leafGridView());

      transportLogger.logAccuracyAfterRefinement(aposteriori_s,
                                                 transportAccuracy);
    }
  }
  return aposteriori_s;
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class SubGridView, class HostGridBasis>
std::vector<typename ASTI<ScatteringKernelApproximation,
                             RHSApproximation>::VectorType>
ASTI<ScatteringKernelApproximation, RHSApproximation>::
interpolate_solutions_to_hostgrid(
      const std::vector<VectorType>& solution,
      SubGridSpaces<SubGridView>& subGridSpaces,
      const HostGridBasis& hostGridBasis)
{
  std::vector<VectorType> solutionHost(solution.size());
  for(size_t i = 0, imax = solution.size(); i < imax; i++) {
    const auto& subGridBasis = subGridSpaces.scatteringSubGridBasis(i);
    const size_t offset = subGridSpaces.scatteringSubGridBasisOffset(i);
    interpolateFromSubGrid(
        subGridBasis, solution[i].begin() + offset,
        hostGridBasis, solutionHost[i]);
  }
  return solutionHost;
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class SubGridView, class HostGridBasis, class Grid, class GridIdSet>
std::vector<typename ASTI<ScatteringKernelApproximation,
                             RHSApproximation>::VectorType>
ASTI<ScatteringKernelApproximation, RHSApproximation>::apply_scattering(
      ScatteringKernelApproximation& kernelApproximation,
      const std::vector<VectorType>& x,
      SubGridSpaces<SubGridView>& subGridSpaces,
      const HostGridBasis& hostGridBasis,
      std::vector<Direction>& sVector,
      std::vector<std::unique_ptr<Grid>>& grids,
      std::vector<GridIdSet>& gridIdSets,
      double accuracy)
{
  sVector = kernelApproximation.setAccuracyAndInputSize(accuracy, x.size());

  const auto xHost =
      interpolate_solutions_to_hostgrid(x, subGridSpaces, hostGridBasis);

  const auto scatteringAssembler =
      make_ApproximateScatteringAssembler(hostGridBasis,
                                          kernelApproximation);
  const size_t numS = sVector.size();
  std::vector<VectorType> rhsFunctional(numS);
  scatteringAssembler.computeScattering(rhsFunctional, xHost);

  create_new_grids(grids, numS);
  save_grids_to_gridIdSets(gridIdSets, grids);
#if ASTI_NORMALIZED_SPACES
  subGridSpaces.create_new_spaces(grids, sVector);
#else
  subGridSpaces.create_new_spaces(grids);
#endif

  return rhsFunctional;
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class Grid>
void
ASTI<ScatteringKernelApproximation, RHSApproximation>::
create_new_grids(
      std::vector<std::unique_ptr<Grid>>& grids,
      size_t numNewGrids)
{
  const size_t numOldGrids = grids.size();
  if(numOldGrids == numNewGrids) {
    // no new grids need to be added
    return;
  }
  std::vector<std::unique_ptr<Grid>> newGrids;
  newGrids.reserve(numNewGrids);
  const size_t numCopies = numNewGrids / numOldGrids;
  const size_t lastOldGrid = numOldGrids-1;
  for(size_t i = 0; i < lastOldGrid; i++) {
    newGrids.push_back(intersectSubGrids(*grids[i], *grids[i+1]));
    const Grid& intersectionGrid = *newGrids.back();
    for(size_t copies = 1; copies < numCopies; ++copies) {
      newGrids.push_back(copySubGrid(intersectionGrid));
    }
  }
  // last one needs to be handled extra, as one of the neighboring grids is
  // the first one (since the direction parameter lives on the unit sphere).
  newGrids.push_back(intersectSubGrids(*grids[lastOldGrid], *grids[0]));
  const Grid& intersectionGrid = *newGrids.back();
  for(size_t copies = 1; copies < numCopies; ++copies) {
    newGrids.push_back(copySubGrid(intersectionGrid));
  }

  std::swap(grids, newGrids);
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class GridIdSet, class Grid>
void
ASTI<ScatteringKernelApproximation, RHSApproximation>::
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

#endif // DUNE_DPG_RADIATIVE_TRANSFER_ASTI_HH
