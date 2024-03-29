// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_ASTI_UNIFORM_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_ASTI_UNIFORM_HH

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <string>
#include <vector>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/io.hh>
#include <dune/istl/umfpack.hh>

#include <dune/functions/functionspacebases/bernsteinbasis.hh>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/dpg/bilinearformfactory.hh>
#include <dune/dpg/innerproductfactory.hh>
#include <dune/dpg/linearformfactory.hh>
#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/errortools.hh>
#include <dune/dpg/functionplotter.hh>
#include <dune/dpg/functions/gridviewfunctions.hh>
#include <dune/dpg/functions/interpolate.hh>
#include <dune/dpg/linearfunctionalterm.hh>
#include <dune/dpg/radiative_transfer/approximate_scattering.hh>
#include <dune/dpg/radiative_transfer/boundary_extension.hh>
#include <dune/dpg/radiative_transfer/passkey.hh>
#include <dune/dpg/radiative_transfer/asti_common.hh>
#include <dune/dpg/rhs_assembler.hh>
#include <dune/dpg/dpg_system_assembler.hh>
#include <dune/dpg/type_traits.hh>

#include <boost/math/constants/constants.hpp>

namespace Dune {

class ASTILogger;
class TransportLogger;

/**
 * This class describes the ASTI algorithm for radiative transfer problems
 *
 * \tparam ScatteringKernelApproximation
 *         specifies the method used to approximate the scattering kernel
 * \tparam RHSApproximation  if right hand side and lifting of boundary
 *                           values are finite element functions, set this
 *                           to FeRHS, otherwise set this to
 *                           ApproximateRHS
 */
template<class ScatteringKernelApproximation, class RHSApproximation>
class ASTI {
  public:

  using Direction = typename ScatteringKernelApproximation::Direction;

  /**
   * Solve a radiative transfer problem using the ASTI algorithm
   *
   * \param grid
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
   */
  template<class Grid, class F, class G, class HB, class Sigma, class Kernel>
  void solve(Grid& grid,
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
             ASTIPlotFlags plotFlags = ASTIPlotFlags::doNotPlot,
             ASTILogFlags logFlags = ASTILogFlags::defaultLog);

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
   * \param spaces
   * \param s the transport direction
   * \param sigma the absorption coefficient
   * \param rhsFunctional the data for the unmodified rhs
   * \param boundary_is_homogeneous
   * \param bvExtension a lifting of the boundary values
   *
   * \return the squared a posteriori error of the solution
   */
  template<class Spaces, class Sigma>
  static double compute_transport_solution(
      VectorType& x,
      const Spaces& spaces,
      const FieldVector<double, 2>& s,
      const Sigma sigma,
      const VectorType& rhsFunctional,
      bool boundary_is_homogeneous,
      const VectorType& bvExtension);

  /**
   * Adaptively compute transport solution
   *
   * This function calls compute_transport_solution iteratively
   * until the given accuracy or the maximal number of iterations
   * is reached.
   */
  template<class Spaces, class Grid, class Sigma, class KernelApproximation>
  static double compute_adaptive_transport_solution(
      std::vector<VectorType>& x,
      Spaces& spaces,
      Grid& grid,
      const std::vector<FieldVector<double, 2>>& sVector,
      const Sigma sigma,
      std::vector<VectorType>& rhsFunctional,
      VectorType& boundaryExtension,
      const std::vector<bool>& boundary_is_homogeneous,
      ASTILogger& logger,
      unsigned int n,
      double transportAccuracy,
      const KernelApproximation& kernelApproximation,
      unsigned int maxNumberOfInnerIterations);

  /**
   * Apply the scattering integral to a solution x
   *
   * This corresponds to [K, u_n, κ_1 * η] in the ASTI algorithm
   * (see Dahmen, Gruber, Mula).
   *
   * \param kernelApproximation an approximation to the scattering kernel
   * \param x  solution to which we want to apply the scattering kernel
   * \param solutionSpaces  tuple of solution spaces
   * \param gridView
   * \param accuracy
   */
  template<class Spaces, class GridView>
  static std::vector<VectorType> apply_scattering(
      ScatteringKernelApproximation& kernelApproximation,
      const std::vector<VectorType>& x,
      const Spaces& spaces,
      std::vector<Direction>& sVector,
      const GridView& gridView,
      double accuracy);
};


class TransportLogger {
  public:
  explicit TransportLogger(std::ofstream& ofs,
                           unsigned int outerIteration,
                           unsigned int direction,
                           ASTILogFlags logFlags,
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
      const size_t numDofs)
  {
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

  private:
  std::ofstream& ofs;
  const unsigned int n;
  const unsigned int i;
};

class ASTILogger {
  public:
  template<class ScatteringKernelApproximation, class RHSApproximation>
  explicit ASTILogger(ASTILogFlags logFlags,
                      std::string filename,
                      PassKey<ASTI<ScatteringKernelApproximation,
                                   RHSApproximation>>)
  : logFlags(logFlags)
  , ofs(filename) {}

  TransportLogger transportLogger(unsigned int outerIteration,
                                  unsigned int direction)
  {
    return TransportLogger(ofs, outerIteration, direction, logFlags, {});
  }

  template<class Kernel, class KernelApproximation>
  void logASTIOverview(
      const double targetAccuracy,
      const Kernel& kernel,
      const KernelApproximation& kernelApproximation,
      const ASTIApproximationParameters& approximationParameters)
  {
    ofs << "uniform ASTI algorithm\n"
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
    if(flagIsSet(logFlags, ASTILogFlags::logDirectionSolveTime)) {
      ofs << "Computing time: "
          << std::chrono::duration_cast<std::chrono::microseconds>
            (endScatteringApproximation - startScatteringApproximation).count()
          << "us\n";
    }
    ofs << std::flush;
  }

  void logInnerIterationsHeader()
  {
    ofs << "\n-----------------------------------\n"
             "Inner iterations (transport solves)\n"
             "-----------------------------------\n";
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
  const ASTILogFlags logFlags;
  std::ofstream ofs;
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

  template<class Spaces, class VectorType>
  void plotSolutions(
      const Spaces& spaces,
      const std::vector<VectorType>& x,
      const unsigned int n,
      const unsigned int numS) const
  {
    if(flagIsSet(plotFlags, ASTIPlotFlags::plotOuterIterations)) {
      //////////////////////////////////////////////////////////////////////
      //  Write result to VTK file
      //  We need to subsample, because VTK cannot natively display
      //  real second-order functions
      //////////////////////////////////////////////////////////////////////
      std::cout << "Print solutions:\n";

      const auto& feBasisInterior = spaces.interiorSolutionSpace();
      const auto& feBasisTrace    = spaces.traceSolutionSpace();


      for(unsigned int i = 0; i < numS; ++i)
      {
        std::cout << "Direction " << i << '\n';

        std::string name = outputfolder
                        + "/u_rad_trans_n"
                        + std::to_string(n)
                        + "_s"
                        + std::to_string(i);
        FunctionPlotter uPlotter(name);
        uPlotter.plot("u", x[i], feBasisInterior, 0, 0);
        name = outputfolder
                        + "/theta_rad_trans_n"
                        + std::to_string(n)
                        + "_s"
                        + std::to_string(i);
        FunctionPlotter thetaPlotter(name);
        thetaPlotter.plot("theta", x[i], feBasisTrace, 2,
                          feBasisInterior.size());
      }
    }
  }

  template<class ScatteringBasis, class VectorType>
  void plotScatteringOnHostGrid(
      const ScatteringBasis& hostGridGlobalBasis,
      const std::vector<VectorType>& scattering,
      const unsigned int n,
      const size_t numS) const
  {
    if(flagIsSet(plotFlags, ASTIPlotFlags::plotScattering)) {
      std::cout << "Plot scattering:\n";

      for(unsigned int i = 0; i < numS; ++i)
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
      const ScatteringBasis& globalBasis,
      const std::vector<VectorType>& solution,
      const KernelApproximation& kernelApproximation,
      const unsigned int n) const
  {
    if(flagIsSet(plotFlags, ASTIPlotFlags::plotIntegratedSolution)) {
      std::cout << "Plot integrated solution of outer iteration "
                << n << ":\n";

      VectorType integratedSolution(globalBasis.size());
      integratedSolution = 0.;

      const auto quadWeightsSubinterval =
          kernelApproximation.quadWeightsOfSubintervalOnCurrentLevel();
      const auto numSperInterval = kernelApproximation.numSperInterval;
      const size_t numS = solution.size();
      for(unsigned int i = 0; i < numS; ++i)
      {
        const auto quadWeight = quadWeightsSubinterval[i % numSperInterval];
        integratedSolution.axpy(quadWeight, solution[i]);
      }

      std::string name = outputfolder
                      + "/integrated_solution_n"
                      + std::to_string(n);
      FunctionPlotter scatteringPlotter(name);
      scatteringPlotter.plot("integrated solution", integratedSolution,
                             globalBasis, 0);
    }
  }

  private:
  const ASTIPlotFlags plotFlags;
  const std::string outputfolder;
};

#ifndef DOXYGEN
namespace detail {
  constexpr size_t dim = 2;
  using VectorType = BlockVector<FieldVector<double,1>>;
  using Direction = FieldVector<double, dim>;

  template<class FEBasisInterior, class Grid, class F>
  inline void approximate_rhs (
      std::vector<VectorType>& rhsFunctional,
      Grid&,
      double,
      FEBasisInterior& feBasisInterior,
      const std::vector<Direction>& sVector,
      const F& f,
      FeRHS) {
    static_assert(!is_RefinedFiniteElement<FEBasisInterior>::value,
        "Functions::interpolate won't work for refined finite elements");
    auto rhsFunction = f(feBasisInterior.gridView());
    for(auto& rhsFunctionalEntry : rhsFunctional)
    {
      VectorType gInterpolation(feBasisInterior.size());
      Functions::interpolate(feBasisInterior, gInterpolation, rhsFunction);

      auto rIt = rhsFunctionalEntry.begin();
      for(const auto& gEntry : gInterpolation) {
        *rIt += gEntry;
        ++rIt;
      }
    }
  }


  // Refine grid until accuracy kappa2*eta is reached for
  // approximation of rhs and boundary values.
  template<class FEBasisInterior, class Grid, class F>
  inline void approximate_rhs (
      std::vector<VectorType>& rhsFunctional,
      Grid& grid,
      double accuracy,
      FEBasisInterior& feBasisInterior,
      const std::vector<Direction>& sVector,
      const F& f,
      ApproximateRHS) {
    using LevelGridView = typename Grid::LevelGridView;
    using FEBasisCoarseInterior =
      changeGridView_t<FEBasisInterior, LevelGridView>;
    const size_t numS = sVector.size();
    auto gridView = grid.leafGridView();
    std::vector<VectorType> boundaryValues(numS);
    const unsigned int maxNumberOfRefinements = 3;
    for(unsigned int refinement = 1; ; refinement++) {
      static_assert(!is_RefinedFiniteElement<FEBasisInterior>::value,
          "Functions::interpolate won't work for refined finite elements");
      auto rhsFunction = f(gridView);
      for(auto& bv : boundaryValues)
      {
        VectorType gInterpolation(feBasisInterior.size());
        Functions::interpolate(feBasisInterior, gInterpolation, rhsFunction);
        std::swap(bv, gInterpolation);
      }

      double rhsError = 0.;
      for(const auto& bv : boundaryValues)
      {
        auto gApprox = Functions::makeDiscreteGlobalBasisFunction<double>(
              feBasisInterior, bv);

        auto localGExact = localFunction(rhsFunction);
        auto localGApprox = localFunction(gApprox);
        auto localView = feBasisInterior.localView();

        double rhsError_i = 0.;
        for(const auto& e : elements(gridView)) {
          localView.bind(e);
          localGExact.bind(e);
          localGApprox.bind(e);

          size_t quadratureOrder = 2*localView.tree().finiteElement()
                                              .localBasis().order()  + 4;
          const Dune::QuadratureRule<double, dim>& quad =
                Dune::QuadratureRules<double, dim>::rule(e.type(),
                                                         quadratureOrder);
          auto geometry = e.geometry();
          double local_error = 0.;
          for (const auto& quadPoint : quad) {
            const FieldVector<double,dim>& quadPos = quadPoint.position();
            const double diff = localGExact(quadPos) - localGApprox(quadPos);
            local_error += diff*diff
                         * geometry.integrationElement(quadPos)
                         * quadPoint.weight();
          }
          rhsError_i += local_error;
        }
        rhsError += std::sqrt(rhsError_i);
      }
      rhsError /= numS;
      if(rhsError <= accuracy || refinement == maxNumberOfRefinements) {
        break;
      } else {
        const auto levelGridView = grid.levelGridView(grid.maxLevel());
        FEBasisCoarseInterior coarseInteriorBasis(levelGridView);

        grid.globalRefine(1);
        feBasisInterior.update(gridView);

        std::vector<VectorType> rhsFunctionalCoarse(numS);
        std::swap(rhsFunctional, rhsFunctionalCoarse);
        for(unsigned int i = 0; i < numS; ++i)
        {
          rhsFunctional[i] = interpolateToUniformlyRefinedGrid(
              coarseInteriorBasis, feBasisInterior,
              rhsFunctionalCoarse[i]);
        }
      }
    }
    for(unsigned int i = 0; i < numS; ++i)
    {
      auto rIt = rhsFunctional[i].begin();
      for(const auto& gEntry : boundaryValues[i]) {
        *rIt += gEntry;
        ++rIt;
      }
    }
  }
} // end namespace detail
#endif

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class Grid, class F, class G, class HB, class Sigma, class Kernel>
void ASTI<ScatteringKernelApproximation, RHSApproximation>::solve(
           Grid& grid,
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
           ASTIPlotFlags plotFlags,
           ASTILogFlags logFlags) {
  static_assert(std::is_same<RHSApproximation, FeRHS>::value
      || std::is_same<RHSApproximation, ApproximateRHS>::value,
      "Unknown type provided for RHSApproximation!\n"
      "Should be either FeRHS or ApproximateRHS.");

  typedef typename Grid::LeafGridView  LeafGridView;
  LeafGridView gridView = grid.leafGridView();

  const unsigned int dim = LeafGridView::dimension;

  /////////////////////////////////////////////
  // To print information in dune-dpg/results/
  /////////////////////////////////////////////

  ASTILogger logger(logFlags, outputfolder+"/output",
      PassKey<ASTI<ScatteringKernelApproximation, RHSApproximation>>{});
  ASTIPlotter plotter(plotFlags, outputfolder);

  ////////////////////////////////////////////
  // Handle directions of discrete ordinates
  ////////////////////////////////////////////
  using Direction = FieldVector<double, dim>;

  using FEBasisTrace = Functions::BernsteinBasis<LeafGridView, 2>;
  using Spaces = TransportSpaces<LeafGridView, FEBasisTrace>;
  Spaces spaces(gridView);

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

  //////////////////////////////////
  //   Solution vectors
  //////////////////////////////////
  std::vector<VectorType> x(numS);

  for(auto& entry : x)
  {
    entry.resize(spaces.interiorSolutionSpace().size()
                 + spaces.traceSolutionSpace().size());
    entry = 0;
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

    auto startScatteringApproximation = std::chrono::steady_clock::now();

    const double kappaNorm = 1.;

    const double accuKernel = approximationParameters.scatteringAccuracy()
                                           // To prevent division by zero.
                            / (kappaNorm * ((uNorm>0.)?uNorm:1e-5));
    std::vector<VectorType> rhsFunctional =
        apply_scattering (
          kernelApproximation, x, spaces,
          sVector, gridView, accuKernel);
    numS = sVector.size();
    x.resize(numS);

    auto endScatteringApproximation = std::chrono::steady_clock::now();

    plotter.plotScatteringOnHostGrid(spaces.interiorSolutionSpace(),
                                     rhsFunctional, n, numS);

    {
      const auto& feBasisInterior = spaces.interiorSolutionSpace();

      detail::approximate_rhs (
          rhsFunctional,
          grid,
          approximationParameters.rhsAccuracy(),
          feBasisInterior,
          sVector,
          f,
          RHSApproximation{});
    }

    std::vector<bool> boundary_is_homogeneous(numS, false);
    std::transform(sVector.cbegin(), sVector.cend(),
                   boundary_is_homogeneous.begin(),
                   [&] (auto s) {
                     // TODO: write a generic test for homogeneous inflow boundary
                     return is_inflow_boundary_homogeneous(s);
                   });
    // get bv contribution to rhs
    VectorType boundaryExtension
        = harmonic_extension_of_boundary_values(g,
            spaces.traceSolutionSpace());

    logger.logKernelApproximationInfo(approximationParameters,
        kernelApproximation, accuKernel, sVector,
        startScatteringApproximation, endScatteringApproximation);

    ////////////////////////////////////////////////////
    // Inner loop
    ////////////////////////////////////////////////////
    logger.logInnerIterationsHeader();
    const double aposterioriTransportGlobal =
        compute_adaptive_transport_solution(
          x, spaces, grid, sVector, sigma, rhsFunctional,
          boundaryExtension, boundary_is_homogeneous, logger, n,
          approximationParameters.transportAccuracy(),
          kernelApproximation, maxNumberOfInnerIterations);

    plotter.plotSolutions(spaces, x, n, numS);

    // Error bound for || u_{n+1} - bar u_{n+1} ||
    deviationOfInexactIterate
      = approximationParameters.errorBetweenExactAndInexactIterate(
                                          deviationOfInexactIterate,
                                          aposterioriTransportGlobal);

    // Error bound for || u - \bar u_n || based on a priori errors
    accuracy = approximationParameters.combinedAccuracy();

    logger.logInnerIterationStats(x, aposterioriTransportGlobal,
        approximationParameters, deviationOfInexactIterate, accuracy, n);

    // compute L2 norm of u
    {
      uNorm = 0.;
      const std::vector<double> quadWeight
        = kernelApproximation.quadWeightsOfSubintervalOnCurrentLevel();
      for(size_t i=0; i<numS; ++i) {
        const double uiNorm =
          ErrorTools::l2norm(spaces.interiorSolutionSpace(), x[i]);
        uNorm += uiNorm * uiNorm
                  * quadWeight[i % kernelApproximation.numSperInterval];
      }
      uNorm = std::sqrt(uNorm);
    }

    approximationParameters.decreaseEta(uNorm);

    plotter.plotIntegratedSolution(spaces.interiorSolutionSpace(),
                                   x, kernelApproximation, n);
  }
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class Spaces, class Sigma>
double
ASTI<ScatteringKernelApproximation, RHSApproximation>::
compute_transport_solution(
    VectorType& x,
    const Spaces& spaces,
    const FieldVector<double, 2>& s,
    const Sigma sigma,
    const VectorType& rhsFunctional,
    bool boundary_is_homogeneous,
    const VectorType& bvExtension)
{
  auto bilinearForm
    = bilinearFormWithSpaces(spaces.testSpacePtr(), spaces.solutionSpacePtr())
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

  // Determine Dirichlet dofs for theta (inflow boundary)
  std::vector<bool> dirichletNodesInflow;
  BoundaryTools::getInflowBoundaryMask(
              spaces.traceSolutionSpace(),
              dirichletNodesInflow,
              s);

  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

  VectorType rhs;
  MatrixType stiffnessMatrix;

  /////////////////////////////////////////////////////////
  //  Assemble the systems
  /////////////////////////////////////////////////////////
  // loop of the discrete ordinates
  if(boundary_is_homogeneous) {
    auto rhsFunction
      = linearFormWithSpace(systemAssembler.getTestSearchSpaces())
        .template addFunctionalTerm<0>(rhsFunctional,
                                       spaces.interiorSolutionSpace())
        .create();
    systemAssembler.assembleSystem(
        stiffnessMatrix, rhs,
        rhsFunction);
  } else {
    assert(bvExtension.size() > 0);
    auto rhsFunction
      = linearFormWithSpace(systemAssembler.getTestSearchSpaces())
        .template addFunctionalTerm<0>
              (rhsFunctional, spaces.interiorSolutionSpace())
        .template addSkeletalFunctionalTerm<0, IntegrationType::normalVector>
              (bvExtension, spaces.traceSolutionSpace(), -1., s)
        .create();
    systemAssembler.assembleSystem(
        stiffnessMatrix, rhs,
        rhsFunction);
  }
  systemAssembler.template applyHomogeneousDirichletBoundary<1>
      (stiffnessMatrix,
       rhs,
       dirichletNodesInflow);
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
    auto rhsFunction
      = linearFormWithSpace(rhsAssemblerEnriched.getTestSpaces())
        .template addFunctionalTerm<0>(rhsFunctional,
                                       spaces.interiorSolutionSpace())
        .create();
    rhsAssemblerEnriched.assembleRhs(rhs, rhsFunction);
  } else {
    auto rhsAssemblerEnriched
      = make_RhsAssembler(spaces.enrichedTestSpacePtr());
    auto rhsFunction
      = linearFormWithSpace(rhsAssemblerEnriched.getTestSpaces())
        .template addFunctionalTerm<0>(rhsFunctional,
                                       spaces.interiorSolutionSpace())
        .template addSkeletalFunctionalTerm<0, IntegrationType::normalVector>
              (bvExtension, spaces.traceSolutionSpace(), -1., s)
        .create();
    rhsAssemblerEnriched.assembleRhs(rhs, rhsFunction);
  }
  // - Computation of the a posteriori error
  const double aposteriori_s
      = ErrorTools::aPosterioriError(bilinearFormEnriched,
                                     innerProductEnriched,
                                     x, rhs);

  return aposteriori_s;
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class Spaces, class Grid, class Sigma, class KernelApproximation>
double
ASTI<ScatteringKernelApproximation, RHSApproximation>::
compute_adaptive_transport_solution(
    std::vector<VectorType>& x,
    Spaces& spaces,
    Grid& grid,
    const std::vector<FieldVector<double, 2>>& sVector,
    const Sigma sigma,
    std::vector<VectorType>& rhsFunctional,
    VectorType& boundaryExtension,
    const std::vector<bool>& boundary_is_homogeneous,
    ASTILogger& logger,
    const unsigned int n,
    const double transportAccuracy,
    const KernelApproximation& kernelApproximation,
    const unsigned int maxNumberOfInnerIterations)
{
  const size_t numS = sVector.size();
  double aposterioriTransportGlobal = 0.;

  for(unsigned int nRefinement = 0; ; )
    // At the end of the loop, we will break if
    // aposterioriTransportGlobal < kapp3*eta
    // or ++nRefinement >= maxNumberOfInnerIterations
    // thus the inner loop terminates eventually.
  {
    aposterioriTransportGlobal = 0.;
    for(unsigned int i = 0; i < numS; ++i)
    {
      auto transportLogger = logger.transportLogger(n, i);
      transportLogger.printCurrentIteration(nRefinement);

      const double aposteriori_s
          = compute_transport_solution(x[i], spaces,
              sVector[i], sigma, rhsFunctional[i],
              boundary_is_homogeneous[i],
              boundaryExtension);

      const std::vector<double> quadWeight
        = kernelApproximation.quadWeightsOfSubintervalOnCurrentLevel();

      aposterioriTransportGlobal
        +=  quadWeight[i % kernelApproximation.numSperInterval]
          * aposteriori_s;

      // TODO: Add (interpolation of) g to theta part of x?

      transportLogger.logSolverStats(nRefinement, std::sqrt(aposteriori_s),
          grid.maxLevel(), x[i].size());
    }
    aposterioriTransportGlobal = std::sqrt(aposterioriTransportGlobal);

    if(++nRefinement >= maxNumberOfInnerIterations
        || aposterioriTransportGlobal <= transportAccuracy) {
      break;
    } else {
      grid.globalRefine(1);
      spaces.update(grid.leafGridView());

      using LevelGridView = typename Grid::LevelGridView;
      const LevelGridView levelGridView
          = grid.levelGridView(grid.maxLevel()-1);

      {
        using FEBasisCoarseInterior
            = changeGridView_t<typename Spaces::FEBasisInterior,
                               LevelGridView>;

        FEBasisCoarseInterior coarseInteriorBasis(levelGridView);
        const auto& feBasisInterior = spaces.interiorSolutionSpace();

        std::vector<VectorType> rhsFunctionalCoarse(numS);
        std::swap(rhsFunctional, rhsFunctionalCoarse);

        std::transform(rhsFunctionalCoarse.cbegin(),
                       rhsFunctionalCoarse.cend(),
                       rhsFunctional.begin(),
                       [&] (const VectorType& rhsCoarse) {
                         return interpolateToUniformlyRefinedGrid(
                                coarseInteriorBasis, feBasisInterior,
                                rhsCoarse);
                       });
      }
      {
        using FEBasisCoarseTrace
            = changeGridView_t<typename Spaces::FEBasisTrace, LevelGridView>;

        FEBasisCoarseTrace coarseTraceBasis(levelGridView);
        const auto& feBasisTrace = spaces.traceSolutionSpace();

        VectorType boundaryExtensionFine
            = interpolateToUniformlyRefinedGrid(
                  coarseTraceBasis, feBasisTrace,
                  boundaryExtension);
        std::swap(boundaryExtension, boundaryExtensionFine);
      }
    }
  }
  return aposterioriTransportGlobal;
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class Spaces, class GridView>
std::vector<Dune::BlockVector<Dune::FieldVector<double, 1> >>
ASTI<ScatteringKernelApproximation, RHSApproximation>::apply_scattering(
      ScatteringKernelApproximation& kernelApproximation,
      const std::vector<VectorType>& x,
      const Spaces& spaces,
      std::vector<Direction>& sVector,
      const GridView& gridView,
      double accuracy) {
  sVector = kernelApproximation.setAccuracyAndInputSize(accuracy, x.size());

  const size_t numS = sVector.size();
  std::vector<VectorType> rhsFunctional(numS);
  const auto& feBasisInterior = spaces.interiorSolutionSpace();

  auto scatteringAssembler =
      make_ApproximateScatteringAssembler(feBasisInterior,
                                          kernelApproximation);
  scatteringAssembler.computeScattering(rhsFunctional, x);

  return rhsFunctional;
}

} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_ASTI_UNIFORM_HH
