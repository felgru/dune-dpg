// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_PERITER_UNIFORM_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_PERITER_UNIFORM_HH

#include <chrono>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <string>
#include <vector>

#include <dune/common/version.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/io.hh>
#include <dune/istl/umfpack.hh>

#include <dune/functions/functionspacebases/bernsteinbasis.hh>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/errortools.hh>
#include <dune/dpg/functionplotter.hh>
#include <dune/dpg/functions/interpolate.hh>
#include <dune/dpg/linearfunctionalterm.hh>
#include <dune/dpg/radiative_transfer/approximate_scattering.hh>
#include <dune/dpg/radiative_transfer/boundary_extension.hh>
#include <dune/dpg/radiative_transfer/passkey.hh>
#include <dune/dpg/radiative_transfer/periter_common.hh>
#include <dune/dpg/rhs_assembler.hh>
#include <dune/dpg/dpg_system_assembler.hh>
#include <dune/dpg/type_traits.hh>

#include <boost/math/constants/constants.hpp>

namespace Dune {

class PeriterLogger;
class TransportLogger;

/**
 * This class describes the Periter algorithm for radiative transfer problems
 *
 * \tparam ScatteringKernelApproximation
 *         specifies the method used to approximate the scattering kernel
 * \tparam RHSApproximation  if right hand side and lifting of boundary
 *                           values are finite element functions, set this
 *                           to FeRHS, otherwise set this to
 *                           ApproximateRHS
 */
template<class ScatteringKernelApproximation, class RHSApproximation>
class Periter {
  public:

  using Direction = typename ScatteringKernelApproximation::Direction;

  /**
   * Solve a radiative transfer problem using the Periter algorithm
   *
   * \param grid
   * \param f  right hand side function
   * \param g  lifting of the boundary values
   * \param is_inflow_boundary_homogeneous
   *            checks if g is 0 on the inflow boundary
   * \param sigma   absorbtion coefficient
   * \param kernel  the scattering kernel, e.g. a Henyey–Greenstein kernel
   * \param rho  the contraction parameter ρ
   * \param CT  the constant C_T from the paper
   * \param cB  the inf-sup constant of the operator B = T - K
   * \param targetAccuracy  periter solves up to this accuracy
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
             double rho,
             double CT,
             double cB,
             double targetAccuracy,
             unsigned int maxNumberOfIterations,
             unsigned int maxNumberOfInnerIterations,
             const std::string& outputfolder,
             PeriterPlotFlags plotFlags = PeriterPlotFlags::doNotPlot);

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
   * \param sigma the absorbtion coefficient
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
      PeriterLogger& logger,
      unsigned int n,
      double transportAccuracy,
      const KernelApproximation& kernelApproximation,
      const std::vector<double>& quadWeight,
      unsigned int maxNumberOfInnerIterations);

  /**
   * Apply the scattering integral to a solution x
   *
   * This corresponds to [K, u_n, κ_1 * η] in the Periter algorithm
   * (see Dahmen, Gruber, Mula).
   *
   * \param kernelApproximation an approximation to the scattering kernel
   * \param x  solution to which we want to apply the scattering kernel
   * \param solutionSpaces  tuple of solution spaces
   * \param gridView
   * \param accuracy
   */
  template<class SolutionSpaces, class GridView>
  static std::vector<VectorType> apply_scattering(
      ScatteringKernelApproximation& kernelApproximation,
      const std::vector<VectorType>& x,
      const SolutionSpaces& solutionSpaces,
      std::vector<Direction>& sVector,
      const GridView& gridView,
      double accuracy);
};


class TransportLogger {
  public:
  explicit TransportLogger(std::ofstream& ofs,
                           unsigned int outerIteration,
                           unsigned int direction,
                           PassKey<PeriterLogger>)
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
        << " for direction " << i << ": \n"
        << "  - A posteriori estimation of || (u,trace u) - (u_fem,theta) || = "
        << aposteriori_s
        << "\n  - Grid level: "     << maxGridLevel
        << "\n  - Number of DOFs: " << numDofs
        << std::endl;

    std::cout << "\nIteration " << n << '.' << nRefinement
        << " for direction " << i << ": \n"
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

class PeriterLogger {
  public:
  template<class ScatteringKernelApproximation, class RHSApproximation>
  explicit PeriterLogger(std::string filename,
                         PassKey<Periter<ScatteringKernelApproximation,
                                         RHSApproximation>>) : ofs(filename) {}

  TransportLogger transportLogger(unsigned int outerIteration,
                                  unsigned int direction)
  {
    return TransportLogger(ofs, outerIteration, direction, {});
  }

  template<class Kernel, class KernelApproximation>
  void logPeriterOverview(
      const double targetAccuracy,
      const Kernel& kernel,
      const KernelApproximation& kernelApproximation,
      const detail::ApproximationParameters& approximationParameters,
      const double CT)
  {
    ofs << "uniform PERITER algorithm\n"
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
      const std::vector<double> singularValues
        = kernelApproximation.getSingularValues();
      ofs << "Singular values of kernel matrix:\n";
      for(size_t i=0; i<singularValues.size(); i++){
        ofs << singularValues[i] << '\n';
      }
    }
  }

  void logOuterIterationHeader(const unsigned int n)
  {
    ofs << "\nIteration n=" << n
        << "\n================\n";
    std::cout << "\nIteration " << n << "\n\n";
  }

  template<class KernelApproximation>
  void logKernelApproximationInfo(
      const detail::ApproximationParameters& approximationParameters,
      const KernelApproximation& kernelApproximation,
      const double accuKernel,
      const std::vector<FieldVector<double, 2>>& sVector,
      const std::chrono::steady_clock::time_point startScatteringApproximation,
      const std::chrono::steady_clock::time_point endScatteringApproximation)
  {
    ofs << "eta_n = rhobar^{-n}: " << approximationParameters.eta() << '\n'
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
        << kernelApproximation.info()      << '\n'
        << "Computing time: "
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

  template<class VectorType>
  void logInnerIterationStats(
      const std::vector<VectorType>& x,
      const double aposterioriTransportGlobal,
      const detail::ApproximationParameters& approximationParameters,
      const std::vector<double>& aposterioriIter,
      const double accuracy,
      const unsigned int n)
  {
    const size_t accumulatedDoFs = std::accumulate(x.cbegin(), x.cend(),
        static_cast<size_t>(0),
        [](size_t acc, auto vec) { return acc + vec.size(); });

    // Error bound for || u_n - \bar u_n || based on a posteriori errors
    const double errorAPosteriori
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
  }

  private:
  std::ofstream ofs;
};

class PeriterPlotter {
  public:
  PeriterPlotter(PeriterPlotFlags plotFlags, std::string outputfolder)
    : plotFlags(plotFlags)
    , outputfolder(outputfolder)
  {
    if((plotFlags & PeriterPlotFlags::plotLastIteration)
        == PeriterPlotFlags::plotLastIteration) {
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
    if(plotOuterIterationsEnabled()) {
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
                        + std::string("/u_rad_trans_n")
                        + std::to_string(n)
                        + std::string("_s")
                        + std::to_string(i);
        FunctionPlotter uPlotter(name);
        uPlotter.plot("u", x[i], feBasisInterior, 0, 0);
        name = outputfolder
                        + std::string("/theta_rad_trans_n")
                        + std::to_string(n)
                        + std::string("_s")
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
    if(plotScatteringEnabled()) {
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
        scatteringPlotter.plot("scattering", scattering[i],
                               hostGridGlobalBasis, 0);
      }
    }
  }

  private:
  bool plotOuterIterationsEnabled() const {
    return (plotFlags & PeriterPlotFlags::plotOuterIterations)
           == PeriterPlotFlags::plotOuterIterations;
  }

  bool plotScatteringEnabled() const {
    return (plotFlags & PeriterPlotFlags::plotScattering)
           == PeriterPlotFlags::plotScattering;
  }

  const PeriterPlotFlags plotFlags;
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
    const size_t numS = sVector.size();
    auto rhsFunction = f(feBasisInterior.gridView());
    for(unsigned int i = 0; i < numS; ++i)
    {
      VectorType gInterpolation(feBasisInterior.size());
      Functions::interpolate(feBasisInterior, gInterpolation, rhsFunction);

      // Add gInterpolate to first feBasisInterior.size() entries of
      // rhsFunctional[i].
      using Iterator = std::decay_t<decltype(rhsFunctional[i].begin())>;
      for(Iterator rIt=rhsFunctional[i].begin(),
                   rEnd=rhsFunctional[i].begin()+feBasisInterior.size(),
                   gIt=gInterpolation.begin(); rIt!=rEnd; ++rIt, ++gIt) {
        *rIt += *gIt;
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
      for(unsigned int i = 0; i < numS; ++i)
      {
        VectorType gInterpolation(feBasisInterior.size());
        Functions::interpolate(feBasisInterior, gInterpolation, rhsFunction);
        std::swap(boundaryValues[i], gInterpolation);
      }

      double rhsError = 0.;
      for(unsigned int i = 0; i < numS; ++i)
      {
        auto gApprox = Functions::makeDiscreteGlobalBasisFunction<double>(
              feBasisInterior, boundaryValues[i]);

        auto localGExact = localFunction(rhsFunction);
        auto localGApprox = localFunction(gApprox);
        auto localView = feBasisInterior.localView();
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
        auto localIndexSet = feBasisInterior.localIndexSet();
#endif

        double rhsError_i = 0.;
        for(const auto& e : elements(gridView)) {
          localView.bind(e);
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
          localIndexSet.bind(localView);
#endif
          localGExact.bind(e);
          localGApprox.bind(e);

          size_t quadratureOrder = 2*localView.tree().finiteElement()
                                              .localBasis().order()  + 4;
          const Dune::QuadratureRule<double, dim>& quad =
                Dune::QuadratureRules<double, dim>::rule(e.type(),
                                                         quadratureOrder);
          auto geometry = e.geometry();
          double local_error = 0.;
          for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {
            const FieldVector<double,dim>& quadPos = quad[pt].position();
            const double diff = localGExact(quadPos) - localGApprox(quadPos);
            local_error += diff*diff
                         * geometry.integrationElement(quadPos)
                         * quad[pt].weight();
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
      FEBasisInterior feBasisInterior(gridView);
      // Add boundaryValues[i] to first feBasisInterior.size() entries of
      // rhsFunctional[i].
      using Iterator = std::decay_t<decltype(rhsFunctional[i].begin())>;
      for(Iterator rIt=rhsFunctional[i].begin(),
                   rEnd=rhsFunctional[i].begin()+feBasisInterior.size(),
                   gIt=boundaryValues[i].begin();
          rIt!=rEnd; ++rIt, ++gIt) {
        *rIt += *gIt;
      }
    }
  }
} // end namespace detail
#endif

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class Grid, class F, class G, class HB, class Sigma, class Kernel>
void Periter<ScatteringKernelApproximation, RHSApproximation>::solve(
           Grid& grid,
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
           PeriterPlotFlags plotFlags) {
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

  PeriterLogger logger(outputfolder+"/output",
      PassKey<Periter<ScatteringKernelApproximation, RHSApproximation>>{});
  PeriterPlotter plotter(plotFlags, outputfolder);

  // TODO: estimate norm of rhs f
  const double fnorm = 1;
  const double err0 = fnorm / cB;
  detail::ApproximationParameters approximationParameters(rho, CT, err0,
                                                          RHSApproximation{});

  ////////////////////////////////////////////
  // Handle directions of discrete ordinates
  ////////////////////////////////////////////
  using Direction = FieldVector<double, dim>;

  using FEBasisTrace = Functions::BernsteinBasis<LeafGridView, 2>;
  using Spaces = TransportSpaces<LeafGridView, FEBasisTrace>;
  Spaces spaces(gridView);

  ScatteringKernelApproximation kernelApproximation(kernel,
      approximationParameters.finalScatteringAccuracy(targetAccuracy));

  logger.logPeriterOverview(targetAccuracy, kernel,
      kernelApproximation, approximationParameters, CT);

  // As the solution u we use for the initial scattering is 0, and the
  // formula for the accuracy contains a 1/\|u\|, we set the initial
  // accuracy to a large enough value.
  std::vector<Direction>
    sVector(kernelApproximation.setAccuracyAndInputSize(1e5, 0));
  size_t numS = sVector.size();

  //////////////////////////////////
  //   right hand side vector type
  //////////////////////////////////
  typedef BlockVector<FieldVector<double,1> > VectorType;

  //////////////////////////////////
  //   Solution vectors
  //////////////////////////////////
  std::vector<VectorType> x(numS);

  for(unsigned int i = 0; i < numS; ++i)
  {
    x[i].resize(spaces.interiorSolutionSpace().size()
                 + spaces.traceSolutionSpace().size());
    x[i] = 0;
  }

  /////////////////////////////////////////////////////////
  //  Fixed-point iterations
  /////////////////////////////////////////////////////////
  double accuracy = err0;

  std::vector<double> aposterioriIter(maxNumberOfIterations, 0.);

  for(unsigned int n = 0; accuracy > targetAccuracy
                          && n < maxNumberOfIterations; ++n)
  {
    logger.logOuterIterationHeader(n);

    auto startScatteringApproximation = std::chrono::steady_clock::now();

    const std::vector<double> quadWeight
      = kernelApproximation.getQuadWeightSubinterval();

    double kappaNorm = 1.;
    double uNorm = 0.;
    for(size_t i=0; i<numS; ++i) {
      const double uiNorm =
        ErrorTools::l2norm(spaces.interiorSolutionSpace(), x[i]);
      uNorm += uiNorm * uiNorm
                / (1 << kernelApproximation.getLevel())
                * quadWeight[i % kernelApproximation.numSperInterval];
    }
    uNorm = std::sqrt(uNorm);
    // To prevent division by zero.
    if(uNorm == 0.) uNorm = 1e-5;

    const double accuKernel = approximationParameters.scatteringAccuracy()
                            / (kappaNorm * uNorm);
    std::vector<VectorType> rhsFunctional =
        apply_scattering (
          kernelApproximation, x, *spaces.solutionSpacePtr(),
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
    for(size_t i = 0; i < numS; i++) {
      boundary_is_homogeneous[i]
          = is_inflow_boundary_homogeneous(sVector[i]);
      // TODO: write a generic test for homogeneous inflow boundary
    }
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
          kernelApproximation, quadWeight,
          maxNumberOfInnerIterations);

    plotter.plotSolutions(spaces, x, n, numS);

    // A posteriori estimation of error ||bar u_n -T^{-1}K bar u_{n-1}||
    aposterioriIter[n] = aposterioriTransportGlobal
                       + CT * approximationParameters.scatteringAccuracy();

    // Error bound for || u - \bar u_n || based on a priori errors
    accuracy = approximationParameters.combinedAccuracy();

    logger.logInnerIterationStats(x, aposterioriTransportGlobal,
        approximationParameters, aposterioriIter, accuracy, n);

    approximationParameters.decreaseEta();
  }
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class Spaces, class Sigma>
double
Periter<ScatteringKernelApproximation, RHSApproximation>::
compute_transport_solution(
    VectorType& x,
    const Spaces& spaces,
    const FieldVector<double, 2>& s,
    const Sigma sigma,
    const VectorType& rhsFunctional,
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


  using LeafGridView = typename  Spaces::GridView;
  using Geometry = typename LeafGridView::template Codim<0>::Geometry;
  using GeometryBuffer_t = GeometryBuffer<Geometry>;

  GeometryBuffer_t geometryBuffer;
  auto systemAssembler =
      make_DPGSystemAssembler(bilinearForm,
                              innerProduct,
                              geometryBuffer);

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
    auto rhsFunction = make_LinearForm(
          systemAssembler.getTestSearchSpaces(),
          std::make_tuple(
            make_LinearFunctionalTerm<0>
              (rhsFunctional, spaces.interiorSolutionSpace())));
    systemAssembler.assembleSystem(
        stiffnessMatrix, rhs,
        rhsFunction);
  } else {
    assert(bvExtension.size() > 0);
    auto rhsFunction = make_LinearForm(
          systemAssembler.getTestSearchSpaces(),
          std::make_tuple(
            make_LinearFunctionalTerm<0>
              (rhsFunctional, spaces.interiorSolutionSpace()),
            make_SkeletalLinearFunctionalTerm
              <0, IntegrationType::normalVector>
              (bvExtension, spaces.traceSolutionSpace(), -1, s)));
    systemAssembler.assembleSystem(
        stiffnessMatrix, rhs,
        rhsFunction);
  }
  systemAssembler.template applyDirichletBoundary<1>
      (stiffnessMatrix,
      rhs,
      dirichletNodesInflow,
      0.);
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
               std::get<0>(*systemAssembler.getSolutionSpaces()))));
    rhsAssemblerEnriched.assembleRhs(rhs, rhsFunction);
  } else {
    auto rhsAssemblerEnriched
      = make_RhsAssembler(spaces.enrichedTestSpacePtr());
    auto rhsFunction = make_LinearForm(
          rhsAssemblerEnriched.getTestSpaces(),
          std::make_tuple(
            make_LinearFunctionalTerm<0>
              (rhsFunctional,
               std::get<0>(*systemAssembler.getSolutionSpaces())),
            make_SkeletalLinearFunctionalTerm
              <0, IntegrationType::normalVector>
              (bvExtension,
               std::get<1>(*systemAssembler.getSolutionSpaces()),
               -1, s)));
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
Periter<ScatteringKernelApproximation, RHSApproximation>::
compute_adaptive_transport_solution(
    std::vector<VectorType>& x,
    Spaces& spaces,
    Grid& grid,
    const std::vector<FieldVector<double, 2>>& sVector,
    const Sigma sigma,
    std::vector<VectorType>& rhsFunctional,
    VectorType& boundaryExtension,
    const std::vector<bool>& boundary_is_homogeneous,
    PeriterLogger& logger,
    const unsigned int n,
    const double transportAccuracy,
    const KernelApproximation& kernelApproximation,
    const std::vector<double>& quadWeight,
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

      aposterioriTransportGlobal
        += 1. / (1 << kernelApproximation.getLevel())
          * quadWeight[i % kernelApproximation.numSperInterval]
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

        for(unsigned int i = 0; i < numS; ++i)
        {
          rhsFunctional[i] = interpolateToUniformlyRefinedGrid(
              coarseInteriorBasis, feBasisInterior,
              rhsFunctionalCoarse[i]);
        }
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
template<class SolutionSpaces, class GridView>
std::vector<Dune::BlockVector<Dune::FieldVector<double, 1> >>
Periter<ScatteringKernelApproximation, RHSApproximation>::apply_scattering(
      ScatteringKernelApproximation& kernelApproximation,
      const std::vector<VectorType>& x,
      const SolutionSpaces& solutionSpaces,
      std::vector<Direction>& sVector,
      const GridView& gridView,
      double accuracy) {
  sVector = kernelApproximation.setAccuracyAndInputSize(accuracy, x.size());

  using FEBasisInterior = std::tuple_element_t<0, SolutionSpaces>;

  const size_t numS = sVector.size();
  std::vector<VectorType> rhsFunctional(numS);
  const FEBasisInterior& feBasisInterior = std::get<0>(solutionSpaces);

  auto scatteringAssembler =
      make_ApproximateScatteringAssembler(feBasisInterior,
                                          kernelApproximation);
  scatteringAssembler.computeScattering(rhsFunctional, x);

  return rhsFunctional;
}

} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_PERITER_UNIFORM_HH
