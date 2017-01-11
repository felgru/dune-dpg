// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_PERITER_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_PERITER_HH

#include <chrono>
#include <vector>

#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/io.hh>
#include <dune/istl/umfpack.hh>

#include <dune/functions/functionspacebases/pqkdgrefineddgnodalbasis.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pqktracenodalbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/pqksubsampleddgbasis.hh>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/errortools.hh>
#include <dune/dpg/functionplotter.hh>
#include <dune/dpg/functions/interpolate.hh>
#include <dune/dpg/linearfunctionalterm.hh>
#include <dune/dpg/radiative_transfer/approximate_scattering.hh>
#include <dune/dpg/rhs_assembler.hh>
#include <dune/dpg/dpg_system_assembler.hh>
#include <dune/dpg/type_traits.hh>

#include <boost/math/constants/constants.hpp>

namespace Dune {

enum class PlotSolutions {
  doNotPlot,
  plotOuterIterations,
  plotLastIteration
};

/**
 * This class describes the Periter algorithm for radiative transfer problems
 *
 * \tparam ScatteringKernelApproximation
 *         specifies the method used to approximate the scattering kernel
 * \tparam RHSApproximation  if right hand side and lifting of boundary
 *                           values are finite element functions, set this
 *                           to FeRHSandBoundary, otherwise set this to
 *                           ApproximateRHSandBoundary
 */
template<class ScatteringKernelApproximation, class RHSApproximation>
class Periter {
  public:

  /**
   * Solve a radiative transfer problem using the Periter algorithm
   *
   * \param createGrids
   * \param f  right hand side function
   * \param g  lifting of the boundary values
   * \param gDeriv  derivative of g in direction s
   * \param sigma
   * \param kernel  the scattering kernel, e.g. a Henyey–Greenstein kernel
   * \param numS  number of directions used in the discretization
   * \param rho  the contraction parameter ρ
   * \param CT  the constant C_T from the paper
   * \param targetAccuracy  periter solves up to this accuracy
   * \param maxNumberOfIterations  ... or up to the given number of iterations
   *                               (whatever comes first)
   * \param plotSolutions  specifies when to create .vtu files for plotting
   *                       the solution
   */
  template<class GridsConstructor, class F, class G, class GDeriv,
           class Kernel>
  void solve(const GridsConstructor& createGrids,
             const F& f,
             const G& g,
             const GDeriv& gDeriv,
             double sigma,
             const Kernel& kernel,
             unsigned int numS,
             double rho,
             double CT,
             double targetAccuracy,
             unsigned int maxNumberOfIterations,
             PlotSolutions plotSolutions = PlotSolutions::doNotPlot);

  private:
  using VectorType = BlockVector<FieldVector<double,1>>;

  /**
   * Apply the scattering integral to a solution x
   *
   * This corresponds to [K, u_n, κ_1 * η] in the Periter algorithm
   * (see Dahmen, Gruber, Mula).
   *
   * \tparam FEBasisInterior  type of the space for inner DOFs
   * \tparam FEBasisTrace     type of the space for trace DOFs
   * \param kernelApproximation an approximation to the scattering kernel
   * \param x  solution to which we want to apply the scattering kernel
   * \param gridViews
   * \param accuracy
   */
  template<class FEBasisInterior, class FEBasisTrace, class GridView>
  static std::vector<VectorType> apply_scattering(
      ScatteringKernelApproximation& kernelApproximation,
      const std::vector<VectorType>& x,
      const std::vector<GridView>& gridViews,
      double accuracy);
};

struct FeRHSandBoundary {};
struct ApproximateRHSandBoundary {};

#ifndef DOXYGEN
namespace detail {
  constexpr size_t dim = 2;
  using VectorType = BlockVector<FieldVector<double,1>>;
  using Direction = FieldVector<double, dim>;

  template<class FEBasisInterior, class Grids,
           class F, class G, class GDeriv>
  inline void approximate_rhs_and_bv (
      std::vector<VectorType>& rhsFunctional,
      Grids& grids,
      double,
      const std::vector<Direction>& sVector,
      F& f, G& g, GDeriv& gDeriv, double sigma,
      FeRHSandBoundary) {
    static_assert(!is_RefinedFiniteElement<FEBasisInterior>::value,
        "Functions::interpolate won't work for refined finite elements");
    size_t numS = sVector.size();
    for(unsigned int i = 0; i < numS; ++i)
    {
      FEBasisInterior feBasisInterior(grids[i]->leafGridView());
      VectorType gInterpolation(feBasisInterior.size());
      const Direction s = sVector[i];
      Functions::interpolate(feBasisInterior, gInterpolation,
          [s,&f,&g,&gDeriv,sigma](const Direction& x)
          { return f(x,s) + gDeriv(x,s) - sigma*g(x,s); });

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
  template<class FEBasisInterior, class Grids,
           class F, class G, class GDeriv>
  inline void approximate_rhs_and_bv (
      std::vector<VectorType>& rhsFunctional,
      Grids& grids,
      double accuracy,
      const std::vector<Direction>& sVector,
      F& f, G& g, GDeriv& gDeriv, double sigma,
      ApproximateRHSandBoundary) {
    using Grid = typename Grids::value_type::element_type;
    using LevelGridView = typename Grid::LevelGridView;
    using FEBasisCoarseInterior =
      changeGridView_t<FEBasisInterior, LevelGridView>;
    size_t numS = sVector.size();

    std::vector<FEBasisInterior> feBasesInterior;
    feBasesInterior.reserve(numS);
    for(unsigned int i = 0; i < numS; ++i)
      feBasesInterior.emplace_back(grids[i]->leafGridView());

    std::vector<VectorType> boundaryValues(numS);
    const unsigned int maxNumberOfRefinements = 3;
    for(unsigned int refinement = 1; ; refinement++) {
      static_assert(!is_RefinedFiniteElement<FEBasisInterior>::value,
          "Functions::interpolate won't work for refined finite elements");
      for(unsigned int i = 0; i < numS; ++i)
      {
        VectorType gInterpolation(feBasesInterior[i].size());
        const Direction s = sVector[i];
        Functions::interpolate(feBasesInterior[i], gInterpolation,
            [s,&f,&g,&gDeriv,sigma](const Direction& x)
            { return f(x,s) + gDeriv(x,s) - sigma*g(x,s); });

        std::swap(boundaryValues[i], gInterpolation);
      }

      double rhsError = 0.;
      for(unsigned int i = 0; i < numS; ++i)
      {
        const Direction s = sVector[i];
        auto gExact = Functions::makeGridViewFunction(
              [s,&f,&g,&gDeriv,sigma](const Direction& x)
              { return f(x,s) + gDeriv(x,s) - sigma*g(x,s); },
              grids[i]->leafGridView());
        auto gApprox = Functions::makeDiscreteGlobalBasisFunction<double>(
              feBasesInterior[i], boundaryValues[i]);

        auto localGExact = localFunction(gExact);
        auto localGApprox = localFunction(gApprox);
        auto localView = feBasesInterior[i].localView();
        auto localIndexSet = feBasesInterior[i].localIndexSet();

        double rhsError_i = 0.;
        for(const auto& e : elements(grids[i]->leafGridView())) {
          localView.bind(e);
          localIndexSet.bind(localView);
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
        std::vector<VectorType> rhsFunctionalCoarse(numS);
        std::swap(rhsFunctional, rhsFunctionalCoarse);
        for(unsigned int i = 0; i < numS; ++i)
        {
          const auto levelGridView =
              grids[i]->levelGridView(grids[i]->maxLevel());
          FEBasisCoarseInterior coarseInteriorBasis(levelGridView);

          grids[i]->globalRefine(1);
          feBasesInterior[i].update(grids[i]->leafGridView());

          rhsFunctional[i] = interpolateToUniformlyRefinedGrid(
              coarseInteriorBasis, feBasesInterior[i],
              rhsFunctionalCoarse[i]);
        }
      }
    }
    for(unsigned int i = 0; i < numS; ++i)
    {
      FEBasisInterior feBasisInterior(grids[i]->leafGridView());
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
template<class GridsConstructor, class F, class G, class GDeriv,
         class Kernel>
void Periter<ScatteringKernelApproximation, RHSApproximation>::solve(
           const GridsConstructor& createGrids,
           const F& f,
           const G& g,
           const GDeriv& gDeriv,
           double sigma,
           const Kernel& kernel,
           unsigned int numS,
           double rho,
           double CT,
           double targetAccuracy,
           unsigned int maxNumberOfIterations,
           PlotSolutions plotSolutions) {
  if(plotSolutions == PlotSolutions::plotLastIteration) {
    std::cerr
        << "Plotting of only the last iteration is not implemented yet!\n";
    std::abort();
  }
  static_assert(std::is_same<RHSApproximation, FeRHSandBoundary>::value
      || std::is_same<RHSApproximation, ApproximateRHSandBoundary>::value,
      "Unknown type provided for RHSApproximation!\n"
      "Should be either FeRHSandBoundary or ApproximateRHSandBoundary.");
  constexpr bool rhsIsFeFunction
      = std::is_same<RHSApproximation, FeRHSandBoundary>::value;

  const unsigned int dim = 2;
  using Direction = FieldVector<double, dim>;
  using Grid = typename
    decltype(createGrids(std::declval<std::vector<Direction>>()))::
      value_type::element_type;
  using LeafGridView  = typename Grid::LeafGridView;
  using LevelGridView = typename Grid::LevelGridView;
  using Geometry = typename LeafGridView::template Codim<0>::Geometry;
  using Domain = typename Geometry::GlobalCoordinate;

  ///////////////////////////////////
  // To print information
  ///////////////////////////////////
  std::ofstream ofs("output_rad_trans");

  ////////////////////////////////////////////
  // Handle directions of discrete ordinates
  ////////////////////////////////////////////
  // Vector of directions: sVector
  std::vector< Direction > sVector(numS);
  for(unsigned int i = 0; i < numS; ++i)
  {
    using namespace boost::math::constants;
    sVector[i] = {cos(2*pi<double>()*i/numS),
                  sin(2*pi<double>()*i/numS)};
  }

  auto grids = createGrids(sVector);
  std::vector<LeafGridView> gridViews;
  gridViews.reserve(grids.size());
  std::for_each(grids.begin(), grids.end(),
      [&gridViews](const std::shared_ptr<Grid>& grid) {
        gridViews.push_back(grid->leafGridView());
      });

  /////////////////////////////////////////////////////////
  //   Choose finite element spaces for the solution
  /////////////////////////////////////////////////////////

  typedef Functions::LagrangeDGBasis<LeafGridView, 1> FEBasisInterior; // u
  typedef Functions::PQkTraceNodalBasis<LeafGridView, 2> FEBasisTrace; // u^

  typedef changeGridView_t<FEBasisInterior, LevelGridView>
          FEBasisCoarseInterior;

  /////////////////////////////////////////////////////////
  //   Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////
  typedef BlockVector<FieldVector<double,1> > VectorType;
  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

  /////////////////////////////////////////////////
  //   Solution vectors
  /////////////////////////////////////////////////
  std::vector<VectorType> x(numS);

  for(unsigned int i = 0; i < numS; ++i)
  {
    FEBasisInterior feBasisInterior(gridViews[i]);
    FEBasisTrace feBasisTrace(gridViews[i]);

    x[i].resize(feBasisTrace.size()
               +feBasisInterior.size());
    x[i] = 0;
  }

  ScatteringKernelApproximation kernelApproximation(kernel, numS);

  /////////////////////////////////////////////////////////
  //  Fixed-point iterations
  /////////////////////////////////////////////////////////
  // TODO: A priori estimate for the accuracy of our solution:
  double accuracy = 1.;
  // η_n:
  double eta = 1;
  // TODO: estimate norm of rhs f
  const double fnorm = 1;
  // ρ̄:
  const double rhobar = (1./rho > 2*rho)? (1./rho) : (2*rho);

  // CT*kappa1 + (1+CT)*kappa2 + 2*kappa3 = 1.
  const double kappa1 = rhsIsFeFunction? 1./(2.*CT) : 1./(3.*CT);
  const double kappa2 = rhsIsFeFunction? 0.         : 1./(3.*(1+CT));
  const double kappa3 = rhsIsFeFunction? 1./4.      : 1./6.;

  ofs << "Periter with " << numS
      << " directions, rho = " << rho << ", CT = " << CT
      << ", kappa1 = " << kappa1
      << ", kappa2 = " << kappa2
      << ", kappa3 = " << kappa3
      << std::endl;

  for(unsigned int n = 0; accuracy > targetAccuracy
                          && n < maxNumberOfIterations; ++n)
  {
    ofs << "\nIteration " << n << std::endl;
    std::cout << "\nIteration " << n << std::endl << std::endl;

    std::chrono::steady_clock::time_point startScatteringApproximation
        = std::chrono::steady_clock::now();

    std::vector<VectorType> rhsFunctional =
        apply_scattering<FEBasisInterior, FEBasisTrace> (
          kernelApproximation, x, gridViews, kappa1*eta);

    std::chrono::steady_clock::time_point endScatteringApproximation
        = std::chrono::steady_clock::now();

    detail::approximate_rhs_and_bv<FEBasisInterior> (
        rhsFunctional,
        grids,
        kappa2*eta,
        sVector,
        f, g, gDeriv, sigma, RHSApproximation{});


    ////////////////////////////////////////////////////
    // Inner loop
    ////////////////////////////////////////////////////
    const unsigned int maxNumberOfInnerIterations = 3;
    double aposterioriErr;
    for(unsigned int nRefinement = 0; ; )
        // At the end of the loop, we will break if
        // aposterioriErr < kapp3*eta
        // or ++nRefinement >= maxNumberOfInnerIterations
        // thus the inner loop terminates eventually.
    {
      std::cout << "Inner iteration " << nRefinement << std::endl;

      using SolutionSpacesPtr
          = std::shared_ptr<std::tuple<FEBasisInterior, FEBasisTrace>>;

      using FEBasisTest = Functions::LagrangeDGBasis<LeafGridView, 4>;
      using FEBasisTestEnriched = FEBasisTest;

      using TestSpacesPtr         = std::shared_ptr<std::tuple<FEBasisTest>>;
      using TestSpacesEnriched    = std::tuple<FEBasisTestEnriched>;
      using TestSpacesEnrichedPtr = std::shared_ptr<TestSpacesEnriched>;

      typedef decltype(make_BilinearForm(
                std::declval<TestSpacesPtr>(), std::declval<SolutionSpacesPtr>(),
                make_tuple(
                  make_IntegralTerm<0,0,IntegrationType::valueValue,
                                        DomainOfIntegration::interior>(0.),
                  make_IntegralTerm<0,0,IntegrationType::gradValue,
                                        DomainOfIntegration::interior>(-1.,
                                           FieldVector<double, dim>{1.,1.}),
                  make_IntegralTerm<0,1,IntegrationType::normalVector,
                                        DomainOfIntegration::face>(1.,
                                           FieldVector<double, dim>{1.,1.}))))
              BilinearForm;
      typedef decltype(make_InnerProduct(
                std::declval<TestSpacesPtr>(),
                make_tuple(
                  make_IntegralTerm<0,0,IntegrationType::gradGrad,
                                        DomainOfIntegration::interior>(1.,
                                           FieldVector<double, dim>{1.,1.}),
                  make_IntegralTerm<0,0,IntegrationType::travelDistanceWeighted,
                                        DomainOfIntegration::face>(1.,
                                           FieldVector<double, dim>{1.,1.}))))
              InnerProduct;

      using BilinearFormEnriched
          = replaceTestSpaces_t<BilinearForm, TestSpacesEnrichedPtr>;
      using InnerProductEnriched
          = replaceTestSpaces_t<InnerProduct, TestSpacesEnrichedPtr>;

      typedef GeometryBuffer<typename LeafGridView::template Codim<0>::Geometry>
          GeometryBuffer_t;

      typedef decltype(make_DPGSystemAssembler(
                  std::declval<BilinearForm&>(),
                  std::declval<InnerProduct&>(),
                  std::declval<GeometryBuffer_t&>()))
              SystemAssembler_t;

      std::vector<SystemAssembler_t> systemAssemblers;
      systemAssemblers.reserve(numS);
      std::vector<GeometryBuffer_t> geometryBuffers(numS);

      /* All the following objects have to be created outside of the
       * following for loop, as the DPGSystemAssembler holds references
       * to them which will otherwise go out of scope. */
      std::vector<BilinearForm> bilinearForms;
      bilinearForms.reserve(numS);
      std::vector<InnerProduct> innerProducts;
      innerProducts.reserve(numS);
      std::vector<BilinearFormEnriched> bilinearFormsEnriched;
      bilinearFormsEnriched.reserve(numS);
      std::vector<InnerProductEnriched> innerProductsEnriched;
      innerProductsEnriched.reserve(numS);

      for(unsigned int i = 0; i < numS; ++i)
      {
        Direction s = sVector[i];

        auto gridView = grids[i]->leafGridView();
        auto solutionSpaces
          = make_space_tuple<FEBasisInterior, FEBasisTrace>(gridView);
        // v enriched
        auto testSpaces = make_space_tuple<FEBasisTest>(gridView);
        auto testSpacesEnriched
          = make_space_tuple<FEBasisTestEnriched>(gridView);

        bilinearForms.emplace_back(
          make_BilinearForm(testSpaces, solutionSpaces,
              make_tuple(
                  make_IntegralTerm<0,0,IntegrationType::valueValue,
                                        DomainOfIntegration::interior>(sigma),
                  make_IntegralTerm<0,0,IntegrationType::gradValue,
                                        DomainOfIntegration::interior>(-1., s),
                  make_IntegralTerm<0,1,IntegrationType::normalVector,
                                        DomainOfIntegration::face>(1., s))));
        bilinearFormsEnriched.emplace_back(
            replaceTestSpaces(bilinearForms[i], testSpacesEnriched));
        innerProducts.emplace_back(
          make_InnerProduct(testSpaces,
              make_tuple(
                  make_IntegralTerm<0,0,IntegrationType::gradGrad,
                                        DomainOfIntegration::interior>(1., s),
                  make_IntegralTerm<0,0,IntegrationType::travelDistanceWeighted,
                                        DomainOfIntegration::face>(1., s))));
        innerProductsEnriched.emplace_back(
            replaceTestSpaces(innerProducts[i], testSpacesEnriched));

        systemAssemblers.emplace_back(
            make_DPGSystemAssembler(bilinearForms[i],
                                    innerProducts[i],
                                    geometryBuffers[i]));
      }

      // Determine Dirichlet dofs for u^ (inflow boundary)
      std::vector<std::vector<bool>> dirichletNodesInflow(numS);
      // Contribution of inflow boundary for the rhs
      std::vector<std::vector<double>> rhsInflowContrib(numS);
      for(unsigned int i = 0; i < numS; ++i)
      {
        Direction s = sVector[i];
        auto solutionSpaces = bilinearForms[i].getSolutionSpaces();
        BoundaryTools::getInflowBoundaryMask(std::get<1>(*solutionSpaces),
                                              dirichletNodesInflow[i],
                                              s);

        auto gSfixed = std::make_tuple([s] (const Domain& x){ return 0.;});
        BoundaryTools::getInflowBoundaryValue(std::get<1>(*solutionSpaces),
                                              rhsInflowContrib[i],
                                              gSfixed);
      }

      std::vector<VectorType> rhs(numS);
      std::vector<MatrixType> stiffnessMatrix(numS);

      /////////////////////////////////////////////////////////
      //  Assemble the systems
      /////////////////////////////////////////////////////////

      // loop of the discrete ordinates
      for(unsigned int i = 0; i < numS; ++i)
      {
        Direction s = sVector[i];
        auto solutionSpaces = systemAssemblers[i].getSolutionSpaces();

        auto rhsFunction = make_LinearForm(
              systemAssemblers[i].getTestSearchSpaces(),
              std::make_tuple(
                make_LinearFunctionalTerm<0, DomainOfIntegration::interior>
                  (rhsFunctional[i], std::get<0>(*solutionSpaces))));
        systemAssemblers[i].assembleSystem(
            stiffnessMatrix[i], rhs[i],
            rhsFunction);
        systemAssemblers[i].template applyDirichletBoundary<1>
            (stiffnessMatrix[i],
             rhs[i],
             dirichletNodesInflow[i],
             rhsInflowContrib[i]);
#if 1
        systemAssemblers[i].template defineCharacteristicFaces<1,dim>(
            stiffnessMatrix[i],
            rhs[i], s);
#endif
      }

      // std::ofstream of("stiffnessNew.dat");
      // printmatrix(of, stiffnessMatrix[0], "stiffnessNew", "--");

      ////////////////////////////////////
      //   Initialize solution vector
      ////////////////////////////////////
      for(unsigned int i = 0; i < numS; ++i)
      {
        auto solutionSpaces = systemAssemblers[i].getSolutionSpaces();
        x[i].resize(std::get<0>(*solutionSpaces).size()
                    + std::get<1>(*solutionSpaces).size());
        x[i] = 0;
      }

      ////////////////////////////
      //   Compute solution
      ////////////////////////////

      std::cout << "rhs size = " << rhs[0].size()
                << " matrix size = " << stiffnessMatrix[0].N()
                            << " x " << stiffnessMatrix[0].M()
                << " solution size = " << x[0].size() << std::endl;


      for(unsigned int i = 0; i < numS; ++i)
      {
        const int verbosity = 0; // 0: not verbose; >0: verbose
        UMFPack<MatrixType> umfPack(stiffnessMatrix[i], verbosity);
        InverseOperatorResult statistics;
        umfPack.apply(x[i], rhs[i], statistics);
      }

      ////////////////////////////////////
      //  A posteriori error
      ////////////////////////////////////
      std::cout << "Compute a posteriori errors\n";
      aposterioriErr = 0.;
      for(unsigned int i = 0; i < numS; ++i)
      {
        const Direction s = sVector[i];

        std::cout << "Direction " << i << '\n';

        // We compute the a posteriori error
        // - We compute the rhs with the enriched test space ("rhs[i]=f(v_i)")
        // -- Contribution of the source term f that has an analytic expression
        // -- Contribution of the scattering term
        auto testSpacesEnriched
          = make_space_tuple<FEBasisTestEnriched>(grids[i]->leafGridView());
        auto rhsAssemblerEnriched = make_RhsAssembler(testSpacesEnriched);
        auto rhsFunction = make_LinearForm(
              rhsAssemblerEnriched.getTestSpaces(),
              std::make_tuple(
                make_LinearFunctionalTerm<0, DomainOfIntegration::interior>
                  (rhsFunctional[i],
                   std::get<0>(*systemAssemblers[i].getSolutionSpaces()))));
        rhsAssemblerEnriched.assembleRhs(rhs[i],
            rhsFunction);
        // - Computation of the a posteriori error
        double aposterioriErr_i = ErrorTools::aPosterioriError(
            bilinearFormsEnriched[i], innerProductsEnriched[i], x[i], rhs[i]);
        aposterioriErr += aposterioriErr_i * aposterioriErr_i;
      }
      aposterioriErr = std::sqrt(aposterioriErr / numS);

      {
        static_assert(!is_RefinedFiniteElement<FEBasisInterior>::value,
            "Functions::interpolate won't work for refined finite elements");
        for(unsigned int i = 0; i < numS; ++i)
        {
          auto feBasisInterior = std::get<0>(
              *systemAssemblers[i].getSolutionSpaces());
          VectorType gInterpolation(feBasisInterior.size());
          const Direction s = sVector[i];
          Functions::interpolate(feBasisInterior, gInterpolation,
              [&g,s](const Direction& x) { return g(x,s); });

          // Add gInterpolation to first feBasisInterior.size() entries of x.
          using Iterator = std::decay_t<decltype(x[i].begin())>;
          for(Iterator xIt=x[i].begin(),
                       xEnd=x[i].begin()+feBasisInterior.size(),
                       gIt=gInterpolation.begin(); xIt!=xEnd; ++xIt, ++gIt) {
            *xIt += *gIt;
          }
          // TODO: Add (interpolation of) g to theta part of x?
        }
      }

      ofs << "Iteration " << n << '.' << nRefinement << ": "
          << "A posteriori estimation of || (u,trace u) - (u_fem,theta) || = "
          << aposterioriErr << ", grid level: " << grids[0]->maxLevel()
          << ", number of DOFs: " << x[0].size() * x.size()
          << ", applying the kernel took "
          << std::chrono::duration_cast<std::chrono::microseconds>
             (endScatteringApproximation - startScatteringApproximation).count()
          << "us, " << kernelApproximation.info()
          << std::endl;
      std::cout << std::endl;

      std::cout << "\nStatistics at end of inner iteration:\n";
      std::cout << "Grid level: " << grids[0]->maxLevel() << '\n';
      std::cout << "A posteriori error: " << aposterioriErr << std::endl;

      if(++nRefinement >= maxNumberOfInnerIterations
          || aposterioriErr <= kappa3*eta) {
        break;
      } else {
        std::vector<VectorType> rhsFunctionalCoarse(numS);
        std::swap(rhsFunctional, rhsFunctionalCoarse);

        for(unsigned int i = 0; i < numS; ++i)
        {
          grids[i]->globalRefine(1);

          const LevelGridView levelGridView
              = grids[i]->levelGridView(grids[i]->maxLevel()-1);
          FEBasisCoarseInterior coarseInteriorBasis(levelGridView);
          FEBasisInterior       feBasisInterior(grids[i]->leafGridView());

          rhsFunctional[i] = interpolateToUniformlyRefinedGrid(
              coarseInteriorBasis, feBasisInterior,
              rhsFunctionalCoarse[i]);
        }
      }
    }

    accuracy = std::pow(rho, n) * CT * fnorm + 2*eta;
    eta /= rhobar;

    if(plotSolutions == PlotSolutions::plotOuterIterations) {
      ////////////////////////////////////////////////////////////////////////
      //  Write result to VTK file
      //  We need to subsample, because VTK cannot natively display
      //  real second-order functions
      ////////////////////////////////////////////////////////////////////////
      std::cout << "Print solutions:\n";

      for(unsigned int i = 0; i < numS; ++i)
      {
        Direction s = sVector[i];

        auto gridView = grids[i]->leafGridView();
        FEBasisInterior feBasisInterior(gridView);
        FEBasisTrace feBasisTrace(gridView);

        std::cout << "Direction " << i << '\n';

        std::string name = std::string("u_rad_trans_n")
                        + std::to_string(n)
                        + std::string("_s")
                        + std::to_string(i);
        FunctionPlotter uPlotter(name);
        uPlotter.plot("u", x[i], feBasisInterior, 0, 0);
        name = std::string("theta_rad_trans_n")
                        + std::to_string(n)
                        + std::string("_s")
                        + std::to_string(i);
        FunctionPlotter thetaPlotter(name);
        thetaPlotter.plot("theta", x[i], feBasisTrace, 2,
                          feBasisInterior.size());
      }
    }
  }
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class FEBasisInterior, class FEBasisTrace, class GridView>
std::vector<Dune::BlockVector<Dune::FieldVector<double, 1> >>
Periter<ScatteringKernelApproximation, RHSApproximation>::apply_scattering(
      ScatteringKernelApproximation& kernelApproximation,
      const std::vector<VectorType>& x,
      const std::vector<GridView>& gridViews,
      double accuracy) {
  kernelApproximation.setAccuracy(accuracy);

  const size_t numS = x.size();
  std::vector<VectorType> rhsFunctional(numS);

  for(unsigned int i = 0; i < numS; ++i) {
    auto solutionSpaces =
        make_space_tuple<FEBasisInterior, FEBasisTrace>(gridViews[i]);
    auto scatteringAssembler_i =
        make_ApproximateScatteringAssembler(solutionSpaces,
                                            kernelApproximation,
                                            i);
    scatteringAssembler_i.template precomputeScattering<0>(
        rhsFunctional[i],
        x);
  }

  return rhsFunctional;
}

} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_PERITER_HH
