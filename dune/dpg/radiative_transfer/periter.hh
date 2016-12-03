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

template<class ScatteringKernelApproximation>
class Periter {
  public:
  template<class Grid, class F, class G, class GDeriv,
           class RHSApproximation, class Kernel>
  void solve(Grid& grid,
             const F& f,
             const G& g,
             const GDeriv& gDeriv,
             RHSApproximation rhsApproximation,
             double sigma,
             const Kernel& kernel,
             unsigned int numS,
             double rho,
             double CT,
             double targetAccuracy,
             unsigned int maxNumberOfIterations,
             PlotSolutions plotSolutions = PlotSolutions::doNotPlot);
};

struct FeRHSandBoundary {};
struct ApproximateRHSandBoundary {};

template<class ScatteringKernelApproximation>
template<class Grid, class F, class G, class GDeriv,
         class RHSApproximation, class Kernel>
void Periter<ScatteringKernelApproximation>::solve(Grid& grid,
           const F& f,
           const G& g,
           const GDeriv& gDeriv,
           RHSApproximation rhsApproximation,
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
  const unsigned int dim = 2;
  constexpr bool rhsIsFeFunction
      = std::is_same<RHSApproximation, FeRHSandBoundary>::value;

  typedef typename Grid::LeafGridView  LeafGridView;
  typedef typename Grid::LevelGridView LevelGridView;
  typedef typename LeafGridView::template Codim<0>::Geometry Geometry;
  LeafGridView gridView = grid.leafGridView();

  ///////////////////////////////////
  // To print information
  ///////////////////////////////////
  std::ofstream ofs("output_rad_trans");

  ////////////////////////////////////////////
  // Handle directions of discrete ordinates
  ////////////////////////////////////////////
  using Domain
    = typename Geometry::GlobalCoordinate;
  using Direction = FieldVector<double, dim>;
  // Vector of directions: sVector
  std::vector< Direction > sVector(numS);
  for(unsigned int i = 0; i < numS; ++i)
  {
    using namespace boost::math::constants;
    sVector[i] = {cos(2*pi<double>()*i/numS),
                  sin(2*pi<double>()*i/numS)};
  }

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
    FEBasisInterior feBasisInterior(gridView);
    FEBasisTrace feBasisTrace(gridView);

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
    kernelApproximation.setAccuracy(kappa1*eta);

    ofs << "\nIteration " << n << std::endl;
    std::cout << "\nIteration " << n << std::endl << std::endl;


    std::chrono::steady_clock::time_point startScatteringApproximation
        = std::chrono::steady_clock::now();
    std::vector<VectorType> rhsFunctional(numS);
    {
      FEBasisInterior feBasisInterior(gridView);
      FEBasisTrace feBasisTrace(gridView);
      auto solutionSpaces =
          std::make_tuple(FEBasisInterior(gridView), FEBasisTrace(gridView));

      for(unsigned int i = 0; i < numS; ++i) {
        auto scatteringAssembler_i =
            make_ApproximateScatteringAssembler(solutionSpaces,
                                                kernelApproximation,
                                                i);
        scatteringAssembler_i.template precomputeScattering<0>(
            rhsFunctional[i],
            x);
      }
    }
    std::chrono::steady_clock::time_point endScatteringApproximation
        = std::chrono::steady_clock::now();

    // Refine grid until accuracy kappa2*eta is reached for
    // approximation of rhs and boundary values.
    if(!rhsIsFeFunction) {
      std::vector<VectorType> boundaryValues(numS);
      const unsigned int maxNumberOfRefinements = 3;
      for(unsigned int refinement = 1; ; refinement++) {
        static_assert(!is_RefinedFiniteElement<FEBasisInterior>::value,
            "Functions::interpolate won't work for refined finite elements");
        FEBasisInterior feBasisInterior(gridView);
        for(unsigned int i = 0; i < numS; ++i)
        {
          VectorType gInterpolation(feBasisInterior.size());
          const Direction s = sVector[i];
          Functions::interpolate(feBasisInterior, gInterpolation,
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
                gridView);
          auto gApprox = Functions::makeDiscreteGlobalBasisFunction<double>(
                feBasisInterior, boundaryValues[i]);

          auto localGExact = localFunction(gExact);
          auto localGApprox = localFunction(gApprox);
          auto localView = feBasisInterior.localView();
          auto localIndexSet = feBasisInterior.localIndexSet();

          double rhsError_i = 0.;
          for(const auto& e : elements(gridView)) {
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
        if(rhsError <= kappa2*eta || refinement == maxNumberOfRefinements) {
          break;
        } else {
          const auto levelGridView = grid.levelGridView(grid.maxLevel());
          FEBasisCoarseInterior coarseInteriorBasis(levelGridView);

          grid.globalRefine(1);
          FEBasisInterior feBasisInterior(gridView);

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
    } else { // rhsIsFeFunction:
      static_assert(!is_RefinedFiniteElement<FEBasisInterior>::value,
          "Functions::interpolate won't work for refined finite elements");
      FEBasisInterior feBasisInterior(gridView);
      for(unsigned int i = 0; i < numS; ++i)
      {
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

      FEBasisInterior feBasisInterior(gridView);
      FEBasisTrace feBasisTrace(gridView);

      auto solutionSpaces
        = std::make_tuple(FEBasisInterior(gridView), FEBasisTrace(gridView));

      // v enriched
      typedef Functions::LagrangeDGBasis<LeafGridView, 4> FEBasisTest;
      auto testSpaces = std::make_tuple(FEBasisTest(gridView));

      typedef FEBasisTest FEBasisTestEnriched;
      auto testSpacesEnriched = std::make_tuple(FEBasisTestEnriched(gridView));
      using TestSpacesEnriched = decltype(testSpacesEnriched);
      auto rhsAssemblerEnriched = make_RhsAssembler(testSpacesEnriched);

      typedef decltype(make_BilinearForm(testSpaces, solutionSpaces,
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
      typedef decltype(make_InnerProduct(testSpaces,
                make_tuple(
                  make_IntegralTerm<0,0,IntegrationType::gradGrad,
                                        DomainOfIntegration::interior>(1.,
                                           FieldVector<double, dim>{1.,1.}),
                  make_IntegralTerm<0,0,IntegrationType::travelDistanceWeighted,
                                        DomainOfIntegration::face>(1.,
                                           FieldVector<double, dim>{1.,1.}))))
              InnerProduct;

      using BilinearFormEnriched
          = replaceTestSpaces_t<BilinearForm, TestSpacesEnriched>;
      using InnerProductEnriched
          = replaceTestSpaces_t<InnerProduct, TestSpacesEnriched>;

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
        BoundaryTools boundaryTools = BoundaryTools();
        boundaryTools.getInflowBoundaryMask(std::get<1>(solutionSpaces),
                                              dirichletNodesInflow[i],
                                              s);

        auto gSfixed = std::make_tuple([s] (const Domain& x){ return 0.;});
        boundaryTools.getInflowBoundaryValue(std::get<1>(solutionSpaces),
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

        auto rhsFunction = make_LinearForm(
              systemAssemblers[i].getTestSearchSpaces(),
              std::make_tuple(
                make_LinearFunctionalTerm<0, DomainOfIntegration::interior>
                  (rhsFunctional[i], std::get<0>(solutionSpaces))));
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
        x[i].resize(feBasisTrace.size()
                   +feBasisInterior.size());
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

        ErrorTools errorTools = ErrorTools();
        // We compute the a posteriori error
        // - We compute the rhs with the enriched test space ("rhs[i]=f(v_i)")
        // -- Contribution of the source term f that has an analytic expression
        // -- Contribution of the scattering term
        auto rhsFunction = make_LinearForm(
              rhsAssemblerEnriched.getTestSpaces(),
              std::make_tuple(
                make_LinearFunctionalTerm<0, DomainOfIntegration::interior>
                  (rhsFunctional[i], std::get<0>(solutionSpaces))));
        rhsAssemblerEnriched.assembleRhs(rhs[i],
            rhsFunction);
        // - Computation of the a posteriori error
        double aposterioriErr_i = errorTools.aPosterioriError(
            bilinearFormsEnriched[i], innerProductsEnriched[i], x[i], rhs[i]);
        aposterioriErr += aposterioriErr_i * aposterioriErr_i;
      }
      aposterioriErr = std::sqrt(aposterioriErr / numS);

      {
        static_assert(!is_RefinedFiniteElement<FEBasisInterior>::value,
            "Functions::interpolate won't work for refined finite elements");
        for(unsigned int i = 0; i < numS; ++i)
        {
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
          << aposterioriErr << ", grid level: " << grid.maxLevel()
          << ", number of DOFs: " << x[0].size() * x.size()
          << ", applying the kernel took "
          << std::chrono::duration_cast<std::chrono::microseconds>
             (endScatteringApproximation - startScatteringApproximation).count()
          << "us, " << kernelApproximation.info()
          << std::endl;
      std::cout << std::endl;

      std::cout << "\nStatistics at end of inner iteration:\n";
      std::cout << "Grid level: " << grid.maxLevel() << '\n';
      std::cout << "A posteriori error: " << aposterioriErr << std::endl;

      if(++nRefinement >= maxNumberOfInnerIterations
          || aposterioriErr <= kappa3*eta) {
        break;
      } else {
        grid.globalRefine(1);

        const LevelGridView levelGridView
            = grid.levelGridView(grid.maxLevel()-1);
        FEBasisCoarseInterior coarseInteriorBasis(levelGridView);
        FEBasisInterior       feBasisInterior(gridView);

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

    accuracy = std::pow(rho, n) * CT * fnorm + 2*eta;
    eta /= rhobar;

    if(plotSolutions == PlotSolutions::plotOuterIterations) {
      ////////////////////////////////////////////////////////////////////////
      //  Write result to VTK file
      //  We need to subsample, because VTK cannot natively display
      //  real second-order functions
      ////////////////////////////////////////////////////////////////////////
      std::cout << "Print solutions:\n";

      FEBasisInterior feBasisInterior(gridView);
      FEBasisTrace feBasisTrace(gridView);


      for(unsigned int i = 0; i < numS; ++i)
      {
        Direction s = sVector[i];

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

} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_PERITER_HH
