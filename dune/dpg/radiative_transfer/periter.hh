// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_PERITER_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_PERITER_HH

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

#include <dune/functions/functionspacebases/optimaltestbasis.hh>
#include <dune/functions/functionspacebases/pqkdgrefineddgnodalbasis.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pqktracenodalbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/pqksubsampleddgbasis.hh>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/errortools.hh>
#include <dune/dpg/functions/interpolate.hh>
#include <dune/dpg/radiative_transfer/approximate_scattering.hh>
#include <dune/dpg/rhs_assembler.hh>
#include <dune/dpg/system_assembler.hh>
#include <dune/dpg/type_traits.hh>

#include <boost/math/constants/constants.hpp>

namespace Dune {

template<class ScatteringKernelApproximation>
class Periter {
  public:
  template<class Grid, class F, class Kernel>
  void solve(Grid& grid,
             const F& f,
             const Kernel& kernel,
             unsigned int numS,
             double targetAccuracy,
             unsigned int numberOfIterations);
};

// Get solution u or theta out of the solution vector x
template<class VectorType>
std::vector<VectorType>
extractSolution(const std::vector<VectorType>& x,
                unsigned int offset,
                unsigned int size
               )
{
  unsigned int numS = x.size();
  std::vector<VectorType> u(numS);
  for(unsigned int iDir=0; iDir<numS; iDir++) {
    u[iDir].resize(size);
    std::copy(x[iDir].begin()+offset,
              x[iDir].begin()+offset+size,
              u[iDir].begin());
  }
  return u;
}

template<class ScatteringKernelApproximation>
template<class Grid, class F, class Kernel>
void Periter<ScatteringKernelApproximation>::solve(Grid& grid,
           const F& f,
           const Kernel& kernel,
           unsigned int numS,
           double targetAccuracy,
           unsigned int maxNumberOfIterations) {
  const unsigned int dim = 2;

  typedef typename Grid::LeafGridView GridView;
  GridView gridView = grid.leafGridView();

  ///////////////////////////////////
  // To print information
  ///////////////////////////////////
  std::ofstream ofs("output_rad_trans");

  ////////////////////////////////////////////
  // Handle directions of discrete ordinates
  ////////////////////////////////////////////
  using Domain
    = typename GridView::template Codim<0>::Geometry::GlobalCoordinate;
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

  typedef Functions::LagrangeDGBasis<GridView, 1> FEBasisInterior; // u
  typedef Functions::PQkTraceNodalBasis<GridView, 2> FEBasisTrace; // u^

  /////////////////////////////////////////////////////////
  //   Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////
  typedef BlockVector<FieldVector<double,1> > VectorType;
  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

  /////////////////////////////////////////////////
  //   Solution vectors
  /////////////////////////////////////////////////
  std::vector<VectorType> x(numS), xPrevious(numS);

  for(unsigned int i = 0; i < numS; ++i)
  {
    FEBasisInterior feBasisInterior(gridView);
    FEBasisTrace feBasisTrace(gridView);

    x[i].resize(feBasisTrace.size()
               +feBasisInterior.size());
    x[i] = 0;
  }

  /////////////////////////////////////////////////////////
  //  Fixed-point iterations
  /////////////////////////////////////////////////////////
  // TODO: Estimate ρ from the paper.
  const double rho = 1.7/5.;
  // The accuracy η_n:
  double accuracy = 1.;
  for(unsigned int n = 0; accuracy > targetAccuracy
                          && n < maxNumberOfIterations; ++n)
  {
    accuracy *= rho/2.;
    ofs << "Iteration " << n << std::endl;
    std::cout << "Iteration " << n << std::endl << std::endl;

    std::vector<VectorType> u, theta;

    /////////////////////////////////////////////////////////
    //  Get previous solutions
    /////////////////////////////////////////////////////////
    std::swap(xPrevious, x);
    size_t numDofsInterior, numDofsTrace;
    {
      FEBasisInterior feBasisInterior(gridView);
      FEBasisTrace feBasisTrace(gridView);

      numDofsInterior = feBasisInterior.size();
      numDofsTrace = feBasisTrace.size();
    }
    std::vector<VectorType> uPrevious
        = extractSolution(xPrevious, 0, numDofsInterior);
    std::vector<VectorType> thetaPrevious
        = extractSolution(xPrevious, numDofsInterior,
                          numDofsTrace);

    const auto levelGridView = grid.levelGridView(grid.maxLevel());
    typedef std::decay_t<decltype(levelGridView)> LevelGridView;

    typedef changeGridView_t<FEBasisInterior, LevelGridView>
            FEBasisCoarseInterior;
    FEBasisCoarseInterior coarseInteriorBasis(levelGridView);

    ////////////////////////////////////////////////////
    // Inner loop
    ////////////////////////////////////////////////////
    const unsigned int maxNumberOfInnerIterations = 10;
    double aposterioriErr;
    for(unsigned int nRefinement = 0;
        // At the end of the loop, we will break if
        // aposterioriErr < accuracy/2.
        nRefinement < maxNumberOfInnerIterations;
        ++nRefinement)
    {
      std::cout << "Inner iteration " << nRefinement << std::endl;

      FEBasisInterior feBasisInterior(gridView);
      FEBasisTrace feBasisTrace(gridView);

      auto solutionSpaces = std::make_tuple(FEBasisInterior(gridView), FEBasisTrace(gridView));

      typedef Functions::LagrangeDGBasis<GridView, 4> FEBasisTest; // v enriched
      auto testSpaces = std::make_tuple(FEBasisTest(gridView));

      auto rhsAssembler = make_RhsAssembler(testSpaces);

      typedef FEBasisTest FEBasisTestEnriched;
      auto testSpacesEnriched = std::make_tuple(FEBasisTestEnriched(gridView));

      // typedef decltype(testSpaces) TestSpaces;
      typedef decltype(solutionSpaces) SolutionSpaces;

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

      typedef Functions::TestspaceCoefficientMatrix<BilinearForm, InnerProduct>
          TestspaceCoefficientMatrix;
      typedef Functions::OptimalTestBasis<TestspaceCoefficientMatrix>
          FEBasisOptimalTest;              // v

      typedef decltype(make_DPG_SystemAssembler(
                  std::declval<std::tuple<FEBasisOptimalTest>>(), solutionSpaces,
                  std::declval<BilinearForm>()))
              SystemAssembler_t;

      std::vector<SystemAssembler_t> systemAssemblers;
      systemAssemblers.reserve(numS);

      ScatteringKernelApproximation kernelApproximation(kernel, numS);

      // Scattering assemblers with optimal test spaces
      std::vector<ApproximateScatteringAssembler
                      <std::tuple<FEBasisOptimalTest>,
                       SolutionSpaces,
                       decltype(kernelApproximation),
                       DPGFormulation>
                 > scatteringAssemblers;
      scatteringAssemblers.reserve(numS);

      // Scattering assembler with enriched test space
      std::vector<ApproximateScatteringAssembler
                      <std::tuple<FEBasisTestEnriched>,
                       SolutionSpaces,
                       decltype(kernelApproximation),
                       DPGFormulation>
                 > scatteringAssemblersEnriched;
      scatteringAssemblersEnriched.reserve(numS);

      /* create an FEBasisOptimalTest for each direction */
      std::vector<std::tuple<FEBasisOptimalTest> > optimalTestSpaces;
      optimalTestSpaces.reserve(numS);
      /* All the following objects have to be created outside of the
       * following for loop, as the optimalTestSpace holds references
       * to them which will otherwise go out of scope. */
      std::vector<BilinearForm> bilinearForms;
      bilinearForms.reserve(numS);
      std::vector<InnerProduct> innerProducts;
      innerProducts.reserve(numS);
      std::vector<TestspaceCoefficientMatrix> coefficientMatrices;
      coefficientMatrices.reserve(numS);

      for(unsigned int i = 0; i < numS; ++i)
      {
        Direction s = sVector[i];

        bilinearForms.emplace_back(
          make_BilinearForm(testSpaces, solutionSpaces,
              make_tuple(
                  make_IntegralTerm<0,0,IntegrationType::valueValue,
                                        DomainOfIntegration::interior>(5.),
                  make_IntegralTerm<0,0,IntegrationType::gradValue,
                                        DomainOfIntegration::interior>(-1., s),
                  make_IntegralTerm<0,1,IntegrationType::normalVector,
                                        DomainOfIntegration::face>(1., s))));
        innerProducts.emplace_back(
          make_InnerProduct(testSpaces,
              make_tuple(
                  make_IntegralTerm<0,0,IntegrationType::gradGrad,
                                        DomainOfIntegration::interior>(1., s),
                  make_IntegralTerm<0,0,IntegrationType::travelDistanceWeighted,
                                        DomainOfIntegration::face>(1., s))));

        coefficientMatrices.emplace_back(bilinearForms[i], innerProducts[i]);

        optimalTestSpaces.emplace_back(
                make_tuple(FEBasisOptimalTest(coefficientMatrices[i])));

        systemAssemblers.emplace_back(
            make_DPG_SystemAssembler(optimalTestSpaces[i], solutionSpaces,
                                     bilinearForms[i]));
        scatteringAssemblers.emplace_back(
            make_DPG_ApproximateScatteringAssembler(optimalTestSpaces[i],
                                                    solutionSpaces,
                                                    kernelApproximation,
                                                    i));
        scatteringAssemblersEnriched.emplace_back(
            make_DPG_ApproximateScatteringAssembler(
                testSpacesEnriched,
                solutionSpaces,
                kernelApproximation,
                i));
      }

      std::vector<VectorType> rhs(numS);
      std::vector<MatrixType> stiffnessMatrix(numS);

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

      /////////////////////////////////////////////////
      //   Initialize solution vector
      /////////////////////////////////////////////////
      for(unsigned int i = 0; i < numS; ++i)
      {
        x[i].resize(feBasisTrace.size()
                   +feBasisInterior.size());
        x[i] = 0;
      }

      std::vector<VectorType> uPreviousFine(numS);
      for(unsigned int i = 0; i < numS; ++i)
      {
        uPreviousFine[i] = interpolateToUniformlyRefinedGrid(
            coarseInteriorBasis, feBasisInterior, uPrevious[i]);
      }

      /////////////////////////////////////////////////////////
      //  Assemble the systems
      /////////////////////////////////////////////////////////
      // using Domain = GridType::template Codim<0>::Geometry::GlobalCoordinate;
      //auto f = [] (const Domain& x, const Direction& s) { return 1.;};

      // TODO: Consider the norm of the transport solver in the accuracy.
      kernelApproximation.setAccuracy(accuracy/2.);

      // loop of the discrete ordinates
      for(unsigned int i = 0; i < numS; ++i)
      {
        Direction s = sVector[i];

        auto rhsFunction = make_DPG_LinearForm(
              systemAssemblers[i].getTestSpaces(),
              std::make_tuple(
                make_LinearIntegralTerm
                  < 0
                  , LinearIntegrationType::valueFunction
                  , DomainOfIntegration::interior>
                  ([s,&f] (const Domain& x) { return f(x,s); })));
        systemAssemblers[i].assembleSystem(
            stiffnessMatrix[i], rhs[i],
            rhsFunction);
        VectorType scattering;
        scatteringAssemblers[i].template assembleScattering<0>(
            scattering,
            uPreviousFine);
        rhs[i] += scattering;
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

      ////////////////////////////
      //   Compute solution
      ////////////////////////////

      std::cout << "rhs size = " << rhs[0].size()
                << " matrix size = " << stiffnessMatrix[0].N()
                            << " x " << stiffnessMatrix[0].M()
                << " solution size = " << x[0].size() << std::endl;


      for(unsigned int i = 0; i < numS; ++i)
      {
        int verbosity = 0; // 0: not verbose; >0: verbose
        UMFPack<MatrixType> umfPack(stiffnessMatrix[i], verbosity);
        InverseOperatorResult statistics;
        umfPack.apply(x[i], rhs[i], statistics);
      }

      u = extractSolution(x, 0, feBasisInterior.size());
      theta = extractSolution(x, feBasisInterior.size(), feBasisTrace.size());

      ////////////////////////////////////
      //  A posteriori error
      ////////////////////////////////////
      kernelApproximation.setAccuracy(0.);
      std::cout << "Compute a posteriori errors\n";
      for(unsigned int i = 0; i < numS; ++i)
      {
        Direction s = sVector[i];

        std::cout << "Direction " << i << '\n';

        ErrorTools errorTools = ErrorTools();
        // We compute the a posteriori error
        // - We compute the rhs with the enriched test space ("rhs[i]=f(v_i)")
        // -- Contribution of the source term f that has an analytic expression
        auto rhsFunction = make_DPG_LinearForm(
              rhsAssembler.getTestSpaces(),
              std::make_tuple(
                make_LinearIntegralTerm
                  < 0
                  , LinearIntegrationType::valueFunction
                  , DomainOfIntegration::interior>
                  ([s,&f] (const Domain& x) { return f(x,s); })));
        rhsAssembler.assembleRhs(rhs[i],
            rhsFunction);
        // -- Contribution of the scattering term
        VectorType scattering;
        scatteringAssemblersEnriched[i]
            .template assembleScattering<0>(scattering, uPreviousFine);
        rhs[i] += scattering;
        // - Computation of the a posteriori error
        aposterioriErr = errorTools.aPosterioriError(
            bilinearForms[i], innerProducts[i], x[i], rhs[i]);
            //change with contribution of scattering rhs[i]
        ofs << "A posteriori estimation of || (u,trace u) - (u_fem,theta) || = " << aposterioriErr << std::endl;
      }
      ofs << std::endl;
      std::cout << std::endl;

      if(aposterioriErr < accuracy/2.) {
        break;
      } else {
        grid.globalRefine(1);
      }
    }

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

      // - Make a discrete function from the FE basis and the coefficient vector
      auto uFunction
          = Dune::Functions::makeDiscreteGlobalBasisFunction<double>
                (feBasisInterior, Dune::TypeTree::hybridTreePath(), u[i]);
      auto localUFunction = localFunction(uFunction);

      auto thetaFunction
          = Dune::Functions::makeDiscreteGlobalBasisFunction<double>
                (feBasisTrace, Dune::TypeTree::hybridTreePath(), theta[i]);
      auto localThetaFunction = localFunction(thetaFunction);
      // - VTK writer
      SubsamplingVTKWriter<GridView> vtkWriterInterior(gridView,0);
      vtkWriterInterior.addVertexData(localUFunction,
                      VTK::FieldInfo("u", VTK::FieldInfo::Type::scalar, 1));
      std::string name = std::string("u_rad_trans_n")
                       + std::to_string(n)
                       + std::string("_s")
                       + std::to_string(i);
      vtkWriterInterior.write(name);

      SubsamplingVTKWriter<GridView> vtkWriterTrace(gridView,2);
      vtkWriterTrace.addVertexData(localThetaFunction, VTK::FieldInfo("theta",VTK::FieldInfo::Type::scalar, 1));
      name = std::string("theta_rad_trans_n")
                       + std::to_string(n)
                       + std::string("_s")
                       + std::to_string(i);
      vtkWriterTrace.write(name);
    }
  }
}

} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_PERITER_HH
