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

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pqktracenodalbasis.hh>
#include <dune/functions/functionspacebases/optimaltestbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/pqksubsampleddgbasis.hh>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/dpg/system_assembler.hh>
#include <dune/dpg/errortools.hh>
#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/rhs_assembler.hh>
#include <dune/dpg/radiative_transfer/approximate_scattering.hh>

#include <boost/math/constants/constants.hpp>

namespace Dune {

// template<class Apply, class Solve>
class Periter {
  public:
  template<class GridView, class F, class Kernel>
  void solve(GridView gridView,
             const F& f,
             const Kernel& kernel,
             unsigned int numS,
             /* TODO: replace with accuracy */
             unsigned int numberOfIterations);
};

// Get solution u or theta out of the solution vector x
template<class FieldVector>
void extractSolution(std::vector< FieldVector >& u,
                     const std::vector< FieldVector >& x,
                     const int offset
                     )
{
  int numS = x.size();
  int i_max = u[0].size();
  for(int iDir=0;iDir<numS;iDir++){
    for(int i=0;i<i_max;i++){
      u[iDir][i] = x[iDir][i+offset];
    }
  }
}

template<class GridView, class F, class Kernel>
void Periter::solve(GridView gridView,
           const F& f,
           const Kernel& kernel,
           unsigned int numS,
           /* TODO: replace with accuracy */
           unsigned int numberOfIterations) {
  const unsigned int dim = 2;

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
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  typedef Functions::LagrangeDGBasis<GridView, 1> FEBasisInterior; // u
  FEBasisInterior feBasisInterior(gridView);

  typedef Functions::PQkTraceNodalBasis<GridView, 2> FEBasisTrace; // u^
  FEBasisTrace feBasisTrace(gridView);

  auto solutionSpaces = std::make_tuple(FEBasisInterior(gridView), FEBasisTrace(gridView));

  typedef Functions::PQkSubsampledDGNodalBasis<GridView, 4, 3> FEBasisTest; // v enriched
  auto testSpaces = std::make_tuple(FEBasisTest(gridView));

  auto rhsAssembler = make_RhsAssembler(testSpaces);

  typedef decltype(testSpaces) TestSpaces;
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
              make_IntegralTerm<0,0,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(1.),
              make_IntegralTerm<0,0,IntegrationType::gradGrad,
                                    DomainOfIntegration::interior>(1.,
                                       FieldVector<double, dim>{1.,1.}))))
          InnerProduct;

  typedef Functions::TestspaceCoefficientMatrix<BilinearForm, InnerProduct>
      TestspaceCoefficientMatrix;
  typedef Functions::OptimalTestBasis<TestspaceCoefficientMatrix>
      FEBasisOptimalTest;              // v

  typedef decltype(make_DPG_SystemAssembler(
#if 1
              std::declval<std::tuple<FEBasisOptimalTest>>(), solutionSpaces,
#else
              testSpaces, solutionSpaces,
#endif
              std::declval<BilinearForm>()))
          SystemAssembler_t;

  std::vector<SystemAssembler_t> systemAssemblers;
  systemAssemblers.reserve(numS);

  ScatteringKernelApproximation::SVD kernelSVD(kernel, numS);

  // Scattering assemblers with optimal test spaces
  std::vector<ApproximateScatteringAssembler
                  <std::tuple<FEBasisOptimalTest>,
                   SolutionSpaces,
                   decltype(kernelSVD),
                   DPGFormulation>
             > scatteringAssemblers;
  scatteringAssemblers.reserve(numS);

  // Scattering assembler with enriched test space
  std::vector<ApproximateScatteringAssembler
                  <std::tuple<FEBasisTest>,
                   SolutionSpaces,
                   decltype(kernelSVD),
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
              make_IntegralTerm<0,0,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(1.),
              make_IntegralTerm<0,0,IntegrationType::gradGrad,
                                    DomainOfIntegration::interior>(1., s))));

    coefficientMatrices.emplace_back(bilinearForms[i], innerProducts[i]);

    optimalTestSpaces.emplace_back(
            make_tuple(FEBasisOptimalTest(coefficientMatrices[i])));

    systemAssemblers.emplace_back(
        make_DPG_SystemAssembler(optimalTestSpaces[i], solutionSpaces,
                                 bilinearForms[i]));
    scatteringAssemblers.emplace_back(
        make_DPG_ApproximateScatteringAssembler(optimalTestSpaces[i],
                                                solutionSpaces,
                                                kernelSVD,
                                                i));
    scatteringAssemblersEnriched.emplace_back(
        make_DPG_ApproximateScatteringAssembler(
            testSpaces,
            solutionSpaces,
            kernelSVD,
            i));
  }

  /////////////////////////////////////////////////////////
  //   Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////
  typedef BlockVector<FieldVector<double,1> > VectorType;
  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

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
  //   Choose an initial iterate
  /////////////////////////////////////////////////
  std::vector<VectorType> x,xPrevious;
  x.reserve(numS);
  xPrevious.reserve(numS);
  for(unsigned int i = 0; i < numS; ++i)
  {
    x.emplace_back(feBasisTrace.indexSet().size()
                   +feBasisInterior.indexSet().size());
    xPrevious.emplace_back(feBasisTrace.indexSet().size()
                   +feBasisInterior.indexSet().size());
    x[i] = 0;
    xPrevious[i] = 0;
  }

  ///////////////////////////////////////////////////
  // Vector to store solution of previous iteration
  // (useful to compute error between two iterates)
  ///////////////////////////////////////////////////
  std::vector<VectorType> u,uPrevious;
  u.reserve(numS);
  uPrevious.reserve(numS);
  for(unsigned int i = 0; i < numS; ++i)
  {
    u.emplace_back(feBasisInterior.indexSet().size());
    u[i] = 0;
    uPrevious.emplace_back(feBasisInterior.indexSet().size());
    uPrevious[i] = 0;
  }

  std::vector<VectorType> theta,thetaPrevious;
  theta.reserve(numS);
  thetaPrevious.reserve(numS);
  for(unsigned int i = 0; i < numS; ++i)
  {
    theta.emplace_back(feBasisTrace.indexSet().size());
    theta[i] = 0;
    thetaPrevious.emplace_back(feBasisTrace.indexSet().size());
    thetaPrevious[i] = 0;
  }

  VectorType diffU(feBasisInterior.indexSet().size());
  diffU = 0;

  VectorType diffTheta(feBasisTrace.indexSet().size());
  diffTheta = 0;

  /////////////////////////////////////////////////////////
  //  Fixed-point iterations
  /////////////////////////////////////////////////////////
  // TODO: Estimate ρ from the paper.
  const double rho = .5;
  // The accuracy η_n:
  double accuracy = 1.;
  for(unsigned int n = 0; n < numberOfIterations; ++n)
  {
    accuracy *= rho/2.;
    ofs << "Iteration " << n << std::endl;
    std::cout << "Iteration " << n << std::endl << std::endl;

    /////////////////////////////////////////////////////////
    //  Update solutions
    /////////////////////////////////////////////////////////
    xPrevious = x;
    extractSolution(uPrevious,xPrevious,0);
    extractSolution(thetaPrevious,xPrevious,feBasisInterior.indexSet().size());

    /////////////////////////////////////////////////////////
    //  Assemble the systems
    /////////////////////////////////////////////////////////
    // using Domain = GridType::template Codim<0>::Geometry::GlobalCoordinate;
    //auto f = [] (const Domain& x, const Direction& s) { return 1.;};

    // TODO: Consider the norm of the transport solver in the accuracy.
    kernelSVD.setAccuracy(accuracy/2.);

    // loop of the discrete ordinates
    for(unsigned int i = 0; i < numS; ++i)
    {
      Direction s = sVector[i];

      systemAssemblers[i].assembleSystem(
          stiffnessMatrix[i], rhs[i],
          std::make_tuple([s,&f] (const Domain& x) { return f(x,s); }));
      VectorType scattering;
      scatteringAssemblers[i].template assembleScattering<0>(
          scattering,
          xPrevious);
      rhs[i] += scattering;
      systemAssemblers[i].template applyDirichletBoundarySolution<1>
          (stiffnessMatrix[i],
           rhs[i],
           dirichletNodesInflow[i],
           rhsInflowContrib[i]);
      systemAssemblers[i].template defineCharacteristicFaces<1,dim>(
          stiffnessMatrix[i],
          rhs[i], s);
    }

    // std::ofstream of("stiffnessNew.dat");
    // printmatrix(of, stiffnessMatrix[0], "stiffnessNew", "--");

    ////////////////////////////
    //   Compute solution
    ////////////////////////////

    std::cout <<"rhs size = "<< rhs[0].size()
              <<" matrix size = " << stiffnessMatrix[0].N()
                         << " x " << stiffnessMatrix[0].M()
              <<" solution size = "<< x[0].size() <<std::endl;


    for(unsigned int i = 0; i < numS; ++i)
    {
      int verbosity = 0; // 0: not verbose; >0: verbose
      UMFPack<MatrixType> umfPack(stiffnessMatrix[i], verbosity);
      InverseOperatorResult statistics;
      umfPack.apply(x[i], rhs[i], statistics);
    }

    extractSolution(u,x,0);
    extractSolution(theta,x,feBasisInterior.indexSet().size());

    ////////////////////////////////////
    //  Error computation and print in VTK file
    ////////////////////////////////////
    for(unsigned int i = 0; i < numS; ++i)
    {
      Direction s = sVector[i];

      std::cout << "Direction " << i << std::endl;

      ////////////////////////////////////////////////////////////////////////
      //  Write result to VTK file
      //  We need to subsample, because VTK cannot natively display
      //  real second-order functions
      ////////////////////////////////////////////////////////////////////////
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

      ////////////////////////////////////
      //  A posteriori error
      ////////////////////////////////////
      ErrorTools errorTools = ErrorTools();
      // We compute the a posteriori error
      // - We compute the rhs with the enriched test space ("rhs[i]=f(v_i)")
      // -- Contribution of the source term f that has an analytic expression
      rhsAssembler.assembleRhs(rhs[i],
          std::make_tuple([s,&f] (const Domain& x) { return f(x,s); }));
      // -- Contribution of the scattering term
      VectorType scattering;
      kernelSVD.setAccuracy(0.);
      scatteringAssemblersEnriched[i]
          .template assembleScattering<0>(scattering, xPrevious);
      rhs[i] += scattering;
      // - Computation of the a posteriori error
      double aposterioriErr = errorTools.aPosterioriError(bilinearForms[i],innerProducts[i],u[i],theta[i],rhs[i]); //change with contribution of scattering rhs[i]
      ofs << "A posteriori estimation of || (u,trace u) - (u_fem,theta) || = " << aposterioriErr << std::endl;

      // We compute the L2 error wrt previous iterate
      for (unsigned int j=0; j<feBasisInterior.indexSet().size(); j++)
        diffU[j] = u[i][j]-uPrevious[i][j];

      for (unsigned int j=0; j<feBasisTrace.indexSet().size(); j++)
        diffTheta[j] = theta[i][j]-thetaPrevious[i][j];

      ofs << "Diff wrt previous iteration: " << std::endl;
      ofs << "  -> || u["<< i << "] - u_previous["<< i <<"] ||_L2 = " << diffU.two_norm() << std::endl;
      ofs << "  -> || theta["<< i << "] - theta_previous["<< i <<"] ||_L2 = " << diffTheta.two_norm() << std::endl << std::endl;

    }
    ofs << std::endl;
    std::cout << std::endl;
  }

  ofs.std::ofstream::close();
}

} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_PERITER_HH
