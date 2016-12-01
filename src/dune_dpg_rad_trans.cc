#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <fstream>
#include <cstdlib> // for std::abort()

#include <vector>

#include <dune/common/exceptions.hh> // We use exceptions

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>



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

#include <dune/dpg/dpg_system_assembler.hh>
#include <dune/dpg/errortools.hh>
#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/rhs_assembler.hh>

#include <dune/dpg/radiative_transfer/scattering.hh>

#include <boost/math/constants/constants.hpp>


using namespace Dune;


// Value of the analytic solution "for the interior of the domain"
template <class Domain,class Direction>
double fInner(const Domain& x,
              const Direction& s)
{
  double c0 = 1.;
  double c1 = 1.;
  double value = std::expm1(c0*x[0])*std::expm1(c1*x[1]);//v pure transport
  //double value = 1-(x[0]-0.5)*(x[0]-0.5)-(x[1]-0.5)*(x[1]-0.5); //v RT
  return value;
}
// Partial derivative of fInner with respect to x[0]
template <class Domain,class Direction>
double fInnerD0(const Domain& x,
                const Direction& s)
{
  double c0 = 1.;
  double c1 = 1.;
  double value = c0*std::exp(c0*x[0])*std::expm1(c1*x[1]);//v pure transport
  //double value = -2*(x[0]-0.5); //v RT
  return value;
}
// Partial derivative of fInner with respect to x[1]
template <class Domain,class Direction>
double fInnerD1(const Domain& x,
                const Direction& s)
{
  double c0 = 1.;
  double c1 = 1.;
  double value = std::expm1(c0*x[0])*std::exp(c1*x[1])*c1;//v pure transport
  //double value = -2*(x[1]-0.5); //v RT
  return value;
}

// This function satifies the zero incoming flux bounday conditions
template <class Domain,class Direction>
double fBoundary(const Domain& x,
                 const Direction& s)
{
  // double value = 1.;//v pure transport
  double value = ( (s[0]>0)*x[0] + (s[0]==0)*1. + (s[0]<0)*(1-x[0]) )*
                 ( (s[1]>0)*x[1] + (s[1]==0)*1. + (s[1]<0)*(1-x[1]) );//v RT
  return value;
}
// Partial derivative of fBoundary with respect to x[0]
template <class Domain,class Direction>
double fBoundaryD0(const Domain& x,
                   const Direction& s)
{
  //double value = 0.;//v pure transport
  double value = ( (s[0]>0)*1 + (s[0]==0)*0. + (s[0]<0)*(-1.) )*
                 ( (s[1]>0)*x[1] + (s[1]==0)*1. + (s[1]<0)*(1-x[1]) );//v RT
  return value;
}
// Partial derivative of fBoundary with respect to x[1]
template <class Domain,class Direction>
double fBoundaryD1(const Domain& x,
                   const Direction& s)
{
  // double value = 0.;//v pure transport
  double value = ( (s[0]>0)*x[0] + (s[0]==0)*1. + (s[0]<0)*(1-x[0]) )*
                 ( (s[1]>0)*1 + (s[1]==0)*0. + (s[1]<0)*(-1.) ); //v RT
  return value;
}

//The analytic solution
template <class Domain,class Direction>
double uAnalytic(const Domain& x,
                 const Direction& s)
{
  return fInner(x,s)*fBoundary(x,s);
}

// Optical parameter: sigma
template <class Domain,class Direction>
double sigma(const Domain& x,
             const Direction& s)
{
  return 5.;
}
// Optical parameter: kernel
template <class Domain,class Direction>
double kernel(const Domain& x,
              const Direction& sIntegration,
              const Direction& s)
{
  return 1.7;
}
// The scattering kernel computed with an analytic solution u
// and a mid-point rule for the quadrature formula
template <class Domain,class Direction,class lambdaExpr>
double collisionTerm(const Domain& x,
                     const Direction& s,
                     const lambdaExpr& u,
                     const std::vector< Direction >& sVector)
{
  int numS = sVector.size();
  double sum = 0.;
  for(int i=0; i<numS; i++)
  {
    sum += kernel(x,sVector[i],s)*u(x,sVector[i]);
  }

  return sum/numS;
}

// The right hand-side
template <class Domain,class Direction,class lambdaExpr>
double f(const Domain& x,
         const Direction& s,
         const lambdaExpr& u,
         const std::vector< Direction >& sVector)
{
  double value = s[0]*( fInnerD0(x,s) * fBoundary(x,s) + fInner(x,s)*fBoundaryD0(x,s)) +
                 s[1]*( fInnerD1(x,s) * fBoundary(x,s) + fInner(x,s)*fBoundaryD1(x,s)) +
                 sigma(x,s)*u(x,s) - collisionTerm(x,s,u,sVector);
  return value;
}


// The values on the incoming boundary
template <class Domain,class Direction>
double gFunction(const Domain& x,
                 const Direction& s)
{
  return 1.5;
}

// Get solution u or theta out of the solution vector x
template<class FieldVector>
void extractSolution(std::vector< FieldVector >& u,
                     const std::vector< FieldVector >& x,
                     unsigned int offset
                     )
{
  unsigned int numS = x.size();
  unsigned int i_max = u[0].size();
  for(unsigned int iDir=0;iDir<numS;iDir++){
    for(unsigned int i=0;i<i_max;i++){
      u[iDir][i] = x[iDir][i+offset];
    }
  }
}


int main(int argc, char** argv)
{
  try{

  ///////////////////////////////////
  // To print information
  ///////////////////////////////////
  std::ofstream ofs("output_rad_trans");

  ///////////////////////////////////
  // Get arguments
  // argv[1]: number of discrete ordinates
  // argv[2]: number of fixed-point iterations
  // argv[3]: size of grid
  ///////////////////////////////////

  if(argc != 4) {
      std::cerr << "Usage: " << argv[0] << " <# of ordinates>"
                << " <# of iterations>"
                << " <size of grid>" << std::endl;
      std::abort();
  }

  int numS = atoi(argv[1]);
  int N = atoi(argv[2]);
  unsigned int sizeGrid = atoi(argv[3]);

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  const int dim = 2;
  typedef UGGrid<dim> GridType;

  FieldVector<double,dim> lower = {0,0};
  FieldVector<double,dim> upper = {1,1};
  array<unsigned int,dim> elements = {sizeGrid,sizeGrid};

  shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createCubeGrid(lower, upper, elements);

  //shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createSimplexGrid(lower, upper, elements);

  //shared_ptr<GridType> grid = shared_ptr<GridType>(GmshReader<GridType>::read("irregular-square.msh"));

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid->leafGridView();

  ///////////////////////////////////
  // Handle directions of discrete ordinates
  ///////////////////////////////////
  using Domain = GridType::template Codim<0>::Geometry::GlobalCoordinate;
  using Direction = FieldVector<double, dim>;
  // Vector of directions: sVector
  std::vector< Direction > sVector(numS);
  for(int i = 0; i < numS; ++i)
  {
    using namespace boost::math::constants;
    sVector[i] = {cos(2*pi<double>()*i/numS),
                  sin(2*pi<double>()*i/numS)};
  }

  ///////////////////////////////////////////////
  // Define the analytical solution (if possible)
  ///////////////////////////////////////////////
  //auto uExact = std::make_tuple(uAnalytic);
  auto uExact = [] (const Domain& x, const Direction& s){ return uAnalytic(x,s);};

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
              make_IntegralTerm<0,0,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(1.),
              make_IntegralTerm<0,0,IntegrationType::gradGrad,
                                    DomainOfIntegration::interior>(1.,
                                       FieldVector<double, dim>{1.,1.}))))
          InnerProduct;

  typedef UnbufferedTestspaceCoefficientMatrix<BilinearForm, InnerProduct>
      TestspaceCoefficientMatrix;
  typedef Functions::OptimalTestBasis<TestspaceCoefficientMatrix>
      FEBasisOptimalTest;              // v

  typedef GeometryBuffer<GridView::template Codim<0>::Geometry>
      GeometryBuffer_t;

  typedef decltype(make_DPGSystemAssembler(
              std::declval<BilinearForm&>(),
              std::declval<InnerProduct&>(),
              std::declval<GeometryBuffer_t&>()))
          SystemAssembler_t;

  std::vector<SystemAssembler_t> systemAssemblers;
  systemAssemblers.reserve(numS);

  // Scattering assemblers with optimal test spaces
  std::vector<ScatteringAssembler<std::tuple<FEBasisOptimalTest>,
                                  SolutionSpaces>
             > scatteringAssemblers;
  scatteringAssemblers.reserve(numS);

  // Scattering assembler with enriched test space
  ScatteringAssembler<std::tuple<FEBasisTest>,
                                  SolutionSpaces
                      > scatteringAssemblerEnriched
                          = make_DPG_ScatteringAssembler(
                                testSpaces,
                                solutionSpaces);

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
  std::vector<GeometryBuffer_t> geometryBuffers(numS);

  for(int i = 0; i < numS; ++i)
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
        make_DPGSystemAssembler(bilinearForms[i],
                                innerProducts[i],
                                geometryBuffers[i]));
    scatteringAssemblers.emplace_back(
        make_DPG_ScatteringAssembler(optimalTestSpaces[i],
                                     solutionSpaces));
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
  for(int i = 0; i < numS; ++i)
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
  for(int i = 0; i < numS; ++i)
  {
    x.emplace_back(feBasisTrace.size()
                   +feBasisInterior.size());
    xPrevious.emplace_back(feBasisTrace.size()
                   +feBasisInterior.size());
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
  for(int i = 0; i < numS; ++i)
  {
    u.emplace_back(feBasisInterior.size());
    u[i] = 0;
    uPrevious.emplace_back(feBasisInterior.size());
    uPrevious[i] = 0;
  }

  std::vector<VectorType> theta,thetaPrevious;
  theta.reserve(numS);
  thetaPrevious.reserve(numS);
  for(int i = 0; i < numS; ++i)
  {
    theta.emplace_back(feBasisTrace.size());
    theta[i] = 0;
    thetaPrevious.emplace_back(feBasisTrace.size());
    thetaPrevious[i] = 0;
  }

  VectorType diffU(feBasisInterior.size());
  diffU = 0;

  VectorType diffTheta(feBasisTrace.size());
  diffTheta = 0;

  /////////////////////////////////////////////////////////
  //  Fixed-point iterations
  /////////////////////////////////////////////////////////
  for(int n = 0; n < N; ++n)
  {
    ofs << "Iteration " << n << std::endl;
    std::cout << "Iteration " << n << std::endl << std::endl;

    /////////////////////////////////////////////////////////
    //  Update solutions
    /////////////////////////////////////////////////////////
    xPrevious = x;
    extractSolution(uPrevious, xPrevious, 0);
    extractSolution(thetaPrevious, xPrevious, feBasisInterior.size());

    /////////////////////////////////////////////////////////
    //  Assemble the systems
    /////////////////////////////////////////////////////////
    // using Domain = GridType::template Codim<0>::Geometry::GlobalCoordinate;
    //auto f = [] (const Domain& x, const Direction& s) { return 1.;};

    // loop of the discrete ordinates
    for(int i = 0; i < numS; ++i)
    {
      Direction s = sVector[i];

      auto g = make_LinearForm(
          systemAssemblers[i].getTestSearchSpaces(),
          std::make_tuple(
            make_LinearIntegralTerm
              < 0
              , LinearIntegrationType::valueFunction
              , DomainOfIntegration::interior>
              ([s,&uExact,&sVector] (const Domain& x)
                  { return f(x,s,uExact,sVector);})));
      //       [s,&f] (const Domain& x) { return f(x,s);}

      auto kernelS = std::make_tuple([s] (const Domain& x, const Direction& sIntegration) {return kernel(x,sIntegration,s);});

      systemAssemblers[i].assembleSystem(stiffnessMatrix[i], rhs[i], g);
      VectorType scattering;
      scatteringAssemblers[i].assembleScattering<0>(scattering, xPrevious, sVector, kernelS);
      rhs[i] += scattering;
      systemAssemblers[i].applyDirichletBoundary<1>
          (stiffnessMatrix[i],
           rhs[i],
           dirichletNodesInflow[i],
           rhsInflowContrib[i]);
      systemAssemblers[i].defineCharacteristicFaces<1,dim>(stiffnessMatrix[i],
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


    for(int i = 0; i < numS; ++i)
    {
      int verbosity = 0; // 0: not verbose; >0: verbose
      UMFPack<MatrixType> umfPack(stiffnessMatrix[i], verbosity);
      InverseOperatorResult statistics;
      umfPack.apply(x[i], rhs[i], statistics);
    }

    extractSolution(u, x, 0);
    extractSolution(theta, x, feBasisInterior.size());

    ////////////////////////////////////
    //  Error computation and print in VTK file
    ////////////////////////////////////
    for(int i = 0; i < numS; ++i)
    {
      Direction s = sVector[i];

      std::cout << "Direction " << i << std::endl;

      // ////////////////////////////////////////////////////////////////////////
      // //  Write result to VTK file
      // //  We need to subsample, because VTK cannot natively display
      // //  real second-order functions
      // ////////////////////////////////////////////////////////////////////////
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
      //  Error wrt exact solution
      ////////////////////////////////////
      //We build an object of type ErrorTools to study errors, residuals and do h-adaptivity
      ErrorTools errorTools = ErrorTools();
      //We compute the L2 error between the exact and the fem solutions
      auto uExactSfixed = std::make_tuple([s] (const Domain& x){ return uAnalytic(x,s);});
      double err = errorTools.computeL2error<1>(std::get<0>(solutionSpaces),u[i],uExactSfixed);
      ofs << "'Exact' error u: || u["<< i << "] - u_fem["<< i <<"] ||_L2 = " << err << std::endl;
      // We compute the a posteriori error
          // - We compute the rhs with the enriched test space ("rhs[i]=f(v_i)")
          // -- Contribution of the source term f that has an analytic expression
      auto g = make_LinearForm(
          rhsAssembler.getTestSpaces(),
          std::make_tuple(
            make_LinearIntegralTerm
              < 0
              , LinearIntegrationType::valueFunction
              , DomainOfIntegration::interior>
              ([s,&uExact,&sVector] (const Domain& x)
                  { return f(x,s,uExact,sVector);})));
      //       [s,&f] (const Domain& x) { return f(x,s);}
      rhsAssembler.assembleRhs(rhs[i], g);
          // -- Contribution of the scattering term
      auto kernelS = std::make_tuple([s] (const Domain& x, const Direction& sIntegration) {return kernel(x,sIntegration,s);});
      VectorType scattering;
      scatteringAssemblerEnriched.assembleScattering<0>(scattering, xPrevious, sVector, kernelS);
      rhs[i] += scattering;
          // - Computation of the a posteriori error
      double aposterioriErr = errorTools.aPosterioriError(
          bilinearForms[i], innerProducts[i], x[i], rhs[i]);
          //change with contribution of scattering rhs[i]
      ofs << "A posteriori estimation of || (u,trace u) - (u_fem,theta) || = " << aposterioriErr << std::endl;

      // We compute the L2 error wrt previous iterate
      for (unsigned int j=0; j<feBasisInterior.size(); j++)
        diffU[j] = u[i][j]-uPrevious[i][j];

      for (unsigned int j=0; j<feBasisTrace.size(); j++)
        diffTheta[j] = theta[i][j]-thetaPrevious[i][j];

      ofs << "Diff wrt previous iteration: " << std::endl;
      ofs << "  -> || u["<< i << "] - u_previous["<< i <<"] ||_L2 = " << diffU.two_norm() << std::endl;
      ofs << "  -> || theta["<< i << "] - theta_previous["<< i <<"] ||_L2 = " << diffTheta.two_norm() << std::endl << std::endl;

    }
    ofs << std::endl;
    std::cout << std::endl;
  }

  ofs.std::ofstream::close();

  return 0;
  }
  catch (Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
