#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <fstream>
#include <cstdlib> // for std::abort()

#include <vector>

#include <dune/common/exceptions.hh> // We use exceptions

#include <dune/grid/common/intersection.hh> //TODO necessary?
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/gmshreader.hh>



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

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discretescalarglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <dune/dpg/system_assembler.hh>
#include <dune/dpg/errortools.hh>
#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/radiative_transfer/scattering.hh>

#include <boost/math/constants/constants.hpp>


using namespace Dune;


// Value of the analytic solution "for the interior of the domain"
template <class Domain,class Direction>
double fInner(const Domain& x,
              const Direction& s)
{
  double value = 1-(x[0]-0.5)*(x[0]-0.5)-(x[1]-0.5)*(x[1]-0.5) ;
  return value ;
}
// Partial derivative of fInner with respect to x[0]
template <class Domain,class Direction>
double fInnerD0(const Domain& x,
                const Direction& s)
{
  double value = -2*(x[0]-0.5) ;
  return value ;
}
// Partial derivative of fInner with respect to x[1]
template <class Domain,class Direction>
double fInnerD1(const Domain& x,
                const Direction& s)
{
  double value = -2*(x[1]-0.5) ;
  return value ;
}

// This function satifies the zero incoming flux bounday conditions
template <class Domain,class Direction>
double fBoundary(const Domain& x,
                 const Direction& s)
{

  double value = ( (s[0]>0)*x[0] + (s[0]==0)*1. + (s[0]<0)*(1-x[0]) )*
                 ( (s[1]>0)*x[1] + (s[1]==0)*1. + (s[1]<0)*(1-x[1]) );
  return value ;
}
// Partial derivative of fBoundary with respect to x[0]
template <class Domain,class Direction>
double fBoundaryD0(const Domain& x,
                   const Direction& s)
{

  double value = ( (s[0]>0)*1 + (s[0]==0)*0. + (s[0]<0)*(-1.) )*
                 ( (s[1]>0)*x[1] + (s[1]==0)*1. + (s[1]<0)*(1-x[1]) );
  return value ;
}
// Partial derivative of fBoundary with respect to x[1]
template <class Domain,class Direction>
double fBoundaryD1(const Domain& x,
                   const Direction& s)
{

  double value = ( (s[0]>0)*x[0] + (s[0]==0)*1. + (s[0]<0)*(1-x[0]) )*
                 ( (s[1]>0)*1 + (s[1]==0)*0. + (s[1]<0)*(-1.) );
  return value ;
}

//The analytic solution
template <class Domain,class Direction>
double uAnalytic(const Domain& x,
                 const Direction& s)
{
  return fInner(x,s)*fBoundary(x,s);
}
// double uAnalytic(const FieldVector<double, 2>& x,
//                  const FieldVector<double, 2>& s)
// {
//   return fInner(x,s)*fBoundary(x,s);
// }

// Optical parameter: sigma
template <class Domain,class Direction>
double sigma(const Domain& x,
             const Direction& s)
{
  return 3.;
}
// Optical parameter: kernel
template <class Domain,class Direction>
double kernel(const Domain& x,
              const Direction& sIntegration,
              const Direction& s)
{
  return 0.;
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
  double sum = 0. ;
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
                 sigma(x,s)*u(x,s) - collisionTerm(x,s,u,sVector) ;
  return value;
}


int main(int argc, char** argv)
{
  try{

  ///////////////////////////////////
  // To print information
  ///////////////////////////////////
  std::ofstream ofs("rad_trans_output");

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  const int dim = 2;
  typedef UGGrid<dim> GridType;

  FieldVector<double,dim> lower = {0,0};
  FieldVector<double,dim> upper = {1,1};
  array<unsigned int,dim> elements = {4,4};

  //shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createCubeGrid(lower, upper, elements);

  shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createSimplexGrid(lower, upper, elements);

  //shared_ptr<GridType> grid = shared_ptr<GridType>(GmshReader<GridType>::read("irregular-square.msh"));

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid->leafGridView();

  ///////////////////////////////////
  // Get number of discrete ordinates
  ///////////////////////////////////

  if(argc != 2) {
      std::cerr << "Usage: " << argv[0] << " <# of ordinates>"
                << std::endl;
      std::abort();
  }
  // number of discrete ordinates
  int numS = atoi(argv[1]);
  // Vector of directions
  using Domain = GridType::template Codim<0>::Geometry::GlobalCoordinate;
  using Direction = FieldVector<double, dim> ;
  std::vector< Direction > sVector(numS) ;
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

  typedef Functions::PQKTraceNodalBasis<GridView, 2> FEBasisTrace; // u^
  FEBasisTrace feBasisTrace(gridView);

  auto solutionSpaces = std::make_tuple(FEBasisInterior(gridView), FEBasisTrace(gridView));

  typedef Functions::LagrangeDGBasis<GridView, 4> FEBasisTest;     // v enriched
  auto testSpaces = std::make_tuple(FEBasisTest(gridView));

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

  typedef decltype(make_SystemAssembler(
#if 1
              std::declval<std::tuple<FEBasisOptimalTest>>(), solutionSpaces,
#else
              testSpaces, solutionSpaces,
#endif
              std::declval<BilinearForm>(),
              std::declval<InnerProduct>(), DPGFormulation()))
          SystemAssembler_t;

  if(argc != 3) {
      std::cerr << "Usage: " << argv[0] << " <# of ordinates>"
                << " <# of iterations>" << std::endl;
      std::abort();
  }
  int numS = atoi(argv[1]);
  int N = atoi(argv[2]);

  std::vector<SystemAssembler_t> systemAssemblers;
  systemAssemblers.reserve(numS);

  std::vector<ScatteringAssembler<std::tuple<FEBasisOptimalTest>,
                                  SolutionSpaces,
                                  DPGFormulation>
             > scatteringAssemblers;
  scatteringAssemblers.reserve(numS);

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

  for(int i = 0; i < numS; ++i)
  {
    using namespace boost::math::constants;
    Direction s = sVector[i];
    // FieldVector<double, dim> s = {cos(2*pi<double>()*i/numS),
    //                               sin(2*pi<double>()*i/numS)};
    bilinearForms.emplace_back(
      make_BilinearForm(testSpaces, solutionSpaces,
          make_tuple(
              make_IntegralTerm<0,0,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(2.),
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
        make_SystemAssembler(optimalTestSpaces[i], solutionSpaces,
                             bilinearForms[i], innerProducts[i],
                             DPGFormulation()));
    scatteringAssemblers.emplace_back(
        make_ScatteringAssembler(optimalTestSpaces[i],
                                 solutionSpaces,
                                 DPGFormulation()));
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
  for(int i = 0; i < numS; ++i)
  {
    Direction s = sVector[i];
    // using namespace boost::math::constants;
    // FieldVector<double, dim> s = {cos(2*pi<double>()*i/numS),
    //                               sin(2*pi<double>()*i/numS)};
    BoundaryTools boundaryTools = BoundaryTools();
    boundaryTools.boundaryTreatmentInflow(std::get<1>(solutionSpaces),
                                          dirichletNodesInflow[i],
                                          s);
  }

  /////////////////////////////////////////////////
  //   Choose an initial iterate
  /////////////////////////////////////////////////
  std::vector<VectorType> x;
  x.reserve(numS);
  for(int i = 0; i < numS; ++i)
  {
    x.emplace_back(feBasisTrace.indexSet().size()
                   +feBasisInterior.indexSet().size());
    x[i] = 0;
  }

  /////////////////////////////////////////////////////////
  //  Fixed-point iterations
  /////////////////////////////////////////////////////////
  for(int n = 0; n < N; ++n)
  {

    /////////////////////////////////////////////////////////
    //  Assemble the systems
    /////////////////////////////////////////////////////////
    // using Domain = GridType::template Codim<0>::Geometry::GlobalCoordinate;
    //auto f = [] (const Domain& x, const Direction& s) { return 1.;};

    // loop of the discrete ordinates
    for(int i = 0; i < numS; ++i)
    {
      Direction s = sVector[i];
      // using namespace boost::math::constants;
      // FieldVector<double, dim> s = {cos(2*pi<double>()*i/numS),
      //                               sin(2*pi<double>()*i/numS)};
      //auto g = std::make_tuple([s,&f] (const Domain& x) { return f(x,s);});
      auto g = std::make_tuple([s,&uExact,&sVector] (const Domain& x) { return f(x,s,uExact,sVector);});

      systemAssemblers[i].assembleSystem(stiffnessMatrix[i], rhs[i], g);
      VectorType scattering;
      scatteringAssemblers[i].assembleScattering<0>(scattering, x);
      rhs[i] += scattering;
      printvector(ofs, scattering, "scatering", "--");
      printvector(ofs, rhs[i], "rhs", "--");
      systemAssemblers[i].applyDirichletBoundarySolution<1>
          (stiffnessMatrix[i],
           rhs[i],
           dirichletNodesInflow[i],
           0.);
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
      UMFPack<MatrixType> umfPack(stiffnessMatrix[i], 2);
      InverseOperatorResult statistics;
      umfPack.apply(x[i], rhs[i], statistics);
    }



    ////////////////////////////////////
    //  Error computation and print in VTK file
    ////////////////////////////////////
    for(int i = 0; i < numS; ++i)
    {

      //////////////////////////////////////////////////////////////////////////
      //  Make a discrete function from the FE basis and the coefficient vector
      //////////////////////////////////////////////////////////////////////////
      VectorType u(feBasisInterior.indexSet().size());
      u=0;
      for (unsigned int j=0; j<feBasisInterior.indexSet().size(); j++)
      {
        u[j] = x[i][j];
      }

      Dune::Functions::DiscreteScalarGlobalBasisFunction
          <decltype(feBasisInterior),decltype(u)>
          uFunction(feBasisInterior,u);
      auto localUFunction = localFunction(uFunction);

      // ////////////////////////////////////////////////////////////////////////
      // //  Write result to VTK file
      // //  We need to subsample, because VTK cannot natively display
      // //  real second-order functions
      // ////////////////////////////////////////////////////////////////////////
      // SubsamplingVTKWriter<GridView> vtkWriter(gridView,2);
      // vtkWriter.addVertexData(localUFunction,
      //                 VTK::FieldInfo("u", VTK::FieldInfo::Type::scalar, 1));
      // std::string name = std::string("solution_rad_trans_n")
      //                  + std::to_string(n)
      //                  + std::string("_s")
      //                  + std::to_string(i);
      // vtkWriter.write(name);

      ////////////////////////////////////
      //  Error wrt exact solution
      ////////////////////////////////////
      // Error tolerance to do h-refinement: I guess we will never do this so remove
      double adaptivityTol = 0.001;
      //We build an object of type ErrorTools to study errors, residuals and do hp-adaptivity
      ErrorTools errorTools = ErrorTools(adaptivityTol);
      //We compute the L2 error between the exact and the fem solutions
      Direction s = sVector[i];
      auto uExactSfixed = std::make_tuple([s] (const Domain& x){ return uAnalytic(x,s);});
      double err = errorTools.computeL2error(std::get<1>(solutionSpaces),u,uExactSfixed);
      ofs << " 'Exact' error u: || u["<< i << "] - u_fem["<< i <<"] ||_L2 = " << err << std::endl ;
    }
    ofs << std::endl ;


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
