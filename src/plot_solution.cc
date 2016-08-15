#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <cstdlib> // for std::abort()

#include <array>
#include <functional>
#include <tuple>
#include <vector>

#include <boost/math/constants/constants.hpp>

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

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/optimaltestbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/pqkdgrefineddgnodalbasis.hh>

#include <dune/dpg/system_assembler.hh>
#include <dune/dpg/rhs_assembler.hh>
#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/errortools.hh>


using namespace Dune;

// The right hand-side
template <class Direction, class Domain = Direction>
std::function<double(const Domain&)> f(const Direction& s)
{
  return [] (const Domain& x) { return 1.;};
}


int main(int argc, char** argv)
{
  try{
  if(argc != 5) {
    std::cerr << "Usage: " << argv[0] << " n c betaX betaY" << std::endl << std::endl
              << "Solves the transport problem $beta . grad(phi) +c phi = 1$ with direction beta=(betaX,betaY) on an nxn grid." << std::endl
              << "Direction beta will be automatically normalized" << std::endl;
    std::abort();
  }

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  const int dim = 2;
  typedef UGGrid<dim> GridType;

  unsigned int nelements = atoi(argv[1]);

  if(nelements==0) {
    std::cerr << "n has to be nonzero." << std::endl;
    std::abort();
  }

  FieldVector<double,dim> lower = {0,0};
  FieldVector<double,dim> upper = {1,1};
  array<unsigned int,dim> elements = {nelements,nelements};

  // shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createCubeGrid(lower, upper, elements);

  shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createSimplexGrid(lower, upper, elements);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid->leafGridView();

  ////////////////////////////////////////////////////
  //   Direction of propagation beta and coefficient c
  ////////////////////////////////////////////////////

  // coefficient c
  double c = atof(argv[2]);

  // direction beta
  double betaX=atof(argv[3]);
  double betaY=atof(argv[4]);
  if(betaX==0. && betaY==0.) {
    std::cerr << "beta=(betaX,betaY) has to be a nonzero vector." << std::endl;
    std::abort();
  }
  double normBeta=std::sqrt(betaX*betaX+betaY*betaY);
  betaX=betaX/normBeta;
  betaY=betaY/normBeta;
  FieldVector<double, dim> beta = {betaX,betaY};

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  // u
  using FEBasisInterior = Functions::LagrangeDGBasis<GridView, 1>;
  FEBasisInterior feBasisInterior(gridView);

  // bulk term corresponding to theta
  using FEBasisTrace = Functions::PQkNodalBasis<GridView, 2>;
  FEBasisTrace feBasisTrace(gridView);

  auto solutionSpaces
      = std::make_tuple(FEBasisInterior(gridView), FEBasisTrace(gridView));

  // v search space
  using FEBasisTest
      = Functions::PQkDGRefinedDGBasis<GridView, 1, 3>;
  auto testSpaces = std::make_tuple(FEBasisTest(gridView));

  auto bilinearForm = make_BilinearForm(testSpaces, solutionSpaces,
          make_tuple(
              make_IntegralTerm<0,0,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(c),
              make_IntegralTerm<0,0,IntegrationType::gradValue,
                                    DomainOfIntegration::interior>(-1., beta),
              make_IntegralTerm<0,1,IntegrationType::normalVector,
                                    DomainOfIntegration::face>(1., beta)));
  auto innerProduct = make_InnerProduct(testSpaces,
          make_tuple(
              make_IntegralTerm<0,0,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(1.),
              make_IntegralTerm<0,0,IntegrationType::gradGrad,
                                    DomainOfIntegration::interior>(1., beta)));

  using BilinearForm = decltype(bilinearForm);
  using InnerProduct = decltype(innerProduct);

  using TestspaceCoefficientMatrix
      = Functions::TestspaceCoefficientMatrix<BilinearForm, InnerProduct>;

  TestspaceCoefficientMatrix
      testspaceCoefficientMatrix(bilinearForm, innerProduct);

  // v
  using FEBasisOptimalTest
      = Functions::OptimalTestBasis<TestspaceCoefficientMatrix>;
  auto optimalTestSpaces
          = make_tuple(FEBasisOptimalTest(testspaceCoefficientMatrix));

  auto systemAssembler
     = make_DPG_SystemAssembler(optimalTestSpaces, solutionSpaces,
                                bilinearForm);

  /////////////////////////////////////////////////////////
  //   Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////


  typedef BlockVector<FieldVector<double,1> > VectorType;
  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

  VectorType rhs;
  MatrixType stiffnessMatrix;

  /////////////////////////////////////////////////////////
  //  Assemble the system
  /////////////////////////////////////////////////////////
  auto rightHandSide = std::make_tuple(f(beta));
  systemAssembler.assembleSystem(stiffnessMatrix, rhs, rightHandSide);

  /////////////////////////////////////////////////
  //   Choose an initial iterate
  /////////////////////////////////////////////////
  VectorType x(feBasisTrace.size()
               +feBasisInterior.size());
  x = 0;

  // Determine Dirichlet dofs for theta (inflow boundary)
  {
    std::vector<bool> dirichletNodesInflow;
    BoundaryTools boundaryTools = BoundaryTools();
    boundaryTools.getInflowBoundaryMask(std::get<1>(solutionSpaces),
                                        dirichletNodesInflow,
                                        beta);
    systemAssembler.applyDirichletBoundarySolution<1>
        (stiffnessMatrix,
         rhs,
         dirichletNodesInflow,
         0.);
  }

  ////////////////////////////
  //   Compute solution
  ////////////////////////////

  std::cout <<"rhs size = "<< rhs.size()
            <<" matrix size = " << stiffnessMatrix.N() <<" x " << stiffnessMatrix.M()
            <<" solution size = "<< x.size() <<std::endl;


  UMFPack<MatrixType> umfPack(stiffnessMatrix, 0);
  InverseOperatorResult statistics;
  umfPack.apply(x, rhs, statistics);


  ////////////////////////////////////////////////////////////////////////////
  //  Make a discrete function from the FE basis and the coefficient vector
  ////////////////////////////////////////////////////////////////////////////

  VectorType u(feBasisInterior.size());
  u=0;
  for (unsigned int i=0; i<feBasisInterior.size(); i++)
  {
    u[i] = x[i];
  }

  VectorType theta(feBasisTrace.size());
  theta=0;
  for (unsigned int i=0; i<feBasisTrace.size(); i++)
  {
    theta[i] = x[i+feBasisInterior.size()];
  }

  auto uFunction
      = Dune::Functions::makeDiscreteGlobalBasisFunction<double>
            (feBasisInterior, Dune::TypeTree::hybridTreePath(), u);

  auto thetaFunction
      = Dune::Functions::makeDiscreteGlobalBasisFunction<double>
            (feBasisTrace, Dune::TypeTree::hybridTreePath(), theta);

  /////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //  We need to subsample, because VTK cannot natively display
  //  real second-order functions
  /////////////////////////////////////////////////////////////////////////
  SubsamplingVTKWriter<GridView> vtkWriter(gridView,2);
  vtkWriter.addVertexData(uFunction,
               VTK::FieldInfo("u", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write("transport_solution");

  SubsamplingVTKWriter<GridView> vtkWriter1(gridView,2);
  vtkWriter1.addVertexData(thetaFunction,
                VTK::FieldInfo("theta", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter1.write("transport_solution_trace");

  std::cout << "Solution of the transport problem" << std::endl
            << "  beta . grad(phi) +c phi = 1 in [0,1]x[0,1]" << std::endl
            << "                      phi = 0 on boundary," << std::endl
            << "with beta=(" << betaX << "," << betaY << ")"
            << " and c=" << c << "."<< std::endl
            << "Mesh size H=1/n=" << 1./nelements << std::endl;


  return 0;
  }
  catch (Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
