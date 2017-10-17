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
#include <dune/functions/functionspacebases/pqktracenodalbasis.hh>
#include <dune/functions/functionspacebases/pqkfacenodalbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/pqkdgrefineddgnodalbasis.hh>

#include <dune/dpg/dpg_system_assembler.hh>
#include <dune/dpg/rhs_assembler.hh>
#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/errortools.hh>


using namespace Dune;

//The analytic solution
template <class Direction, class Domain = Direction>
std::function<double(const Domain&)> uAnalytic(const Direction& s)
{
  return [s] (const Domain& x) -> double
    { double crossproduct = s[0]*x[1]-s[1]*x[0];
      // return distance to inflow boundary along s
      if(crossproduct > 0)
        return sqrt(s[1]*s[1]/(s[0]*s[0])+1)*x[0];
      else
        return sqrt(s[0]*s[0]/(s[1]*s[1])+1)*x[1];
    };
}

// The right hand-side
template <class Direction, class Domain = Direction>
std::function<double(const Domain&)> f(const Direction& s)
{
  return [] (const Domain& x) { return 1.;};
}


int main(int argc, char** argv)
{
  try{
  if(argc != 2) {
    std::cerr << "Usage: " << argv[0] << " n" << std::endl << std::endl
              << "Solves the transport problem on an nxn grid." << std::endl;
    std::abort();
  }
  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  const int dim = 2;
  typedef UGGrid<dim> GridType;

  unsigned int nelements = atoi(argv[1]);

  FieldVector<double,dim> lower = {0,0};
  FieldVector<double,dim> upper = {1,1};
  std::array<unsigned int,dim> elements = {nelements,nelements};

  // std::shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createCubeGrid(lower, upper, elements);

  std::shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createSimplexGrid(lower, upper, elements);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid->leafGridView();

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
    = make_space_tuple<FEBasisInterior, FEBasisTrace>(gridView);

  // v search space
#if LEVEL_SEARCH>0
  using FEBasisTest
      = Functions::PQkDGRefinedDGBasis<GridView, LEVEL_SEARCH, K_SEARCH>;
#else
  using FEBasisTest
      = Functions::LagrangeDGBasis<GridView, K_SEARCH>;
#endif
  auto testSpaces = make_space_tuple<FEBasisTest>(gridView);

  // enriched test space for error estimation
  using FEBasisTest_aposteriori
      = Functions::PQkDGRefinedDGBasis<GridView, LEVEL_APOSTERIORI,
                                                 K_APOSTERIORI>;
  auto testSpaces_aposteriori
      = make_space_tuple<FEBasisTest_aposteriori>(gridView);

  FieldVector<double, dim> beta
             = {cos(boost::math::constants::pi<double>()/8),
                sin(boost::math::constants::pi<double>()/8)};
  double c = 0;

  auto bilinearForm = make_BilinearForm(testSpaces, solutionSpaces,
          make_tuple(
              make_IntegralTerm<0,0,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(c),
              make_IntegralTerm<0,0,IntegrationType::gradValue,
                                    DomainOfIntegration::interior>(-1., beta),
              make_IntegralTerm<0,1,IntegrationType::normalVector,
                                    DomainOfIntegration::face>(1., beta)));
  auto bilinearForm_aposteriori
      = replaceTestSpaces(bilinearForm, testSpaces_aposteriori);
  auto innerProduct = make_InnerProduct(testSpaces,
          make_tuple(
              make_IntegralTerm<0,0,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(1.),
              make_IntegralTerm<0,0,IntegrationType::gradGrad,
                                    DomainOfIntegration::interior>(1., beta)));
  auto innerProduct_aposteriori
      = replaceTestSpaces(innerProduct, testSpaces_aposteriori);

  auto systemAssembler
     = make_DPGSystemAssembler(bilinearForm, innerProduct);

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

  auto rightHandSide
    = make_LinearForm(testSpaces,
                      std::make_tuple(make_LinearIntegralTerm<0,
                                            LinearIntegrationType::valueFunction,
                                            DomainOfIntegration::interior>(
                                 f(beta))));
  systemAssembler.assembleSystem(stiffnessMatrix, rhs, rightHandSide);

  /////////////////////////////////////////////////
  //   Choose an initial iterate
  /////////////////////////////////////////////////
  VectorType x(feBasisTrace.size()
               +feBasisInterior.size());
  x = 0;

#if LEVEL_SEARCH==0
  double delta = 1e-8;
  systemAssembler.defineCharacteristicFaces<1>
                    (stiffnessMatrix,
                     rhs,
                     beta,
                     delta);
#endif

  // Determine Dirichlet dofs for theta (inflow boundary)
  {
    std::vector<bool> dirichletNodesInflow;
    BoundaryTools::getInflowBoundaryMask(std::get<1>(*solutionSpaces),
                                         dirichletNodesInflow,
                                         beta);
    systemAssembler.applyDirichletBoundary<1>
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

  double err = ErrorTools::computeL2error<1>(std::get<0>(*solutionSpaces),
                                         u, std::make_tuple(uAnalytic(beta)));

  // We compute the a posteriori error
  auto rhsAssembler_aposteriori
    = make_RhsAssembler(testSpaces_aposteriori);
  auto rightHandSide_aposteriori
    = replaceTestSpaces(rightHandSide,
        rhsAssembler_aposteriori.getTestSpaces());
  rhsAssembler_aposteriori.assembleRhs(rhs, rightHandSide_aposteriori);

  double aposterioriErr
    = ErrorTools::aPosterioriError(bilinearForm_aposteriori,
                                  innerProduct_aposteriori,
                                  x, rhs);

  std::cout << "test search space: level=" << LEVEL_SEARCH
            << ", polynomial degree=" << K_SEARCH << std::endl;
  std::cout << "aposteriori search space: level=" << LEVEL_APOSTERIORI
            << ", polynomial degree=" << K_APOSTERIORI << std::endl;
  std::cout << "grid size n=" << nelements << std::endl;
  std::cout << "L^2 error of numerical solution: " << err << std::endl;
  std::cout << "A posteriori estimation of || (u,trace u) - (u_fem,theta) || = "
            << aposterioriErr << std::endl;

  std::string filename
      = std::string("convergence_error_ls")
      + std::to_string(LEVEL_SEARCH)
      + std::string("_ks") + std::to_string(K_SEARCH)
      + std::string("_la") + std::to_string(LEVEL_APOSTERIORI)
      + std::string("_ka") + std::to_string(K_APOSTERIORI);

  std::ofstream ofs(filename, std::ios_base::app);
  ofs << LEVEL_SEARCH << " " << K_SEARCH << " "
      << LEVEL_APOSTERIORI << " " << K_APOSTERIORI << " "
      << nelements << " "
      << std::scientific << std::setprecision(12) << err << " "
      << aposterioriErr << std::endl;


  return 0;
  }
  catch (Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
