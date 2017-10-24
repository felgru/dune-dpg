#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib> // for std::exit()
#include <functional>
#include <iostream>
#include <tuple>
#include <vector>

#include <boost/math/constants/constants.hpp>

#include <dune/common/exceptions.hh> // We use exceptions

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
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/pqkdgrefineddgnodalbasis.hh>

#include <dune/dpg/dpg_system_assembler.hh>
#include <dune/dpg/boundarytools.hh>


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
  if(argc != 5 && argc != 2) {
    std::cerr << "Usage: " << argv[0] << " n [c βx βy]\n\n"
              << "Solves the transport problem β.∇ϕ + c ϕ = 1"
                 " with direction β=(βx, βy) on an nxn grid.\n"
              << "Direction β will be automatically normalized.\n\n"
              << "When unspecified, c=0 and β=(cos(π/8), sin(π/8))."
              << std::endl;
    std::exit(1);
  }

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  const int dim = 2;
  typedef UGGrid<dim> GridType;

  unsigned int nelements = atoi(argv[1]);

  if(nelements==0) {
    std::cerr << "n has to be nonzero." << std::endl;
    std::exit(1);
  }

  FieldVector<double,dim> lower = {0,0};
  FieldVector<double,dim> upper = {1,1};
  std::array<unsigned int,dim> elements = {nelements,nelements};

  // std::shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createCubeGrid(lower, upper, elements);

  std::shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createSimplexGrid(lower, upper, elements);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid->leafGridView();

  ////////////////////////////////////////////////////
  //   Direction of propagation beta and coefficient c
  ////////////////////////////////////////////////////

  double c = 0;
  FieldVector<double, dim> beta
               = {cos(boost::math::constants::pi<double>()/8),
                  sin(boost::math::constants::pi<double>()/8)};

  if(argc==5) {
  // coefficient c
    c = atof(argv[2]);

    // direction beta
    double betaX = atof(argv[3]);
    double betaY = atof(argv[4]);
    if(betaX==0. && betaY==0.) {
      std::cerr << "β=(βx, βy) has to be a nonzero vector."
                << std::endl;
      std::exit(1);
    }
    const double normBeta = std::sqrt(betaX*betaX+betaY*betaY);
    betaX = betaX/normBeta;
    betaY = betaY/normBeta;
    beta = {betaX, betaY};
  }

  ////////////////////////////////////////////////////////////////////
  //   Choose finite element spaces and weak formulation of problem
  ////////////////////////////////////////////////////////////////////

  // phi
  using FEBasisInterior = Functions::LagrangeDGBasis<GridView, 1>;
  using Phi = FEBasisInterior;

  // w
  using FEBasisTraceLifting = Functions::PQkNodalBasis<GridView, 2>;
  using W = FEBasisTraceLifting;

  auto solutionSpaces
    = make_space_tuple<FEBasisInterior, FEBasisTraceLifting>(gridView);

#if LEVEL_SEARCH>0
  using FEBasisTest
      = Functions::PQkDGRefinedDGBasis<GridView, LEVEL_SEARCH, K_SEARCH>;
#else
  using FEBasisTest
      = Functions::LagrangeDGBasis<GridView, K_SEARCH>;
#endif
  auto testSearchSpaces = make_space_tuple<FEBasisTest>(gridView);

  auto bilinearForm = make_BilinearForm(testSearchSpaces, solutionSpaces,
          make_tuple(
              make_IntegralTerm<0,0,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(c),
              make_IntegralTerm<0,0,IntegrationType::gradValue,
                                    DomainOfIntegration::interior>(-1., beta),
              make_IntegralTerm<0,1,IntegrationType::normalVector,
                                    DomainOfIntegration::face>(1., beta)));
  auto innerProduct = make_InnerProduct(testSearchSpaces,
          make_tuple(
              make_IntegralTerm<0,0,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(1.),
              make_IntegralTerm<0,0,IntegrationType::gradGrad,
                                    DomainOfIntegration::interior>(1., beta)));

  typedef GeometryBuffer<GridView::template Codim<0>::Geometry> GeometryBuffer;
  GeometryBuffer geometryBuffer;

  auto bufferedSystemAssembler
      = make_DPGSystemAssembler(bilinearForm, innerProduct, geometryBuffer);

  auto unbufferedSystemAssembler
      = make_DPGSystemAssembler(bilinearForm, innerProduct);

  /////////////////////////////////////////////////////////
  //  Assemble the system
  /////////////////////////////////////////////////////////

  typedef BlockVector<FieldVector<double,1> > VectorType;
  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

  auto rhsFunctions
    = make_LinearForm(testSearchSpaces,
          std::make_tuple(make_LinearIntegralTerm<0,
                                LinearIntegrationType::valueFunction,
                                DomainOfIntegration::interior>(f(beta))));

  {
    std::chrono::steady_clock::time_point startsystemassembler
      = std::chrono::steady_clock::now();

    VectorType rhsVector;
    MatrixType stiffnessMatrix;

    bufferedSystemAssembler
      .assembleSystem(stiffnessMatrix, rhsVector, rhsFunctions);

    std::chrono::steady_clock::time_point endsystemassembler
      = std::chrono::steady_clock::now();
    std::cout << "The   buffered system assembler took "
              << std::chrono::duration_cast<std::chrono::microseconds>
                    (endsystemassembler - startsystemassembler).count()
              << "us.\n";
  }

  {
    std::chrono::steady_clock::time_point startsystemassembler
      = std::chrono::steady_clock::now();

    VectorType rhsVector;
    MatrixType stiffnessMatrix;

    unbufferedSystemAssembler
      .assembleSystem(stiffnessMatrix, rhsVector, rhsFunctions);

    std::chrono::steady_clock::time_point endsystemassembler
      = std::chrono::steady_clock::now();
    std::cout << "The unbuffered system assembler took "
              << std::chrono::duration_cast<std::chrono::microseconds>
                    (endsystemassembler - startsystemassembler).count()
              << "us.\n";

#if LEVEL_SEARCH==0
    double delta = 1e-8;
    unbufferedSystemAssembler.defineCharacteristicFaces<1>
                      (stiffnessMatrix,
                       rhsVector,
                       beta,
                       delta);
#endif

    // Determine Dirichlet dofs for w (inflow boundary)
    {
      std::vector<bool> dirichletNodesInflow;
      BoundaryTools::getInflowBoundaryMask(std::get<1>(*solutionSpaces),
                                           dirichletNodesInflow,
                                           beta);
      unbufferedSystemAssembler.applyDirichletBoundary<1>
          (stiffnessMatrix,
           rhsVector,
           dirichletNodesInflow,
           0.);
    }

    ////////////////////////////
    //   Compute solution
    ////////////////////////////

    VectorType x(std::get<W>(*solutionSpaces).size()
                 + std::get<Phi>(*solutionSpaces).size());
    x = 0;

    std::chrono::steady_clock::time_point startsolver
      = std::chrono::steady_clock::now();

    UMFPack<MatrixType> umfPack(stiffnessMatrix, 0);
    InverseOperatorResult statistics;
    umfPack.apply(x, rhsVector, statistics);

    std::chrono::steady_clock::time_point endsolver
      = std::chrono::steady_clock::now();
    std::cout << "Solving the system with UMFPACK took "
              << std::chrono::duration_cast<std::chrono::microseconds>
                    (endsolver - startsolver).count()
              << "us.\n";
  }


  return 0;
  }
  catch (Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
