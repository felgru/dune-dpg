#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib> // for std::exit()
#include <iostream>
#include <memory>
#include <tuple>
#include <vector>

#include <boost/math/constants/constants.hpp>

#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/io.hh>
#include <dune/istl/umfpack.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/pqkdgrefineddgnodalbasis.hh>

#include <dune/dpg/bilinearformfactory.hh>
#include <dune/dpg/innerproductfactory.hh>
#include <dune/dpg/linearformfactory.hh>
#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/dpg_system_assembler.hh>


using namespace Dune;

int main(int argc, char** argv)
{
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

  constexpr int dim = 2;
  typedef UGGrid<dim> GridType;

  const unsigned int nelements = atoi(argv[1]);

  if(nelements==0) {
    std::cerr << "n has to be nonzero." << std::endl;
    std::exit(1);
  }

  const FieldVector<double,dim> lower = {0, 0};
  const FieldVector<double,dim> upper = {1, 1};
  const std::array<unsigned int,dim> elements = {nelements, nelements};

  // std::unique_ptr<GridType> grid = StructuredGridFactory<GridType>::createCubeGrid(lower, upper, elements);

  std::unique_ptr<GridType> grid = StructuredGridFactory<GridType>::createSimplexGrid(lower, upper, elements);

  typedef GridType::LeafGridView GridView;
  const GridView gridView = grid->leafGridView();

  ////////////////////////////////////////////////////
  //   Direction of propagation beta and coefficient c
  ////////////////////////////////////////////////////

  double c = 0;
  FieldVector<double, dim> beta
               = {std::cos(boost::math::constants::pi<double>()/8),
                  std::sin(boost::math::constants::pi<double>()/8)};

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
  using FEBasisTraceLifting = Functions::LagrangeBasis<GridView, 2>;
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

  auto bilinearForm
    = bilinearFormWithSpaces(testSearchSpaces, solutionSpaces)
      .addIntegralTerm<0,0,IntegrationType::valueValue,
                           DomainOfIntegration::interior>(c)
      .addIntegralTerm<0,0,IntegrationType::gradValue,
                           DomainOfIntegration::interior>(-1., beta)
      .addIntegralTerm<0,1,IntegrationType::normalVector,
                           DomainOfIntegration::face>(1., beta)
      .create();
  auto innerProduct
    = innerProductWithSpace(testSearchSpaces)
      .addIntegralTerm<0,0,IntegrationType::valueValue,
                           DomainOfIntegration::interior>(1.)
      .addIntegralTerm<0,0,IntegrationType::gradGrad,
                           DomainOfIntegration::interior>(1., beta)
      .create();

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
    = linearFormWithSpace(testSearchSpaces)
      .addIntegralTerm<0,LinearIntegrationType::valueFunction,
                         DomainOfIntegration::interior>(1.)
      .create();

  {
    const auto startsystemassembler = std::chrono::steady_clock::now();

    VectorType rhsVector;
    MatrixType stiffnessMatrix;

    bufferedSystemAssembler
      .assembleSystem(stiffnessMatrix, rhsVector, rhsFunctions);

    const auto endsystemassembler = std::chrono::steady_clock::now();
    std::cout << "The   buffered system assembler took "
              << std::chrono::duration_cast<std::chrono::microseconds>
                    (endsystemassembler - startsystemassembler).count()
              << "us.\n";
  }

  {
    const auto startsystemassembler = std::chrono::steady_clock::now();

    VectorType rhsVector;
    MatrixType stiffnessMatrix;

    unbufferedSystemAssembler
      .assembleSystem(stiffnessMatrix, rhsVector, rhsFunctions);

    const auto endsystemassembler = std::chrono::steady_clock::now();
    std::cout << "The unbuffered system assembler took "
              << std::chrono::duration_cast<std::chrono::microseconds>
                    (endsystemassembler - startsystemassembler).count()
              << "us.\n";

#if LEVEL_SEARCH==0
    const double delta = 1e-8;
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

    const auto startsolver = std::chrono::steady_clock::now();

    UMFPack<MatrixType> umfPack(stiffnessMatrix, 0);
    InverseOperatorResult statistics;
    umfPack.apply(x, rhsVector, statistics);

    const auto endsolver = std::chrono::steady_clock::now();
    std::cout << "Solving the system with UMFPACK took "
              << std::chrono::duration_cast<std::chrono::microseconds>
                    (endsolver - startsolver).count()
              << "us.\n";
  }


  return 0;
}
