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

#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/dpg_system_assembler.hh>
#include <dune/dpg/errortools.hh>
#include <dune/dpg/functionplotter.hh>
#include <dune/dpg/rhs_assembler.hh>


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
  std::array<unsigned int,dim> elements = {nelements,nelements};

  // shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createCubeGrid(lower, upper, elements);

  shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createSimplexGrid(lower, upper, elements);

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
      std::abort();
    }
    const double normBeta = std::sqrt(betaX*betaX+betaY*betaY);
    betaX = betaX/normBeta;
    betaY = betaY/normBeta;
    beta = {betaX, betaY};
  }

  std::cout << "Computing solution of the transport problem" << std::endl
            << "  β.∇ϕ + c ϕ = 1 in [0,1]x[0,1]" << std::endl
            << "           ϕ = 0 on boundary," << std::endl
            << "with β=(" << beta[0] << ", " << beta[1] << ")"
            << " and c=" << c << "."<< std::endl
            << "Mesh size H=1/n=" << 1./nelements << std::endl;

  ////////////////////////////////////////////////////////////////////
  //   Choose finite element spaces and weak formulation of problem
  ////////////////////////////////////////////////////////////////////

  using FEBasisInterior = Functions::LagrangeDGBasis<GridView, 1>;
  FEBasisInterior spacePhi(gridView);

  using FEBasisTraceLifting = Functions::PQkNodalBasis<GridView, 2>;
  FEBasisTraceLifting spaceW(gridView);

  auto solutionSpaces
      = make_space_tuple<FEBasisInterior, FEBasisTraceLifting>(gridView);

  using FEBasisTest
      = Functions::PQkDGRefinedDGBasis<GridView, 1, 3>;
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

  auto systemAssembler
      = make_DPGSystemAssembler(bilinearForm, innerProduct, geometryBuffer);

  /////////////////////////////////////////////////////////
  //  Assemble the system
  /////////////////////////////////////////////////////////

  typedef BlockVector<FieldVector<double,1> > VectorType;
  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

  VectorType rhsVector;
  MatrixType stiffnessMatrix;

  auto rhsFunctions
    = make_LinearForm(testSearchSpaces,
          std::make_tuple(make_LinearIntegralTerm<0,
                                LinearIntegrationType::valueFunction,
                                DomainOfIntegration::interior>(f(beta))));
  systemAssembler.assembleSystem(stiffnessMatrix, rhsVector, rhsFunctions);

  // Determine Dirichlet dofs for w (inflow boundary)
  {
    std::vector<bool> dirichletNodesInflow;
    BoundaryTools boundaryTools = BoundaryTools();
    boundaryTools.getInflowBoundaryMask(std::get<1>(*solutionSpaces),
                                        dirichletNodesInflow,
                                        beta);
    systemAssembler.applyDirichletBoundary<1>
        (stiffnessMatrix,
         rhsVector,
         dirichletNodesInflow,
         0.);
  }

  ////////////////////////////
  //   Compute solution
  ////////////////////////////

  VectorType x(spaceW.size()
               +spacePhi.size());
  x = 0;

  std::cout << "rhs size = " << rhsVector.size()
            << " matrix size = " << stiffnessMatrix.N()
                        << " x " << stiffnessMatrix.M()
            << " solution size = " << x.size() << std::endl;


  UMFPack<MatrixType> umfPack(stiffnessMatrix, 0);
  InverseOperatorResult statistics;
  umfPack.apply(x, rhsVector, statistics);


  //////////////////////////////////////////////////////////////////
  //  Write result to VTK files
  //////////////////////////////////////////////////////////////////
  FunctionPlotter phiPlotter("transport_solution");
  FunctionPlotter wPlotter("transport_solution_trace");
  phiPlotter.plot("phi", x, spacePhi, 0, 0);
  wPlotter.plot("w", x, spaceW, 2, spacePhi.size());


  return 0;
  }
  catch (Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
