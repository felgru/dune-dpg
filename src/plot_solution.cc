#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <cmath>
#include <cstdlib> // for std::exit()
#include <iostream>

#include <array>
#include <chrono>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/math/constants/constants.hpp>

#include <dune/common/exceptions.hh> // We use exceptions

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/io.hh>
#include <dune/istl/umfpack.hh>

#include <dune/functions/functionspacebases/hangingnodep2nodalbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/pqkdgrefineddgnodalbasis.hh>

#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/dpg_system_assembler.hh>
#include <dune/dpg/errorplotter.hh>
#include <dune/dpg/errortools.hh>
#include <dune/dpg/rhs_assembler.hh>
#include <dune/dpg/functionplotter.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <dune/subgrid/subgrid.hh>
#pragma GCC diagnostic pop

using namespace Dune;

// The right hand-side
template <class Direction, class Domain = Direction>
auto f(const Direction& s)
{
  return [] (const Domain& x) { return 1.;};
}


int main(int argc, char** argv)
{
  try{
  if(argc != 2) {
    std::cerr << "Usage: " << argv[0] << " n" << std::endl << std::endl
              << "Solves the transport problem on an nxn grid." << std::endl;
    std::exit(1);
  }
  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  constexpr int dim = 2;
  using HostGrid = UGGrid<dim>;
  using Grid = SubGrid<dim, HostGrid, false>;

  const unsigned int nelements = atoi(argv[1]);

  const FieldVector<double,dim> lower = {0, 0};
  const FieldVector<double,dim> upper = {1, 1};
  const std::array<unsigned int,dim> elements = {nelements, nelements};

  // std::shared_ptr<HostGrid> hostGrid = StructuredGridFactory<HostGrid>
  //                                 ::createCubeGrid(lower, upper, elements);

  std::shared_ptr<HostGrid> hostGrid = StructuredGridFactory<HostGrid>
                                  ::createSimplexGrid(lower, upper, elements);
  hostGrid->setClosureType(HostGrid::NONE);

  // We use a SubGrid as it will automatically make sure that we do
  // not have more than difference 1 in the levels of neighboring
  // elements. This is necessary since HangingNodeP2NodalBasis does
  // not implement higher order hanging nodes constraints.
  std::unique_ptr<Grid> grid = std::make_unique<Grid>(*hostGrid);
  {
    grid->createBegin();
    grid->insertLevel(hostGrid->maxLevel());
    grid->createEnd();
    grid->setMaxLevelDifference(1);
  }

  using GridView = Grid::LeafGridView;
  //TODO how to define this without GridView?
  using GeometryBuffer = GeometryBuffer<GridView::template Codim<0>::Geometry>;
  GeometryBuffer geometryBuffer;

  double err = 1.;
  constexpr double tol = 1e-10;
  for(unsigned int i = 0; err > tol && i < 200; ++i)
  {
    const auto startiteration = std::chrono::steady_clock::now();
    GridView gridView = grid->leafGridView();

    /////////////////////////////////////////////////////////
    //   Choose finite element spaces
    /////////////////////////////////////////////////////////

    // u
    using FEBasisInterior = Functions::LagrangeDGBasis<GridView, 1>;

    // bulk term corresponding to theta
    using FEBasisTrace = Functions::HangingNodeP2NodalBasis<GridView>;

    auto solutionSpaces
      = make_space_tuple<FEBasisInterior, FEBasisTrace>(gridView);

    // v search space
    using FEBasisTest = Functions::PQkDGRefinedDGBasis<GridView, 1, 3>;
    auto testSpaces = make_space_tuple<FEBasisTest>(gridView);

    // enriched test space for error estimation
    using FEBasisTest_aposteriori
        = Functions::PQkDGRefinedDGBasis<GridView, 1, 4>;
    auto testSpaces_aposteriori
        = make_space_tuple<FEBasisTest_aposteriori>(gridView);

    FieldVector<double, dim> beta
               = {std::cos(boost::math::constants::pi<double>()/8),
                  std::sin(boost::math::constants::pi<double>()/8)};
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
              make_IntegralTerm<0,0,IntegrationType::gradGrad,
                                    DomainOfIntegration::interior>(1., beta),
              make_IntegralTerm<0,0,IntegrationType::travelDistanceWeighted,
                                    DomainOfIntegration::face>(1., beta)));
     auto innerProduct_aposteriori
        = replaceTestSpaces(innerProduct, testSpaces_aposteriori);

    //  System assembler without geometry buffer
    //auto systemAssembler
    //   = make_DPGSystemAssembler(bilinearForm, innerProduct);

    //  System assembler with geometry buffer
    auto systemAssembler
     = make_DPGSystemAssembler(bilinearForm, innerProduct, geometryBuffer);
    /////////////////////////////////////////////////////////
    //   Stiffness matrix and right hand side vector
    /////////////////////////////////////////////////////////


    typedef BlockVector<FieldVector<double,1> > VectorType;
    typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

    VectorType rhs, rhs1;
    MatrixType stiffnessMatrix, stiffnessMatrix1;

    /////////////////////////////////////////////////////////
    //  Assemble the system
    /////////////////////////////////////////////////////////
    auto rightHandSide
      = make_LinearForm(testSpaces,
                        std::make_tuple(make_LinearIntegralTerm<0,
                                            LinearIntegrationType::valueFunction,
                                            DomainOfIntegration::interior>(
                                   f(beta))));

    const auto startsystemassembler = std::chrono::steady_clock::now();
    systemAssembler.assembleSystem(stiffnessMatrix, rhs, rightHandSide);

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
    const auto endsystemassembler = std::chrono::steady_clock::now();
    std::cout << "The system assembler took "
              << std::chrono::duration_cast<std::chrono::microseconds>
                 (endsystemassembler - startsystemassembler).count()
              << "us.\n";


    /////////////////////////////////////////////////
    //   Choose an initial iterate
    /////////////////////////////////////////////////
    VectorType x(std::get<FEBasisTrace>(*solutionSpaces).size()
                 + std::get<FEBasisInterior>(*solutionSpaces).size());
    x = 0;
    ////////////////////////////
    //   Compute solution
    ////////////////////////////

    std::cout <<"rhs size = "<< rhs.size()
              <<" matrix size = " << stiffnessMatrix.N() <<" x " << stiffnessMatrix.M()
              <<" solution size = "<< x.size() <<std::endl;

    const auto startsolve = std::chrono::steady_clock::now();

    UMFPack<MatrixType> umfPack(stiffnessMatrix, 0);
    InverseOperatorResult statistics;
    umfPack.apply(x, rhs, statistics);

    const auto endsolve = std::chrono::steady_clock::now();
    std::cout << "The solution took "
              << std::chrono::duration_cast<std::chrono::microseconds>
                 (endsolve - startsolve).count()
              << "us.\n";

#if 1
    const auto startresults = std::chrono::steady_clock::now();
    //////////////////////////////////////////////////////////////////
    //  Write result to VTK file
    //////////////////////////////////////////////////////////////////
    FunctionPlotter uPlotter("transport_solution_"
                            + std::to_string(nelements)
                            + "_" + std::to_string(i));
    FunctionPlotter thetaPlotter("transport_solution_trace_"
                                + std::to_string(nelements)
                                + "_" + std::to_string(i));
    uPlotter.plot("u", x, std::get<FEBasisInterior>(*solutionSpaces), 0, 0);
    thetaPlotter.plot("theta", x, std::get<FEBasisTrace>(*solutionSpaces),
                      2, std::get<FEBasisInterior>(*solutionSpaces).size());

    const auto endresults = std::chrono::steady_clock::now();
    std::cout << "Saving the results took "
              << std::chrono::duration_cast<std::chrono::microseconds>
                 (endresults - startresults).count()
              << "us.\n";
#endif

    ////////////////////////////////////////////////////
    // Estimate a posteriori error and refine
    ////////////////////////////////////////////////////

    auto rhsAssembler_aposteriori = make_RhsAssembler(testSpaces_aposteriori);
    auto rightHandSide_aposteriori
      = replaceTestSpaces(rightHandSide, testSpaces_aposteriori);
    rhsAssembler_aposteriori.assembleRhs(rhs, rightHandSide_aposteriori);

    const auto starterror = std::chrono::steady_clock::now();
    const double ratio = .2;
    auto errorEstimates = ErrorTools::squaredCellwiseResidual(
                                     bilinearForm_aposteriori,
                                     innerProduct_aposteriori,
                                     x, rhs);
    ErrorPlotter errPlotter("transport_error_"
                            + std::to_string(nelements)
                            + "_" + std::to_string(i));
    errPlotter.plot("errors", errorEstimates, gridView);
    err = std::sqrt(
        ErrorTools::DoerflerMarking(*grid, ratio, std::move(errorEstimates)));

    std::cout << "A posteriori error in iteration " << i << ": "
              << err << std::endl;

    const auto enderror = std::chrono::steady_clock::now();
    std::cout << "The error computation took "
              << std::chrono::duration_cast<std::chrono::microseconds>
                 (enderror - starterror).count()
              << "us.\n";

    grid->preAdapt();
    grid->adapt();
    grid->postAdapt();

    const auto endwholeiteration = std::chrono::steady_clock::now();
    std::cout << "The whole iteration took "
              << std::chrono::duration_cast<std::chrono::microseconds>
                 (endwholeiteration - startiteration).count()
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
