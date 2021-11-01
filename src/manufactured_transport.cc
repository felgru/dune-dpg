#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <cmath>
#include <cstdlib> // for std::exit()
#include <iostream>
#include <unistd.h>

#include <array>
#include <chrono>
#include <memory>
#include <tuple>
#include <vector>

#include <boost/math/constants/constants.hpp>

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

#include <dune/functions/functionspacebases/bernsteindgrefineddgnodalbasis.hh>
#include <dune/functions/functionspacebases/hangingnodebernsteinp2basis.hh>
#include <dune/functions/functionspacebases/bernsteindgbasis.hh>

#include <dune/dpg/functions/analyticgridviewfunction.hh>

#include <dune/dpg/bilinearformfactory.hh>
#include <dune/dpg/innerproductfactory.hh>
#include <dune/dpg/linearformfactory.hh>
#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/dpg_system_assembler.hh>
#include <dune/dpg/errorplotter.hh>
#include <dune/dpg/errortools.hh>
#include <dune/dpg/rhs_assembler.hh>
#include <dune/dpg/functionplotter.hh>
#include <dune/dpg/functions/normalizedspaces.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <dune/subgrid/subgrid.hh>
#pragma GCC diagnostic pop

using namespace Dune;

// Value of the analytic solution "for the interior of the domain"
template <class Domain,class Direction>
double fInner(const Domain& x,
              const Direction& s)
{
  const FieldVector<double,2> c{1., 1.};
  return std::expm1(c[0]*x[0])*std::expm1(c[1]*x[1]);
}
// Partial derivative of fInner with respect to x[0]
template <class Domain,class Direction>
double fInnerD0(const Domain& x,
                const Direction& s)
{
  const FieldVector<double,2> c{1., 1.};
  return c[0]*std::exp(c[0]*x[0])*std::expm1(c[1]*x[1]);
}
// Partial derivative of fInner with respect to x[1]
template <class Domain,class Direction>
double fInnerD1(const Domain& x,
                const Direction& s)
{
  const FieldVector<double,2> c{1., 1.};
  return std::expm1(c[0]*x[0])*std::exp(c[1]*x[1])*c[1];
}

// This function satisfies the zero incoming flux boundary condition
template <class Domain,class Direction>
double fBoundary(const Domain& x,
                 const Direction& s)
{
  return 1.;
}

// Partial derivative of fBoundary with respect to x[0]
template <class Domain,class Direction>
double fBoundaryD0(const Domain& x,
                   const Direction& s)
{
  return 0.;
}

// Partial derivative of fBoundary with respect to x[1]
template <class Domain,class Direction>
double fBoundaryD1(const Domain& x,
                   const Direction& s)
{
  return 0.;
}

// Optical parameter: sigma
static const double sigma = 5.;

//The analytic solution
template <class Domain, class Direction>
double uAnalytic(const Domain& x,
                 const Direction& s)
{
  return fInner(x,s)*fBoundary(x,s);
}

//The analytic solution
template <class Domain, class Direction>
double sGradUAnalytic(const Domain& x,
                      const Direction& s)
{
  return s[0]*(fInnerD0(x,s)*fBoundary(x,s) + fInner(x,s)*fBoundaryD0(x,s)) +
         s[1]*(fInnerD1(x,s)*fBoundary(x,s) + fInner(x,s)*fBoundaryD1(x,s));
}

// The right hand-side
template <class Domain, class Direction>
double f(const Domain& x,
         const Direction& s)
{
  return sGradUAnalytic(x,s) + sigma*uAnalytic(x,s);
}

template<typename FEBasisInterior, typename FEBasisTrace>
auto make_solution_spaces(const typename FEBasisInterior::GridView& gridView)
{
  auto interiorSpace = make_space_tuple<FEBasisInterior>(gridView);
  auto l2InnerProduct
    = innerProductWithSpace(interiorSpace)
      .template addIntegralTerm<0,0,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(1.)
      .create();
  auto normedSpace = make_normalized_space(l2InnerProduct);
  using NormedSpace = decltype(normedSpace);

  return std::make_shared<std::tuple<NormedSpace, FEBasisTrace>>(
      std::make_tuple(std::move(normedSpace), FEBasisTrace(gridView)));
}

template<typename FEBasisTest>
auto
make_test_spaces(const typename FEBasisTest::GridView& gridView,
                 const FieldVector<double, 2> direction)
{
  auto unnormalizedTestSpaces = make_space_tuple<FEBasisTest>(gridView);
  auto unnormalizedInnerProduct
    = innerProductWithSpace(unnormalizedTestSpaces)
      .template addIntegralTerm<0,0,IntegrationType::gradGrad,
                                    DomainOfIntegration::interior>
                               (1., direction)
      .template addIntegralTerm<0,0,IntegrationType::travelDistanceWeighted,
                                    DomainOfIntegration::face>(1., direction)
      .create();
  return make_normalized_space_tuple(unnormalizedInnerProduct);
}

void printHelp(const char* name) {
  std::cerr << "Usage: " << name << " [-p] <n>\n"
            << "Solves the transport problem on an nxn grid.\n\n"
            << "Options:\n"
            << " -p: plot solutions and error estimates\n";
  std::exit(0);
}

int main(int argc, char** argv)
{
  bool plot = false;
  {
    int opt;
    while ((opt = getopt(argc,argv,"ph")) != EOF)
      switch(opt)
      {
        case 'p': plot = true; break;
        default:
        case '?':
        case 'h':
          printHelp(argv[0]);
      }
    if(optind != argc-1) {
      printHelp(argv[0]);
    }
  }
  const unsigned int nelements = atoi(argv[optind]);

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  constexpr int dim = 2;
  using HostGrid = UGGrid<dim>;
  using Grid = SubGrid<dim, HostGrid, false>;

  const FieldVector<double,dim> lower = {0,0};
  const FieldVector<double,dim> upper = {1,1};
  const std::array<unsigned int,dim> elements = {nelements,nelements};

  std::unique_ptr<HostGrid> hostGrid = StructuredGridFactory<HostGrid>
                                  ::createSimplexGrid(lower, upper, elements);
  hostGrid->setClosureType(HostGrid::NONE);

  // We use a SubGrid as it will automatically make sure that we do
  // not have more than difference 1 in the levels of neighboring
  // elements. This is necessary since HangingNodeLagrangeP2Basis does
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
  const double tol = 1e-10;
  for(unsigned int i = 0; err > tol && i < 200; ++i)
  {
    const auto startiteration = std::chrono::steady_clock::now();
    GridView gridView = grid->leafGridView();

    const FieldVector<double, dim> beta
               = {std::cos(boost::math::constants::pi<double>()/8),
                  std::sin(boost::math::constants::pi<double>()/8)};
    const double c = sigma;

    /////////////////////////////////////////////////////////
    //   Choose finite element spaces
    /////////////////////////////////////////////////////////

    // u
    using FEBasisInterior = Functions::BernsteinDGBasis<GridView, 1>;

    // bulk term corresponding to theta
    using FEBasisTrace = Functions::HangingNodeBernsteinP2Basis<GridView>;

    auto solutionSpaces
      = make_solution_spaces<FEBasisInterior, FEBasisTrace>(gridView);

    // v search space
    using FEBasisTest = Functions::BernsteinDGRefinedDGBasis<GridView, 1, 3>;
    auto testSpaces = make_test_spaces<FEBasisTest>(gridView, beta);

    // enriched test space for error estimation
    using FEBasisTest_aposteriori
        = Functions::BernsteinDGRefinedDGBasis<GridView, 1, 4>;
    auto testSpaces_aposteriori
        = make_test_spaces<FEBasisTest_aposteriori>(gridView, beta);

    auto innerProduct = replaceTestSpaces(
        std::get<0>(*testSpaces).preBasis().innerProduct(),
        testSpaces);
    auto innerProduct_aposteriori = replaceTestSpaces(
        std::get<0>(*testSpaces_aposteriori).preBasis().innerProduct(),
        testSpaces_aposteriori);

    auto bilinearForm
      = bilinearFormWithSpaces(testSpaces, solutionSpaces)
        .addIntegralTerm<0,0,IntegrationType::valueValue,
                             DomainOfIntegration::interior>(c)
        .addIntegralTerm<0,0,IntegrationType::gradValue,
                             DomainOfIntegration::interior>(-1., beta)
        .addIntegralTerm<0,1,IntegrationType::normalVector,
                             DomainOfIntegration::face>(1., beta)
        .create();
    auto bilinearForm_aposteriori
        = replaceTestSpaces(bilinearForm, testSpaces_aposteriori);

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
    auto rhsLambda
      = [beta](const FieldVector<double, 2>& x){ return f(x, beta); };
    auto rhsFunc
      = Functions::makeAnalyticGridViewFunctionWithQuadratureOrder<4>
                                                  (rhsLambda, gridView);
    auto rightHandSide
      = linearFormWithSpace(testSpaces)
        .addIntegralTerm<0,LinearIntegrationType::valueFunction,
                           DomainOfIntegration::interior>(rhsFunc)
        .create();

    const auto startsystemassembler = std::chrono::steady_clock::now();
    systemAssembler.assembleSystem(stiffnessMatrix, rhs, rightHandSide);

    // Determine Dirichlet dofs for theta (inflow boundary)
    {
      std::vector<bool> dirichletNodesInflow;
      BoundaryTools::getInflowBoundaryMask(std::get<1>(*solutionSpaces),
                                           dirichletNodesInflow,
                                           beta);
      systemAssembler.applyHomogeneousDirichletBoundary<1>
          (stiffnessMatrix,
           rhs,
           dirichletNodesInflow);
    }
    const auto endsystemassembler = std::chrono::steady_clock::now();
    std::cout << "The system assembler took "
              << std::chrono::duration_cast<std::chrono::microseconds>
                 (endsystemassembler - startsystemassembler).count()
              << "us.\n";


    /////////////////////////////////////////////////
    //   Choose an initial iterate
    /////////////////////////////////////////////////
    VectorType x(std::get<1>(*solutionSpaces).size()
                 + std::get<0>(*solutionSpaces).size());
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

    if(plot) {
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
      uPlotter.plot("u", x, std::get<0>(*solutionSpaces), 0, 0);
      thetaPlotter.plot("theta", x, std::get<1>(*solutionSpaces),
                        2, std::get<0>(*solutionSpaces).size());

      const auto endresults = std::chrono::steady_clock::now();
      std::cout << "Saving the results took "
                << std::chrono::duration_cast<std::chrono::microseconds>
                   (endresults - startresults).count()
                << "us.\n";
    }

    ////////////////////////////////////////////////////
    // Estimate a posteriori error and refine
    ////////////////////////////////////////////////////

    auto rhsAssembler_aposteriori = make_RhsAssembler(testSpaces_aposteriori);
    auto rightHandSide_aposteriori
      = replaceTestSpaces(rightHandSide, testSpaces_aposteriori);
    rhsAssembler_aposteriori.assembleRhs(rhs, rightHandSide_aposteriori);

    const auto starterror = std::chrono::steady_clock::now();
    const double ratio = .2;
    auto squaredErrorEstimates = ErrorTools::squaredCellwiseResidual(
                                     bilinearForm_aposteriori,
                                     innerProduct_aposteriori,
                                     x, rhs);
    if(plot) {
      ErrorPlotter errPlotter("transport_error_"
                              + std::to_string(nelements)
                              + "_" + std::to_string(i));
      errPlotter.plot("errors", squaredErrorEstimates, gridView);
    }
    err = std::sqrt(
        ErrorTools::DoerflerMarking(*grid, ratio,
                                    std::move(squaredErrorEstimates)));

    // Error with respect to exact solution
    const double l2err
      = ErrorTools::computeL2error<2>(std::get<0>(*solutionSpaces), x,
          [beta](const FieldVector<double, dim>& x)
          { return uAnalytic(x, beta); });

    std::cout << "A posteriori error in iteration " << i << ": "
              << err << std::endl;

    const auto enderror = std::chrono::steady_clock::now();
    std::cout << "The error computation took "
              << std::chrono::duration_cast<std::chrono::microseconds>
                 (enderror - starterror).count()
              << "us.\n";

    std::cout <<   "exact L2 error:     " << l2err
              << "\na posteriori error: " << err
              << "\nL2 / a posteriori:  " << l2err / err << '\n';

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
