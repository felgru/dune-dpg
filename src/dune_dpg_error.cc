#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>

#include <cmath>

#include <algorithm>
#include <array>
#include <memory>
#include <tuple>
#include <vector>

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

#include <dune/functions/functionspacebases/hangingnodelagrangep2basis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgrefineddgbasis.hh>
#include <dune/dpg/functions/analyticgridviewfunction.hh>

#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/bilinearformfactory.hh>
#include <dune/dpg/innerproductfactory.hh>
#include <dune/dpg/linearformfactory.hh>
#include <dune/dpg/dpg_system_assembler.hh>
#include <dune/dpg/errortools.hh>
#include <dune/dpg/functionplotter.hh>
#include <dune/dpg/rhs_assembler.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <dune/subgrid/subgrid.hh>
#pragma GCC diagnostic pop


using namespace Dune;


// The right-hand side explicit expression
double fieldRHS(const Dune::FieldVector<double, 2>& x) {

  double beta0 = -1.;
  double beta1 = -1.;

  double c0 = 1.;
  double c1 = 1.;

  return beta0*c0*std::exp(c0*x[0])*(std::exp(c1*x[1])-std::exp(1.))
        +beta1*c1*std::exp(c1*x[1])*(std::exp(c0*x[0])-std::exp(1.))
        +(std::exp(c0*x[0])-std::exp(1.))*(std::exp(c1*x[1])-std::exp(1.));
  //return 3*std::exp(x[0]+x[1])-2*(std::exp(x[0])+std::exp(x[1]))+1;
}

// The exact transport solution
double fieldExact(const Dune::FieldVector<double, 2>& x) {

  double c0 = 1.;
  double c1 = 1.;

  return (std::exp(c0*x[0])-std::exp(1.))*(std::exp(c1*x[1])-std::exp(1.));
}

int main()
{
  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  constexpr int dim = 2;

  using HostGrid = UGGrid<dim>;
  using Grid = SubGrid<dim, HostGrid, false>;
  const FieldVector<double,dim> lower = {0, 0};
  const FieldVector<double,dim> upper = {1, 1};
  const std::array<unsigned int,dim> elements = {1, 1};

  // Square mesh
  //std::unique_ptr<HostGrid> hostGrid = StructuredGridFactory<HostGrid>::createCubeGrid(lower, upper, elements);
  // Triangular mesh
  std::unique_ptr<HostGrid> hostGrid = StructuredGridFactory<HostGrid>::createSimplexGrid(lower, upper, elements);
  // Read mesh from an input file
  // std::unique_ptr<HostGrid> hostGrid = GmshReader<HostGrid>::read("irregular-square.msh"); // for an irregular mesh square
  // std::unique_ptr<HostGrid> hostGrid = GmshReader<HostGrid>::read("circle.msh"); // for a circle-shaped mesh
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

  double err = 1.;
  constexpr double tol = 1e-10;
  for(unsigned int i = 0; err > tol && i < 20; ++i)
  {
    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();


    /////////////////////////////////////////////////////////
    //   Choose finite element spaces
    /////////////////////////////////////////////////////////

    using FEBasisInterior = Functions::LagrangeDGBasis<GridView, 1>;
    using FEBasisTrace = Functions::HangingNodeLagrangeP2Basis<GridView>;
    auto solutionSpaces
      = make_space_tuple<FEBasisInterior, FEBasisTrace>(gridView);

    using FEBasisTest = Functions::LagrangeDGRefinedDGBasis<GridView, 1, 3>;
    auto testSpaces = make_space_tuple<FEBasisTest>(gridView);

    // enriched test space for error estimation
    using FEBasisTest_aposteriori
        = Functions::LagrangeDGRefinedDGBasis<GridView, 1, 4>;
    auto testSpaces_aposteriori
        = make_space_tuple<FEBasisTest_aposteriori>(gridView);

    /////////////////////////////////////////////////////////
    //   Choose a bilinear form
    /////////////////////////////////////////////////////////

    const double c = 1.;
    const FieldVector<double, dim> beta = {-1., -1.};

    auto bilinearForm
      = bilinearFormWithSpaces(testSpaces, solutionSpaces)
        .addIntegralTerm<0,0,IntegrationType::valueValue,
                             DomainOfIntegration::interior>(c)
        .addIntegralTerm<0,0,IntegrationType::gradValue,
                             DomainOfIntegration::interior>(-1., beta)
        .addIntegralTerm<0,1,IntegrationType::normalVector,
                             DomainOfIntegration::face>(1., beta)
        .create();
    auto innerProduct
      = innerProductWithSpace(testSpaces)
        .addIntegralTerm<0,0,IntegrationType::gradGrad,
                             DomainOfIntegration::interior>(1., beta)
        .addIntegralTerm<0,0,IntegrationType::travelDistanceWeighted,
                             DomainOfIntegration::face>(1., beta)
        .create();

    auto bilinearForm_aposteriori
        = replaceTestSpaces(bilinearForm, testSpaces_aposteriori);
    auto innerProduct_aposteriori
        = replaceTestSpaces(innerProduct, testSpaces_aposteriori);

    auto aPosterioriInnerProduct
      = innerProductWithSpace(solutionSpaces)
        .addIntegralTerm<0,0,IntegrationType::valueValue,            // (u,u)
                             DomainOfIntegration::interior>(1.)
        .addIntegralTerm<1,1,IntegrationType::valueValue,            // (w,w)
                             DomainOfIntegration::interior>(1.)
        .addIntegralTerm<0,1,IntegrationType::valueValue,            // -2(u,w)
                             DomainOfIntegration::interior>(-2.)
        .addIntegralTerm<1,1,IntegrationType::gradGrad,              // (beta grad w,beta grad w)
                             DomainOfIntegration::interior>(1., beta)
        .addIntegralTerm<1,1,IntegrationType::valueValue,            // (cw,cw)
                             DomainOfIntegration::interior>(c*c)
        .addIntegralTerm<1,1,IntegrationType::gradValue,             // 2(beta grad w, cw)
                             DomainOfIntegration::interior>(2.*c, beta)
        .create();
    auto m2rhsFunc
      = Functions::makeAnalyticGridViewFunctionWithQuadratureOrder<4>(
        [](const FieldVector<double, dim>& x){return (-2)*fieldRHS(x);},
        gridView);
    auto m2crhsFunc
      = Functions::makeAnalyticGridViewFunctionWithQuadratureOrder<4>(
        [c](const FieldVector<double, dim>& x){return (-2)*c*fieldRHS(x);},
        gridView);
    auto aPosterioriLinearForm
      = linearFormWithSpace(solutionSpaces)
        .addIntegralTerm<1,LinearIntegrationType::gradFunction,// -2(beta grad w,f)
                           DomainOfIntegration::interior>(m2rhsFunc, beta)
        .addIntegralTerm<1,LinearIntegrationType::valueFunction,    // -2(cw,f)
                           DomainOfIntegration::interior>(m2crhsFunc)
        .create();

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

  auto rhsFunc
    = Functions::makeAnalyticGridViewFunctionWithQuadratureOrder<4>
                                                (fieldRHS, gridView);
    auto rightHandSide
      = linearFormWithSpace(testSpaces)
        .addIntegralTerm<0,LinearIntegrationType::valueFunction,
                           DomainOfIntegration::interior>(rhsFunc)
        .create();

    systemAssembler.assembleSystem(stiffnessMatrix, rhs, rightHandSide);

    /////////////////////////////////////////////////
    //   Choose an initial iterate
    /////////////////////////////////////////////////
    /* TODO: compute the correct size from the .sizes of the FE spaces. */
    VectorType x(rhs.size());
    x = 0;

    // Determine Dirichlet dofs for u^ (inflow boundary)
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

    ////////////////////////////
    //   Compute solution
    ////////////////////////////

    std::cout <<"rhs size = "<< rhs.size()
              <<" matrix size = " << stiffnessMatrix.N() <<" x " << stiffnessMatrix.M()
              <<" solution size = "<< x.size() <<std::endl;

    //writeMatrixToMatlab(stiffnessMatrix, "TestMatrix1cell");

    UMFPack<MatrixType> umfPack(stiffnessMatrix, 0);
    InverseOperatorResult statistics;
    umfPack.apply(x, rhs, statistics);

    ////////////////////////////////////////////////////////////////////////////
    //  Make a discrete function from the FE basis and the coefficient vector
    ////////////////////////////////////////////////////////////////////////////

    const size_t nFace = std::get<1>(*solutionSpaces).size();
    const size_t nInner = std::get<0>(*solutionSpaces).size();

    // We extract the solution vector u
    VectorType u(nInner);
    std::copy_n(x.begin(), nInner, u.begin());

    // We extract the solution vector theta of the faces
    VectorType theta(nFace);
    std::copy_n(x.begin() + nInner, nFace, theta.begin());

    auto innerSpace = std::get<0>(*solutionSpaces);
    auto feBasisTrace = std::get<1>(*solutionSpaces);

    ////////////////////////////////////////////////////////////////////////////
    //  Error evaluation
    ////////////////////////////////////////////////////////////////////////////

    std::cout << std::endl << "******** Computation of errors *************" << std::endl;

    // We build an object of type ErrorTools to study errors, residuals
    // and do hp-adaptivity
    // We compute the L2 error between the exact and the fem solutions
    err = ErrorTools::computeL2error<1>(innerSpace, u, fieldExact);
    std::cout << "'Exact' error u: || u - u_fem ||_L2 = " << err << std::endl;

    // A posteriori error
    // We compute the rhs in the form given by the projection approach
    auto rhsAssembler_aposteriori = make_RhsAssembler(testSpaces_aposteriori);
    auto rightHandSide_aposteriori
      = replaceTestSpaces(rightHandSide, testSpaces_aposteriori);
    rhsAssembler_aposteriori.assembleRhs(rhs, rightHandSide_aposteriori);

    const double aposterioriErr
        = ErrorTools::aPosterioriError(bilinearForm_aposteriori,
                                      innerProduct_aposteriori, x, rhs);
    std::cout << "A posteriori error: || (u,trace u) - (u_fem,theta) || = "
              << aposterioriErr << std::endl;
    const double aposterioriL2Err
        = ErrorTools::aPosterioriL2Error(aPosterioriInnerProduct,
                                        aPosterioriLinearForm, fieldRHS, x);
    std::cout << "A posteriori L2 error: || (u,trace u) - (u_fem,theta) || = "
              << aposterioriL2Err << std::endl;

    ///////////////////////////////////
    //  Write result to VTK file
    ///////////////////////////////////
    FunctionPlotter uPlotter("solution_transport"+std::to_string(i));
    uPlotter.plot("u", u, innerSpace, 0);
    FunctionPlotter thetaPlotter("solution_trace"+std::to_string(i));
    thetaPlotter.plot("theta", theta, feBasisTrace, 2);

    ////////////////
    // Refine
    ////////////////
    const double ratio = .2;
    ErrorTools::DoerflerMarking(*grid, ratio,
                            ErrorTools::squaredCellwiseResidual(
                               bilinearForm_aposteriori,
                               innerProduct_aposteriori,
                               aPosterioriInnerProduct,
                               aPosterioriLinearForm, fieldRHS,
                               x, rhs, 0));   // the last parameter is in [0,1] and
                                              // determines which error indicator
                                              // is used
                                              // 1 = residuum
                                              // 0 = |u-lifting(u^)|_L_2^2
                                              //     + conforming residuum^2
    grid->preAdapt();
    grid->adapt();
    grid->postAdapt();
  }

  return 0;
}
