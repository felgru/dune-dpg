#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <cmath>
#include <cstdlib> // for std::exit()
#include <iostream>

#include <array>
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

#include <dune/functions/functionspacebases/bernsteindgrefineddgnodalbasis.hh>
#include <dune/functions/functionspacebases/bernsteinbasis.hh>
#include <dune/functions/functionspacebases/bernsteindgbasis.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>

#include <dune/dpg/bilinearformfactory.hh>
#include <dune/dpg/innerproductfactory.hh>
#include <dune/dpg/linearformfactory.hh>
#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/dpg_system_assembler.hh>
#include <dune/dpg/errortools.hh>
#include <dune/dpg/functionplotter.hh>
#include <dune/dpg/functions/normedspaces.hh>
#include <dune/dpg/rhs_assembler.hh>


using namespace Dune;

// Value of the analytic solution "for the interior of the domain"
template <class Domain,class Direction>
double fInner(const Domain& x,
              const Direction& s)
{
  const FieldVector<double,2> c{1., 1.};
  return std::expm1(c[0]*x[0])*std::expm1(c[1]*x[1]); //v pure transport
  // return 1-(x[0]-0.5)*(x[0]-0.5)-(x[1]-0.5)*(x[1]-0.5); //v RT
}
// Partial derivative of fInner with respect to x[0]
template <class Domain,class Direction>
double fInnerD0(const Domain& x,
                const Direction& s)
{
  const FieldVector<double,2> c{1., 1.};
  return c[0]*std::exp(c[0]*x[0])*std::expm1(c[1]*x[1]); //v pure transport
  // return -2*(x[0]-0.5); //v RT
}
// Partial derivative of fInner with respect to x[1]
template <class Domain,class Direction>
double fInnerD1(const Domain& x,
                const Direction& s)
{
  const FieldVector<double,2> c{1., 1.};
  return std::expm1(c[0]*x[0])*std::exp(c[1]*x[1])*c[1]; //v pure transport
  // return -2*(x[1]-0.5); //v RT
}

// This function satifies the zero incoming flux bounday conditions
template <class Domain,class Direction>
double fBoundary(const Domain& x,
                 const Direction& s)
{
  return 1.; //v pure transport
  // return ( (s[0]>0)*x[0] + (s[0]==0)*1. + (s[0]<0)*(1-x[0]) ) *
  //        ( (s[1]>0)*x[1] + (s[1]==0)*1. + (s[1]<0)*(1-x[1]) ); //v RT
}

// Partial derivative of fBoundary with respect to x[0]
template <class Domain,class Direction>
double fBoundaryD0(const Domain& x,
                   const Direction& s)
{
  return 0.; //v pure transport
  // return ( (s[0]>0)*1 + (s[0]==0)*0. + (s[0]<0)*(-1.) ) *
  //        ( (s[1]>0)*x[1] + (s[1]==0)*1. + (s[1]<0)*(1-x[1]) ); //v RT
}

// Partial derivative of fBoundary with respect to x[1]
template <class Domain,class Direction>
double fBoundaryD1(const Domain& x,
                   const Direction& s)
{
  return 0.; //v pure transport
  // return ( (s[0]>0)*x[0] + (s[0]==0)*1. + (s[0]<0)*(1-x[0]) )*
  //        ( (s[1]>0)*1 + (s[1]==0)*0. + (s[1]<0)*(-1.) ); //v RT
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
  using InnerProduct = decltype(l2InnerProduct);
  using WrappedSpaces = typename InnerProduct::TestSpaces;
  using NormedSpace = std::conditional_t<
    is_RefinedFiniteElement<std::tuple_element_t<0, WrappedSpaces>>::value,
    Functions::NormalizedRefinedBasis<InnerProduct>,
    Functions::NormalizedBasis<InnerProduct>>;

  return std::make_shared<std::tuple<NormedSpace, FEBasisTrace>>(
      std::make_tuple(NormedSpace(l2InnerProduct), FEBasisTrace(gridView)));
}

int main(int argc, char** argv)
{
  if(argc != 4 && argc != 2) {
    std::cerr << "Usage: " << argv[0] << " n [βx βy]\n\n"
              << "Solves the transport problem β.∇ϕ + c ϕ = 1"
                 " with direction β=(βx, βy) on an nxn grid.\n"
              << "Direction β will be automatically normalized.\n\n"
              << "When unspecified, β=(cos(π/8), sin(π/8))."
              << std::endl;
    std::exit(1);
  }

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  constexpr int dim = 2;
  typedef UGGrid<dim> GridType;

  unsigned int nelements = atoi(argv[1]);

  if(nelements==0) {
    std::cerr << "n has to be nonzero." << std::endl;
    std::exit(1);
  }

  FieldVector<double,dim> lower = {0,0};
  FieldVector<double,dim> upper = {1,1};
  std::array<unsigned int,dim> elements = {nelements,nelements};

  // std::unique_ptr<GridType> grid = StructuredGridFactory<GridType>::createCubeGrid(lower, upper, elements);

  std::unique_ptr<GridType> grid = StructuredGridFactory<GridType>::createSimplexGrid(lower, upper, elements);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid->leafGridView();

  ////////////////////////////////////////////////////
  //   Direction of propagation beta and coefficient c
  ////////////////////////////////////////////////////

  const double c = sigma;
  FieldVector<double, dim> beta
               = {std::cos(boost::math::constants::pi<double>()/8),
                  std::sin(boost::math::constants::pi<double>()/8)};

  if(argc==4) {
    // direction beta
    double betaX = atof(argv[2]);
    double betaY = atof(argv[3]);
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

  std::cout << "Computing solution of the transport problem" << std::endl
            << "  β.∇ϕ + c ϕ = f in [0,1]x[0,1]" << std::endl
            << "           ϕ = 0 on boundary," << std::endl
            << "with β=(" << beta[0] << ", " << beta[1] << ")"
            << " and c=" << c << "."<< std::endl
            << "Mesh size H=1/n=" << 1./nelements << std::endl;

  const size_t maxNumRefinements = 4;
  for(size_t n = 0; n < maxNumRefinements; n++) {
    ////////////////////////////////////////////////////////////////////
    //   Choose finite element spaces and weak formulation of problem
    ////////////////////////////////////////////////////////////////////

    using FEBasisInterior = Functions::BernsteinDGBasis<GridView, 1>;
    using FEBasisTraceLifting = Functions::BernsteinBasis<GridView, 2>;

    auto solutionSpaces
      = make_solution_spaces<FEBasisInterior, FEBasisTraceLifting>(gridView);

    using FEBasisTest
        = Functions::BernsteinDGRefinedDGBasis<GridView, 2, 3>;
    auto unnormalizedTestSearchSpaces = make_space_tuple<FEBasisTest>(gridView);

    // enriched test space for error estimation
    using FEBasisTest_aposteriori
        = Functions::BernsteinDGRefinedDGBasis<GridView, 2, 4>;
    auto unnormalizedTestSpaces_aposteriori
        = make_space_tuple<FEBasisTest_aposteriori>(gridView);

    auto unnormalizedInnerProduct
      = innerProductWithSpace(unnormalizedTestSearchSpaces)
        .addIntegralTerm<0,0,IntegrationType::valueValue,
                             DomainOfIntegration::interior>(1.)
        .addIntegralTerm<0,0,IntegrationType::gradGrad,
                             DomainOfIntegration::interior>(1., beta)
        .create();
    auto unnormalizedInnerProduct_aposteriori
        = replaceTestSpaces(unnormalizedInnerProduct,
                            unnormalizedTestSpaces_aposteriori);

    auto testSearchSpaces
        = make_normalized_space_tuple(unnormalizedInnerProduct);
    auto testSpaces_aposteriori
        = make_normalized_space_tuple(unnormalizedInnerProduct_aposteriori);

    auto innerProduct
        = replaceTestSpaces(unnormalizedInnerProduct, testSearchSpaces);

    auto bilinearForm
      = bilinearFormWithSpaces(testSearchSpaces, solutionSpaces)
        .addIntegralTerm<0,0,IntegrationType::valueValue,
                             DomainOfIntegration::interior>(c)
        .addIntegralTerm<0,0,IntegrationType::gradValue,
                             DomainOfIntegration::interior>(-1., beta)
        .addIntegralTerm<0,1,IntegrationType::normalVector,
                             DomainOfIntegration::face>(1., beta)
        .create();

    auto bilinearForm_aposteriori
        = replaceTestSpaces(bilinearForm, testSpaces_aposteriori);
    auto innerProduct_aposteriori
        = replaceTestSpaces(innerProduct, testSpaces_aposteriori);

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

    auto rhsLambda
      = [beta](const FieldVector<double, 2>& x){ return f(x, beta); };
    auto rhsFunc
      = Functions::makeAnalyticGridViewFunction(rhsLambda, gridView);
    auto rhsFunctions
      = linearFormWithSpace(testSearchSpaces)
        .addIntegralTerm<0,LinearIntegrationType::valueFunction,
                           DomainOfIntegration::interior>(rhsFunc)
        .create();
    systemAssembler.assembleSystem(stiffnessMatrix, rhsVector, rhsFunctions);

    // Determine Dirichlet dofs for w (inflow boundary)
    {
      std::vector<bool> dirichletNodesInflow;
      BoundaryTools::getInflowBoundaryMask(std::get<1>(*solutionSpaces),
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

    VectorType x(std::get<0>(*solutionSpaces).size()
                 + std::get<1>(*solutionSpaces).size());
    x = 0;

    std::cout << "rhs size = " << rhsVector.size()
              << " matrix size = " << stiffnessMatrix.N()
                          << " x " << stiffnessMatrix.M()
              << " solution size = " << x.size() << std::endl;


    UMFPack<MatrixType> umfPack(stiffnessMatrix, 0);
    InverseOperatorResult statistics;
    umfPack.apply(x, rhsVector, statistics);

    ///////////////////
    // Compute errors
    ///////////////////

    // Error with respect to exact solution
    const double l2err
      = ErrorTools::computeL2error<2>(std::get<0>(*solutionSpaces), x,
          [beta](const FieldVector<double, dim>& x)
          { return uAnalytic(x, beta); });

    // A posteriori error
    // We compute the rhsVector in the form given by the projection approach
    auto rhsAssembler_aposteriori = make_RhsAssembler(testSpaces_aposteriori);
    auto rightHandSide_aposteriori
      = replaceTestSpaces(rhsFunctions, testSpaces_aposteriori);
    rhsAssembler_aposteriori.assembleRhs(rhsVector, rightHandSide_aposteriori);

    const double aposterioriErr
        = ErrorTools::aPosterioriError(bilinearForm_aposteriori,
                                       innerProduct_aposteriori, x, rhsVector);

    std::cout <<   "exact L2 error:     " << l2err
              << "\na posteriori error: " << aposterioriErr
              << "\nL2 / a posteriori:  " << l2err / aposterioriErr << '\n';

    //////////////////////////////////////////////////////////////////
    //  Write result to VTK files
    //////////////////////////////////////////////////////////////////
    // auto& spacePhi = std::get<0>(*solutionSpaces);
    // auto& spaceW = std::get<1>(*solutionSpaces);
    // FunctionPlotter phiPlotter("transport_solution");
    // FunctionPlotter wPlotter("transport_solution_trace");
    // phiPlotter.plot("phi", x, spacePhi, 0, 0);
    // wPlotter.plot("w", x, spaceW, 2, spacePhi.size());

    grid->globalRefine(1);
  }
}
