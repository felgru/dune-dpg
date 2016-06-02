#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>

#include <cmath>

#include <vector>
#include <array>
#include <tuple>

#include <dune/common/exceptions.hh> // We use exceptions

#include <dune/grid/io/file/gmshreader.hh>
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
#include <dune/functions/functionspacebases/optimaltestbasis.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pqktracenodalbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dune/dpg/system_assembler.hh>
#include <dune/dpg/errortools.hh>
#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/rhs_assembler.hh>


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
  try{

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  const int dim = 2;

  typedef UGGrid<dim> GridType;
  FieldVector<double,dim> lower = {0,0};
  FieldVector<double,dim> upper = {1,1};
  array<unsigned int,dim> elements = {1,1};

  // Square mesh
  //std::shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createCubeGrid(lower, upper, elements);
  // Triangular mesh
  std::shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createSimplexGrid(lower, upper, elements);
  // Read mesh from an input file
  // std::shared_ptr<GridType> grid(GmshReader<GridType>::read("irregular-square.msh")); // for an irregular mesh square
  // std::shared_ptr<GridType> grid(GmshReader<GridType>::read("circle.msh")); // for a circle-shaped mesh

  double err = 1.;
  const double tol = 1e-10;
  for(unsigned int i = 0; err > tol && i < 20; ++i)
  {
    typedef GridType::LeafGridView GridView;
    GridView gridView = grid->leafGridView();


    /////////////////////////////////////////////////////////
    //   Choose finite element spaces
    /////////////////////////////////////////////////////////

    typedef Functions::PQkTraceNodalBasis<GridView, 2> FEBasisTrace; // u^
    typedef Functions::LagrangeDGBasis<GridView, 1> FEBasisInterior; // u
    auto solutionSpaces = std::make_tuple(FEBasisInterior(gridView),
                                          FEBasisTrace(gridView));

    typedef Functions::LagrangeDGBasis<GridView, 4> FEBasisTest;     // v
    auto testSpaces = std::make_tuple(FEBasisTest(gridView));

    /////////////////////////////////////////////////////////
    //   Choose a bilinear form
    /////////////////////////////////////////////////////////

    double beta0 = -1.0;
    double beta1 = -1.0;
    double c = 1.;
    FieldVector<double, dim> beta = {beta0,beta1};
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
                                      DomainOfIntegration::interior>(1., beta),
                make_IntegralTerm<0,0,IntegrationType::travelDistanceWeighted,
                                      DomainOfIntegration::face>(1., beta)));

    auto minInnerProduct = make_InnerProduct(solutionSpaces,
            make_tuple(
                make_IntegralTerm<1,1,IntegrationType::valueValue,              // (u^,u^)
                                      DomainOfIntegration::interior>(1),
                make_IntegralTerm<1,1,IntegrationType::gradGrad,                // (beta grad u^,beta grad u^)
                                      DomainOfIntegration::interior>(1, beta)
          ));
    auto aPosterioriInnerProduct = make_InnerProduct(solutionSpaces,
            make_tuple(
                make_IntegralTerm<0,0,IntegrationType::valueValue,              // (u,u)
                                      DomainOfIntegration::interior>(1),
                make_IntegralTerm<1,1,IntegrationType::valueValue,              // (w,w)
                                      DomainOfIntegration::interior>(1),
                make_IntegralTerm<0,1,IntegrationType::valueValue,              // -2(u,w)
                                      DomainOfIntegration::interior>(-2),
                make_IntegralTerm<1,1,IntegrationType::gradGrad,              // (beta grad w,beta grad w)
                                      DomainOfIntegration::interior>(1, beta),
                make_IntegralTerm<1,1,IntegrationType::valueValue,              // (cw,cw)
                                      DomainOfIntegration::interior>(c*c),
                make_IntegralTerm<1,1,IntegrationType::gradValue,              // 2(beta grad w, cw)
                                      DomainOfIntegration::interior>(2*c, beta)
          ));
    auto aPosterioriLinearForm = make_LinearForm(solutionSpaces,
            make_tuple(
                make_LinearIntegralTerm<1,LinearIntegrationType::gradFunction,// -2(beta grad w,f)
                                      DomainOfIntegration::interior>([](const FieldVector<double, dim>& x){return (-2)*fieldRHS(x);}, beta),
                make_LinearIntegralTerm<1,LinearIntegrationType::valueFunction,    // -2(cw,f)
                                      DomainOfIntegration::interior>([c](const FieldVector<double, dim>& x){return (-2)*c*fieldRHS(x);})
          ));
    auto rhsAssembler = make_RhsAssembler(testSpaces);

    typedef decltype(bilinearForm) BilinearForm;
    typedef decltype(innerProduct) InnerProduct;
    typedef Functions::TestspaceCoefficientMatrix<BilinearForm, InnerProduct> TestspaceCoefficientMatrix;

    TestspaceCoefficientMatrix testspaceCoefficientMatrix(bilinearForm, innerProduct);

    typedef Functions::OptimalTestBasis<TestspaceCoefficientMatrix> FEBasisOptimalTest;              // v
    FEBasisOptimalTest feBasisTest(testspaceCoefficientMatrix);
    auto optimalTestSpaces = make_tuple(feBasisTest);

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

    auto rightHandSide
      = make_DPG_LinearForm(systemAssembler.getTestSpaces(),
                    std::make_tuple(make_LinearIntegralTerm<0,
                                        LinearIntegrationType::valueFunction,
                                        DomainOfIntegration::interior>(fieldRHS)));

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
      std::vector<bool> dirichletNodesInflowErrorTools;

      BoundaryTools boundaryTools = BoundaryTools();
      boundaryTools.getInflowBoundaryMask(std::get<1>(solutionSpaces),
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

    //writeMatrixToMatlab(stiffnessMatrix, "TestMatrix1cell");

    UMFPack<MatrixType> umfPack(stiffnessMatrix, 0);
    InverseOperatorResult statistics;
    umfPack.apply(x, rhs, statistics);

    ////////////////////////////////////////////////////////////////////////////
    //  Make a discrete function from the FE basis and the coefficient vector
    ////////////////////////////////////////////////////////////////////////////

    size_t nFace = std::get<1>(solutionSpaces).size();
    size_t nInner = std::get<0>(solutionSpaces).size();
    VectorType u(nInner);
    VectorType theta(nFace);
    u=0;
    theta=0;

    // We extract the solution vector u
    for (size_t i=0; i<nInner; i++)
    {
      u[i] = x[i];
    }

    // We extract the solution vector theta of the faces
    for (size_t i=0; i<nFace; i++)
    {
      theta[i] = x[nInner+i];
    }

    auto innerSpace = std::get<0>(solutionSpaces);
    auto uFunction
        = Dune::Functions::makeDiscreteGlobalBasisFunction<double>
              (innerSpace, Dune::TypeTree::hybridTreePath(), u);
    auto localUFunction = localFunction(uFunction);

    auto feBasisTrace = std::get<1>(solutionSpaces);
    auto thetaFunction
        = Dune::Functions::makeDiscreteGlobalBasisFunction<double>
              (feBasisTrace, Dune::TypeTree::hybridTreePath(), theta);
    auto localThetaFunction = localFunction(thetaFunction);


    ////////////////////////////////////////////////////////////////////////////
    //  Error evaluation
    ////////////////////////////////////////////////////////////////////////////

    std::cout << std::endl << "******** Computation of errors *************" << std::endl;
    // The exact solution against which we are comparing our FEM solution
    auto uExact = std::make_tuple(fieldExact);
    // to retrive the value out of this:
    // std::get<0>(uExact)(x)

    // We build an object of type ErrorTools to study errors, residuals and do hp-adaptivity
    ErrorTools errorTools = ErrorTools();

    // We compute the L2 error between the exact and the fem solutions
    err = errorTools.computeL2error<1>(innerSpace,u,uExact);
    std::cout << "'Exact' error u: || u - u_fem ||_L2 = " << err << std::endl;

    // A posteriori error
    // We compute the rhs in the form given by the projection approach
    auto rightHandSideEnriched
      = make_LinearForm(rhsAssembler.getTestSpaces(),
                    std::make_tuple(make_LinearIntegralTerm<0,
                                        LinearIntegrationType::valueFunction,
                                        DomainOfIntegration::interior>(fieldRHS)));
    rhsAssembler.assembleRhs(rhs, rightHandSideEnriched);
    // It is necessary to provide rhs in the above form to call this aPosterioriError method
    double aposterioriErr
        = errorTools.aPosterioriError(bilinearForm, innerProduct, x, rhs);
    std::cout << "A posteriori error: || (u,trace u) - (u_fem,theta) || = " << aposterioriErr << std::endl;
    double aposterioriL2Err
        = errorTools.aPosterioriL2Error(aPosterioriInnerProduct, aPosterioriLinearForm, fieldRHS, x);
    std::cout << "A posteriori L2 error: || (u,trace u) - (u_fem,theta) || = " << aposterioriL2Err << std::endl;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //  Write result to VTK file
    //  We need to subsample, because VTK cannot natively display real second-order functions
    //////////////////////////////////////////////////////////////////////////////////////////////
    SubsamplingVTKWriter<GridView> vtkWriter(gridView,0);
    vtkWriter.addVertexData(localUFunction, VTK::FieldInfo("u", VTK::FieldInfo::Type::scalar, 1));
    vtkWriter.write(std::string{"solution_transport_"}+std::to_string(i));

    SubsamplingVTKWriter<GridView> vtkWriter1(gridView,2);
    vtkWriter1.addVertexData(localThetaFunction, VTK::FieldInfo("theta",VTK::FieldInfo::Type::scalar, 1));
    vtkWriter1.write(std::string{"solution_trace_"}+std::to_string(i));

    ////////////////
    // Refine
    ////////////////
    const double ratio = .2;
    errorTools.DoerflerMarking(*grid, ratio,
                               bilinearForm, innerProduct,
                               aPosterioriInnerProduct,
                               aPosterioriLinearForm, fieldRHS,
                               x, rhs, 0);    // the last parameter is in [0,1] and
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
  catch (Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
