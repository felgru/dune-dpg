#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>

#include <cmath>

#include <vector>
#include <array>
#include <tuple>

#include <dune/common/exceptions.hh> // We use exceptions

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/io.hh>
#include <dune/istl/umfpack.hh>

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pqktracenodalbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/optimaltestbasis.hh>
#include <dune/functions/gridfunctions/discretescalarglobalbasisfunction.hh>

#include <dune/dpg/system_assembler.hh>
#include <dune/dpg/errortools.hh>
#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/rhs_assembler.hh>

#include <dune/grid/uggrid.hh>   // for triangular meshes that are locally adaptive
#include <dune/grid/utility/structuredgridfactory.hh> // for triangular meshes that are locally adaptive
#include <dune/grid/io/file/gmshreader.hh> // to read a mesh from an input file of type gmsh


using namespace Dune;



// This method marks all vertices on the boundary of the grid.
// In our problem these are precisely the Dirichlet nodes.
// The result can be found in the 'dirichletNodes' variable.  There, a bit
// is set precisely when the corresponding vertex is on the grid boundary.
template <class FEBasis>
void boundaryTreatmentInflow (const FEBasis& feBasis,
                        std::vector<bool>& dirichletNodes )
{
  const int dim = FEBasis::GridView::dimension;

  // Interpolating the identity function wrt to a Lagrange basis
  // yields the positions of the Lagrange nodes

  // TODO: We are hacking our way around the fact that interpolation
  // of vector-value functions is not supported yet.
  BlockVector<FieldVector<double,dim> > lagrangeNodes;
  interpolate(feBasis, lagrangeNodes, [](FieldVector<double,dim> x){ return x; });

  dirichletNodes.resize(lagrangeNodes.size());

  //std::cout<< lagrangeNodes.size() << std::endl;

  // Mark all Lagrange nodes on the bounding box as Dirichlet
  for (size_t i=0; i<lagrangeNodes.size(); i++)
  {
    bool isBoundary = false;
    for (int j=0; j<dim; j++)
    {
      //std::cout << "lagrangeNodes["<<i <<"]["<<j <<"]="<< lagrangeNodes[i][j] << std::endl;
      isBoundary = isBoundary || lagrangeNodes[i][j] < 1e-8;
    }

    if (isBoundary)
      dirichletNodes[i] = true;

    std::cout << dirichletNodes[i] << std::endl;
  }

}


// The right-hand side explicit expression
double fieldRHS(const Dune::FieldVector<double, 2>& x) {

  double beta0 = -1.;
  double beta1 = -1.;

  double c0 = 1.;
  double c1 = 1.;

  return beta0*c0*std::exp(c0*x[0])*(std::exp(c1*x[1])-std::exp(1.))
        +beta1*c1*std::exp(c1*x[1])*(std::exp(c0*x[0])-std::exp(1.))
        +(std::exp(c0*x[0])-std::exp(1.))*(std::exp(c1*x[1])-std::exp(1.)) ;
  //return 3*std::exp(x[0]+x[1])-2*(std::exp(x[0])+std::exp(x[1]))+1  ;
}

// The exact transport solution
double fieldExact(const Dune::FieldVector<double, 2>& x) {

  double c0 = 1.;
  double c1 = 1.;

  return (std::exp(c0*x[0])-std::exp(1.))*(std::exp(c1*x[1])-std::exp(1.)) ;
}

int main(int argc, char** argv)
{
  try{

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

   const int dim = 2;

   typedef UGGrid<dim> GridType;
   FieldVector<double,dim> lower = {0,0};
   FieldVector<double,dim> upper = {1,1};
   array<unsigned int,dim> elements = {10,10};

   // Square mesh
   //shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createCubeGrid(lower, upper, elements);
   // Triangular mesh
   shared_ptr<GridType> grid  =  StructuredGridFactory<GridType>::createSimplexGrid(lower, upper, elements);
   // Read mesh from an input file
   // shared_ptr<GridType> grid = shared_ptr<GridType>(GmshReader<GridType>::read("irregular-square.msh")); // for an irregular mesh square
   // shared_ptr<GridType> grid = shared_ptr<GridType>(GmshReader<GridType>::read("circle.msh")); // for an irregular mesh square

    typedef GridType::LeafGridView GridView;
    GridView gridView = grid->leafGridView();


  /////////////////////////////////////////////////////////
  //   Choose finite element spaces
  /////////////////////////////////////////////////////////

  typedef Functions::PQKTraceNodalBasis<GridView, 2> FEBasisTrace; // u^
  typedef Functions::LagrangeDGBasis<GridView, 1> FEBasisInterior; // u
  auto solutionSpaces = std::make_tuple(FEBasisTrace(gridView),
                                        FEBasisInterior(gridView));

  typedef Functions::LagrangeDGBasis<GridView, 4> FEBasisTest;     // v
  auto testSpaces = std::make_tuple(FEBasisTest(gridView));

  /////////////////////////////////////////////////////////
  //   Choose a bilinear form
  /////////////////////////////////////////////////////////

  double beta0 = -1.0;
  double beta1 = -1.0;
  FieldVector<double, dim> beta = {beta0,beta1};
  auto bilinearForm = make_BilinearForm(testSpaces, solutionSpaces,
          make_tuple(
              make_IntegralTerm<0,1,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(1.),
              make_IntegralTerm<0,1,IntegrationType::gradValue,
                                    DomainOfIntegration::interior>(-1., beta),
              make_IntegralTerm<0,0,IntegrationType::normalVector,
                                    DomainOfIntegration::face>(1., beta)));
  auto innerProduct = make_InnerProduct(testSpaces,
          make_tuple(
              make_IntegralTerm<0,0,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(1.),
              make_IntegralTerm<0,0,IntegrationType::gradGrad,
                                    DomainOfIntegration::interior>(1., beta)));

  auto rhsAssembler = make_RhsAssembler(testSpaces);

  typedef decltype(bilinearForm) BilinearForm;
  typedef decltype(innerProduct) InnerProduct;
  typedef Functions::TestspaceCoefficientMatrix<BilinearForm, InnerProduct> TestspaceCoefficientMatrix;

  TestspaceCoefficientMatrix testspaceCoefficientMatrix(bilinearForm, innerProduct);

  typedef Functions::OptimalTestBasis<TestspaceCoefficientMatrix> FEBasisOptimalTest;              // v
  FEBasisOptimalTest feBasisTest(testspaceCoefficientMatrix);
  auto optimalTestSpaces = make_tuple(feBasisTest);

  auto systemAssembler = make_SystemAssembler(optimalTestSpaces, solutionSpaces,
          bilinearForm, innerProduct, DPGFormulation());

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

  using Domain = GridType::template Codim<0>::Geometry::GlobalCoordinate;

  auto rightHandSide = std::make_tuple(fieldRHS);
  // how to retrive the value out of this:
  // std::get<i>(rightHandSide)(x)

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

    // boundaryTreatmentInflow(std::get<0>(solutionSpaces),
    //                         dirichletNodesInflow);

    BoundaryTools boundaryTools = BoundaryTools();
    boundaryTools.boundaryTreatmentInflow(std::get<0>(solutionSpaces),
                                          dirichletNodesInflow,
                                          beta);


    systemAssembler.applyDirichletBoundarySolution<0>
        (stiffnessMatrix,
         rhs,
         dirichletNodesInflow,
         0.);
  }

#if 0
  // Determine Dirichlet dofs for v (inflow boundary)
  {
    std::vector<bool> dirichletNodesInflowTest;
    boundaryTreatmentInflow(std::get<0>(testSpaces),
                            dirichletNodesInflowTest);
    systemAssembler.applyDirichletBoundary<SpaceType::test,0,double>
        (stiffnessMatrix,
         rhs,
         dirichletNodesInflowTest,
         0);
  }
#endif

  ////////////////////////////
  //   Compute solution
  ////////////////////////////

  std::cout <<"rhs size = "<< rhs.size()
            <<" matrix size = " << stiffnessMatrix.N() <<" x " << stiffnessMatrix.M()
            <<" solution size = "<< x.size() <<std::endl;

  //writeMatrixToMatlab(stiffnessMatrix, "TestMatrix1cell");

  UMFPack<MatrixType> umfPack(stiffnessMatrix, 2);
  InverseOperatorResult statistics;
  umfPack.apply(x, rhs, statistics);

  ////////////////////////////////////////////////////////////////////////////
  //  Make a discrete function from the FE basis and the coefficient vector
  ////////////////////////////////////////////////////////////////////////////

  size_t nTest = std::get<0>(testSpaces).indexSet().size();
  size_t nFace = std::get<0>(solutionSpaces).indexSet().size();
  size_t nInner = std::get<1>(solutionSpaces).indexSet().size();
  VectorType u(nInner);
  VectorType theta(nFace);
  u=0;
  theta=0;

  //We extract the solution vector u
  for (size_t i=0; i<nInner; i++)
  {
    /* TODO: use a precomputed solution space offset. */
    u[i] = x[nFace+i];
  }

  //We extract the solution vector theta of the faces
  for (size_t i=0; i<nFace; i++)
  {
    theta[i] = x[i];
  }

  auto innerSpace = std::get<1>(solutionSpaces);
  Dune::Functions::DiscreteScalarGlobalBasisFunction
      <typename std::remove_reference<decltype(innerSpace)>::type, decltype(u)>
          uFunction(innerSpace,u);
  auto localUFunction = localFunction(uFunction);


  ////////////////////////////////////////////////////////////////////////////
  //  Error evaluation
  ////////////////////////////////////////////////////////////////////////////

  std::cout << std::endl << "******** Computation of errors *************" << std::endl;
  // The exact solution against which we are comparing our FEM solution
  auto uExact = std::make_tuple(fieldExact);
  // to retrive the value out of this:
  // std::get<0>(uExact)(x)

  // Error tolerance to do h-refinement
  double adaptivityTol = 0.001;

  //We build an object of type ErrorTools to study errors, residuals and do hp-adaptivity
  ErrorTools errorTools = ErrorTools(adaptivityTol);

  //We compute the L2 error between the exact and the fem solutions
  double err = errorTools.computeL2error(innerSpace,u,uExact) ;
  std::cout << "'Exact' error u: || u - u_fem ||_L2 = " << err << std::endl ;

  //// todo: h-refinement
  //errorTools->hRefinement(grid);
  //// todo: p-refinement

  // A posteriori error
      //We compute the rhs in the form given by the projection approach
  rhsAssembler.assembleRhs(rhs, rightHandSide);
      //It is necessary to provide rhs in the above form to call this aPosterioriError method
  double aposterioriErr = errorTools.aPosterioriError(bilinearForm,innerProduct,u,theta,rhs) ;
  std::cout << "A posteriori error: || (u,trace u) - (u_fem,theta) || = " << aposterioriErr << std::endl ;

  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //  We need to subsample, because VTK cannot natively display real second-order functions
  //////////////////////////////////////////////////////////////////////////////////////////////
  SubsamplingVTKWriter<GridView> vtkWriter(gridView,0);
  vtkWriter.addVertexData(localUFunction, VTK::FieldInfo("u", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write("solution_transport");

    return 0;
  }
  catch (Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
