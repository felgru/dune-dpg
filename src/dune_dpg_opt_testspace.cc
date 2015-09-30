#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <cstdlib> // for std::abort()

#include <vector>
#include <array>
#include <tuple>

#include <dune/common/exceptions.hh> // We use exceptions

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/gmshreader.hh>



#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/io.hh>
#include <dune/istl/umfpack.hh>


#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pqktracenodalbasis.hh>
#include <dune/functions/functionspacebases/pqkfacenodalbasis.hh>
#include <dune/functions/functionspacebases/optimaltestbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discretescalarglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <dune/dpg/system_assembler.hh>
#include <dune/dpg/boundarytools.hh>


using namespace Dune;



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
  array<unsigned int,dim> elements = {nelements,nelements};

  //shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createCubeGrid(lower, upper, elements);

  shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createSimplexGrid(lower, upper, elements);

  //shared_ptr<GridType> grid = shared_ptr<GridType>(GmshReader<GridType>::read("irregular-square.msh"));

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid->leafGridView();

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  typedef Functions::LagrangeDGBasis<GridView, 1> FEBasisInterior; // u
  FEBasisInterior feBasisInterior(gridView);

  typedef Functions::PQkNodalBasis<GridView, 2> FEBasisTrace; // bulk term corresponding to u^
  FEBasisTrace feBasisTrace(gridView);

  auto solutionSpaces = std::make_tuple(FEBasisInterior(gridView), FEBasisTrace(gridView));

  typedef Functions::LagrangeDGBasis<GridView, 5> FEBasisTest;     // v enriched
  auto testSpaces = std::make_tuple(FEBasisTest(gridView));

    typedef decltype(testSpaces) TestSpaces;
    typedef decltype(solutionSpaces) SolutionSpaces;
  //atof(argv[2])
  FieldVector<double, dim> beta = {1,1};
  double c = 1;
  double delta = 0.0001;

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
                                    DomainOfIntegration::interior>(1., beta)));

  auto minInnerProduct = make_InnerProduct(solutionSpaces,
          make_tuple(
              make_IntegralTerm<1,1,IntegrationType::valueValue,              // (u^,u^)
                                    DomainOfIntegration::interior>(1),
              make_IntegralTerm<1,1,IntegrationType::gradGrad,                // (beta grad u^,beta grad u^)
                                    DomainOfIntegration::interior>(1, beta)
          ));

  typedef decltype(bilinearForm) BilinearForm;
  typedef decltype(innerProduct) InnerProduct;
  typedef decltype(minInnerProduct) MinInnerProduct;
  typedef Functions::TestspaceCoefficientMatrix<BilinearForm, InnerProduct> TestspaceCoefficientMatrix;

  TestspaceCoefficientMatrix testspaceCoefficientMatrix(bilinearForm, innerProduct);

  typedef Functions::OptimalTestBasis<TestspaceCoefficientMatrix> FEBasisOptimalTest;              // v
  auto optimalTestSpaces
          = make_tuple(FEBasisOptimalTest(testspaceCoefficientMatrix));

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
//#if 0
  using Domain = GridType::template Codim<0>::Geometry::GlobalCoordinate;

  auto rightHandSide = std::make_tuple([] (const Domain& x) { return 1.;});
  systemAssembler.assembleSystem(stiffnessMatrix, rhs, rightHandSide);

  /////////////////////////////////////////////////
  //   Choose an initial iterate
  /////////////////////////////////////////////////
  VectorType x(feBasisTrace.indexSet().size()
               +feBasisInterior.indexSet().size());
  x = 0;
  //printmatrix(std::cout , stiffnessMatrix, "stiffnessMatrixWithoutBV", "--");
  //std::cout <<"rhsWithoutBV" <<std::endl;
  //for (unsigned int i=0; i<rhs.size(); i++)
  //{
  //  std::cout <<rhs[i] <<std::endl;
  //}

  // Add minimization property for u^ on (near-)characteristic boundary if epsilon is closed to zero
  systemAssembler.applyMinimization<1, MinInnerProduct,2>
                    (stiffnessMatrix,
                     minInnerProduct,
                     beta,
                     delta);

  // Determine Dirichlet dofs for u^ (inflow boundary)
  {
    std::vector<bool> dirichletNodesInflow;
    BoundaryTools boundaryTools = BoundaryTools();
    boundaryTools.getInflowBoundaryMask(std::get<1>(solutionSpaces),
                                        dirichletNodesInflow,
                                        beta);
    systemAssembler.applyDirichletBoundarySolution<1>
        (stiffnessMatrix,
         rhs,
         dirichletNodesInflow,
         0.);
  }

  writeMatrixToMatlab(stiffnessMatrix, "transportMatrix_"+ std::to_string(nelements));

 // printmatrix(std::cout , stiffnessMatrix, "stiffnessMatrix", "--");
//  std::cout <<"rhs" <<std::endl;
//  for (unsigned int i=0; i<rhs.size(); i++)
//  {
//    std::cout <<rhs[i] <<std::endl;
//  }
  ////////////////////////////
  //   Compute solution
  ////////////////////////////

  std::cout <<"rhs size = "<< rhs.size()
            <<" matrix size = " << stiffnessMatrix.N() <<" x " << stiffnessMatrix.M()
            <<" solution size = "<< x.size() <<std::endl;


//#if 0
  UMFPack<MatrixType> umfPack(stiffnessMatrix, 2);
  InverseOperatorResult statistics;
  umfPack.apply(x, rhs, statistics);
//#endif

#if 0
  // Technicality:  turn the matrix into a linear operator
  MatrixAdapter<MatrixType,VectorType,VectorType> op(stiffnessMatrix);

  // Sequential incomplete LU decomposition as the preconditioner
  SeqILU0<MatrixType,VectorType,VectorType> ilu0(stiffnessMatrix,1.0);

  // Preconditioned conjugate-gradient solver
  CGSolver<VectorType> cg(op,
                          ilu0, // preconditioner
                          1e-4, // desired residual reduction factor
                          50,   // maximum number of iterations
                          2);   // verbosity of the solver

  // Object storing some statistics about the solving process
  InverseOperatorResult statistics;

  // Solve!
  cg.apply(x, rhs, statistics);
#endif



  ////////////////////////////////////////////////////////////////////////////
  //  Make a discrete function from the FE basis and the coefficient vector
  ////////////////////////////////////////////////////////////////////////////



  VectorType u(feBasisInterior.indexSet().size());
  u=0;
  for (unsigned int i=0; i<feBasisInterior.indexSet().size(); i++)
  {
    u[i] = x[i];
  }

  VectorType uhat(feBasisTrace.indexSet().size());
  uhat=0;
  for (unsigned int i=0; i<feBasisTrace.indexSet().size(); i++)
  {
    uhat[i] = x[i+feBasisInterior.indexSet().size()];
  }

//  std::cout << u << std::endl;

  Dune::Functions::DiscreteScalarGlobalBasisFunction<decltype(feBasisInterior),decltype(u)> uFunction(feBasisInterior,u);
  auto localUFunction = localFunction(uFunction);

  Dune::Functions::DiscreteScalarGlobalBasisFunction<decltype(feBasisTrace),decltype(uhat)> uhatFunction(feBasisTrace,uhat);
  auto localUhatFunction = localFunction(uhatFunction);

  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //  We need to subsample, because VTK cannot natively display real second-order functions
  //////////////////////////////////////////////////////////////////////////////////////////////
  SubsamplingVTKWriter<GridView> vtkWriter(gridView,2);
  vtkWriter.addVertexData(localUFunction, VTK::FieldInfo("u", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write("transport_simplex_"+std::to_string(nelements) +"_"+ std::to_string(beta[0]) + "_" + std::to_string(beta[1]));

  SubsamplingVTKWriter<GridView> vtkWriter1(gridView,2);
  vtkWriter1.addVertexData(localUhatFunction, VTK::FieldInfo("uhat", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter1.write("transport_simplex_trace_"+std::to_string(nelements) +"_"+ std::to_string(beta[0]) + "_" + std::to_string(beta[1]));

//#endif

# if 0
  VectorType sol(feBasisTest.indexSet().size());
  for(unsigned int idx = 0; idx<feBasisTest.indexSet().size(); idx++)
  {
    sol=0;
  //int idx = 1;
  //int idx =atoi(argv[1]);
    sol[idx]=1;

  //std::cout << u << std::endl;
    Dune::Functions::DiscreteScalarGlobalBasisFunction<decltype(feBasisTest),decltype(sol)> uFunction(feBasisTest,sol);
    auto localUFunction = localFunction(uFunction);
  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //  We need to subsample, because VTK cannot natively display real second-order functions
  //////////////////////////////////////////////////////////////////////////////////////////////
    SubsamplingVTKWriter<GridView> vtkWriter(gridView,2);
    vtkWriter.addVertexData(localUFunction, VTK::FieldInfo("testfkt", VTK::FieldInfo::Type::scalar, 1));
    vtkWriter.write("optimal_testfunction_irregGrid_"+ std::to_string(idx));
  }
# endif

    return 0;
  }
  catch (Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
