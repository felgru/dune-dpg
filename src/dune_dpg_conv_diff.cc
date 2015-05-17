#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>

#include <vector>
#include <array>
#include <tuple>
#include <ctime>

#define FUSION_MAX_VECTOR_SIZE 15

#include <dune/common/exceptions.hh> // We use exceptions

#include <dune/grid/common/intersection.hh> //TODO necessary?
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



using namespace Dune;



// This method marks all vertices on the boundary of the grid.
// In our problem these are precisely the Dirichlet nodes.
// The result can be found in the 'dirichletNodes' variable.  There, a bit
// is set precisely when the corresponding vertex is on the grid boundary.
template <class FEBasis>
void boundaryTreatmentInflow (const FEBasis& feBasis,
                        std::vector<bool>& dirichletNodes )
{
  static const int dim = FEBasis::GridView::dimension;

  // Interpolating the identity function wrt to a Lagrange basis
  // yields the positions of the Lagrange nodes

  // TODO: We are hacking our way around the fact that interpolation
  // of vector-value functions is not supported yet.
  BlockVector<FieldVector<double,dim> > lagrangeNodes;
  interpolate(feBasis, lagrangeNodes, [](FieldVector<double,dim> x){ return x; });

  dirichletNodes.resize(lagrangeNodes.size());

  // Mark all Lagrange nodes on the bounding box as Dirichlet
  for (size_t i=0; i<lagrangeNodes.size(); i++)
  {
    bool isBoundary = false;
    for (int j=0; j<dim; j++)
      isBoundary = isBoundary || lagrangeNodes[i][j] < 1e-8;

    if (isBoundary)

      dirichletNodes[i] = true;
  }
}

/*template <class FEBasis>
void boundaryTreatment (const FEBasis& feBasis,
                        std::vector<bool>& dirichletNodes )
{
  static const int dim = FEBasis::GridView::dimension;

  // Interpolating the identity function wrt to a Lagrange basis
  // yields the positions of the Lagrange nodes

  // TODO: We are hacking our way around the fact that interpolation
  // of vector-value functions is not supported yet.
  BlockVector<FieldVector<double,dim> > lagrangeNodes;
  interpolate(feBasis, lagrangeNodes, [](FieldVector<double,dim> x){ return x; });

  dirichletNodes.resize(lagrangeNodes.size());

  // Mark all Lagrange nodes on the bounding box as Dirichlet
  for (size_t i=0; i<lagrangeNodes.size(); i++)
  {
    bool isBoundary = false;
    for (int j=0; j<dim; j++)
      isBoundary = isBoundary || lagrangeNodes[i][j] < 1e-8 || lagrangeNodes[i][j] > 1-1e-8;

    if (isBoundary)

      dirichletNodes[i] = true;
  }
}*/


int main(int argc, char** argv)
{
  try{
  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////
  time_t tstart, tassembled, tsolved;
  tstart = time(0);

  const int dim = 2;
  typedef UGGrid<dim> GridType;

  int nelements =atoi(argv[1]);

  FieldVector<double,dim> lower = {0,0};
  FieldVector<double,dim> upper = {1,1};
  array<unsigned int,dim> elements = {nelements,nelements};

  shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createCubeGrid(lower, upper, elements);

  //shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createSimplexGrid(lower, upper, elements);

  //shared_ptr<GridType> grid = shared_ptr<GridType>(GmshReader<GridType>::read("irregular-square.msh"));

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid->leafGridView();

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  typedef Functions::LagrangeDGBasis<GridView, 1> FEBasisInterior; // u
  FEBasisInterior feBasisInterior(gridView);

  typedef Functions::LagrangeDGBasis<GridView, 1> FEBasisSigma; // sigma 1 und 2
//  FEBasisSigma feBasisSigma1(gridView);
//  FEBasisSigma feBasisSigma2(gridView);

  typedef Functions::PQKTraceNodalBasis<GridView, 2> FEBasisTrace; // u^
  FEBasisTrace feBasisTrace(gridView);

  typedef Functions::PQKFaceNodalBasis<GridView, 2> FEBasisFace; // sigma_n^
//  FEBasisFace feBasisFace(gridView);

  auto solutionSpaces = std::make_tuple(FEBasisInterior(gridView),
                                        FEBasisSigma(gridView),
                                        FEBasisSigma(gridView),
                                        FEBasisTrace(gridView),
                                        FEBasisFace(gridView));

  typedef Functions::LagrangeDGBasis<GridView, 5> FEBasisTestV;   // v enriched

  typedef Functions::LagrangeDGBasis<GridView, 5> FEBasisTestTau; // tau enriched
  auto testSpaces = std::make_tuple(FEBasisTestV(gridView),
                                    FEBasisTestTau(gridView),
                                    FEBasisTestTau(gridView));

    typedef decltype(testSpaces) TestSpaces;          // testSpaces = (v, tau1, tau2)
    typedef decltype(solutionSpaces) SolutionSpaces;  // solutionSpaces =
                                                      // (u, sigma1, sigma2, u^, sigma_n^)

  FieldVector<double, dim> beta = {1,1};
  double c = 1;
  double epsilon = 5*1e-2;
  double sqrtepsilon = std::sqrt(epsilon);
  double mu = 1;
  std::cout <<"c = " << c <<" epsilon = " <<epsilon <<" sqrtepsilon = " << sqrtepsilon << " beta = [" << beta[0] <<"," << beta[1] <<"] mu = " <<mu <<std::endl;

  FieldVector<double, dim> firstcomponent = {1,0};
  FieldVector<double, dim> secondcomponent = {0,1};

  auto bilinearForm = make_BilinearForm(testSpaces, solutionSpaces,
          make_tuple(
              make_IntegralTerm<0,0,IntegrationType::valueValue,              // (cu, v))
                                    DomainOfIntegration::interior>(c),
              make_IntegralTerm<0,0,IntegrationType::gradValue,               // -(u, beta Grad v)
                                    DomainOfIntegration::interior>(-1, beta),
              make_IntegralTerm<1,0,IntegrationType::gradValue,               // sqrtepsilon (u, dx_1 tau1)
                                    DomainOfIntegration::interior>(sqrtepsilon, firstcomponent),
              make_IntegralTerm<2,0,IntegrationType::gradValue,               // sqrtepsilon (u, dx_2 tau2)
                                    DomainOfIntegration::interior>(sqrtepsilon, secondcomponent),
              make_IntegralTerm<1,1,IntegrationType::valueValue,              // (sigma1, tau1)
                                    DomainOfIntegration::interior>(1),
              make_IntegralTerm<2,2,IntegrationType::valueValue,              // (sigma2, tau2)
                                    DomainOfIntegration::interior>(1),
              make_IntegralTerm<0,1,IntegrationType::gradValue,               // sqrtepsilon (sigma1, dx_1 v)
                                    DomainOfIntegration::interior>(sqrtepsilon, firstcomponent),
              make_IntegralTerm<0,2,IntegrationType::gradValue,               // sqrtepsilon (sigma2, dx_2 v)
                                    DomainOfIntegration::interior>(sqrtepsilon, secondcomponent),
              make_IntegralTerm<0,3,IntegrationType::normalVector,              // <u^, beta n v>
                                    DomainOfIntegration::face>(1, beta),
              make_IntegralTerm<1,3,IntegrationType::normalVector,              // <u^, n_1 tau1>
                                    DomainOfIntegration::face>((-1*sqrtepsilon), firstcomponent),
              make_IntegralTerm<2,3,IntegrationType::normalVector,              // <u^, n_2 tau2>
                                    DomainOfIntegration::face>((-1*sqrtepsilon), secondcomponent),
              make_IntegralTerm<0,4,IntegrationType::normalSign,              // <sigma^ sgn(n), v>
                                    DomainOfIntegration::face>((-1*sqrtepsilon))
          ));
  auto innerProduct = make_InnerProduct(testSpaces,
          make_tuple(
              make_IntegralTerm<0,0,IntegrationType::valueValue,              // (v,v)
                                    DomainOfIntegration::interior>(1),
              make_IntegralTerm<0,0,IntegrationType::gradGrad,                // (beta grad v,beta grad v)
                                    DomainOfIntegration::interior>(1, beta),
              make_IntegralTerm<0,0,IntegrationType::gradGrad,                // (dx_1 v,dx_1 v)
                                    DomainOfIntegration::interior>(epsilon, firstcomponent),
              make_IntegralTerm<0,0,IntegrationType::gradGrad,                // (dx_2 v,dx_2 v)
                                    DomainOfIntegration::interior>(epsilon, secondcomponent),
              make_IntegralTerm<1,1,IntegrationType::valueValue,              // (tau1,tau1)
                                    DomainOfIntegration::interior>(1),
              make_IntegralTerm<2,2,IntegrationType::valueValue,              // (tau2,tau2)
                                    DomainOfIntegration::interior>(1),
              make_IntegralTerm<1,1,IntegrationType::gradGrad,              // epsilon (dx_1 tau1,dx_1 tau1)
                                    DomainOfIntegration::interior>(epsilon, firstcomponent),
              make_IntegralTerm<2,2,IntegrationType::gradGrad,              // epsilon (dx_2 tau2,dx_2 tau2)
                                    DomainOfIntegration::interior>(epsilon, secondcomponent),
              make_IntegralTerm<1,2,IntegrationType::gradGrad,              // epsilon (dx_1 tau1, dx_2 tau2)
                                    DomainOfIntegration::interior>(epsilon,
                                                                   firstcomponent,
                                                                   secondcomponent),
              make_IntegralTerm<2,1,IntegrationType::gradGrad,              // epsilon (dx_2 tau2, dx_1 tau1)
                                    DomainOfIntegration::interior>(epsilon,
                                                                   secondcomponent,
                                                                   firstcomponent)
          ));

    typedef decltype(bilinearForm) BilinearForm;
    typedef decltype(innerProduct) InnerProduct;

  typedef Functions::OptimalTestBasis<BilinearForm, InnerProduct, 0> FEBasisOptimalTest0;              // v
  FEBasisOptimalTest0 feBasisTest0(bilinearForm, innerProduct);
  typedef Functions::OptimalTestBasis<BilinearForm, InnerProduct, 1> FEBasisOptimalTest1;              // tau1
  FEBasisOptimalTest1 feBasisTest1(bilinearForm, innerProduct);
  typedef Functions::OptimalTestBasis<BilinearForm, InnerProduct, 2> FEBasisOptimalTest2;              // tau2
  FEBasisOptimalTest2 feBasisTest2(bilinearForm, innerProduct);

  auto optimalTestSpaces = make_tuple(feBasisTest0, feBasisTest1, feBasisTest2);

  auto systemAssembler = make_SystemAssembler(optimalTestSpaces, solutionSpaces,    //DPG
          bilinearForm, innerProduct, DPGFormulation());
//  auto systemAssembler = make_SystemAssembler(testSpaces, solutionSpaces,    //Saddlepoint
//          bilinearForm, innerProduct, SaddlepointFormulation());
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

  auto rightHandSide = std::make_tuple([] (const Domain& x) { return 1;},
                                       [] (const Domain& x) { return 0;},
                                       [] (const Domain& x) { return 0;});
  systemAssembler.assembleSystem(stiffnessMatrix, rhs, rightHandSide);

/*  std::ofstream filesymwobc("matrixsymmetrywithoutbc_test.txt");
  filesymwobc <<"Stiffness Matrix (without BC) Symmetry:" <<std::endl;
  for (unsigned int i=0; i<stiffnessMatrix.N(); i++)
  {
    for (unsigned int j=0; j<stiffnessMatrix.M(); j++)
    {
      if (stiffnessMatrix.exists(i,j))
      {
        if ((stiffnessMatrix.entry(i,j)-stiffnessMatrix.entry(j,i))>1e-8 or
            (stiffnessMatrix.entry(i,j)-stiffnessMatrix.entry(j,i))<-1e-8)
        {
          filesymwobc << " 1 ";
        }
        else
        {
          if (stiffnessMatrix.entry(i,j)>1e-8 or stiffnessMatrix.entry(i,j)<-1e-8)
          {
            filesymwobc << " 0 ";
          }
          else
          {
            filesymwobc << " . ";
          }
        }
      }
      else
      {
        filesymwobc << "   ";
      }
    }
    filesymwobc <<std::endl;
  }*/

//std::ofstream file2("matrixwithoutbc_test.txt");
//printmatrix(file2 , stiffnessMatrix, "matrix", "--");
  // Set weak boundary conditions for u^ (outflow boundary)
  systemAssembler.applyWeakBoundaryCondition<3,2>
                    (stiffnessMatrix,
                     beta,
                     mu);

//std::ofstream file1("matrixweakbc.txt");
//printmatrix(file1 , stiffnessMatrix, "matrix", "--");

  // Determine Dirichlet dofs for u^ (inflow boundary)
  {
    std::vector<bool> dirichletNodesInflow;
    boundaryTreatmentInflow(std::get<3>(solutionSpaces),
                            dirichletNodesInflow);
    systemAssembler.applyDirichletBoundarySolution<3,double>
        (stiffnessMatrix,
         rhs,
         dirichletNodesInflow,
         0.);
  }

  tassembled = time(0);

/*  bool issymmetric = true;
  std::ofstream filesym("matrixsymmetry.txt");
  filesym <<"Stiffness Matrix Symmetry:" <<std::endl;
  for (unsigned int i=0; i<stiffnessMatrix.N(); i++)
  {
    for (unsigned int j=0; j<stiffnessMatrix.M(); j++)
    {
      if (stiffnessMatrix.exists(i,j))
      {
        if ((stiffnessMatrix.entry(i,j)-stiffnessMatrix.entry(j,i))>1e-8 or
            (stiffnessMatrix.entry(i,j)-stiffnessMatrix.entry(j,i))<-1e-8)
        {
          issymmetric = false;
          filesym << " 1 ";
        }
        else
        {
          if (stiffnessMatrix.entry(i,j)>1e-8 or stiffnessMatrix.entry(i,j)<-1e-8)
          {
            filesym << " 0 ";
          }
          else
          {
            filesym << " . ";
          }
        }
      }
      else
      {
        filesym << "   ";
      }
    }
    filesym <<std::endl;
  }
  std::cout <<"is symmetric = " <<issymmetric <<std::endl; */

//std::ofstream file("matrix.txt");
//printmatrix(file , stiffnessMatrix, "matrix", "--");

//printmatrix(std::cout , stiffnessMatrix, "matrix", "--");

//file <<"rhs = " <<std::endl;
//for (unsigned int i=0; i<rhs.size(); i++)
//  file <<rhs[i] <<std::endl;

//  writeMatrixToMatlab(stiffnessMatrix, "stiffnessMatrix");

  ////////////////////////////
  //   Compute solution
  ////////////////////////////
  VectorType x(rhs.size());
  x = 0;

  std::cout <<"rhs size = "<< rhs.size()
            <<" matrix size = " << stiffnessMatrix.N() <<" x " << stiffnessMatrix.M()
            <<" solution size = "<< x.size() <<std::endl;


//#if 0
  UMFPack<MatrixType> umfPack(stiffnessMatrix, 2);
  InverseOperatorResult statistics;
  umfPack.apply(x, rhs, statistics);

  tsolved = time(0);

  std::cout <<"It took " << difftime(tassembled, tstart) <<" seconds to assemble the matrix and " << difftime (tsolved, tassembled) << " seconds to solve it." <<std::endl;
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

  unsigned int shift = 0;
  //shift+= std::get<0>(testSpaces).indexSet().size();  //Saddlepoint TODO here for 3 testSpaces
  //shift+= std::get<1>(testSpaces).indexSet().size();  //Saddlepoint TODO here for 3 testSpaces
  //shift+= std::get<2>(testSpaces).indexSet().size();  //Saddlepoint TODO here for 3 testSpaces

  VectorType u(std::get<0>(solutionSpaces).indexSet().size());
  u=0;
  for (unsigned int i=0; i<std::get<0>(solutionSpaces).indexSet().size(); i++)
  {
    u[i] = x[shift+i];
  }

//  VectorType uhat(feBasisInterior.indexSet().size());
//  uhat=0;
//  for (unsigned int i=0; i<feBasisTrace.indexSet().size(); i++)
//  {
//    uhat[i] = x[i+feBasisInterior.indexSet().size()];
//  }

//  std::cout << u << std::endl;

  Dune::Functions::DiscreteScalarGlobalBasisFunction<decltype(feBasisInterior),decltype(u)> uFunction(feBasisInterior,u);
  auto localUFunction = localFunction(uFunction);

//  Dune::Functions::DiscreteScalarGlobalBasisFunction<decltype(feBasisTrace),decltype(uhat)> uhatFunction(feBasisTrace,uhat);
//  auto localUhatFunction = localFunction(uhatFunction);

  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //  We need to subsample, because VTK cannot natively display real second-order functions
  //////////////////////////////////////////////////////////////////////////////////////////////
  SubsamplingVTKWriter<GridView> vtkWriter(gridView,2);
  vtkWriter.addVertexData(localUFunction, VTK::FieldInfo("u", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write("convdiff_cube_"+ std::to_string(nelements));

//  SubsamplingVTKWriter<GridView> vtkWriter1(gridView,2);
//  vtkWriter1.addVertexData(localUhatFunction, VTK::FieldInfo("uhat", VTK::FieldInfo::Type::scalar, 1));
//  vtkWriter1.write("solution_transport_trace");

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

# if 0
  VectorType sol(feBasisTrace.indexSet().size());
  for(unsigned int idx = 0; idx<feBasisTrace.indexSet().size(); idx++)
  {
    sol=0;
  //int idx = 1;
  //int idx =atoi(argv[1]);
    sol[idx]=1;

  //std::cout << u << std::endl;
    Dune::Functions::DiscreteScalarGlobalBasisFunction<decltype(feBasisTrace),decltype(sol)> uhatbasisFunction(feBasisTrace,sol);
    auto localUhatbasisFunction = localFunction(uhatbasisFunction);
  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //  We need to subsample, because VTK cannot natively display real second-order functions
  //////////////////////////////////////////////////////////////////////////////////////////////
    SubsamplingVTKWriter<GridView> vtkWriterUhatbasis(gridView,2);
    vtkWriterUhatbasis.addVertexData(localUhatbasisFunction, VTK::FieldInfo("uhatbasisfkt", VTK::FieldInfo::Type::scalar, 1));
    vtkWriterUhatbasis.write("basisfunction_uhat_"+ std::to_string(idx));
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
