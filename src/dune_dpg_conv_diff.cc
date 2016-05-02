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

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk.hh>
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
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pqktracenodalbasis.hh>
#include <dune/functions/functionspacebases/pqkfacenodalbasis.hh>
#include <dune/functions/functionspacebases/optimaltestbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dune/dpg/system_assembler.hh>
#include <dune/dpg/errortools.hh>
#include <dune/dpg/rhs_assembler.hh>

#include <boost/math/constants/constants.hpp>

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

  // Mark all Lagrange nodes on the bounding box as Dirichlet
  for (size_t i=0; i<lagrangeNodes.size(); i++)
  {
    bool isBoundary = false;
    //for (int j=0; j<dim; j++)
    int j = 0;
      isBoundary = isBoundary || lagrangeNodes[i][j] < 1e-8;

    if (isBoundary)

      dirichletNodes[i] = true;
  }
}

/*template <class FEBasis>
void boundaryTreatment (const FEBasis& feBasis,
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

// The right-hand side explicit expression

double epsilon;   //TODO classes for fieldRHS and fielsExact with epsilon as variable

double fieldRHS(const Dune::FieldVector<double, 2>& x) {
  const double pi = boost::math::constants::pi<double>();
  return (-(-pi*pi*epsilon*x[0] + pi*pi*epsilon - 1) * std::sin(pi*x[1]));
}

// The exact transport solution
double fieldExact(const Dune::FieldVector<double, 2>& x) {
  const double pi = boost::math::constants::pi<double>();
  const double r1 = (-1+std::sqrt(1+4*epsilon*epsilon*pi*pi))/(-2*epsilon);
  const double r2 = (-1-std::sqrt(1+4*epsilon*epsilon*pi*pi))/(-2*epsilon);
  return (((std::exp(r1*(x[0]-1))-std::exp(r2*(x[0]-1)))/(std::exp(-r1)-std::exp(-r2)) - (1-x[0]))* std::sin(pi*x[1]) );
}


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

  unsigned int nelements = atoi(argv[1]);
  int epsinv = atoi(argv[2]);
//  int muinv = atoi(argv[3]);
  if (epsinv == 0)
  {
    epsilon = 0;
  }
  else
  {
    epsilon = 1/double(epsinv);
  }
  //epsilon = 0;

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
  typedef Functions::LagrangeDGBasis<GridView, 1> FEBasisSigma; // sigma 1 und 2
  typedef Functions::PQkNodalBasis<GridView, 2> FEBasisTrace; // u^
//  typedef Functions::PQkTraceNodalBasis<GridView, 2> FEBasisTrace; // u^
  typedef Functions::PQkFaceNodalBasis<GridView, 2> FEBasisFace; // sigma_n^
  auto solutionSpaces = std::make_tuple(FEBasisInterior(gridView),
                                        FEBasisSigma(gridView),
                                        FEBasisSigma(gridView),
                                        FEBasisTrace(gridView),
                                        FEBasisFace(gridView));

  typedef Functions::LagrangeDGBasis<GridView, 4> FEBasisTestV;   // v search space
  typedef Functions::LagrangeDGBasis<GridView, 4> FEBasisTestTau; // tau search space
  auto testSpaces = std::make_tuple(FEBasisTestV(gridView),
                                    FEBasisTestTau(gridView),
                                    FEBasisTestTau(gridView));

  typedef decltype(testSpaces) TestSpaces;          // testSpaces = (v, tau1, tau2)
  typedef decltype(solutionSpaces) SolutionSpaces;  // solutionSpaces =
                                                      // (u, sigma1, sigma2, u^, sigma_n^)

//  FieldVector<double, dim> beta = {1,0};
  FieldVector<double, dim> beta = {1,0};
  const double c(0);
//  const double epsilon(0.01);
  const double sqrtepsilon (std::sqrt(epsilon));
  const double mu(epsilon*epsilon);   //1000./double(muinv)
  const double delta(1e-2);

  double sqrtEpsilonCutof = sqrtepsilon;
  if (sqrtepsilon < delta)
  {
//    sqrtEpsilonCutof = 0;
  }

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
//
              make_IntegralTerm<0,3,IntegrationType::normalVector,            // <u^, beta n v>
                                    DomainOfIntegration::face>(1, beta),
              make_IntegralTerm<1,3,IntegrationType::normalVector,            // <u^, n_1 tau1>
                                    DomainOfIntegration::face>((-1*sqrtEpsilonCutof), firstcomponent),
              make_IntegralTerm<2,3,IntegrationType::normalVector,            // <u^, n_2 tau2>
                                    DomainOfIntegration::face>((-1*sqrtEpsilonCutof), secondcomponent),
//
              make_IntegralTerm<0,4,IntegrationType::normalSign,              // -sqrtepsilon <sigma^ sgn(n), v>
                                    DomainOfIntegration::face>((-1))//((-1*sqrtepsilon) (-1) TODO which formulation makes more sense?
          ));

  auto innerProduct = make_InnerProduct(testSpaces,
          make_tuple(
              make_IntegralTerm<0,0,IntegrationType::valueValue,              // (v,v)
                                    DomainOfIntegration::interior>(1),
              make_IntegralTerm<0,0,IntegrationType::gradGrad,                // (beta grad v,beta grad v)
                                    DomainOfIntegration::interior>(1, beta),
              make_IntegralTerm<0,0,IntegrationType::gradGrad,                // epsilon (dx_1 v,dx_1 v)
                                    DomainOfIntegration::interior>(epsilon, firstcomponent),
              make_IntegralTerm<0,0,IntegrationType::gradGrad,                // epsilon (dx_2 v,dx_2 v)
                                    DomainOfIntegration::interior>(epsilon, secondcomponent),
              make_IntegralTerm<1,1,IntegrationType::valueValue,              // (tau1,tau1)
                                    DomainOfIntegration::interior>(1),
              make_IntegralTerm<2,2,IntegrationType::valueValue,              // (tau2,tau2)
                                    DomainOfIntegration::interior>(1),
              make_IntegralTerm<1,1,IntegrationType::gradGrad,                // epsilon (dx_1 tau1,dx_1 tau1)
                                    DomainOfIntegration::interior>(epsilon, firstcomponent),
              make_IntegralTerm<2,2,IntegrationType::gradGrad,                // epsilon (dx_2 tau2,dx_2 tau2)
                                    DomainOfIntegration::interior>(epsilon, secondcomponent),
              make_IntegralTerm<1,2,IntegrationType::gradGrad,                // epsilon (dx_1 tau1, dx_2 tau2)
                                    DomainOfIntegration::interior>(epsilon,
                                                                   firstcomponent,
                                                                   secondcomponent),
              make_IntegralTerm<2,1,IntegrationType::gradGrad,                // epsilon (dx_2 tau2, dx_1 tau1)
                                    DomainOfIntegration::interior>(epsilon,
                                                                   secondcomponent,
                                                                   firstcomponent)
          ));

  auto minInnerProduct = make_InnerProduct(solutionSpaces,
          make_tuple(
              make_IntegralTerm<3,3,IntegrationType::valueValue,              // (u^,u^)
                                    DomainOfIntegration::interior>(1),
              make_IntegralTerm<3,3,IntegrationType::gradGrad,                // (beta grad u^,beta grad u^)
                                    DomainOfIntegration::interior>(1, beta),
              make_IntegralTerm<3,3,IntegrationType::gradGrad,                // epsilon (dx_1 u^,dx_1 u^)
                                    DomainOfIntegration::interior>(epsilon, firstcomponent),
              make_IntegralTerm<3,3,IntegrationType::gradGrad,                // epsilon (dx_2 u^,dx_2 u^)
                                    DomainOfIntegration::interior>(epsilon, secondcomponent)
          ));

  typedef decltype(bilinearForm) BilinearForm;
  typedef decltype(innerProduct) InnerProduct;
  typedef decltype(minInnerProduct) MinInnerProduct;
  typedef Functions::TestspaceCoefficientMatrix<BilinearForm, InnerProduct> TestspaceCoefficientMatrix;

  TestspaceCoefficientMatrix testspaceCoefficientMatrix(bilinearForm, innerProduct);

  typedef Functions::OptimalTestBasis<TestspaceCoefficientMatrix, 0> FEBasisOptimalTest0;              // v
  FEBasisOptimalTest0 feBasisTest0(testspaceCoefficientMatrix);
  typedef Functions::OptimalTestBasis<TestspaceCoefficientMatrix, 1> FEBasisOptimalTest1;              // tau1
  FEBasisOptimalTest1 feBasisTest1(testspaceCoefficientMatrix);
  typedef Functions::OptimalTestBasis<TestspaceCoefficientMatrix, 2> FEBasisOptimalTest2;              // tau2
  FEBasisOptimalTest2 feBasisTest2(testspaceCoefficientMatrix);


  auto optimalTestSpaces = make_tuple(feBasisTest0, feBasisTest1, feBasisTest2);

  auto systemAssembler
     = make_DPG_SystemAssembler(optimalTestSpaces, solutionSpaces,
                                bilinearForm);
//  auto systemAssembler
//     = make_Saddlepoint_SystemAssembler(testSpaces, solutionSpaces,
//                                        bilinearForm, innerProduct);
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

  auto rightHandSide
    = make_DPG_LinearForm(systemAssembler.getTestSpaces(),
    // make_Saddlepoint_LinearForm(
    //  systemAssembler.getTestSpaces(),
    //  systemAssembler.getSolutionSpaces(),
        std::make_tuple(
            make_LinearIntegralTerm<0>(fieldRHS)
          , make_LinearIntegralTerm<1>([] (const Domain& x) { return 0.;})
          , make_LinearIntegralTerm<2>([] (const Domain& x) { return 0.;})
          ));
  systemAssembler.assembleSystem(stiffnessMatrix, rhs, rightHandSide);

  MatrixType testMatrix(stiffnessMatrix);

  for (unsigned int i=0; i<testMatrix.N(); i++)
  {
    for (unsigned int j=0; j<testMatrix.M(); j++)
    {
      if (testMatrix.exists(i,j))
      {
        testMatrix[i][j]=0;
      }
    }
  }
  // Add minimization property for u^ on (near-)characteristic boundary if epsilon is closed to zero
  systemAssembler.applyMinimization<3, MinInnerProduct,2>
                    (testMatrix,
                     minInnerProduct,
                     beta,
                     delta,
                     sqrtepsilon);  //TODO is this really sqrtepsilon or something similar?

std::ofstream file("testmatrix.txt");
printmatrix(file , testMatrix, "testmatrix", "--");


  // Add minimization property for u^ on (near-)characteristic boundary if epsilon is closed to zero
  systemAssembler.applyMinimization<3, MinInnerProduct,2>
                    (stiffnessMatrix,
                     minInnerProduct,
                     beta,
                     delta,
                     sqrtepsilon);  //TODO is this really sqrtepsilon or something similar?


  // Set weak zero-boundary conditions for u^ (outflow boundary)
  systemAssembler.applyWeakBoundaryCondition<3, 2>
                    (stiffnessMatrix,
                     beta,
                     mu);

  // Determine Dirichlet dofs for u^ (inflow boundary) and set them to zero
  {
    std::vector<bool> dirichletNodesInflow;
    boundaryTreatmentInflow(std::get<3>(solutionSpaces),
                            dirichletNodesInflow);
    systemAssembler.applyDirichletBoundary<3,double>
        (stiffnessMatrix,
         rhs,
         dirichletNodesInflow,
         0.);
  }

  tassembled = time(0);

// computation of condition number of stiffnessMatrix (does not work)

//  const bool verbose = true;
//  const unsigned int arppp_a_verbosity_level = 2;
//  const unsigned int pia_verbosity_level = 1;
//  MatrixInfo<MatrixType> matrixInfo
//      (stiffnessMatrix,verbose,arppp_a_verbosity_level,pia_verbosity_level);

//  MatrixType::field_type condition = matrixInfo.getCond2(false);
//  std::cout << "Condition number in 2 Norm = " << condition <<std::endl;


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

//printmatrix(std::cout , stiffnessMatrix, "matrix", "--");

std::ofstream file1("matrix.txt");
printmatrix(file1 , stiffnessMatrix, "matrix", "--");

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

  size_t nu = std::get<0>(solutionSpaces).size();
  size_t nsigma1 = std::get<1>(solutionSpaces).size();
  size_t nsigma2 = std::get<2>(solutionSpaces).size();
  size_t nuhat = std::get<3>(solutionSpaces).size();
  size_t nsigmahat = std::get<4>(solutionSpaces).size();

  size_t ntestv = std::get<0>(testSpaces).size();
  size_t ntesttau1 = std::get<1>(testSpaces).size();
  size_t ntesttau2 = std::get<2>(testSpaces).size();
  unsigned int shift = 0;

  VectorType u(nu);
  VectorType sigma1(nsigma1);
  VectorType sigma2(nsigma2);
  VectorType uhat(nuhat);
  VectorType sigmahat(nsigmahat);

  u=0;
  sigma1=0;
  sigma2=0;
  u=0;
  sigmahat=0;

  for (unsigned int i=0; i<nu; i++)
  {
    u[i] = x[shift+i];
  }
  shift += nu;

  for (unsigned int i=0; i<nsigma1; i++)
  {
    sigma1[i] = x[shift+i];
  }
  shift += nsigma1;

  for (unsigned int i=0; i<nsigma2; i++)
  {
    sigma2[i] = x[shift+i];
  }
  shift += nsigma2;

  for (unsigned int i=0; i<nuhat; i++)
  {
    uhat[i] = x[shift+i];
  }
  shift += nuhat;

  for (unsigned int i=0; i<nsigmahat; i++)
  {
    sigmahat[i] = x[shift+i];
  }
  shift += nsigmahat;

  auto feBasisInterior = std::get<0>(solutionSpaces);
  auto feBasisSigma1 = std::get<1>(solutionSpaces);
  auto feBasisSigma2 = std::get<2>(solutionSpaces);
  auto feBasisTrace = std::get<3>(solutionSpaces);
  auto feBasisFace = std::get<4>(solutionSpaces);


  ////////////////////////////////////////////////////////////////////////////
  //  Error evaluation
  ////////////////////////////////////////////////////////////////////////////


  std::cout << std::endl << "******** Computation of errors *************" << std::endl;
  // The exact solution against which we are comparing our FEM solution
  auto uExact = std::make_tuple(fieldExact);
  // to retrive the value out of this:
  // std::get<0>(uExact)(x)

  //We build an object of type ErrorTools to study errors, residuals and do hp-adaptivity
  ErrorTools errorTools = ErrorTools();

  //We compute the L2 error between the exact and the fem solutions
  double err = errorTools.computeL2error<1>(feBasisInterior,u,uExact);
  std::cout << "'Exact' error u: || u - u_fem ||_L2 = " << err << std::endl;

  //// todo: h-refinement
  //errorTools->hRefinement(grid);
  //// todo: p-refinement

/*  // A posteriori error
      //We compute the rhs in the form given by the projection approach
  rhsAssembler.assembleRhs(rhs, rightHandSide);
      //It is necessary to provide rhs in the above form to call this aPosterioriError method
  double aposterioriErr = errorTools.aPosterioriError(bilinearForm,innerProduct,x,rhs);
  std::cout << "A posteriori error: || (u,trace u) - (u_fem,theta) || = " << aposterioriErr << std::endl; */


  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //  We need to subsample, because VTK cannot natively display real second-order functions
  //////////////////////////////////////////////////////////////////////////////////////////////

  auto uFunction
      = Dune::Functions::makeDiscreteGlobalBasisFunction<double>
            (feBasisInterior, Dune::TypeTree::hybridTreePath(), u);
  auto localUFunction = localFunction(uFunction);

  auto uhatFunction
      = Dune::Functions::makeDiscreteGlobalBasisFunction<double>
            (feBasisTrace, Dune::TypeTree::hybridTreePath(), uhat);
  auto localUhatFunction = localFunction(uhatFunction);

//  auto sigma1Function
//      = Dune::Functions::makeDiscreteGlobalBasisFunction<double>
//            (feBasisSigma1, Dune::TypeTree::hybridTreePath(), sigma1);
//  auto localSigma1Function = localFunction(sigma1Function);

//  auto sigma2Function
//      = Dune::Functions::makeDiscreteGlobalBasisFunction<double>
//            (feBasisSigma2, Dune::TypeTree::hybridTreePath(), sigma2);
//  auto localSigma2Function = localFunction(sigma2Function);

  auto sigmahatFunction
      = Dune::Functions::makeDiscreteGlobalBasisFunction<double>
            (feBasisFace, Dune::TypeTree::hybridTreePath(), sigmahat);
  auto localSigmahatFunction = localFunction(sigmahatFunction);

  SubsamplingVTKWriter<GridView> vtkWriter(gridView,2);
  vtkWriter.addVertexData(localUFunction, VTK::FieldInfo("u", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write("convdiff_cube_scaled_min_"+ std::to_string(nelements)+"_"+ std::to_string(epsinv));

//  SubsamplingVTKWriter<GridView> vtkWriter1(gridView,2);
//  vtkWriter1.addVertexData(localUhatFunction, VTK::FieldInfo("uhat", VTK::FieldInfo::Type::scalar, 1));
//  vtkWriter1.write("convdiff_cube_trace_"+ std::to_string(nelements)+"_"+ std::to_string(epsinv));

//  SubsamplingVTKWriter<GridView> vtkWriter2(gridView,2);
//  vtkWriter2.addVertexData(localSigma1Function, VTK::FieldInfo("sigma1", VTK::FieldInfo::Type::scalar, 1));
//  vtkWriter2.write("convdiff_cube_sigma1_"+ std::to_string(nelements));

//  SubsamplingVTKWriter<GridView> vtkWriter3(gridView,2);
//  vtkWriter3.addVertexData(localSigma2Function, VTK::FieldInfo("sigma2", VTK::FieldInfo::Type::scalar, 1));
//  vtkWriter3.write("convdiff_cube_sigma2_"+ std::to_string(nelements));

//  SubsamplingVTKWriter<GridView> vtkWriter4(gridView,2);
//  vtkWriter4.addVertexData(localSigmahatFunction, VTK::FieldInfo("sigmahat", VTK::FieldInfo::Type::scalar, 1));
//  vtkWriter4.write("convdiff_cube_sigmahat_"+ std::to_string(nelements)+"_"+ std::to_string(epsinv));

//#endif

# if 0
  VectorType sol(feBasisTest.size());
  for(unsigned int idx = 0; idx<feBasisTest.size(); idx++)
  {
    sol=0;
  //int idx = 1;
  //int idx =atoi(argv[1]);
    sol[idx]=1;

  //std::cout << u << std::endl;
    auto uFunction
        = Dune::Functions::makeDiscreteGlobalBasisFunction<double>
              (feBasisTest, Dune::TypeTree::hybridTreePath(), sol);
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
  VectorType sol(feBasisTrace.size());
  for(unsigned int idx = 0; idx<feBasisTrace.size(); idx++)
  {
    sol=0;
  //int idx = 1;
  //int idx =atoi(argv[1]);
    sol[idx]=1;

  //std::cout << u << std::endl;
    auto uhatbasisFunction
        = Dune::Functions::makeDiscreteGlobalBasisFunction<double>
              (feBasisTrace, Dune::TypeTree::hybridTreePath(), sol);
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
