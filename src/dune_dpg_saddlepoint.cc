#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>

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
#include <dune/functions/functionspacebases/pqktransportnodalbasis.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discretescalarglobalbasisfunction.hh>

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
      isBoundary = isBoundary || lagrangeNodes[i][j] < 1e-8;

    if (isBoundary)

      dirichletNodes[i] = true;
  }
}




int main(int argc, char** argv)
{
  try{
  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {5, 5};
  GridType grid(l,elements);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid.leafGridView();
  FieldVector<double,dim> beta={1.,0.5};

  /////////////////////////////////////////////////////////
  //   Choose finite element spaces
  /////////////////////////////////////////////////////////

  typedef Functions::PQkTraceNodalBasis<GridView, 2> FEBasisTrace; // u^
  typedef Functions::LagrangeDGBasis<GridView, 1> FEBasisInterior; // u
  auto solutionSpaces = std::make_tuple(FEBasisTrace(gridView),
                                        FEBasisInterior(gridView));

  //typedef Functions::LagrangeDGBasis<GridView, 4> FEBasisTest;     // v
  //auto testSpaces = std::make_tuple(FEBasisTest(gridView));

  typedef Functions::PQkTransportBasis<GridView,4> FEBasisTest;              // v
  auto testSpaces = std::make_tuple(FEBasisTest(gridView,beta));

  /////////////////////////////////////////////////////////
  //   Choose a bilinear form
  /////////////////////////////////////////////////////////


  auto bilinearForm = make_BilinearForm(testSpaces, solutionSpaces,
          make_tuple(
              make_IntegralTerm<0,1,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(0.),
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
  auto systemAssembler = make_SystemAssembler(testSpaces, solutionSpaces,
          bilinearForm, innerProduct, SaddlepointFormulation());

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

  auto rightHandSide = std::make_tuple([] (const Domain& x) { return 1.;});
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
    boundaryTreatmentInflow(std::get<0>(solutionSpaces),
                            dirichletNodesInflow);
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
    systemAssembler.applyDirichletBoundaryTest<0>
        (stiffnessMatrix,
         rhs,
         dirichletNodesInflowTest,
         0.);
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
  u=0;
  for (size_t i=0; i<nInner; i++)
  {
    /* TODO: use a precomputed solution space offset. */
    u[i] = x[nTest+nFace+i];
  }

  auto innerSpace = std::get<1>(solutionSpaces);
  Dune::Functions::DiscreteScalarGlobalBasisFunction
      <typename std::remove_reference<decltype(innerSpace)>::type, decltype(u)>
          uFunction(innerSpace,u);
  auto localUFunction = localFunction(uFunction);

  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //  We need to subsample, because VTK cannot natively display real second-order functions
  //////////////////////////////////////////////////////////////////////////////////////////////
  SubsamplingVTKWriter<GridView> vtkWriter(gridView,2);
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
