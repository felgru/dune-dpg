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

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pqktracenodalbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/saddlepoint_system_assembler.hh>
#include <dune/dpg/functionplotter.hh>



using namespace Dune;



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
  auto solutionSpaces
    = make_space_tuple<FEBasisTrace, FEBasisInterior>(gridView);

  typedef Functions::LagrangeDGBasis<GridView, 4> FEBasisTest;     // v
  auto testSpaces = make_space_tuple<FEBasisTest>(gridView);

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
  auto systemAssembler
     = make_SaddlepointSystemAssembler(bilinearForm, innerProduct);

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

  auto rightHandSide
    = make_Saddlepoint_LinearForm(
        systemAssembler.getTestSpaces(),
        systemAssembler.getSolutionSpaces(),
        std::make_tuple(
            make_LinearIntegralTerm<0,
                                    LinearIntegrationType::valueFunction,
                                    DomainOfIntegration::interior>
                                   ([] (const Domain& x) { return 1.;})
          ));
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
    BoundaryTools::getInflowBoundaryMask(std::get<0>(*solutionSpaces),
                                         dirichletNodesInflow,
                                         beta);
    systemAssembler.applyDirichletBoundary<0>
        (stiffnessMatrix,
         rhs,
         dirichletNodesInflow,
         0.);
  }

#if 0
  // Determine Dirichlet dofs for v (inflow boundary)
  {
    std::vector<bool> dirichletNodesInflowTest;
    BoundaryTools::getInflowBoundaryMask(std::get<0>(*testSpaces),
                                         dirichletNodesInflowTest,
                                         beta);
    /* TODO: applyDirichletBoundaryTest has been removed */
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

  //////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //////////////////////////////////////////////////////////////////
  size_t nTest = std::get<0>(*testSpaces).size();
  size_t nFace = std::get<0>(*solutionSpaces).size();

  FunctionPlotter uPlotter("solution_transport");
  uPlotter.plot("u", x, std::get<1>(*solutionSpaces), 2, nTest+nFace);

  return 0;
  }
  catch (Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
