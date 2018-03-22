#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>

#include <array>
#include <tuple>
#include <vector>

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

#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/saddlepoint_system_assembler.hh>
#include <dune/dpg/functionplotter.hh>
#include <dune/dpg/functions/gridviewfunctions.hh>



using namespace Dune;



int main()
{
  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  constexpr int dim = 2;
  using GridType = YaspGrid<dim>;
  const FieldVector<double,dim> l(1);
  const std::array<int,dim> elements = {5, 5};
  const GridType grid(l, elements);

  using GridView = GridType::LeafGridView;
  const GridView gridView = grid.leafGridView();
  const FieldVector<double,dim> beta = {1., 0.5};

  /////////////////////////////////////////////////////////
  //   Choose finite element spaces
  /////////////////////////////////////////////////////////

  typedef Functions::PQkTraceNodalBasis<GridView, 2> FEBasisTrace; // u^
  typedef Functions::LagrangeDGBasis<GridView, 1> FEBasisInterior; // u

  typedef Functions::LagrangeDGBasis<GridView, 4> FEBasisTest;     // v

  auto spaces = make_space_tuple<FEBasisTest,
                                 FEBasisTrace, FEBasisInterior>(gridView);

  auto testSpaces = make_space_tuple_view<0, 1>(spaces);
  auto solutionSpaces = make_space_tuple_view<1, 2>(spaces);

  /////////////////////////////////////////////////////////
  //   Choose a bilinear form
  /////////////////////////////////////////////////////////

  auto cFunc = Functions::makeConstantGridViewFunction(0., gridView);
  auto betaFunc = Functions::makeConstantGridViewFunction(beta, gridView);
  auto oneFunc = Functions::makeConstantGridViewFunction(1., gridView);
  auto minusOneFunc = Functions::makeConstantGridViewFunction(-1., gridView);


  auto bilinearForm = make_BilinearForm(testSpaces, solutionSpaces,
          make_tuple(
              make_IntegralTerm<0,1,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(cFunc),
              make_IntegralTerm<0,1,IntegrationType::gradValue,
                                    DomainOfIntegration::interior>
                                (minusOneFunc, betaFunc),
              make_IntegralTerm<0,0,IntegrationType::normalVector,
                                    DomainOfIntegration::face>
                                (oneFunc, betaFunc)));
  auto innerProduct = make_InnerProduct(testSpaces,
          make_tuple(
              make_IntegralTerm<0,0,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(oneFunc),
              make_IntegralTerm<0,0,IntegrationType::gradGrad,
                                    DomainOfIntegration::interior>
                                (oneFunc, betaFunc)));
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

  auto rightHandSide
    = make_LinearForm(
        spaces,
        std::make_tuple(
            make_LinearIntegralTerm<0,
                                    LinearIntegrationType::valueFunction,
                                    DomainOfIntegration::interior>
                                   (oneFunc)
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
  const size_t nTest = std::get<0>(*testSpaces).size();
  const size_t nFace = std::get<0>(*solutionSpaces).size();

  FunctionPlotter uPlotter("solution_transport");
  uPlotter.plot("u", x, std::get<1>(*solutionSpaces), 2, nTest+nFace);

  return 0;
}
