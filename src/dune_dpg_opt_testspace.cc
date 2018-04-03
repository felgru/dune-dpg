#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <cstdlib> // for std::exit()

#include <array>
#include <memory>
#include <tuple>
#include <vector>

#include <boost/math/constants/constants.hpp>

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

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pqktracenodalbasis.hh>
#include <dune/functions/functionspacebases/pqkfacenodalbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/pqkdgrefineddgnodalbasis.hh>

#include <dune/dpg/dpg_system_assembler.hh>
#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/functionplotter.hh>

#include <chrono>

using namespace Dune;



int main(int argc, char** argv)
{
  if(argc != 2) {
    std::cerr << "Usage: " << argv[0] << " n" << std::endl << std::endl
              << "Solves the transport problem on an nxn grid." << std::endl;
    std::exit(1);
  }
  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  constexpr int dim = 2;
  typedef UGGrid<dim> GridType;

  const unsigned int nelements = atoi(argv[1]);

  const FieldVector<double,dim> lower = {0, 0};
  const FieldVector<double,dim> upper = {1, 1};
  const std::array<unsigned int,dim> elements = {nelements, nelements};

  //std::shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createCubeGrid(lower, upper, elements);

  std::shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createSimplexGrid(lower, upper, elements);

  //std::shared_ptr<GridType> grid = std::shared_ptr<GridType>(GmshReader<GridType>::read("irregular-square.msh"));

  typedef GridType::LeafGridView GridView;
  const GridView gridView = grid->leafGridView();

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  // u
  typedef Functions::LagrangeDGBasis<GridView, 1> FEBasisInterior;
  FEBasisInterior feBasisInterior(gridView);

  // bulk term corresponding to u^
  typedef Functions::PQkNodalBasis<GridView, 2> FEBasisTrace;
  FEBasisTrace feBasisTrace(gridView);

  auto solutionSpaces
    = make_space_tuple<FEBasisInterior, FEBasisTrace>(gridView);

  // v search space
  typedef Functions::PQkDGRefinedDGBasis<GridView, 1, 3> FEBasisTest;
  //typedef Functions::LagrangeDGBasis<GridView, 3> FEBasisTest;
  auto testSpaces = make_space_tuple<FEBasisTest>(gridView);

  const FieldVector<double, dim> beta
             = {cos(boost::math::constants::pi<double>()/8),
                sin(boost::math::constants::pi<double>()/8)};
  const double c = 2;

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
              make_IntegralTerm<0,0,IntegrationType::gradGrad,
                                    DomainOfIntegration::interior>(1., beta),
              make_IntegralTerm<0,0,IntegrationType::travelDistanceWeighted,
                                    DomainOfIntegration::face>(1., beta)));

  /////////////////////////////////////////////////////////
  //   Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////

  typedef BlockVector<FieldVector<double,1> > VectorType;
  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

  VectorType rhs;
  MatrixType stiffnessMatrix;

  typedef decltype(std::declval<typename GridView::template Codim<0>::Entity>().geometry()) Geometry;
  GeometryBuffer<Geometry> geometryBuffer;

  auto systemAssembler
     = make_DPGSystemAssembler(bilinearForm, innerProduct, geometryBuffer);

  /////////////////////////////////////////////////////////
  //  Assemble the system
  /////////////////////////////////////////////////////////
  using Domain = GridType::template Codim<0>::Geometry::GlobalCoordinate;

  auto rightHandSide
    = make_LinearForm(testSpaces,
                      std::make_tuple(make_LinearIntegralTerm<0,
                                            LinearIntegrationType::valueFunction,
                                            DomainOfIntegration::interior>(
                                 [] (const Domain& x) { return 1.;})));

  systemAssembler.assembleSystem(stiffnessMatrix, rhs, rightHandSide);

  /////////////////////////////////////////////////
  //   Choose an initial iterate
  /////////////////////////////////////////////////
  VectorType x(feBasisTrace.size()
               +feBasisInterior.size());
  x = 0;

#if 0
  const double delta = 1e-8;
  systemAssembler.defineCharacteristicFaces<1>
                    (stiffnessMatrix,
                     rhs,
                     beta,
                     delta);
#endif
  // Determine Dirichlet dofs for u^ (inflow boundary)
  {
    std::vector<bool> dirichletNodesInflow;
    BoundaryTools::getInflowBoundaryMask(std::get<1>(*solutionSpaces),
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

  UMFPack<MatrixType> umfPack(stiffnessMatrix, 2);
  InverseOperatorResult statistics;
  umfPack.apply(x, rhs, statistics);

  //////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //////////////////////////////////////////////////////////////////
  FunctionPlotter uPlotter("transport_simplex_"
                          + std::to_string(nelements) + "_"
                          + std::to_string(beta[0]) + "_"
                          + std::to_string(beta[1]));
  FunctionPlotter uhatPlotter("transport_simplex_trace_"
                             + std::to_string(nelements) + "_"
                             + std::to_string(beta[0]) + "_"
                             + std::to_string(beta[1]));
  uPlotter.plot("u", x, feBasisInterior, 2, 0);
  uhatPlotter.plot("uhat", x, feBasisTrace, 2,
                   feBasisInterior.size());

  return 0;
}
