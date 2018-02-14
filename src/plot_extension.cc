#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <cstdlib> // for std::exit()

#include <array>
#include <tuple>
#include <vector>

#include <boost/math/constants/constants.hpp>

#include <dune/common/exceptions.hh> // We use exceptions

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

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/pqkdgrefineddgnodalbasis.hh>

#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/dpg_system_assembler.hh>
#include <dune/dpg/errortools.hh>
#include <dune/dpg/radiative_transfer/boundary_extension.hh>
#include <dune/dpg/rhs_assembler.hh>
#include <dune/dpg/functionplotter.hh>

#include <chrono>

using namespace Dune;

int main(int argc, char** argv)
{
  try{
  if(argc != 2) {
    std::cerr << "Usage: " << argv[0] << " n" << std::endl << std::endl
              << "Plot harmonic extension of boundary values on an nxn grid.\n";
    std::exit(1);
  }
  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  constexpr int dim = 2;
  typedef UGGrid<dim> GridType;

  unsigned int nelements = atoi(argv[1]);

  FieldVector<double,dim> lower = {0,0};
  FieldVector<double,dim> upper = {1,1};
  std::array<unsigned int,dim> elements = {nelements,nelements};

  std::shared_ptr<GridType> grid
      = StructuredGridFactory<GridType>::createSimplexGrid(lower, upper,
                                                           elements);

  auto gridView = grid->leafGridView();

  /////////////////////////////////////////////////////////
  //   Choose finite element space
  /////////////////////////////////////////////////////////

  using GridView = GridType::LeafGridView;
  using FEBasis = Functions::PQkNodalBasis<GridView, 2>;
  FEBasis feBasis(gridView);

  using Direction = FieldVector<double, 2>;
  using VectorType = BlockVector<FieldVector<double, 1>>;

  auto dirichletValueFunction =
      [](const Direction& x) {
        if(x[0] < .1) {
          const double xProj = x[1];
          if(xProj>=.5+.125 || xProj<=.5-.125) {
            return 0.;
          } else if(xProj<=.5) {
            return 8*(xProj-(.5-.125));
          } else {
            return 1-8*(xProj-.5);
          }
        } else {
          return 0.;
        }
      };

  VectorType extension = harmonic_extension_of_boundary_values(
      dirichletValueFunction,
      feBasis);

  //////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //////////////////////////////////////////////////////////////////
  FunctionPlotter plotter("boundary_extension_"
                          + std::to_string(nelements));
  plotter.plot("extension", extension, feBasis, 2);


  return 0;
  }
  catch (Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
