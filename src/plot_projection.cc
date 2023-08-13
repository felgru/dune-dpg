#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <cstdlib> // for std::exit()

#include <array>
#include <cmath>
#include <vector>

#include <boost/math/constants/constants.hpp>

#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dune/dpg/radiative_transfer/projection.hh>
#include <dune/dpg/functionplotter.hh>

#include <dune/subgrid/subgrid.hh>

using namespace Dune;

int main(int argc, char** argv)
{
  if(argc != 2) {
    std::cerr << "Usage: " << argv[0] << " n" << std::endl << std::endl
              << "Plot projection on an nxn grid.\n";
    std::exit(1);
  }
  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  constexpr int dim = 2;
  using HostGrid = UGGrid<dim>;
  using Grid = SubGrid<dim, HostGrid, false>;

  unsigned int nelements = atoi(argv[1]);

  FieldVector<double,dim> lower = {0,0};
  FieldVector<double,dim> upper = {1,1};
  std::array<unsigned int,dim> elements = {nelements,nelements};

  std::unique_ptr<HostGrid> hostGrid = StructuredGridFactory<HostGrid>
                                  ::createSimplexGrid(lower, upper, elements);
  hostGrid->setClosureType(HostGrid::NONE);

  // We use a SubGrid as it will automatically make sure that we do
  // not have more than difference 1 in the levels of neighboring
  // elements. This is necessary since HangingNodeBernsteinP2Basis does
  // not implement higher order hanging nodes constraints.
  std::unique_ptr<Grid> grid = std::make_unique<Grid>(*hostGrid);
  {
    grid->createBegin();
    grid->insertLevel(hostGrid->maxLevel());
    grid->createEnd();
    grid->setMaxLevelDifference(1);
  }

  auto gridView = grid->leafGridView();

  /////////////////////////////////////////////////////////
  //   Choose finite element space
  /////////////////////////////////////////////////////////

  using GridView = Grid::LeafGridView;
  using FEBasis = Functions::LagrangeDGBasis<GridView, 1>;
  FEBasis feBasis(gridView);

  using Direction = FieldVector<double, 2>;

  auto targetFunction =
      [](const Direction& x) {
        return std::sin(x[0] * 2 * boost::math::constants::pi<double>());
      };

  const std::vector<double> projection = adaptivelyProjectFunctionToGrid(
      *grid,
      feBasis,
      targetFunction,
      0.4,
      0.00001);

  //////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //////////////////////////////////////////////////////////////////
  FunctionPlotter plotter("projection_" + std::to_string(nelements));
  plotter.plot("projection", projection, feBasis, 0);


  return 0;
}
