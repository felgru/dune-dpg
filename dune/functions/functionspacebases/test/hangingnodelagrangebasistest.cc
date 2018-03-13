// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/dpg/functions/localindexsetiteration.hh>
#include <dune/functions/functionspacebases/hangingnodep2nodalbasis.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include "hangingnodebasistest.hh"

using namespace Dune;

int main() try
{
  constexpr int dim = 2;
  using HostGrid = UGGrid<dim>;
  const FieldVector<double,dim> lower = {0, 0};
  const FieldVector<double,dim> upper = {1, 1};
  const std::array<unsigned int,dim> numElements = {1, 1};

  std::shared_ptr<HostGrid> hostGrid = StructuredGridFactory<HostGrid>
                                ::createSimplexGrid(lower, upper, numElements);
  hostGrid->setClosureType(HostGrid::NONE);

  // We use a SubGrid as it will automatically make sure that we do
  // not have more than difference 1 in the levels of neighboring
  // elements. This is necessary since HangingNodeP2NodalBasis does
  // not implement higher order hanging nodes constraints.

  const auto grid = createHangingNodeRefinedSubGrid(*hostGrid);
  using Grid = decltype(grid)::element_type;
  using GridView = Grid::LeafGridView;
  const GridView gridView = grid->leafGridView();

  bool success = true;

  std::cout << "Testing HangingNodeP2NodalBasis on 2d"
            << " triangular grid." << std::endl;

  Functions::HangingNodeP2NodalBasis<GridView> hangingNodeBasis(gridView);

  success &= globalIndicesFormConsecutiveRange(hangingNodeBasis);
  success &= constraintsFulfillContinuityEquation(hangingNodeBasis);

  return success ? 0 : 1;
}
catch (Dune::Exception e)
{
  std::cout << e << std::endl;
  return 1;
}
