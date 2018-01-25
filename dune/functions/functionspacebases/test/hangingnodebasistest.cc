// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/functions/functionspacebases/hangingnodep2nodalbasis.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <dune/subgrid/subgrid.hh>
#pragma GCC diagnostic pop

using namespace Dune;

template<class GlobalBasis>
bool globalIndicesFormConsecutiveRange(const GlobalBasis& feBasis)
{
  std::vector<bool> seen(feBasis.size(), false);

  auto localView = feBasis.localView();
  auto localIndexSet = feBasis.localIndexSet();
  auto gridView = feBasis.gridView();

  for (const auto& element : elements(gridView))
  {
    localView.bind(element);
    localIndexSet.bind(localView);

    const auto& indicesLocalGlobal = localIndexSet.indicesLocalGlobal();
    for (const auto index : indicesLocalGlobal)
    {
      if (index[0] >= seen.size()) {
        std::cout << "Global index " << index
                  << " is larger than allowed.\n";
        return false;
      }

      seen[index[0]] = true;
    }
  }

  for (size_t i=0; i<seen.size(); i++)
    if (! seen[i]) {
      std::cout << "Global index " << i << " does not exist in index set.\n";
      return false;
    }

  return true;
}

template<class HostGrid>
std::unique_ptr<SubGrid<HostGrid::dimension, HostGrid, false>>
createHangingNodeRefinedSubGrid(HostGrid& hostGrid)
{
  constexpr int dim = HostGrid::dimension;
  using Grid = SubGrid<dim, HostGrid, false>;
  std::unique_ptr<Grid> grid = std::make_unique<Grid>(hostGrid);
  {
    grid->createBegin();
    grid->insertLevel(hostGrid.maxLevel());
    grid->createEnd();
    grid->setMaxLevelDifference(1);
  }
  // Refine the first cell in grid to obtain a hanging node.
  const auto gridView = grid->leafGridView();
  grid->mark(1, *gridView.template begin<0>());

  grid->preAdapt();
  grid->adapt();
  grid->postAdapt();

  return grid;
}

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

  std::cout << "Testing HangingNodeP2DNodalBasis on 2d"
            << " triangular grid." << std::endl;

  Functions::HangingNodeP2NodalBasis<GridView> hangingNodeBasis(gridView);

  success &= globalIndicesFormConsecutiveRange(hangingNodeBasis);

  return success ? 0 : 1;
}
catch (Dune::Exception e)
{
  std::cout << e << std::endl;
  return 1;
}
