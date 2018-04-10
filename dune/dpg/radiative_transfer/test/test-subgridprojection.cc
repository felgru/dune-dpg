// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <memory>

#include "../subgridprojection.hh"
#include <dune/functions/functionspacebases/bernsteindgbasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <dune/subgrid/subgrid.hh>
#pragma GCC diagnostic pop

using namespace Dune;

template<class Function1, class Function2>
bool functions_are_identical
      (const Function1& function1, const Function2& function2,
       double tolerance = 1e-8)
{
  auto localFunction1 = localFunction(function1);
  auto localFunction2 = localFunction(function2);

  const auto entitySet = function1.entitySet();
  constexpr int dim = decltype(entitySet)::Element::dimension;
  bool identical = true;
  for(const auto& e : entitySet) {
    localFunction1.bind(e);
    localFunction2.bind(e);
    const unsigned int quadratureOrder = 2;
    const Dune::QuadratureRule<double, dim>& quad =
          Dune::QuadratureRules<double, dim>::rule(e.type(),
                                                   quadratureOrder);
    for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {
      const FieldVector<double,dim>& quadPos = quad[pt].position();
      const auto difference
        = std::fabs(localFunction1(quadPos)-localFunction2(quadPos));
      if(difference > tolerance)
      {
        std::cout << "On element with center " << e.geometry().center()
          << ", function evaluations differ by "
          << difference << '\n';
        identical = false;
      }
    }
  }
  return identical;
}


template<class FEVector, class FEBasis, class Function>
bool finite_element_data_matches_function
      (const FEVector& feData, const FEBasis& feBasis,
       const Function& function, double tolerance = 1e-8)
{
  auto feFunction
      = Dune::Functions::makeDiscreteGlobalBasisFunction<double>
            (feBasis, Dune::TypeTree::hybridTreePath(), feData);
  auto analyticFunction = Functions::makeAnalyticGridViewFunction(
                                function, feBasis.gridView());
  return functions_are_identical(feFunction, analyticFunction, tolerance);
}

template<class HostGrid>
std::unique_ptr<SubGrid<HostGrid::dimension, HostGrid, false>>
createSubGrid(HostGrid& hostGrid)
{
  using SubGrid = SubGrid<HostGrid::dimension, HostGrid, false>;
  auto subGrid = std::make_unique<SubGrid>(hostGrid);
  subGrid->createBegin();
  subGrid->insertLevel(0);
  subGrid->createEnd();
  subGrid->setMaxLevelDifference(1);

  return subGrid;
}

template<class Grid, class Coordinate>
void refineNear(Grid& grid, const Coordinate& coordinate)
{
  auto leafGridView = grid.leafGridView();
  for(const auto& e : elements(leafGridView)) {
    Coordinate eCenter = e.geometry().center();
    if((eCenter - coordinate).two_norm() < 1e-5) {
      grid.mark(1, e);
    }
  }
  grid.preAdapt();
  grid.adapt();
  grid.postAdapt();
}

int main() {
  constexpr int dim = 2;
  using HostGrid = UGGrid<dim>;
  const FieldVector<double,dim> lower = {0, 0};
  const FieldVector<double,dim> upper = {1, 1};
  const std::array<unsigned int,dim> numElements = {1, 1};

  std::unique_ptr<HostGrid> hostGrid = StructuredGridFactory<HostGrid>
                                ::createSimplexGrid(lower, upper, numElements);
  hostGrid->setClosureType(HostGrid::NONE);
  refineNear(*hostGrid, FieldVector<double, dim>{1./3., 2./3.});

  const auto subGrid = createSubGrid(*hostGrid);
  refineNear(*subGrid, FieldVector<double, dim>{1./3., 2./3.});
  refineNear(*subGrid, FieldVector<double, dim>{2./3., 1./3.});

  using SubGrid = decltype(subGrid)::element_type;
  using SubGridView = SubGrid::LeafGridView;
  using HostGridView = typename HostGrid::LeafGridView;
  using FEBasisSubGrid = Functions::BernsteinDGBasis<SubGridView, 2>;
  using FEBasisHostGrid = changeGridView_t<FEBasisSubGrid, HostGridView>;
  using FEBasisProjecteeSubGrid = Functions::BernsteinDGBasis<SubGridView, 3>;

  using FEVector = BlockVector<FieldVector<double,1> >;

  std::cout << "Testing SubGrid projection on 2d triangular grid.\n";

  bool success = true;

  auto testFunction =
    [](const FieldVector<double, dim>& pos)
      { return pos[0] * pos[1] + pos[0] * pos[0]; };

  FEBasisSubGrid feBasisSubGrid(subGrid->leafGridView());

  FEVector x(feBasisSubGrid.size());
  Functions::interpolate(feBasisSubGrid, x, testFunction);

  FEBasisHostGrid feBasisHostGrid(hostGrid->leafGridView());

  FEVector xHost(feBasisHostGrid.size());

  // interpolateFromSubGrid(feBasisSubGrid, x, feBasisHostGrid, xHost);
  Functions::interpolate(feBasisHostGrid, xHost, testFunction);
  success = success && finite_element_data_matches_function(xHost,
                                        feBasisHostGrid, testFunction);

  std::cout << "Testing projection without refinement.\n";

  FEBasisProjecteeSubGrid feBasisProjecteeSubGrid(subGrid->leafGridView());
  auto subGridData = attachDataToSubGrid(
            feBasisProjecteeSubGrid,
            feBasisHostGrid,
            xHost);

  auto restoredX
    = subGridData.restoreDataToRefinedSubGrid(feBasisProjecteeSubGrid);
  success = success && finite_element_data_matches_function(restoredX,
                                        feBasisProjecteeSubGrid, testFunction);

  std::cout << "Testing projection with refinement.\n";

  refineNear(*subGrid, FieldVector<double, dim>{1./3., 2./3.});
  refineNear(*subGrid, FieldVector<double, dim>{2./3., 1./3.});
  feBasisSubGrid.update(subGrid->leafGridView());
  feBasisProjecteeSubGrid.update(subGrid->leafGridView());

  restoredX
    = subGridData.restoreDataToRefinedSubGrid(feBasisProjecteeSubGrid);
  success = success && finite_element_data_matches_function(restoredX,
                                        feBasisProjecteeSubGrid, testFunction);

  return success ? 0 : 1;
}
