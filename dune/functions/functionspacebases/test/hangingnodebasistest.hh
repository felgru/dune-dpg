// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_HANGINGNODEBASISTEST_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_HANGINGNODEBASISTEST_HH

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/dpg/functions/localindexsetiteration.hh>
#include <dune/geometry/quadraturerules.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <dune/subgrid/subgrid.hh>
#pragma GCC diagnostic pop

namespace Dune {

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

template<class GlobalBasis>
bool constraintsFulfillContinuityEquation(const GlobalBasis& feBasis)
{
  auto localView = feBasis.localView();
  auto localIndexSet = feBasis.localIndexSet();
  auto dominatedElementLocalView = feBasis.localView();
  auto dominatedElementLocalIndexSet = feBasis.localIndexSet();
  auto gridView = feBasis.gridView();

  bool success = true;
  constexpr double eps = 1e-10;

  for (const auto& element : elements(gridView))
  {
    localView.bind(element);
    localIndexSet.bind(localView);
    for(auto&& intersection : intersections(gridView, element))
    {
      if(!intersection.conforming())
        if(intersection.inside().level() < intersection.outside().level())
          // inside dominates outside (with one level difference)
        {
          dominatedElementLocalView.bind(intersection.outside());
          dominatedElementLocalIndexSet.bind(dominatedElementLocalView);

          const auto geometryInDominatingElement
              = intersection.geometryInInside();
          const auto geometryInDominatedElement
              = intersection.geometryInOutside();
          const auto& quad // TODO: replace 3 with degree of basis
              = QuadratureRules<double, 1>::rule(GeometryTypes::line, 3);
          for(size_t pt=0; pt < quad.size(); pt++)
          {
            // point-wise check of continuity condition
            const FieldVector<double,1>& quadPos = quad[pt].position();
            const auto quadPosDominated
                = geometryInDominatedElement.global(quadPos);
            const auto quadPosDominating
                = geometryInDominatingElement.global(quadPos);

            std::vector<FieldVector<double,1> > valuesDominating;
            localView.tree().finiteElement().localBasis()
              .evaluateFunction(quadPosDominating, valuesDominating);
            std::vector<FieldVector<double,1> > valuesDominated;
            dominatedElementLocalView.tree().finiteElement().localBasis()
              .evaluateFunction(quadPosDominated, valuesDominated);
            const auto& localCoefficientsDominating
              = localView.tree().finiteElement().localCoefficients();
            iterateOverLocalIndexSet(
                localIndexSet,
                [&](size_t j, auto gj)
                {
                  if(std::fabs(valuesDominating[j]) > eps) {
                    double valueDominated = 0.;
                    std::vector<std::pair<double, double>> constraints;
                    iterateOverLocalIndexSet(
                        dominatedElementLocalIndexSet,
                        [&](size_t i, auto gi)
                        {
                          if(gi == gj) {
                            valueDominated += valuesDominated[i];
                            constraints.emplace_back(1., valuesDominated[i]);
                          }
                        },
                        [](size_t /* i */) {},
                        [&](size_t i, auto gi, double wi)
                        {
                          if(gi == gj) {
                            valueDominated += wi*valuesDominated[i];
                            constraints.emplace_back(wi, valuesDominated[i]);
                          }
                        });
                    if(std::fabs(valuesDominating[j] - valueDominated) > eps)
                    {
                      std::cout << "Dominating FE at local index " << j
                        << " and position " << quadPosDominating
                        << " (" << quadPos << " locally)"
                        << " evaluates to " << valuesDominating[j]
                        << " but dominated basis evaluations at "
                        << quadPosDominated << " have weighted sum "
                        << valueDominated << "!\n";
                      std::cout << "It is a DoF of codim "
                        << localCoefficientsDominating.localKey(j).codim()
                        << '\n';
                      std::cout << "weights * dominated basis evaluations:\n";
                      for(const auto& c: constraints) {
                        std::cout << c.first << " * " << c.second << '\n';
                      }
                      success = false;
                    }
                  }
                },
                [](size_t /* j */) {},
                [&](size_t j, auto /* gj */, double /* wj */)
                {
                  if(std::fabs(valuesDominating[j]) > eps) {
                    DUNE_THROW(Dune::InvalidStateException,
                      "On this edge, the indices of the dominating"
                      " element should not be constrained.");
                  }
                });
          }
        }
    }
  }

  return success;
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

}

#endif // defined(DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_HANGINGNODEBASISTEST_HH)
