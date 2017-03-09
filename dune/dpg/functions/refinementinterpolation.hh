// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_FUNCTIONS_REFINEMENTINTERPOLATION_HH
#define DUNE_DPG_FUNCTIONS_REFINEMENTINTERPOLATION_HH

#include <utility>
#include <vector>

#include <dune/dpg/functions/localindexsetiteration.hh>

namespace Dune {

template<class GlobalBasis, class Vector>
std::pair<std::vector<typename GlobalBasis::LocalView::Element::EntitySeed>,
          std::vector<FieldVector<double, 1>>>
attachDataToGrid(
    const GlobalBasis& globalBasis,
    const Vector& data)
{
  using EntitySeed = typename GlobalBasis::LocalView::Element::EntitySeed;

  auto gridView = globalBasis.gridView();
  auto localView = globalBasis.localView();
  auto localIndexSet = globalBasis.localIndexSet();

  size_t maxEntityDoFs = localView.maxSize();

  std::vector<EntitySeed> entitySeeds(gridView.size(0));
  std::vector<typename Vector::value_type>
      entityData(gridView.size(0) * maxEntityDoFs);

  auto localEntitySeed = entitySeeds.begin();
  auto localData = entityData.begin();
  for(const auto& e : elements(gridView))
  {
    localView.bind(e);
    localIndexSet.bind(localView);

    *(localEntitySeed++) = e.seed();

    iterateOverLocalIndexSet(
      localIndexSet,
      [&](size_t i, auto gi)
      {
        localData[i] = data[gi[0]];
      },
      [&](size_t i){ localData[i] = 0; },
      [&](size_t i, auto gi, double wi)
      {
        localData[i] += wi * data[gi[0]];
      }
    );
    localData += maxEntityDoFs;
  }
  assert(localData == entityData.end());
  assert(localEntitySeed == entitySeeds.end());

  return std::make_pair(entitySeeds, entityData);
}

namespace detail {
  template<class GlobalBasis, class Geometry, class Iterator>
  struct RestoreDataToRefinedGridFunction {
    using FiniteElement = typename GlobalBasis::LocalView::Tree::FiniteElement;
    using Domain = typename Geometry::LocalCoordinate;
    using Range = FieldVector<double, 1>;

    RestoreDataToRefinedGridFunction(
        const FiniteElement& finiteElement,
        const Geometry& childInElementGeometry,
        const Iterator& elementData)
      : finiteElement(finiteElement),
        childInElementGeometry(childInElementGeometry),
        elementData(elementData)
    {}

    void evaluate(const Domain& x, Range& y) const {
      const Domain xElement = childInElementGeometry.global(x);

      y = Range(0);

      auto&& localBasis = finiteElement.localBasis();

      shapeFunctionValues.resize(localBasis.size());
      localBasis.evaluateFunction(xElement, shapeFunctionValues);

      for(size_t i = 0; i < localBasis.size(); i++) {
        y += elementData[i] * shapeFunctionValues[i];
      }
    }

    const FiniteElement& finiteElement;
    const Geometry childInElementGeometry;
    const Iterator& elementData;
    mutable std::vector<Range> shapeFunctionValues;
  };
}

/**
 * interpolate data saved with attachDataToGrid to refined grid
 *
 * \note This function assumes that the grid only includes one entity type.
 * \note Additionally it assumes that the finite elements before and after
 *       refinement are the same on all entities.
 */
template<class GlobalBasis>
std::vector<FieldVector<double, 1>>
restoreDataToRefinedGrid(
    const GlobalBasis& globalBasis,
    const std::pair<std::vector<typename GlobalBasis::LocalView::
                                              Element::EntitySeed>,
                    std::vector<FieldVector<double, 1>>>& storedData)
{
  using Element = typename GlobalBasis::LocalView::Element;

  auto gridView = globalBasis.gridView();
  auto& grid = gridView.grid();
  auto localView = globalBasis.localView();
  auto localIndexSet = globalBasis.localIndexSet();
  auto node = globalBasis.nodeFactory().node(Dune::TypeTree::hybridTreePath());

  size_t maxEntityDoFs = localView.maxSize();

  std::vector<FieldVector<double, 1>> data(globalBasis.size());

  for(size_t eIdx = 0, eMax = storedData.first.size(); eIdx < eMax; eIdx++)
  {
    auto e = grid.entity(storedData.first[eIdx]);
    localView.bind(e);
    localIndexSet.bind(localView);

    auto localData = storedData.second.begin() + eIdx * maxEntityDoFs;
    if(e.isLeaf()) {
      // directly copy cell data
      iterateOverLocalIndexSet(
        localIndexSet,
        [&](size_t i, auto gi)
        {
          data[gi[0]] = localData[i];
        },
        [](size_t i){},
        [](size_t i, auto gi, double wi)
        { /* will be copied from neighboring element */ }
      );
    } else {
      std::vector<FieldVector<double, 1>> childLocalData;
      for (const auto& child : descendantElements(e, grid.maxLevel()))
      {
        node.bind(child);
        localView.bind(child);
        localIndexSet.bind(localView);

        // This assumes that e and child share the same finite element
        // and thus the same entity type.
        auto&& localFiniteElement = node.finiteElement();
        assert(child.father() == e);
        auto oldGridFunction = detail::RestoreDataToRefinedGridFunction
          <GlobalBasis, typename Element::LocalGeometry, decltype(localData)>(
              localFiniteElement,
              child.geometryInFather(),
              localData);
        localFiniteElement.localInterpolation().interpolate(oldGridFunction,
            childLocalData);

        iterateOverLocalIndexSet(
          localIndexSet,
          [&](size_t i, auto gi)
          {
            data[gi[0]] = childLocalData[i];
          },
          [](size_t i){},
          [](size_t i, auto gi, double wi)
          { /* will be copied from neighboring element */ }
        );
      }
    }
  }
  return data;
}

} // end namespace Dune

#endif
