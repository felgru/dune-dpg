// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_FUNCTIONS_REFINEMENTINTERPOLATION_HH
#define DUNE_DPG_FUNCTIONS_REFINEMENTINTERPOLATION_HH

#include <numeric>
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

  const size_t maxEntityDoFs = localView.maxSize();

  std::vector<EntitySeed> entitySeeds(gridView.size(0));
  std::vector<typename Vector::value_type>
      entityData(gridView.size(0) * maxEntityDoFs);

  auto localEntitySeed = entitySeeds.begin();
  auto localData = entityData.begin();
  for(const auto& e : elements(gridView))
  {
    localView.bind(e);

    *(localEntitySeed++) = e.seed();

    copyToLocalVector(data, localData, localView);
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

    struct Traits
    {
      using DomainType = Domain;
      using RangeType  = Range;
    };


    RestoreDataToRefinedGridFunction(
        const FiniteElement& finiteElement,
        const Geometry& childInElementGeometry,
        const Iterator elementData)
      : finiteElement(finiteElement),
        childInElementGeometry(childInElementGeometry),
        elementData(elementData)
    {}

    Range operator()(const Domain& x) const {
      const Domain xElement = childInElementGeometry.global(x);

      auto&& localBasis = finiteElement.localBasis();

      shapeFunctionValues.resize(localBasis.size());
      localBasis.evaluateFunction(xElement, shapeFunctionValues);

      return std::inner_product(shapeFunctionValues.cbegin(),
                                shapeFunctionValues.cend(),
                                elementData, Range(0));
    }

    const FiniteElement& finiteElement;
    const Geometry childInElementGeometry;
    const Iterator elementData;
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

  auto& [entitySeeds, entityData] = storedData;
  const auto& grid = globalBasis.gridView().grid();
  auto localView = globalBasis.localView();

  const size_t maxEntityDoFs = localView.maxSize();

  std::vector<FieldVector<double, 1>> data(globalBasis.size());

  auto localData = entityData.cbegin();
  for(auto entitySeed = entitySeeds.cbegin(),
           entitySeedEnd = entitySeeds.cend();
      entitySeed != entitySeedEnd;
      ++entitySeed, localData += maxEntityDoFs)
  {
    const auto e = grid.entity(*entitySeed);
    localView.bind(e);

    if(e.isLeaf()) {
      // directly copy cell data
      iterateOverLocalIndices(
        localView,
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
        localView.bind(child);

        // This assumes that e and child share the same finite element
        // and thus the same entity type.
        auto&& localFiniteElement = localView.tree().finiteElement();
        assert(child.father() == e);
        auto oldGridFunction = detail::RestoreDataToRefinedGridFunction
          <GlobalBasis, typename Element::LocalGeometry, decltype(localData)>(
              localFiniteElement,
              child.geometryInFather(),
              localData);
        localFiniteElement.localInterpolation().interpolate(oldGridFunction,
            childLocalData);

        iterateOverLocalIndices(
          localView,
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
