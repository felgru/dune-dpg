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

} // end namespace Dune

#endif
