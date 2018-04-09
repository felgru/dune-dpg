// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_SUBGRID_WORKAROUNDS_HH
#define DUNE_DPG_SUBGRID_WORKAROUNDS_HH

#include <memory>

namespace Dune {

template<class SubGrid>
std::unique_ptr<SubGrid> copySubGrid(const SubGrid& subGrid) {
  std::unique_ptr<SubGrid> gr
      = std::make_unique<SubGrid>(subGrid.getHostGrid());
  gr->createBegin();
  const auto gridView = subGrid.leafGridView();
  for(const auto& e : elements(gridView)) {
    gr->insert(subGrid.template getHostEntity<0>(e));
  }
  gr->createEnd();
  gr->setMaxLevelDifference(1);
  return gr;
}

} // end namespace Dune

#endif // DUNE_DPG_SUBGRID_WORKAROUNDS_HH
