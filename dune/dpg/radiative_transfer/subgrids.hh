// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_SUBGRIDS_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_SUBGRIDS_HH

#include <algorithm>
#include <memory>
#include <set>
#include <type_traits>

#include <dune/grid/common/rangegenerators.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <dune/subgrid/subgrid.hh>
#pragma GCC diagnostic pop

namespace Dune {

template<class SubGrid, class HostGrid>
std::unique_ptr<SubGrid> fullSubGrid(HostGrid& hostGrid) {
  std::unique_ptr<SubGrid> gr
      = std::make_unique<SubGrid>(hostGrid);
  gr->createBegin();
  const auto gridView = hostGrid.leafGridView();
  for(const auto& e : elements(gridView)) {
    gr->insert(e);
  }
  gr->createEnd();
  gr->setMaxLevelDifference(1);
  return gr;
}

/**
 * compute intersection of two SubGrids
 */
template<class SubGrid>
std::unique_ptr<SubGrid> intersectSubGrids(const SubGrid& subGrid1,
                                           const SubGrid& subGrid2)
{
  // assert(subGrid1.getHostGrid() == subGrid2.getHostGrid());
  std::unique_ptr<SubGrid> gr
      = std::make_unique<SubGrid>(subGrid1.getHostGrid());
  gr->createBegin();
  const int maxLevel = std::min(subGrid1.maxLevel(), subGrid2.maxLevel());
  for(int level=0; level <= maxLevel; ++level)
  {
    const auto levelGridView = subGrid1.levelGridView(level);
    for (auto&& e : elements(levelGridView))
    {
      const auto eHost = subGrid1.template getHostEntity<0>(e);
      if(subGrid2.template contains<0>(eHost))
      {
        gr->insertRaw(eHost);
      }
    }
  }
  gr->createEnd();
  gr->setMaxLevelDifference(1);
  return gr;
}

template<class SubGrid>
std::set<typename SubGrid::HostGridType::GlobalIdSet::IdType>
saveSubGridToIdSet(const SubGrid& subGrid)
{
  auto& idSet = subGrid.getHostGrid().globalIdSet();
  auto subGridView = subGrid.leafGridView();
  std::set<typename std::decay_t<decltype(idSet)>::IdType> subGridElements;
  for(const auto& e : elements(subGridView)) {
    subGridElements.insert(idSet.id(subGrid.template getHostEntity<0>(e)));
  }
  return subGridElements;
}

template<class SubGrid>
std::unique_ptr<SubGrid>
restoreSubGridFromIdSet(
    std::set<typename SubGrid::HostGridType::GlobalIdSet::IdType>&& idSet,
    typename SubGrid::HostGridType& hostGrid)
{
  const int numberOfLeafElements = idSet.size();
  auto subGrid = std::make_unique<SubGrid>(hostGrid);
  subGrid->createBegin();
  subGrid->insertSetPartial(idSet);
  subGrid->createEnd();
  subGrid->setMaxLevelDifference(1);
  if(subGrid->leafGridView().size(0) != numberOfLeafElements) {
    DUNE_THROW(Dune::Exception, "restored SubGrid does not match idSet!");
  }
  return subGrid;
}

template<class SubGrid>
std::unique_ptr<SubGrid>
restoreSubGridFromIdSet(
    const std::set<typename SubGrid::HostGridType::GlobalIdSet::IdType>& idSet,
    typename SubGrid::HostGridType& hostGrid)
{
  std::set<typename SubGrid::HostGridType::GlobalIdSet::IdType>
    idSetCopy(idSet);
  return restoreSubGridFromIdSet<SubGrid>(std::move(idSetCopy), hostGrid);
}

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
  if(gr->leafGridView().size(0) != subGrid.leafGridView().size(0)) {
    DUNE_THROW(Dune::Exception, "copied SubGrid does not match original one!");
  }
  return gr;
}

} // end namespace Dune

#endif
