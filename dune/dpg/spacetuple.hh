// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_SPACETUPLE_HH
#define DUNE_DPG_SPACETUPLE_HH

#include <memory>
#include <tuple>
#include <type_traits>

#include <dune/dpg/functions/concepts.hh>

namespace Dune {

  template<class... Spaces>
  using SpaceTuplePtr = std::shared_ptr<std::tuple<Spaces...>>;

  template<class T>
  struct is_SpaceTuplePtr : std::false_type {};

  template<class... Spaces>
  struct is_SpaceTuplePtr<SpaceTuplePtr<Spaces...>> : std::true_type {};

  /**
   * Create a shared_ptr to a tuple of spaces over the same GridView
   */
  template<class... Spaces, class GridView>
  SpaceTuplePtr<Spaces...>
  make_space_tuple(const GridView& gridView)
  {
    static_assert(Concept::tupleEntriesModel<
        Functions::Concept::GeneralizedGlobalBasis<GridView>,
        std::tuple<Spaces...>>(),
        "Spaces need to model the GeneralizedGlobalBasis concept.");

    return std::make_shared<typename std::tuple<Spaces...>>(
        std::make_tuple(Spaces(gridView)...));
  }

} // end namespace Dune

#endif
