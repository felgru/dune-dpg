// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_SPACETUPLE_HH
#define DUNE_DPG_SPACETUPLE_HH

#include <memory>
#include <tuple>
#include <type_traits>

namespace Dune {

  template<class... Spaces>
  using SpaceTuplePtr = std::shared_ptr<std::tuple<Spaces...>>;

  template<class T>
  struct is_SpaceTuplePtr : std::false_type {};

  template<class... Spaces>
  struct is_SpaceTuplePtr<SpaceTuplePtr<Spaces...>> : std::true_type {};

} // end namespace Dune

#endif
