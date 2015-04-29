// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/* Convenience interface for enable_if, taken from
 * http://flamingdangerzone.com/cxx11/2012/06/01/almost-static-if.html */
#ifndef DUNE_DPG_TRAITS_HH
#define DUNE_DPG_TRAITS_HH

#include<type_traits>

namespace Dune {

namespace detail {
  enum class enabler {};
}

template <typename Condition>
using EnableIf =
  typename std::enable_if<Condition::value, detail::enabler>::type;

} // end namespace Dune

#endif // DUNE_DPG_TRAITS_HH
