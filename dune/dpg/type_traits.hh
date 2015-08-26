// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_TYPE_TRAITS_HH
#define DUNE_DPG_TYPE_TRAITS_HH

#include<type_traits>

namespace Dune {

namespace Functions {
  template<class GV, int s, int k>
  class PQKSubsampledDGBasis;
}

/* Convenience interface for enable_if, taken from
 * http://flamingdangerzone.com/cxx11/2012/06/01/almost-static-if.html */

namespace detail {
  enum class enabler {};
}

template <typename Condition>
using EnableIf =
  typename std::enable_if<Condition::value, detail::enabler>::type;


template <typename FiniteElement>
struct is_SubsampledFiniteElement : std::false_type {};

template<class GV, int s, int k>
struct is_SubsampledFiniteElement<Functions::PQKSubsampledDGBasis
                                  <GV, s, k> > : std::true_type {};

template <typename FiniteElement>
struct numberOfSamples : std::integral_constant<int, 1> {};

template<class GV, int s, int k>
struct numberOfSamples<Functions::PQKSubsampledDGBasis<GV, s, k> >
                : std::integral_constant<int, s> {};

} // end namespace Dune

#endif // DUNE_DPG_TYPE_TRAITS_HH
