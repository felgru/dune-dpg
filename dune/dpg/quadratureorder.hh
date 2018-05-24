// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_QUADRATUREORDER_HH
#define DUNE_DPG_QUADRATUREORDER_HH

#include <dune/common/std/type_traits.hh>

namespace Dune {

template<class LocalFunction>
struct requiredQuadratureOrder;

namespace detail {
  template<class T>
  using has_requiredQuadratureOrderValue_t
    = decltype(requiredQuadratureOrder<T>::value);
}

template<class LocalFunction>
struct hasRequiredQuadratureOrder :
  Std::is_detected<detail::has_requiredQuadratureOrderValue_t, LocalFunction>
{};

} // End namespace Dune

#endif
