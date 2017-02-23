// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_HENYEY_GREENSTEIN_SCATTERING_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_HENYEY_GREENSTEIN_SCATTERING_HH

#include <cmath>
#include <functional>

#include <dune/common/exceptions.hh>
#include <dune/common/deprecated.hh>
DUNE_NO_DEPRECATED_BEGIN
// Older Boost versions spill a lot of deprecation warnings as they
// still use std::auto_ptr.
#include <boost/math/constants/constants.hpp>
DUNE_NO_DEPRECATED_END

namespace Dune {

template<class Direction>
std::function<double(const Direction&, const Direction&)>
HenyeyGreensteinScattering(double gamma) {
  if (Direction().dim() == 2) {
    return [gamma](const Direction& s1, const Direction& s2) {
      double scalarProduct = s1 * s2;
      if (scalarProduct > 1) scalarProduct = 1;
      if (scalarProduct < 0) scalarProduct = 0;
      return (1-gamma*gamma)/(1+gamma*gamma-2*gamma*scalarProduct);
    };
  } else if (Direction().dim() == 3) {
    return [gamma](const Direction& s1, const Direction& s2) {
      double scalarProduct = s1 * s2;
      if (scalarProduct > 1) scalarProduct = 1;
      if (scalarProduct < 0) scalarProduct = 0;
      return (1-gamma*gamma)/pow(1+gamma*gamma-2*gamma*scalarProduct,3./2.);
    };
  } else {
    DUNE_THROW(Dune::NotImplemented,
        "Henyey-Greenstein scattering only implemented in 2d and 3d.");
  }
}

}

#endif // DUNE_DPG_RADIATIVE_TRANSFER_HENYEY_GREENSTEIN_SCATTERING_HH
