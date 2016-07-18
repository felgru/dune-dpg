// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_HENYEY_GREENSTEIN_SCATTERING_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_HENYEY_GREENSTEIN_SCATTERING_HH

#include <cmath>
#include <functional>
#include <boost/math/constants/constants.hpp>

namespace Dune {

template<class Direction>
std::function<double(const Direction& s1, const Direction& s2)>
HenyeyGreensteinScattering(double gamma) {
  return [gamma](const Direction& s1, const Direction& s2) {
    double scalarProduct = s1 * s2;
    if (scalarProduct > 1) scalarProduct = 1;
    if (scalarProduct < 0) scalarProduct = 0;
    using namespace boost::math::constants;
    return 1./(2*pi<double>())
           *(1-gamma*gamma)/(1+gamma*gamma-2*gamma*scalarProduct);
  };
}

}

#endif // DUNE_DPG_RADIATIVE_TRANSFER_HENYEY_GREENSTEIN_SCATTERING_HH
