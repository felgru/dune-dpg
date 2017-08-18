// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_HENYEY_GREENSTEIN_SCATTERING_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_HENYEY_GREENSTEIN_SCATTERING_HH

#include <cmath>
#include <functional>
#include <boost/math/constants/constants.hpp>
#include <dune/common/exceptions.hh>

namespace Dune {

inline auto
HenyeyGreensteinScattering(double gamma) {
  return [gamma](double angle) {
    return (1-gamma*gamma)/(2 * boost::math::constants::pi<double>()
                              * (1+gamma*gamma-2*gamma*cos(angle)));
  };
}

}

#endif // DUNE_DPG_RADIATIVE_TRANSFER_HENYEY_GREENSTEIN_SCATTERING_HH
