// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_HENYEY_GREENSTEIN_SCATTERING_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_HENYEY_GREENSTEIN_SCATTERING_HH

#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <string>

namespace Dune {

struct HenyeyGreensteinScattering {
  HenyeyGreensteinScattering(double gamma) : gamma(gamma) {}

  double operator()(double angle) const {
      return (1-gamma*gamma)/(2 * boost::math::constants::pi<double>()
                              * (1+gamma*gamma-2*gamma*std::cos(angle)));
  }

  std::string info() const {
    return "Henyey Greenstein kernel with gamma = " + std::to_string(gamma);
  }

private:
  const double gamma;
};

}
#endif // DUNE_DPG_RADIATIVE_TRANSFER_HENYEY_GREENSTEIN_SCATTERING_HH
