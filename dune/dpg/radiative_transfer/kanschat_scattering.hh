// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_KANSCHAT_SCATTERING_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_KANSCHAT_SCATTERING_HH

#include <cmath>
#include <functional>

namespace Dune {

template<class Direction>
std::function<double(const Direction& s1, const Direction& s2)>
KanschatScattering(std::vector<double>&& a) {
  return [a](const Direction& s1, const Direction& s2) {
    double scalarProduct = s1 * s2;
    if (scalarProduct > 1) scalarProduct = 1;
    if (scalarProduct < 0) scalarProduct = 0;
    const double phi = acos(scalarProduct);
    double scattering = 0;
    for (unsigned int i = 0, asize = a.size(); i < asize; ++i) {
      scattering += a[i]*cos(i*phi);
    }
    return scattering;
  };
}

template<class Direction>
std::function<double(const Direction& s1, const Direction& s2)>
ACP2011Scattering() {
  return [](const Direction& s1, const Direction& s2) {
    return 1 + 0.2*(s1 * s2);
  };
}

}

#endif // DUNE_DPG_RADIATIVE_TRANSFER_KANSCHAT_SCATTERING_HH
