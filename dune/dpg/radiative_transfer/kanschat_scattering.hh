// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_KANSCHAT_SCATTERING_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_KANSCHAT_SCATTERING_HH

#include <cmath>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace Dune {

struct KanschatScattering {
  KanschatScattering(std::vector<double>&& a) : a(std::move(a)) {}

  double operator()(double angle) const {
    double scattering = 0;
    for (unsigned int i = 0, asize = a.size(); i < asize; ++i) {
      scattering += a[i]*std::cos(i*angle);
    }
    return scattering;
  };

  std::string info() const {
    std::stringstream buf;
    buf << "Kanschat scattering kernel with a = {" << a[0];
    for(size_t i = 1, imax = a.size(); i < imax; i++)
      buf << ", " << a[i];
    buf << '}';
    return buf.str();
  }

private:
  const std::vector<double> a;
};

struct ACP2011Scattering {
  double operator() (double angle) const {
    return 1 + 0.2*std::cos(angle);
  };

  std::string info() const {
    return "scattering kernel from Avila, Codina and Principe 2011";
  }
};

}

#endif // DUNE_DPG_RADIATIVE_TRANSFER_KANSCHAT_SCATTERING_HH
