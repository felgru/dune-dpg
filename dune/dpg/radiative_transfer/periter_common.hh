// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_PERITER_COMMON_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_PERITER_COMMON_HH

#include <type_traits>

namespace Dune {

enum class PlotSolutions : unsigned char {
  doNotPlot = 0,
  plotOuterIterations = 1 << 0,
  plotLastIteration = 1 << 1,
  plotScattering = 1 << 2
};

constexpr inline PlotSolutions operator|(PlotSolutions a, PlotSolutions b) {
  using T = std::underlying_type_t<PlotSolutions>;
  return static_cast<PlotSolutions>(static_cast<T>(a) | static_cast<T>(b));
}

constexpr inline PlotSolutions operator&(PlotSolutions a, PlotSolutions b) {
  using T = std::underlying_type_t<PlotSolutions>;
  return static_cast<PlotSolutions>(static_cast<T>(a) & static_cast<T>(b));
}

} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_PERITER_COMMON_HH
