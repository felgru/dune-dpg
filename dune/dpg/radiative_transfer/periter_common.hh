// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_PERITER_COMMON_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_PERITER_COMMON_HH

#include <ostream>
#include <type_traits>

#include <dune/dpg/assemble_helper.hh>
#include <dune/dpg/spacetuple.hh>

#include <dune/functions/functionspacebases/bernsteindgrefineddgnodalbasis.hh>
#include <dune/functions/functionspacebases/bernsteindgbasis.hh>

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

constexpr inline PlotSolutions& operator|=(PlotSolutions& a, PlotSolutions b) {
  return a = (a | b);
}

struct FeRHS {};
struct ApproximateRHS {};

template<class GV, class TraceBasis>
class TransportSpaces {
  public:
  using GridView = GV;

  using FEBasisInterior = Functions::BernsteinDGBasis<GridView, 1>;
  using FEBasisTrace = TraceBasis;
  static_assert(std::is_same<typename TraceBasis::GridView, GridView>::value,
                "GridViews of transport spaces don't match!");

  using FEBasisTest = Functions::BernsteinDGRefinedDGBasis<GridView, 1, 3>;
  using FEBasisEnrichedTest = FEBasisTest;

  using SolutionSpacePtr = decltype(make_space_tuple<FEBasisInterior,
                                    FEBasisTrace>(std::declval<GridView>()));
  using TestSpacePtr = decltype(make_space_tuple<FEBasisTest>
                                (std::declval<GridView>()));
  using EnrichedTestSpacePtr = decltype(make_space_tuple<FEBasisEnrichedTest>
                                        (std::declval<GridView>()));

  TransportSpaces(const GridView& gridView)
    : solutionSpace_(make_space_tuple<FEBasisInterior, FEBasisTrace>(gridView))
    , testSpace_(make_space_tuple<FEBasisTest>(gridView))
    , enrichedTestSpace_(make_space_tuple<FEBasisEnrichedTest>(gridView))
  {}

  void update(const GridView& gridView) {
    detail::updateSpaces(*solutionSpace_,     gridView);
    detail::updateSpaces(*testSpace_,         gridView);
    detail::updateSpaces(*enrichedTestSpace_, gridView);
  }

  const FEBasisInterior& interiorSolutionSpace() const {
    return std::get<0>(*solutionSpace_);
  }

  const FEBasisTrace& traceSolutionSpace() const {
    return std::get<1>(*solutionSpace_);
  }

  const FEBasisTest& testSpace() const {
    return std::get<0>(*testSpace_);
  }

  const SolutionSpacePtr& solutionSpacePtr() const {
    return solutionSpace_;
  }

  const TestSpacePtr& testSpacePtr() const {
    return testSpace_;
  }

  const EnrichedTestSpacePtr& enrichedTestSpacePtr() const {
    return enrichedTestSpace_;
  }

  private:

  SolutionSpacePtr solutionSpace_;
  TestSpacePtr testSpace_;
  EnrichedTestSpacePtr enrichedTestSpace_;
};

namespace detail {
  class ApproximationParameters
  {
    unsigned int n = 0;
    // η_n:
    double eta_ = 1;
    const double rho;
    const double rhobar;
    const double err0;
    // CT*kappa1 + CT*kappa2 + 2*kappa3 = 1.
    const double kappa1;
    const double kappa2;
    const double kappa3;

    friend std::ostream& operator<<(std::ostream& os,
                                    const ApproximationParameters& params);

    public:

    ApproximationParameters(double rho, double CT, double err0, FeRHS)
      : rho(rho)
      , rhobar(2./rho)
      , err0(err0)
      , kappa1(1./(2.*CT))
      , kappa2(0.)
      , kappa3(1./4.)
    {}

    ApproximationParameters(double rho, double CT, double err0, ApproximateRHS)
      : rho(rho)
      , rhobar(2./rho)
      , err0(err0)
      , kappa1(1./(3.*CT))
      , kappa2(1./(3.*CT))
      , kappa3(1./6.)
    {}

    // TODO: The accuracy also depends on the kappa from K = κG and on \|u\|.
    //       Adding a factor 1/4. to compensate for that.
    double finalScatteringAccuracy(double targetAccuracy) const {
      return kappa1*targetAccuracy/4.;
    }

    double scatteringAccuracy() const {
      return kappa1 * eta_;
    }

    double rhsAccuracy() const {
      return kappa2 * eta_;
    }

    double transportAccuracy() const {
      return kappa3 * eta_;
    }

    double combinedAccuracy() const {
      return (rho*err0 + 2) * std::pow(rho,n);
    }

    double aPosterioriError(const std::vector<double>& aposterioriIter) const {
      double errorAPosteriori = 0.;
      for(size_t j=0; j < n+1; j++) {
        errorAPosteriori += std::pow(rho,j)*aposterioriIter[n-j];
      }
      return errorAPosteriori;
    }

    double eta() const {
      return eta_;
    }

    void decreaseEta() {
      eta_ /= rhobar;
      n++;
    }
  };

  std::ostream& operator<<(std::ostream& os,
                           const ApproximationParameters& params)
  {
    os << "rho = "    << params.rho    << '\n'
       << "rhobar = " << params.rhobar << '\n'
       << "kappa1 = " << params.kappa1 << '\n'
       << "kappa2 = " << params.kappa2 << '\n'
       << "kappa3 = " << params.kappa3 << '\n';
    return os;
  }
} // end namespace detail
} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_PERITER_COMMON_HH
