// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_PERITER_COMMON_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_PERITER_COMMON_HH

#include <ostream>
#include <stdexcept>
#include <type_traits>

#include <dune/dpg/assemble_helper.hh>
#include <dune/dpg/spacetuple.hh>

#include <dune/functions/functionspacebases/bernsteindgrefineddgnodalbasis.hh>
#include <dune/functions/functionspacebases/bernsteindgbasis.hh>

namespace Dune {

enum class PeriterPlotFlags : unsigned char {
  doNotPlot = 0,
  plotOuterIterations = 1 << 0,
  plotLastIteration = 1 << 1,
  plotScattering = 1 << 2,
  plotRhs = 1 << 3,
  plotIntegratedSolution = 1 << 4
};

constexpr inline
PeriterPlotFlags operator|(PeriterPlotFlags a, PeriterPlotFlags b) {
  using T = std::underlying_type_t<PeriterPlotFlags>;
  return static_cast<PeriterPlotFlags>(static_cast<T>(a) | static_cast<T>(b));
}

constexpr inline
PeriterPlotFlags operator&(PeriterPlotFlags a, PeriterPlotFlags b) {
  using T = std::underlying_type_t<PeriterPlotFlags>;
  return static_cast<PeriterPlotFlags>(static_cast<T>(a) & static_cast<T>(b));
}

constexpr inline
PeriterPlotFlags& operator|=(PeriterPlotFlags& a, PeriterPlotFlags b) {
  return a = (a | b);
}

struct FeRHS {};
struct ApproximateRHS {};

template<class GV, class TraceBasis>
class TransportSpaces {
  static_assert(std::is_same<typename TraceBasis::GridView, GV>::value,
                "GridViews of transport spaces don't match!");

  public:
  using GridView = GV;

  using FEBasisInterior = Functions::BernsteinDGBasis<GridView, 1>;
  using FEBasisTrace = TraceBasis;

  using FEBasisTest = Functions::BernsteinDGRefinedDGBasis<GridView, 1, 3>;
  using FEBasisEnrichedTest
      = Functions::BernsteinDGRefinedDGBasis<GridView, 1, 4>;

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

  GridView gridView() const {
    return interiorSolutionSpace().gridView();
  }

  private:

  SolutionSpacePtr solutionSpace_;
  TestSpacePtr testSpace_;
  EnrichedTestSpacePtr enrichedTestSpace_;
};

  class PeriterApproximationParameters
  {
    unsigned int n = 0;
    // η_n:
    double eta_ = 1;
    const double rho;
    const double rhobar;
    const double CT;
    const double err0;
    // CT*kappa1 + CT*kappa2 + 2*kappa3 = 1.
    const double kappa1;
    const double kappa2;
    const double kappa3;

    friend std::ostream& operator<<
        (std::ostream& os, const PeriterApproximationParameters& params);

    public:

    /**
     * \param accuracyRatio a value in (0,1) where larger values mean
     *        more accuracy for the transport solver and smaller values
     *        mean more accuracy for the kernel approximation
     */
    PeriterApproximationParameters(double accuracyRatio,
                                   double rho, double CT, double err0, FeRHS)
      : rho(rho)
      , rhobar(2./rho)
      , CT(CT)
      , err0(err0)
      , kappa1(accuracyRatio/CT)
      , kappa2(0.)
      , kappa3((1.-accuracyRatio)/2.)
    {
      if(accuracyRatio < 0. || accuracyRatio > 1.)
        throw std::domain_error("accuracyRatio needs to be in (0,1).");
    }

    PeriterApproximationParameters(double accuracyRatio,
                                   double rho, double CT, double err0,
                                   ApproximateRHS)
      : rho(rho)
      , rhobar(2./rho)
      , CT(CT)
      , err0(err0)
      , kappa1(accuracyRatio/CT)
      , kappa2((1.-accuracyRatio)/(2.*CT))
      , kappa3((1.-accuracyRatio)/4.)
    {}

    // TODO: The accuracy also depends on the kappa from K = κG and on \|u\|.
    //       Adding a factor 1/4. to compensate for that.
    double finalScatteringAccuracy(double targetAccuracy) const {
      const int m = maxOuterIterationsForTargetAccuracy(targetAccuracy);
      return kappa1*std::pow(rhobar, -m)/4.;
    }

    double aPrioriAccuracy() const {
      return err0;
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

    //! a priori estimate for $\|u - \bar u_{n+1}\|$
    double combinedAccuracy() const {
      return (rho*err0 + 2) * std::pow(rho,n);
    }

    unsigned int maxOuterIterationsForTargetAccuracy(double target) const
    {
      const double eps2 = target/(rho*err0+2);
      const int m = static_cast<int>(std::ceil(std::log(eps2) / std::log(rho)));
      return static_cast<unsigned int>(std::max(m, 0));
    }

    double aPosterioriErrorInLastOuterIteration(
        double aposterioriTransportGlobal) const
    {
      return aposterioriTransportGlobal + CT * scatteringAccuracy();
    }

    //! a posteriori estimate for $\|u_{n+1} - \bar u_{n+1}\|$
    double errorBetweenExactAndInexactIterate(
        const std::vector<double>& aposterioriIter) const
    {
      double errorAPosteriori = 0.;
      for(size_t j=0; j < n+1; j++) {
        errorAPosteriori += std::pow(rho,j)*aposterioriIter[n-j];
      }
      return errorAPosteriori;
    }

    //! a posteriori estimate for $\|u - \bar u_{n+1}\|$
    double combinedAPosterioriError(
        const std::vector<double>& aposterioriIter) const
    {
      return std::pow(rho,n+1)*err0
              + errorBetweenExactAndInexactIterate(aposterioriIter);
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
                           const PeriterApproximationParameters& params)
  {
    os << "rho = "    << params.rho    << '\n'
       << "rhobar = " << params.rhobar << '\n'
       << "kappa1 = " << params.kappa1 << '\n'
       << "kappa2 = " << params.kappa2 << '\n'
       << "kappa3 = " << params.kappa3 << '\n'
       << "CT = "     << params.CT     << '\n';
    return os;
  }
} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_PERITER_COMMON_HH
