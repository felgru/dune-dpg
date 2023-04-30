// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_ASTI_COMMON_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_ASTI_COMMON_HH

#include <cmath>
#include <ostream>
#include <stdexcept>
#include <type_traits>

#include <dune/dpg/assemble_helper.hh>
#if ASTI_NORMALIZED_SPACES
#  include <dune/dpg/functions/normalizedspaces.hh>
#  include <dune/dpg/integralterm.hh>
#endif
#include <dune/dpg/spacetuple.hh>

#include <dune/functions/functionspacebases/bernsteindgrefineddgnodalbasis.hh>
#include <dune/functions/functionspacebases/bernsteindgbasis.hh>

namespace Dune {

enum class ASTIPlotFlags : unsigned char {
  doNotPlot = 0,
  plotOuterIterations = 1 << 0,
  plotLastIteration = 1 << 1,
  plotScattering = 1 << 2,
  plotRhs = 1 << 3,
  plotIntegratedSolution = 1 << 4
};

constexpr inline
ASTIPlotFlags operator|(ASTIPlotFlags a, ASTIPlotFlags b) {
  using T = std::underlying_type_t<ASTIPlotFlags>;
  return static_cast<ASTIPlotFlags>(static_cast<T>(a) | static_cast<T>(b));
}

constexpr inline
ASTIPlotFlags operator&(ASTIPlotFlags a, ASTIPlotFlags b) {
  using T = std::underlying_type_t<ASTIPlotFlags>;
  return static_cast<ASTIPlotFlags>(static_cast<T>(a) & static_cast<T>(b));
}

constexpr inline
ASTIPlotFlags& operator|=(ASTIPlotFlags& a, ASTIPlotFlags b) {
  return a = (a | b);
}

constexpr inline
bool flagIsSet(ASTIPlotFlags flags, ASTIPlotFlags flag) {
  using T = std::underlying_type_t<ASTIPlotFlags>;
  return static_cast<T>(flags & flag) != T{0};
}

enum class ASTILogFlags : unsigned char {
  defaultLog = 0,
  logDirectionSolveTime = 1 << 0,
  logKernelApproximationTime = 1 << 1,
  logTimings = 0b11,
};

constexpr inline
ASTILogFlags operator|(ASTILogFlags a, ASTILogFlags b) {
  using T = std::underlying_type_t<ASTILogFlags>;
  return static_cast<ASTILogFlags>(static_cast<T>(a) | static_cast<T>(b));
}

constexpr inline
ASTILogFlags operator&(ASTILogFlags a, ASTILogFlags b) {
  using T = std::underlying_type_t<ASTILogFlags>;
  return static_cast<ASTILogFlags>(static_cast<T>(a) & static_cast<T>(b));
}

constexpr inline
ASTILogFlags& operator|=(ASTILogFlags& a, ASTILogFlags b) {
  return a = (a | b);
}

constexpr inline
bool flagIsSet(ASTILogFlags flags, ASTILogFlags flag) {
  using T = std::underlying_type_t<ASTILogFlags>;
  return static_cast<T>(flags & flag) != T{0};
}

struct FeRHS {};
struct ApproximateRHS {};

template<class GV, class TraceBasis>
class TransportSpaces {
  static_assert(std::is_same<typename TraceBasis::GridView, GV>::value,
                "GridViews of transport spaces don't match!");

  template<typename TestSpacesPtr>
  static auto
  make_test_inner_product(const TestSpacesPtr& testSpaces,
                          const FieldVector<double, 2> direction)
  {
    auto innerProduct
      = innerProductWithSpace(testSpaces)
        .template addIntegralTerm<0,0,IntegrationType::gradGrad,
                                      DomainOfIntegration::interior>
                                 (1., direction)
        .template addIntegralTerm<0,0,IntegrationType::valueValue,
                                      DomainOfIntegration::interior>(1.)
        .create();

    return innerProduct;
  }

#if ASTI_NORMALIZED_SPACES
  template<typename FEBasisTest>
  static auto
  make_test_spaces(const typename FEBasisTest::GridView& gridView,
                   const FieldVector<double, 2> direction)
  {
    auto unnormalizedTestSpaces = make_space_tuple<FEBasisTest>(gridView);
    auto innerProduct
      = make_test_inner_product(unnormalizedTestSpaces, direction);

    return make_normalized_space_tuple(innerProduct);
  }
#endif

  template<typename FEBasisInterior, typename FEBasisTrace>
  static auto
  make_solution_spaces(const typename FEBasisInterior::GridView& gridView)
  {
#if ASTI_NORMALIZED_SPACES
    auto interiorSpace = make_space_tuple<FEBasisInterior>(gridView);
    auto l2InnerProduct
      = innerProductWithSpace(interiorSpace)
        .template addIntegralTerm<0,0,IntegrationType::valueValue,
                                      DomainOfIntegration::interior>(1.)
        .create();
    auto normedSpace = make_normalized_space(l2InnerProduct);
    using NormedSpace = decltype(normedSpace);

    return std::make_shared<std::tuple<NormedSpace, FEBasisTrace>>(
        std::make_tuple(std::move(normedSpace), FEBasisTrace(gridView)));
#else
    return make_space_tuple<FEBasisInterior, FEBasisTrace>(gridView);
#endif
  }

  using UnnormalizedFEBasisInterior = Functions::BernsteinDGBasis<GV, 1>;
  using UnnormalizedFEBasisTrace = TraceBasis;

  using UnnormalizedFEBasisTest
      = Functions::BernsteinDGRefinedDGBasis<GV, 1, 3>;
  using UnnormalizedFEBasisEnrichedTest
      = Functions::BernsteinDGRefinedDGBasis<GV, 1, 4>;

  public:
  using GridView = GV;

  using SolutionSpacePtr
    = decltype(make_solution_spaces<UnnormalizedFEBasisInterior,
                                    UnnormalizedFEBasisTrace>
                                   (std::declval<GridView>()));
#if ASTI_NORMALIZED_SPACES
  using TestSpacePtr = decltype(make_test_spaces<UnnormalizedFEBasisTest>
                                (std::declval<GridView>(),
                                 std::declval<FieldVector<double, 2>>()));
  using EnrichedTestSpacePtr
      = decltype(make_test_spaces<UnnormalizedFEBasisEnrichedTest>
                 (std::declval<GridView>(),
                  std::declval<FieldVector<double, 2>>()));
#else
  using TestSpacePtr = decltype(make_space_tuple<UnnormalizedFEBasisTest>
                                (std::declval<GridView>()));
  using EnrichedTestSpacePtr
      = decltype(make_space_tuple<UnnormalizedFEBasisEnrichedTest>
                 (std::declval<GridView>()));
#endif

  using FEBasisInterior
      = std::tuple_element_t<0, typename SolutionSpacePtr::element_type>;
  using FEBasisTrace
      = std::tuple_element_t<1, typename SolutionSpacePtr::element_type>;
  using FEBasisTest
      = std::tuple_element_t<0, typename TestSpacePtr::element_type>;
  using FEBasisEnrichedTest
      = std::tuple_element_t<0, typename EnrichedTestSpacePtr::element_type>;

#if ASTI_NORMALIZED_SPACES
  TransportSpaces(const GridView& gridView, FieldVector<double, 2> direction)
#else
  TransportSpaces(const GridView& gridView)
#endif
    : solutionSpace_(make_solution_spaces<UnnormalizedFEBasisInterior,
                                          UnnormalizedFEBasisTrace>(gridView))
#if ASTI_NORMALIZED_SPACES
    , testSpace_(make_test_spaces<UnnormalizedFEBasisTest>(gridView, direction))
    , enrichedTestSpace_(make_test_spaces<UnnormalizedFEBasisEnrichedTest>
                                         (gridView, direction))
#else
    , testSpace_(make_space_tuple<FEBasisTest>(gridView))
    , enrichedTestSpace_(make_space_tuple<FEBasisEnrichedTest>(gridView))
#endif
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

  auto testInnerProduct(const FieldVector<double, 2> direction) const {
    return make_test_inner_product(testSpace_, direction);
  }

  private:

  SolutionSpacePtr solutionSpace_;
  TestSpacePtr testSpace_;
  EnrichedTestSpacePtr enrichedTestSpace_;
};

  class ASTIApproximationParameters
  {
    unsigned int n = 0;
    // η_n:
    double eta_ = 1;
    const double beta = 1.5;
    const double rho;
    const double CT;
    double uBound; // b(u) in our paper
    // CT*kappa1 + CT*kappa2 + 2*kappa3 = 1.
    const double kappa1;
    const double kappa2;
    const double kappa3;
    const size_t scatteringTruthLevel_;

    friend std::ostream& operator<<
        (std::ostream& os, const ASTIApproximationParameters& params);

    double etaInStep(unsigned int m) const {
      return std::pow(1+m, -beta) * std::pow(rho, m);
    }

    public:

    /**
     * \param accuracyRatio a value in (0,1) where larger values mean
     *        more accuracy for the transport solver and smaller values
     *        mean more accuracy for the kernel approximation
     */
    ASTIApproximationParameters(double accuracyRatio,
                                   double rho, double CT, double err0,
                                   size_t scatteringTruthLevel, FeRHS)
      : rho(rho)
      , CT(CT)
      , uBound(err0)
      , kappa1(accuracyRatio/CT)
      , kappa2(0.)
      , kappa3(1.-accuracyRatio)
      , scatteringTruthLevel_(scatteringTruthLevel)
    {
      if(accuracyRatio < 0. || accuracyRatio > 1.)
        throw std::domain_error("accuracyRatio needs to be in (0,1).");
    }

    ASTIApproximationParameters(double accuracyRatio,
                                   double rho, double CT, double err0,
                                   size_t scatteringTruthLevel,
                                   ApproximateRHS)
      : rho(rho)
      , CT(CT)
      , uBound(err0)
      , kappa1(accuracyRatio/CT)
      , kappa2((1.-accuracyRatio)/(2.*CT))
      , kappa3((1.-accuracyRatio)/2.)
      , scatteringTruthLevel_(scatteringTruthLevel)
    {}

    // TODO: The accuracy also depends on the kappa from K = κG and on \|u\|.
    //       Adding a factor 1/4. to compensate for that.
    double finalScatteringAccuracy(double targetAccuracy) const {
      const int m = maxOuterIterationsForTargetAccuracy(targetAccuracy);
      return kappa1*etaInStep(m)/4.;
    }

    size_t scatteringTruthLevel() const {
      return scatteringTruthLevel_;
    }

    double aPrioriAccuracy() const {
      return uBound;
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
      // TODO: here we can replace riemann_zeta(beta) with the sum
      //       over the previous eta_j
      return (rho*uBound + std::riemann_zeta(beta)) * std::pow(rho,n);
    }

    unsigned int maxOuterIterationsForTargetAccuracy(double target) const
    {
      const double eps2 = target/(rho*uBound+std::riemann_zeta(beta));
      const int m = static_cast<int>(std::ceil(std::log(eps2) / std::log(rho)));
      return static_cast<unsigned int>(std::max(m, 0));
    }

    //! a posteriori estimate for $\|u_{n+1} - \bar u_{n+1}\|$
    double errorBetweenExactAndInexactIterate(
        double previousDeviationOfInexactIterate,
        double transportError) const
    {
      if(n > 0) {
        return rho * previousDeviationOfInexactIterate
             + CT * (kappa1 + kappa2) * eta_
             + transportError;
      } else {
        // In the first iteration we have K u_0 = u_0 = 0
        return rho * previousDeviationOfInexactIterate
             + CT * kappa2 * eta_
             + transportError;
      }
    }

    //! a posteriori estimate for $\|u - \bar u_{n+1}\|$
    double combinedAPosterioriError(
        double deviationOfInexactIterate) const
    {
      return std::pow(rho,n+1)*uBound + deviationOfInexactIterate;
    }

    double eta() const {
      return eta_;
    }

    void decreaseEta(double uNorm) {
      updateUBound(uNorm);
      n++;
      eta_ = etaInStep(n);
    }

    private:

    void updateUBound(const double uNorm) {
      // TODO: here we can probably replace riemann_zeta(beta) with the sum
      //       over the previous eta_j
      const double uBoundUpdate = uNorm
          + (rho*uBound + std::riemann_zeta(beta)) * std::pow(rho,n-1);
      uBound = std::min(uBound, uBoundUpdate);
    }
  };

  std::ostream& operator<<(std::ostream& os,
                           const ASTIApproximationParameters& params)
  {
    os << "rho = "    << params.rho    << '\n'
       << "kappa1 = " << params.kappa1 << '\n'
       << "kappa2 = " << params.kappa2 << '\n'
       << "kappa3 = " << params.kappa3 << '\n'
       << "CT = "     << params.CT     << '\n';
    return os;
  }
} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_ASTI_COMMON_HH
