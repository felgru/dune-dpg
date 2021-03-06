// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_LOCALCOEFFICIENTS_HH
#define DUNE_DPG_LOCALCOEFFICIENTS_HH

#include <type_traits>

namespace Dune {
namespace Functions {
  template<class Range, class GridView>
  class ConstantGridViewFunction;
}

namespace detail {
  namespace LocalCoefficients {

    template<class F>
    class OnlyFactor {
      F factor_;
      typename F::LocalFunction localFactor_;

    public:
      using Factor = F;
      using LocalFactor = typename F::LocalFunction;
      using Element = typename F::Element;

      OnlyFactor(F f) : factor_(f), localFactor_(localFunction(factor_)) {};

      void bind(const Element& element)
      {
        localFactor_.bind(element);
      }

      const LocalFactor& localFactor() const {
        return localFactor_;
      }
    };

    template<class F, class D>
    class FactorAndDirection {
      F factor_;
      typename F::LocalFunction localFactor_;
      D direction_;
      typename D::LocalFunction localDirection_;

    public:
      using Factor = F;
      using LocalFactor = typename F::LocalFunction;
      using Direction = D;
      using LocalDirection = typename D::LocalFunction;
      using Element = typename F::Element;

      FactorAndDirection(F f, D d)
        : factor_(f), localFactor_(localFunction(factor_)),
          direction_(d), localDirection_(localFunction(direction_)) {};

      void bind(const Element& element)
      {
        localFactor_.bind(element);
        localDirection_.bind(element);
      }

      const LocalFactor& localFactor() const {
        return localFactor_;
      }

      const LocalDirection& localDirection() const {
        return localDirection_;
      }

      const LocalDirection& localSecondDirection() const {
        return localDirection_;
      }
    };

    template<class F, class D1, class D2>
    class FactorAndTwoDirections {
      F factor_;
      typename F::LocalFunction localFactor_;
      D1 direction1_;
      typename D1::LocalFunction localDirection1_;
      D2 direction2_;
      typename D2::LocalFunction localDirection2_;

    public:
      using Factor = F;
      using LocalFactor = typename F::LocalFunction;
      using Direction = D1;
      using LocalDirection = typename D1::LocalFunction;
      using SecondDirection = D2;
      using LocalSecondDirection = typename D2::LocalFunction;
      using Element = typename F::Element;

      FactorAndTwoDirections(F f, D1 d1, D2 d2)
        : factor_(f), localFactor_(localFunction(factor_)),
          direction1_(d1), localDirection1_(localFunction(direction1_)),
          direction2_(d2), localDirection2_(localFunction(direction2_)) {};

      void bind(const Element& element)
      {
        localFactor_.bind(element);
        localDirection1_.bind(element);
        localDirection2_.bind(element);
      }

      const LocalFactor& localFactor() const {
        return localFactor_;
      }

      const LocalDirection& localDirection() const {
        return localDirection1_;
      }

      const LocalSecondDirection& localSecondDirection() const {
        return localDirection2_;
      }
    };

    template<class Coeffs>
    struct uses_only_constant_coefficients : std::false_type {};

    template<class Function>
    struct is_ConstantGridViewFunction : std::false_type {};

    template<class Range, class GridView>
    struct is_ConstantGridViewFunction
      <Dune::Functions::ConstantGridViewFunction<Range, GridView>>
      : std::true_type {};

    template<class F>
    struct uses_only_constant_coefficients<OnlyFactor<F>>
      : is_ConstantGridViewFunction<F> {};

    template<class F, class D>
    struct uses_only_constant_coefficients<FactorAndDirection<F, D>>
      : std::conditional<
                 is_ConstantGridViewFunction<F>::value
              && is_ConstantGridViewFunction<D>::value,
             std::true_type, std::false_type>::type {};

    template<class F, class D1, class D2>
    struct uses_only_constant_coefficients<FactorAndTwoDirections<F, D1, D2>>
      : std::conditional<
                 is_ConstantGridViewFunction<F>::value
              && is_ConstantGridViewFunction<D1>::value
              && is_ConstantGridViewFunction<D2>::value,
             std::true_type, std::false_type>::type {};
  } // namespace LocalCoefficients
} // namespace detail
} // namespace Dune
#endif
