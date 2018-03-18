// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_LOCALCOEFFICIENTS_HH
#define DUNE_DPG_LOCALCOEFFICIENTS_HH

namespace Dune {
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

    public:
      using Factor = F;
      using LocalFactor = typename F::LocalFunction;
      using Direction = D;
      using Element = typename F::Element;

      FactorAndDirection(F f, D d)
        : factor_(f), localFactor_(localFunction(factor_)), direction_(d) {};

      void bind(const Element& element)
      {
        localFactor_.bind(element);
      }

      const LocalFactor& localFactor() const {
        return localFactor_;
      }

      const Direction& direction() const {
        return direction_;
      }

      const Direction& secondDirection() const {
        return direction_;
      }
    };

    template<class F, class D1, class D2>
    class FactorAndTwoDirections {
      F factor_;
      typename F::LocalFunction localFactor_;
      D1 direction1_;
      D2 direction2_;

    public:
      using Factor = F;
      using LocalFactor = typename F::LocalFunction;
      using Direction = D1;
      using SecondDirection = D2;
      using Element = typename F::Element;

      FactorAndTwoDirections(F f, D1 d1, D2 d2)
        : factor_(f), localFactor_(localFunction(factor_)),
          direction1_(d1), direction2_(d2) {};

      void bind(const Element& element)
      {
        localFactor_.bind(element);
      }

      const LocalFactor& localFactor() const {
        return localFactor_;
      }

      const Direction& direction() const {
        return direction1_;
      }

      const Direction& secondDirection() const {
        return direction2_;
      }
    };

  } // namespace LocalCoefficients
} // namespace detail
} // namespace Dune
#endif
