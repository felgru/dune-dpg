// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NORMALIZEDBASISADAPTOR_FINITEELEMENT_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NORMALIZEDBASISADAPTOR_FINITEELEMENT_HH

#include <algorithm>
#include <vector>

#include <dune/common/version.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>

namespace Dune
{
  template<class LocalBasis>
  class ScaledLocalBasis
  {
  public:

    using Traits = typename LocalBasis::Traits;

    ScaledLocalBasis () = delete;

    ScaledLocalBasis (const std::vector<double>& scalingWeights) :
      wrappedBasis_(nullptr),
      scalingWeights(scalingWeights) {}

    //! \brief number of shape functions
    unsigned int size () const
    {
      return wrappedBasis_->size();
    }

    inline void evaluateFunction (const typename Traits::DomainType& x,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      wrappedBasis_->evaluateFunction(x, out);
      auto weight = scalingWeights.begin();
      std::for_each(std::begin(out), std::end(out),
          [&](typename Traits::RangeType& o)
          {
            o *= *weight;
            ++weight;
          });
    }

    inline void
    evaluateJacobian (const typename Traits::DomainType& x,
                      std::vector<typename Traits::JacobianType>& out) const
    {
      wrappedBasis_->evaluateJacobian(x, out);
      auto weight = scalingWeights.begin();
      std::for_each(std::begin(out), std::end(out),
          [&](typename Traits::JacobianType& o)
          {
            o *= *weight;
            ++weight;
          });
    }

    void partial(const std::array<unsigned int,2>& order,
                 const typename Traits::DomainType& in,
                 std::vector<typename Traits::RangeType>& out) const
    {
      wrappedBasis_->partial(order, in, out);
      auto weight = scalingWeights.begin();
      std::for_each(std::begin(out), std::end(out),
          [&](typename Traits::RangeType& o)
          {
            o *= *weight;
            ++weight;
          });
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return wrappedBasis_->order();
    }

    void setWrappedBasis (const LocalBasis& wrappedBasis)
    {
      wrappedBasis_ = &wrappedBasis;
    }

  private:
    const LocalBasis* wrappedBasis_;
    const std::vector<double>& scalingWeights;
  };

  template<class LocalInterpolation>
  class ScaledLocalInterpolation
  {
  public:

    ScaledLocalInterpolation() = delete;
    ScaledLocalInterpolation(const std::vector<double>& scalingWeights) :
      wrappedInterpolation_(nullptr),
      scalingWeights(scalingWeights) {};

    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      wrappedInterpolation_->interpolate(f, out);
      auto weight = scalingWeights.begin();
      std::for_each(std::begin(out), std::end(out),
          [&](C& o)
          {
            o /= *weight;
            ++weight;
          });
    }

    void setWrappedInterpolation(const LocalInterpolation& wrappedInterpolation)
    {
      wrappedInterpolation_ = &wrappedInterpolation;
    }

  private:
    const LocalInterpolation* wrappedInterpolation_;
    const std::vector<double>& scalingWeights;
  };

  template<class LocalFiniteElement>
  class ScaledLocalFiniteElement
  {
  public:
    using Traits = LocalFiniteElementTraits<
        ScaledLocalBasis<typename LocalFiniteElement::Traits::LocalBasisType>,
        typename LocalFiniteElement::Traits::LocalCoefficientsType,
        ScaledLocalInterpolation<
          typename LocalFiniteElement::Traits::LocalInterpolationType>>;

    ScaledLocalFiniteElement () = delete;

    // TODO: replace double with more appropriate type
    ScaledLocalFiniteElement (const std::vector<double>& scalingWeights) :
      wrappedFiniteElement_(nullptr),
      scalingWeights_(scalingWeights),
      localBasis_{scalingWeights_},
      localInterpolation_{scalingWeights_}
    {}

    void setWrappedFiniteElement (const LocalFiniteElement& newFiniteElement)
    {
      wrappedFiniteElement_ = &newFiniteElement;
      localBasis_.setWrappedBasis(wrappedFiniteElement_->localBasis());
      localInterpolation_.setWrappedInterpolation(
          wrappedFiniteElement_->localInterpolation());
    }

    const typename Traits::LocalBasisType& localBasis () const
    {
      return localBasis_;
    }

    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return wrappedFiniteElement_->localCoefficients();
    }

    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return localInterpolation_;
    }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return wrappedFiniteElement_->size();
    }

    GeometryType type () const
    {
      return wrappedFiniteElement_->type();
    }

  private:
    const LocalFiniteElement* wrappedFiniteElement_;
    const std::vector<double>& scalingWeights_;
    typename Traits::LocalBasisType localBasis_;
    typename Traits::LocalInterpolationType localInterpolation_;
  };

}

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NORMALIZEDBASISADAPTOR_FINITEELEMENT_HH
