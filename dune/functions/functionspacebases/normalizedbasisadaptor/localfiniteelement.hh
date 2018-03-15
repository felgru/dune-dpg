// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NORMALIZEDBASISADAPTOR_FINITEELEMENT_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NORMALIZEDBASISADAPTOR_FINITEELEMENT_HH

#include <algorithm>
#include <vector>

#include <dune/common/version.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>

namespace Dune
{
  template<class LocalBasis>
  class ScaledLocalBasis
  {
  public:

    using Traits = typename LocalBasis::Traits;

    ScaledLocalBasis () :
      wrappedBasis_(nullptr),
      weightsBegin_() {}

    //! \brief number of shape functions
    unsigned int size () const
    {
      return wrappedBasis_->size();
    }

    inline void evaluateFunction (const typename Traits::DomainType& x,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      wrappedBasis_->evaluateFunction(x, out);
      auto weight = weightsBegin_;
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
      auto weight = weightsBegin_;
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
      auto weight = weightsBegin_;
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

    void setWrappedBasisAndWeights(
        const LocalBasis& wrappedBasis,
        typename std::vector<double>::const_iterator weightsBegin)
    {
      wrappedBasis_ = &wrappedBasis;
      weightsBegin_ = std::move(weightsBegin);
    }

  private:
    const LocalBasis* wrappedBasis_;
    typename std::vector<double>::const_iterator weightsBegin_;
  };

  template<class LocalInterpolation>
  class ScaledLocalInterpolation
  {
  public:

    ScaledLocalInterpolation() :
      wrappedInterpolation_(nullptr),
      weightsBegin_() {};

    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      wrappedInterpolation_->interpolate(f, out);
      auto weight = weightsBegin_;
      std::for_each(std::begin(out), std::end(out),
          [&](C& o)
          {
            o /= *weight;
            ++weight;
          });
    }

    void setWrappedInterpolationAndWeights(
        const LocalInterpolation& wrappedInterpolation,
        typename std::vector<double>::const_iterator weightsBegin)
    {
      wrappedInterpolation_ = &wrappedInterpolation;
      weightsBegin_ = std::move(weightsBegin);
    }

  private:
    const LocalInterpolation* wrappedInterpolation_;
    typename std::vector<double>::const_iterator weightsBegin_;
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

    ScaledLocalFiniteElement () :
      wrappedFiniteElement_(nullptr),
      localBasis_(),
      localInterpolation_()
    {}

    // TODO: replace double with more appropriate type
    void setWrappedFiniteElementAndWeights(
        const LocalFiniteElement& newFiniteElement,
        const typename std::vector<double>::const_iterator weightsBegin)
    {
      wrappedFiniteElement_ = &newFiniteElement;
      localBasis_.setWrappedBasisAndWeights
          (wrappedFiniteElement_->localBasis(), weightsBegin);
      localInterpolation_.setWrappedInterpolationAndWeights(
          wrappedFiniteElement_->localInterpolation(), weightsBegin);
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
    typename Traits::LocalBasisType localBasis_;
    typename Traits::LocalInterpolationType localInterpolation_;
  };

}

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NORMALIZEDBASISADAPTOR_FINITEELEMENT_HH
