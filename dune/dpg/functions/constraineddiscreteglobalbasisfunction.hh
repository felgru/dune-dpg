// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_DPG_CONSTRAINEDDISCRETEGLOBALBASISFUNCTIONS_HH
#define DUNE_FUNCTIONS_DPG_CONSTRAINEDDISCRETEGLOBALBASISFUNCTIONS_HH

#include <numeric>
#include <type_traits>
#include <vector>

#include <dune/dpg/functions/localindexsetiteration.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>
#include <dune/functions/gridfunctions/gridfunction.hh>

namespace Dune {
namespace Functions {



template<typename B, typename V,
         typename R = typename V::value_type>
class ConstrainedDiscreteGlobalBasisFunction
{
public:
  using Basis = B;
  using Vector = V;

  using GridView = typename Basis::GridView;
  using EntitySet = GridViewEntitySet<GridView, 0>;
  using Node = typename Basis::LocalView::Tree;

  using Domain = typename EntitySet::GlobalCoordinate;
  using Range = R;

  using LocalDomain = typename EntitySet::LocalCoordinate;
  using Element = typename EntitySet::Element;

  using Traits = Imp::GridFunctionTraits<Range(Domain), EntitySet,
                                         DefaultDerivativeTraits, 16>;

  class LocalFunction
  {
    using LocalBasisView = typename Basis::LocalView;
    using size_type = typename Basis::size_type;

    using LocalBasisRange = typename
        Node::FiniteElement::Traits::LocalBasisType::Traits::RangeType;

    using NodeData = std::vector<LocalBasisRange>;

  public:

    using GlobalFunction = ConstrainedDiscreteGlobalBasisFunction;
    using Domain = LocalDomain;
    using Range = GlobalFunction::Range;
    using Element = GlobalFunction::Element;

    LocalFunction(const ConstrainedDiscreteGlobalBasisFunction& globalFunction)
      : globalFunction_(globalFunction)
      , localBasisView_(globalFunction.basis().localView())
      , node_(&localBasisView_.tree())
    {
      shapeFunctionValues_.reserve(localBasisView_.maxSize());
      localDoFs_.reserve(localBasisView_.maxSize());
    }

    LocalFunction(const LocalFunction& other)
      : LocalFunction(other.globalFunction_) {}

    LocalFunction operator=(const LocalFunction& other)
    {
      globalFunction_ = other.globalFunction_;
      localBasisView_ = other.localBasisView_;
      node_ = &localBasisView_.tree();

      shapeFunctionValues_.reserve(localBasisView_.maxSize());
      localDoFs_.reserve(localBasisView_.maxSize());
    }

    /**
     * \brief Bind LocalFunction to grid element.
     *
     * You must call this method before evaluation
     * and after changes to the coefficient vector.
     */
    void bind(const Element& element)
    {
      localBasisView_.bind(element);
      localDoFs_.resize(localBasisView_.size());
      copyToLocalVector(globalFunction_.dofs(), localDoFs_, localBasisView_);
    }

    void unbind()
    {
      localBasisView_.unbind();
    }

    /**
     * \brief Evaluate LocalFunction at bound element.
     *
     * The result of this method is undefined if you did
     * not call bind() beforehand or changed the coefficient
     * vector after the last call to bind(). In the latter case
     * you have to call bind() again in order to make operator()
     * usable.
     */
    Range operator()(const Domain& x) const
    {
      auto&& fe = node_->finiteElement();
      auto&& localBasis = fe.localBasis();

      shapeFunctionValues_.resize(localBasis.size());
      localBasis.evaluateFunction(x, shapeFunctionValues_);

      Range y = std::inner_product(shapeFunctionValues_.cbegin(),
                                   shapeFunctionValues_.cend(),
                                   localDoFs_.cbegin(), Range(0));
      return y;
    }

    const Element& localContext() const
    {
      return localBasisView_.element();
    }

    friend typename Traits::LocalFunctionTraits::DerivativeInterface derivative(const LocalFunction& t)
    {
      DUNE_THROW(NotImplemented,"not implemented");
    }

  private:

    const ConstrainedDiscreteGlobalBasisFunction& globalFunction_;
    LocalBasisView localBasisView_;

    mutable std::vector<typename V::value_type> shapeFunctionValues_;
    std::vector<typename V::value_type> localDoFs_;
    const Node* node_;
  };

  ConstrainedDiscreteGlobalBasisFunction(const Basis & basis, const V & coefficients) :
    entitySet_(basis.gridView()),
    basis_(basis),
    coefficients_(coefficients)
  {}

  const Basis& basis() const
  {
    return basis_;
  }

  const V& dofs() const
  {
    return coefficients_;
  }

  // TODO: Implement this using hierarchic search
  Range operator() (const Domain& x) const
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  friend typename Traits::DerivativeInterface derivative(const ConstrainedDiscreteGlobalBasisFunction& t)
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  friend LocalFunction localFunction(const ConstrainedDiscreteGlobalBasisFunction& t)
  {
    return LocalFunction(t);
  }

  /**
   * \brief Get associated EntitySet
   */
  const EntitySet& entitySet() const
  {
    return entitySet_;
  }

private:

  EntitySet entitySet_;
  const Basis& basis_;
  const V& coefficients_;
};



template<typename R, typename B, typename V>
auto makeConstrainedDiscreteGlobalBasisFunction(const B& basis, const V& vector)
{
  using Basis = std::decay_t<B>;
  using Vector = std::decay_t<V>;
  return ConstrainedDiscreteGlobalBasisFunction<Basis, Vector, R>
      (basis, vector);
}



} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_DPG_CONSTRAINEDDISCRETEGLOBALBASISFUNCTIONS_HH
