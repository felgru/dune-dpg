// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_FUNCTIONS_FUNCTIONSPACEBASES_CONCEPTS_HH
#define DUNE_DPG_FUNCTIONS_FUNCTIONSPACEBASES_CONCEPTS_HH


#include <dune/common/concept.hh>

#include <dune/functions/common/utility.hh>

#include <dune/functions/functionspacebases/nodes.hh>

#include <dune/functions/functionspacebases/concepts.hh>


namespace Dune {
namespace Functions {
namespace Concept {

using namespace Dune::Concept;


template<class PreBasis>
struct ConstrainedNodeIndexSet
{
  template<class I>
  auto require(const I& indexSet) -> decltype(
    requireType<typename I::size_type>(),
    requireType<typename I::MultiIndex>(),
    requireType<typename I::PreBasis>(),
    requireType<typename I::Node>(),
    requireSameType<typename I::PreBasis, PreBasis>(),
    const_cast<I&>(indexSet).bind(std::declval<typename I::Node>()),
    const_cast<I&>(indexSet).unbind(),
    requireConvertible<typename I::size_type>(indexSet.size()),
    requireConvertible<const std::vector<typename I::MultiIndex>&>
        (indexSet.indicesLocalGlobal()),
    requireConvertible<typename I::size_type>(indexSet.constraintsSize()),
    requireConvertible<typename I::size_type>(
        indexSet.constraintOffset(std::declval<typename I::size_type>())),
    requireConvertible<const typename I::ConstraintWeights&>(
        indexSet.constraintWeights(std::declval<typename I::size_type>()))
  );
};

template<class GridView>
struct ConstrainedPreBasis
{
  template<class F>
  auto require(const F& factory) -> decltype(
    requireType<typename F::GridView>(),
    requireType<typename F::size_type>(),
    requireType<typename F::MultiIndex>(),
    requireType<typename F::SizePrefix>(),
    requireType<typename F::Node>(),
    requireType<typename F::IndexSet>(),
    requireSameType<typename F::GridView, GridView>(),
    const_cast<F&>(factory).initializeIndices(),
    requireConvertible<typename F::GridView>(factory.gridView()),
    requireConvertible<typename F::Node>(factory.makeNode()),
    requireConvertible<typename F::IndexSet>(factory.makeIndexSet()),
    requireConvertible<typename F::size_type>(factory.size()),
    requireConvertible<typename F::size_type>(factory.size(std::declval<typename F::SizePrefix>())),
    requireConvertible<typename F::size_type>(factory.dimension()),
    requireConvertible<typename F::size_type>(factory.maxNodeSize()),
    requireConcept<BasisTree<typename F::GridView>>(factory.makeNode()),
    requireConcept<ConstrainedNodeIndexSet<F>>(factory.makeIndexSet())
  );
};

template<class GlobalBasis>
struct ConstrainedLocalView
{
  template<class V>
  auto require(const V& localView) -> decltype(
    requireType<typename V::size_type>(),
    requireType<typename V::MultiIndex>(),
    requireType<typename V::GlobalBasis>(),
    requireType<typename V::Tree>(),
    requireType<typename V::GridView>(),
    requireType<typename V::Element>(),
    requireSameType<typename V::GlobalBasis, GlobalBasis>(),
    requireSameType<typename V::GridView, typename GlobalBasis::GridView>(),
    requireSameType<typename V::size_type, typename GlobalBasis::size_type>(),
    requireSameType<typename V::Element, typename GlobalBasis::GridView::template Codim<0>::Entity>(),
    const_cast<V&>(localView).bind(std::declval<typename V::Element>()),
    const_cast<V&>(localView).unbind(),
    requireConvertible<typename V::Tree>(localView.tree()),
    requireConvertible<typename V::size_type>(localView.size()),
    requireConvertible<typename V::size_type>(localView.maxSize()),
    requireConvertible<const std::vector<typename V::MultiIndex>&>
        (localView.indicesLocalGlobal()),
    requireConvertible<typename V::size_type>(localView.constraintsSize()),
    requireConvertible<typename V::size_type>(
        localView.constraintOffset(std::declval<typename V::size_type>())),
    requireConvertible<const typename V::ConstraintWeights&>(
        localView.constraintWeights(std::declval<typename V::size_type>())),
    requireConvertible<typename V::GlobalBasis>(localView.globalBasis()),
    requireConcept<BasisTree<typename V::GridView>>(localView.tree()),
    0
  );
};

template<class GridView>
struct ConstrainedGlobalBasis
{
  template<class B>
  auto require(const B& basis) -> decltype(
    requireType<typename B::GridView>(),
    requireType<typename B::size_type>(),
    requireType<typename B::MultiIndex>(),
    requireType<typename B::SizePrefix>(),
    requireType<typename B::LocalView>(),
    requireSameType<typename B::GridView, GridView>(),
    requireConvertible<typename B::GridView>(basis.gridView()),
    requireConvertible<typename B::LocalView>(basis.localView()),
    requireConvertible<typename B::size_type>(basis.size()),
    requireConvertible<typename B::size_type>(basis.size(std::declval<typename B::SizePrefix>())),
    requireConvertible<typename B::size_type>(basis.dimension()),
    requireConcept<ConstrainedLocalView<B>>(basis.localView())
  );
};

template<class GridView>
struct GeneralizedGlobalBasis
{
  template<class B>
  auto require(const B& basis) -> decltype(
    requireTrue<models<GlobalBasis<GridView>, B>() or
                models<ConstrainedGlobalBasis<GridView>, B>()>()
  );
};


} // namespace Dune::Functions::Concept
} // namespace Dune::Functions
} // namespace Dune


#endif // DUNE_DPG_FUNCTIONS_FUNCTIONSPACEBASES_CONCEPTS_HH
