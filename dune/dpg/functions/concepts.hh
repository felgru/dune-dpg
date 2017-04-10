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


template<class NodeFactory>
struct ConstrainedNodeIndexSet
{
  template<class I>
  auto require(const I& indexSet) -> decltype(
    requireType<typename I::size_type>(),
    requireType<typename I::MultiIndex>(),
    requireType<typename I::NodeFactory>(),
    requireType<typename I::Node>(),
    requireSameType<typename I::NodeFactory, NodeFactory>(),
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

template<class LocalView>
struct ConstrainedLocalIndexSet
{
  template<class I>
  auto require(const I& indexSet) -> decltype(
    requireType<typename I::size_type>(),
    requireType<typename I::MultiIndex>(),
    requireType<typename I::LocalView>(),
    requireSameType<typename I::LocalView, LocalView>(),
    requireConvertible<typename I::LocalView>(indexSet.localView()),
    const_cast<I&>(indexSet).bind(std::declval<typename I::LocalView>()),
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


// Concept for a NodeFactory
template<class GridView>
struct ConstrainedNodeFactory
{
  using RootTreePath = decltype(TypeTree::hybridTreePath());

  template<class F>
  auto require(const F& factory) -> decltype(
    requireType<typename F::GridView>(),
    requireType<typename F::size_type>(),
    requireType<typename F::MultiIndex>(),
    requireType<typename F::SizePrefix>(),
    requireType<typename F::template Node<RootTreePath>>(),
    requireType<typename F::template IndexSet<RootTreePath>>(),
    requireSameType<typename F::GridView, GridView>(),
    const_cast<F&>(factory).initializeIndices(),
    requireConvertible<typename F::GridView>(factory.gridView()),
    requireConvertible<typename F::template Node<RootTreePath>>(factory.node(RootTreePath())),
    requireConvertible<typename F::template IndexSet<RootTreePath>>(factory.template indexSet<RootTreePath>()),
    requireConvertible<typename F::size_type>(factory.size()),
    requireConvertible<typename F::size_type>(factory.size(std::declval<typename F::SizePrefix>())),
    requireConvertible<typename F::size_type>(factory.dimension()),
    requireConvertible<typename F::size_type>(factory.maxNodeSize()),
    requireConcept<BasisTree<typename F::GridView>>(factory.node(RootTreePath())),
    requireConcept<ConstrainedNodeIndexSet<F>>(factory.template indexSet<RootTreePath>())
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
    requireType<typename B::LocalIndexSet>(),
    requireType<typename B::LocalView>(),
    requireSameType<typename B::GridView, GridView>(),
    requireConvertible<typename B::GridView>(basis.gridView()),
    requireConvertible<typename B::LocalIndexSet>(basis.localIndexSet()),
    requireConvertible<typename B::LocalView>(basis.localView()),
    requireConvertible<typename B::size_type>(basis.size()),
    requireConvertible<typename B::size_type>(basis.size(std::declval<typename B::SizePrefix>())),
    requireConvertible<typename B::size_type>(basis.dimension()),
    requireConcept<ConstrainedLocalIndexSet<typename B::LocalView>>
        (basis.localIndexSet()),
    requireConcept<LocalView<B>>(basis.localView())
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
