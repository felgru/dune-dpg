// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NORMALIZEDREFINEDBASISADAPTOR_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NORMALIZEDREFINEDBASISADAPTOR_HH

#include <array>
#include <cmath>

#include <dune/dpg/assemble_helper.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/refinednode.hh>

#include <dune/istl/matrix.hh>

#include <dune/typetree/leafnode.hh>

#include "normalizedbasisadaptor/localfiniteelement.hh"


namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   NormalizedRefinedPreBasis
//   NormalizedRefinedNodeIndexSet
//   NormalizedRefinedNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename InnerProduct, typename TP>
class NormalizedRefinedNode;

template<typename Basis, class TP>
class NormalizedRefinedNodeIndexSet;


template<typename InnerProduct>
class NormalizedRefinedPreBasis
{
  using Basis = std::tuple_element_t<0, typename InnerProduct::TestSpaces>;
  using WrappedPreBasis = typename Basis::PreBasis;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = typename Basis::GridView;
  using size_type = std::size_t;

  using WrappedBasis = Basis;

  template<class TP>
  using Node = NormalizedRefinedNode<InnerProduct, TP>;

  template<class TP>
  using IndexSet = NormalizedRefinedNodeIndexSet<
                      NormalizedRefinedPreBasis<InnerProduct>, TP>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = typename Basis::MultiIndex;

  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  /** \brief Constructor for a given grid view object */
  NormalizedRefinedPreBasis(const InnerProduct& ip) :
    wrappedPreBasis_(std::get<0>(*ip.getTestSpaces()).preBasis()),
    innerProduct_(ip)
  {}


  void initializeIndices()
  {
    wrappedPreBasis_.initializeIndices();
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return wrappedPreBasis_.gridView();
  }

  void update(const GridView& gv)
  {
    wrappedPreBasis_.update(gv);
    Dune::detail::updateSpaces(*innerProduct_.getTestSpaces(), gv);
  }

  template<class TP>
  Node<TP> node(const TP& tp) const
  {
    return Node<TP>{innerProduct_, tp};
  }

  template<class TP>
  IndexSet<TP> indexSet() const
  {
    return IndexSet<TP>{*this};
  }

  size_type size() const
  {
    return wrappedPreBasis_.size();
  }

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix prefix) const
  {
    return wrappedPreBasis_.size(prefix);
  }

  /** \todo This method has been added to the interface without prior discussion. */
  size_type dimension() const
  {
    return size();
  }

  size_type maxNodeSize() const
  {
    return wrappedPreBasis_.maxNodeSize();
  }

  const WrappedPreBasis& wrappedPreBasis() const
  {
    return wrappedPreBasis_;
  }

  const InnerProduct& innerProduct() const
  {
    return innerProduct_;
  }

protected:
  WrappedPreBasis wrappedPreBasis_;
  InnerProduct innerProduct_;
};



template<typename InnerProduct, typename TP>
class NormalizedRefinedNode :
  public LeafBasisNode<std::size_t, TP>
{
  using Basis = std::tuple_element_t<0, typename InnerProduct::TestSpaces>;
  using LocalViews
      = Dune::detail::getLocalViews_t<typename InnerProduct::TestSpaces>;

  using Base = LeafBasisNode<std::size_t, TP>;
  using WrappedPreBasis = typename Basis::PreBasis;
  using WrappedNode = typename WrappedPreBasis::template Node<TP>;
  using PreBasis = NormalizedRefinedPreBasis<InnerProduct>;

public:

  using size_type = std::size_t;
  using TreePath = TP;
  using Element = typename WrappedNode::Element;
  using FiniteElement
      = ScaledLocalFiniteElement<typename WrappedNode::FiniteElement>;
  using RefinementGrid = typename WrappedNode::RefinementGrid;
  using RefinementGridView = typename WrappedNode::RefinementGridView;

  NormalizedRefinedNode(const PreBasis& preBasis,
                        const TreePath& treePath) :
    wrappedNode_(treePath),
    innerProduct_(preBasis.innerProduct()),
    localViews_(Dune::detail::getLocalViews(*innerProduct_.getTestSpaces())),
    scalingWeights_(),
    finiteElement_()
  {
    scalingWeights_.reserve(preBasis.maxNodeSize());
  }

  const Element& element() const
  {
    return wrappedNode_.element();
  }

  const FiniteElement& finiteElement() const
  {
    return finiteElement_;
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    wrappedNode_.bind(e);
    this->setSize(wrappedNode_.size());
    Dune::detail::bindLocalViews(localViews_, e);
    innerProduct_.bind(localViews_);

    Matrix<FieldMatrix<double,1,1> > localGramian;
    innerProduct_.getLocalMatrix(localGramian);
    scalingWeights_.resize(localGramian.N());
    for(size_t i = 0, size = scalingWeights_.size(); i < size; i++)
      scalingWeights_[i] = 1./std::sqrt(localGramian[i][i]);

    finiteElement_.setWrappedFiniteElementAndWeights(
        wrappedNode_.finiteElement(), scalingWeights_.begin());
  }

  const RefinementGridView refinedReferenceElementGridView() const
  {
    return wrappedNode_.refinedReferenceElementGridView();
  }

  const WrappedNode& wrappedBasisNode() const
  {
    return wrappedNode_;
  }

protected:

  WrappedNode wrappedNode_;
  InnerProduct innerProduct_;
  LocalViews localViews_;
  std::vector<double> scalingWeights_;
  FiniteElement finiteElement_;
};



template<typename PB, class TP>
class NormalizedRefinedNodeIndexSet
{
  using Basis = typename PB::WrappedBasis;

  using WrappedIndexSet = typename Basis::PreBasis::template IndexSet<TP>;

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = typename Basis::MultiIndex;

  using PreBasis = PB;

  using Node = typename PreBasis::template Node<TP>;

  NormalizedRefinedNodeIndexSet(const PreBasis& preBasis) :
    wrappedIndexSet_(preBasis.wrappedPreBasis())
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Node& node)
  {
    wrappedIndexSet_.bind(node.wrappedBasisNode());
  }

  /** \brief Unbind the view
   */
  void unbind()
  {
    wrappedIndexSet_.unbind();
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const
  {
    return wrappedIndexSet_.size();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  template<typename It>
  It indices(It it) const
  {
    return wrappedIndexSet_.indices(it);
  }

protected:
  WrappedIndexSet wrappedIndexSet_;
};



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Basis adaptor to normalize a finite element space
 *         wrt given inner product
 */
template<typename InnerProduct>
using NormalizedRefinedBasis
    = DefaultGlobalBasis<NormalizedRefinedPreBasis<InnerProduct>>;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NORMALIZEDREFINEDBASISADAPTOR_HH
