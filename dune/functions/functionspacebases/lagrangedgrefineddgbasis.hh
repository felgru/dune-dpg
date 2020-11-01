// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGREFINEDDGBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGREFINEDDGBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/power.hh>
#include <dune/common/version.hh>

#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/refinedglobalbasis.hh>




namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   LagrangeDGRefinedDGPreBasis
//   LagrangeDGRefinedDGNodeIndexSet
//   LagrangeDGRefinedDGNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
template<typename GV, int level, int k>
class LagrangeDGRefinedDGNode;

template<typename GV, int level, int k, class MI>
class LagrangeDGRefinedDGNodeIndexSet;
#else
template<typename GV, int level, int k, typename TP>
class LagrangeDGRefinedDGNode;

template<typename GV, int level, int k, class MI, class TP>
class LagrangeDGRefinedDGNodeIndexSet;
#endif


template<typename GV, int level, int k, class MI>
class LagrangeDGRefinedDGPreBasis
  : public DGRefinedPreBasisConstants<GV::dimension, level, k>
{
  static constexpr int dim = GV::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = std::size_t;

  using RefinementConstants = DGRefinedPreBasisConstants<dim, level, k>;

  // Precompute the number of dofs per entity type
  constexpr static int dofsPerEdge
      = RefinementConstants::numberOfSubEdges
      * RefinementConstants::dofsPerSubEdge;
  constexpr static int dofsPerTriangle
      = RefinementConstants::numberOfSubTriangles
      * RefinementConstants::dofsPerSubTriangle;
  constexpr static int dofsPerQuad
      = RefinementConstants::numberOfSubQuads
      * RefinementConstants::dofsPerSubQuad;


#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
  using Node = LagrangeDGRefinedDGNode<GV, level, k>;

  using IndexSet = LagrangeDGRefinedDGNodeIndexSet<GV, level, k, MI>;
#else
  template<class TP>
  using Node = LagrangeDGRefinedDGNode<GV, level, k, TP>;

  template<class TP>
  using IndexSet = LagrangeDGRefinedDGNodeIndexSet<GV, level, k, MI, TP>;
#endif

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  /** \brief Constructor for a given grid view object */
  LagrangeDGRefinedDGPreBasis(const GridView& gv) :
    gridView_(gv)
  {}


  void initializeIndices()
  {
    if constexpr (dim == 1) {
      return;
    } else if constexpr (dim == 2) {
      quadrilateralOffset_ = dofsPerTriangle
                             * gridView_.size(GeometryTypes::triangle);
    } else {
      static_assert(dim >= 1 && dim <= 2,
                      "LagrangeDGRefinedDGPreBasis not implmented for grids of this dimension.");
    }
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return gridView_;
  }

  void update (const GridView& gv)
  {
    gridView_ = gv;
  }

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
  Node makeNode() const
  {
    return Node{};
  }

  IndexSet makeIndexSet() const
  {
    return IndexSet{*this};
  }
#else
  template<class TP>
  Node<TP> node(const TP& tp) const
  {
    return Node<TP>{tp};
  }

  template<class TP>
  IndexSet<TP> indexSet() const
  {
    return IndexSet<TP>{*this};
  }
#endif

  size_type size() const
  {
    if constexpr (dim == 1) {
      return dofsPerEdge * gridView_.size(0);
    } else if constexpr (dim == 2) {
      return dofsPerTriangle * gridView_.size(GeometryTypes::triangle)
           + dofsPerQuad * gridView_.size(GeometryTypes::quadrilateral);
    } else {
      static_assert(dim >= 1 && dim <= 2,
                    "No size method implemented for this dimension!");
    }
  }

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix prefix) const
  {
    assert(prefix.size() == 0 || prefix.size() == 1);
    return (prefix.size() == 0) ? size() : 0;
  }

  /** \todo This method has been added to the interface without prior discussion. */
  size_type dimension() const
  {
    return size();
  }

  size_type maxNodeSize() const
  {
    return StaticPower<4,level>::power*StaticPower<(k+1),GV::dimension>::power;
  }

//protected:
  GridView gridView_;

  size_type quadrilateralOffset_;
};



#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
template<typename GV, int level, int k>
class LagrangeDGRefinedDGNode :
  public LeafBasisNode,
#else
template<typename GV, int level, int k, typename TP>
class LagrangeDGRefinedDGNode :
  public LeafBasisNode<std::size_t, TP>,
#endif
  public RefinedNode < typename GV::template Codim<0>::Entity
                     , typename GV::ctype, GV::dimension, level>
{
  static constexpr int dim = GV::dimension;

#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,7)
  using Base = LeafBasisNode<std::size_t, TP>;
#endif
  using RefinedNodeBase =
          RefinedNode < typename GV::template Codim<0>::Entity
                      , typename GV::ctype, dim, level>;
  using FiniteElementCache = typename Dune::PQkLocalFiniteElementCache<typename GV::ctype, double, dim, k>;

public:

  using size_type = std::size_t;
#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,7)
  using TreePath = TP;
#endif
  using Element = typename GV::template Codim<0>::Entity;
  using SubElement = typename RefinedNodeBase::SubElement;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
  LagrangeDGRefinedDGNode() :
#else
  LagrangeDGRefinedDGNode(const TreePath& treePath) :
    Base(treePath),
#endif
    RefinedNodeBase(),
    finiteElement_(nullptr)
  {}

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement() const
  {
    return *finiteElement_;
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    // TODO: We assume that all subElements are of the same geometry type
    //       which is the same as for e.
    this->element_ = &e;
    finiteElement_ = &(feCache_.get(this->element_->type()));
    using PreBasis = LagrangeDGRefinedDGPreBasis<GV, level, k, void>;
    size_type numberOfSubElements;
    if(e.type().isTriangle()) {
      numberOfSubElements = PreBasis::numberOfSubTriangles;
    } else if(e.type().isQuadrilateral()) {
      numberOfSubElements = PreBasis::numberOfSubQuads;
    } else {
      DUNE_THROW(Dune::NotImplemented,
                 "LagrangeDGRefinedNode::bind() not implemented for element type "
                 << e.type().id());
    }
    this->setSize(numberOfSubElements*finiteElement_->size());
  }

  void bindSubElement(const SubElement& se)
  {
    this->subElement_ = &se;
    assert(this->element_->type() == this->subElement_->type());
  }

protected:

  FiniteElementCache feCache_;
  const FiniteElement* finiteElement_;
};



#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
template<typename GV, int level, int k, class MI>
#else
template<typename GV, int level, int k, class MI, class TP>
#endif
class LagrangeDGRefinedDGNodeIndexSet
{
  static constexpr int dim = GV::dimension;

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = LagrangeDGRefinedDGPreBasis<GV, level, k, MI>;

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
  using Node = typename PreBasis::Node;
#else
  using Node = typename PreBasis::template Node<TP>;
#endif

  LagrangeDGRefinedDGNodeIndexSet(const PreBasis& preBasis) :
    preBasis_(&preBasis)
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Node& node)
  {
    node_ = &node;
  }

  /** \brief Unbind the view
   */
  void unbind()
  {
    node_ = nullptr;
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const
  {
    assert(node_ != nullptr);
    return node_->size();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  template<typename It>
  It indices(It it) const
  {
    assert(node_ != nullptr);
    const auto& gridIndexSet = preBasis_->gridView().indexSet();
    const auto& element = node_->element();

    for (size_type i = 0, end = this->size(); i < end; ++it, ++i)
    {
      if constexpr (dim == 1) {
        *it = {{ preBasis_->dofsPerEdge
                 * gridIndexSet.subIndex(element,0,0) + i }};
        continue;
      } else if constexpr (dim == 2) {
        if (element.type().isTriangle())
        {
          *it = {{ preBasis_->dofsPerTriangle
                   * gridIndexSet.subIndex(element,0,0) + i }};
          continue;
        }
        else if (element.type().isQuadrilateral())
        {
          *it = {{ preBasis_->quadrilateralOffset_
                   + preBasis_->dofsPerQuad
                     * gridIndexSet.subIndex(element,0,0) + i }};
          continue;
        }
        else
          DUNE_THROW(Dune::NotImplemented,
              "2d elements have to be triangles or quadrilaterals");
      } else {
        static_assert(dim >= 1 && dim <= 2,
                      "The index method is not yet implemented for grids of this dimension!");
      }
    }
    return it;
  }

protected:
  const PreBasis* preBasis_;

  const Node* node_;
};



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Basis of a scalar k-th-order Lagrangean-DG finite element space
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int level, int k>
using LagrangeDGRefinedDGBasis = RefinedGlobalBasis<LagrangeDGRefinedDGPreBasis<GV, level, k, FlatMultiIndex<std::size_t> > >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGREFINEDDGBASIS_HH
