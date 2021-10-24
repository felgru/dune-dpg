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

template<typename GV, int level, int k, typename R=double>
class LagrangeDGRefinedDGNode;

#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,8)
template<typename GV, int level, int k, class MI, typename R=double>
class LagrangeDGRefinedDGNodeIndexSet;
#endif


#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
template<typename GV, int level, int k, typename R=double>
#else
template<typename GV, int level, int k, class MI, typename R=double>
#endif
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


  using Node = LagrangeDGRefinedDGNode<GV, level, k, R>;

#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,8)
  using IndexSet = LagrangeDGRefinedDGNodeIndexSet<GV, level, k, MI, R>;
#endif

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
  static constexpr size_type maxMultiIndexSize = 1;
  static constexpr size_type minMultiIndexSize = 1;
  static constexpr size_type multiIndexBufferSize = 1;
#else
  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 1>;
#endif

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

  Node makeNode() const
  {
    return Node{};
  }

#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,8)
  IndexSet makeIndexSet() const
  {
    return IndexSet{*this};
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
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
  template<class SizePrefix>
  size_type size(const SizePrefix& prefix) const
#else
  size_type size(const SizePrefix prefix) const
#endif
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

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  template<typename It>
  It indices(const Node& node, It it) const
  {
    const auto& gridIndexSet = gridView().indexSet();
    const auto& element = node.element();

    for (size_type i = 0, end = node.size(); i < end; ++it, ++i)
    {
      if constexpr (dim == 1) {
        *it = {{ dofsPerEdge * gridIndexSet.subIndex(element,0,0) + i }};
        continue;
      } else if constexpr (dim == 2) {
        if (element.type().isTriangle())
        {
          *it = {{ dofsPerTriangle * gridIndexSet.subIndex(element,0,0) + i }};
          continue;
        }
        else if (element.type().isQuadrilateral())
        {
          *it = {{ quadrilateralOffset_
                   + dofsPerQuad * gridIndexSet.subIndex(element,0,0) + i }};
          continue;
        }
        else
          DUNE_THROW(Dune::NotImplemented,
              "2d elements have to be triangles or quadrilaterals");
      } else {
        static_assert(dim >= 1 && dim <= 2,
                      "The indices method is not yet implemented for grids of this dimension!");
      }
    }
    return it;
  }

//protected:
  GridView gridView_;

  size_type quadrilateralOffset_;
};



template<typename GV, int level, int k, typename R>
class LagrangeDGRefinedDGNode :
  public LeafBasisNode,
  public RefinedNode < typename GV::template Codim<0>::Entity
                     , typename GV::ctype, GV::dimension, level>
{
  static constexpr int dim = GV::dimension;

  using RefinedNodeBase =
          RefinedNode < typename GV::template Codim<0>::Entity
                      , typename GV::ctype, dim, level>;
  using FiniteElementCache = typename Dune::PQkLocalFiniteElementCache<typename GV::ctype, R, dim, k>;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using SubElement = typename RefinedNodeBase::SubElement;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  LagrangeDGRefinedDGNode() :
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



#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,8)
template<typename GV, int level, int k, class MI, typename R>
class LagrangeDGRefinedDGNodeIndexSet
{
  static constexpr int dim = GV::dimension;

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = LagrangeDGRefinedDGPreBasis<GV, level, k, MI, R>;

  using Node = typename PreBasis::Node;

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
    return preBasis_->indices(*node_, it);
  }

protected:
  const PreBasis* preBasis_;

  const Node* node_;
};
#endif



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Basis of a scalar k-th-order Lagrangean-DG finite element space
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 * \tparam R The range type of the local basis
 */
template<typename GV, int level, int k, typename R=double>
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
using LagrangeDGRefinedDGBasis = RefinedGlobalBasis<LagrangeDGRefinedDGPreBasis<GV, level, k, R>>;
#else
using LagrangeDGRefinedDGBasis = RefinedGlobalBasis<LagrangeDGRefinedDGPreBasis<GV, level, k, FlatMultiIndex<std::size_t>, R> >;
#endif



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGREFINEDDGBASIS_HH
