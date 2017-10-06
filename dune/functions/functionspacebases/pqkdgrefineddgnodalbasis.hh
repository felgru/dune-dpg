// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKDGREFINEDDGBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKDGREFINEDDGBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/version.hh>

#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/refinednode.hh>




namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   PQkDGRefinedDGNodeFactory
//   PQkDGRefinedDGNodeIndexSet
//   PQkDGRefinedDGNode
//
// The factory allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, int level, int k, typename TP>
class PQkDGRefinedDGNode;

template<typename GV, int level, int k, class MI, class TP>
class PQkDGRefinedDGNodeIndexSet;


template<typename GV, int level, int k, class MI>
class PQkDGRefinedDGNodeFactory
  : public DGRefinedNodeFactoryConstants<GV::dimension, level, k>
{
  static const int dim = GV::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = std::size_t;

  using RefinementConstants = DGRefinedNodeFactoryConstants<dim, level, k>;

  // Precompute the number of dofs per entity type
  const static int dofsPerEdge     = RefinementConstants::numberOfSubEdges
                                   * RefinementConstants::dofsPerSubEdge;
  const static int dofsPerTriangle = RefinementConstants::numberOfSubTriangles
                                   * RefinementConstants::dofsPerSubTriangle;
  const static int dofsPerQuad     = RefinementConstants::numberOfSubQuads
                                   * RefinementConstants::dofsPerSubQuad;


  template<class TP>
  using Node = PQkDGRefinedDGNode<GV, level, k, TP>;

  template<class TP>
  using IndexSet = PQkDGRefinedDGNodeIndexSet<GV, level, k, MI, TP>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  /** \brief Constructor for a given grid view object */
  PQkDGRefinedDGNodeFactory(const GridView& gv) :
    gridView_(gv)
  {}


  void initializeIndices()
  {
    switch (dim)
    {
      case 1:
      {
        break;
      }
      case 2:
      {
        quadrilateralOffset_ = dofsPerTriangle
                               * gridView_.size(GeometryTypes::triangle);
        break;
      }
      case 3:
      {
        DUNE_THROW(Dune::NotImplemented,
                   "PQkDGRefinedDGNodeFactory not implmented in 3d.");
      }
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

  size_type size() const
  {
    switch (dim)
    {
      case 1:
        return dofsPerEdge * gridView_.size(0);
      case 2:
      {
        return dofsPerTriangle * gridView_.size(GeometryTypes::triangle)
             + dofsPerQuad * gridView_.size(GeometryTypes::quadrilateral);
      }
    }
    DUNE_THROW(Dune::NotImplemented, "No size method for " << dim
                                     << "d grids available yet!");
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



template<typename GV, int level, int k, typename TP>
class PQkDGRefinedDGNode :
  public LeafBasisNode<std::size_t, TP>,
  public RefinedNode < typename GV::template Codim<0>::Entity
                     , typename GV::ctype, GV::dimension, level>
{
  static const int dim = GV::dimension;
  static const int maxSize = StaticPower<(k+1),GV::dimension>::power;

  using Base = LeafBasisNode<std::size_t, TP>;
  using RefinedNodeBase =
          RefinedNode < typename GV::template Codim<0>::Entity
                      , typename GV::ctype, dim, level>;
  using FiniteElementCache = typename Dune::PQkLocalFiniteElementCache<typename GV::ctype, double, dim, k>;

public:

  using size_type = std::size_t;
  using TreePath = TP;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  PQkDGRefinedDGNode(const TreePath& treePath) :
    Base(treePath),
    RefinedNodeBase(),
    finiteElement_(nullptr)
  {}

  //! Return current element, throw if unbound
  const Element& element() const
  {
    return *(this->element_);
  }

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
    this->element_ = &e;
    finiteElement_ = &(feCache_.get(this->element_->type()));
    using Factory = PQkDGRefinedDGNodeFactory<GV, level, k, void>;
    size_type numberOfSubElements;
    if(e.type().isTriangle()) {
      numberOfSubElements = Factory::numberOfSubTriangles;
    } else if(e.type().isQuadrilateral()) {
      numberOfSubElements = Factory::numberOfSubQuads;
    } else {
      DUNE_THROW(Dune::NotImplemented,
                 "PQkDGRefinedNode::bind() not implemented for element type "
                 << e.type().id());
    }
    this->setSize(numberOfSubElements*finiteElement_->size());
  }

protected:

  FiniteElementCache feCache_;
  const FiniteElement* finiteElement_;
};



template<typename GV, int level, int k, class MI, class TP>
class PQkDGRefinedDGNodeIndexSet
{
  // Cannot be an enum -- otherwise the switch statement below produces compiler warnings
  static const int dim = GV::dimension;

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = PQkDGRefinedDGNodeFactory<GV, level, k, MI>;

  using Node = typename NodeFactory::template Node<TP>;

  PQkDGRefinedDGNodeIndexSet(const NodeFactory& nodeFactory) :
    nodeFactory_(&nodeFactory)
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
#if DUNE_VERSION_NEWER(DUNE_GRID,2,6)
  template<typename It>
  It indices(It it) const
  {
    assert(node_ != nullptr);
    const auto& gridIndexSet = nodeFactory_->gridView().indexSet();
    const auto& element = node_->element();

    for (size_type i = 0, end = this->size(); i < end; ++it, ++i)
    {
      switch (dim)
      {
        case 1:
        {
          *it = {{ nodeFactory_->dofsPerEdge
                   * gridIndexSet.subIndex(element,0,0) + i }};
          continue;
        }
        case 2:
        {
          if (element.type().isTriangle())
          {
            *it = {{ nodeFactory_->dofsPerTriangle
                     * gridIndexSet.subIndex(element,0,0) + i }};
            continue;
          }
          else if (element.type().isQuadrilateral())
          {
            *it = {{ nodeFactory_->quadrilateralOffset_
                     + nodeFactory_->dofsPerQuad
                       * gridIndexSet.subIndex(element,0,0) + i }};
            continue;
          }
          else
            DUNE_THROW(Dune::NotImplemented,
                "2d elements have to be triangles or quadrilaterals");
        }
      }
      DUNE_THROW(Dune::NotImplemented, "No index method for " << dim << "d grids available yet!");
    }
    return it;
  }
#else
  MultiIndex index(size_type i) const
  {
    const auto& gridIndexSet = nodeFactory_->gridView().indexSet();
    const auto& element = node_->element();

    switch (dim)
    {
      case 1:
      {
        return {nodeFactory_->dofsPerEdge*gridIndexSet.subIndex(element,0,0) + i};
      }
      case 2:
      {
        if (element.type().isTriangle())
        {
          return {nodeFactory_->dofsPerTriangle*gridIndexSet.subIndex(element,0,0) + i};
        }
        else if (element.type().isQuadrilateral())
        {
          return { nodeFactory_->quadrilateralOffset_ + nodeFactory_->dofsPerQuad*gridIndexSet.subIndex(element,0,0) + i};
        }
        else
          DUNE_THROW(Dune::NotImplemented, "2d elements have to be triangles or quadrilaterals");
      }
    }
    DUNE_THROW(Dune::NotImplemented, "No index method for " << dim << "d grids available yet!");
  }
#endif

protected:
  const NodeFactory* nodeFactory_;

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
using PQkDGRefinedDGBasis = DefaultGlobalBasis<PQkDGRefinedDGNodeFactory<GV, level, k, FlatMultiIndex<std::size_t> > >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKDGREFINEDDGBASIS_HH
