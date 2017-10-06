// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKDGSUBSAMPLEDDGBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKDGSUBSAMPLEDDGBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/version.hh>

#include <dune/localfunctions/lagrange/pqkdgsubsampledfactory.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>


namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   PQkDGSubsampledDGNodeFactory
//   PQkDGSubsampledDGNodeIndexSet
//   PQkDGSubsampledDGNode
//
// The factory allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, int s, int k, typename TP>
class PQkDGSubsampledDGNode;

template<typename GV, int s, int k, class MI, class TP>
class PQkDGSubsampledDGNodeIndexSet;

template<typename GV, int s, int k, class MI>
class PQkDGSubsampledDGNodeFactory;



template<typename GV, int s, int k, class MI>
class PQkDGSubsampledDGNodeFactory
{
  static const int dim = GV::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = std::size_t;


  // Precompute the number of dofs per entity type
  const int dofsPerEdge        = s*(k+1);
  const int dofsPerTriangle    = s*s*((k+1)*(k+2))/2;
  const int dofsPerQuad        = s*s*(k+1)*(k+1);


  template<class TP>
  using Node = PQkDGSubsampledDGNode<GV, s, k, TP>;

  template<class TP>
  using IndexSet = PQkDGSubsampledDGNodeIndexSet<GV, s, k, MI, TP>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  /** \brief Constructor for a given grid view object */
  PQkDGSubsampledDGNodeFactory(const GridView& gv) :
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
        return dofsPerEdge*gridView_.size(0);
      case 2:
      {
        return dofsPerTriangle * gridView_.size(GeometryTypes::triangle)
                + dofsPerQuad * gridView_.size(GeometryTypes::quadrilateral);
      }
    }

    DUNE_THROW(Dune::NotImplemented, "No size method for " << dim << "d grids available yet!");
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
    return StaticPower<s,GV::dimension>::power
           * StaticPower<k+1,GV::dimension>::power;
  }

//protected:
  GridView gridView_;

  size_type quadrilateralOffset_;

};



template<typename GV, int s, int k, typename TP>
class PQkDGSubsampledDGNode :
  public LeafBasisNode<std::size_t, TP>
{
  static const int dim = GV::dimension;
  static const int maxSize = StaticPower<s,dim>::power
                             * StaticPower<k+1,dim>::power;

  using Base = LeafBasisNode<std::size_t, TP>;
  using FiniteElementCache
      = typename Dune::PQkDGSubsampledLocalFiniteElementCache
                        <typename GV::ctype, double, dim, s, k>;

public:

  using size_type = std::size_t;
  using TreePath = TP;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  PQkDGSubsampledDGNode(const TreePath& treePath) :
    Base(treePath),
    finiteElement_(nullptr),
    element_(nullptr)
  {}

  //! Return current element, throw if unbound
  const Element& element() const
  {
    return *element_;
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
    element_ = &e;
    finiteElement_ = &(cache_.get(element_->type()));
    this->setSize(finiteElement_->size());
  }

protected:

  FiniteElementCache cache_;
  const FiniteElement* finiteElement_;
  const Element* element_;
};



template<typename GV, int s, int k, class MI, class TP>
class PQkDGSubsampledDGNodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = PQkDGSubsampledDGNodeFactory<GV, s, k, MI>;

  using Node = typename NodeFactory::template Node<TP>;

  PQkDGSubsampledDGNodeIndexSet(const NodeFactory& nodeFactory) :
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
    return node_->finiteElement().size();
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
      DUNE_THROW(Dune::NotImplemented,
                 "No index method for " << dim << "d grids available yet!");
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
        return {nodeFactory_->dofsPerEdge
                * gridIndexSet.subIndex(element,0,0) + i};
      }
      case 2:
      {
        if (element.type().isTriangle())
        {
          return {nodeFactory_->dofsPerTriangle
                  * gridIndexSet.subIndex(element,0,0) + i};
        }
        else if (element.type().isQuadrilateral())
        {
          return {nodeFactory_->quadrilateralOffset_
                  + nodeFactory_->dofsPerQuad
                    * gridIndexSet.subIndex(element,0,0) + i};
        }
        else
          DUNE_THROW(Dune::NotImplemented,
                     "2d elements have to be triangles or quadrilaterals");
      }
    }
    DUNE_THROW(Dune::NotImplemented,
               "No index method for " << dim << "d grids available yet!");
  }
#endif

protected:
  const NodeFactory* nodeFactory_;

  const Node* node_;
};



namespace BasisBuilder {

namespace Imp {

template<std::size_t s, std::size_t k>
struct PQkDGSubsampledDGNodeFactoryBuilder
{
  static const std::size_t requiredMultiIndexSize=1;

  template<class MultiIndex, class GridView>
  auto build(const GridView& gridView)
    -> PQkDGSubsampledDGNodeFactory<GridView, s, k, MultiIndex>
  {
    return {gridView};
  }
};

} // end namespace BasisBuilder::Imp

template<std::size_t s, std::size_t k>
Imp::PQkDGSubsampledDGNodeFactoryBuilder<s, k>
pqDGSubsampledDG()
{
  return{};
}

} // end namespace BasisBuilder



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Nodal basis of a scalar k-th-order discontinuous Lagrangean finite element space
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int s, int k>
using PQkDGSubsampledDGNodalBasis = DefaultGlobalBasis<PQkDGSubsampledDGNodeFactory<GV, s, k, FlatMultiIndex<std::size_t> > >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKDGSUBSAMPLEDDGBASIS_HH
