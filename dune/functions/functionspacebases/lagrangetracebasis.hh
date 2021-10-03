// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGETRACEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGETRACEBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/power.hh>
#include <dune/common/version.hh>

#include <dune/localfunctions/lagrange/pqktracefactory.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>


namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   LagrangeTracePreBasis
//   LagrangeTraceNodeIndexSet
//   LagrangeTraceNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
template<typename GV, int k>
class LagrangeTraceNode;

template<typename GV, int k, class MI>
class LagrangeTraceNodeIndexSet;
#else
template<typename GV, int k, typename TP>
class LagrangeTraceNode;

template<typename GV, int k, class MI, class TP>
class LagrangeTraceNodeIndexSet;
#endif

template<typename GV, int k, class MI>
class LagrangeTracePreBasis;



template<typename GV, int k, class MI>
class LagrangeTracePreBasis
{
  static constexpr int dim = GV::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = std::size_t;


  // Precompute the number of dofs per entity type
  static constexpr int dofsPerEdge     = k-1;
  static constexpr int dofsPerTriangle = (k-1)*(k-2)/2;
  static constexpr int dofsPerQuad     = (k-1)*(k-1);


#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
  using Node = LagrangeTraceNode<GV, k>;

  using IndexSet = LagrangeTraceNodeIndexSet<GV, k, MI>;
#else
  template<class TP>
  using Node = LagrangeTraceNode<GV, k, TP>;

  template<class TP>
  using IndexSet = LagrangeTraceNodeIndexSet<GV, k, MI, TP>;
#endif

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  /** \brief Constructor for a given grid view object */
  LagrangeTracePreBasis(const GridView& gv) :
    gridView_(gv)
  {}


  void initializeIndices()
  {
    vertexOffset_        = 0;
    edgeOffset_          = vertexOffset_ + gridView_.size(dim);
    if constexpr (dim==3)
    {
      triangleOffset_      = edgeOffset_
                           + dofsPerEdge * gridView_.size(dim-1);

      quadrilateralOffset_ = triangleOffset_
                             + dofsPerTriangle
                               * gridView_.size(GeometryTypes::triangle);
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
      return gridView_.size(dim);
    } else if constexpr (dim == 2) {
      return gridView_.size(dim) + dofsPerEdge * gridView_.size(1);
    } else if constexpr (dim == 3) {
      return gridView_.size(dim) + dofsPerEdge * gridView_.size(2)
           + dofsPerTriangle * gridView_.size(GeometryTypes::triangle)
           + dofsPerQuad * gridView_.size(GeometryTypes::quadrilateral);
    } else {
      static_assert(dim >= 1 && dim <= 3,
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
    return StaticPower<(k+1),GV::dimension>::power
           - StaticPower<(k-1),GV::dimension>::power;
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  template<typename It>
  It indices(const Node& node, It it) const
  {
    const auto& gridIndexSet = gridView().indexSet();
    const auto& element = node.element();
    const auto& finiteElement = node.finiteElement();

    for (size_type i = 0, end = finiteElement.size(); i < end; ++it, ++i)
    {
      const Dune::LocalKey localKey
          = finiteElement.localCoefficients().localKey(i);
      // The dimension of the entity that the current dof is related to
      const size_t dofDim = dim - localKey.codim();
      if (dofDim==0) {  // vertex dof
        *it = {{ gridIndexSet.subIndex(element,localKey.subEntity(),dim) }};
        continue;
      }

      if (dofDim==1)
      {  // edge dof
        if constexpr (dim==1)   // element dof -- any local numbering is fine
        {
          DUNE_THROW(Dune::NotImplemented,
              "traces have no elements of codimension 0");
        }
        else
        {
          const auto refElement
              = Dune::referenceElement<double,dim>(element.type());

          // we have to reverse the numbering if the local triangle edge is
          // not aligned with the global edge
          size_t v0 = gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),0,dim),dim);
          size_t v1 = gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),1,dim),dim);
          bool flip = (v0 > v1);
          *it = {{ (flip)
                   ? edgeOffset_
                     + (k-1)*gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim())
                     + (k-2)-localKey.index()
                   : edgeOffset_
                     + (k-1)*gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim())
                     + localKey.index() }};
          continue;
        }
      }

      if (dofDim==2)
      {
        if constexpr (dim==2)   // element dof -- any local numbering is fine
        {
          DUNE_THROW(Dune::NotImplemented,
              "traces have no elements of codimension 0");
        } else
        {
          const auto refElement
              = Dune::referenceElement<double,dim>(element.type());

          if (! refElement.type(localKey.subEntity(), localKey.codim()).isTriangle()
              or k>3)
            DUNE_THROW(Dune::NotImplemented, "LagrangeTraceNodalBasis for 3D grids is only implemented if k<=3 and if the grid is a simplex grid");

          *it = {{ triangleOffset_
                   + gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()) }};
          continue;
        }
      }
      DUNE_THROW(Dune::NotImplemented,
          "Grid contains elements not supported for the LagrangeTraceNodalBasis");
    }
    return it;
  }

//protected:
  GridView gridView_;

  size_type vertexOffset_;
  size_type edgeOffset_;
  size_type triangleOffset_;
  size_type quadrilateralOffset_;

};



#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
template<typename GV, int k>
class LagrangeTraceNode :
  public LeafBasisNode
#else
template<typename GV, int k, typename TP>
class LagrangeTraceNode :
  public LeafBasisNode<std::size_t, TP>
#endif
{
  static constexpr int dim = GV::dimension;

#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,7)
  using Base = LeafBasisNode<std::size_t, TP>;
#endif
  using FiniteElementCache = typename Dune::PQkTraceLocalFiniteElementCache
                                        <typename GV::ctype, double, dim, k>;

public:

  using size_type = std::size_t;
#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,7)
  using TreePath = TP;
#endif
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
  LagrangeTraceNode() :
#else
  LagrangeTraceNode(const TreePath& treePath) :
    Base(treePath),
#endif
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



#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
template<typename GV, int k, class MI>
#else
template<typename GV, int k, class MI, class TP>
#endif
class LagrangeTraceNodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = LagrangeTracePreBasis<GV, k, MI>;

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
  using Node = typename PreBasis::Node;
#else
  using Node = typename PreBasis::template Node<TP>;
#endif

  LagrangeTraceNodeIndexSet(const PreBasis& preBasis) :
    preBasis_(&preBasis),
    node_(nullptr)
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



namespace BasisFactory {

namespace Imp {

template<std::size_t k>
struct LagrangeTracePreBasisFactory
{
  static const std::size_t requiredMultiIndexSize = 1;

  template<class MultiIndex, class GridView>
  auto makePreBasis(const GridView& gridView) const
  {
    return LagrangeTracePreBasis<GridView, k, MultiIndex>(gridView);
  }
};

} // end namespace BasisFactory::Imp

template<std::size_t k>
auto lagrangeTrace()
{
  return Imp::LagrangeTracePreBasisFactory<k>();
}

} // end namespace BasisFactory



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Nodal basis of a scalar k-th-order Lagrangean finite element space
 *         without inner dofs
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int k>
using LagrangeTraceBasis = DefaultGlobalBasis<LagrangeTracePreBasis<GV, k, FlatMultiIndex<std::size_t> > >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGETRACEBASIS_HH
