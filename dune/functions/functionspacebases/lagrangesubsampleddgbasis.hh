// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGESUBSAMPLEDDGBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGESUBSAMPLEDDGBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/power.hh>
#include <dune/common/version.hh>

#include <dune/localfunctions/lagrange/pqksubsampledfactory.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>


namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   LagrangeSubsampledDGPreBasis
//   LagrangeSubsampledDGNodeIndexSet
//   LagrangeSubsampledDGNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
template<typename GV, int s, int k>
class LagrangeSubsampledDGNode;

template<typename GV, int s, int k, class MI>
class LagrangeSubsampledDGNodeIndexSet;
#else
template<typename GV, int s, int k, typename TP>
class LagrangeSubsampledDGNode;

template<typename GV, int s, int k, class MI, class TP>
class LagrangeSubsampledDGNodeIndexSet;
#endif

template<typename GV, int s, int k, class MI>
class LagrangeSubsampledDGPreBasis;



template<typename GV, int s, int k, class MI>
class LagrangeSubsampledDGPreBasis
{
  static constexpr int dim = GV::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = std::size_t;


  // Precompute the number of dofs per entity type
  static constexpr int dofsPerEdge        = s*k+1;
  static constexpr int dofsPerTriangle    = ((s*k+1)*(s*k+2))/2;
  static constexpr int dofsPerQuad        = (s*k+1)*(s*k+1);
  static constexpr int dofsPerTetrahedron = (s*k+1)*(s*k+2)*(s*k+3)/6;
  static constexpr int dofsPerPrism       = ((s*k+1)*(s*k+1)*(s*k+2))/2;
  static constexpr int dofsPerHexahedron  = (s*k+1)*(s*k+1)*(s*k+1);
  static constexpr int dofsPerPyramid     = ((s*k+1)*(s*k+2)*(2*s*k+3))/6;


#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
  using Node = LagrangeSubsampledDGNode<GV, s, k>;

  using IndexSet = LagrangeSubsampledDGNodeIndexSet<GV, s, k, MI>;
#else
  template<class TP>
  using Node = LagrangeSubsampledDGNode<GV, s, k, TP>;

  template<class TP>
  using IndexSet = LagrangeSubsampledDGNodeIndexSet<GV, s, k, MI, TP>;
#endif

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  /** \brief Constructor for a given grid view object */
  LagrangeSubsampledDGPreBasis(const GridView& gv) :
    gridView_(gv)
  {}


  void initializeIndices()
  {
    if constexpr (dim == 1) {
      return;
    } else if constexpr (dim == 2) {
      quadrilateralOffset_ = dofsPerTriangle
                             * gridView_.size(GeometryTypes::triangle);
    } else if constexpr (dim == 3) {
      prismOffset_         = dofsPerTetrahedron
                             * gridView_.size(GeometryTypes::tetrahedron);

      hexahedronOffset_    = prismOffset_
                             + dofsPerPrism
                               * gridView_.size(GeometryTypes::prism);

      pyramidOffset_       = hexahedronOffset_
                             + dofsPerHexahedron
                               * gridView_.size(GeometryTypes::hexahedron);
    } else {
      static_assert(dim >= 1 && dim <= 3,
          "LagrangeSubsampledDGPreBasis not implemented for Grids of this dimension.");
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
      return dofsPerEdge*gridView_.size(0);
    } else if constexpr (dim == 2) {
      return dofsPerTriangle * gridView_.size(GeometryTypes::triangle)
           + dofsPerQuad * gridView_.size(GeometryTypes::quadrilateral);
    } else if constexpr (dim == 3) {
      return dofsPerTetrahedron * gridView_.size(GeometryTypes::tetrahedron)
           + dofsPerPyramid * gridView_.size(GeometryTypes::pyramid)
           + dofsPerPrism * gridView_.size(GeometryTypes::prism)
           + dofsPerHexahedron * gridView_.size(GeometryTypes::hexahedron);
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
    return StaticPower<(s*k+1),GV::dimension>::power;
  }

//protected:
  GridView gridView_;

  size_type quadrilateralOffset_;
  size_type pyramidOffset_;
  size_type prismOffset_;
  size_type hexahedronOffset_;

};



#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
template<typename GV, int s, int k>
class LagrangeSubsampledDGNode :
  public LeafBasisNode
#else
template<typename GV, int s, int k, typename TP>
class LagrangeSubsampledDGNode :
  public LeafBasisNode<std::size_t, TP>
#endif
{
  static constexpr int dim = GV::dimension;

#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,7)
  using Base = LeafBasisNode<std::size_t, TP>;
#endif
  using FiniteElementCache
      = typename Dune::PQkSubsampledLocalFiniteElementCache
                        <typename GV::ctype, double, dim, s, k>;

public:

  using size_type = std::size_t;
#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,7)
  using TreePath = TP;
#endif
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
  LagrangeSubsampledDGNode() :
#else
  LagrangeSubsampledDGNode(const TreePath& treePath) :
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
template<typename GV, int s, int k, class MI>
#else
template<typename GV, int s, int k, class MI, class TP>
#endif
class LagrangeSubsampledDGNodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = LagrangeSubsampledDGPreBasis<GV, s, k, MI>;

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
  using Node = typename PreBasis::Node;
#else
  using Node = typename PreBasis::template Node<TP>;
#endif

  LagrangeSubsampledDGNodeIndexSet(const PreBasis& preBasis) :
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
    return node_->finiteElement().size();
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
          DUNE_THROW(Dune::NotImplemented, "2d elements have to be triangles or quadrilaterals");
      } else if constexpr (dim == 3) {
        if (element.type().isTetrahedron())
        {
          *it = {{ preBasis_->dofsPerTetrahedron
                   * gridIndexSet.subIndex(element,0,0) + i }};
          continue;
        }
        else if (element.type().isPrism())
        {
          *it = {{ preBasis_->prismOffset_
                   + preBasis_->dofsPerPrism
                     * gridIndexSet.subIndex(element,0,0) + i }};
          continue;
        }
        else if (element.type().isHexahedron())
        {
          *it = {{ preBasis_->hexahedronOffset_
                   + preBasis_->dofsPerHexahedron
                     * gridIndexSet.subIndex(element,0,0) + i }};
          continue;
        }
        else if (element.type().isPyramid())
        {
          *it = {{ preBasis_->pyramidOffset_
                   + preBasis_->dofsPerPyramid
                     * gridIndexSet.subIndex(element,0,0) + i }};
          continue;
        }
        else
          DUNE_THROW(Dune::NotImplemented, "3d elements have to be tetrahedrons, prisms, hexahedrons or pyramids");
      } else {
        static_assert(dim >= 1 && dim <= 3,
            "The index method is not yet implemented for grids of this dimension!");
      }
    }
    return it;
  }

protected:
  const PreBasis* preBasis_;

  const Node* node_;
};



namespace BasisFactory {

namespace Imp {

template<std::size_t s, std::size_t k>
struct LagrangeSubsampledDGPreBasisFactory
{
  static const std::size_t requiredMultiIndexSize = 1;

  template<class MultiIndex, class GridView>
  auto makePreBasis(const GridView& gridView) const
  {
    return LagrangeSubsampledDGPreBasis<GridView, s, k, MultiIndex>(gridView);
  }
};

} // end namespace BasisFactory::Imp

template<std::size_t s, std::size_t k>
auto pqSubsampledDG()
{
  return Imp::LagrangeSubsampledDGPreBasisFactory<s, k>();
}

} // end namespace BasisFactory



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
using LagrangeSubsampledDGBasis = DefaultGlobalBasis<LagrangeSubsampledDGPreBasis<GV, s, k, FlatMultiIndex<std::size_t> > >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGESUBSAMPLEDDGBASIS_HH
