// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BERNSTEINDGBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BERNSTEINDGBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/power.hh>
#include <dune/common/version.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/bernsteinbasis.hh>




namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   BernsteinDGPreBasis
//   BernsteinDGNodeIndexSet
//   BernsteinDGNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
template<typename GV, int k>
using BernsteinDGNode = BernsteinNode<GV, k>;

template<typename GV, int k, class MI>
class BernsteinDGNodeIndexSet;
#else
template<typename GV, int k, typename TP>
using BernsteinDGNode = BernsteinNode<GV, k, TP>;

template<typename GV, int k, class MI, class TP>
class BernsteinDGNodeIndexSet;
#endif


template<typename GV, int k, class MI>
class BernsteinDGPreBasis
{
  static constexpr int dim = GV::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = std::size_t;


  // Precompute the number of dofs per entity type
  constexpr static int dofsPerEdge        = k+1;
  constexpr static int dofsPerTriangle    = (k+1)*(k+2)/2;
  constexpr static int dofsPerQuad        = (k+1)*(k+1);
  constexpr static int dofsPerTetrahedron = (k+1)*(k+2)*(k+3)/6;
  constexpr static int dofsPerPrism       = (k+1)*(k+1)*(k+2)/2;
  constexpr static int dofsPerHexahedron  = (k+1)*(k+1)*(k+1);
  constexpr static int dofsPerPyramid     = (k+1)*(k+2)*(2*k+3)/6;


#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
  using Node = BernsteinDGNode<GV, k>;

  using IndexSet = BernsteinDGNodeIndexSet<GV, k, MI>;
#else
  template<class TP>
  using Node = BernsteinDGNode<GV, k, TP>;

  template<class TP>
  using IndexSet = BernsteinDGNodeIndexSet<GV, k, MI, TP>;
#endif

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  /** \brief Constructor for a given grid view object */
  BernsteinDGPreBasis(const GridView& gv) :
    gridView_(gv)
  {}


  void initializeIndices()
  {
    if constexpr (dim == 1) {
      return;
    } else if constexpr (dim == 2) {
      quadrilateralOffset_ = dofsPerTriangle * gridView_.size(Dune::GeometryTypes::triangle);
    } else if constexpr (dim == 3) {
      prismOffset_         = dofsPerTetrahedron * gridView_.size(Dune::GeometryTypes::tetrahedron);
      hexahedronOffset_    = prismOffset_         +   dofsPerPrism * gridView_.size(Dune::GeometryTypes::prism);
      pyramidOffset_       = hexahedronOffset_    +   dofsPerHexahedron * gridView_.size(Dune::GeometryTypes::hexahedron);
    }
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return gridView_;
  }

  void update(const GridView& gv)
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
      return dofsPerTriangle * gridView_.size(Dune::GeometryTypes::triangle)
           + dofsPerQuad * gridView_.size(Dune::GeometryTypes::quadrilateral);
    } else if constexpr (dim == 3) {
      return dofsPerTetrahedron*gridView_.size(Dune::GeometryTypes::tetrahedron)
           + dofsPerPyramid*gridView_.size(Dune::GeometryTypes::pyramid)
           + dofsPerPrism*gridView_.size(Dune::GeometryTypes::prism)
           + dofsPerHexahedron*gridView_.size(Dune::GeometryTypes::hexahedron);
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
    return StaticPower<(k+1),GV::dimension>::power;
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  template<typename It>
  It indices(const Node& node, It it) const
  {
    const auto& gridIndexSet = gridView().indexSet();
    const auto& element = node.element();

    for (size_type i = 0, end = size(); i < end; ++i, ++it) {
      if constexpr (dim == 1) {
        *it = {dofsPerEdge*gridIndexSet.subIndex(element,0,0) + i};
        continue;
      } else if constexpr (dim == 2) {
        if (element.type().isTriangle()) {
          *it = {dofsPerTriangle*gridIndexSet.subIndex(element,0,0) + i};
          continue;
        } else if (element.type().isQuadrilateral()) {
          *it = {quadrilateralOffset_
                 + dofsPerQuad*gridIndexSet.subIndex(element,0,0) + i};
          continue;
        } else {
          DUNE_THROW(Dune::NotImplemented, "2d elements have to be triangles or quadrilaterals");
        }
      } else if constexpr (dim == 3) {
        if (element.type().isTetrahedron()) {
          *it = {dofsPerTetrahedron*gridIndexSet.subIndex(element,0,0) + i};
          continue;
        } else if (element.type().isPrism()) {
          *it = {prismOffset_
                 + dofsPerPrism*gridIndexSet.subIndex(element,0,0) + i};
          continue;
        } else if (element.type().isHexahedron()) {
          *it = {hexahedronOffset_
                 + dofsPerHexahedron*gridIndexSet.subIndex(element,0,0) + i};
          continue;
        } else if (element.type().isPyramid()) {
          *it = {pyramidOffset_
                 + dofsPerPyramid*gridIndexSet.subIndex(element,0,0) + i};
          continue;
        } else {
          DUNE_THROW(Dune::NotImplemented, "3d elements have to be tetrahedrons, prisms, hexahedrons or pyramids");
        }
      } else {
        static_assert(dim >= 1 && dim <= 3,
                      "The index method is not yet implemented for grids of this dimension!");
      }
    }
    return it;
  }

//protected:
  GridView gridView_;

  size_t quadrilateralOffset_;
  size_t pyramidOffset_;
  size_t prismOffset_;
  size_t hexahedronOffset_;
};



#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
template<typename GV, int k, class MI>
#else
template<typename GV, int k, class MI, class TP>
#endif
class BernsteinDGNodeIndexSet
{
  static constexpr int dim = GV::dimension;

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = BernsteinDGPreBasis<GV, k, MI>;

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
  using Node = typename PreBasis::Node;
#else
  using Node = typename PreBasis::template Node<TP>;
#endif

  BernsteinDGNodeIndexSet(const PreBasis& preBasis) :
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



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Basis of a scalar k-th-order Bernstein-DG finite element space
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int k>
using BernsteinDGBasis = DefaultGlobalBasis<BernsteinDGPreBasis<GV, k, FlatMultiIndex<std::size_t>> >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BERNSTEINDGBASIS_HH
