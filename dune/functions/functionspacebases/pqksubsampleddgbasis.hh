// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKSUBSAMPLEDDGBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKSUBSAMPLEDDGBASIS_HH

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
//   PQkSubsampledDGPreBasis
//   PQkSubsampledDGNodeIndexSet
//   PQkSubsampledDGNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, int s, int k, typename TP>
class PQkSubsampledDGNode;

template<typename GV, int s, int k, class MI, class TP>
class PQkSubsampledDGNodeIndexSet;

template<typename GV, int s, int k, class MI>
class PQkSubsampledDGPreBasis;



template<typename GV, int s, int k, class MI>
class PQkSubsampledDGPreBasis
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


  template<class TP>
  using Node = PQkSubsampledDGNode<GV, s, k, TP>;

  template<class TP>
  using IndexSet = PQkSubsampledDGNodeIndexSet<GV, s, k, MI, TP>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  /** \brief Constructor for a given grid view object */
  PQkSubsampledDGPreBasis(const GridView& gv) :
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
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
        quadrilateralOffset_ = dofsPerTriangle
                               * gridView_.size(GeometryTypes::triangle);
#else
        GeometryType triangle;
        triangle.makeTriangle();
        quadrilateralOffset_ = dofsPerTriangle * gridView_.size(triangle);
#endif
        break;
      }
      case 3:
      {
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
        prismOffset_         = dofsPerTetrahedron
                               * gridView_.size(GeometryTypes::tetrahedron);

        hexahedronOffset_    = prismOffset_
                               + dofsPerPrism
                                 * gridView_.size(GeometryTypes::prism);

        pyramidOffset_       = hexahedronOffset_
                               + dofsPerHexahedron
                                 * gridView_.size(GeometryTypes::hexahedron);
#else
        GeometryType tetrahedron;
        tetrahedron.makeSimplex(3);
        prismOffset_         = dofsPerTetrahedron * gridView_.size(tetrahedron);

        GeometryType prism;
        prism.makePrism();
        hexahedronOffset_    = prismOffset_
                             + dofsPerPrism * gridView_.size(prism);

        GeometryType hexahedron;
        hexahedron.makeCube(3);
        pyramidOffset_       = hexahedronOffset_
                             + dofsPerHexahedron * gridView_.size(hexahedron);
#endif
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
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
        return dofsPerTriangle * gridView_.size(GeometryTypes::triangle)
             + dofsPerQuad * gridView_.size(GeometryTypes::quadrilateral);
#else
        GeometryType triangle, quad;
        triangle.makeTriangle();
        quad.makeQuadrilateral();
        return dofsPerTriangle * gridView_.size(triangle)
             + dofsPerQuad * gridView_.size(quad);
#endif
      }
      case 3:
      {
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
        return dofsPerTetrahedron * gridView_.size(GeometryTypes::tetrahedron)
             + dofsPerPyramid * gridView_.size(GeometryTypes::pyramid)
             + dofsPerPrism * gridView_.size(GeometryTypes::prism)
             + dofsPerHexahedron * gridView_.size(GeometryTypes::hexahedron);
#else
        GeometryType tetrahedron, pyramid, prism, hexahedron;
        tetrahedron.makeTetrahedron();
        pyramid.makePyramid();
        prism.makePrism();
        hexahedron.makeCube(3);
        return dofsPerTetrahedron * gridView_.size(tetrahedron)
             + dofsPerPyramid * gridView_.size(pyramid)
             + dofsPerPrism * gridView_.size(prism)
             + dofsPerHexahedron * gridView_.size(hexahedron);
#endif
      }

    }

    DUNE_THROW(Dune::NotImplemented,
               "No size method for " << dim << "d grids available yet!");
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



template<typename GV, int s, int k, typename TP>
class PQkSubsampledDGNode :
  public LeafBasisNode<std::size_t, TP>
{
  static constexpr int dim = GV::dimension;

  using Base = LeafBasisNode<std::size_t, TP>;
  using FiniteElementCache
      = typename Dune::PQkSubsampledLocalFiniteElementCache
                        <typename GV::ctype, double, dim, s, k>;

public:

  using size_type = std::size_t;
  using TreePath = TP;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  PQkSubsampledDGNode(const TreePath& treePath) :
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
class PQkSubsampledDGNodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = PQkSubsampledDGPreBasis<GV, s, k, MI>;

  using Node = typename PreBasis::template Node<TP>;

  PQkSubsampledDGNodeIndexSet(const PreBasis& preBasis) :
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
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,6)
  template<typename It>
  It indices(It it) const
  {
    assert(node_ != nullptr);
    const auto& gridIndexSet = preBasis_->gridView().indexSet();
    const auto& element = node_->element();

    for (size_type i = 0, end = this->size(); i < end; ++it, ++i)
    {
      switch (dim)
      {
        case 1:
        {
          *it = {{ preBasis_->dofsPerEdge
                   * gridIndexSet.subIndex(element,0,0) + i }};
          continue;
        }
        case 2:
        {
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
        }
        case 3:
        {
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
        }
      }
      DUNE_THROW(Dune::NotImplemented, "No index method for " << dim << "d grids available yet!");
    }
    return it;
  }
#else
  MultiIndex index(size_type i) const
  {
    const auto& gridIndexSet = preBasis_->gridView().indexSet();
    const auto& element = node_->element();

    switch (dim)
    {
      case 1:
      {
        return {preBasis_->dofsPerEdge*gridIndexSet.subIndex(element,0,0) + i};
      }
      case 2:
      {
        if (element.type().isTriangle())
        {
          return {preBasis_->dofsPerTriangle*gridIndexSet.subIndex(element,0,0) + i};
        }
        else if (element.type().isQuadrilateral())
        {
          return { preBasis_->quadrilateralOffset_ + preBasis_->dofsPerQuad*gridIndexSet.subIndex(element,0,0) + i};
        }
        else
          DUNE_THROW(Dune::NotImplemented, "2d elements have to be triangles or quadrilaterals");
      }
      case 3:
      {
        if (element.type().isTetrahedron())
        {
          return {preBasis_->dofsPerTetrahedron*gridIndexSet.subIndex(element,0,0) + i};
        }
        else if (element.type().isPrism())
        {
          return { preBasis_->prismOffset_ + preBasis_->dofsPerPrism*gridIndexSet.subIndex(element,0,0) + i};
        }
        else if (element.type().isHexahedron())
        {
          return { preBasis_->hexahedronOffset_ + preBasis_->dofsPerHexahedron*gridIndexSet.subIndex(element,0,0) + i};
        }
        else if (element.type().isPyramid())
        {
          return { preBasis_->pyramidOffset_ + preBasis_->dofsPerPyramid*gridIndexSet.subIndex(element,0,0) + i};
        }
        else
          DUNE_THROW(Dune::NotImplemented, "3d elements have to be tetrahedrons, prisms, hexahedrons or pyramids");
      }
    }
    DUNE_THROW(Dune::NotImplemented, "No index method for " << dim << "d grids available yet!");
  }
#endif

protected:
  const PreBasis* preBasis_;

  const Node* node_;
};



namespace BasisBuilder {

namespace Imp {

template<std::size_t s, std::size_t k>
struct PQkSubsampledDGPreBasisFactory
{
  static const std::size_t requiredMultiIndexSize = 1;

  template<class MultiIndex, class GridView>
  auto makePreBasis(const GridView& gridView) const
  {
    return PQkSubsampledDGPreBasis<GridView, s, k, MultiIndex>(gridView);
  }
};

} // end namespace BasisBuilder::Imp

template<std::size_t s, std::size_t k>
auto pqSubsampledDG()
{
  return Imp::PQkSubsampledDGPreBasisFactory<s, k>();
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
using PQkSubsampledDGNodalBasis = DefaultGlobalBasis<PQkSubsampledDGPreBasis<GV, s, k, FlatMultiIndex<std::size_t> > >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKSUBSAMPLEDDGBASIS_HH
