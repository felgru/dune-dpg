// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGESUBSAMPLEDDGBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGESUBSAMPLEDDGBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/math.hh>

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

template<typename GV, int s, int k, typename R=double>
class LagrangeSubsampledDGNode;

template<typename GV, int s, int k, class MI, typename R=double>
class LagrangeSubsampledDGNodeIndexSet;

template<typename GV, int s, int k, class MI, typename R=double>
class LagrangeSubsampledDGPreBasis;



template<typename GV, int s, int k, class MI, typename R>
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


  using Node = LagrangeSubsampledDGNode<GV, s, k, R>;

  using IndexSet = LagrangeSubsampledDGNodeIndexSet<GV, s, k, MI, R>;

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

  Node makeNode() const
  {
    return Node{};
  }

  IndexSet makeIndexSet() const
  {
    return IndexSet{*this};
  }

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
    return power(s*k+1, dim);
  }

//protected:
  GridView gridView_;

  size_type quadrilateralOffset_;
  size_type pyramidOffset_;
  size_type prismOffset_;
  size_type hexahedronOffset_;

};



template<typename GV, int s, int k, typename R>
class LagrangeSubsampledDGNode :
  public LeafBasisNode
{
  static constexpr int dim = GV::dimension;

  using FiniteElementCache
      = typename Dune::PQkSubsampledLocalFiniteElementCache
                        <typename GV::ctype, R, dim, s, k>;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  LagrangeSubsampledDGNode() :
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



template<typename GV, int s, int k, class MI, typename R>
class LagrangeSubsampledDGNodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = LagrangeSubsampledDGPreBasis<GV, s, k, MI, R>;

  using Node = typename PreBasis::Node;

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

template<std::size_t s, std::size_t k, typename R=double>
struct LagrangeSubsampledDGPreBasisFactory
{
  static const std::size_t requiredMultiIndexSize = 1;

  template<class MultiIndex, class GridView>
  auto makePreBasis(const GridView& gridView) const
  {
    return LagrangeSubsampledDGPreBasis<GridView, s, k, MultiIndex, R>(gridView);
  }
};

} // end namespace BasisFactory::Imp

template<std::size_t s, std::size_t k, typename R=double>
auto pqSubsampledDG()
{
  return Imp::LagrangeSubsampledDGPreBasisFactory<s, k, R>();
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
 * \tparam R The range type of the local basis
 */
template<typename GV, int s, int k>
using LagrangeSubsampledDGBasis = DefaultGlobalBasis<LagrangeSubsampledDGPreBasis<GV, s, k, FlatMultiIndex<std::size_t>, R> >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGESUBSAMPLEDDGBASIS_HH
