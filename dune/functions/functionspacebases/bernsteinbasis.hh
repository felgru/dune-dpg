// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BERNSTEINBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BERNSTEINBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/power.hh>
#include <dune/common/version.hh>

#include <dune/localfunctions/bernstein/pqkfactory.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>


namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   BernsteinPreBasis
//   BernsteinNodeIndexSet
//   BernsteinNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k, typename R=double>
class BernsteinNode;

#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,8)
template<typename GV, int k, class MI, typename R=double>
class BernsteinNodeIndexSet;
#endif

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
template<typename GV, int k, typename R=double>
#else
template<typename GV, int k, class MI, typename R=double>
#endif
class BernsteinPreBasis;



/**
 * \brief A pre-basis for PQ Bernstein bases with given order
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV  The grid view that the FE basis is defined on
 * \tparam k   The polynomial order of ansatz functions
 * \tparam MI  Type to be used for multi-indices
 * \tparam R   Range type used for shape function values
 *
 * \note This only works for 2d simplex grids.
 */
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
template<typename GV, int k, typename R>
#else
template<typename GV, int k, class MI, typename R>
#endif
class BernsteinPreBasis
{
  static constexpr int dim = GV::dimension;

public:

  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

private:

#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,8)
  template<typename, int, class, typename>
  friend class BernsteinNodeIndexSet;
#endif

  // Precompute the number of dofs per entity type
  constexpr static size_type dofsPerVertex =
      k == 0 ? (dim == 0 ? 1 : 0) : 1;
  constexpr static size_type dofsPerEdge =
      k == 0 ? (dim == 1 ? 1 : 0) : k-1;
  constexpr static size_type dofsPerTriangle =
      k == 0 ? (dim == 2 ? 1 : 0) : (k-1)*(k-2)/2;
  constexpr static size_type dofsPerQuad =
      k == 0 ? (dim == 2 ? 1 : 0) : (k-1)*(k-1);
  constexpr static size_type dofsPerTetrahedron =
      k == 0 ? (dim == 3 ? 1 : 0) : (k-3)*(k-2)*(k-1)/6;
  constexpr static size_type dofsPerPrism =
      k == 0 ? (dim == 3 ? 1 : 0) : (k-1)*(k-1)*(k-2)/2;
  constexpr static size_type dofsPerHexahedron =
      k == 0 ? (dim == 3 ? 1 : 0) : (k-1)*(k-1)*(k-1);
  constexpr static size_type dofsPerPyramid =
      k == 0 ? (dim == 3 ? 1 : 0) : (k-2)*(k-1)*(2*k-3)/6;

public:

  using Node = BernsteinNode<GV, k, R>;

#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,8)
  using IndexSet = BernsteinNodeIndexSet<GV, k, MI, R>;
#endif

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
  static constexpr size_type maxMultiIndexSize = 1;
  static constexpr size_type minMultiIndexSize = 1;
  static constexpr size_type multiIndexBufferSize = 1;
#else
  //! Type used for global numbering of the basis vectors
  using MultiIndex = MI;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = Dune::ReservedVector<size_type, 1>;
#endif

  //! Constructor for a given grid view object
  BernsteinPreBasis(const GridView& gv) :
    gridView_(gv)
  {}

  //! Initialize the global indices
  void initializeIndices()
  {
    vertexOffset_        = 0;
    edgeOffset_          = vertexOffset_          + dofsPerVertex * static_cast<size_type>(gridView_.size(dim));
    triangleOffset_      = edgeOffset_            + dofsPerEdge * static_cast<size_type>(gridView_.size(dim-1));

    quadrilateralOffset_ = triangleOffset_        + dofsPerTriangle * static_cast<size_type>(gridView_.size(Dune::GeometryTypes::triangle));

    if constexpr (dim==3) {
      tetrahedronOffset_   = quadrilateralOffset_ + dofsPerQuad * static_cast<size_type>(gridView_.size(Dune::GeometryTypes::quadrilateral));

      prismOffset_         = tetrahedronOffset_   +   dofsPerTetrahedron * static_cast<size_type>(gridView_.size(Dune::GeometryTypes::tetrahedron));

      hexahedronOffset_    = prismOffset_         +   dofsPerPrism * static_cast<size_type>(gridView_.size(Dune::GeometryTypes::prism));

      pyramidOffset_       = hexahedronOffset_    +   dofsPerHexahedron * static_cast<size_type>(gridView_.size(Dune::GeometryTypes::hexahedron));
    }
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView() const
  {
    return gridView_;
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update (const GridView& gv)
  {
    gridView_ = gv;
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    return Node{};
  }

  /**
   * \brief Create tree node index set
   *
   * Create an index set suitable for the tree node obtained
   * by makeNode().
   */
#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,8)
  IndexSet makeIndexSet() const
  {
    return IndexSet{*this};
  }
#endif

  //! Same as size(prefix) with empty prefix
  size_type size() const
  {
    if constexpr (dim == 1) {
      return dofsPerVertex * static_cast<size_type>(gridView_.size(dim))
        + dofsPerEdge * static_cast<size_type>(gridView_.size(dim-1));
    } else if constexpr (dim == 2) {
      return dofsPerVertex * static_cast<size_type>(gridView_.size(dim))
        + dofsPerEdge * static_cast<size_type>(gridView_.size(dim-1))
        + dofsPerTriangle * static_cast<size_type>(gridView_.size(Dune::GeometryTypes::triangle))
        + dofsPerQuad * static_cast<size_type>(gridView_.size(Dune::GeometryTypes::quadrilateral));
    } else if constexpr (dim == 3) {
      return dofsPerVertex * static_cast<size_type>(gridView_.size(dim))
        + dofsPerEdge * static_cast<size_type>(gridView_.size(dim-1))
        + dofsPerTriangle * static_cast<size_type>(gridView_.size(Dune::GeometryTypes::triangle))
        + dofsPerQuad * static_cast<size_type>(gridView_.size(Dune::GeometryTypes::quadrilateral))
        + dofsPerTetrahedron * static_cast<size_type>(gridView_.size(Dune::GeometryTypes::tetrahedron))
        + dofsPerPyramid * static_cast<size_type>(gridView_.size(Dune::GeometryTypes::pyramid))
        + dofsPerPrism * static_cast<size_type>(gridView_.size(Dune::GeometryTypes::prism))
        + dofsPerHexahedron * static_cast<size_type>(gridView_.size(Dune::GeometryTypes::hexahedron));
    } else {
      static_assert(dim >= 1 && dim <= 3,
                    "No size method implemented for this dimension!");
    }
  }

  //! Return number of possible values for next position in multi index
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

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    return size();
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize() const
  {
    return StaticPower<(k+1),GV::dimension>::power;
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  template<typename It>
  It indices(const Node& node, It it) const
  {
    for (size_type i = 0, end = node.finiteElement().size() ; i < end ; ++it, ++i)
      {
        Dune::LocalKey localKey = node.finiteElement().localCoefficients().localKey(i);
        const auto& gridIndexSet = gridView().indexSet();
        const auto& element = node.element();

        // The dimension of the entity that the current dof is related to
        auto dofDim = dim - localKey.codim();

        // The test for k==1 is redundant, but having it here allows the
        // compiler to conclude at compile-time that the dofDim==0 case
        // is the only one that will ever happen.
        // See https://gitlab.dune-project.org/staging/dune-functions/issues/30
        if (k==1 || dofDim==0) {  // vertex dof
          *it = {{ static_cast<size_type>(gridIndexSet.subIndex(element,localKey.subEntity(),dim)) }};
          continue;
        }

        if (dofDim==1)
          {  // edge dof
            if (dim==1)  // element dof -- any local numbering is fine
              {
                *it = {{ edgeOffset_
                         + dofsPerEdge * static_cast<size_type>(gridIndexSet.subIndex(element,0,0))
                         + localKey.index() }};
                continue;
              }
            else
              {
                const auto refElement
                  = Dune::referenceElement<double,dim>(element.type());

                // we have to reverse the numbering if the local triangle edge is
                // not aligned with the global edge
                const auto v0 = static_cast<size_type>(gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),0,dim),dim));
                const auto v1 = static_cast<size_type>(gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),1,dim),dim));
                const bool flip = (v0 > v1);
                *it = {{ (flip)
                         ? edgeOffset_
                         + dofsPerEdge*static_cast<size_type>(gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()))
                         + (dofsPerEdge-1)-localKey.index()
                         : edgeOffset_
                         + dofsPerEdge*static_cast<size_type>(gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()))
                         + localKey.index() }};
                continue;
              }
          }

        if (dofDim==2)
          {
            if (dim==2)   // element dof -- any local numbering is fine
              {
                if (element.type().isTriangle())
                  {
                    const int interiorLagrangeNodesPerTriangle = (k-1)*(k-2)/2;
                    *it = {{ triangleOffset_ + interiorLagrangeNodesPerTriangle*static_cast<size_type>(gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                    continue;
                  }
                else if (element.type().isQuadrilateral())
                  {
                    const int interiorLagrangeNodesPerQuadrilateral = (k-1)*(k-1);
                    *it = {{ quadrilateralOffset_ + interiorLagrangeNodesPerQuadrilateral*static_cast<size_type>(gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                    continue;
                  }
                else
                  DUNE_THROW(Dune::NotImplemented, "2d elements have to be triangles or quadrilaterals");
              } else
              {
                const auto refElement
                  = Dune::referenceElement<double,dim>(element.type());

                if (k>3)
                  DUNE_THROW(Dune::NotImplemented, "BernsteinBasis for 3D grids is only implemented if k<=3");

                if (k==3 and !refElement.type(localKey.subEntity(), localKey.codim()).isTriangle())
                  DUNE_THROW(Dune::NotImplemented, "BernsteinBasis for 3D grids with k==3 is only implemented if the grid is a simplex grid");

                *it = {{ triangleOffset_ + static_cast<size_type>(gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim())) }};
                continue;
              }
          }

        if (dofDim==3)
          {
            if (dim==3)   // element dof -- any local numbering is fine
              {
                if (element.type().isTetrahedron())
                  {
                    *it = {{ tetrahedronOffset_ + dofsPerTetrahedron*static_cast<size_type>(gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                    continue;
                  }
                else if (element.type().isHexahedron())
                  {
                    *it = {{ hexahedronOffset_ + dofsPerHexahedron*static_cast<size_type>(gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                    continue;
                  }
                else if (element.type().isPrism())
                  {
                    *it = {{ prismOffset_ + dofsPerPrism*static_cast<size_type>(gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                    continue;
                  }
                else if (element.type().isPyramid())
                  {
                    *it = {{ pyramidOffset_ + dofsPerPyramid*static_cast<size_type>(gridIndexSet.subIndex(element,0,0)) + localKey.index() }};
                    continue;
                  }
                else
                  DUNE_THROW(Dune::NotImplemented, "3d elements have to be tetrahedra, hexahedra, prisms, or pyramids");
              } else
              DUNE_THROW(Dune::NotImplemented, "Grids of dimension larger than 3 are no supported");
          }
        DUNE_THROW(Dune::NotImplemented, "Grid contains elements not supported for the BernsteinBasis");
      }
    return it;
  }

protected:
  GridView gridView_;

  size_type vertexOffset_;
  size_type edgeOffset_;
  size_type triangleOffset_;
  size_type quadrilateralOffset_;
  size_type tetrahedronOffset_;
  size_type pyramidOffset_;
  size_type prismOffset_;
  size_type hexahedronOffset_;

};



template<typename GV, int k, typename R>
class BernsteinNode :
  public LeafBasisNode
{
  static constexpr int dim = GV::dimension;
  static const int maxSize = StaticPower<(k+1),GV::dimension>::power;

  using FiniteElementCache = typename Dune::BernsteinLocalFiniteElementCache<typename GV::ctype, R, dim, k>;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  BernsteinNode() :
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



#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,8)
template<typename GV, int k, class MI, typename R>
class BernsteinNodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = BernsteinPreBasis<GV, k, MI, R>;

  using Node = typename PreBasis::Node;

  BernsteinNodeIndexSet(const PreBasis& preBasis) :
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
#endif



namespace BasisFactory {

#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,9)
namespace Imp {

template<std::size_t k, typename R=double>
struct BernsteinPreBasisFactory
{
  static const std::size_t requiredMultiIndexSize = 1;

  template<class MultiIndex, class GridView>
  auto makePreBasis(const GridView& gridView) const
  {
    return BernsteinPreBasis<GridView, k, MultiIndex, R>(gridView);
  }
};

} // end namespace BasisFactory::Imp
#endif



/**
 * \brief Create a pre-basis factory that can build a Bernstein pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam k   The polynomial order of ansatz functions
 * \tparam R   The range type of the local basis
 */
template<std::size_t k, typename R=double>
auto bernstein()
{
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
  return [](const auto& gridView) {
    return BernsteinPreBasis<std::decay_t<decltype(gridView)>, k, R>(gridView);
  };
#else
  return Imp::BernsteinPreBasisFactory<k, R>();
#endif
}

} // end namespace BasisFactory



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Nodal basis of a scalar k-th-order Bernstein finite element space
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 *
 * All arguments passed to the constructor will be forwarded to the constructor
 * of BernsteinPreBasis.
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 * \tparam R The range type of the local basis
 */
template<typename GV, int k, typename R=double>
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
using BernsteinBasis = DefaultGlobalBasis<BernsteinPreBasis<GV, k, R>>;
#else
using BernsteinBasis = DefaultGlobalBasis<BernsteinPreBasis<GV, k, FlatMultiIndex<std::size_t>, R> >;
#endif



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BERNSTEINBASIS_HH
