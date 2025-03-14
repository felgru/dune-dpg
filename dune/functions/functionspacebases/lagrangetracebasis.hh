// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGETRACEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGETRACEBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/math.hh>

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

template<typename GV, int k, typename R=double>
class LagrangeTraceNode;

template<typename GV, int k, typename R=double>
class LagrangeTracePreBasis;



template<typename GV, int k, typename R>
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


  using Node = LagrangeTraceNode<GV, k, R>;

  static constexpr size_type maxMultiIndexSize = 1;
  static constexpr size_type minMultiIndexSize = 1;
  static constexpr size_type multiIndexBufferSize = 1;

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

  Node makeNode() const
  {
    return Node{};
  }

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
  template<class SizePrefix>
  size_type size(const SizePrefix& prefix) const
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
    return power(k+1, dim) - power(k-1, dim);
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



template<typename GV, int k, typename R>
class LagrangeTraceNode :
  public LeafBasisNode
{
  static constexpr int dim = GV::dimension;

  using FiniteElementCache = typename Dune::PQkTraceLocalFiniteElementCache
                                        <typename GV::ctype, R, dim, k>;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  LagrangeTraceNode() :
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



namespace BasisFactory {

template<std::size_t k, typename R=double>
auto lagrangeTrace()
{
  return [](const auto& gridView) {
    return LagrangeTracePreBasis<std::decay_t<decltype(gridView)>, k, R>(gridView);
  };
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
 * \tparam R The range type of the local basis
 */
template<typename GV, int k, typename R=double>
using LagrangeTraceBasis = DefaultGlobalBasis<LagrangeTracePreBasis<GV, k, R>>;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGETRACEBASIS_HH
