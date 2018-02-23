// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKFACENODALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKFACENODALBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/version.hh>

#include <dune/localfunctions/lagrange/pqkfacefactory.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>


namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   PQkFacePreBasis
//   PQkFaceNodeIndexSet
//   PQkFaceNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k, typename TP>
class PQkFaceNode;

template<typename GV, int k, class MI, class TP>
class PQkFaceNodeIndexSet;

template<typename GV, int k, class MI>
class PQkFacePreBasis;



template<typename GV, int k, class MI>
class PQkFacePreBasis
{
  static constexpr int dim = GV::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = std::size_t;


  // Precompute the number of dofs per entity type
  static constexpr int dofsPerEdge     = k+1;
  static constexpr int dofsPerTriangle = (k+1)*(k+2)/2;
  static constexpr int dofsPerQuad     = (k+1)*(k+1);


  template<class TP>
  using Node = PQkFaceNode<GV, k, TP>;

  template<class TP>
  using IndexSet = PQkFaceNodeIndexSet<GV, k, MI, TP>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  /** \brief Constructor for a given grid view object */
  PQkFacePreBasis(const GridView& gv) :
    gridView_(gv)
  {}


  void initializeIndices()
  {
    edgeOffset_          = 0;
    if (dim==3)
    {
      DUNE_THROW(Dune::NotImplemented, "PQkFaceNodalBasis for 3D grids is not implemented");
      triangleOffset_      = 0;
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
      quadrilateralOffset_ = triangleOffset_
                           + dofsPerTriangle
                             * gridView_.size(GeometryTypes::triangle);
#else
      GeometryType triangle;
      triangle.makeTriangle();
      quadrilateralOffset_ = triangleOffset_
                           + dofsPerTriangle * gridView_.size(triangle);
#endif
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
        DUNE_THROW(Dune::NotImplemented,
                   "PQkFaceNodalBasis for 1D grids is not implemented");
        return 2*gridView_.size(dim)-2;
      case 2:
        return dofsPerEdge*gridView_.size(1);
      case 3:
      {
        DUNE_THROW(Dune::NotImplemented,
                   "PQkFaceNodalBasis for 3D grids is not implemented");
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
    return 4*(k+1);
  }

//protected:
  GridView gridView_;

  size_type edgeOffset_;
  size_type triangleOffset_;
  size_type quadrilateralOffset_;

};



template<typename GV, int k, typename TP>
class PQkFaceNode :
  public LeafBasisNode<std::size_t, TP>
{
  static constexpr int dim = GV::dimension;

  using Base = LeafBasisNode<std::size_t, TP>;
  using FiniteElementCache = typename Dune::PQkFaceLocalFiniteElementCache
                                        <typename GV::ctype, double, dim, k>;

public:

  using size_type = std::size_t;
  using TreePath = TP;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  PQkFaceNode(const TreePath& treePath) :
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



template<typename GV, int k, class MI, class TP>
class PQkFaceNodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = PQkFacePreBasis<GV, k, MI>;
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,6))
  using NodeFactory = PreBasis;
#endif

  using Node = typename PreBasis::template Node<TP>;

  PQkFaceNodeIndexSet(const PreBasis& preBasis) :
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
      const Dune::LocalKey localKey
          = node_->finiteElement().localCoefficients().localKey(i);
      // The dimension of the entity that the current dof is related to
      const size_t dofDim = dim - localKey.codim();
      if (dofDim==0) {  // vertex dof
        *it = {{ gridIndexSet.subIndex(element,localKey.subEntity(),dim) }};
        continue;
      }

      if (dofDim==1)
      {  // edge dof
        if (dim==1)   // element dof -- any local numbering is fine
        {
          DUNE_THROW(Dune::NotImplemented,
              "faces have no elements of codimension 0");
        }
        else
        {
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
          const Dune::ReferenceElement<double,dim> refElement
              = Dune::referenceElement<double,dim>(element.type());
#else
          const Dune::ReferenceElement<double,dim>& refElement
              = Dune::ReferenceElements<double,dim>::general(element.type());
#endif

          // we have to reverse the numbering if the local triangle edge is
          // not aligned with the global edge
          size_t v0 = gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),0,dim),dim);
          size_t v1 = gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),1,dim),dim);
          bool flip = (v0 > v1);
          *it = {{ (flip)
                   ? preBasis_->edgeOffset_
                     + (k+1)*gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim())
                     + k-localKey.index()
                   : preBasis_->edgeOffset_
                     + (k+1)*gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim())
                     + localKey.index() }};
          continue;
        }
      }

      if (dofDim==2)
      {
        if (dim==2)   // element dof -- any local numbering is fine
        {
          DUNE_THROW(Dune::NotImplemented,
                     "faces have no elements of codimension 0");
        } else
        {
          DUNE_THROW(Dune::NotImplemented,
                     "PQkFaceNodalBasis for 3D grids is not implemented");
        }
      }
      DUNE_THROW(Dune::NotImplemented,
          "Grid contains elements not supported for the PQkFaceNodalBasis");
    }
    return it;
  }
#else
  MultiIndex index(size_type i) const
  {
    Dune::LocalKey localKey = node_->finiteElement().localCoefficients().localKey(i);
    const auto& gridIndexSet = preBasis_->gridView().indexSet();
    const auto& element = node_->element();

    // The dimension of the entity that the current dof is related to
    size_t dofDim = dim - localKey.codim();
    if (dofDim==0) {  // vertex dof
      return { gridIndexSet.subIndex(element,localKey.subEntity(),dim) };
    }

    if (dofDim==1)
    {  // edge dof
      if (dim==1)   // element dof -- any local numbering is fine
      {
        DUNE_THROW(Dune::NotImplemented, "faces have no elements of codimension 0");
      }
      else
      {
        const Dune::ReferenceElement<double,dim>& refElement
            = Dune::ReferenceElements<double,dim>::general(element.type());

        // we have to reverse the numbering if the local triangle edge is
        // not aligned with the global edge
        size_t v0 = gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),0,dim),dim);
        size_t v1 = gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),1,dim),dim);
        bool flip = (v0 > v1);
        return { (flip)
          ? preBasis_->edgeOffset_ + (k+1)*gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()) + k-localKey.index()
              : preBasis_->edgeOffset_ + (k+1)*gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()) + localKey.index() };
      }
    }

    if (dofDim==2)
    {
      if (dim==2)   // element dof -- any local numbering is fine
      {
        DUNE_THROW(Dune::NotImplemented,
                   "faces have no elements of codimension 0");
      } else
      {
        DUNE_THROW(Dune::NotImplemented,
                   "PQkFaceNodalBasis for 3D grids is not implemented");
      }
    }
    DUNE_THROW(Dune::NotImplemented,
        "Grid contains elements not supported for the PQkFaceNodalBasis");
  }
#endif

protected:
  const PreBasis* preBasis_;

  const Node* node_;
};



namespace BasisBuilder {

namespace Imp {

template<std::size_t k>
struct PQkFacePreBasisFactory
{
  static const std::size_t requiredMultiIndexSize = 1;

  template<class MultiIndex, class GridView, class size_type=std::size_t>
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,6)
  auto makePreBasis(const GridView& gridView) const
#else
  auto build(const GridView& gridView) const
#endif
  {
    return PQkFacePreBasis<GridView, k, MultiIndex>(gridView);
  }
};

} // end namespace BasisBuilder::Imp

template<std::size_t k>
auto pqFace()
{
  return Imp::PQkFacePreBasisFactory<k>();
}

} // end namespace BasisBuilder



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Nodal basis of a scalar k-th-order Lagrangean finite element space
 *         defined on the faces of a cell
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
using PQkFaceNodalBasis = DefaultGlobalBasis<PQkFacePreBasis<GV, k, FlatMultiIndex<std::size_t> > >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKFACENODALBASIS_HH
