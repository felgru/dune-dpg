// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKFACENODALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKFACENODALBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/std/final.hh>

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
//   PQkFaceNodeFactory
//   PQkFaceNodeIndexSet
//   PQkFaceNode
//
// The factory allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k, typename ST, typename TP>
class PQkFaceNode;

template<typename GV, int k, class MI, class TP, class ST>
class PQkFaceNodeIndexSet;

template<typename GV, int k, class MI, class ST>
class PQkFaceNodeFactory;



template<typename GV, int k, class MI, class ST>
class PQkFaceNodeFactory
{
  static const int dim = GV::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = ST;


  // Precompute the number of dofs per entity type
  const static int dofsPerEdge        = k+1;
  const static int dofsPerTriangle    = (k+1)*(k+2)/2;
  const static int dofsPerQuad        = (k+1)*(k+1);


  template<class TP>
  using Node = PQkFaceNode<GV, k, size_type, TP>;

  template<class TP>
  using IndexSet = PQkFaceNodeIndexSet<GV, k, MI, TP, ST>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 2>;

  /** \brief Constructor for a given grid view object */
  PQkFaceNodeFactory(const GridView& gv) :
    gridView_(gv)
  {}


  void initializeIndices()
  {
    edgeOffset_          = 0;
    if (dim==3)
    {
      DUNE_THROW(Dune::NotImplemented, "PQkFaceNodalBasis for 3D grids is not implemented");
      triangleOffset_      = 0;
      GeometryType triangle;
      triangle.makeTriangle();
      quadrilateralOffset_ = triangleOffset_        + dofsPerTriangle * gridView_.size(triangle);
    }
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return gridView_;
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
        DUNE_THROW(Dune::NotImplemented, "PQkFaceNodalBasis for 1D grids is not implemented");
        return 2*gridView_.size(dim)-2;
      case 2:
        return dofsPerEdge*gridView_.size(1);
      case 3:
      {
        DUNE_THROW(Dune::NotImplemented, "PQkFaceNodalBasis for 3D grids is not implemented");
        GeometryType triangle, quad;
        triangle.makeTriangle();
        quad.makeQuadrilateral();
        return dofsPerTriangle*gridView_.size(triangle) + dofsPerQuad*gridView_.size(quad);
      }
    }
    DUNE_THROW(Dune::NotImplemented, "No size method for " << dim << "d grids available yet!");
  }

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix prefix) const
  {
    if (prefix.size() == 0)
      return size();
    assert(false);
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
  const GridView gridView_;

  size_t edgeOffset_;
  size_t triangleOffset_;
  size_t quadrilateralOffset_;

};



template<typename GV, int k, typename ST, typename TP>
class PQkFaceNode :
  public GridFunctionSpaceBasisLeafNodeInterface<
    typename GV::template Codim<0>::Entity,
    typename Dune::PQkFaceLocalFiniteElementCache
             < typename GV::ctype
             , double, GV::dimension, k >::FiniteElementType,
    ST,
    TP>
{
  static const int dim = GV::dimension;
  static const int maxSize = 4*(k+1);

  typedef typename GV::template Codim<0>::Entity E;
  typedef typename Dune::PQkFaceLocalFiniteElementCache
                   <typename GV::ctype, double, dim, k> FiniteElementCache;
  typedef typename FiniteElementCache::FiniteElementType FE;

public:
  typedef GridFunctionSpaceBasisLeafNodeInterface<E,FE,ST,TP> Interface;
  typedef typename Interface::size_type size_type;
  typedef typename Interface::Element Element;
  typedef typename Interface::FiniteElement FiniteElement;
  typedef typename Interface::TreePath TreePath;

  PQkFaceNode(const TreePath& treePath) :
    Interface(treePath),
    finiteElement_(nullptr),
    element_(nullptr)
  {}

  //! Return current element, throw if unbound
  const Element& element() const DUNE_FINAL
  {
    return *element_;
  }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement() const DUNE_FINAL
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



template<typename GV, int k, class MI, class TP, class ST>
class PQkFaceNodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = ST;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = PQkFaceNodeFactory<GV, k, MI, ST>;

  using Node = typename NodeFactory::template Node<TP>;

  PQkFaceNodeIndexSet(const NodeFactory& nodeFactory) :
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
    return node_->finiteElement().size();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  const MultiIndex index(size_type i) const
  {
    Dune::LocalKey localKey = node_->finiteElement().localCoefficients().localKey(i);
    const auto& gridIndexSet = nodeFactory_->gridView().indexSet();
    const auto& element = node_->element();

    // The dimension of the entity that the current dof is related to
    size_t dofDim = dim - localKey.codim();
    if (dofDim==0) {  // vertex dof
      return {{ gridIndexSet.subIndex(element,localKey.subEntity(),dim) }};
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
        return {{ (flip)
          ? nodeFactory_->edgeOffset_ + (k+1)*gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()) + k-localKey.index()
              : nodeFactory_->edgeOffset_ + (k+1)*gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()) + localKey.index() }} ;
      }
    }

    if (dofDim==2)
    {
      if (dim==2)   // element dof -- any local numbering is fine
      {
        DUNE_THROW(Dune::NotImplemented, "faces have no elements of codimension 0");
      } else
      {
        DUNE_THROW(Dune::NotImplemented, "PQkFaceNodalBasis for 3D grids is not implemented");
      }
    }
    DUNE_THROW(Dune::NotImplemented, "Grid contains elements not supported for the PQkFaceNodalBasis");
  }

protected:
  const NodeFactory* nodeFactory_;

  const Node* node_;
};



namespace BasisBuilder {

namespace Imp {

template<std::size_t k>
struct PQkFaceNodeFactoryBuilder
{
  static const std::size_t requiredMultiIndexSize=1;

  template<class MultiIndex, class GridView, class size_type=std::size_t>
  auto build(const GridView& gridView)
    -> PQkFaceNodeFactory<GridView, k, MultiIndex, size_type>
  {
    return {gridView};
  }
};

} // end namespace BasisBuilder::Imp

template<std::size_t k>
Imp::PQkFaceNodeFactoryBuilder<k> pqFace()
{
  return{};
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
template<typename GV, int k, class ST = std::size_t>
using PQkFaceNodalBasis = DefaultGlobalBasis<PQkFaceNodeFactory<GV, k, FlatMultiIndex<ST>, ST> >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKFACENODALBASIS_HH
