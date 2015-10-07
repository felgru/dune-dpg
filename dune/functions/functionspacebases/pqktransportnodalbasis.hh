// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKTRANSPORTBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKTRANSPORTBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/std/final.hh>

#include <dune/localfunctions/lagrange/pk2d.hh>
#include <dune/localfunctions/lagrange/qk.hh>
#include <dune/localfunctions/transportfunctions/pqktransportfactory.hh>
#include <dune/localfunctions/transportfunctions/qktransport.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/gridviewfunctionspacebasis.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>


namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   PQkTransportFactory
//   PQkTransportIndexSet
//   PQkTransport
//
// The factory allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k, typename ST, typename TP>
class PQkTransportNode;

template<typename GV, int k, class MI, class TP, class ST>
class PQkTransportIndexSet;

template<typename GV, int k, class MI, class ST>
class PQkTransportFactory;



template<typename GV, int k, class MI, class ST>
class PQkTransportFactory
{
  static const int dim = GV::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = ST;


  template<class TP>
  using Node = PQkTransportNode<GV, k, size_type, TP>;

  template<class TP>
  using IndexSet = PQkTransportIndexSet<GV, k, MI, TP, ST>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  /** \brief Constructor for a given grid view object */
  PQkTransportFactory(const GridView& gv,
                      const FieldVector<double,dim>& transportDirection) :
    gridView_(gv),
    beta_(transportDirection)
  {}


  void initializeIndices()
  {}

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return gridView_;
  }

  template<class TP>
  Node<TP> node(const TP& tp) const
  {
    return Node<TP>{tp, beta_};
  }

  template<class TP>
  IndexSet<TP> indexSet() const
  {
    return IndexSet<TP>{*this};
  }

  size_type size() const
  {
    GeometryType quad;
    quad.makeQuadrilateral();
    if(   beta_[0] == 0
       || beta_[1] == 0
       || beta_[0] == beta_[1]
      )
      return sizeQ * gridView_.size(quad);
    else
      return dofsPerQuad * gridView_.size(quad);
    //DUNE_THROW(Dune::NotImplemented, "No size method for " << dim << "d grids available yet!");
  }

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix prefix) const
  {
    if (prefix.size() == 0)
      return size();
    if (prefix.size() == 1)
      return 0;
    assert(false);
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

  FieldVector<double,dim> transportDirection() const
  {
    return beta_;
  }

//protected:
  const GridView gridView_;

  Pk2DLocalFiniteElement<double,double,k> localFiniteElementP;
  QkLocalFiniteElement<double,double,dim,k> localFiniteElementQ;
  // Precompute the number of dofs per entity type
  const unsigned int sizeP = localFiniteElementP.size();
  const unsigned int sizeQ = localFiniteElementQ.size();
  const int dofsPerQuad    = (sizeP+sizeQ-(k+1));
  FieldVector<double,dim> beta_;

};



template<typename GV, int k, typename ST, typename TP>
class PQkTransportNode :
  public GridFunctionSpaceBasisLeafNodeInterface<
    typename GV::template Codim<0>::Entity,
    typename Dune::PQkTransportLocalFiniteElementCache<typename GV::ctype,
                                   double, GV::dimension,k>::FiniteElementType,
    ST,
    TP>
{
  static const int dim = GV::dimension;

  typedef typename GV::template Codim<0>::Entity E;
  typedef typename Dune::PQkTransportLocalFiniteElementCache<typename GV::ctype, double, dim,k> FiniteElementCache;
  typedef typename FiniteElementCache::FiniteElementType FE;

public:
  typedef GridFunctionSpaceBasisLeafNodeInterface<E,FE,ST,TP> Interface;
  typedef typename Interface::size_type size_type;
  typedef typename Interface::Element Element;
  typedef typename Interface::FiniteElement FiniteElement;
  typedef typename Interface::TreePath TreePath;

  PQkTransportNode(const TreePath& treePath,
                   const FieldVector<double,dim>& transportDirection) :
    Interface(treePath),
    finiteElement_(nullptr),
    element_(nullptr),
    beta_(transportDirection)
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
    FieldVector<typename GV::ctype,dim> normalized, unnormalized;
    unnormalized=beta_;
    if(unnormalized[0]!=0 && unnormalized[1]!=0){

      double factor=unnormalized[1]/unnormalized[0];
      if(factor<=1){
        normalized[0]=1.;
        normalized[1]=factor;
      }
      else if(factor>1){
        normalized[0]=1./factor;
        normalized[1]=1.;
      }

    }
    if(unnormalized[0]==0){
      normalized[0]=0.;
      normalized[1]=1.;
    }
    if(unnormalized[1]==0){
      normalized[0]=1.;
      normalized[1]=0.;
    }
    finiteElement_ = &(cache_.get(element_->type(),normalized));
    this->setSize(finiteElement_->size());
  }

protected:

  FiniteElementCache cache_;
  const FiniteElement* finiteElement_;
  const Element* element_;
  FieldVector<double,dim> beta_;
};



template<typename GV, int k, class MI, class TP, class ST>
class PQkTransportIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = ST;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = PQkTransportFactory<GV, k, MI, ST>;

  using Node = typename NodeFactory::template Node<TP>;

  PQkTransportIndexSet(const NodeFactory& nodeFactory) :
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
  MultiIndex index(size_type i) const
  {
    const auto& gridIndexSet = nodeFactory_->gridView().indexSet();
    const auto& element = node_->element();

    const auto& beta = nodeFactory_->transportDirection();


    if(   beta[0] == 0
       || beta[1] == 0
       || beta[0] == beta[1]
      )
      return {nodeFactory_->sizeQ*gridIndexSet.subIndex(element,0,0) + i};
    else
      return {nodeFactory_->dofsPerQuad*gridIndexSet.subIndex(element,0,0) + i};
    //DUNE_THROW(Dune::NotImplemented, "2d elements have to be triangles or quadrilaterals");

    //DUNE_THROW(Dune::NotImplemented, "No index method for " << dim << "d grids available yet!");
  }

protected:
  const NodeFactory* nodeFactory_;

  const Node* node_;
};



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Nodal basis of a scalar k-th-order Lagrangean finite element space
 *         aligned with the transport direction
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - Grids must be 1d or 2d
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int k, class ST = std::size_t>
using PQkTransportBasis = DefaultGlobalBasis<PQkTransportFactory<GV, k, FlatMultiIndex<ST>, ST> >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKTRANSPORTBASIS_HH
