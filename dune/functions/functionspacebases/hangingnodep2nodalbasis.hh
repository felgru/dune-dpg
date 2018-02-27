// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HANGINGNODEP2NODALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HANGINGNODEP2NODALBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/power.hh>
#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
#include <dune/common/std/optional.hh>
#else
#include <dune/functions/common/optional.hh>
#endif

#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/constrainedglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/nodes.hh>


namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   HangingNodeP2PreBasis
//   HangingNodeP2NodeIndexSet
//   HangingNodeP2Node
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, typename TP>
class HangingNodeP2Node;

template<typename GV, class MI, class TP>
class HangingNodeP2NodeIndexSet;

template<typename GV, class MI>
class HangingNodeP2PreBasis;



/**
 * \brief A pre-basis for PQ-lagrange bases of order 2 with hanging nodes
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV  The grid view that the FE basis is defined on
 * \tparam MI  Type to be used for multi-indices
 *
 * \note This only works on 2d grids
 */
template<typename GV, class MI>
class HangingNodeP2PreBasis
{
  static constexpr int dim = GV::dimension;

public:

  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

private:

  template<typename, class, class>
  friend class HangingNodeP2NodeIndexSet;

  using FiniteElementCache = typename Dune::PQkLocalFiniteElementCache<typename GV::ctype, double, dim, 2>;

public:

  //! Template mapping root tree path to type of created tree node
  template<class TP>
  using Node = HangingNodeP2Node<GV, TP>;

  //! Template mapping root tree path to type of created tree node index set
  template<class TP>
  using IndexSet = HangingNodeP2NodeIndexSet<GV, MI, TP>;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = MI;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  //! Constructor for a given grid view object
  HangingNodeP2PreBasis(const GridView& gv) :
    gridView_(gv),
    constraint{3./8, 3./4, -1./8}
  {
    static_assert(dim==2, "HangingNodeP2GlobalBasis only defined in 2D!");
  }

  //! Initialize the global indices
  void initializeIndices()
  {
    const auto& gridIndexSet = gridView_.indexSet();
    std::vector<size_t> edgeDofs(gridView_.size(1), SIZE_MAX);
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
    std::vector<Dune::Std::optional<std::array<size_t,3>>>
#else
    std::vector<Optional<std::array<size_t,3>>>
#endif
        edgeConstraints(edgeDofs.size());
    size_t nextEdgeDof = gridView_.size(dim); // edges start after vertices
    for(const auto& e : elements(gridView_)) {
      for (auto&& intersection : intersections(gridView_, e)) {
        if (!intersection.conforming()) {
          if(intersection.inside().level() < intersection.outside().level()) {
            // inside dominates outside with one level difference
            assert(intersection.outside().level()
                 - intersection.inside().level() == 1);
            const auto edgeInside = intersection.indexInInside();
            const auto edgeOutside = intersection.indexInOutside();
            std::array<size_t, 2> verticesInside {0, 0};
            switch(edgeInside) {
              case 0:
                verticesInside
                  = {gridIndexSet.subIndex(e, 0, dim),
                     gridIndexSet.subIndex(e, 1, dim)};
                break;
              case 1:
                verticesInside
                  = {gridIndexSet.subIndex(e, 0, dim),
                     gridIndexSet.subIndex(e, 2, dim)};
                break;
              case 2:
                verticesInside
                  = {gridIndexSet.subIndex(e, 1, dim),
                     gridIndexSet.subIndex(e, 2, dim)};
                break;
            }
            const auto outsideElement = intersection.outside();
            std::array<size_t, 2> verticesOutside {0, 0};
            switch(edgeOutside) {
              case 0:
                verticesOutside
                  = {gridIndexSet.subIndex(outsideElement, 0, dim),
                     gridIndexSet.subIndex(outsideElement, 1, dim)};
                break;
              case 1:
                verticesOutside
                  = {gridIndexSet.subIndex(outsideElement, 0, dim),
                     gridIndexSet.subIndex(outsideElement, 2, dim)};
                break;
              case 2:
                verticesOutside
                  = {gridIndexSet.subIndex(outsideElement, 1, dim),
                     gridIndexSet.subIndex(outsideElement, 2, dim)};
                break;
            }
            size_t commonVertex, midVertex, farVertex;
            if(verticesInside[0] == verticesOutside[0]) {
              commonVertex = verticesInside[0];
              midVertex = verticesOutside[1];
              farVertex = verticesInside[1];
            } else if(verticesInside[0] == verticesOutside[1]) {
              commonVertex = verticesInside[0];
              midVertex = verticesOutside[0];
              farVertex = verticesInside[1];
            } else if(verticesInside[1] == verticesOutside[0]) {
              commonVertex = verticesInside[1];
              midVertex = verticesOutside[1];
              farVertex = verticesInside[0];
            } else { // if(verticesInside[1] == verticesOutside[1])
              commonVertex = verticesInside[1];
              midVertex = verticesOutside[0];
              farVertex = verticesInside[0];
            }
            edgeDofs[gridIndexSet.subIndex(e, edgeInside, 1)] = midVertex;
            edgeConstraints[gridIndexSet.subIndex(outsideElement,edgeOutside,1)]
                = std::array<size_t,3>{commonVertex, midVertex, farVertex};
          } else {
            // Nothing to be done here, as the edge DoF from this face
            // is constrained by the larger neighboring element.
          }
        } else {
          size_t& edgeDof
            = edgeDofs[gridIndexSet.subIndex(e,intersection.indexInInside(),1)];
          if(edgeDof == SIZE_MAX)
            edgeDof = nextEdgeDof++;
        }
      }
    }

    numberOfDofs = nextEdgeDof;

    indicesLocalGlobal.resize(gridView_.size(1));
    constraintIndicator.resize(gridView_.size(1));
    FiniteElementCache feCache;
    for(const auto& e : elements(gridView_)) {
      auto& finiteElement = feCache.get(e.type());
      if(!e.type().isTriangle ()) {
        DUNE_THROW(Dune::NotImplemented,
                   "HangingNodeP2PreBasis only implemented on triangles.");
      }
      const unsigned short numDofs = finiteElement.size();
      std::vector<MultiIndex> localToGlobal;
      localToGlobal.reserve(numDofs);
      std::vector<size_t> constraintOffsets;
      unsigned short preceedingUnconstrainedIndices = 0;
      for(unsigned short i = 0; i < numDofs; i++) {
        LocalKey localKey = finiteElement.localCoefficients().localKey(i);
        const auto dofCodim = localKey.codim();
        const size_type subIndex =
            gridIndexSet.subIndex(e, localKey.subEntity(), dofCodim);
        if(dofCodim == 2) {
          localToGlobal.emplace_back(MultiIndex{subIndex});
          ++preceedingUnconstrainedIndices;
        } else { // dofCodim == 1, since P2 does not have inner dofs
          if(edgeConstraints[subIndex]) {
            constraintOffsets.push_back(preceedingUnconstrainedIndices);
            preceedingUnconstrainedIndices = 0;
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
            for(size_type idx : *edgeConstraints[subIndex])
#else
            for(size_type idx : edgeConstraints[subIndex].value())
#endif
            {
              localToGlobal.emplace_back(MultiIndex{idx});
            }
          } else {
            assert(edgeDofs[subIndex] != SIZE_MAX);
            localToGlobal.emplace_back(MultiIndex{edgeDofs[subIndex]});
            ++preceedingUnconstrainedIndices;
          }
        }
      }
      localToGlobal.shrink_to_fit();
      constraintOffsets.shrink_to_fit();
      size_t elementIndex = gridIndexSet.subIndex(e, 0, 0);
      std::swap(indicesLocalGlobal[elementIndex], localToGlobal);
      std::swap(constraintIndicator[elementIndex], constraintOffsets);
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
   * \brief Create tree node with given root tree path
   *
   * \tparam TP Type of root tree path
   * \param tp Root tree path
   *
   * By passing a non-trivial root tree path this can be used
   * to create a node suitable for being placed in a tree at
   * the position specified by the root tree path.
   */
  template<class TP>
  Node<TP> node(const TP& tp) const
  {
    return Node<TP>{tp};
  }

  /**
   * \brief Create tree node index set with given root tree path
   *
   * \tparam TP Type of root tree path
   * \param tp Root tree path
   *
   * Create an index set suitable for the tree node obtained
   * by node(tp).
   */
  template<class TP>
  IndexSet<TP> indexSet() const
  {
    return IndexSet<TP>{*this};
  }

  //! Same as size(prefix) with empty prefix
  size_type size() const
  {
    return numberOfDofs;
  }

  //! Return number of possible values for next position in multi index
  size_type size(const SizePrefix prefix) const
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
    return StaticPower<(2+1),GV::dimension>::power;
  }

protected:
  GridView gridView_;

  std::vector<std::vector<MultiIndex>> indicesLocalGlobal;
  std::vector<std::vector<size_type>> constraintIndicator;
  const std::array<double, 3> constraint;

  size_type numberOfDofs;
};



template<typename GV, typename TP>
class HangingNodeP2Node :
  public LeafBasisNode<std::size_t, TP>
{
  static constexpr int dim = GV::dimension;

  using Base = LeafBasisNode<std::size_t,TP>;
  using FiniteElementCache = typename Dune::PQkLocalFiniteElementCache<typename GV::ctype, double, dim, 2>;

public:

  using size_type = std::size_t;
  using TreePath = TP;
  using Element = typename GV::template Codim<0>::Entity;
  using ElementSeed = typename Element::EntitySeed;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  HangingNodeP2Node(const TreePath& treePath) :
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


template<typename GV, class MI, class TP>
class HangingNodeP2NodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = HangingNodeP2PreBasis<GV, MI>;
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,6))
  using NodeFactory = PreBasis;
#endif

  using Node = typename PreBasis::template Node<TP>;

  using ConstraintWeights = std::array<double, 3>;

  HangingNodeP2NodeIndexSet(const PreBasis& preBasis) :
    preBasis_(&preBasis),
    node_(nullptr),
    indicesLocalGlobal_(nullptr)
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Node& node)
  {
    node_ = &node;
    const auto& gridIndexSet = preBasis_->gridView().indexSet();
    size_t elementIndex = gridIndexSet.subIndex(node_->element(), 0, 0);
    indicesLocalGlobal_ = &preBasis_->indicesLocalGlobal[elementIndex];
    constraintIndicator_ = &preBasis_->constraintIndicator[elementIndex];
  }

  /** \brief Unbind the view
   */
  void unbind()
  {
    node_ = nullptr;
    indicesLocalGlobal_ = nullptr;
    constraintIndicator_ = nullptr;
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const
  {
    assert(node_ != nullptr);
    return node_->finiteElement().size();
  }

  // TODO: That might not be the most user friendly interface!
  const std::vector<MultiIndex>& indicesLocalGlobal() const
  {
    assert(node_ != nullptr);
    return *indicesLocalGlobal_;
  }

  size_type constraintsSize() const
  {
    assert(node_ != nullptr);
    return constraintIndicator_->size();
  }

  size_type constraintOffset(size_type i) const
  {
    assert(node_ != nullptr);
    return (*constraintIndicator_)[i];
  }

  const ConstraintWeights& constraintWeights(size_type i) const
  {
    assert(node_ != nullptr);
    return preBasis_->constraint;
  }

protected:
  const PreBasis* preBasis_;

  const Node* node_;

  const std::vector<MultiIndex>* indicesLocalGlobal_;
  const std::vector<size_t>* constraintIndicator_;
};



namespace BasisBuilder {

namespace Imp {

struct HangingNodeP2PreBasisFactory
{
  static const std::size_t requiredMultiIndexSize = 1;

  template<class MultiIndex, class GridView>
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,6)
  auto makePreBasis(const GridView& gridView) const
#else
  auto build(const GridView& gridView) const
#endif
  {
    return HangingNodeP2PreBasis<GridView, MultiIndex>(gridView);
  }
};

} // end namespace BasisBuilder::Imp



/**
 * \brief Create a pre-basis builder that can build a hanging node P_2 pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 */
auto hangingNodeP2()
{
  return Imp::HangingNodeP2PreBasisFactory();
}

} // end namespace BasisBuilder



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Nodal basis of a scalar 2nd-order Lagrangean finite element space
 *         with hanging nodes
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \note This only works for 2d grids.
 *
 * All arguments passed to the constructor will be forwarded to the constructor
 * of HangingNodeP2PreBasis.
 *
 * \tparam GV The GridView that the space is defined on
 */
template<typename GV>
using HangingNodeP2NodalBasis = ConstrainedGlobalBasis<HangingNodeP2PreBasis<GV, FlatMultiIndex<std::size_t>> >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HANGINGNODEP2NODALBASIS_HH
