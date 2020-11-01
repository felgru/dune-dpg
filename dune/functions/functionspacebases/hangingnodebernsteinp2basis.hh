// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HANGINGNODEBERNSTEINP2BASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HANGINGNODEBERNSTEINP2BASIS_HH

#include <array>
#include <optional>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/power.hh>
#include <dune/common/version.hh>

#include <dune/localfunctions/bernstein/pqkfactory.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/constrainedglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/nodes.hh>


namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   HangingNodeBernsteinP2PreBasis
//   HangingNodeBernsteinP2NodeIndexSet
//   HangingNodeBernsteinP2Node
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
template<typename GV>
class HangingNodeBernsteinP2Node;

template<typename GV, class MI>
class HangingNodeBernsteinP2NodeIndexSet;
#else
template<typename GV, typename TP>
class HangingNodeBernsteinP2Node;

template<typename GV, class MI, class TP>
class HangingNodeBernsteinP2NodeIndexSet;
#endif

template<typename GV, class MI>
class HangingNodeBernsteinP2PreBasis;



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
class HangingNodeBernsteinP2PreBasis
{
  static constexpr int dim = GV::dimension;

public:

  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

private:

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
  template<typename, class>
#else
  template<typename, class, class>
#endif
  friend class HangingNodeBernsteinP2NodeIndexSet;

  using FiniteElementCache = typename Dune::BernsteinLocalFiniteElementCache<typename GV::ctype, double, dim, 2>;

public:

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
  using Node = HangingNodeBernsteinP2Node<GV>;

  using IndexSet = HangingNodeBernsteinP2NodeIndexSet<GV, MI>;
#else
  //! Template mapping root tree path to type of created tree node
  template<class TP>
  using Node = HangingNodeBernsteinP2Node<GV, TP>;

  //! Template mapping root tree path to type of created tree node index set
  template<class TP>
  using IndexSet = HangingNodeBernsteinP2NodeIndexSet<GV, MI, TP>;
#endif

  //! Type used for global numbering of the basis vectors
  using MultiIndex = MI;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  //! Constructor for a given grid view object
  HangingNodeBernsteinP2PreBasis(const GridView& gv) :
    gridView_(gv),
    constraintWeights{{{   1.,    0.,    0.},
                       {1./2., 1./2.,    0.},
                       {1./4., 1./2., 1./4.}}}
  {
    static_assert(dim==2,
        "HangingNodeBernsteinP2GlobalBasis only defined in 2D!");
  }

  //! Initialize the global indices
  void initializeIndices()
  {
    const auto& gridIndexSet = gridView_.indexSet();
    std::vector<size_t> edgeDofs(gridView_.size(1), SIZE_MAX);
    std::vector<std::optional<std::array<size_t,3>>>
        edgeConstraints(edgeDofs.size());
    std::vector<std::optional<std::array<size_t,3>>>
        vertexConstraints(gridView_.size(2));
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
            vertexConstraints[midVertex]
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
                   "HangingNodeBernsteinP2PreBasis only implemented on triangles.");
      }
      const unsigned short numDofs = finiteElement.size();
      std::vector<MultiIndex> localToGlobal;
      localToGlobal.reserve(numDofs);
      std::vector<std::pair<size_t,size_t>> constraintOffsets;
      unsigned short preceedingUnconstrainedIndices = 0;
      for(unsigned short i = 0; i < numDofs; i++) {
        LocalKey localKey = finiteElement.localCoefficients().localKey(i);
        const auto dofCodim = localKey.codim();
        const size_type subIndex =
            gridIndexSet.subIndex(e, localKey.subEntity(), dofCodim);
        if(dofCodim == 2) { // dof is vertex
          if(vertexConstraints[subIndex]) {
            constraintOffsets.emplace_back(preceedingUnconstrainedIndices, 2);
            preceedingUnconstrainedIndices = 0;
            for(size_type idx : *vertexConstraints[subIndex])
            {
              localToGlobal.push_back(MultiIndex{idx});
            }
          } else {
            localToGlobal.push_back(MultiIndex{subIndex});
            ++preceedingUnconstrainedIndices;
          }
        } else { // dofCodim == 1, i.e. edge, since P2 does not have inner dofs
          if(edgeConstraints[subIndex]) {
            constraintOffsets.emplace_back(preceedingUnconstrainedIndices, 1);
            preceedingUnconstrainedIndices = 0;
            for(size_type idx : *edgeConstraints[subIndex])
            {
              localToGlobal.push_back(MultiIndex{idx});
            }
          } else {
            assert(edgeDofs[subIndex] != SIZE_MAX);
            localToGlobal.push_back(MultiIndex{edgeDofs[subIndex]});
            ++preceedingUnconstrainedIndices;
          }
        }
      }
      localToGlobal.shrink_to_fit();
      constraintOffsets.shrink_to_fit();
      const size_t elementIndex = gridIndexSet.index(e);
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

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
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
  IndexSet makeIndexSet() const
  {
    return IndexSet{*this};
  }
#else
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
#endif

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
  // First size_type is the offset of the constrained indices,
  // second size_type is index of constraintsWeights.
  std::vector<std::vector<std::pair<size_type, size_type>>> constraintIndicator;
  const std::array<std::array<double, 3>, 3> constraintWeights;

  size_type numberOfDofs;
};



#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
template<typename GV>
class HangingNodeBernsteinP2Node :
  public LeafBasisNode
#else
template<typename GV, typename TP>
class HangingNodeBernsteinP2Node :
  public LeafBasisNode<std::size_t, TP>
#endif
{
  static constexpr int dim = GV::dimension;

#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,7)
  using Base = LeafBasisNode<std::size_t,TP>;
#endif
  using FiniteElementCache = typename Dune::BernsteinLocalFiniteElementCache<typename GV::ctype, double, dim, 2>;

public:

  using size_type = std::size_t;
#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,7)
  using TreePath = TP;
#endif
  using Element = typename GV::template Codim<0>::Entity;
  using ElementSeed = typename Element::EntitySeed;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
  HangingNodeBernsteinP2Node() :
#else
  HangingNodeBernsteinP2Node(const TreePath& treePath) :
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
template<typename GV, class MI>
#else
template<typename GV, class MI, class TP>
#endif
class HangingNodeBernsteinP2NodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = HangingNodeBernsteinP2PreBasis<GV, MI>;

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
  using Node = typename PreBasis::Node;
#else
  using Node = typename PreBasis::template Node<TP>;
#endif

  using ConstraintWeights = std::array<double, 3>;

  HangingNodeBernsteinP2NodeIndexSet(const PreBasis& preBasis) :
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
    const size_t elementIndex = gridIndexSet.index(node_->element());
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
    return (*constraintIndicator_)[i].first;
  }

  const ConstraintWeights& constraintWeights(size_type i) const
  {
    assert(node_ != nullptr);
    return preBasis_->constraintWeights[(*constraintIndicator_)[i].second];
  }

protected:
  const PreBasis* preBasis_;

  const Node* node_;

  const std::vector<MultiIndex>* indicesLocalGlobal_;
  const std::vector<std::pair<size_t,size_t>>* constraintIndicator_;
};



namespace BasisFactory {

namespace Imp {

struct HangingNodeBernsteinP2PreBasisFactory
{
  static const std::size_t requiredMultiIndexSize = 1;

  template<class MultiIndex, class GridView>
  auto makePreBasis(const GridView& gridView) const
  {
    return HangingNodeBernsteinP2PreBasis<GridView, MultiIndex>(gridView);
  }
};

} // end namespace BasisFactory::Imp



/**
 * \brief Create a pre-basis builder that can build a hanging node P_2 pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 */
auto hangingNodeP2()
{
  return Imp::HangingNodeBernsteinP2PreBasisFactory();
}

} // end namespace BasisFactory



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Nodal basis of a scalar 2nd-order Bernstein finite element space
 *         with hanging nodes
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \note This only works for 2d grids.
 *
 * All arguments passed to the constructor will be forwarded to the constructor
 * of HangingNodeBernsteinP2PreBasis.
 *
 * \tparam GV The GridView that the space is defined on
 */
template<typename GV>
using HangingNodeBernsteinP2Basis = ConstrainedGlobalBasis<HangingNodeBernsteinP2PreBasis<GV, FlatMultiIndex<std::size_t>> >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HANGINGNODEBERNSTEINP2BASIS_HH
