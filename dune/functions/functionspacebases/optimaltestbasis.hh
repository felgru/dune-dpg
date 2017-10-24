// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_OPTIMALTESTBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_OPTIMALTESTBASIS_HH

#include <array>
#include <tuple>
#include <functional>
#include <memory>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/version.hh>

#include <dune/localfunctions/optimaltestfunctions/optimaltest.hh>
#include <dune/localfunctions/optimaltestfunctions/refinedoptimaltest.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/dpg/assemble_helper.hh>
#include <dune/dpg/testspace_coefficient_matrix.hh>
#include <dune/dpg/type_traits.hh>


namespace Dune {
namespace Functions {




//////////////////////////////////////////////////////////////////////////////
// Helper begin
//////////////////////////////////////////////////////////////////////////////


struct computeIndex
{
    computeIndex(size_t& space_index, size_t& index_result, bool& index_found)
    : space_index(space_index),
      index_result(index_result),
      index_found(index_found)
    {}

    template<class LIS>
    void operator()(LIS& localIndexSet) const
    {
      if (!index_found)
      {
        if (localIndexSet.size() > index_result)
        {
          index_found  = true;
          index_result = (localIndexSet.index(index_result))[0];
        }
        else
        {
          space_index  += 1;
          index_result -= localIndexSet.size();
        }
      }
    }

private:
    size_t& space_index;
    size_t& index_result;
    bool&   index_found;
};


//////////////////////////////////////////////////////////////////////////////
// Helper end
//////////////////////////////////////////////////////////////////////////////

/* some forward declarations for refined bases */

template<typename Element, typename ctype, int dim, int level>
class RefinedNode;

template<typename Element>
class UnrefinedNode
{
public:

  UnrefinedNode() :
    element_(nullptr)
  {}

protected:

  const Element* element_;
};

struct EmptyNodeFactoryConstants {};


// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   OptimalTestBasisNodeFactory
//   OptimalTestBasisNodeIndexSet
//   OptimalTestBasisNode
//
// The factory allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename TestspaceCoefficientMatrix, std::size_t testIndex, typename TP>
class OptimalTestBasisNode;

template<typename TestspaceCoefficientMatrix, std::size_t testIndex, class MI, class TP>
class OptimalTestBasisNodeIndexSet;

template<typename TestspaceCoefficientMatrix, std::size_t testIndex, class MI>
class OptimalTestBasisNodeFactory;


template<typename Space>
struct RefinementConstants : public EmptyNodeFactoryConstants {};

template<class GV, int level, int k>
struct RefinementConstants
       < DefaultGlobalBasis<
             PQkDGRefinedDGNodeFactory
             <GV, level, k, FlatMultiIndex<std::size_t> > > >
  : public PQkDGRefinedDGNodeFactory<GV, level, k, FlatMultiIndex<std::size_t> >
           ::RefinementConstants {};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex, class MI>
class OptimalTestBasisNodeFactory
  : public RefinementConstants<
               typename std::tuple_element<testIndex,
                     typename TestspaceCoefficientMatrix::TestSpaces
               >::type>
{
  static const int dim =
      TestspaceCoefficientMatrix::GridView::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = typename TestspaceCoefficientMatrix::GridView;
  using size_type = std::size_t;


  template<class TP>
  using Node = OptimalTestBasisNode<TestspaceCoefficientMatrix, testIndex, TP>;

  template<class TP>
  using IndexSet = OptimalTestBasisNodeIndexSet<TestspaceCoefficientMatrix, testIndex, MI, TP>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  using TestSearchSpaces = typename TestspaceCoefficientMatrix::TestSpaces;
  using SolutionSpaces = typename TestspaceCoefficientMatrix::SolutionSpaces;


  /** \brief Constructor for a given test coefficient matrix */
  OptimalTestBasisNodeFactory(TestspaceCoefficientMatrix& testCoeffMat) :
    testspaceCoefficientMatrix_(testCoeffMat)
  {}


  void initializeIndices()
  {
    Dune::detail::computeOffsets(
             globalOffsets,
             *testspaceCoefficientMatrix_.bilinearForm().getSolutionSpaces());
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return testspaceCoefficientMatrix_.gridView();
  }

  /**
   * \note We don't update anything here, as we expect the spaces referred
   * to in testspaceCoefficientMatrix_ to be updated seperately.
   */
  void update (const GridView&)
  {}

  template<class TP>
  Node<TP> node(const TP& tp) const
  {
    return Node<TP>{tp, testspaceCoefficientMatrix_};
  }

  template<class TP>
  IndexSet<TP> indexSet() const
  {
    return IndexSet<TP>{*this};
  }

  size_type size() const
  {
    return Hybrid::accumulate(*testspaceCoefficientMatrix_.bilinearForm()
                              .getSolutionSpaces(), 0,
                              [&](size_type acc, const auto& s) {
                                return acc + s.size();
                              });
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
    return Hybrid::accumulate(*testspaceCoefficientMatrix_.bilinearForm()
                              .getSolutionSpaces(), 0,
                              [&](size_type acc, const auto& s) {
                                return acc + s.nodeFactory().maxNodeSize();
                              });
  }

//protected:
  // TODO: store testspaceCoefficientMatrix_ by reference or by value?
  TestspaceCoefficientMatrix& testspaceCoefficientMatrix_;

  size_t globalOffsets[std::tuple_size<SolutionSpaces>::value];

};



template<typename TestspaceCoefficientMatrix, std::size_t testIndex, typename TP>
class OptimalTestBasisNode :
  public LeafBasisNode<std::size_t, TP>,
  public std::conditional
           < Dune::is_RefinedFiniteElement<
               typename std::tuple_element<testIndex,
                     typename TestspaceCoefficientMatrix::TestSpaces
               >::type>::value
           , RefinedNode< typename TestspaceCoefficientMatrix::GridView::
                                   template Codim<0>::Entity
                        , typename TestspaceCoefficientMatrix::GridView::ctype
                        , TestspaceCoefficientMatrix::GridView::dimension
                        , levelOfFE<typename std::tuple_element<testIndex,
                              typename TestspaceCoefficientMatrix::TestSpaces
                            >::type
                          >::value >
           , UnrefinedNode<typename TestspaceCoefficientMatrix::GridView::
                                   template Codim<0>::Entity >
           >::type
{
public:
  using SolutionSpaces = typename TestspaceCoefficientMatrix::SolutionSpaces;
  using TestSearchSpaces = typename TestspaceCoefficientMatrix::TestSpaces;

  using TestSearchSpace
          = typename std::tuple_element<testIndex,TestSearchSpaces>::type;

private:
  using GV = typename TestspaceCoefficientMatrix::GridView;
  static const int dim = GV::dimension;

  using Base = LeafBasisNode<std::size_t, TP>;
  using TestSearchFiniteElement
      = typename TestSearchSpace::LocalView::Tree::FiniteElement;

  using SolutionLocalViews = Dune::detail::getLocalViews_t<SolutionSpaces>;
  using TestLocalViews = Dune::detail::getLocalViews_t<TestSearchSpaces>;

  static const bool testSearchSpaceIsRefined
    = is_RefinedFiniteElement<TestSearchSpace>::value;

public:

  using size_type = std::size_t;
  using TreePath = TP;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename std::conditional
                          < testSearchSpaceIsRefined
                          , Dune::RefinedOptimalTestLocalFiniteElement
                                  <typename GV::ctype,
                                   double,
                                   GV::dimension,
                                   levelOfFE<TestSearchSpace>::value,
                                   RefinementConstants<TestSearchSpace>,
                                   TestSearchFiniteElement>
                          , Dune::OptimalTestLocalFiniteElement
                                  <typename GV::ctype,
                                   double,
                                   GV::dimension,
                                   TestSearchFiniteElement>
                          >::type;

  OptimalTestBasisNode(const TreePath& treePath,
                       TestspaceCoefficientMatrix& testCoeffMat) :
    Base(treePath),
    testspaceCoefficientMatrix(testCoeffMat),
    finiteElement_(nullptr),
    testSearchSpace_(nullptr),
    localViewsSolution_(Dune::detail::getLocalViews(*testCoeffMat.bilinearForm()
                                                .getSolutionSpaces())),
    localViewsTest(Dune::detail::getLocalViews(*testCoeffMat.bilinearForm()
                                                    .getTestSpaces()))
  {}

  //! Return current element, throw if unbound
  const Element& element() const
  {
    return *this->element_;
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
    using namespace Dune::detail;

    this->element_ = &e;
    bindLocalViews(localViewsTest, e);
    testSearchSpace_ =
        &(std::get<testIndex>(localViewsTest).tree().finiteElement());
    bindLocalViews(localViewsSolution_, e);

    testspaceCoefficientMatrix.bind(e);

    const size_t offset = computeOffset<testIndex>(localViewsTest);

    finiteElement_ = std::make_unique<FiniteElement>
                        (testspaceCoefficientMatrix.coefficientMatrix(),
                         *testSearchSpace_, offset);
    this->setSize(finiteElement_->size());
  }

  const SolutionLocalViews& localViewsSolution() const
  {
    return localViewsSolution_;
  }
protected:

  TestspaceCoefficientMatrix& testspaceCoefficientMatrix;
  //FiniteElementCache cache_;
  std::unique_ptr<FiniteElement> finiteElement_;
  const TestSearchFiniteElement* testSearchSpace_;
  SolutionLocalViews localViewsSolution_;
  TestLocalViews localViewsTest;
};



template<typename TestspaceCoefficientMatrix, std::size_t testIndex, class MI, class TP>
class OptimalTestBasisNodeIndexSet
{
  using GV = typename TestspaceCoefficientMatrix::GridView;
  enum {dim = GV::dimension};

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = OptimalTestBasisNodeFactory<TestspaceCoefficientMatrix, testIndex, MI>;

  using Node = typename NodeFactory::template Node<TP>;

  typedef typename TestspaceCoefficientMatrix::SolutionSpaces SolutionSpaces;
  typedef Dune::detail::getLocalIndexSets_t<SolutionSpaces>
      SolutionLocalIndexSets;

  OptimalTestBasisNodeIndexSet(const NodeFactory& nodeFactory) :
    nodeFactory_(&nodeFactory),
    solutionLocalIndexSets_(Dune::detail::getLocalIndexSets(
                      *nodeFactory.testspaceCoefficientMatrix_
                          .bilinearForm().getSolutionSpaces()))
  {}

  constexpr OptimalTestBasisNodeIndexSet(const OptimalTestBasisNodeIndexSet&)
          = default;

  OptimalTestBasisNodeIndexSet(OptimalTestBasisNodeIndexSet&& indexSet)
          = default;

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Node& node)
  {
    node_ = &node;

    Dune::detail::bindLocalIndexSets(solutionLocalIndexSets_,
                                     node.localViewsSolution());
  }

  /** \brief Unbind the view
   */
  void unbind()
  {
    node_ = nullptr;
    Hybrid::forEach(solutionLocalIndexSets_, [](auto& lis) {
          lis.unbind();
        });
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const
  {
    assert(node_ != nullptr);
    return node_->finiteElement().size();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
#if DUNE_VERSION_NEWER(DUNE_GRID,2,6)
  template<typename It>
  It indices(It it) const
  {
    assert(node_ != nullptr);
    for (size_type i = 0, end = this->size(); i < end; ++it, ++i)
    {
      size_t space_index = 0;
      size_t index_result = i;
      bool index_found=false;

      Hybrid::forEach(solutionLocalIndexSets_,
                      computeIndex(space_index, index_result, index_found));

      *it = {{ nodeFactory_->globalOffsets[space_index] + index_result }};
    }
    return it;
  }
#else
  MultiIndex index(size_type i) const
  {
    size_t space_index=0;
    size_t index_result=i;
    bool index_found=false;

    Hybrid::forEach(solutionLocalIndexSets_,
                    computeIndex(space_index, index_result, index_found));

    MultiIndex result;
    result[0]=(nodeFactory_->globalOffsets[space_index]+index_result);
    return result;
  }
#endif

protected:
  const NodeFactory* nodeFactory_;

  const Node* node_;

  SolutionLocalIndexSets solutionLocalIndexSets_;
};



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Nodal basis of a scalar k-th-order Lagrangean finite element space
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 *
 * \tparam TestspaceCoefficentMatrix  caches the computation of local
 *                                    optimal test spaces
 * \tparam testIndex  index of the optimal test space in the test space tuple
 */
template<typename TestspaceCoefficientMatrix, std::size_t testIndex = 0>
using OptimalTestBasis = DefaultGlobalBasis<OptimalTestBasisNodeFactory<TestspaceCoefficientMatrix, testIndex, FlatMultiIndex<std::size_t> > >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_OPTIMALTESTBASIS_HH
