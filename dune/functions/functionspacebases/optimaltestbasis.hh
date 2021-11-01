// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_OPTIMALTESTBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_OPTIMALTESTBASIS_HH

#include <array>
#include <tuple>
#include <memory>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/version.hh>

#include <dune/localfunctions/optimaltestfunctions/optimaltest.hh>
#include <dune/localfunctions/optimaltestfunctions/refinedoptimaltest.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/dpg/assemble_helper.hh>
#include <dune/dpg/testspace_coefficient_matrix.hh>
#include <dune/dpg/type_traits.hh>

#include <boost/hana.hpp>


namespace Dune {
namespace Functions {

/* some forward declarations for refined bases */

template<class PB>
class DefaultGlobalBasis;

template<class PB>
class RefinedGlobalBasis;

template<typename Element, typename ctype, int dim, int level>
class RefinedNode;

template<typename Element>
class UnrefinedNode
{
public:

  UnrefinedNode() :
    element_(nullptr)
  {}

  //! Return current element, throw if unbound
  const Element& element() const
  {
    return *element_;
  }

protected:

  const Element* element_;
};

struct EmptyPreBasisConstants {};


// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   OptimalTestBasisPreBasis
//   OptimalTestBasisNodeIndexSet
//   OptimalTestBasisNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
class OptimalTestBasisNode;

#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,8)
template<typename TestspaceCoefficientMatrix, std::size_t testIndex, class MI>
class OptimalTestBasisNodeIndexSet;
#endif

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
#else
template<typename TestspaceCoefficientMatrix, std::size_t testIndex, class MI>
#endif
class OptimalTestBasisPreBasis;


template<typename Space>
struct RefinementConstants : public EmptyPreBasisConstants {};

template<class GV, int level, int k>
struct RefinementConstants
       < DefaultGlobalBasis<
             LagrangeDGRefinedDGPreBasis
             <GV, level, k, FlatMultiIndex<std::size_t> > > >
  : public LagrangeDGRefinedDGPreBasis<GV, level, k, FlatMultiIndex<std::size_t> >
           ::RefinementConstants {};

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
#else
template<typename TestspaceCoefficientMatrix, std::size_t testIndex, class MI>
#endif
class OptimalTestBasisPreBasis
  : public RefinementConstants<
               typename std::tuple_element<testIndex,
                     typename TestspaceCoefficientMatrix::TestSpaces
               >::type>
{
  static constexpr int dim =
      TestspaceCoefficientMatrix::GridView::dimension;

  using SolutionSpaces = typename TestspaceCoefficientMatrix::SolutionSpaces;
  using SolutionLocalViews = Dune::detail::getLocalViews_t<SolutionSpaces>;

  template<typename It, typename LocalView>
  static It computeIndices(It it, const LocalView& localView,
                           size_t globalOffset)
  {
    for (size_type i = 0, end = localView.size(); i < end; ++it, ++i)
    {
      *it = {{ globalOffset + (localView.index(i))[0] }};
    }
    return it;
  }

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = typename TestspaceCoefficientMatrix::GridView;
  using size_type = std::size_t;


  using Node = OptimalTestBasisNode<TestspaceCoefficientMatrix, testIndex>;

#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,8)
  using IndexSet = OptimalTestBasisNodeIndexSet<TestspaceCoefficientMatrix,
                                                testIndex, MI>;
#endif

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
  static constexpr size_type maxMultiIndexSize = 1;
  static constexpr size_type minMultiIndexSize = 1;
  static constexpr size_type multiIndexBufferSize = 1;
#else
  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 1>;
#endif

  using TestSearchSpaces = typename TestspaceCoefficientMatrix::TestSpaces;
  using SolutionSpaces = typename TestspaceCoefficientMatrix::SolutionSpaces;


  /** \brief Constructor for a given test coefficient matrix */
  OptimalTestBasisPreBasis(TestspaceCoefficientMatrix& testCoeffMat) :
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
   * to in testspaceCoefficientMatrix_ to be updated separately.
   */
  void update (const GridView&)
  {}

  Node makeNode() const
  {
    return Node{testspaceCoefficientMatrix_};
  }

#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,8)
  IndexSet makeIndexSet() const
  {
    return IndexSet{*this};
  }
#endif

  size_type size() const
  {
    return Hybrid::accumulate(*testspaceCoefficientMatrix_.bilinearForm()
                              .getSolutionSpaces(), 0,
                              [&](size_type acc, const auto& s) {
                                return acc + s.size();
                              });
  }

  //! Return number possible values for next position in multi index
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
                                return acc + s.preBasis().maxNodeSize();
                              });
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  template<typename It>
  It indices(const Node& node, It it) const
  {
    using namespace boost::hana;
    return fold_left(
              make_range(int_c<0>,
                int_c<std::tuple_size<SolutionLocalViews>::value>),
              it,
              [&](It it, auto i) {
                return computeIndices(it,
                  std::get<i>(node.localViewsSolution()),
                  globalOffsets[i]);
              });
  }

//protected:
  // TODO: store testspaceCoefficientMatrix_ by reference or by value?
  TestspaceCoefficientMatrix& testspaceCoefficientMatrix_;

  std::array<size_t,std::tuple_size<SolutionSpaces>::value> globalOffsets;

};



template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
class OptimalTestBasisNode :
  public LeafBasisNode,
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
  static constexpr int dim = GV::dimension;

  using TestSearchFiniteElement
      = typename TestSearchSpace::LocalView::Tree::FiniteElement;

  using SolutionLocalViews = Dune::detail::getLocalViews_t<SolutionSpaces>;
  using TestLocalViews = Dune::detail::getLocalViews_t<TestSearchSpaces>;

  static const bool testSearchSpaceIsRefined
    = is_RefinedFiniteElement<TestSearchSpace>::value;

public:

  using size_type = std::size_t;
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

  OptimalTestBasisNode(TestspaceCoefficientMatrix& testCoeffMat) :
    testspaceCoefficientMatrix(testCoeffMat),
    finiteElement_(nullptr),
    testSearchSpace_(nullptr),
    localViewsSolution_(Dune::detail::getLocalViews(*testCoeffMat.bilinearForm()
                                                .getSolutionSpaces())),
    localViewsTest(Dune::detail::getLocalViews(*testCoeffMat.bilinearForm()
                                                    .getTestSpaces()))
  {}

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



#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,8)
template<typename TestspaceCoefficientMatrix, std::size_t testIndex, class MI>
class OptimalTestBasisNodeIndexSet
{
public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = OptimalTestBasisPreBasis<TestspaceCoefficientMatrix, testIndex, MI>;

  using Node = typename PreBasis::Node;

  OptimalTestBasisNodeIndexSet(const PreBasis& preBasis) :
    preBasis_(&preBasis)
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
using OptimalTestBasis = typename std::conditional
           < Dune::is_RefinedFiniteElement<
               typename std::tuple_element<testIndex,
                     typename TestspaceCoefficientMatrix::TestSpaces
               >::type>::value
           , RefinedGlobalBasis<OptimalTestBasisPreBasis<
                                     TestspaceCoefficientMatrix,
                                     testIndex, FlatMultiIndex<std::size_t>>>
           , DefaultGlobalBasis<OptimalTestBasisPreBasis<
                                     TestspaceCoefficientMatrix,
                                     testIndex, FlatMultiIndex<std::size_t>>>
           >::type;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_OPTIMALTESTBASIS_HH
