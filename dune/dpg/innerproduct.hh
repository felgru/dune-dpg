// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_INNERPRODUCT_HH
#define DUNE_DPG_INNERPRODUCT_HH

#include <tuple>
#include <vector>
#include <functional>
#include <memory>
#include <type_traits>

#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/fusion/adapted/array.hpp>
#include <boost/fusion/container/vector/convert.hpp>
#include <boost/fusion/container/set/convert.hpp>
#include <boost/fusion/sequence/intrinsic/size.hpp>
#include <boost/fusion/algorithm/auxiliary/copy.hpp>
#include <boost/fusion/algorithm/transformation/transform.hpp>
#include <boost/fusion/algorithm/transformation/zip.hpp>
#include <boost/fusion/algorithm/iteration/accumulate.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/functional/generation/make_fused_procedure.hpp>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include "assemble_helper.hh"
#include "assemble_types.hh"
#include "integralterm.hh"

namespace Dune {

  /**
   * \brief This class describes an inner product.
   *
   * \tparam TSpaces            tuple of test spaces
   * \tparam InnerProductTerms  tuple of IntegralTerm
   */
  template<class TSpaces, class InnerProductTerms>
  class InnerProduct
  {
  public:
    //! tuple type of test spaces
    typedef TSpaces TestSpaces;
    //! tuple type for the local views of the test spaces
    typedef typename boost::fusion::result_of::as_vector<
        typename boost::fusion::result_of::
        transform<TestSpaces, detail::getLocalView>::type>::type TestLocalView;

    InnerProduct () = delete;
    /**
     * \brief constructor for InnerProduct
     *
     * \note For your convenience, use make_InnerProduct() instead.
     */
    constexpr InnerProduct (const TestSpaces&        testSpaces,
                            const InnerProductTerms& terms)
               : testSpaces(testSpaces),
                 terms(terms),
                 testLocalView(nullptr)
    { }

    /**
     * \brief Compute the stiffness matrix for a single element.
     *
     * \pre bind() has to be called on the current element before.
     *
     * \param[out] elementMatrix will hold the local system matrix
     */
    template <class MatrixType>
    void getLocalMatrix(MatrixType& elementMatrix) const
    {
      // Set all matrix entries to zero
      elementMatrix.setSize(localTotalTestSize,
                            localTotalTestSize);
      elementMatrix = 0;

      boost::fusion::for_each(terms,
              detail::getLocalMatrixHelper
                      <MatrixType,
                       TestSpaces,
                       TestSpaces>
                      (*testLocalView,
                       *testLocalView,
                       elementMatrix,
                       localTestSpaceOffsets,
                       localTestSpaceOffsets));
    }

    /**
     * \brief Creates the occupation pattern for the system matrix.
     *
     * This creates an index set of all those entries of the system
     * matrix that might be non-zero.
     *
     * \param nb  the matrix index set of the non-zero entries
     */
    void getOccupationPattern(MatrixIndexSet& nb) const;

    /**
     * \brief Bind the InnerProduct to the given localViews.
     *
     * \pre The given localViews have to be bound to the same element.
     */
    void bind(const TestLocalView& tlv)
    {
      testLocalView = std::addressof(tlv);

      using namespace boost::fusion;
      using namespace Dune::detail;

      /* set up local offsets */
      localTotalTestSize =
          fold(zip(localTestSpaceOffsets, tlv), (size_t)0, offsetHelper());
    }

    /**
     * \brief Does exactly what it says on the tin.
     */
    const TestSpaces& getTestSpaces() const
    { return testSpaces; }

    /**
     * \brief Does exactly what it says on the tin.
     */
    const InnerProductTerms& getTerms() const
    {
      return terms;
    }

  private:
    TestSpaces         testSpaces;
    InnerProductTerms  terms;
    size_t localTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
    size_t localTotalTestSize;
    const TestLocalView* testLocalView;
  };

/**
 * \brief Creates an InnerProduct,
 *        deducing the target type from the types of arguments.
 *
 * \param testSpaces  a tuple of test spaces
 * \param terms       a tuple of IntegralTerm
 */
template<class TestSpaces, class InnerProductTerms>
auto make_InnerProduct(TestSpaces        testSpaces,
                       InnerProductTerms terms)
    -> InnerProduct<TestSpaces, InnerProductTerms>
{
  return InnerProduct<TestSpaces, InnerProductTerms> (testSpaces, terms);
}

template<class TestSpaces, class InnerProductTerms>
void InnerProduct<TestSpaces, InnerProductTerms>::
getOccupationPattern(MatrixIndexSet& nb) const
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  /* set up global offsets */
  size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
  fold(zip(globalTestSpaceOffsets, testSpaces),
       (size_t)0, globalOffsetHelper());

  auto testLocalView = as_vector(transform(testSpaces,
                                           getLocalView()));
  auto testLocalIndexSet = as_vector(transform(testSpaces,
                                               getLocalIndexSet()));

  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(testSpaces).gridView();

  /* create set of index pairs from innerProductTerms to loop over. */
  typedef typename boost::mpl::fold<
      typename boost::mpl::transform<
          /* This as_vector is probably not needed for boost::fusion 1.58
           * or higher. */
          typename result_of::as_vector<InnerProductTerms>::type
        , mpl::firstTwo<boost::mpl::_1>
        >::type
    , boost::mpl::set0<>
    , boost::mpl::insert<boost::mpl::_1,boost::mpl::_2>
    >::type IndexPairs;
  auto indexPairs = IndexPairs{};

  for(const auto& e : elements(gridView))
  {
    for_each(testLocalView, applyBind<decltype(e)>(e));

    for_each(zip(testLocalIndexSet, testLocalView),
             make_fused_procedure(bindLocalIndexSet()));

    auto gOPH = getOccupationPatternHelper<decltype(testLocalView),
                                           decltype(testLocalView),
                                           decltype(testLocalIndexSet),
                                           decltype(testLocalIndexSet),
                                           false>
                        (testLocalView,
                         testLocalView,
                         testLocalIndexSet,
                         testLocalIndexSet,
                         globalTestSpaceOffsets,
                         globalTestSpaceOffsets,
                         nb);
    for_each(indexPairs,
        std::ref(gOPH));
  }
}

} // end namespace Dune

#endif // DUNE_DPG_INNERPRODUCT_HH
