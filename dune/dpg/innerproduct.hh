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

#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>

#include <dune/common/std/final.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discretescalarglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include "assemble_helper.hh"
#include "assemble_types.hh"
#include "integralterm.hh"

namespace Dune {

  /**
   * class InnerProduct
   *
   * \tparam TSpaces         tuple of test spaces
   * \tparam InnerProductTerms  tuple of IntegralTerm
   */
  template<class TSpaces, class InnerProductTerms>
  class InnerProduct
  {
  public:
    typedef TSpaces TestSpaces;
    typedef typename boost::fusion::result_of::as_vector<
        typename boost::fusion::result_of::
        transform<TestSpaces, getLocalView>::type>::type TestLocalView;

    InnerProduct () = delete;
    constexpr InnerProduct (TestSpaces        testSpaces,
                            InnerProductTerms terms)
               : testSpaces(testSpaces),
                 terms(terms),
                 testLocalView(nullptr)
    { };

    /** Compute the stiffness matrix for a single element
     *  \pre bind has to be called on the current element before. */
    template <class MatrixType>
    void getLocalMatrix(MatrixType& elementMatrix) const
    {
      // Set all matrix entries to zero
      elementMatrix.setSize(localTotalTestSize,
                            localTotalTestSize);
      elementMatrix = 0;      // fills the entire matrix with zeroes

      boost::fusion::for_each(terms,
              getLocalMatrixHelper
                      <MatrixType,
                       typename std::remove_pointer
                                <decltype(testLocalView)>::type,
                       typename std::remove_pointer
                                <decltype(testLocalView)>::type>
                      (*testLocalView,
                       *testLocalView,
                       elementMatrix,
                       localTestSpaceOffsets,
                       localTestSpaceOffsets));
    };

    void getOccupationPattern(MatrixIndexSet& nb) const;

    void bind(const TestLocalView& tlv)
    {
      testLocalView = const_cast<TestLocalView*>(std::addressof(tlv));

      using namespace boost::fusion;
      using TestSize = typename result_of::size<TestLocalView>::type;

      /* set up local offsets */
      fold(zip(localTestSpaceOffsets, tlv), 0, offsetHelper());
      localTotalTestSize =
          localTestSpaceOffsets[TestSize::value-1]
          + at_c<TestSize::value-1>(tlv)->size();
    };

    const TestSpaces& getTestSpaces() const
    { return testSpaces; };

    const InnerProductTerms& getTerms() const
    {
      return terms;
    };

  private:
    TestSpaces         testSpaces;
    InnerProductTerms  terms;
    size_t localTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
    size_t localTotalTestSize;
    TestLocalView*     testLocalView;
  };

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

  // Total number of degrees of freedom
  auto testBasisIndexSet = as_vector(transform(testSpaces, getIndexSet()));

  /* set up global offsets */
  size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
  fold(zip(globalTestSpaceOffsets, testBasisIndexSet), 0, globalOffsetHelper());
  size_t globalTotalTestSize =
      globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value-1]
      + at_c<std::tuple_size<TestSpaces>::value-1>(testBasisIndexSet).size();

  // A view on the FE basis on a single element
  auto testLocalView = as_vector(transform(testSpaces,
                                           localViewFromFEBasis()));

  auto testLocalIndexSet = as_vector(transform(testBasisIndexSet,
                                               getLocalIndexSet()));

  // Get the grid view from the finite element basis
  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(testSpaces).gridView();

  /* create set of index pairs from innerProductTerms to loop over. */
  typedef typename boost::mpl::fold<
      typename boost::mpl::transform<
          /* This as_vector is probably not needed for boost::fusion 1.58
           * or higher. */
          typename result_of::as_vector<InnerProductTerms>::type
        , firstTwo::mpl<boost::mpl::_1>
        >::type
    , boost::mpl::set0<>
    , boost::mpl::insert<boost::mpl::_1,boost::mpl::_2>
    >::type IndexPairs;
  auto indexPairs = IndexPairs{};

  // Loop over all leaf elements
  for(const auto& e : elements(gridView))
  {
    // Bind the local FE basis view to the current element
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

  /* free memory handled by raw pointers */
  for_each(testLocalIndexSet,     default_deleter());
  for_each(testLocalView,         default_deleter());
}

} // end namespace Dune

#endif // DUNE_DPG_INNERPRODUCT_HH
