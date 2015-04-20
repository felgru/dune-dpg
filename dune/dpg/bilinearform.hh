// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_BILINEARFORM_HH
#define DUNE_DPG_BILINEARFORM_HH

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
   * class BilinearForm
   *
   * \tparam TestSpaces     tuple of test spaces
   * \tparam SolutionSpaces tuple of solution spaces
   * \tparam BilinearTerms  tuple of IntegralTerm
   */
  template<class TestSpaces, class SolutionSpaces, class BilinearTerms>
  class BilinearForm
  {
  public:
    typedef typename boost::fusion::result_of::as_vector<
        typename boost::fusion::result_of::
        transform<TestSpaces, getLocalView>::type>::type TestLocalView;
    typedef typename boost::fusion::result_of::as_vector<
        typename boost::fusion::result_of::
        transform<SolutionSpaces, getLocalView>::type>::type SolutionLocalView;

    BilinearForm () = delete;
    constexpr BilinearForm (TestSpaces     testSpaces,
                            SolutionSpaces solutionSpaces,
                            BilinearTerms  terms)
               : testSpaces(testSpaces),
                 solutionSpaces(solutionSpaces),
                 terms(terms),
                 testLocalView(nullptr),
                 solutionLocalView(nullptr)
    { };

    /** Compute the stiffness matrix for a single element
     *  \pre bind has to be called on the current element before. */
    template <class MatrixType>
    void getLocalMatrix(MatrixType& elementMatrix) const
    {
      // Set all matrix entries to zero
      elementMatrix.setSize(localTotalTestSize,
                            localTotalSolutionSize);
      elementMatrix = 0;      // fills the entire matrix with zeroes

      boost::fusion::for_each(terms,
              getLocalMatrixHelper
                      <MatrixType,
                       typename std::remove_pointer
                                <decltype(testLocalView)>::type,
                       typename std::remove_pointer
                                <decltype(solutionLocalView)>::type>
                      (*testLocalView,
                       *solutionLocalView,
                       elementMatrix,
                       localTestSpaceOffsets,
                       localSolutionSpaceOffsets));
    };

    void getOccupationPattern(MatrixIndexSet& nb) const;

    void bind(const TestLocalView& tlv, const SolutionLocalView& slv)
    {
      testLocalView     = const_cast<TestLocalView*>
                              (std::addressof(tlv));
      solutionLocalView = const_cast<SolutionLocalView*>
                              (std::addressof(slv));

      using namespace boost::fusion;
      using TestSize     = typename result_of::size<TestLocalView>::type;
      using SolutionSize = typename result_of::size<SolutionLocalView>::type;

      /* set up local offsets */
      fold(zip(localTestSpaceOffsets, tlv), 0, offsetHelper());
      localTotalTestSize =
          localTestSpaceOffsets[TestSize::value-1]
          + at_c<TestSize::value-1>(tlv)->size();

      fold(zip(localSolutionSpaceOffsets, slv),
           0, offsetHelper());
      localTotalSolutionSize =
          localSolutionSpaceOffsets[SolutionSize::value-1]
          + at_c<SolutionSize::value-1>(slv)->size();
    };

    const TestSpaces& getTestSpaces() const
    { return testSpaces; };

    const SolutionSpaces& getSolutionSpaces() const
    { return solutionSpaces; };

    const BilinearTerms& getTerms() const
    {
      return terms;
    };

    using TestSpaceIndexArray = size_t[std::tuple_size<TestSpaces>::value];
    using SolutionSpaceIndexArray
        = size_t[std::tuple_size<SolutionSpaces>::value];

    const TestSpaceIndexArray& getLocalTestSpaceOffsets() const
    { return localTestSpaceOffsets; };

    const SolutionSpaceIndexArray& getLocalSolutionSpaceOffsets() const
    { return localSolutionSpaceOffsets; };

  private:
    TestSpaces     testSpaces;
    SolutionSpaces solutionSpaces;
    BilinearTerms  terms;

    size_t localTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
    size_t localSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
    size_t localTotalTestSize,
           localTotalSolutionSize;

    TestLocalView*     testLocalView;
    SolutionLocalView* solutionLocalView;
  };

template<class TestSpaces, class SolutionSpaces, class BilinearTerms>
auto make_BilinearForm(TestSpaces     testSpaces,
                       SolutionSpaces solutionSpaces,
                       BilinearTerms  terms)
    -> BilinearForm<TestSpaces, SolutionSpaces, BilinearTerms>
{
  return BilinearForm<TestSpaces, SolutionSpaces, BilinearTerms>
                      (testSpaces,
                       solutionSpaces,
                       terms);
}

template<class TestSpaces, class SolutionSpaces, class BilinearTerms>
void BilinearForm<TestSpaces, SolutionSpaces, BilinearTerms>::
getOccupationPattern(MatrixIndexSet& nb) const
{
  using namespace boost::fusion;

  // Total number of degrees of freedom
  auto testBasisIndexSet = as_vector(transform(testSpaces, getIndexSet()));
  auto solutionBasisIndexSet = as_vector(transform(solutionSpaces,
                                                   getIndexSet()));

  /* set up global offsets */
  size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
  fold(zip(globalTestSpaceOffsets, testBasisIndexSet), 0, globalOffsetHelper());
  size_t globalTotalTestSize =
      globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value-1]
      + at_c<std::tuple_size<TestSpaces>::value-1>(testBasisIndexSet).size();

  fold(zip(globalSolutionSpaceOffsets, solutionBasisIndexSet),
       globalTotalTestSize, globalOffsetHelper());
  size_t globalTotalSolutionSize =
      globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value-1]
      + at_c<std::tuple_size<SolutionSpaces>::value-1>
            (solutionBasisIndexSet).size()
      - globalTotalTestSize;

  // A view on the FE basis on a single element
  auto solutionLocalView = as_vector(transform(solutionSpaces,
                                               localViewFromFEBasis()));
  auto testLocalView     = as_vector(transform(testSpaces,
                                               localViewFromFEBasis()));

  auto solutionLocalIndexSet = as_vector(transform(solutionBasisIndexSet,
                                                   getLocalIndexSet()));
  auto testLocalIndexSet     = as_vector(transform(testBasisIndexSet,
                                                   getLocalIndexSet()));

  // Get the grid view from the finite element basis
  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(testSpaces).gridView();

  /* create set of index pairs from bilinearTerms to loop over. */
  typedef typename boost::mpl::fold<
      typename boost::mpl::transform<
          /* This as_vector is probably not needed for boost::fusion 1.58
           * or higher. */
          typename result_of::as_vector<BilinearTerms>::type
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
    for_each(solutionLocalView, applyBind<decltype(e)>(e));
    for_each(testLocalView, applyBind<decltype(e)>(e));

    for_each(zip(solutionLocalIndexSet, solutionLocalView),
             make_fused_procedure(bindLocalIndexSet()));
    for_each(zip(testLocalIndexSet, testLocalView),
             make_fused_procedure(bindLocalIndexSet()));

    /* TODO: This is specific to the saddlepoint formulation. */
    auto gOPH = getOccupationPatternHelper<decltype(testLocalView),
                                           decltype(solutionLocalView),
                                           decltype(testLocalIndexSet),
                                           decltype(solutionLocalIndexSet),
                                           true>
                        (testLocalView,
                         solutionLocalView,
                         testLocalIndexSet,
                         solutionLocalIndexSet,
                         globalTestSpaceOffsets,
                         globalSolutionSpaceOffsets,
                         nb);
    for_each(indexPairs,
        std::ref(gOPH));
  }

  /* free memory handled by raw pointers */
  for_each(testLocalIndexSet,     default_deleter());
  for_each(solutionLocalIndexSet, default_deleter());
  for_each(testLocalView,         default_deleter());
  for_each(solutionLocalView,     default_deleter());
}

} // end namespace Dune

#endif // DUNE_DPG_BILINEARFORM_HH