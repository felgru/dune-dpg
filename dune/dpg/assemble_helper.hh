// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_ASSEMBLE_HELPER_HH
#define DUNE_DPG_ASSEMBLE_HELPER_HH

#include <array>
#include <memory>
#include <utility>
#include <tuple>
#include <dune/common/hybridutilities.hh>
#include <dune/common/std/memory.hh>
#include <dune/common/tupleutility.hh>
#include <dune/dpg/assemble_types.hh>
#include <dune/dpg/functions/localindexsetiteration.hh>
#include <dune/dpg/spacetuple.hh>
#include <dune/istl/matrixindexset.hh>
#include <boost/hana.hpp>

namespace Dune {

namespace detail {

struct getLocalViewFunctor
{
  template<class T>
  struct TypeEvaluator
  {
    typedef typename T::LocalView Type;
  };

  template<class T>
  typename TypeEvaluator<T>::Type operator()(const T& t) const
  {
    return t.localView();
  }
};

template<class Spaces>
using getLocalViews_t
    = typename ForEachType<detail::getLocalViewFunctor::TypeEvaluator,
                           Spaces>::Type;

template<class Spaces>
inline getLocalViews_t<Spaces> getLocalViews(const Spaces& spaces) {
  return genericTransformSpaceTuple(spaces, getLocalViewFunctor());
}

template<class LocalViews, class Element>
inline void bindLocalViews(LocalViews& localViews,
                           const Element& e)
{
  Hybrid::forEach(localViews, [&](auto& lv) { lv.bind(e); });
}

template<class Spaces, class GridView>
inline void updateSpaces(Spaces& spaces,
                         const GridView& gridView)
{
  Hybrid::forEach(spaces, [&](auto& s) { s.update(gridView); });
}

template <class MatrixType,
          class LhsSpaces,
          class RhsSpaces>
struct getLocalMatrixHelper
{
  template<class T, class Tuple>
  using array_of_same_size = std::array<T,std::tuple_size<Tuple>::value>;

  //! tuple type for the local views of the test spaces
  using LhsLocalViews = getLocalViews_t<LhsSpaces>;
  //! tuple type for the local views of the solution spaces
  using RhsLocalViews = getLocalViews_t<RhsSpaces>;

  getLocalMatrixHelper(LhsLocalViews& lhsLocalViews,
                       RhsLocalViews& rhsLocalViews,
                       MatrixType& elementMatrix,
                       const array_of_same_size<size_t, LhsLocalViews>&
                           lhsLocalSpaceOffsets,
                       const array_of_same_size<size_t, RhsLocalViews>&
                           rhsLocalSpaceOffsets)
      : lhsLocalViews(lhsLocalViews),
        rhsLocalViews(rhsLocalViews),
        elementMatrix(elementMatrix),
        lhsLocalSpaceOffsets(lhsLocalSpaceOffsets),
        rhsLocalSpaceOffsets(rhsLocalSpaceOffsets)
  {}

  /**
   * \tparam Term an IntegralTerm
   */
  template <class lhsSpaceIndex,
            class rhsSpaceIndex,
            class BilinearTerm>
  void operator()
         (const BilinearTermWithIndices<
          lhsSpaceIndex,
          rhsSpaceIndex,
          BilinearTerm>& termWithIndices) const
  {
    auto& lhsLV =
        std::get<lhsSpaceIndex::value>(lhsLocalViews);
    auto& rhsLV =
        std::get<rhsSpaceIndex::value>(rhsLocalViews);
    const size_t lhsLocalSpaceOffset =
        lhsLocalSpaceOffsets[lhsSpaceIndex::value];
    const size_t rhsLocalSpaceOffset =
        rhsLocalSpaceOffsets[rhsSpaceIndex::value];

    termWithIndices.term
                   .getLocalMatrix(lhsLV,
                                   rhsLV,
                                   elementMatrix,
                                   lhsLocalSpaceOffset,
                                   rhsLocalSpaceOffset);
  }

private:
  LhsLocalViews& lhsLocalViews;
  RhsLocalViews& rhsLocalViews;
  MatrixType& elementMatrix;
  const array_of_same_size<size_t, LhsLocalViews>&
      lhsLocalSpaceOffsets;
  const array_of_same_size<size_t, RhsLocalViews>&
      rhsLocalSpaceOffsets;
};

template <class VectorType,
          class TestSpaces>
struct getLocalVectorHelper
{
  template<class T, class Tuple>
  using array_of_same_size = std::array<T,std::tuple_size<Tuple>::value>;

  //! tuple type for the local views of the test spaces
  using TestLocalViews = getLocalViews_t<TestSpaces>;

  getLocalVectorHelper(TestLocalViews& testLocalViews,
                       VectorType& elementVector,
                       const array_of_same_size<size_t, TestLocalViews>&
                           localTestSpaceOffsets)
      : testLocalViews(testLocalViews),
        elementVector(elementVector),
        localTestSpaceOffsets(localTestSpaceOffsets)
  {}

  /**
   * \tparam Term a LinearIntegralTerm
   */
  template <class testSpaceIndex,
            class Term>
  void operator()
         (const LinearTermWithIndex<testSpaceIndex, Term>& termWithIndex) const
  {
    auto& testLV = std::get<testSpaceIndex::value>(testLocalViews);
    const size_t localTestSpaceOffset =
        localTestSpaceOffsets[testSpaceIndex::value];

    termWithIndex.term
                 .getLocalVector(testLV,
                                 elementVector,
                                 localTestSpaceOffset);
  }

private:
  TestLocalViews& testLocalViews;
  VectorType& elementVector;
  const array_of_same_size<size_t, TestLocalViews>&
      localTestSpaceOffsets;
};

template <class TestLocalViews,
          class SolutionLocalViews,
          bool mirror = false>
struct getOccupationPatternHelper
{
  template<class T, class Tuple>
  using array_of_same_size = std::array<T,std::tuple_size<Tuple>::value>;

  getOccupationPatternHelper(
      const TestLocalViews& testLocalViews,
      const SolutionLocalViews& solutionLocalViews,
      const array_of_same_size<size_t, TestLocalViews>&
          globalTestSpaceOffsets,
      const array_of_same_size<size_t, SolutionLocalViews>&
          globalSolutionSpaceOffsets,
      Dune::MatrixIndexSet& nb)
      : testLocalViews(testLocalViews),
        solutionLocalViews(solutionLocalViews),
        globalTestSpaceOffsets(globalTestSpaceOffsets),
        globalSolutionSpaceOffsets(globalSolutionSpaceOffsets),
        nb(nb)
  {}

  template <class testSpaceIndex,
            class solutionSpaceIndex>
  void operator()
         (const boost::hana::tuple<testSpaceIndex,
                                   solutionSpaceIndex>&)
  {
    const auto& testLIS =
        std::get<testSpaceIndex::value>(testLocalViews);
    const auto& solutionLIS =
        std::get<solutionSpaceIndex::value>(solutionLocalViews);
    const size_t globalTestSpaceOffset =
        globalTestSpaceOffsets[testSpaceIndex::value];
    const size_t globalSolutionSpaceOffset =
        globalSolutionSpaceOffsets[solutionSpaceIndex::value];
    using TestMultiIndex
        = typename std::decay_t<decltype(testLIS)>::MultiIndex;
    using SolutionMultiIndex
        = typename std::decay_t<decltype(solutionLIS)>::MultiIndex;

    auto fillOccupationPattern = [&](size_t /* i */, TestMultiIndex gi)
        {
          auto fillOccupationPatternInner
            = [&](size_t /* j */, SolutionMultiIndex gj)
              {
                // Add a nonzero entry to the matrix
                nb.add(gi[0]+globalTestSpaceOffset,
                       gj[0]+globalSolutionSpaceOffset);
                if(mirror) {
                    nb.add(gj[0]+globalSolutionSpaceOffset,
                           gi[0]+globalTestSpaceOffset);
                }
              };
          iterateOverLocalIndices(
              solutionLIS,
              fillOccupationPatternInner,
              [](size_t /* j */){},
              [&](size_t j, SolutionMultiIndex gj, double /* wj */)
              {
                fillOccupationPatternInner(j, gj);
              }
          );
        };
    iterateOverLocalIndices(
        testLIS,
        fillOccupationPattern,
        [](size_t /* i */){},
        [&](size_t i, TestMultiIndex gi, double /* wi */)
        {
          fillOccupationPattern(i, gi);
        }
    );
  }

private:
  const TestLocalViews& testLocalViews;
  const SolutionLocalViews& solutionLocalViews;
  const array_of_same_size<size_t, TestLocalViews>&
      globalTestSpaceOffsets;
  const array_of_same_size<size_t, SolutionLocalViews>&
      globalSolutionSpaceOffsets;
  Dune::MatrixIndexSet& nb;
};

template<class Indices,
         bool mirror,
         class LeftLocalViews, class RightLocalViews,
         class LeftOffsets, class RightOffsets>
inline void getOccupationPattern(
      const LeftLocalViews& leftLocalViews,
      const RightLocalViews& rightLocalViews,
      const LeftOffsets& leftGlobalOffsets,
      const RightOffsets& rightGlobalOffsets,
      MatrixIndexSet& occupationPattern) {
  auto gOPH = getOccupationPatternHelper<LeftLocalViews,
                                         RightLocalViews,
                                         mirror>
                      (leftLocalViews,
                       rightLocalViews,
                       leftGlobalOffsets,
                       rightGlobalOffsets,
                       occupationPattern);
  boost::hana::for_each(Indices{}, [&](const auto& index) { gOPH(index); });
}

template<class LocalMatrix, class GlobalMatrix,
         class TestLocalViews, class SolutionLocalViews,
         class TestOffsets, class SolutionOffsets,
         bool mirror = false>
struct localToGlobalCopier
{

  localToGlobalCopier
         (const LocalMatrix& lm, GlobalMatrix& gm,
          TestLocalViews const & testLocalViews,
          TestOffsets const & testLocalOffsets,
          TestOffsets const & testGlobalOffsets,
          SolutionLocalViews const & solutionLocalViews,
          SolutionOffsets const & solutionLocalOffsets,
          SolutionOffsets const & solutionGlobalOffsets)
      : elementMatrix(lm), matrix(gm),
        testLocalViews(testLocalViews),
        testLocalOffsets(testLocalOffsets),
        testGlobalOffsets(testGlobalOffsets),
        solutionLocalViews(solutionLocalViews),
        solutionLocalOffsets(solutionLocalOffsets),
        solutionGlobalOffsets(solutionGlobalOffsets) {}

  template <class testSpaceIndex,
            class solutionSpaceIndex>
  void operator()
         (const boost::hana::tuple<testSpaceIndex,
                                   solutionSpaceIndex>&)
  {
    const auto& testLocalView
        = std::get<testSpaceIndex::value>(testLocalViews);
    const size_t testLocalOffset = testLocalOffsets[testSpaceIndex::value];
    const size_t testGlobalOffset = testGlobalOffsets[testSpaceIndex::value];
    const auto& solutionLocalView
        = std::get<solutionSpaceIndex::value>(solutionLocalViews);
    const size_t solutionLocalOffset
        = solutionLocalOffsets[solutionSpaceIndex::value];
    const size_t solutionGlobalOffset
        = solutionGlobalOffsets[solutionSpaceIndex::value];

    using TestMultiIndex
        = typename std::decay_t<decltype(testLocalView)>::MultiIndex;
    using SolutionMultiIndex
        = typename std::decay_t<decltype(solutionLocalView)>::MultiIndex;

    addToGlobalMatrix(
        testLocalView,
        solutionLocalView,
        [&](size_t i, size_t j) -> double
        {
          return elementMatrix[i+testLocalOffset][j+solutionLocalOffset];
        },
        [&](TestMultiIndex gi, SolutionMultiIndex gj)
        -> typename GlobalMatrix::block_type&
        {
          const auto row = gi[0] + testGlobalOffset;
          const auto col = gj[0] + solutionGlobalOffset;
          return matrix[row][col];
        }
    );
    if(mirror) {
      addToGlobalMatrix(
          testLocalView,
          solutionLocalView,
          [&](size_t i, size_t j) -> double
          {
            return elementMatrix[i+testLocalOffset][j+solutionLocalOffset];
          },
          [&](TestMultiIndex gi, SolutionMultiIndex gj)
          -> typename GlobalMatrix::block_type&
          {
            const auto row = gi[0] + testGlobalOffset;
            const auto col = gj[0] + solutionGlobalOffset;
            return matrix[col][row];
          }
      );
    }
  }

private:
  const LocalMatrix& elementMatrix;
  GlobalMatrix& matrix;

  const TestLocalViews & testLocalViews;
  const TestOffsets & testLocalOffsets;
  const TestOffsets & testGlobalOffsets;
  const SolutionLocalViews & solutionLocalViews;
  const SolutionOffsets & solutionLocalOffsets;
  const SolutionOffsets & solutionGlobalOffsets;
};

template<class Indices,
         class LocalMatrix, class GlobalMatrix,
         class TestLocalViews, class SolutionLocalViews,
         class TestOffsets, class SolutionOffsets,
         bool mirror = false>
inline void copyLocalToGlobalMatrix(
    const LocalMatrix&        elementMatrix,
    GlobalMatrix&             matrix,
    const TestLocalViews&     testLocalViews,
    const TestOffsets&        localTestSpaceOffsets,
    const TestOffsets&        globalTestSpaceOffsets,
    const SolutionLocalViews& solutionLocalViews,
    const SolutionOffsets&    localSolutionSpaceOffsets,
    const SolutionOffsets&    globalSolutionSpaceOffsets) {

  auto cpMatrix = localToGlobalCopier<LocalMatrix, GlobalMatrix,
                 TestLocalViews, SolutionLocalViews,
                 TestOffsets, SolutionOffsets,
                 mirror> (elementMatrix, matrix,
                          testLocalViews,
                          localTestSpaceOffsets,
                          globalTestSpaceOffsets,
                          solutionLocalViews,
                          localSolutionSpaceOffsets,
                          globalSolutionSpaceOffsets);

  /* copy every local submatrix indexed by a pair of indices from
   * Indices exactly once. */
  boost::hana::for_each(Indices{}, cpMatrix);
}

template<class Indices,
         class LocalMatrix, class GlobalMatrix,
         class TestLocalViews, class SolutionLocalViews,
         class TestOffsets, class SolutionOffsets>
inline void copyLocalToGlobalMatrixSymmetric(
    const LocalMatrix&        elementMatrix,
    GlobalMatrix&             matrix,
    const TestLocalViews&     testLocalViews,
    const TestOffsets&        localTestSpaceOffsets,
    const TestOffsets&        globalTestSpaceOffsets,
    const SolutionLocalViews& solutionLocalViews,
    const SolutionOffsets&    localSolutionSpaceOffsets,
    const SolutionOffsets&    globalSolutionSpaceOffsets) {
  copyLocalToGlobalMatrix<Indices, LocalMatrix, GlobalMatrix,
                          TestLocalViews, SolutionLocalViews,
                          TestOffsets, SolutionOffsets,
                          true>
                         ( elementMatrix,
                           matrix,
                           testLocalViews,
                           localTestSpaceOffsets,
                           globalTestSpaceOffsets,
                           solutionLocalViews,
                           localSolutionSpaceOffsets,
                           globalSolutionSpaceOffsets);
}

template<class LocalVector, class GlobalVector,
         class TestLocalViews, class TestOffsets>
struct localToGlobalRHSCopier
{

  localToGlobalRHSCopier
         (const LocalVector& lv, GlobalVector& gv,
          TestLocalViews const & testLocalViews,
          TestOffsets const & testLocalOffsets,
          TestOffsets const & testGlobalOffsets)
      : localRhs(lv), rhs(gv),
        testLocalViews(testLocalViews),
        testLocalOffsets(testLocalOffsets),
        testGlobalOffsets(testGlobalOffsets) {}

  template <class TestSpaceIndex>
  void operator() (const TestSpaceIndex&) const
  {
    constexpr auto testSpaceIndex = TestSpaceIndex::type::value;
    const auto& testLocalView
        = std::get<testSpaceIndex>(testLocalViews);
    const size_t testLocalOffset = testLocalOffsets[testSpaceIndex];
    const size_t testGlobalOffset = testGlobalOffsets[testSpaceIndex];

    using MultiIndex
        = typename std::decay_t<decltype(testLocalView)>::MultiIndex;

    iterateOverLocalIndices(
        testLocalView,
        [&](size_t i, MultiIndex gi)
        {
          rhs[gi[0] + testGlobalOffset]
              += localRhs[i+testLocalOffset];
        },
        [](size_t /* i */){},
        [&](size_t i, MultiIndex gi, double wi)
        {
          rhs[gi[0] + testGlobalOffset]
              += wi * localRhs[i+testLocalOffset];
        }
    );
  }

private:
  const LocalVector& localRhs;
  GlobalVector& rhs;

  const TestLocalViews & testLocalViews;
  const TestOffsets & testLocalOffsets;
  const TestOffsets & testGlobalOffsets;
};

template<class Indices,
         class LocalVector, class GlobalVector,
         class TestLocalViews,
         class TestOffsets>
inline void copyLocalToGlobalVector(
    const LocalVector&        elementVector,
    GlobalVector&             vector,
    const TestLocalViews&     testLocalViews,
    const TestOffsets&        localTestSpaceOffsets,
    const TestOffsets&        globalTestSpaceOffsets) {

  auto cpRhs = localToGlobalRHSCopier<LocalVector, GlobalVector,
                    TestLocalViews, TestOffsets>
                         (elementVector, vector,
                          testLocalViews,
                          localTestSpaceOffsets,
                          globalTestSpaceOffsets);

  /* copy every local subvector indexed by an index from
   * Indices exactly once. */
  boost::hana::for_each(Indices{}, cpRhs);
}

template<class SpacesOrLocalViews, class Offsets>
inline size_t computeOffsets(Offsets& offsets, const SpacesOrLocalViews& s,
                             size_t start = 0)
{
  size_t numDofs = start;

  Hybrid::forEach(
      Std::make_index_sequence<
          std::tuple_size<SpacesOrLocalViews>::value>{},
      [&](auto i) {
        offsets[i] = numDofs;
        numDofs += std::get<i>(s).size();
      });

  return numDofs;
}

template<size_t spaceIndex, class SpacesOrLocalViews>
inline size_t computeOffset(const SpacesOrLocalViews& s,
                            size_t start = 0)
{
  size_t numDofs = start;

  Hybrid::forEach(
      Std::make_index_sequence<spaceIndex>{},
      [&](auto i) {
        numDofs += std::get<i>(s).size();
      });

  return numDofs;
}

} } // end namespace Dune::detail

#endif // DUNE_DPG_ASSEMBLE_HELPER_HH
