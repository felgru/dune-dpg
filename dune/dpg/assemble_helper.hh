// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_ASSEMBLE_HELPER_HH
#define DUNE_DPG_ASSEMBLE_HELPER_HH

#include <functional>
#include <utility>
#include <tuple>
#include <dune/common/hybridutilities.hh>
#include <dune/common/std/memory.hh>
#include <dune/common/tupleutility.hh>
#include <dune/dpg/functions/concepts.hh>
#include <dune/dpg/functions/localindexsetiteration.hh>
#include <dune/istl/matrixindexset.hh>
#include <boost/hana.hpp>

namespace Dune {

namespace detail {

struct getLocalIndexSetFunctor
{
  template<class T>
  struct TypeEvaluator
  {
    typedef typename T::LocalIndexSet Type;
  };

  template<class T>
  typename TypeEvaluator<T>::Type operator()(const T& t) const
  {
    return t.localIndexSet();
  }
};

template<class Spaces>
using getLocalIndexSets_t
    = typename ForEachType<detail::getLocalIndexSetFunctor::TypeEvaluator,
                           Spaces>::Type;

template<class Spaces>
inline getLocalIndexSets_t<Spaces> getLocalIndexSets(const Spaces& spaces) {
  return genericTransformTuple(spaces, getLocalIndexSetFunctor());
}

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
  return genericTransformTuple(spaces, getLocalViewFunctor());
}

template<class LocalViews, class Element>
inline void bindLocalViews(LocalViews& localViews,
                           const Element& e)
{
  Hybrid::forEach(localViews, [&](auto& lv) { lv.bind(e); });
}

template<class LocalIndexSets, class LocalViews>
inline void bindLocalIndexSets(LocalIndexSets&   lis,
                               const LocalViews& lvs)
{
  Hybrid::forEach(
      Std::make_index_sequence<
          std::tuple_size<LocalIndexSets>::value>{},
      [&](auto i) {
        std::get<i>(lis).bind(std::get<i>(lvs));
      });
}

template <class MatrixType,
          class TestSpaces,
          class SolutionSpaces>
struct getLocalMatrixHelper
{
  template<class T, class Tuple>
  using array_of_same_size =
      T[std::tuple_size<Tuple>::value];

  //! tuple type for the local views of the test spaces
  using TestLocalViews = getLocalViews_t<TestSpaces>;
  //! tuple type for the local views of the solution spaces
  using SolutionLocalViews = getLocalViews_t<SolutionSpaces>;

  getLocalMatrixHelper(const TestLocalViews& testLocalViews,
                       const SolutionLocalViews& solutionLocalViews,
                       MatrixType& elementMatrix,
                       const array_of_same_size<size_t, TestLocalViews>&
                           localTestSpaceOffsets,
                       const array_of_same_size<size_t, SolutionLocalViews>&
                           localSolutionSpaceOffsets)
      : testLocalViews(testLocalViews),
        solutionLocalViews(solutionLocalViews),
        elementMatrix(elementMatrix),
        localTestSpaceOffsets(localTestSpaceOffsets),
        localSolutionSpaceOffsets(localSolutionSpaceOffsets)
  {}

  /**
   * \tparam Term an IntegralTerm
   */
  template <class solutionSpaceIndex,
            class testSpaceIndex,
            class Term>
  void operator()
         (const std::tuple<
          testSpaceIndex,
          solutionSpaceIndex,
          Term>& termTuple) const
  {
    const auto& term = std::get<2>(termTuple);

    const auto& testLV =
        std::get<testSpaceIndex::value>(testLocalViews);
    const auto& solutionLV =
        std::get<solutionSpaceIndex::value>(solutionLocalViews);
    const size_t localTestSpaceOffset =
        localTestSpaceOffsets[testSpaceIndex::value];
    const size_t localSolutionSpaceOffset =
        localSolutionSpaceOffsets[solutionSpaceIndex::value];

    term.getLocalMatrix(testLV,
                        solutionLV,
                        elementMatrix,
                        localTestSpaceOffset,
                        localSolutionSpaceOffset);
  }

private:
  const TestLocalViews& testLocalViews;
  const SolutionLocalViews& solutionLocalViews;
  MatrixType& elementMatrix;
  const array_of_same_size<size_t, TestLocalViews>&
      localTestSpaceOffsets;
  const array_of_same_size<size_t, SolutionLocalViews>&
      localSolutionSpaceOffsets;
};

template <class VectorType,
          class TestSpaces>
struct getLocalVectorHelper
{
  template<class T, class Tuple>
  using array_of_same_size =
      T[std::tuple_size<Tuple>::value];

  //! tuple type for the local views of the test spaces
  using TestLocalViews = getLocalViews_t<TestSpaces>;

  getLocalVectorHelper(const TestLocalViews& testLocalViews,
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
         (const std::tuple<testSpaceIndex, Term>& termTuple) const
  {
    const auto& term = std::get<1>(termTuple);

    const auto& testLV = std::get<testSpaceIndex::value>(testLocalViews);
    const size_t localTestSpaceOffset =
        localTestSpaceOffsets[testSpaceIndex::value];

    term.getLocalVector(testLV,
                        elementVector,
                        localTestSpaceOffset);
  }

private:
  const TestLocalViews& testLocalViews;
  VectorType& elementVector;
  const array_of_same_size<size_t, TestLocalViews>&
      localTestSpaceOffsets;
};

template <class TestLocalIndexSets,
          class SolutionLocalIndexSets,
          bool mirror = false>
struct getOccupationPatternHelper
{
  template<class T, class Tuple>
  using array_of_same_size = T[std::tuple_size<Tuple>::value];

  getOccupationPatternHelper(
      const TestLocalIndexSets& testLocalIndexSets,
      const SolutionLocalIndexSets& solutionLocalIndexSets,
      const array_of_same_size<size_t, TestLocalIndexSets>&
          globalTestSpaceOffsets,
      const array_of_same_size<size_t, SolutionLocalIndexSets>&
          globalSolutionSpaceOffsets,
      Dune::MatrixIndexSet& nb)
      : testLocalIndexSets(testLocalIndexSets),
        solutionLocalIndexSets(solutionLocalIndexSets),
        globalTestSpaceOffsets(globalTestSpaceOffsets),
        globalSolutionSpaceOffsets(globalSolutionSpaceOffsets),
        nb(nb)
  {}

  template <class testSpaceIndex,
            class solutionSpaceIndex>
  void operator()
         (const boost::hana::tuple<testSpaceIndex,
                                   solutionSpaceIndex>& indexTuple)
  {
    const auto& testLIS =
        std::get<testSpaceIndex::value>(testLocalIndexSets);
    const auto& solutionLIS =
        std::get<solutionSpaceIndex::value>(solutionLocalIndexSets);
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
          iterateOverLocalIndexSet(
              solutionLIS,
              fillOccupationPatternInner,
              [](size_t j){},
              [&](size_t j, SolutionMultiIndex gj, double /* wj */)
              {
                fillOccupationPatternInner(j, gj);
              }
          );
        };
    iterateOverLocalIndexSet(
        testLIS,
        fillOccupationPattern,
        [](size_t i){},
        [&](size_t i, TestMultiIndex gi, double /* wi */)
        {
          fillOccupationPattern(i, gi);
        }
    );
  }

private:
  const TestLocalIndexSets& testLocalIndexSets;
  const SolutionLocalIndexSets& solutionLocalIndexSets;
  const array_of_same_size<size_t, TestLocalIndexSets>&
      globalTestSpaceOffsets;
  const array_of_same_size<size_t, SolutionLocalIndexSets>&
      globalSolutionSpaceOffsets;
  Dune::MatrixIndexSet& nb;
};

template<class Indices,
         bool mirror,
         class LeftLocalIndexSets, class RightLocalIndexSets,
         class LeftOffsets, class RightOffsets>
inline void getOccupationPattern(
      const LeftLocalIndexSets& leftLocalIndexSets,
      const RightLocalIndexSets& rightLocalIndexSets,
      const LeftOffsets& leftGlobalOffsets,
      const RightOffsets& rightGlobalOffsets,
      MatrixIndexSet& occupationPattern) {
  auto gOPH = getOccupationPatternHelper<LeftLocalIndexSets,
                                         RightLocalIndexSets,
                                         mirror>
                      (leftLocalIndexSets,
                       rightLocalIndexSets,
                       leftGlobalOffsets,
                       rightGlobalOffsets,
                       occupationPattern);
  boost::hana::for_each(Indices{}, [&](const auto& index) { gOPH(index); });
}

template<class LocalMatrix, class GlobalMatrix,
         class TestLocalIndexSets, class SolutionLocalIndexSets,
         class TestOffsets, class SolutionOffsets,
         bool mirror = false>
struct localToGlobalCopier
{

  localToGlobalCopier
         (const LocalMatrix& lm, GlobalMatrix& gm,
          TestLocalIndexSets const & testLocalIndexSets,
          TestOffsets const & testLocalOffsets,
          TestOffsets const & testGlobalOffsets,
          SolutionLocalIndexSets const & solutionLocalIndexSets,
          SolutionOffsets const & solutionLocalOffsets,
          SolutionOffsets const & solutionGlobalOffsets)
      : elementMatrix(lm), matrix(gm),
        testLocalIndexSets(testLocalIndexSets),
        testLocalOffsets(testLocalOffsets),
        testGlobalOffsets(testGlobalOffsets),
        solutionLocalIndexSets(solutionLocalIndexSets),
        solutionLocalOffsets(solutionLocalOffsets),
        solutionGlobalOffsets(solutionGlobalOffsets) {}

  template <class testSpaceIndex,
            class solutionSpaceIndex>
  void operator()
         (const boost::hana::tuple<testSpaceIndex,
                                   solutionSpaceIndex>& indexTuple)
  {
    const auto& testLocalIndexSet
        = std::get<testSpaceIndex::value>(testLocalIndexSets);
    const size_t testLocalOffset = testLocalOffsets[testSpaceIndex::value];
    const size_t testGlobalOffset = testGlobalOffsets[testSpaceIndex::value];
    const auto& solutionLocalIndexSet
        = std::get<solutionSpaceIndex::value>(solutionLocalIndexSets);
    const size_t solutionLocalOffset
        = solutionLocalOffsets[solutionSpaceIndex::value];
    const size_t solutionGlobalOffset
        = solutionGlobalOffsets[solutionSpaceIndex::value];

    using TestMultiIndex
        = typename std::decay_t<decltype(testLocalIndexSet)>::MultiIndex;
    using SolutionMultiIndex
        = typename std::decay_t<decltype(solutionLocalIndexSet)>::MultiIndex;

    addToGlobalMatrix(
        testLocalIndexSet,
        solutionLocalIndexSet,
        [&](size_t i, size_t j) -> double
        {
          return elementMatrix[i+testLocalOffset][j+solutionLocalOffset];
        },
        [&](TestMultiIndex gi, SolutionMultiIndex gj)
        -> FieldMatrix<double, 1, 1>&
        {
          const auto row = gi[0] + testGlobalOffset;
          const auto col = gj[0] + solutionGlobalOffset;
          return matrix[row][col];
        }
    );
    if(mirror) {
      addToGlobalMatrix(
          testLocalIndexSet,
          solutionLocalIndexSet,
          [&](size_t i, size_t j) -> double
          {
            return elementMatrix[i+testLocalOffset][j+solutionLocalOffset];
          },
          [&](TestMultiIndex gi, SolutionMultiIndex gj)
          -> FieldMatrix<double, 1, 1>&
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

  const TestLocalIndexSets & testLocalIndexSets;
  const TestOffsets & testLocalOffsets;
  const TestOffsets & testGlobalOffsets;
  const SolutionLocalIndexSets & solutionLocalIndexSets;
  const SolutionOffsets & solutionLocalOffsets;
  const SolutionOffsets & solutionGlobalOffsets;
};

template<class Indices,
         class LocalMatrix, class GlobalMatrix,
         class TestLocalIndexSets, class SolutionLocalIndexSets,
         class TestOffsets, class SolutionOffsets,
         bool mirror = false>
inline void copyLocalToGlobalMatrix(
    const LocalMatrix&            elementMatrix,
    GlobalMatrix&                 matrix,
    const TestLocalIndexSets&     testLocalIndexSets,
    const TestOffsets&            localTestSpaceOffsets,
    const TestOffsets&            globalTestSpaceOffsets,
    const SolutionLocalIndexSets& solutionLocalIndexSets,
    const SolutionOffsets&        localSolutionSpaceOffsets,
    const SolutionOffsets&        globalSolutionSpaceOffsets) {

  auto cpMatrix = localToGlobalCopier<LocalMatrix, GlobalMatrix,
                 TestLocalIndexSets, SolutionLocalIndexSets,
                 TestOffsets, SolutionOffsets,
                 mirror> (elementMatrix, matrix,
                          testLocalIndexSets,
                          localTestSpaceOffsets,
                          globalTestSpaceOffsets,
                          solutionLocalIndexSets,
                          localSolutionSpaceOffsets,
                          globalSolutionSpaceOffsets);

  /* copy every local submatrix indexed by a pair of indices from
   * Indices exactly once. */
  boost::hana::for_each(Indices{}, cpMatrix);
}

template<class Indices,
         class LocalMatrix, class GlobalMatrix,
         class TestLocalIndexSets, class SolutionLocalIndexSets,
         class TestOffsets, class SolutionOffsets>
inline void copyLocalToGlobalMatrixSymmetric(
    const LocalMatrix&            elementMatrix,
    GlobalMatrix&                 matrix,
    const TestLocalIndexSets&     testLocalIndexSets,
    const TestOffsets&            localTestSpaceOffsets,
    const TestOffsets&            globalTestSpaceOffsets,
    const SolutionLocalIndexSets& solutionLocalIndexSets,
    const SolutionOffsets&        localSolutionSpaceOffsets,
    const SolutionOffsets&        globalSolutionSpaceOffsets) {
  copyLocalToGlobalMatrix<Indices, LocalMatrix, GlobalMatrix,
                          TestLocalIndexSets, SolutionLocalIndexSets,
                          TestOffsets, SolutionOffsets,
                          true>
                         ( elementMatrix,
                           matrix,
                           testLocalIndexSets,
                           localTestSpaceOffsets,
                           globalTestSpaceOffsets,
                           solutionLocalIndexSets,
                           localSolutionSpaceOffsets,
                           globalSolutionSpaceOffsets);
}

template<class LocalVector, class GlobalVector,
         class TestLocalIndexSets, class TestOffsets>
struct localToGlobalRHSCopier
{

  localToGlobalRHSCopier
         (const LocalVector& lv, GlobalVector& gv,
          TestLocalIndexSets const & testLocalIndexSets,
          TestOffsets const & testLocalOffsets,
          TestOffsets const & testGlobalOffsets)
      : localRhs(lv), rhs(gv),
        testLocalIndexSets(testLocalIndexSets),
        testLocalOffsets(testLocalOffsets),
        testGlobalOffsets(testGlobalOffsets) {}

  template <class TestSpaceIndex>
  void operator() (const TestSpaceIndex& index) const
  {
    constexpr auto testSpaceIndex = TestSpaceIndex::type::value;
    const auto& testLocalIndexSet
        = std::get<testSpaceIndex>(testLocalIndexSets);
    const size_t testLocalOffset = testLocalOffsets[testSpaceIndex];
    const size_t testGlobalOffset = testGlobalOffsets[testSpaceIndex];

    using MultiIndex
        = typename std::decay_t<decltype(testLocalIndexSet)>::MultiIndex;

    iterateOverLocalIndexSet(
        testLocalIndexSet,
        [&](size_t i, MultiIndex gi)
        {
          rhs[gi[0] + testGlobalOffset]
              += localRhs[i+testLocalOffset];
        },
        [](size_t i){},
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

  const TestLocalIndexSets & testLocalIndexSets;
  const TestOffsets & testLocalOffsets;
  const TestOffsets & testGlobalOffsets;
};

template<class Indices,
         class LocalVector, class GlobalVector,
         class TestLocalIndexSets,
         class TestOffsets>
inline void copyLocalToGlobalVector(
    const LocalVector&            elementVector,
    GlobalVector&                 vector,
    const TestLocalIndexSets&     testLocalIndexSets,
    const TestOffsets&            localTestSpaceOffsets,
    const TestOffsets&            globalTestSpaceOffsets) {

  auto cpRhs = localToGlobalRHSCopier<LocalVector, GlobalVector,
                    TestLocalIndexSets, TestOffsets>
                         (elementVector, vector,
                          testLocalIndexSets,
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
