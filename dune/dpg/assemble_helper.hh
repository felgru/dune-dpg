// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_ASSEMBLE_HELPER_HH
#define DUNE_DPG_ASSEMBLE_HELPER_HH

#include <functional>
#include <utility>
#include <tuple>
#include <dune/common/hybridutilities.hh>
#include <dune/common/std/memory.hh>
#include <dune/istl/matrixindexset.hh>
#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/fusion/adapted/array.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/algorithm/transformation/zip.hpp>
#include <boost/fusion/sequence/intrinsic/at.hpp>

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

template <class TestLocalViews,
          class SolutionLocalViews,
          class TestLocalIndexSets,
          class SolutionLocalIndexSets,
          bool mirror = false>
struct getOccupationPatternHelper
{
  template<class T, class Tuple>
  using array_of_same_size = T[std::tuple_size<Tuple>::value];

  getOccupationPatternHelper(
                       const TestLocalViews& testLocalViews,
                       const SolutionLocalViews& solutionLocalViews,
                       const TestLocalIndexSets& testLocalIndexSets,
                       const SolutionLocalIndexSets& solutionLocalIndexSets,
                       const array_of_same_size<size_t, TestLocalViews>&
                           globalTestSpaceOffsets,
                       const array_of_same_size<size_t, SolutionLocalViews>&
                           globalSolutionSpaceOffsets,
                       Dune::MatrixIndexSet& nb)
      : testLocalViews(testLocalViews),
        solutionLocalViews(solutionLocalViews),
        testLocalIndexSets(testLocalIndexSets),
        solutionLocalIndexSets(solutionLocalIndexSets),
        globalTestSpaceOffsets(globalTestSpaceOffsets),
        globalSolutionSpaceOffsets(globalSolutionSpaceOffsets),
        nb(nb)
  {}

  template <class testSpaceIndex,
            class solutionSpaceIndex>
  void operator()
         (const std::tuple<testSpaceIndex, solutionSpaceIndex>& indexTuple)
  {
    const auto& testLV =
        std::get<testSpaceIndex::value>(testLocalViews);
    const auto& solutionLV =
        std::get<solutionSpaceIndex::value>(solutionLocalViews);
    const auto& testLIS =
        std::get<testSpaceIndex::value>(testLocalIndexSets);
    const auto& solutionLIS =
        std::get<solutionSpaceIndex::value>(solutionLocalIndexSets);
    const size_t globalTestSpaceOffset =
        globalTestSpaceOffsets[testSpaceIndex::value];
    const size_t globalSolutionSpaceOffset =
        globalSolutionSpaceOffsets[solutionSpaceIndex::value];

    for (size_t i=0, i_max=testLV.size(); i<i_max; i++) {

      const auto iIdx = testLIS.index(i)[0];

      for (size_t j=0, j_max=solutionLV.size(); j<j_max; j++) {

        const auto jIdx = solutionLIS.index(j)[0];

        // Add a nonzero entry to the matrix
        nb.add(iIdx+globalTestSpaceOffset,
               jIdx+globalSolutionSpaceOffset);
        if(mirror) {
            nb.add(jIdx+globalSolutionSpaceOffset,
                   iIdx+globalTestSpaceOffset);
        }

      }
    }
  }

private:
  const TestLocalViews& testLocalViews;
  const SolutionLocalViews& solutionLocalViews;
  const TestLocalIndexSets& testLocalIndexSets;
  const SolutionLocalIndexSets& solutionLocalIndexSets;
  const array_of_same_size<size_t, TestLocalViews>&
      globalTestSpaceOffsets;
  const array_of_same_size<size_t, SolutionLocalViews>&
      globalSolutionSpaceOffsets;
  Dune::MatrixIndexSet& nb;
};

template<class Indices,
         bool mirror,
         class LeftLocalViews, class RightLocalViews,
         class LeftLocalIndexSets, class RightLocalIndexSets,
         class LeftOffsets, class RightOffsets>
inline void getOccupationPattern(
      const LeftLocalViews& leftLocalViews,
      const RightLocalViews& rightLocalViews,
      const LeftLocalIndexSets& leftLocalIndexSets,
      const RightLocalIndexSets& rightLocalIndexSets,
      const LeftOffsets& leftGlobalOffsets,
      const RightOffsets& rightGlobalOffsets,
      MatrixIndexSet& occupationPattern) {
  auto gOPH = getOccupationPatternHelper<LeftLocalViews,
                                         RightLocalViews,
                                         LeftLocalIndexSets,
                                         RightLocalIndexSets,
                                         mirror>
                      (leftLocalViews,
                       rightLocalViews,
                       leftLocalIndexSets,
                       rightLocalIndexSets,
                       leftGlobalOffsets,
                       rightGlobalOffsets,
                       occupationPattern);
  boost::fusion::for_each(Indices{},
      std::ref(gOPH));
}

template<class LocalMatrix, class GlobalMatrix,
         class TestLocalViews, class SolutionLocalViews,
         class TestLocalIndexSets, class SolutionLocalIndexSets,
         class TestOffsets, class SolutionOffsets,
         bool mirror = false>
struct localToGlobalCopier
{

  localToGlobalCopier
         (const LocalMatrix& lm, GlobalMatrix& gm,
          TestLocalViews const & testLocalViews,
          TestLocalIndexSets const & testLocalIndexSets,
          TestOffsets const & testLocalOffsets,
          TestOffsets const & testGlobalOffsets,
          SolutionLocalViews const & solutionLocalViews,
          SolutionLocalIndexSets const & solutionLocalIndexSets,
          SolutionOffsets const & solutionLocalOffsets,
          SolutionOffsets const & solutionGlobalOffsets)
      : elementMatrix(lm), matrix(gm),
        testLocalViews(testLocalViews),
        testLocalIndexSets(testLocalIndexSets),
        testLocalOffsets(testLocalOffsets),
        testGlobalOffsets(testGlobalOffsets),
        solutionLocalViews(solutionLocalViews),
        solutionLocalIndexSets(solutionLocalIndexSets),
        solutionLocalOffsets(solutionLocalOffsets),
        solutionGlobalOffsets(solutionGlobalOffsets) {}

  template <class testSpaceIndex,
            class solutionSpaceIndex>
  void operator()
         (const std::tuple<
          testSpaceIndex,
          solutionSpaceIndex>& indexTuple) const
  {
    const auto& testLocalView
        = std::get<testSpaceIndex::value>(testLocalViews);
    const auto& testLocalIndexSet
        = std::get<testSpaceIndex::value>(testLocalIndexSets);
    const size_t testLocalOffset = testLocalOffsets[testSpaceIndex::value];
    const size_t testGlobalOffset = testGlobalOffsets[testSpaceIndex::value];
    const auto& solutionLocalView
        = std::get<solutionSpaceIndex::value>(solutionLocalViews);
    const auto& solutionLocalIndexSet
        = std::get<solutionSpaceIndex::value>(solutionLocalIndexSets);
    const size_t solutionLocalOffset
        = solutionLocalOffsets[solutionSpaceIndex::value];
    const size_t solutionGlobalOffset
        = solutionGlobalOffsets[solutionSpaceIndex::value];

    const size_t nTest(testLocalView.size());
    const size_t nSolution(solutionLocalView.size());

    for (size_t i=0; i<nTest; i++)
    {
      const auto row = testLocalIndexSet.index(i)[0] + testGlobalOffset;

      for (size_t j=0; j<nSolution; j++)
      {
        const auto col = solutionLocalIndexSet.index(j)[0]
                         + solutionGlobalOffset;
        matrix[row][col] += elementMatrix[i+testLocalOffset]
                                         [j+solutionLocalOffset];
        if(mirror) {
          matrix[col][row] += elementMatrix[i+testLocalOffset]
                                           [j+solutionLocalOffset];
        }
      }
    }
  }

private:
  const LocalMatrix& elementMatrix;
  GlobalMatrix& matrix;

  const TestLocalViews & testLocalViews;
  const TestLocalIndexSets & testLocalIndexSets;
  const TestOffsets & testLocalOffsets;
  const TestOffsets & testGlobalOffsets;
  const SolutionLocalViews & solutionLocalViews;
  const SolutionLocalIndexSets & solutionLocalIndexSets;
  const SolutionOffsets & solutionLocalOffsets;
  const SolutionOffsets & solutionGlobalOffsets;
};

template<class Indices,
         class LocalMatrix, class GlobalMatrix,
         class TestLocalViews, class SolutionLocalViews,
         class TestLocalIndexSets, class SolutionLocalIndexSets,
         class TestOffsets, class SolutionOffsets,
         bool mirror = false>
inline void copyLocalToGlobalMatrix(
    const LocalMatrix&            elementMatrix,
    GlobalMatrix&                 matrix,
    const TestLocalViews&         testLocalViews,
    const TestLocalIndexSets&     testLocalIndexSets,
    const TestOffsets&            localTestSpaceOffsets,
    const TestOffsets&            globalTestSpaceOffsets,
    const SolutionLocalViews&     solutionLocalViews,
    const SolutionLocalIndexSets& solutionLocalIndexSets,
    const SolutionOffsets&        localSolutionSpaceOffsets,
    const SolutionOffsets&        globalSolutionSpaceOffsets) {

  auto cpMatrix = localToGlobalCopier<LocalMatrix, GlobalMatrix,
                 TestLocalViews, SolutionLocalViews,
                 TestLocalIndexSets, SolutionLocalIndexSets,
                 TestOffsets, SolutionOffsets,
                 mirror> (elementMatrix, matrix,
                          testLocalViews,
                          testLocalIndexSets,
                          localTestSpaceOffsets,
                          globalTestSpaceOffsets,
                          solutionLocalViews,
                          solutionLocalIndexSets,
                          localSolutionSpaceOffsets,
                          globalSolutionSpaceOffsets);

  /* copy every local submatrix indexed by a pair of indices from
   * Indices exactly once. */
  boost::fusion::for_each(Indices{}, cpMatrix);
}

template<class Indices,
         class LocalMatrix, class GlobalMatrix,
         class TestLocalViews, class SolutionLocalViews,
         class TestLocalIndexSets, class SolutionLocalIndexSets,
         class TestOffsets, class SolutionOffsets>
inline void copyLocalToGlobalMatrixSymmetric(
    const LocalMatrix&            elementMatrix,
    GlobalMatrix&                 matrix,
    const TestLocalViews&         testLocalViews,
    const TestLocalIndexSets&     testLocalIndexSets,
    const TestOffsets&            localTestSpaceOffsets,
    const TestOffsets&            globalTestSpaceOffsets,
    const SolutionLocalViews&     solutionLocalViews,
    const SolutionLocalIndexSets& solutionLocalIndexSets,
    const SolutionOffsets&        localSolutionSpaceOffsets,
    const SolutionOffsets&        globalSolutionSpaceOffsets) {
  copyLocalToGlobalMatrix<Indices, LocalMatrix, GlobalMatrix,
                          TestLocalViews, SolutionLocalViews,
                          TestLocalIndexSets, SolutionLocalIndexSets,
                          TestOffsets, SolutionOffsets,
                          true>
                         ( elementMatrix,
                           matrix,
                           testLocalViews,
                           testLocalIndexSets,
                           localTestSpaceOffsets,
                           globalTestSpaceOffsets,
                           solutionLocalViews,
                           solutionLocalIndexSets,
                           localSolutionSpaceOffsets,
                           globalSolutionSpaceOffsets);
}

template<class LocalVector, class GlobalVector,
         class TestLocalViews, class TestLocalIndexSets,
         class TestOffsets>
struct localToGlobalRHSCopier
{

  localToGlobalRHSCopier
         (const LocalVector& lv, GlobalVector& gv,
          TestLocalViews const & testLocalViews,
          TestLocalIndexSets const & testLocalIndexSets,
          TestOffsets const & testLocalOffsets,
          TestOffsets const & testGlobalOffsets)
      : localRhs(lv), rhs(gv),
        testLocalViews(testLocalViews),
        testLocalIndexSets(testLocalIndexSets),
        testLocalOffsets(testLocalOffsets),
        testGlobalOffsets(testGlobalOffsets) {}

  template <class TestSpaceIndex>
  void operator() (const TestSpaceIndex& index) const
  {
    const auto& testLocalView
        = std::get<TestSpaceIndex::value>(testLocalViews);
    const auto& testLocalIndexSet
        = std::get<TestSpaceIndex::value>(testLocalIndexSets);
    const size_t testLocalOffset = testLocalOffsets[TestSpaceIndex::value];
    const size_t testGlobalOffset = testGlobalOffsets[TestSpaceIndex::value];

    const size_t nTest(testLocalView.size());

    for (size_t i=0; i<nTest; i++) {
      const auto row = testLocalIndexSet.index(i)[0] + testGlobalOffset;
      rhs[row] += localRhs[i+testLocalOffset];
    }
  }

private:
  const LocalVector& localRhs;
  GlobalVector& rhs;

  const TestLocalViews & testLocalViews;
  const TestLocalIndexSets & testLocalIndexSets;
  const TestOffsets & testLocalOffsets;
  const TestOffsets & testGlobalOffsets;
};

template<class Indices,
         class LocalVector, class GlobalVector,
         class TestLocalViews,
         class TestLocalIndexSets,
         class TestOffsets>
inline void copyLocalToGlobalVector(
    const LocalVector&            elementVector,
    GlobalVector&                 vector,
    const TestLocalViews&         testLocalViews,
    const TestLocalIndexSets&     testLocalIndexSets,
    const TestOffsets&            localTestSpaceOffsets,
    const TestOffsets&            globalTestSpaceOffsets) {

  auto cpRhs = localToGlobalRHSCopier<LocalVector, GlobalVector,
                    TestLocalViews, TestLocalIndexSets,
                    TestOffsets>
                         (elementVector, vector,
                          testLocalViews,
                          testLocalIndexSets,
                          localTestSpaceOffsets,
                          globalTestSpaceOffsets);

  /* copy every local subvector indexed by an index from
   * Indices exactly once. */
  boost::fusion::for_each(Indices{}, cpRhs);
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


namespace mpl {
  template<class Seq>
  struct first
  {
    typedef typename std::remove_reference<
      typename boost::fusion::result_of::at_c<Seq,0>::type>::type
        type;
  };

  template<class Seq>
  struct second
  {
    typedef typename std::remove_reference<
      typename boost::fusion::result_of::at_c<Seq,1>::type>::type
        type;
  };

  template<class Seq>
  struct firstTwo
  {
    template <size_t i>
    using type_at = typename std::remove_reference<
        typename boost::fusion::result_of::at_c<Seq,i>::type>::type;

    typedef std::tuple<type_at<0>, type_at<1> > type;
  };

  template<class T>
  struct tupleOf0And
  {
    typedef std::tuple<std::integral_constant<size_t,0>, T> type;
  };

  template<class I, class J>
  struct tupleOfIAndJ
  {
    typedef std::tuple<I, J> type;
  };

  template<class Seq, class I>
  struct prefixSequenceWith
  {
    typedef
        typename boost::mpl::transform<
            Seq
          , tupleOfIAndJ<I, boost::mpl::_1>
          >::type type;
  };
}

} } // end namespace Dune::detail

#endif // DUNE_DPG_ASSEMBLE_HELPER_HH
