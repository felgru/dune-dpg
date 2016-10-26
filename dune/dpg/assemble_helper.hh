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
#include <boost/fusion/algorithm/transformation/join.hpp>
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

template<class E>
struct applyBind
{
  applyBind(const E& e) : e(e) {}

  template<class T>
  void operator()(T& t) const
  {
    t.bind(e);
  }

private:
  const E& e;
};

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
  typedef typename ForEachType<detail::getLocalViewFunctor::TypeEvaluator,
                               TestSpaces>::Type  TestLocalViews;
  //! tuple type for the local views of the solution spaces
  typedef typename ForEachType<detail::getLocalViewFunctor::TypeEvaluator,
                               SolutionSpaces>::Type  SolutionLocalViews;

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
    size_t localTestSpaceOffset =
        localTestSpaceOffsets[testSpaceIndex::value];
    size_t localSolutionSpaceOffset =
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
  typedef typename ForEachType<detail::getLocalViewFunctor::TypeEvaluator,
                               TestSpaces>::Type  TestLocalViews;

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
         (const std::tuple<
          testSpaceIndex,
          solutionSpaceIndex>& indexTuple)
  {
    const auto& testLV =
        std::get<testSpaceIndex::value>(testLocalViews);
    const auto& solutionLV =
        std::get<solutionSpaceIndex::value>(solutionLocalViews);
    const auto& testLIS =
        std::get<testSpaceIndex::value>(testLocalIndexSets);
    const auto& solutionLIS =
        std::get<solutionSpaceIndex::value>(solutionLocalIndexSets);
    size_t globalTestSpaceOffset =
        globalTestSpaceOffsets[testSpaceIndex::value];
    size_t globalSolutionSpaceOffset =
        globalSolutionSpaceOffsets[solutionSpaceIndex::value];

    for (size_t i=0, i_max=testLV.size(); i<i_max; i++) {

      auto iIdx = testLIS.index(i)[0];

      for (size_t j=0, j_max=solutionLV.size(); j<j_max; j++) {

        auto jIdx = solutionLIS.index(j)[0];

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

template<class LocalMatrix, class GlobalMatrix,
         bool mirror = false>
struct localToGlobalCopier
{

  localToGlobalCopier(const LocalMatrix& lm, GlobalMatrix& gm) :
      elementMatrix(lm), matrix(gm) {}

  template<class TestLocalView, class SolutionLocalView,
           class TestLocalIndexSet, class SolutionLocalIndexSet>
  void operator()
         (TestLocalView const & testLocalView,
          TestLocalIndexSet const & testLocalIndexSet,
          size_t testLocalOffset, size_t testGlobalOffset,
          SolutionLocalView const & solutionLocalView,
          SolutionLocalIndexSet const & solutionLocalIndexSet,
          size_t solutionLocalOffset, size_t solutionGlobalOffset
         )
  {
    const size_t nTest(testLocalView.size());
    const size_t nSolution(solutionLocalView.size());

    for (size_t i=0; i<nTest; i++)
    {
      auto row = testLocalIndexSet.index(i)[0]+testGlobalOffset;

      for (size_t j=0; j<nSolution; j++)
      {
        auto col = solutionLocalIndexSet.index(j)[0]
                    +solutionGlobalOffset;
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
};

template <class SolutionZip,
          class TestZip,
          class LocalToGlobalCopier>
struct localToGlobalCopyHelper
{
  localToGlobalCopyHelper(
                       const SolutionZip& solutionZip,
                       const TestZip& testZip,
                       LocalToGlobalCopier& localToGlobalCopier)
      : solutionZip(solutionZip),
        testZip(testZip),
        localToGlobalCopier(localToGlobalCopier)
  {}

  template <class testSpaceIndex,
            class solutionSpaceIndex>
  void operator()
         (const std::tuple<
          testSpaceIndex,
          solutionSpaceIndex>& indexTuple) const
  {
    using namespace boost::fusion;

    const auto& solutionData =
        at_c<solutionSpaceIndex::value>(solutionZip);
    const auto& testData =
        at_c<testSpaceIndex::value>(testZip);

    localToGlobalCopier(join(testData,solutionData));
  }

private:
  const SolutionZip& solutionZip;
  const TestZip& testZip;

  LocalToGlobalCopier& localToGlobalCopier;
};

template<class LocalVector, class GlobalVector>
struct localToGlobalRHSCopier
{

  localToGlobalRHSCopier(const LocalVector& lv, GlobalVector& gv) :
      localRhs(lv), rhs(gv) {}

  template<class TestLocalView, class TestLocalIndexSet>
  void operator()
         (TestLocalView const & testLocalView,
          TestLocalIndexSet const & testLocalIndexSet,
          size_t testLocalOffset, size_t testGlobalOffset
         )
  {
    const size_t nTest(testLocalView.size());

    for (size_t i=0; i<nTest; i++) {
      auto row = testLocalIndexSet.index(i)[0] + testGlobalOffset;
      rhs[row] += localRhs[i+testLocalOffset];
    }
  }

private:
  const LocalVector& localRhs;
  GlobalVector& rhs;
};

template <class TestZip,
          class LocalToGlobalRHSCopier>
struct localToGlobalRHSCopyHelper
{
  localToGlobalRHSCopyHelper(
                       const TestZip& testZip,
                       LocalToGlobalRHSCopier& localToGlobalRHSCopier)
      : testZip(testZip),
        localToGlobalRHSCopier(localToGlobalRHSCopier)
  {}

  template <class TestSpaceIndex>
  void operator()
         (const TestSpaceIndex& index) const
  {
    using namespace boost::fusion;

    const auto& testData =
        at_c<TestSpaceIndex::value>(testZip);

    localToGlobalRHSCopier(testData);
  }

private:
  const TestZip& testZip;

  LocalToGlobalRHSCopier& localToGlobalRHSCopier;
};

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

template<class VectorType1, class VectorType2>
struct getLocalCoefficients
{
  getLocalCoefficients(const VectorType1& sol, VectorType2& solEl) :
    solution(sol),
    solutionElement(solEl)
  {}

  template<class T>
  void operator()(const T& t) const
  {
    using namespace boost::fusion;

    auto const & localView = at_c<0>(t);
    auto const & localIndexSet = at_c<1>(t);
    size_t dofElement = localView.size();
    for (size_t i=0; i<dofElement; i++)
    {
      solutionElement[i+at_c<2>(t)] = solution[at_c<3>(t)+localIndexSet.index(i)[0] ];
    }
  }
private:
  const VectorType1& solution;
  VectorType2& solutionElement;
};



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
