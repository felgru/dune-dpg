// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_ASSEMBLE_HELPER_HH
#define DUNE_DPG_ASSEMBLE_HELPER_HH

#include <functional>
#include <utility>
#include <tuple>
#include <dune/common/std/memory.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <boost/fusion/container/vector/convert.hpp>
#include <boost/fusion/sequence/intrinsic/size.hpp>
#include <boost/fusion/algorithm/transformation/transform.hpp>

namespace Dune {

namespace detail {

struct getLocalIndexSet
{
  template<class T>
  struct result;

  template<class T>
  struct result<getLocalIndexSet(T)>
  {
    typedef typename T::LocalIndexSet type;
  };

  template<class T>
  typename result<getLocalIndexSet(T)>::type operator()(const T& t) const
  {
    return t.localIndexSet();
  }
};

struct getLocalView
{
  template<class T>
  struct result;

  template<class T>
  struct result<getLocalView(T)>
  {
    typedef typename T::LocalView type;
  };

  template<class T>
  typename result<getLocalView(T)>::type operator()(const T& t) const
  {
    return t.localView();
  }
};

struct getSize
{
  template<class T>
  struct result;

  template<class T>
  struct result<getSize(T)>
  {
    typedef size_t type;
  };

  template<class T>
  size_t operator()(const T& t) const
  {
    return t.size();
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

struct bindLocalIndexSet
{
  template<class LIS, class LV>
  void operator()(const LIS& lis, const LV& lv) const
  {
    /* TODO: I feel uncomfortable casting away the const, but
     * I do not know how else to work around the fact that many
     * boost::fusion functions only take const sequences. */
    const_cast<LIS&>(lis).bind(lv);
  }
};

template <class MatrixType,
          class TestSpaces,
          class SolutionSpaces>
struct getLocalMatrixHelper
{
  template<class T, class Seq>
  using array_of_same_size =
      T[boost::fusion::result_of::size<Seq>::type::value];

  //! tuple type for the local views of the test spaces
  typedef typename boost::fusion::result_of::as_vector<
      typename boost::fusion::result_of::
      transform<TestSpaces, detail::getLocalView>::type>::type TestLocalViews;
  //! tuple type for the local views of the solution spaces
  typedef typename boost::fusion::result_of::as_vector<
      typename boost::fusion::result_of::
      transform<SolutionSpaces, detail::getLocalView>::type
      >::type SolutionLocalViews;

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
    using namespace boost::fusion;

    const auto& term = std::get<2>(termTuple);

    const auto& testLV =
        at_c<testSpaceIndex::value>(testLocalViews);
    const auto& solutionLV =
        at_c<solutionSpaceIndex::value>(solutionLocalViews);
    size_t localTestSpaceOffset =
        at_c<testSpaceIndex::value>(localTestSpaceOffsets);
    size_t localSolutionSpaceOffset =
        at_c<solutionSpaceIndex::value>(localSolutionSpaceOffsets);

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

template <class TestLocalViews,
          class SolutionLocalViews,
          class TestLocalIndexSets,
          class SolutionLocalIndexSets,
          bool mirror = false>
struct getOccupationPatternHelper
{
  template<class T, class Seq>
  using array_of_same_size =
      T[boost::fusion::result_of::size<Seq>::type::value];

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
    using namespace boost::fusion;

    const auto& testLV =
        at_c<testSpaceIndex::value>(testLocalViews);
    const auto& solutionLV =
        at_c<solutionSpaceIndex::value>(solutionLocalViews);
    const auto& testLIS =
        at_c<testSpaceIndex::value>(testLocalIndexSets);
    const auto& solutionLIS =
        at_c<solutionSpaceIndex::value>(solutionLocalIndexSets);
    size_t globalTestSpaceOffset =
        at_c<testSpaceIndex::value>(globalTestSpaceOffsets);
    size_t globalSolutionSpaceOffset =
        at_c<solutionSpaceIndex::value>(globalSolutionSpaceOffsets);

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

template<class GlobalVector>
struct localToGlobalRHSCopier
{

  localToGlobalRHSCopier(GlobalVector& gv) :
      rhs(gv) {}

  template<class LocalVector, class TestLocalIndexSet>
  void operator()
         (LocalVector& localRhs,
          TestLocalIndexSet const & testLocalIndexSet,
          size_t testGlobalOffset
         )
  {
    for (size_t i=0, i_max=localRhs.size(); i<i_max; i++) {
      // The global index of the i-th vertex of the element 'it'
      auto row = testLocalIndexSet.index(i)[0]
                 + testGlobalOffset;
      rhs[row] += localRhs[i];
    }
  }

private:
  GlobalVector& rhs;
};

struct offsetHelper
{
  template<class T>
  struct result;

  template<class T>
  struct result<offsetHelper(const size_t&, T)>
  {
    typedef size_t type;
  };

  template<class T>
  size_t operator()(const size_t& s, const T& t) const
  {
    using namespace boost::fusion;

    // TODO: böser const_cast!
    //       Can we put offset into a reference_wrappers to fix this?
    size_t & offset = const_cast<size_t&>(at_c<0>(t));
    auto const & localView = at_c<1>(t);
    offset = s;
    return s + localView.size();
  }
};

struct globalOffsetHelper
{
  template<class T>
  struct result;

  template<class T>
  struct result<globalOffsetHelper(const size_t&, T)>
  {
    typedef size_t type;
  };

  template<class T>
  size_t operator()(const size_t& s, const T& t) const
  {
    using namespace boost::fusion;

    // TODO: böser const_cast!
    //       Can we put offset into a reference_wrappers to fix this?
    size_t & offset = const_cast<size_t&>(at_c<0>(t));
    auto const & space = at_c<1>(t);
    offset = s;
    return s + space.size();
  }
};


struct getMaxNodeSize
{
    template<class T>
    struct result;

    template<class T>
    struct result<getMaxNodeSize(T)>
    {
        typedef std::size_t type;
    };
    template<class T>
    std::size_t operator()(T t) const
    {
        return t.nodeFactory().maxNodeSize();
    }
};



struct applyUnbind
{
    template<class T>
    void operator()(T t) const
    {
        t.unbind();
    }
};


struct getLocalFiniteElement
{
    template<class T>
    struct result;

    template<class T>
    struct result<getLocalFiniteElement(T)>
    {
        typedef const typename T::Tree::FiniteElement& type;
    };

    template<class T>
    typename result<getLocalFiniteElement(T)>::type operator()(const T& t) const
    {
        return t.tree().finiteElement();
    }
};


namespace mpl {
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
}

} } // end namespace Dune::detail

#endif // DUNE_DPG_ASSEMBLE_HELPER_HH
