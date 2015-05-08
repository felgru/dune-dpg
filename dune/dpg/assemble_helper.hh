// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_ASSEMBLE_HELPER_HH
#define DUNE_DPG_ASSEMBLE_HELPER_HH

#include <functional>
#include <utility>
#include <dune/common/std/memory.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <boost/fusion/sequence/intrinsic/size.hpp>
#include <boost/fusion/algorithm/transformation/pop_back.hpp>

namespace Dune {

namespace detail {

struct getIndexSet
{
  template<class T>
  struct result;

  template<class T>
  struct result<getIndexSet(T)>
  {
    typedef typename T::IndexSet type;
  };

  template<class T>
  typename T::IndexSet operator()(T t) const
  {
    return t.indexSet();
  };
};

struct getLocalIndexSet
{
  template<class T>
  struct result;

  template<class T>
  struct result<getLocalIndexSet(T)>
  {
    typedef typename std::add_pointer<typename T::LocalIndexSet>::type type;
  };

  template<class T>
  typename result<getLocalIndexSet(T)>::type operator()(const T& t) const
  {
    const typename T::LocalIndexSet& lis =
        const_cast<T&>(t).localIndexSet();
    typename T::LocalIndexSet& result =
        const_cast<typename T::LocalIndexSet&>(lis);
    return new typename T::LocalIndexSet(std::move(result));
  };
};

struct getLocalView
{
  template<class T>
  struct result;

  template<class T>
  struct result<getLocalView(T)>
  {
    typedef typename T::LocalView* type;
  };

  template<class T>
  typename result<getLocalView(T)>::type operator()(const T& t) const
  {
    const typename T::LocalView& lv = const_cast<T&>(t).localView();
    typename T::LocalView& result = const_cast<typename T::LocalView&>(lv);
    return new typename T::LocalView(std::move(result));
  };
};

struct localViewFromFEBasis
{
  template<class R>
  struct result;

  template<class B>
  struct result<localViewFromFEBasis(B)>
  {
    typedef typename B::LocalView* type;
  };

  template<class B>
  typename result<localViewFromFEBasis(B)>::type
  operator()(const B& b) const
  {
    const typename B::LocalView& lv = b.localView();
    typename B::LocalView& result =
        const_cast<typename B::LocalView&>(lv);
    return new typename B::LocalView(std::move(result));
  };
};

template<class GridView>
struct getLocalVolumeTerm
{
  getLocalVolumeTerm(const GridView& gridView) : gridView(gridView) {};

  template<class T>
  struct result;

  template<class T>
  struct result<getLocalVolumeTerm(T)>
  {
    typedef decltype(localFunction(
                     Functions::makeGridViewFunction(
                             std::declval<T>(),
                             std::declval<GridView>()))) type;
  };

  template<class T>
  typename result<getLocalVolumeTerm(T)>::type operator()(const T& t) const
  {
    /* TODO: correctly forward t to makeGridViewFunction.
     *       This also needs adaptation at the calling site */
    auto localVolumeTerm =
        localFunction(Functions::makeGridViewFunction(t, gridView));
    return localVolumeTerm;
  };

  private:
  const GridView& gridView;
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
  };
};

struct default_deleter
{
  template<class T>
  void operator()(T* t) const
  {
    delete t;
  }
};

template<class E>
struct applyBind
{
  applyBind(const E& e) : e(e) {};

  /* TODO: Is the T& more expensive than T for pointer types T? */
  template<class T>
  void operator()(T& t) const
  {
    bind_impl(t, std::is_pointer<T>());
  }

  template<class T>
  void bind_impl(T t, std::true_type) const
  { t->bind(e); };

  template<class T>
  void bind_impl(T& t, std::false_type) const
  { t.bind(e); };

private:
  const E& e;
};

struct bindLocalIndexSet
{
  template<class LIS, class LV>
  void operator()(const LIS& lis, const LV* lv) const
  {
    /* TODO: I feel uncomfortable casting away the const, but
     * I do not know how else to work around the fact that many
     * boost::fusion functions only take const sequences. */
    const_cast<LIS&>(lis)->bind(*lv);
  }
};

template <class MatrixType,
          class TestLocalView,
          class SolutionLocalView>
struct getLocalMatrixHelper
{
  template<class T, class Seq>
  using array_of_same_size =
      T[boost::fusion::result_of::size<Seq>::type::value];

  getLocalMatrixHelper(const TestLocalView& testLocalView,
                       const SolutionLocalView& solutionLocalView,
                       MatrixType& elementMatrix,
                       const array_of_same_size<size_t, TestLocalView>&
                           localTestSpaceOffsets,
                       const array_of_same_size<size_t, SolutionLocalView>&
                           localSolutionSpaceOffsets)
      : solutionLocalView(solutionLocalView),
        testLocalView(testLocalView),
        elementMatrix(elementMatrix),
        localSolutionSpaceOffsets(localSolutionSpaceOffsets),
        localTestSpaceOffsets(localTestSpaceOffsets)
  {};

  /**
   * \tparam Term either a BilinearTerm or an InnerProductTerm
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
        at_c<testSpaceIndex::value>(testLocalView);
    const auto& solutionLV =
        at_c<solutionSpaceIndex::value>(solutionLocalView);
    size_t localTestSpaceOffset =
        at_c<testSpaceIndex::value>(localTestSpaceOffsets);
    size_t localSolutionSpaceOffset =
        at_c<solutionSpaceIndex::value>(localSolutionSpaceOffsets);
    term.getLocalMatrix(testLV,
                        solutionLV,
                        elementMatrix,
                        localTestSpaceOffset,
                        localSolutionSpaceOffset);
  };

private:
  const TestLocalView& testLocalView;
  const SolutionLocalView& solutionLocalView;
  MatrixType& elementMatrix;
  const array_of_same_size<size_t, TestLocalView>&
      localTestSpaceOffsets;
  const array_of_same_size<size_t, SolutionLocalView>&
      localSolutionSpaceOffsets;
};

template <class TestLocalView,
          class SolutionLocalView,
          class TestLocalIndexSet,
          class SolutionLocalIndexSet,
          bool mirror = false>
struct getOccupationPatternHelper
{
  template<class T, class Seq>
  using array_of_same_size =
      T[boost::fusion::result_of::size<Seq>::type::value];

  getOccupationPatternHelper(
                       const TestLocalView& testLocalView,
                       const SolutionLocalView& solutionLocalView,
                       const TestLocalIndexSet& testLocalIndexSet,
                       const SolutionLocalIndexSet& solutionLocalIndexSet,
                       const array_of_same_size<size_t, TestLocalView>&
                           globalTestSpaceOffsets,
                       const array_of_same_size<size_t, SolutionLocalView>&
                           globalSolutionSpaceOffsets,
                       Dune::MatrixIndexSet& nb)
      : testLocalView(testLocalView),
        solutionLocalView(solutionLocalView),
        testLocalIndexSet(testLocalIndexSet),
        solutionLocalIndexSet(solutionLocalIndexSet),
        globalTestSpaceOffsets(globalTestSpaceOffsets),
        globalSolutionSpaceOffsets(globalSolutionSpaceOffsets),
        nb(nb)
  {};

  template <class testSpaceIndex,
            class solutionSpaceIndex>
  void operator()
         (const std::tuple<
          testSpaceIndex,
          solutionSpaceIndex>& indexTuple)
  {
    using namespace boost::fusion;

    const auto& testLV =
        at_c<testSpaceIndex::value>(testLocalView);
    const auto& solutionLV =
        at_c<solutionSpaceIndex::value>(solutionLocalView);
    const auto& testLIS =
        at_c<testSpaceIndex::value>(testLocalIndexSet);
    const auto& solutionLIS =
        at_c<solutionSpaceIndex::value>(solutionLocalIndexSet);
    size_t globalTestSpaceOffset =
        at_c<testSpaceIndex::value>(globalTestSpaceOffsets);
    size_t globalSolutionSpaceOffset =
        at_c<solutionSpaceIndex::value>(globalSolutionSpaceOffsets);

    for (size_t i=0, i_max=testLV->size(); i<i_max; i++) {

      auto iIdx = testLIS->index(i)[0];

      for (size_t j=0, j_max=solutionLV->size(); j<j_max; j++) {

        auto jIdx = solutionLIS->index(j)[0];

        // Add a nonzero entry to the matrix
        nb.add(iIdx+globalTestSpaceOffset,
               jIdx+globalSolutionSpaceOffset);
        if(mirror) {
            nb.add(jIdx+globalSolutionSpaceOffset,
                   iIdx+globalTestSpaceOffset);
        }

      }
    }
  };

private:
  const TestLocalView& testLocalView;
  const SolutionLocalView& solutionLocalView;
  const TestLocalIndexSet& testLocalIndexSet;
  const SolutionLocalIndexSet& solutionLocalIndexSet;
  const array_of_same_size<size_t, TestLocalView>&
      globalTestSpaceOffsets;
  const array_of_same_size<size_t, SolutionLocalView>&
      globalSolutionSpaceOffsets;
  Dune::MatrixIndexSet& nb;
};

template<class LocalMatrix, class GlobalMatrix,
         bool mirror = false>
struct localToGlobalCopier
{

  localToGlobalCopier(const LocalMatrix& lm, GlobalMatrix& gm) :
      elementMatrix(lm), matrix(gm) {};

  template<class TestLocalView, class SolutionLocalView,
           class TestLocalIndexSet, class SolutionLocalIndexSet>
  void operator()
         (TestLocalView const * testLocalView,
          TestLocalIndexSet const * testLocalIndexSet,
          size_t testLocalOffset, size_t testGlobalOffset,
          SolutionLocalView const * solutionLocalView,
          SolutionLocalIndexSet const * solutionLocalIndexSet,
          size_t solutionLocalOffset, size_t solutionGlobalOffset
         )
  {
    const size_t nTest(testLocalView->size());
    const size_t nSolution(solutionLocalView->size());

    for (size_t i=0; i<nTest; i++)
    {
      auto row = testLocalIndexSet->index(i)[0]+testGlobalOffset;

      for (size_t j=0; j<nSolution; j++)
      {
        auto col = solutionLocalIndexSet->index(j)[0]
                    +solutionGlobalOffset;
        matrix[row][col] += elementMatrix[i+testLocalOffset]
                                         [j+solutionLocalOffset];
        if(mirror) {
          matrix[col][row] += elementMatrix[i+testLocalOffset]
                                           [j+solutionLocalOffset];
        }
      }
    }
  };

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
                       const std::reference_wrapper<LocalToGlobalCopier>&
                           localToGlobalCopier)
      : solutionZip(solutionZip),
        testZip(testZip),
        localToGlobalCopier(localToGlobalCopier)
  {};

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
  };

private:
  const SolutionZip& solutionZip;
  const TestZip& testZip;

  const std::reference_wrapper<LocalToGlobalCopier>& localToGlobalCopier;
};

template<class GlobalVector>
struct localToGlobalRHSCopier
{

  localToGlobalRHSCopier(GlobalVector& gv) :
      rhs(gv) {};

  template<class LocalVector, class TestLocalIndexSet>
  void operator()
         (LocalVector& localRhs,
          TestLocalIndexSet const * testLocalIndexSet,
          size_t testGlobalOffset
         )
  {
    for (size_t i=0, i_max=localRhs.size(); i<i_max; i++) {
      // The global index of the i-th vertex of the element 'it'
      auto row = testLocalIndexSet->index(i)[0]
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
  struct result<offsetHelper(size_t,T)>
  {
    typedef size_t type;
  };

  template<class T>
  size_t operator()(size_t s, const T& t) const
  {
    using namespace boost::fusion;

    /* offset and localView are assumed to be reference_wrappers */
    size_t & offset = const_cast<size_t&>(at_c<0>(t));
    auto const & localView = at_c<1>(t);
    offset = s;
    return s + localView->size();
  };
};

struct globalOffsetHelper
{
  template<class T>
  struct result;

  template<class T>
  struct result<offsetHelper(size_t,T)>
  {
    typedef size_t type;
  };

  template<class T>
  size_t operator()(size_t s, const T& t) const
  {
    using namespace boost::fusion;

    /* offset and localView are assumed to be reference_wrappers */
    size_t & offset = const_cast<size_t&>(at_c<0>(t));
    auto const & indexSet = at_c<1>(t);
    offset = s;
    return s + indexSet.size();
  };
};


struct getLocalFeSize
{
    template<class T>
    struct result;

    template<class T>
    struct result<getLocalFeSize(T)>
    {
        typedef std::size_t type;
    };
    template<class T>
    std::size_t operator()(T t) const
    {
        return t->tree().finiteElement().size();
    };
};



struct getIndexSetSize
{
    template<class T>
    struct result;

    template<class T>
    struct result<getIndexSetSize(T)>
    {
        typedef std::size_t type;
    };
    template<class T>
    std::size_t operator()(T t) const
    {
        return t.indexSet().size();
    };
};



struct getLocalViewMaxSize
{
    template<class T>
    struct result;

    template<class T>
    struct result<getLocalViewMaxSize(T)>
    {
        typedef std::size_t type;
    };
    template<class T>
    std::size_t operator()(T t) const
    {
        return t.localView().maxSize();
    };
};



struct applyUnbind
{
    template<class T>
    void operator()(T t) const
    {
        t->unbind();
    }
};


struct getLocalFiniteElement
{
    template<class T>
    struct result;

    template<class T>
    struct result<getLocalFiniteElement(T)>
    {
        typedef const typename std::remove_pointer<T>::type::Tree::FiniteElement& type;
    };

    template<class T>
    typename result<getLocalFiniteElement(T)>::type operator()(const T& t) const
    {
        return t->tree().finiteElement();
    };
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
