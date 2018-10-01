// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_UNREFINEDBASISTEST_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_UNREFINEDBASISTEST_HH

#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/version.hh>

#include <dune/localfunctions/test/test-localfe.hh>

namespace Dune {

template <typename Basis>
void testLocalFeForEachElement(const Basis& feBasis)
{
  typedef typename Basis::GridView GridView;
  const GridView gridView = feBasis.gridView();

  typename Basis::LocalView localView(feBasis);


  // Test the LocalFiniteElement
  for (const auto& element : elements(gridView))
  {
    localView.bind(element);

    // The general LocalFiniteElement unit test from
    // dune/localfunctions/test/test-localfe.hh
    const auto& lFE = localView.tree().finiteElement();
    testFE(lFE);

    if (lFE.size() != localView.size())
      DUNE_THROW(Exception,
          "Size of leaf node and finite element do not coincide");
  }
}

template <typename Basis>
void checkRangeOfGlobalIndices(const Basis& feBasis)
{
  /////////////////////////////////////////////////////////////////////
  //  Check whether the global indices are in the correct range,
  //  and whether each global index appears at least once.
  /////////////////////////////////////////////////////////////////////

  const auto gridView = feBasis.gridView();

  auto localView = feBasis.localView();

  std::vector<bool> seen(feBasis.size(), false);

  // Loop over all leaf elements
  for (const auto& element : elements(gridView))
  {
    localView.bind(element);

    for (size_t i=0; i<localView.tree().size(); i++)
    {
      if (localView.index(i)[0] < 0)
        DUNE_THROW(Exception, "Index is negative, which is not allowed");

      if (localView.index(i)[0] >= seen.size())
        DUNE_THROW(Exception, "Local index " << i
                           << " is mapped to global index "
                           << localView.index(i)
                           << ", which is larger than allowed");

      seen[localView.index(i)[0]] = true;
    }
  }

  for (size_t i=0; i<seen.size(); i++)
    if (! seen[i])
      DUNE_THROW(Exception,
          "Index [" << i << "] does not exist as global basis vector");
}

template <typename Basis>
void checkConsistencyOfLocalViewAndIndexSet(const Basis& feBasis)
{
  const auto gridView = feBasis.gridView();
  typename Basis::LocalView localView(feBasis);
  auto localView2 = feBasis.localView();

  for (const auto& element : elements(gridView))
  {
    localView.bind(element);
    localView2.bind(element);

    // paranoia checks
    assert(&(localView.globalBasis()) == &(feBasis));

    assert(localView.size() == localView2.size());
    for (size_t i=0; i<localView.size(); i++)
      assert(localView.index(i) == localView2.index(i));

    typedef typename Basis::LocalView::Tree Tree;
    const Tree& tree = localView.tree();

    // we have a flat tree...
    assert(localView.size() == tree.size());
    assert(localView.size() == tree.finiteElement().localBasis().size());

    localView.unbind();
  }
}

template <typename Basis>
void testScalarBasis(const Basis& feBasis)
{
  testLocalFeForEachElement(feBasis);

  // Check whether the basis exports a type 'MultiIndex'
  typedef typename Basis::MultiIndex MultiIndex;

  // And this type must be indexable
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,7)
  static_assert(IsIndexable<MultiIndex>(),
      "MultiIndex must support operator[]");
#else
  static_assert(is_indexable<MultiIndex>(),
      "MultiIndex must support operator[]");
#endif

  checkRangeOfGlobalIndices(feBasis);
  checkConsistencyOfLocalViewAndIndexSet(feBasis);
}

}
#endif
