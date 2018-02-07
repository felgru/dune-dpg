// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_REFINEDBASISTEST_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_REFINEDBASISTEST_HH

#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>

#include <dune/localfunctions/test/test-localfe.hh>

namespace Dune {

template <typename Basis>
void testLocalFeForEachElement(const Basis& feBasis)
{
  typedef typename Basis::GridView GridView;
  const GridView gridView = feBasis.gridView();

  typename Basis::LocalView localView(feBasis);


  for (const auto& element : elements(gridView))
  {
    localView.bind(element);

    // The general LocalFiniteElement unit test from
    // dune/localfunctions/test/test-localfe.hh
    const auto& lFE = localView.tree().finiteElement();
    testFE(lFE);

    const size_t numSubElements =
        localView.tree().refinedReferenceElementGridView().size(0);
    if (lFE.size() * numSubElements != localView.size())
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

  typename Basis::LocalView localView(feBasis);

  std::vector<bool> seen(feBasis.size(), false);

  auto localIndexSet = feBasis.localIndexSet();

  // Loop over all leaf elements
  for (const auto& element : elements(gridView))
  {
    localView.bind(element);
    localIndexSet.bind(localView);

    for (size_t i=0; i<localView.tree().size(); i++)
    {
      if (localIndexSet.index(i)[0] < 0)
        DUNE_THROW(Exception, "Index is negative, which is not allowed");

      if (localIndexSet.index(i)[0] >= seen.size())
        DUNE_THROW(Exception, "Local index " << i
                           << " is mapped to global index "
                           << localIndexSet.index(i)
                           << ", which is larger than allowed");

      seen[localIndexSet.index(i)[0]] = true;
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
  auto localIndexSet = feBasis.localIndexSet();

  // Objects required in the local context
  auto localIndexSet2 = feBasis.localIndexSet();

  for (const auto& element : elements(gridView))
  {
    localView.bind(element);
    localIndexSet.bind(localView);
    localIndexSet2.bind(localView);

    // paranoia checks
    assert(localView.size() == localIndexSet.size());
    assert(&(localView.globalBasis()) == &(feBasis));
    assert(&(localIndexSet.localView()) == &(localView));

    assert(localIndexSet.size() == localIndexSet2.size());
    for (size_t i=0; i<localIndexSet.size(); i++)
      assert(localIndexSet.index(i) == localIndexSet2.index(i));

    typedef typename Basis::LocalView::Tree Tree;
    const Tree& tree = localView.tree();

    // we have a flat tree...
    assert(localView.size() == tree.size());

    localIndexSet.unbind();
    localView.unbind();
  }
}

template <typename Basis>
void testRefinedScalarBasis(const Basis& feBasis)
{
  testLocalFeForEachElement(feBasis);

  // Check whether the basis exports a type 'MultiIndex'
  typedef typename Basis::MultiIndex MultiIndex;

  // And this type must be indexable
  static_assert(is_indexable<MultiIndex>(),
      "MultiIndex must support operator[]");

  checkRangeOfGlobalIndices(feBasis);
  checkConsistencyOfLocalViewAndIndexSet(feBasis);
}

}
#endif
