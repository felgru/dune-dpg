// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <cmath>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>

#include <dune/localfunctions/test/test-localfe.hh>

#include <dune/dpg/assemble_helper.hh>
#include <dune/dpg/assemble_types.hh>
#include <dune/dpg/innerproduct.hh>
#include <dune/dpg/integralterm.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/normalizedbasisadaptor.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

using namespace Dune;
using namespace Dune::Functions;

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
    assert(localView.size() == tree.finiteElement().localBasis().size());

    localIndexSet.unbind();
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
  static_assert(is_indexable<MultiIndex>(),
      "MultiIndex must support operator[]");

  checkRangeOfGlobalIndices(feBasis);
  checkConsistencyOfLocalViewAndIndexSet(feBasis);
}


template<class Basis>
void testNormedAdaptorOn(const typename Basis::GridView gridView)
{
  auto testSpaces = make_space_tuple<Basis>(gridView);

  const FieldVector<double, 2> beta{1./std::sqrt(2.), 1./std::sqrt(2.)};

  auto innerProduct = make_InnerProduct(testSpaces,
          make_tuple(
              make_IntegralTerm<0,0,IntegrationType::gradGrad,
                                    DomainOfIntegration::interior>(1., beta),
              make_IntegralTerm<0,0,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(1.)));
  using InnerProduct = decltype(innerProduct);

  NormalizedBasis<InnerProduct> normedBasis(innerProduct);
  testScalarBasis(normedBasis);
}

int main () try
{
  // Generate grid for testing
  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  const FieldVector<double,dim> l(1);
  const std::array<int,dim> elements = {{10, 10}};
  const GridType grid(l,elements);

  typedef GridType::LeafGridView GridView;
  const GridView gridView = grid.leafGridView();

  testNormedAdaptorOn<LagrangeDGBasis<GridView, 2>>(gridView);

  return 0;

} catch ( Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
catch(...)
{
  std::cerr << "Unknown exception thrown!" << std::endl;
  return 1;
}
