// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/localfunctions/test/test-localfe.hh>

#include <dune/functions/functionspacebases/pqktracenodalbasis.hh>

using namespace Dune;
using namespace Dune::Functions;

template <typename Basis>
void testScalarBasis(const Basis& feBasis)
{
  /////////////////////////////////////////////////////////////////////
  //  Run the dune-localfunctions test for the LocalFiniteElement of
  //  each grid element
  /////////////////////////////////////////////////////////////////////

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


  // Check whether the basis exports a type 'MultiIndex'
  typedef typename Basis::MultiIndex MultiIndex;

  // And this type must be indexable
  static_assert(is_indexable<MultiIndex>(),
      "MultiIndex must support operator[]");

  /////////////////////////////////////////////////////////////////////
  //  Check whether the global indices are in the correct range,
  //  and whether each global index appears at least once.
  /////////////////////////////////////////////////////////////////////

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
      DUNE_THROW(Exception, "Index [" << i << "] does not exist as global basis vector");

  // Sample the function f(x,y) = x on the grid vertices
  // If we use that as the coefficients of a finite element function,
  // we know its integral and can check whether quadrature returns
  // the correct result
  std::vector<double> x(feBasis.size());

  // TODO: Implement interpolation properly using the global basis.
  const int dim = Basis::GridView::dimension;
  for (const auto& element : elements(gridView))
    x[gridView.indexSet().index(element)] = element.geometry().corner(0)[0];

  // Objects required in the local context
  auto localIndexSet2 = feBasis.localIndexSet();
  std::vector<double> coefficients(localView.maxSize());

  // Loop over elements and integrate over the function
  double integral = 0;
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

    // copy data from global vector
    coefficients.resize(localIndexSet.size());
    for (size_t i=0; i<localIndexSet.size(); i++)
    {
      coefficients[i] = x[localIndexSet.index(i)[0]];
    }

    typedef typename Basis::LocalView::Tree Tree;
    const Tree& tree = localView.tree();

    auto& localFiniteElement = tree.finiteElement();

    // we have a flat tree...
    assert(localView.size() == tree.size());
    assert(localView.size() == tree.finiteElement().localBasis().size());

    const QuadratureRule<double, dim>& quad
        = QuadratureRules<double, dim>::rule(element.type(), 1);

    for ( size_t pt=0; pt < quad.size(); pt++ ) {
      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();

      // The multiplicative factor in the integral transformation formula
      const double integrationElement
          = element.geometry().integrationElement(quadPos);

      std::vector<FieldVector<double,1> > shapeFunctionValues;
      localFiniteElement.localBasis().evaluateFunction(quadPos, shapeFunctionValues);

      // Actually compute the vector entries
      for (size_t i=0; i<localFiniteElement.localBasis().size(); i++)
      {
        integral += coefficients[tree.localIndex(i)] * shapeFunctionValues[i]
                    * quad[pt].weight() * integrationElement;
      }
    }

    localIndexSet.unbind();
    localView.unbind();
  }

  std::cout << "Computed integral is " << integral << std::endl;
  if (std::abs(integral-0.5) > 1e-10)
    std::cerr << "Warning: integral value is wrong!" << std::endl;
}

int main (int argc, char* argv[]) try
{
  // Generate grid for testing
  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  const FieldVector<double,dim> l(1);
  const std::array<int,dim> elements = {{10, 10}};
  const GridType grid(l,elements);

  typedef GridType::LeafGridView GridView;
  const GridView gridView = grid.leafGridView();

  PQkTraceNodalBasis<GridView, 1> pq1Basis(gridView);
  testScalarBasis(pq1Basis);

  PQkTraceNodalBasis<GridView, 2> pq2Basis(gridView);
  testScalarBasis(pq2Basis);

  PQkTraceNodalBasis<GridView, 3> pq3Basis(gridView);
  testScalarBasis(pq3Basis);

  PQkTraceNodalBasis<GridView, 4> pq4Basis(gridView);
  testScalarBasis(pq4Basis);

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
