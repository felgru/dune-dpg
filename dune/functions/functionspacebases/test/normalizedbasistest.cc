// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <cmath>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/dpg/innerproductfactory.hh>
#include <dune/dpg/spacetuple.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/normalizedbasisadaptor.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include "unrefinedbasistest.hh"

using namespace Dune;
using namespace Dune::Functions;

template<class Basis>
void testNormedAdaptorOn(const typename Basis::GridView gridView)
{
  auto testSpaces = make_space_tuple<Basis>(gridView);

  const FieldVector<double, 2> beta{1./std::sqrt(2.), 1./std::sqrt(2.)};

  auto innerProduct
    = innerProductWithSpace(testSpaces)
      .template addIntegralTerm<0,0,IntegrationType::gradGrad,
                                    DomainOfIntegration::interior>(1., beta)
      .template addIntegralTerm<0,0,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(1.)
      .create();
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
