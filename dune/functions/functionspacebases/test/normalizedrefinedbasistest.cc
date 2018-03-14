// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <cmath>
#include <iostream>

#include <dune/common/exceptions.hh>

#include <dune/dpg/assemble_helper.hh>
#include <dune/dpg/assemble_types.hh>
#include <dune/dpg/innerproduct.hh>
#include <dune/dpg/integralterm.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/normalizedrefinedbasisadaptor.hh>
#include <dune/functions/functionspacebases/pqkdgrefineddgnodalbasis.hh>

#include "refinedbasistest.hh"

using namespace Dune;
using namespace Dune::Functions;

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

  NormalizedRefinedBasis<InnerProduct> normedBasis(innerProduct);
  testRefinedScalarBasis(normedBasis);
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

  typedef Functions::PQkDGRefinedDGBasis<GridView, 1, 3> FEBasisTest;
  //typedef Functions::LagrangeDGBasis<GridView, 3> FEBasisTest;
  auto testSpaces = make_space_tuple<FEBasisTest>(gridView);

  testNormedAdaptorOn<PQkDGRefinedDGBasis<GridView, 1, 2>>(gridView);

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
