// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/pqkdgrefineddgnodalbasis.hh>

#include "refinedbasistest.hh"

using namespace Dune;
using namespace Dune::Functions;

int main () try
{
  // Generate grid for testing
  constexpr int dim = 2;
  typedef YaspGrid<dim> GridType;
  const FieldVector<double,dim> l(1);
  const std::array<int,dim> elements = {{10, 10}};
  const GridType grid(l,elements);

  typedef GridType::LeafGridView GridView;
  const GridView gridView = grid.leafGridView();

  Functions::PQkDGRefinedDGBasis<GridView, 1, 2> dgrefinedbasis(gridView);
  testRefinedScalarBasis(dgrefinedbasis);

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
