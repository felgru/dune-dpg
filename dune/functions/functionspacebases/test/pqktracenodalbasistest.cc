// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/pqktracenodalbasis.hh>

#include "unrefinedbasistest.hh"

using namespace Dune;
using namespace Dune::Functions;

int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);

  // Generate grid for testing
  constexpr int dim = 2;
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
