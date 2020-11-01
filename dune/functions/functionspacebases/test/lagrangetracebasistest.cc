// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/lagrangetracebasis.hh>

#include "unrefinedbasistest.hh"

#include <dune/common/version.hh>

using namespace Dune;
using namespace Dune::Functions;

int main(int argc, char* argv[])
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

  TestSuite t;

  LagrangeTraceBasis<GridView, 1> pq1Basis(gridView);
  t.subTest(testScalarBasis(pq1Basis));

#if DUNE_VERSION_GTE(DUNE_LOCALFUNCTIONS,2,8)
  LagrangeTraceBasis<GridView, 2> pq2Basis(gridView);
  t.subTest(testScalarBasis(pq2Basis, DisableRepresentConstants));

  LagrangeTraceBasis<GridView, 3> pq3Basis(gridView);
  t.subTest(testScalarBasis(pq3Basis, DisableRepresentConstants));

  LagrangeTraceBasis<GridView, 4> pq4Basis(gridView);
  t.subTest(testScalarBasis(pq4Basis, DisableRepresentConstants));
#else
  LagrangeTraceBasis<GridView, 2> pq2Basis(gridView);
  t.subTest(testScalarBasis(pq2Basis));

  LagrangeTraceBasis<GridView, 3> pq3Basis(gridView);
  t.subTest(testScalarBasis(pq3Basis));

  LagrangeTraceBasis<GridView, 4> pq4Basis(gridView);
  t.subTest(testScalarBasis(pq4Basis));
#endif

  return t.exit();
}
