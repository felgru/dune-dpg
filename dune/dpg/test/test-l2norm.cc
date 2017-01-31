// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cmath>
#include <iostream>
#include <dune/grid/yaspgrid.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include "../errortools.hh"

int main() try {
  std::cout << "Testing L2 norm of constant FE function on unform 2D grid."
            << std::endl;

  Dune::YaspGrid<2> grid(Dune::FieldVector<double, 2>(1.),
                         std::array<int, 2>{2, 2});
  auto gridView = grid.leafGridView();
  using GridView = decltype(gridView);
  Dune::Functions::PQkNodalBasis<GridView, 2> feBasis(gridView);

  Dune::BlockVector<Dune::FieldVector<double,1> >
      feCoefficients(feBasis.size());
  feCoefficients = 1.;

  const double l2norm = Dune::ErrorTools::l2norm(feBasis, feCoefficients);
  const double error = 1. - l2norm;

  if(std::fabs(error) > 1e-6) {
    std::cout << "L2 norm is " << l2norm << ", but should be 1!" << std::endl;
    return 1;
  } else {
    return 0;
  }
} catch (Dune::Exception e) {
  std::cout << e << std::endl;
  return 1;
}
