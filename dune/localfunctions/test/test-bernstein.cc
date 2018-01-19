// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include "../bernstein/pk2d.hh"
#include <dune/localfunctions/test/test-localfe.hh>

template<unsigned int order>
bool testBernstein2D()
{
  bool success = true;

  std::cout << "order : " << order << std::endl;
  Dune::BernsteinPk2DLocalFiniteElement<double,double,order>
      bernsteinTriangle;
  TEST_FE(bernsteinTriangle);
  return success;
}

int main(int argc, char** argv) try
{
  bool success = true;

  std::cout << "Testing BernsteinPk2DLocalFiniteElement on 2d"
            << " triangular elements with double precision" << std::endl;

  success &= testBernstein2D<0>();
  success &= testBernstein2D<1>();
  success &= testBernstein2D<2>();
  success &= testBernstein2D<3>();

  return success ? 0 : 1;
}
catch (Dune::Exception e)
{
  std::cout << e << std::endl;
  return 1;
}
