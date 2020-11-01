// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include "../lagrange/pkdgsubsampled2d.hh"
#include <dune/localfunctions/test/test-localfe.hh>

template<unsigned int order, unsigned int samples>
bool testDGSubsampled2D()
{
  bool success = true;

  std::cout << "order : " << order
            << ", samples : " << samples << std::endl;
  Dune::PkDGSubsampled2DLocalFiniteElement<double,double,samples,order>
      lagrangeSubsampledTria;
  TEST_FE2(lagrangeSubsampledTria,
           // Disable all tests that assume a standard Lagrange FE
           DisableLocalInterpolation
           | DisableEvaluate
           | DisableRepresentConstants
           );
  return success;
}

int main() try
{
  bool success = true;

  std::cout << "Testing PkDGSubsampledLocalFiniteElement on 2d"
            << " triangular elements with double precision" << std::endl;

  success &= testDGSubsampled2D<2, 1>();
  success &= testDGSubsampled2D<1, 2>();
  success &= testDGSubsampled2D<2, 2>();

  return success ? 0 : 1;
}
catch (Dune::Exception& e)
{
  std::cout << e << std::endl;
  return 1;
}
