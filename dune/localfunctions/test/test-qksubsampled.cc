// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include "../lagrange/qksubsampled.hh"
#include <dune/localfunctions/test/test-localfe.hh>

int main(int argc, char** argv) try
{
  bool success = true;

  const int dim = 2;

  std::cout << "Testing QkDGLocalFiniteElement on 2d"
            << " quadrilateral elements with double precision" << std::endl;
  {
    const unsigned int order=2;
    const unsigned int samples=1;
    std::cout << "order : " << order
              << ", samples : " << samples << std::endl;
    Dune::QkSubsampledLocalFiniteElement<double,double,dim,samples,order>
        lagrangeSubsampledQuad;
    TEST_FE2(lagrangeSubsampledQuad,
            // Disable all tests that assume a standard Lagrange FE
            DisableEvaluate);
  }
  {
    const unsigned int order=3;
    const unsigned int samples=3;
    std::cout << "order : " << order
              << ", samples : " << samples << std::endl;
    Dune::QkSubsampledLocalFiniteElement<double,double,dim,samples,order>
        lagrangeSubsampledQuad;
    //TEST_FE(lagrangeSubsampledQuad);
    TEST_FE2(lagrangeSubsampledQuad,
            DisableEvaluate);
  }

  return success ? 0 : 1;
}
catch (Dune::Exception e)
{
  std::cout << e << std::endl;
  return 1;
}
