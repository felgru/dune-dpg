// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include "../lagrange/qktrace.hh"
#include <dune/localfunctions/test/test-localfe.hh>

#include <dune/common/version.hh>

int main() try
{
  bool success = true;

  constexpr int dim = 3;

  std::cout << "Testing QkTraceLocalFiniteElement on 3d"
            << " hexagonal elements with double precision" << std::endl;
  {
    const unsigned int order=1;
    std::cout << "order : " << order << std::endl;
    Dune::QkTraceLocalFiniteElement<double,double,dim,order> lagrangeFaceQuad;
    success = success && testFE(lagrangeFaceQuad);
  }
  {
    const unsigned int order=2;
    std::cout << "order : " << order << std::endl;
    Dune::QkTraceLocalFiniteElement<double,double,dim,order> lagrangeFaceQuad;
    success = success && testFE(lagrangeFaceQuad,
#if DUNE_VERSION_GTE(DUNE_LOCALFUNCTIONS,2,8)
                                DisableRepresentConstants
#endif
                                );
  }

  return success ? 0 : 1;
}
catch (Dune::Exception& e)
{
  std::cout << e << std::endl;
  return 1;
}
