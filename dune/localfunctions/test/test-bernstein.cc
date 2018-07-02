// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>

#include "../bernstein/pk2d.hh"
#include <dune/localfunctions/test/test-localfe.hh>

template<class FE>
bool testPartitionOfUnity(const FE& fe, unsigned int order = 2)
{
  using LB = typename FE::Traits::LocalBasisType;
  using RangeType = typename LB::Traits::RangeType;
  const auto eps = 1e-10;

  bool success = true;

  // A set of test points
  const Dune::QuadratureRule<double,LB::Traits::dimDomain> quad =
    Dune::QuadratureRules<double,LB::Traits::dimDomain>::rule(fe.type(),order);

  std::vector<RangeType> evaluations;
  for(const auto& qp : quad) {
    const auto& testPoint = qp.position();
    fe.localBasis().evaluateFunction(testPoint, evaluations);

    const auto sumOfBasisFunctionsAtTestPoint
      = std::accumulate(evaluations.begin(), evaluations.end(), RangeType{0});
    if(std::fabs(sumOfBasisFunctionsAtTestPoint - 1) > eps) {
      std::cout << "Test for partition of unity for " << Dune::className(fe)
                << " failed at position " << testPoint
                << " with sum " << sumOfBasisFunctionsAtTestPoint << std::endl;
      std::cout << "Basis function evaluations were: \n";
      for(auto& e : evaluations) {
        std::cout << e << '\n';
      }
      success = false;
    }
  }

  if(success) {
    std::cout << "Test for partition of unity succeeded.\n";
  }
  return success;
}

template<class FE>
bool testJacobianOfPartitionOfUnity(const FE& fe, unsigned int order = 2)
{
  using LB = typename FE::Traits::LocalBasisType;
  using JacobianType = typename LB::Traits::JacobianType;
  using RangeFieldType = typename LB::Traits::RangeFieldType;
  const auto eps = std::numeric_limits<RangeFieldType>::epsilon();

  bool success = true;

  // A set of test points
  const Dune::QuadratureRule<double,LB::Traits::dimDomain> quad =
    Dune::QuadratureRules<double,LB::Traits::dimDomain>::rule(fe.type(),order);

  std::vector<JacobianType> evaluations;
  for(const auto& qp : quad) {
    const auto& testPoint = qp.position();
    fe.localBasis().evaluateJacobian(testPoint, evaluations);

    const auto sumOfBasisFunctionsAtTestPoint
      = std::accumulate(evaluations.begin(), evaluations.end(),
                        JacobianType{0},
                        // why does JacobianType not implement operator+ ?
                        [](const JacobianType& a, const JacobianType& b) {
                          JacobianType sum = a;
                          sum += b;
                          return sum;
                        });
    JacobianType difference{0};
    difference -= sumOfBasisFunctionsAtTestPoint;
    for (int direction=0; direction < LB::Traits::dimDomain; direction++) {
      for(int component=0; component < LB::Traits::dimRange; ++component) {
        if(std::fabs(difference[component][direction]) > eps) {
          std::cout << "Test for Jacobian of partition of unity for "
                    << Dune::className(fe)
                    << " failed at position " << testPoint
                    << ", component " << component
                    << ", direction " << direction
                    << " with sum " << sumOfBasisFunctionsAtTestPoint << '\n';
          std::cout << "Basis function evaluations were: \n";
          for(auto& e : evaluations) {
            std::cout << e << '\n';
          }
          success = false;
        }
      }
    }
  }

  if(success) {
    std::cout << "Test for Jacobian of partition of unity succeeded.\n";
  }
  return success;
}

template<unsigned int order>
bool testBernstein2D()
{
  bool success = true;

  std::cout << "degree: " << order << std::endl;
  Dune::BernsteinPk2DLocalFiniteElement<double,double,order>
      bernsteinTriangle;
  success &= testPartitionOfUnity(bernsteinTriangle);
  success &= testJacobianOfPartitionOfUnity(bernsteinTriangle);
  TEST_FE(bernsteinTriangle);
  return success;
}

int main() try
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
catch (Dune::Exception& e)
{
  std::cout << e << std::endl;
  return 1;
}
