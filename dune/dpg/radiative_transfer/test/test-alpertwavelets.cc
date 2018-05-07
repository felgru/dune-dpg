// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cmath>
#include <iostream>

#include <boost/math/constants/constants.hpp>

#include <Eigen/Core>

#include "../waveletkernelapproximation.hh"

bool sameVector(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
  const double eps = 1e-10;
  for(size_t i=0, imax=x.size(); i < imax; i++) {
    if(std::fabs(x[i]-y[i]) > eps) {
      std::cout << "Comparison failed at index " << i
                << ", value should be " << x[i]
                << " but was " << y[i] << '\n';
      return false;
    }
  }
  return true;
}

int main()
{
  using namespace Dune::ScatteringKernelApproximation;
  bool success = true;

  //****** Alpert wavelets ***********************
  // Number of vanishing moments (polynomials of degree <=L-1 cancel out)
  const size_t L = 1;
  // Starting level of the space V_J
  const size_t J = 1;
  // Gauss-Legendre quadrature order to use in internal routines
  // Remark: quadOrder = 2*nQuad -1 in Gauss-Legendre
  const size_t quadOrder = 20;
  // Interval [-r,r]
  const double r = boost::math::constants::pi<double>();
  // Input data
  // Given a function f in [-r,r], we project it to the space V_J of
  // piecewise polynomials of degree L-1 in each interval I_{j,k}.
  // This will be our input data for DWT.
  auto flambda = [](auto a,auto xmin,auto xmax) { return a; };
  Eigen::VectorXd data = AlpertWavelet::ProjectOntoVJ(flambda,r,J,L,quadOrder);
  // // Other trivial input data
  // // First option
  // Eigen::VectorXd data(static_cast<int>(L*pow(2,J)));
  // data << 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16;
  // // Second option
  // Eigen::VectorXd data(16);
  // for(size_t i=0, imax=data.size(); i < imax; i++) data(i) = 1.;

  // Compute DWT
  auto y = AlpertWavelet::DWT(data,L,J,quadOrder);
  // Compute IDWT
  Eigen::VectorXd dataIDWT = AlpertWavelet::IDWT(y,L,J,quadOrder);
  success = sameVector(dataIDWT,data);
  if(success)
    std::cout << "IDWT(DWT(data)) successful." << std::endl;
  else
    std::cout << "IDWT(DWT(data)) not successful." << std::endl;


  // //****** Felix *********************
  // std::cout << "testing DWT and inverse DWT" << std::endl;
  // Eigen::VectorXd x(16);
  // for(size_t i=0, imax=x.size(); i < imax; i++) x[i] = 1.;
  // Eigen::VectorXd y = x;
  // HaarWavelet::DWT(y);
  // Eigen::VectorXd wlt(x.size());
  // for(size_t i=0, imax=wlt.size(); i < imax; i++) wlt[i] = 0;
  // wlt(0) = std::sqrt(2*boost::math::constants::pi<double>());
  // std::cout << "Wavelet transform of constant 1 function\n";
  // success &= sameVector(wlt, y);

  // std::cout << "IDWT(DWT(1)) == 1\n";
  // HaarWavelet::IDWT(y);
  // success &= sameVector(x, y);

  // std::cout << "Wavelet transform of constant 0 function\n";
  // for(size_t i=0, imax=x.size(); i < imax; i++) x[i] = 0.;
  // y = x;
  // HaarWavelet::DWT(y);
  // success &= sameVector(x, y);

  // std::cout << "IDWT(DWT(0)) == 0\n";
  // HaarWavelet::IDWT(y);
  // success &= sameVector(x, y);

  // std::cout << "IDWT(DWT(x)) == x\n";
  // for(size_t i=0, imax=x.size(); i < imax; i++) x[i] = (double)i;
  // y = x;
  // HaarWavelet::DWT(y); HaarWavelet::IDWT(y);
  // success &= sameVector(x, y);

  // std::cout << "Wavelet transform of f = 1 for x < 0 and f = -1 for x >= 0\n";
  // for(size_t i=0, imax=x.size()/2; i < imax; i++)        x[i] =  1.;
  // for(size_t i=x.size()/2, imax=x.size(); i < imax; i++) x[i] = -1.;
  // y = x;
  // HaarWavelet::DWT(y);
  // for(size_t i=0, imax=wlt.size(); i < imax; i++) wlt[i] = 0;
  // wlt(1) = std::sqrt(boost::math::constants::pi<double>());
  // success &= sameVector(wlt, y);
  // std::cout << "IDWT(DWT(f)) == f\n";
  // HaarWavelet::IDWT(y);
  // success &= sameVector(x, y);

  return success ? 0 : 1;
}
