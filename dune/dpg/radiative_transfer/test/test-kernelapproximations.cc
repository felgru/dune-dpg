// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include "../svdkernelapproximation.hh"
#include "../waveletkernelapproximation.hh"
#include "../henyey_greenstein_scattering.hh"

bool sameVector(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
  const double eps = 1e-10;
  for(size_t i=0, imax=x.size(); i < imax; i++) {
    if(fabs(x[i]-y[i]) > eps) {
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
  bool success = true;

  std::cout << "Testing if SVD and wavelet approximation of"
            << " HenyeyGreenstein kernel are identical" << std::endl;
  auto kernel = Dune::HenyeyGreensteinScattering(0.9);
  const double accuracy = 0.01;

  using namespace Dune::ScatteringKernelApproximation;
  HaarWavelet::MatrixCompression
      matrixCompressionApproximation(kernel, accuracy);
  HaarWavelet::SVD waveletSVDApproximation(kernel, accuracy);
  const size_t num_s = waveletSVDApproximation.maxNumS();
  SVD svdApproximation(kernel, num_s);

  Eigen::VectorXd x(num_s);
  for(size_t i=0; i < num_s; i++) x[i] = 1.;
  Eigen::VectorXd y = x;
  Eigen::VectorXd z = x;

  waveletSVDApproximation.applyToVector(z);
  matrixCompressionApproximation.applyToVector(y);
  svdApproximation.applyToVector(x);
  success &= sameVector(x, y);
  success &= sameVector(x, z);

  for(size_t i=0; i < num_s; i++) x[i] = (double)(i*i);
  y = x;
  z = x;

  waveletSVDApproximation.applyToVector(z);
  matrixCompressionApproximation.applyToVector(y);
  svdApproximation.applyToVector(x);
  success &= sameVector(x, y);
  success &= sameVector(x, z);

  return success ? 0 : 1;
}
