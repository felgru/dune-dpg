#include <boost/math/constants/constants.hpp>

#include <dune/common/fvector.hh>
#include <dune/dpg/radiative_transfer/kanschat_scattering.hh>

#include <Eigen/Core>
#include <Eigen/SVD>

#include <cmath>
#include <iostream>
#include <vector>

std::vector<double> exponentialDecay(size_t n) {
  std::vector<double> a(n);
  for(std::ptrdiff_t i=0; (size_t)i<n; ++i) {
    a[i] = exp(-i);
  }
  return a;
}

std::vector<double> polynomialDecay(size_t n) {
  std::vector<double> a(n);
  for(size_t i=0; i<n; ++i) {
    a[i] = 1./(i+1);
  }
  return a;
}

template<class Direction, class Kernel>
Eigen::MatrixXd kernelMatrix(const std::vector<Direction>& sVector, Kernel&& kernel) {
  size_t numS = sVector.size();
  Eigen::MatrixXd kernelMatrix(numS, numS);
  for(size_t j = 0; j < numS; ++j) {
    Direction s_j = sVector[j];
    for(size_t i = 0; i < numS; ++i) {
      Direction s_i = sVector[i];
      // TODO: maybe use a higher order quadrature
      kernelMatrix(i,j) = kernel(s_i, s_j)/(numS*numS);
    }
  }
  return kernelMatrix;
}

int main() {
  using namespace Dune;

  const unsigned int dim = 2;
  using Direction = FieldVector<double, dim>;

  const size_t numS = 100;
  std::vector<Direction> sVector(numS);
  for(size_t i = 0; i < numS; ++i)
  {
    using namespace boost::math::constants;
    sVector[i] = {cos(2*pi<double>()*i/numS),
                  sin(2*pi<double>()*i/numS)};
  }

  auto exponentialKernel = KanschatScattering<Direction>(exponentialDecay(50));
  auto polynomialKernel = KanschatScattering<Direction>(polynomialDecay(50));
  auto acpKernel = ACP2011Scattering<Direction>();
  using namespace Eigen;
  MatrixXd m = kernelMatrix(sVector, acpKernel);

  Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::NoQRPreconditioner>
      kernelSVD(numS, numS, Eigen::ComputeThinU | Eigen::ComputeThinV);
  kernelSVD.compute(m);

  VectorXd singularValues = kernelSVD.singularValues();
  for(size_t i=0, size=singularValues.size(); i<size; ++i)
    std::cout << singularValues[i] << '\n';
}
