#include <boost/math/constants/constants.hpp>

#include <dune/common/fvector.hh>
#include <dune/dpg/radiative_transfer/kanschat_scattering.hh>

#include <Eigen/Core>
#include <Eigen/SVD>

#include <cmath>
#include <iostream>
#include <vector>

#include <stdlib.h>
#include <unistd.h>

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
    a[i] = 1./((i+1)*(i+1));
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
      kernelMatrix(i,j) = kernel(s_i, s_j)/numS;
    }
  }
  return kernelMatrix;
}

void printHelp(char* name) {
  std::cerr << "Usage: " << name << " [-e n|-p n|-a] s\n"
               " -e n: exponential decay with n terms\n"
               " -p n: polynomial decay with n terms\n"
               " -a: kernel from Avila et al. 2011\n"
               " s: number of directions\n";
  exit(0);
}

int main(int argc, char *argv[]) {
  enum class kernelType {
    exponential,
    polynomial,
    ACP2011
  };

  kernelType kt = kernelType::exponential;
  size_t kernelTerms = 50;
  size_t numS = 100;

  int opt;
  while ((opt = getopt(argc,argv,"e:p:a")) != EOF)
    switch(opt)
    {
      case 'e': kt = kernelType::exponential; kernelTerms = atoi(optarg); break;
      case 'p': kt = kernelType::polynomial; kernelTerms = atoi(optarg); break;
      case 'a': kt = kernelType::ACP2011; break;
      default:
      case '?':
        printHelp(argv[0]);
    }
  if(optind != argc-1) {
    printHelp(argv[0]);
  }
  numS = atoi(argv[optind]);

  using namespace Dune;

  const unsigned int dim = 2;
  using Direction = FieldVector<double, dim>;

  std::vector<Direction> sVector(numS);
  for(size_t i = 0; i < numS; ++i)
  {
    using namespace boost::math::constants;
    sVector[i] = {cos(2*pi<double>()*i/numS),
                  sin(2*pi<double>()*i/numS)};
  }

  Eigen::MatrixXd m;
  switch(kt) {
    case kernelType::exponential:
      {
        auto exponentialKernel
          = KanschatScattering<Direction>(exponentialDecay(kernelTerms));
        m = kernelMatrix(sVector, exponentialKernel);
      }
      break;
    case kernelType::polynomial:
      {
        auto polynomialKernel
          = KanschatScattering<Direction>(polynomialDecay(kernelTerms));
        m = kernelMatrix(sVector, polynomialKernel);
      }
      break;
    case kernelType::ACP2011:
      {
        auto acpKernel = ACP2011Scattering<Direction>();
        m = kernelMatrix(sVector, acpKernel);
      }
  }

  Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::NoQRPreconditioner>
      kernelSVD(numS, numS, Eigen::ComputeThinU | Eigen::ComputeThinV);
  kernelSVD.compute(m);

  Eigen::VectorXd singularValues = kernelSVD.singularValues();
  for(size_t i=0, size=singularValues.size(); i<size; ++i)
    std::cout << singularValues[i] << '\n';
}
