#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <boost/math/constants/constants.hpp>

#include <dune/common/fvector.hh>
#include <dune/dpg/radiative_transfer/kanschat_scattering.hh>
#include <dune/dpg/radiative_transfer/henyey_greenstein_scattering.hh>
#include <dune/dpg/radiative_transfer/waveletkernelapproximation.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>

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
    a[i] = std::exp(-i);
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

std::vector<double> angles(unsigned int wltOrder,
                           unsigned int level)
{
  using namespace Dune;

  // Get Gauss-Legendre quadrature in [0,1]
  const int quadOrder = 2*wltOrder+1;
  const size_t dim = 1;
  const QuadratureRule<double, dim>& quad =
    QuadratureRules<double, dim>::rule(GeometryTypes::line,
      quadOrder, QuadratureType::GaussLegendre);
  using namespace boost::math::constants;
  std::vector<double> angles;
  angles.reserve((1<<level)*quad.size());
  for(int k = 0; k < (1<<level); ++k)
  {
    for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {
      const double angle = 2*pi<double>()*(k+quad[pt].position())
                           / (1<<level) - pi<double>(); // Angle in [-pi,pi]
      angles.push_back(angle);
    }
  }
  return angles;
}

#if 0
template<class Kernel>
Eigen::MatrixXd kernelMatrix(const std::vector<double>& angles,
                             Kernel&& kernel)
{
  size_t numS = angles.size();
  Eigen::MatrixXd kernelMatrix(numS, numS);
  for(size_t j = 0; j < numS; ++j) {
    const double s_j = angles[j];
    for(size_t i = 0; i < numS; ++i) {
      const double s_i = angles[i];
      // TODO: maybe use a higher order quadrature
      kernelMatrix(i,j) = kernel(s_i - s_j);
    }
  }
  return kernelMatrix;
}
#endif

void printHelp(char* name) {
  std::cerr << "Usage: " << name << " [-e n|-p n|-a] [-w] l\n"
               " -e n: exponential decay with n terms\n"
               " -p n: polynomial decay with n terms\n"
               " -g gamma: Henyey-Greenstein kernel with given gamma\n"
               " -a: kernel from Avila et al. 2011\n"
               " -w: use Alpert wavelet representation of kernel matrix\n"
               " l: maximal wavelet level\n";
  exit(0);
}

int main(int argc, char *argv[]) {
  enum class kernelType {
    exponential,
    polynomial,
    henyey_greenstein,
    ACP2011
  };

  kernelType kt = kernelType::exponential;
  size_t kernelTerms = 50;
  double gamma = 0.;
  size_t maxLevel = 5;
  bool wavelet = false;

  int opt;
  while ((opt = getopt(argc,argv,"e:p:g:aw")) != EOF)
    switch(opt)
    {
      case 'e': kt = kernelType::exponential; kernelTerms = atoi(optarg); break;
      case 'p': kt = kernelType::polynomial; kernelTerms = atoi(optarg); break;
      case 'g': kt = kernelType::henyey_greenstein; gamma = atof(optarg); break;
      case 'a': kt = kernelType::ACP2011; break;
      case 'w': wavelet = true; break;
      default:
      case '?':
        printHelp(argv[0]);
    }
  if(optind != argc-1) {
    printHelp(argv[0]);
  }
  if(!wavelet) {
    std::cerr << "computation of kernel SVDs is currently only implemented"
                 " for wavelet representation of kernels.\n";
    exit(0);
  }
  maxLevel = atoi(argv[optind]);
  const size_t wltOrder = 3;

  using namespace Dune;

  using namespace Dune::ScatteringKernelApproximation::AlpertWavelet;

  Eigen::MatrixXd m;
  switch(kt) {
    case kernelType::exponential:
      {
        auto exponentialKernel
          = KanschatScattering(exponentialDecay(kernelTerms));
        std::cout << "Singular values of " << exponentialKernel.info() << '\n';
        if(wavelet) {
          m = waveletKernelMatrix(exponentialKernel, wltOrder, maxLevel, 16);
#if 0
        } else {
          m = kernelMatrix(angles(wltOrder, maxLevel), exponentialKernel);
#endif
        }
      }
      break;
    case kernelType::polynomial:
      {
        auto polynomialKernel
          = KanschatScattering(polynomialDecay(kernelTerms));
        std::cout << "Singular values of " << polynomialKernel.info() << '\n';
        if(wavelet) {
          m = waveletKernelMatrix(polynomialKernel, wltOrder, maxLevel, 16);
#if 0
        } else {
          m = kernelMatrix(angles(wltOrder, maxLevel), polynomialKernel);
#endif
        }
      }
      break;
    case kernelType::henyey_greenstein:
      {
        auto henyeyGreensteinKernel = HenyeyGreensteinScattering(gamma);
        std::cout << "Singular values of "
                  << henyeyGreensteinKernel.info() << '\n';
        if(wavelet) {
          m = waveletKernelMatrix(henyeyGreensteinKernel, wltOrder,
                                  maxLevel, 16);
#if 0
        } else {
          m = kernelMatrix(angles(wltOrder, maxLevel), henyeyGreensteinKernel);
#endif
        }
      }
      break;
    case kernelType::ACP2011:
      {
        auto acpKernel = ACP2011Scattering();
        std::cout << "Singular values of " << acpKernel.info() << '\n';
        if(wavelet) {
          m = waveletKernelMatrix(acpKernel, wltOrder, maxLevel, 16);
#if 0
        } else {
          m = kernelMatrix(angles(wltOrder, maxLevel), acpKernel);
#endif
        }
      }
  }

  const size_t numS = m.rows();
  Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::NoQRPreconditioner>
      kernelSVD(numS, numS, Eigen::ComputeThinU | Eigen::ComputeThinV);
  kernelSVD.compute(m);

  std::cout << kernelSVD.singularValues() << '\n';
}
