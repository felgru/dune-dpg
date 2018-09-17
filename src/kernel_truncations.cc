#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <dune/dpg/radiative_transfer/kanschat_scattering.hh>
#include <dune/dpg/radiative_transfer/henyey_greenstein_scattering.hh>
#include <dune/dpg/radiative_transfer/waveletkernelapproximation.hh>

#include <Eigen/Core>

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

template<size_t wltOrder, class Kernel>
std::vector<double>
computeTruncationErrors(const Kernel& kernel,
                        size_t minLevel, size_t maxLevel)
{
  std::cout << "truncation errors of " << kernel.info() << '\n';
  using namespace Dune::ScatteringKernelApproximation::AlpertWavelet;
  return SVD<wltOrder>(kernel, minLevel, maxLevel).getTruncationErrors();
}

void printHelp(char* name) {
  std::cerr << "Print truncation errors of Alpert wavelet truth matrix\n"
               "Usage: " << name << " [-e n|-p n|-a] [-w] l\n"
               " -e n: exponential decay with n terms\n"
               " -p n: polynomial decay with n terms\n"
               " -g gamma: Henyey-Greenstein kernel with given gamma\n"
               " -a: kernel from Avila et al. 2011\n"
               " l: wavelet level of truth matrix\n";
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

  int opt;
  while ((opt = getopt(argc,argv,"e:p:g:aw")) != EOF)
    switch(opt)
    {
      case 'e': kt = kernelType::exponential; kernelTerms = atoi(optarg); break;
      case 'p': kt = kernelType::polynomial; kernelTerms = atoi(optarg); break;
      case 'g': kt = kernelType::henyey_greenstein; gamma = atof(optarg); break;
      case 'a': kt = kernelType::ACP2011; break;
      default:
      case '?':
        printHelp(argv[0]);
    }
  if(optind != argc-1) {
    printHelp(argv[0]);
  }
  maxLevel = atoi(argv[optind]);
  const size_t minLevel = 2;
  const size_t wltOrder = 3;

  using namespace Dune;

  using namespace Dune::ScatteringKernelApproximation::AlpertWavelet;

  std::vector<double> truncationErrors;
  switch(kt) {
    case kernelType::exponential:
      {
        auto exponentialKernel
          = KanschatScattering(exponentialDecay(kernelTerms));
        truncationErrors = computeTruncationErrors<wltOrder>
                           (exponentialKernel, minLevel, maxLevel);
      }
      break;
    case kernelType::polynomial:
      {
        auto polynomialKernel
          = KanschatScattering(polynomialDecay(kernelTerms));
        truncationErrors = computeTruncationErrors<wltOrder>
                           (polynomialKernel, minLevel, maxLevel);
      }
      break;
    case kernelType::henyey_greenstein:
      {
        auto henyeyGreensteinKernel = HenyeyGreensteinScattering(gamma);
        truncationErrors = computeTruncationErrors<wltOrder>
                           (henyeyGreensteinKernel, minLevel, maxLevel);
      }
      break;
    case kernelType::ACP2011:
      {
        auto acpKernel = ACP2011Scattering();
        truncationErrors = computeTruncationErrors<wltOrder>
                           (acpKernel, minLevel, maxLevel);
      }
  }

  for(auto&& error : truncationErrors) {
    std::cout << error << '\n';
  }
}
