// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_WAVELET_KERNEL_APPROXIMATION_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_WAVELET_KERNEL_APPROXIMATION_HH

#include <algorithm>
#include <cmath>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <boost/math/constants/constants.hpp>

#include <Eigen/Core>

namespace Dune {

namespace ScatteringKernelApproximation {
  class HaarWavelet {
    private:
      // Scaling function $\phi_{0,0}$ on the interval [-r,r)
      static Eigen::VectorXd sf(const Eigen::VectorXd& x, double r) {
        double normfactor = 1./sqrt(2*r);
        Eigen::VectorXd res(x.size());
        for(size_t i=0, imax=x.size(); i<imax; i++) {
          double val = x(i);
          res(i) = normfactor * (val >= -r and val < r);
        }
        return res;
      }

      // Wavelet $\psi_{n,k}$ on the interval [-r,r)
      // Its support is in [xmin,xmax] where
      // 	xmin = r*k*pow(2,-j+1)-r
      // 	xmax = r*(k+1)*pow(2,-j+1)-r
      static Eigen::VectorXd wlt(size_t j, size_t k,
                                 const Eigen::VectorXd& x, double r) {
        double xmin = r *   k   * std::exp2(-(double)j+1) - r;
        double xmax = r * (k+1) * std::exp2(-(double)j+1) - r;
        double xmiddle = (xmin + xmax)/2;
        double normfactor = 1./sqrt(xmax-xmin);

        Eigen::VectorXd res(x.size());
        for(size_t i=0, imax=x.size(); i<imax; i++) {
          double val = x(i);
          double positivePart =     (val >= xmin and val < xmiddle);
          double negativePart = -1.*(val >= xmiddle and val <= xmax);
          res(i) = normfactor * (positivePart + negativePart);
        }
        return res;
      }

      template<class KernelFunction>
      static void evalKernel(
          Eigen::MatrixXd& kernelMatrix, KernelFunction& kernelFunction,
          size_t xIndex, size_t yIndex,
          const Eigen::VectorXd& x, const Eigen::VectorXd& xValues,
          const Eigen::VectorXd& y, const Eigen::VectorXd& yValues
          ) {
        // evaluate the kernel function
        size_t nquad = x.size();
        double xmin = x[0], xmax = x[nquad-1];
        double ymin = y[0], ymax = y[nquad-1];
        double quadweight = (ymax-ymin) * (xmax-xmin)
                          / ((nquad-1)*(nquad-1));
        double eval = 0.;
        for(size_t i = 0; i < nquad; i++) {
          Direction s_i = {cos(x[i]), sin(x[i])};
          for(size_t j = 0; j < nquad; j++) {
            Direction s_j = {cos(y[j]), sin(y[j])};
            // Integral over [-pi,pi]x[-pi,pi]
            eval += kernelFunction(s_i, s_j)
                  * quadweight * xValues[i] * yValues[j];
          }
        }
        kernelMatrix(xIndex, yIndex) = eval;
      }

      static Eigen::VectorXd computeQuadPoints(
          size_t jx, size_t kx, double r, size_t nquad) {
        double xmin = r *   kx   * std::exp2(-(double)jx+1) - r;
        double xmax = r * (kx+1) * std::exp2(-(double)jx+1) - r;
        Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(nquad, xmin, xmax);
        return x;
      }

    public:
      // Discrete Haar wavelet transform
      // data is used for in- and output and is assumed to be of size 2^n
      // for a natural number n.
      static void DWT(Eigen::VectorXd& data) {
        Eigen::VectorXd tmp(data.size());
        for(size_t len = data.size() >> 1; len > 0; len >>=1) {
          for(size_t i = 0; i < len; i++) {
            double sum  = data[2*i] + data[2*i+1];
            double diff = data[2*i] - data[2*i+1];
            tmp[i]     = .5 * sum;
            tmp[len+i] = .5 * diff;
          }
          data.head(2*len) = tmp.head(2*len);
        }
      }

      // Inverse discrete Haar wavelet transform
      // data is used for in- and output and is assumed to be of size 2^n
      // for a natural number n.
      static void IDWT(Eigen::VectorXd& data) {
        Eigen::VectorXd tmp(data.size());
        size_t max_len = data.size();
        for(size_t len = 1; len < max_len; len <<=1) {
          for(size_t i = 0; i < len; i++) {
            double sum  = data[i] + data[len+i];
            double diff = data[i] - data[len+i];
            tmp[2*i]   = sum;
            tmp[2*i+1] = diff;
          }
          data.head(2*len) = tmp.head(2*len);
        }
      }

      enum : unsigned int { dim = 2 };
      using Direction = FieldVector<double, dim>;

      HaarWavelet() = delete;
      HaarWavelet(const HaarWavelet&) = delete;

      template<class Function>
      HaarWavelet(const Function& kernelFunction, size_t num_s)
        : kernelMatrix(num_s, num_s),
          maxLevel(std::ilogb(num_s)),
          level(maxLevel) {
        using namespace Eigen;
        using namespace boost::math::constants;
        if((1u << maxLevel) != num_s)
          DUNE_THROW(MathError, "You are using " << num_s
              << " directions, but only powers of 2 are supported.");

        /* compute wavelet transform of the kernel matrix */
        double r = pi<double>();
        size_t nquad = 100;

        // jy, ky = 0
        {
          VectorXd y = computeQuadPoints(0, 0, r, nquad);
          VectorXd sfy = sf(y, r);
          // jx, kx = 0
          {
            VectorXd x = computeQuadPoints(0, 0, r, nquad);
            VectorXd sfx = sf(x, r);
            evalKernel(kernelMatrix, kernelFunction,
                0, 0, x, sfx, y, sfy);
          }
          for(size_t jx = 0; jx < maxLevel; jx++) {
            for(size_t kx = 0, kx_max = 1 << jx; kx < kx_max; kx++) {
              VectorXd x = computeQuadPoints(jx, kx, r, nquad);
              VectorXd wltx = wlt(jx, kx, x, r);
              size_t i = (1 << jx) + kx;
              evalKernel(kernelMatrix, kernelFunction,
                  i, 0, x, wltx, y, sfy);
            }
          }
        }
        for(size_t jy = 0; jy < maxLevel; jy++) {
          for(size_t ky = 0, ky_max = 1 << jy; ky < ky_max; ky++) {
            VectorXd y = computeQuadPoints(jy, ky, r, nquad);
            VectorXd wlty = wlt(jy, ky, y, r);
            // jx, kx = 0
            {
              VectorXd x = computeQuadPoints(0, 0, r, nquad);
              VectorXd sfx = sf(x, r);
              size_t j = (1 << jy) + ky;
              evalKernel(kernelMatrix, kernelFunction,
                  0, j, x, sfx, y, wlty);
            }
            for(size_t jx = 0; jx < maxLevel; jx++) {
              for(size_t kx = 0, kx_max = 1 << jx; kx < kx_max; kx++) {
                VectorXd x = computeQuadPoints(jx, kx, r, nquad);
                VectorXd wltx = wlt(jx, kx, x, r);
                size_t i = (1 << jx) + kx;
                size_t j = (1 << jy) + ky;
                evalKernel(kernelMatrix, kernelFunction,
                    i, j, x, wltx, y, wlty);
              }
            }
          }
        }
      }

      void applyToVector(Eigen::VectorXd& v) const {
        DWT(v);
        // TODO: trunkate matrix (and vector?) to level
        v = kernelMatrix * v;
        IDWT(v);
      }

      void setAccuracy(double accuracy) {
        // TODO: properly set level dependent on accuracy
        level = maxLevel;
      }

    private:
      Eigen::MatrixXd kernelMatrix;
      size_t maxLevel;
      size_t level;
  };
}

}

#endif // defined DUNE_DPG_RADIATIVE_TRANSFER_WAVELET_KERNEL_APPROXIMATION_HH
