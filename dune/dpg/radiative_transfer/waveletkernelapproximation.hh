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
  namespace HaarWavelet {
      // Scaling function $\phi_{0,0}$ on the interval [-r,r)
      Eigen::VectorXd sf(const Eigen::VectorXd& x, double r) {
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
      Eigen::VectorXd wlt(size_t j, size_t k,
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
      double evalKernel(
          KernelFunction& kernelFunction,
          const Eigen::VectorXd& x, const Eigen::VectorXd& xValues,
          double quadweightx,
          const Eigen::VectorXd& y, const Eigen::VectorXd& yValues,
          double quadweighty
          ) {
        enum : unsigned int { dim = 2 };
        using Direction = FieldVector<double, dim>;

        // evaluate the kernel function
        const size_t nquadx = x.size();
        const size_t nquady = y.size();
        const double quadweight = quadweighty * quadweightx;
        double eval = 0.;
        for(size_t i = 0; i < nquadx; i++) {
          Direction s_i = {cos(x[i]), sin(x[i])};
          for(size_t j = 0; j < nquady; j++) {
            Direction s_j = {cos(y[j]), sin(y[j])};
            // Integral over [-pi,pi]x[-pi,pi]
            eval += kernelFunction(s_i, s_j)
                  * quadweight * xValues[i] * yValues[j];
          }
        }
        // scale integral to be a probability measure
        // TODO: maybe the scaling factor should be included in the
        //       kernel, as in the original definition of
        //       Henyey and Greenstein.
        return eval / (2*boost::math::constants::pi<double>());
      }

      std::tuple<Eigen::VectorXd, double> computeQuadPoints(
          size_t jx, size_t kx, double r, size_t maxLevel) {
        double xmin = r *   kx   * std::exp2(-(double)jx+1) - r;
        double xmax = r * (kx+1) * std::exp2(-(double)jx+1) - r;
        size_t nquad = 1 << (maxLevel - jx);
        Eigen::VectorXd x(nquad);
        for(size_t i=0; i < nquad; i++) {
          x[i] = xmin + (2*i+1)/(2.*nquad) * (xmax - xmin);
        }
        double quadweight = (xmax - xmin) / nquad;
        return std::make_tuple(x, quadweight);
      }

      // Discrete Haar wavelet transform
      // data is used for in- and output and is assumed to be of size 2^n
      // for a natural number n.
      void DWT(Eigen::VectorXd& data) {
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
        double scaling_factor = sqrt(2*boost::math::constants::pi<double>());
        data.segment(0, 1) *= scaling_factor;
        for(size_t len=1, max_len=data.size()>>1; len < max_len; len <<=1) {
          scaling_factor /= sqrt(2);
          data.segment(len, len) *= scaling_factor;
        }
      }

      // Inverse discrete Haar wavelet transform
      // data is used for in- and output and is assumed to be of size 2^n
      // for a natural number n.
      void IDWT(Eigen::VectorXd& data) {
        Eigen::VectorXd tmp(data.size());
        double scaling_factor = 1./sqrt(2*boost::math::constants::pi<double>());
        data.segment(0, 1) *= scaling_factor;
        for(size_t len=1, max_len=data.size()>>1; len < max_len; len <<=1) {
          scaling_factor *= sqrt(2);
          data.segment(len, len) *= scaling_factor;
        }
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

      template<class Function>
      Eigen::MatrixXd waveletKernelMatrix(const Function& kernelFunction,
                                          size_t maxLevel)
      {
        using namespace Eigen;
        using namespace boost::math::constants;

        const size_t num_s = 1ul << maxLevel;
        MatrixXd kernelMatrix(num_s, num_s);

        /* compute wavelet transform of the kernel matrix */
        const double r = pi<double>();

        // jy, ky = 0
        {
          VectorXd y; double qwy;
          std::tie(y, qwy) = computeQuadPoints(0, 0, r, maxLevel);
          VectorXd sfy = sf(y, r);
          // jx, kx = 0
          {
            VectorXd x; double qwx;
            std::tie(x, qwx) = computeQuadPoints(0, 0, r, maxLevel);
            VectorXd sfx = sf(x, r);
            kernelMatrix(0, 0)
                = evalKernel(kernelFunction, x, sfx, qwx, y, sfy, qwy);
          }
          for(size_t jx = 0; jx < maxLevel; jx++) {
            for(size_t kx = 0, kx_max = 1 << jx; kx < kx_max; kx++) {
              VectorXd x; double qwx;
              std::tie(x, qwx) = computeQuadPoints(jx, kx, r, maxLevel);
              VectorXd wltx = wlt(jx, kx, x, r);
              const size_t i = (1 << jx) + kx;
              kernelMatrix(i, 0)
                  = evalKernel(kernelFunction, x, wltx, qwx, y, sfy, qwy);
            }
          }
        }
        for(size_t jy = 0; jy < maxLevel; jy++) {
          for(size_t ky = 0, ky_max = 1 << jy; ky < ky_max; ky++) {
            VectorXd y; double qwy;
            std::tie(y, qwy) = computeQuadPoints(jy, ky, r, maxLevel);
            VectorXd wlty = wlt(jy, ky, y, r);
            // jx, kx = 0
            {
              VectorXd x; double qwx;
              std::tie(x, qwx) = computeQuadPoints(0, 0, r, maxLevel);
              VectorXd sfx = sf(x, r);
              const size_t j = (1 << jy) + ky;
              kernelMatrix(0, j)
                  = evalKernel(kernelFunction, x, sfx, qwx, y, wlty, qwy);
            }
            for(size_t jx = 0; jx < maxLevel; jx++) {
              for(size_t kx = 0, kx_max = 1 << jx; kx < kx_max; kx++) {
                VectorXd x; double qwx;
                std::tie(x, qwx) = computeQuadPoints(jx, kx, r, maxLevel);
                VectorXd wltx = wlt(jx, kx, x, r);
                const size_t i = (1 << jx) + kx;
                const size_t j = (1 << jy) + ky;
                kernelMatrix(i, j)
                    = evalKernel(kernelFunction, x, wltx, qwx, y, wlty, qwy);
              }
            }
          }
        }
        return kernelMatrix;
      }

      class MatrixCompression {
      public:
        enum : unsigned int { dim = 2 };
        using Direction = FieldVector<double, dim>;

        MatrixCompression() = delete;
        MatrixCompression(const MatrixCompression&) = delete;

        template<class Function>
        MatrixCompression(const Function& kernelFunction, size_t num_s)
          : kernelMatrix(waveletKernelMatrix(kernelFunction,
                                             std::ilogb(num_s))),
            maxLevel(std::ilogb(num_s)),
            level(maxLevel) {
          if((1u << maxLevel) != num_s)
            DUNE_THROW(MathError, "You are using " << num_s
                << " directions, but only powers of 2 are supported.");
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

}

#endif // defined DUNE_DPG_RADIATIVE_TRANSFER_WAVELET_KERNEL_APPROXIMATION_HH
