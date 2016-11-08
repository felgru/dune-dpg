// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_SVD_KERNEL_APPROXIMATION_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_SVD_KERNEL_APPROXIMATION_HH

#include <string>

#include <dune/common/fvector.hh>

#include <boost/math/constants/constants.hpp>

#include <Eigen/Core>
#include <Eigen/SVD>

namespace Dune {

namespace ScatteringKernelApproximation {
  class SVD {
    public:
      enum : unsigned int { dim = 2 };
      using Direction = FieldVector<double, dim>;

      SVD() = delete;
      SVD(const SVD&) = delete;

      template<class Function>
      SVD(const Function& kernel, size_t num_s)
        : kernelSVD(num_s, num_s, Eigen::ComputeThinU | Eigen::ComputeThinV),
          rank(num_s) {
        using namespace Eigen;
        using namespace boost::math::constants;
        MatrixXd kernelMatrix(num_s, num_s);
        for(size_t j = 0; j < num_s; ++j) {
          Direction s_j = {cos(2*pi<double>()*j/num_s),
                           sin(2*pi<double>()*j/num_s)};
          for(size_t i = 0; i < num_s; ++i) {
            Direction s_i = {cos(2*pi<double>()*i/num_s),
                             sin(2*pi<double>()*i/num_s)};
            // TODO: maybe use a higher order quadrature
            kernelMatrix(i,j) = kernel(s_i, s_j)/num_s;
          }
        }
        /* initialize SVD of kernel (using Eigen) */
        kernelSVD.compute(kernelMatrix);
      }

      void applyToVector(Eigen::VectorXd& v) const {
        v = kernelSVD.matrixU().leftCols(rank)
          * kernelSVD.singularValues().head(rank).asDiagonal()
          * kernelSVD.matrixV().leftCols(rank).adjoint() * v;
      }

      void setAccuracy(double accuracy) {
        using namespace Eigen;
        VectorXd singularValues = kernelSVD.singularValues();
        size_t i = singularValues.size() - 1;
        double err = 0,
               rank_err = singularValues(i) * singularValues(i);
        accuracy = accuracy * accuracy;
        while (err + rank_err < accuracy && i > 0) {
          err += rank_err;
          i -= 1;
          rank_err = singularValues(i) * singularValues(i);
        }
        rank = i+1;
        // TODO: If accuracy is low enough to allow rank = 0,
        //       this gives rank = 1.
      }

      std::string info() const {
        std::string s = "SVD approximation with rank "
                        + std::to_string(rank);
        return s;
      }

    private:
      Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::NoQRPreconditioner> kernelSVD;
      size_t rank;
  };
}

}

#endif // defined DUNE_DPG_RADIATIVE_TRANSFER_SVD_KERNEL_APPROXIMATION_HH
