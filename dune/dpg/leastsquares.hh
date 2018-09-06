// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LEASTSQUARES_HH
#define DUNE_LEASTSQUARES_HH

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/QR>

#include "eigen_conversions.hh"

namespace Dune {

/**
 * solves least squares problem ||AX-B|| → min
 *
 * Solve a least-squares problem for X
 * to make sure that it is well defined even when A
 * is severly ill conditioned.
 */
template<class Matrix>
void solveLeastSquares(const Matrix& A, const Matrix& b, Matrix& x)
{
  EigenRowMatrix eigenA = toEigenMatrix(A);
  Eigen::ColPivHouseholderQR<Eigen::Ref<EigenRowMatrix>>
      decomposition(eigenA);
  EigenRowMatrix eigenB = toEigenMatrix(b);
  const auto solution = decomposition.solve(eigenB);
  copyEigenMatrixToDuneMatrix(solution, x);
}

/**
 * solves least squares problem ||Ax-b|| → min
 *
 * Solve a least-squares problem for x
 * to make sure that it is well defined even when A
 * is severly ill conditioned.
 *
 * \note This destroys the input b and overwrites it with the
 *       least squares solution.
 */
template<class Matrix, class Vector>
void solveLeastSquares(const Matrix& A, Vector& b)
{
  EigenRowMatrix eigenA = toEigenMatrix(A);
  Eigen::ColPivHouseholderQR<Eigen::Ref<EigenRowMatrix>>
      decomposition(eigenA);
  Eigen::VectorXd eigenB = toEigenVector(b);
  const auto solution = decomposition.solve(eigenB);
  copyEigenVectorToDuneVector(solution, b);
}

} // end namespace Dune

#endif // DUNE_LEASTSQUARES_HH
