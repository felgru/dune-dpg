// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LEASTSQUARES_HH
#define DUNE_LEASTSQUARES_HH

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/QR>

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
      using Eigen::Dynamic;
      using RowMatrix
          = Eigen::Matrix<double, Dynamic, Dynamic, Eigen::RowMajor>;
      RowMatrix eigenA(A.N(), A.M());
      for(size_t i = 0; i < A.N(); i++)
        for(size_t j = 0; j < A.M(); j++)
          eigenA(i, j) = A[i][j][0];
      Eigen::ColPivHouseholderQR<Eigen::Ref<RowMatrix>>
        decomposition(eigenA);
      RowMatrix eigenB(b.N(), b.M());
      for(size_t i = 0; i < b.N(); i++)
        for(size_t j = 0; j < b.M(); j++)
          eigenB(i, j) = b[i][j][0];
      auto solution = decomposition.solve(eigenB);
      x.setSize(eigenA.rows(), eigenB.cols());
      for(size_t i = 0; i < x.N(); i++)
        for(size_t j = 0; j < x.M(); j++)
          x[i][j][0] = solution(i, j);
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
      using Eigen::Dynamic;
      using RowMatrix
          = Eigen::Matrix<double, Dynamic, Dynamic, Eigen::RowMajor>;
      RowMatrix eigenA(A.N(), A.M());
      for(size_t i = 0; i < A.N(); i++)
        for(size_t j = 0; j < A.M(); j++)
          eigenA(i, j) = A[i][j][0];
      Eigen::ColPivHouseholderQR<Eigen::Ref<RowMatrix>>
        decomposition(eigenA);
      Eigen::VectorXd eigenB(b.N());
      for(size_t i = 0; i < b.N(); i++)
        eigenB(i) = b[i][0];
      auto solution = decomposition.solve(eigenB);
      b.resize(eigenA.rows(), false);
      for(size_t i = 0; i < b.N(); i++)
        b[i][0] = solution(i);
}

} // end namespace Dune

#endif // DUNE_LEASTSQUARES_HH
