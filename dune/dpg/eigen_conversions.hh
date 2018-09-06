// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_EIGEN_CONVERSIONS_HH
#define DUNE_EIGEN_CONVERSIONS_HH

#include <Eigen/Core>
#include <Eigen/Dense>

namespace Dune {

  using EigenRowMatrix
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

  template<class Matrix>
  EigenRowMatrix
  toEigenMatrix(const Matrix& matrix) {
    EigenRowMatrix eigenMatrix(matrix.N(), matrix.M());
    for(size_t i = 0; i < matrix.N(); i++)
      for(size_t j = 0; j < matrix.M(); j++)
        eigenMatrix(i, j) = matrix[i][j][0];
    return eigenMatrix;
  }

  template<class Matrix>
  void copyEigenMatrixToDuneMatrix(const EigenRowMatrix& eigenMatrix,
                                   Matrix& matrix) {
    matrix.setSize(eigenMatrix.rows(), eigenMatrix.cols());
    for(size_t i = 0; i < matrix.N(); i++)
      for(size_t j = 0; j < matrix.M(); j++)
        matrix[i][j][0] = eigenMatrix(i, j);
  }

  template<class Vector>
  Eigen::VectorXd
  toEigenVector(const Vector& v) {
    Eigen::VectorXd eigenV(v.N());
    for(size_t i = 0; i < v.N(); i++)
      eigenV(i) = v[i][0];
    return eigenV;
  }

  template<class Vector>
  void copyEigenVectorToDuneVector(const Eigen::VectorXd eigenV,
                                   Vector& v) {
    v.resize(eigenV.rows());
    for(size_t i = 0; i < v.N(); i++)
      v[i][0] = eigenV(i);
  }

}
#endif
