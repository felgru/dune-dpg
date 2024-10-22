// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_CHOLESKY_HH
#define DUNE_DPG_CHOLESKY_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/io.hh>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "eigen_conversions.hh"

namespace Dune {

template<typename MatrixType>
class Cholesky
{
public:
  Cholesky(const MatrixType& matrix)
    : eigenMatrix(toEigenMatrix(matrix))
    , ldltDecomposition(eigenMatrix)
  {}

  void apply(MatrixType& rhsMatrix) const
  {
    const EigenRowMatrix eigenRhs = toEigenMatrix(rhsMatrix);
    const EigenRowMatrix solution = ldltDecomposition.solve(eigenRhs);
    copyEigenMatrixToDuneMatrix(solution, rhsMatrix);
  }

  void apply(const MatrixType& rhsMatrix, MatrixType& solutionMatrix) const
  {
    const EigenRowMatrix eigenRhs = toEigenMatrix(rhsMatrix);
    const EigenRowMatrix solution = ldltDecomposition.solve(eigenRhs);
    copyEigenMatrixToDuneMatrix(solution, solutionMatrix);
  }


  void apply(BlockVector<FieldVector<double,1> >& rhsVector) const
  {
    const Eigen::VectorXd eigenRhs = toEigenVector(rhsVector);
    const Eigen::VectorXd solution = ldltDecomposition.solve(eigenRhs);
    copyEigenVectorToDuneVector(solution, rhsVector);
  }

private:
  EigenRowMatrix eigenMatrix;
  Eigen::LDLT<Eigen::Ref<EigenRowMatrix>> ldltDecomposition;
};

} // end namespace Dune

#endif // DUNE_DPG_CHOLESKY_HH
