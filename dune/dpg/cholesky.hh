// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_CHOLESKY_HH
#define DUNE_DPG_CHOLESKY_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/io.hh>

namespace Dune {

template<typename MatrixType>
class Cholesky
{
public:
  Cholesky<MatrixType>(MatrixType& matrix)
    : matrix(matrix)
  {
    const unsigned int n = matrix.N();
    if(n != matrix.M())
    {
      DUNE_THROW(Dune::Exception,
          "Cholesky decomposition only possible for square-matrices.");
    }
    for (unsigned int k=0; k<n; k++)
    {
      double sum = matrix[k][k];
      for (unsigned int i=0; i<k; ++i)
      {
        sum -= matrix[k][i] * matrix[k][i] * matrix[i][i];
      }
      if (sum < 1e-20)
      {
        std::cout << "sum = " << sum << std::endl;
        std::ofstream ofs("matrix");
        printmatrix(ofs, matrix, "Problematic Matrix in Cholesky", "--");
        DUNE_THROW(Dune::Exception,
            "Matrix for Cholesky decomposition not positive semidefinite.");
      }
      matrix[k][k] = sum;
      for (unsigned int j=k+1; j<n; ++j)
      {
        sum = matrix[j][k];
        for (unsigned int i=0; i<k; ++i)
        {
          // TODO: More efficient to cache matrix[k][i]*matrix[i][i]?
          sum -= matrix[j][i] * matrix[k][i] * matrix[i][i];
        }
        matrix[j][k] = sum / matrix[k][k];
      }
    }
  }

  void apply(MatrixType& rhsMatrix) const
  {
    const unsigned int m = rhsMatrix.M();
    const unsigned int n = matrix.N();
    if(n != rhsMatrix.N())
    {
      DUNE_THROW(Dune::Exception,
          "rhsMatrix in Cholesky has wrong number of rows.");
    }
    for (unsigned int im=0; im<m; im++)
    {
      // solve LY = B
      for (unsigned int j=0; j<n; j++)
      {
        double sum = 0;
        for (unsigned int i=0; i<j; i++)
        {
          sum += matrix[j][i] * rhsMatrix[i][im];
        }
        rhsMatrix[j][im] = rhsMatrix[j][im] - sum;
      }
      // solve L^TX = D^(-1)B
      for (int j=n-1; j>-1; j--)
      {
        rhsMatrix[j][im] = rhsMatrix[j][im]/matrix[j][j];
        double sum = 0;
        for (unsigned int i=j+1; i<n; i++)
        {
          sum += matrix[i][j] * rhsMatrix[i][im];
        }
        rhsMatrix[j][im] = rhsMatrix[j][im] - sum;
      }
    }
  }

  void apply(const MatrixType& rhsMatrix, MatrixType& solutionMatrix) const
  {
    const unsigned int m = rhsMatrix.M();
    const unsigned int n = matrix.N();
    solutionMatrix.setSize(n,m);
    if(n != rhsMatrix.N())
    {
      DUNE_THROW(Dune::Exception,
          "rhsMatrix in Cholesky has wrong number of rows.");
    }
    for (unsigned int im=0; im<m; im++)
    {
      // solve LY = B
      for (unsigned int j=0; j<n; j++)
      {
        double sum = 0;
        for (unsigned int i=0; i<j; i++)
        {
          sum += matrix[j][i] * solutionMatrix[i][im];
        }
        solutionMatrix[j][im] = rhsMatrix[j][im] - sum;
      }

      // solve L^TX = D^(-1)B
      for (int j=n-1; j>-1; j--)
      {
        solutionMatrix[j][im] = solutionMatrix[j][im]/matrix[j][j];
        double sum = 0;
        for (unsigned int i=j+1; i<n; i++)
        {
          sum += matrix[i][j] * solutionMatrix[i][im];
        }
        solutionMatrix[j][im] = solutionMatrix[j][im] - sum;
      }
    }
  }


  void apply(BlockVector<FieldVector<double,1> >& rhsVector) const
  {
    const unsigned int n = rhsVector.size();
    if(n != matrix.N())
    {
      DUNE_THROW(Dune::Exception,
          "rhsVector in Cholesky has wrong number of rows.");
    }
      // solve LY = B
      for (unsigned int j=0; j<n; j++)
      {
        double sum = 0;
        for (unsigned int i=0; i<j; i++)
        {
          sum += matrix[j][i] * rhsVector[i];
        }
        rhsVector[j] = rhsVector[j] - sum;
      }
      // solve L^TX = D^(-1)B
      for (int j=n-1; j>-1; j--)
      {
        rhsVector[j] = rhsVector[j]/matrix[j][j];
        double sum = 0;
        for (unsigned int i=j+1; i<n; i++)
        {
          sum += matrix[i][j] * rhsVector[i];
        }
        rhsVector[j] = rhsVector[j] - sum;
      }
  }

private:
  MatrixType& matrix;
};

} // end namespace Dune

#endif // DUNE_DPG_CHOLESKY_HH
