// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_CHOLESKY_HH
#define DUNE_DPG_CHOLESKY_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

namespace Dune {

template<typename MatrixType>
class Cholesky
{
public:
  Cholesky<MatrixType>(MatrixType& matrix)
    : matrix(matrix)
  {
    unsigned int n = matrix.N();
    if(!(n==matrix.M()))
    {
      DUNE_THROW(Dune::Exception, "Cholesky only possible for square-matrices");
    }
    for (unsigned int k=0; k<n; k++)
    {
      double summe = matrix[k][k];
      for (unsigned int i=0; i<k; ++i)
      {
        summe -= (matrix[k][i]*matrix[k][i]*matrix[i][i]);
      }
      if (summe<1e-20)
      {
        std::cout << "summe = " << summe << std::endl;
        std::ofstream ofs("matrix");
        printmatrix(ofs, matrix, "Problematic Matrix in Cholesky", "--");
        DUNE_THROW(Dune::Exception, "Matrix for Cholesky not positive semidefinite");
      }
      matrix[k][k]=summe;
      for (unsigned int j=k+1; j<n; ++j)
      {
        summe = matrix[j][k];
        for (unsigned int i=0; i<k; ++i)
        {
          summe -= (matrix[j][i]*matrix[k][i]*matrix[i][i]); //TODO matrix[k][i]*matrix[i][i] zwischenspeichern guenstiger?
        }
        matrix[j][k] = (summe/matrix[k][k]);
      }
    }
  }

  void apply(MatrixType& rhsMatrix) const
  {
    const unsigned int m = rhsMatrix.M();
    const unsigned int n = matrix.N();
    if(!(n==rhsMatrix.N()))
    {
      DUNE_THROW(Dune::Exception, "rhsMatrix in Cholesky has wrong number of rows");
    }
    for (unsigned int im=0; im<m; im++)
    {
    // solve LY = B
      for (unsigned int j=0; j<n; j++)
      {
        double summe = 0;
        for (unsigned int i=0; i<j; i++)
        {
          summe+=matrix[j][i]*rhsMatrix[i][im];
        }
        rhsMatrix[j][im] = (rhsMatrix[j][im]-summe);
      }
     // solve L^TX = D^(-1)B
      for (int j=n-1; j>-1; j--)
      {
        rhsMatrix[j][im] = rhsMatrix[j][im]/matrix[j][j];
        double summe = 0;
        for (unsigned int i=j+1; i<n; i++)
        {
          summe+=matrix[i][j]*rhsMatrix[i][im];
        }
        rhsMatrix[j][im] = (rhsMatrix[j][im]-summe);
      }
    }
  }

  void apply(const MatrixType& rhsMatrix, MatrixType& solutionMatrix) const
  {
    const unsigned int m = rhsMatrix.M();
    const unsigned int n = matrix.N();
    solutionMatrix.setSize(n,m);
    if(!(n==rhsMatrix.N()))
    {
      DUNE_THROW(Dune::Exception, "rhsMatrix in Cholesky has wrong number of rows");
    }
    for (unsigned int im=0; im<m; im++)
    {
    // solve LY = B
      for (unsigned int j=0; j<n; j++)
      {
        double summe = 0;
        for (unsigned int i=0; i<j; i++)
        {
          summe+=matrix[j][i]*solutionMatrix[i][im];
        }
        solutionMatrix[j][im] = (rhsMatrix[j][im]-summe);
      }

     // solve L^TX = D^(-1)B
      for (int j=n-1; j>-1; j--)
      {
        solutionMatrix[j][im] = solutionMatrix[j][im]/matrix[j][j];
        double summe = 0;
        for (unsigned int i=j+1; i<n; i++)
        {
          summe+=matrix[i][j]*solutionMatrix[i][im];
        }
        solutionMatrix[j][im] = (solutionMatrix[j][im]-summe);
      }
    }
  }


  void apply(BlockVector<FieldVector<double,1> >& rhsVector) const
  {
    const unsigned int n = rhsVector.size();
    if(!(n==matrix.N()))
    {
      DUNE_THROW(Dune::Exception, "rhsVector in Cholesky has wrong number of rows");
    }
    // solve LY = B
      for (unsigned int j=0; j<n; j++)
      {
        double summe = 0;
        for (unsigned int i=0; i<j; i++)
        {
          summe+=matrix[j][i]*rhsVector[i];
        }
        rhsVector[j] = (rhsVector[j]-summe);
      }
     // solve L^TX = D^(-1)B
      for (int j=n-1; j>-1; j--)
      {
        rhsVector[j] = rhsVector[j]/matrix[j][j];
        double summe = 0;
        for (unsigned int i=j+1; i<n; i++)
        {
          summe+=matrix[i][j]*rhsVector[i];
        }
        rhsVector[j] = (rhsVector[j]-summe);
      }
  }

private:
  MatrixType& matrix;
};

} // end namespace Dune

#endif // DUNE_DPG_CHOLESKY_HH
