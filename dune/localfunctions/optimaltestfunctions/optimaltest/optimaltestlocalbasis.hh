// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_OPTIMALTESTLOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_OPTIMALTESTLOCALBASIS_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/power.hh>

#include <dune/geometry/type.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/io.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>


namespace Dune
{
  /**@ingroup LocalBasisImplementation
     \brief Lagrange shape functions of order k on the reference cube.

     Also known as \f$Q^k\f$.

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.
     \tparam k Polynomial degree
     \tparam d Dimension of the cube

     \nosubgrouping
   */
  template<class D, class R, int d, class TestSearchSpace>
  class OptimalTestLocalBasis
  {
    typedef Matrix<FieldMatrix<double,1,1> > MatrixType;

  public:
    typedef LocalBasisTraits<D,d,Dune::FieldVector<D,d>,R,1,Dune::FieldVector<R,1>,Dune::FieldMatrix<R,1,d> > Traits;

    OptimalTestLocalBasis (TestSearchSpace* testSearchSpace,
                           MatrixType* coeffMat, size_t offset, size_t k)
        : testSearchSpace(testSearchSpace),
          coefficientMatrix(coeffMat),
          offset(offset),
          k(k) { }

    //! \brief number of shape functions
    unsigned int size () const
    {
      return coefficientMatrix->M();
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());
      std::vector<FieldVector<double,1> > valuesTestSearchSpace(k);
      testSearchSpace->localBasis().evaluateFunction(in, valuesTestSearchSpace);
      for (size_t i=0; i<size(); i++)
      {
        out[i] = 0;
        for (size_t j=0; j<k; ++j)
        {
          out[i]+=((*coefficientMatrix)[offset+j][i][0]*valuesTestSearchSpace[j][0]);
        }
      }
    }

    /** \brief Evaluate Jacobian of all shape functions
     * \param in position where to evaluate
     * \param out The return value
     */
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,
                      std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());
      std::vector<typename Traits::JacobianType> JacobianTestSearchSpace(k);
      testSearchSpace->localBasis().evaluateJacobian(in, JacobianTestSearchSpace);
      // Loop over all shape functions
      for (size_t i=0; i<size(); i++)
      {
        // Loop over all coordinate directions
        for (int b=0; b<d; b++)
        {
          // Initialize: the overall expression is a product
          // if j-th bit of i is set to -1, else 1
          out[i][0][b] = 0;
          for (size_t j=0; j<k; ++j)
          {
            out[i][0][b]+=((*coefficientMatrix)[offset+j][i][0]*JacobianTestSearchSpace[j][0][b]);
          }
        }
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return testSearchSpace->localBasis().order();
    }
  private:
    TestSearchSpace* testSearchSpace;
    MatrixType* coefficientMatrix;
    size_t offset;
    size_t k;
  };
}

#endif
