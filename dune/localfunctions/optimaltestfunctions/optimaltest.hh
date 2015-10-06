// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_OPTIMALTEST_LOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_OPTIMALTEST_LOCALFINITEELEMENT_HH

#include "optimaltest/optimaltestlocalinterpolation.hh"
#include "optimaltest/optimaltestlocalbasis.hh"
#include "optimaltest/optimaltestlocalcoefficients.hh"

namespace Dune
{
  /** \brief optimaltestfunctions
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam d dimension of the reference element
   * \tparam k polynomial order
   */
  template<class D, class R, int d, class TestSearchSpace>
  class OptimalTestLocalFiniteElement {

    typedef OptimalTestLocalBasis<D,R,d,TestSearchSpace> LocalBasis;
    typedef OptimalTestLocalCoefficients<d> LocalCoefficients;
    typedef OptimalTestLocalInterpolation<d,LocalBasis> LocalInterpolation;
    typedef Matrix<FieldMatrix<double,1,1> > MatrixType;

  public:

    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<LocalBasis,OptimalTestLocalCoefficients<d>,LocalInterpolation> Traits;

    /** \todo Please doc me !
     */
    OptimalTestLocalFiniteElement () = delete;

    OptimalTestLocalFiniteElement (MatrixType* coeffMat,
                                   TestSearchSpace* testSearchSpace,
                                   std::vector<LocalKey>* localKeyList)
    :testSearchSpace(testSearchSpace),
    basis(testSearchSpace, coeffMat, 0, coeffMat->N()),
    coefficients(localKeyList)
    {
      gt=testSearchSpace->type();
    }

    OptimalTestLocalFiniteElement (MatrixType* coeffMat,
                                   TestSearchSpace* testSearchSpace)
    :testSearchSpace(testSearchSpace),
    basis(testSearchSpace, coeffMat, 0, coeffMat->N())
    {
      gt=testSearchSpace->type();
    }

    OptimalTestLocalFiniteElement (MatrixType* coeffMat,
                                   TestSearchSpace* testSearchSpace,
                                   size_t offset,
                                   size_t k)
    :testSearchSpace(testSearchSpace),
    basis(testSearchSpace, coeffMat, offset, k)
    {
      gt=testSearchSpace->type();
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation;
    }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return basis.size();
    }

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      return gt;
    }

    OptimalTestLocalFiniteElement* clone () const
    {
      return new OptimalTestLocalFiniteElement(*this);
    }

  private:
    TestSearchSpace* testSearchSpace;
    MatrixType* coefficientMatrix;
    LocalBasis basis;
    LocalCoefficients coefficients;
    LocalInterpolation interpolation;
    GeometryType gt;
  };

}

#endif
