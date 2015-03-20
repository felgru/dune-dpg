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
  template<class D, class R, int d, class EnrichedTestspace>
  class OptimalTestLocalFiniteElement {

    typedef OptimalTestLocalBasis<D,R,d,EnrichedTestspace> LocalBasis;
    typedef OptimalTestLocalCoefficients<d> LocalCoefficients;
    typedef OptimalTestLocalInterpolation<d,LocalBasis> LocalInterpolation;
    typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

  public:

    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<LocalBasis,OptimalTestLocalCoefficients<d>,LocalInterpolation> Traits;

    /** \todo Please doc me !
     */
    OptimalTestLocalFiniteElement () = delete;

    OptimalTestLocalFiniteElement (MatrixType* coeffMat, std::vector<LocalKey>* localKeyList)
    :enrichedTestspace(),
    basis(&enrichedTestspace, coeffMat),
    coefficients(localKeyList)
    {
      gt=enrichedTestspace.type();
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
    EnrichedTestspace enrichedTestspace;
    MatrixType* coefficientMatrix;
    LocalBasis basis;
    LocalCoefficients coefficients;
    LocalInterpolation interpolation;
    GeometryType gt;
  };

}

#endif
