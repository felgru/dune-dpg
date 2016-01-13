// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_REFINEDOPTIMALTEST_LOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_REFINEDOPTIMALTEST_LOCALFINITEELEMENT_HH

#include "refinedoptimaltest/refinedoptimaltestlocalinterpolation.hh"
#include "refinedoptimaltest/refinedoptimaltestlocalbasis.hh"
#include "refinedoptimaltest/refinedoptimaltestlocalcoefficients.hh"

namespace Dune
{
  /** \brief refined optimal test functions
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam d dimension of the reference element
   * \tparam level level of the refinement
   * \tparam RefinementConstants
   * \tparam TestSearchSpace
   */
  template<class D, class R, int d,
           int level, class RefinementConstants, class TestSearchSpace>
  class RefinedOptimalTestLocalFiniteElement {

    typedef RefinedOptimalTestLocalBasis<D,R,d,level,RefinementConstants,
              typename TestSearchSpace::Traits::LocalBasisType> LocalBasis;
    typedef RefinedOptimalTestLocalCoefficients<d,level> LocalCoefficients;
    typedef RefinedOptimalTestLocalInterpolation<d,level,LocalBasis> LocalInterpolation;
    typedef Matrix<FieldMatrix<double,1,1> > MatrixType;

  public:

    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<LocalBasis,LocalCoefficients,LocalInterpolation> Traits;

    /** \todo Please doc me !
     */
    RefinedOptimalTestLocalFiniteElement () = delete;

    RefinedOptimalTestLocalFiniteElement (MatrixType& coeffMat,
                                   const TestSearchSpace& testSearchSpace,
                                   std::vector<LocalKey>& localKeyList)
    : basis(testSearchSpace.localBasis(), testSearchSpace.type(),
            coeffMat, 0),
      coefficients(localKeyList),
      gt(testSearchSpace.type())
    { }

    RefinedOptimalTestLocalFiniteElement (MatrixType& coeffMat,
                                   const TestSearchSpace& testSearchSpace,
                                   size_t offset = 0)
    : basis(testSearchSpace.localBasis(), testSearchSpace.type(),
            coeffMat, offset),
      gt(testSearchSpace.type())
    { }

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

    RefinedOptimalTestLocalFiniteElement* clone () const
    {
      return new RefinedOptimalTestLocalFiniteElement(*this);
    }

  private:
    LocalBasis basis;
    LocalCoefficients coefficients;
    LocalInterpolation interpolation;
    GeometryType gt;
  };

}

#endif
