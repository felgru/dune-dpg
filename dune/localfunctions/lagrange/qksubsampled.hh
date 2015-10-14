// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_QK_SUBSAMPLED_LOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_QK_SUBSAMPLED_LOCALFINITEELEMENT_HH

#include "qksubsampled/qksubsampledlocalinterpolation.hh"
#include "qksubsampled/qksubsampledlocalbasis.hh"
#include "qksubsampled/qksubsampledlocalcoefficients.hh"

namespace Dune
{
  /** \brief Subsampled Lagrange finite element for cubes with arbitrary
   *         dimension and polynomial order
   *   \note The general class QkSubsampledLocalCoefficients is available
   *         for k>0 in dimensions 2 and 3 only
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam d dimension of the reference element
   * \tparam s number of samples
   * \tparam k polynomial order
   */
  template<class D, class R, int d, int s, int k>
  class QkSubsampledLocalFiniteElement {

    typedef QkSubsampledLocalBasis<D,R,s,k,d> LocalBasis;
    typedef QkSubsampledLocalCoefficients<s,k,d> LocalCoefficients;
    typedef QkSubsampledLocalInterpolation<s,k,d,LocalBasis> LocalInterpolation;

  public:

    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<LocalBasis,
                                     LocalCoefficients,
                                     LocalInterpolation>  Traits;

    /** \todo Please doc me !
     */
    QkSubsampledLocalFiniteElement ()
    {
      gt.makeCube(d);
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

    QkSubsampledLocalFiniteElement* clone () const
    {
      return new QkSubsampledLocalFiniteElement(*this);
    }

  private:
    LocalBasis basis;
    LocalCoefficients coefficients;
    LocalInterpolation interpolation;
    GeometryType gt;
  };

}

#endif
