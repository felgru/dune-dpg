// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_QKTRACE_LOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_QKTRACE_LOCALFINITEELEMENT_HH

#include <dune/common/version.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include "qktrace/qktracelocalinterpolation.hh"
#include "qktrace/qktracelocalbasis.hh"
#include "qktrace/qktracelocalcoefficients.hh"

namespace Dune
{
  /** \brief General Lagrange finite element suitable for continuous traces for cubes with arbitrary dimension and polynomial order
   *   \note The general class QkLocalCoefficients is available for k>0 in dimensions 2 and 3 only
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam d dimension of the reference element
   * \tparam k polynomial order
   */
  template<class D, class R, int d, int k>
  class QkTraceLocalFiniteElement {

    typedef QkTraceLocalBasis<D,R,k,d> LocalBasis;
    typedef QkTraceLocalCoefficients<k,d> LocalCoefficients;
    typedef QkTraceLocalInterpolation<k,d,LocalBasis> LocalInterpolation;

  public:

    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<LocalBasis,QkTraceLocalCoefficients<k,d>,LocalInterpolation> Traits;

    /** \todo Please doc me !
     */
    QkTraceLocalFiniteElement ()
#if DUNE_VERSION_NEWER(DUNE_GRID,2,6)
    { }
#else
    {
      gt.makeCube(d);
    }
#endif

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
#if DUNE_VERSION_NEWER(DUNE_GRID,2,6)
    static constexpr GeometryType type ()
    {
      return GeometryTypes::cube(d);
    }
#else
    GeometryType type () const
    {
      return gt;
    }
#endif

    QkTraceLocalFiniteElement* clone () const
    {
      return new QkTraceLocalFiniteElement(*this);
    }

  private:
    LocalBasis basis;
    LocalCoefficients coefficients;
    LocalInterpolation interpolation;
#if not(DUNE_VERSION_NEWER(DUNE_GRID,2,6))
    GeometryType gt;
#endif
  };


  /** \brief General Lagrange finite element suitable for functions living only on the faces for cubes with arbitrary dimension and polynomial order
   *   \note The general class QkLocalCoefficients is available for k>0 in dimensions 2 and 3 only
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam d dimension of the reference element
   * \tparam k polynomial order
   */

}

#endif
