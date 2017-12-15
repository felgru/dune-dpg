// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_QKFACE2D_LOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_QKFACE2D_LOCALFINITEELEMENT_HH

#include <dune/common/version.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include "qkface2d/qkface2dlocalinterpolation.hh"
#include "qkface2d/qkface2dlocalbasis.hh"
#include "qkface2d/qkface2dlocalcoefficients.hh"

namespace Dune
{
  /** \brief General Lagrange finite element suitable for discontinuous faces for cubes with dimension 2 and polynomial order
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam k polynomial order
   */
  template<class D, class R, int k>
  class QkFace2DLocalFiniteElement {

    typedef QkFace2DLocalBasis<D,R,k> LocalBasis;
    typedef QkFace2DLocalCoefficients<k> LocalCoefficients;
    typedef QkFace2DLocalInterpolation<k,LocalBasis> LocalInterpolation;

  public:

    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<LocalBasis,QkFace2DLocalCoefficients<k>,LocalInterpolation> Traits;

    /** \todo Please doc me !
     */
    QkFace2DLocalFiniteElement ()
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
    { }
#else
    {
      gt.makeCube(2);
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
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
    static constexpr GeometryType type ()
    {
      return GeometryTypes::quadrilateral;
    }
#else
    GeometryType type () const
    {
      return gt;
    }
#endif

  private:
    LocalBasis basis;
    LocalCoefficients coefficients;
    LocalInterpolation interpolation;
#if not(DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6))
    GeometryType gt;
#endif
  };

}

#endif
