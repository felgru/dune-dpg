// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_QKFACE2D_LOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_QKFACE2D_LOCALFINITEELEMENT_HH

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
      : gt(GeometryTypes::quadrilateral)
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

    QkFace2DLocalFiniteElement* clone () const
    {
      return new QkFace2DLocalFiniteElement(*this);
    }

  private:
    LocalBasis basis;
    LocalCoefficients coefficients;
    LocalInterpolation interpolation;
    GeometryType gt;
  };

}

#endif
