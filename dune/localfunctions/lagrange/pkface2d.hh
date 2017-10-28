// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PKFACE2DLOCALFINITEELEMENT_HH
#define DUNE_PKFACE2DLOCALFINITEELEMENT_HH

#include <cstddef>

#include <dune/common/version.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include "pkface2d/pkface2dlocalbasis.hh"
#include "pkface2d/pkface2dlocalcoefficients.hh"
#include "pkface2d/pkface2dlocalinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R, unsigned int k>
  class PkFace2DLocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<PkFace2DLocalBasis<D,R,k>,
        PkFace2DLocalCoefficients<k>,
        PkFace2DLocalInterpolation<PkFace2DLocalBasis<D,R,k> > > Traits;

    /** \todo Please doc me !
     */
    PkFace2DLocalFiniteElement ()
#if DUNE_VERSION_NEWER(DUNE_GRID,2,6)
    { }
#else
    {
      gt.makeTriangle();
    }
#endif

    /** \todo Please doc me !
     */
    PkFace2DLocalFiniteElement (int variant)
      : coefficients(variant)
#if DUNE_VERSION_NEWER(DUNE_GRID,2,6)
    { }
#else
    {
      gt.makeTriangle();
    }
#endif

    /** Constructor for six variants with permuted vertices.

        \param vertexmap The permutation of the vertices.  This
        can for instance be generated from the global indices of
        the vertices by reducing those to the integers 0...2
     */
    PkFace2DLocalFiniteElement (const unsigned int vertexmap[3])
      : coefficients(vertexmap)
#if DUNE_VERSION_NEWER(DUNE_GRID,2,6)
    { }
#else
    {
      gt.makeTriangle();
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
      return GeometryTypes::triangle;
    }
#else
    GeometryType type () const
    {
      return gt;
    }
#endif

  private:
    PkFace2DLocalBasis<D,R,k> basis;
    PkFace2DLocalCoefficients<k> coefficients;
    PkFace2DLocalInterpolation<PkFace2DLocalBasis<D,R,k> > interpolation;
#if not(DUNE_VERSION_NEWER(DUNE_GRID,2,6))
    GeometryType gt;
#endif
  };

}

#endif
