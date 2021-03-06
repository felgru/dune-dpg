// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PKDGSUBSAMPLED2D_HH
#define DUNE_PKDGSUBSAMPLED2D_HH

#include <cstddef>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include "pkdgsubsampled2d/pkdgsubsampled2dlocalbasis.hh"
#include "pkdgsubsampled2d/pkdgsubsampled2dlocalcoefficients.hh"
#include "pkdgsubsampled2d/pkdgsubsampled2dlocalinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R, unsigned int s, unsigned int k>
  class PkDGSubsampled2DLocalFiniteElement
  {

    typedef PkDGSubsampled2DLocalBasis<D,R,s,k>     LocalBasis;
    typedef PkDGSubsampled2DLocalCoefficients<s,k>  LocalCoefficients;
    typedef PkDGSubsampled2DLocalInterpolation<s,k,
            PkDGSubsampled2DLocalBasis<D,R,s,k> >   LocalInterpolation;

  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<LocalBasis,
                                     LocalCoefficients,
                                     LocalInterpolation>  Traits;

    /** \todo Please doc me !
     */
    PkDGSubsampled2DLocalFiniteElement ()
    { }

    /** \todo Please doc me !
     */
    PkDGSubsampled2DLocalFiniteElement (int variant)
      : coefficients(variant)
    { }

    /** Constructor for six variants with permuted vertices.

        \param vertexmap The permutation of the vertices.  This
        can for instance be generated from the global indices of
        the vertices by reducing those to the integers 0...2
     */
    PkDGSubsampled2DLocalFiniteElement (const unsigned int vertexmap[3])
      : coefficients(vertexmap)
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
    static constexpr GeometryType type ()
    {
      return GeometryTypes::triangle;
    }

  private:
    LocalBasis basis;
    LocalCoefficients coefficients;
    LocalInterpolation interpolation;
  };

}

#endif
