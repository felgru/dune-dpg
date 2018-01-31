// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BERNSTEIN_PK2DLOCALFINITEELEMENT_HH
#define DUNE_BERNSTEIN_PK2DLOCALFINITEELEMENT_HH

#include <cstddef>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include "pk2d/pk2dlocalbasis.hh"
#include "pk2d/pk2dlocalcoefficients.hh"
#include "pk2d/pk2dlocalinterpolation.hh"

namespace Dune
{
  template<class D, class R, unsigned int k>
  class BernsteinPk2DLocalFiniteElement
  {
  public:
    typedef LocalFiniteElementTraits<BernsteinPk2DLocalBasis<D,R,k>,
        BernsteinPk2DLocalCoefficients<k>,
        BernsteinPk2DLocalInterpolation<BernsteinPk2DLocalBasis<D,R,k>>> Traits;

    BernsteinPk2DLocalFiniteElement ()
    {}

    BernsteinPk2DLocalFiniteElement (int variant) :
      coefficients(variant)
    {}

    /** Constructor for six variants with permuted vertices.

        \param vertexmap The permutation of the vertices.  This
        can for instance be generated from the global indices of
        the vertices by reducing those to the integers 0...2
     */
    BernsteinPk2DLocalFiniteElement (const unsigned int vertexmap[3]) :
      coefficients(vertexmap)
    {}

    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis;
    }

    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients;
    }

    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation;
    }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return basis.size();
    }

    static constexpr GeometryType type ()
    {
      return GeometryTypes::triangle;
    }

  private:
    BernsteinPk2DLocalBasis<D,R,k> basis;
    BernsteinPk2DLocalCoefficients<k> coefficients;
    BernsteinPk2DLocalInterpolation<BernsteinPk2DLocalBasis<D,R,k> > interpolation;
  };
}

#endif
