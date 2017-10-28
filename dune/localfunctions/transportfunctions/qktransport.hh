// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:



#ifndef DUNE_QK_TRANSPORTLOCALFINITEELEMENT_HH
#define DUNE_QK_TRANSPORTLOCALFINITEELEMENT_HH

#include <dune/common/version.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include "qktransport/qktransportlocalbasis.hh"
#include "qktransport/qktransportlocalcoefficients.hh"
#include "qktransport/qktransportlocalinterpolation.hh"


namespace Dune
{

/** \brief The local QkTransport finite element on cubes
      \tparam D Domain data type
      \tparam R Range data type
      \tparam dim Dimension of the simplex
 */
template<class D,class R, int dim,int k>
class QkTransportLocalFiniteElement
{
public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<QkTransportLocalBasis<D,R,dim,k>,QkTransportLocalCoefficients<dim,k>,
            QkTransportLocalInterpolation<dim,QkTransportLocalBasis<D,R,dim,k>,k >> Traits;


    QkTransportLocalFiniteElement (const QkTransportLocalFiniteElement & o) = default;

    QkTransportLocalFiniteElement (FieldVector<D,dim> transport)
      : basis(transport)
      , coefficients(transport)
      , interpolation(transport)
#if DUNE_VERSION_NEWER(DUNE_GRID,2,6)
    {
        static_assert(dim==2, "QkTransportLocalFiniteElement is only implemented in 2D.");
    }
#else
    {
        gt.makeCube(dim);
        static_assert(dim==2, "QkTransportLocalFiniteElement is only implemented in 2D.");
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
      return GeometryTypes::cube(dim);
    }
#else
    GeometryType type () const
    {
      return gt;
    }
#endif

private:
    QkTransportLocalBasis<D,R,dim,k> basis;
    QkTransportLocalCoefficients<dim,k> coefficients;
    QkTransportLocalInterpolation<dim,QkTransportLocalBasis<D,R,dim,k>,k>
                                        interpolation;
#if not(DUNE_VERSION_NEWER(DUNE_GRID,2,6))
    GeometryType gt;
#endif

};

}

#endif
