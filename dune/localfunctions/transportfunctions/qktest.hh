// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:



#ifndef DUNE_Qk_TESTLOCALFINITEELEMENT_HH
#define DUNE_Qk_TESTLOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include "qktest/qktestlocalbasis.hh"
#include "qktest/qktestlocalcoefficients.hh"
#include "qktest/qktestlocalinterpolation.hh"


namespace Dune
{

/** \brief The local QkTest finite element on cubes
      \tparam D Domain data type
      \tparam R Range data type
      \tparam dim Dimension of the simplex
 */
template<class D,class R, int dim,int k>
class QkTestLocalFiniteElement
{
public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<QkTestLocalBasis<D,R,dim,k>,QkTestLocalCoefficients<dim,k>,
            QkTestLocalInterpolation<dim,QkTestLocalBasis<D,R,dim,k>,k >> Traits;


    QkTestLocalFiniteElement (const QkTestLocalFiniteElement & o) : basis(o.basis),coefficients(o.coefficients),interpolation(o.interpolation),gt(o.gt)
    {}

    QkTestLocalFiniteElement (FieldVector<D,dim> transport): basis(transport),coefficients(transport),interpolation(transport)
    {
        gt.makeCube(dim);
        assert(dim==2);
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


    QkTestLocalFiniteElement* clone () const
    {
        return new QkTestLocalFiniteElement(*this);
    }

private:
    QkTestLocalBasis<D,R,dim,k> basis;
    QkTestLocalCoefficients<dim,k> coefficients;
    QkTestLocalInterpolation<dim,QkTestLocalBasis<D,R,dim,k>,k> interpolation;
    GeometryType gt;

};

}

#endif
