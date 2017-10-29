/*
 * pqktransportfactory.hh
 *
 *  Created on: May 12, 2015
 *      Author: koenig
 */

#ifndef DUNE_PQKTRANSPORTFACTORY_HH_
#define DUNE_PQKTRANSPORTFACTORY_HH_

// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <map>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include <dune/localfunctions/lagrange/p0.hh>
#include <dune/localfunctions/lagrange/pk.hh>
#include <dune/localfunctions/lagrange/q1.hh>
#include <dune/localfunctions/lagrange/qk.hh>
#include <dune/localfunctions/lagrange/prismp1.hh>
#include <dune/localfunctions/lagrange/prismp2.hh>
#include <dune/localfunctions/lagrange/pyramidp1.hh>
#include <dune/localfunctions/lagrange/pyramidp2.hh>
#include <dune/localfunctions/transportfunctions/qktransport.hh>

namespace Dune
{

/** \brief Factory that only creates dimension specific local finite elements
 *
 * Empty default implementation
 */
template<class D, class R, int d, int k>
struct DimSpecificPQkTransportLocalFiniteElementFactory
{
    typedef typename
      P0LocalFiniteElement<D,R,d>::Traits::LocalBasisType::Traits T;

    //! create finite element for given GeometryType
    static LocalFiniteElementVirtualInterface<T>* create(const GeometryType& gt)
            {
        return 0;
            }
};

/** \brief Factory that only creates dimension specific local finite elements
 *
 * Specialization for dim=3
 */
template<class D, class R, int k>
struct DimSpecificPQkTransportLocalFiniteElementFactory<D,R,3,k>
{
    typedef typename
      P0LocalFiniteElement<D,R,3>::Traits::LocalBasisType::Traits T;
    typedef PrismP1LocalFiniteElement<D,R> PrismP1;
    typedef PrismP2LocalFiniteElement<D,R> PrismP2;
    typedef PyramidP1LocalFiniteElement<D,R> PyramidP1;
    typedef PyramidP2LocalFiniteElement<D,R> PyramidP2;

    //! create finite element for given GeometryType
    static LocalFiniteElementVirtualInterface<T>* create(const GeometryType& gt)
            {
        if ((gt.isPrism())and (k==1))
            return new LocalFiniteElementVirtualImp<PrismP1>(PrismP1());
        if ((gt.isPrism())and (k==2))
            return new LocalFiniteElementVirtualImp<PrismP2>(PrismP2());
        if ((gt.isPyramid())and (k==1))
            return new LocalFiniteElementVirtualImp<PyramidP1>(PyramidP1());
        if ((gt.isPyramid())and (k==2))
            return new LocalFiniteElementVirtualImp<PyramidP2>(PyramidP2());
        return 0;
            }
};


/** \brief Factory to create any kind of QkTransport in 2D like element wrapped for the virtual interface
 *
 */
template<class D, class R, int dim, int k>
struct PQkTransportLocalFiniteElementFactory
{
    typedef typename
      P0LocalFiniteElement<D,R,dim>::Traits::LocalBasisType::Traits T;
    typedef LocalFiniteElementVirtualInterface<T> FiniteElementType;
    typedef P0LocalFiniteElement<D,R,dim> P0;
    typedef PkLocalFiniteElement<D,R,dim,k> Pk;
    typedef QkLocalFiniteElement<D,R,dim,k> Qk;
    typedef QkTransportLocalFiniteElement <D,R,dim,k> QTransport;


    //! create finite element for given GeometryType
    static FiniteElementType* create(const GeometryType& gt, FieldVector<D,dim> transport)
    {
        if (k==0)
            return new LocalFiniteElementVirtualImp<P0>(P0(gt));

        if (gt.isSimplex())
            return new LocalFiniteElementVirtualImp<Pk>(Pk());

        //NEW ELEMENT with transport direction
        if (gt.isCube())
            return new LocalFiniteElementVirtualImp<QTransport>(QTransport(transport));

        return DimSpecificPQkTransportLocalFiniteElementFactory<D,R,dim,k>::create(gt);
    }
};



/** \brief A cache that stores all available Pk/Qk like local finite elements for the given dimension and order
 *
 * An interface for dealing with different vertex orders is currently missing.
 * So you can in general only use this for order=1,2 or with global DG spaces
 *
 * \tparam D Type used for domain coordinates
 * \tparam R Type used for shape function values
 * \tparam dim Element dimension
 * \tparam k Element order
 */
template<class D, class R, int dim, int k>
class PQkTransportLocalFiniteElementCache
{
protected:
    typedef typename
      P0LocalFiniteElement<D,R,dim>::Traits::LocalBasisType::Traits T;
    typedef LocalFiniteElementVirtualInterface<T> FE;
    typedef typename std::map<GeometryType,FE*> FEMap;

public:
    /** \brief Type of the finite elements stored in this cache */
    typedef FE FiniteElementType;


    /** \brief Default constructor */
    PQkTransportLocalFiniteElementCache() {}


    /** \brief Copy constructor */
    PQkTransportLocalFiniteElementCache(const PQkTransportLocalFiniteElementCache& other)
    {
        typename FEMap::iterator it = other.cache_.begin();
        typename FEMap::iterator end = other.cache_.end();
        for(; it!=end; ++it)
            cache_[it->first] = (it->second)->clone();

    }

    ~PQkTransportLocalFiniteElementCache()
    {
        typename FEMap::iterator it = cache_.begin();
        typename FEMap::iterator end = cache_.end();
        for(; it!=end; ++it)
            delete it->second;
    }


    //! Get local finite element for given GeometryType
    const FiniteElementType& get(const GeometryType& gt, FieldVector<double,dim> transport) const
    {
        typename FEMap::const_iterator it = cache_.find(gt);
        if (it==cache_.end())
        {
            FiniteElementType* fe = PQkTransportLocalFiniteElementFactory<D,R,dim,k>::create(gt,transport);
            if (fe==0)
                DUNE_THROW(Dune::NotImplemented,"No Pk/Qk like local finite element available for geometry type " << gt << " and order " << k);

            cache_[gt] = fe;
            return *fe;
        }
        return *(it->second);
    }

protected:
    mutable FEMap cache_;

};

}







#endif /* DUNE_PQKTRANSPORTFACTORY_HH_ */
