// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BERNSTEIN_PQK_FACTORY_HH
#define DUNE_BERNSTEIN_PQK_FACTORY_HH

#include <map>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include <dune/localfunctions/bernstein/pk.hh>

namespace Dune
{
  /** \brief Factory to create any kind of Pk/Qk like element wrapped for the virtual interface
   *
   */
  template<class D, class R, int dim, int k>
  struct BernsteinLocalFiniteElementFactory
  {
    typedef typename BernsteinPkLocalFiniteElement<D,R,dim,k>
        ::Traits::LocalBasisType::Traits T;
    typedef LocalFiniteElementVirtualInterface<T> FiniteElementType;
    typedef BernsteinPkLocalFiniteElement<D,R,dim,k> Pk;


    //! create finite element for given GeometryType
    static FiniteElementType* create(const GeometryType& gt)
    {
      if (gt.isSimplex())
        return new LocalFiniteElementVirtualImp<Pk>(Pk());

      DUNE_THROW(NotImplemented,
          "BernsteinLocalFiniteElement is only implemented for Simplices.");
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
  class BernsteinLocalFiniteElementCache
  {
  protected:
    typedef typename BernsteinPkLocalFiniteElement<D,R,dim,k>
        ::Traits::LocalBasisType::Traits T;
    typedef LocalFiniteElementVirtualInterface<T> FE;
    typedef typename std::map<GeometryType,FE*> FEMap;

  public:
    /** \brief Type of the finite elements stored in this cache */
    typedef FE FiniteElementType;

    /** \brief Default constructor */
    BernsteinLocalFiniteElementCache() {}

    /** \brief Copy constructor */
    BernsteinLocalFiniteElementCache(const BernsteinLocalFiniteElementCache& other)
    {
      for(const auto& entry : other.cache_)
        cache_[entry.first] = (entry.second)->clone();
    }

    ~BernsteinLocalFiniteElementCache()
    {
      for(auto&& entry : cache_)
        delete entry.second;
    }

    //! Get local finite element for given GeometryType
    const FiniteElementType& get(const GeometryType& gt) const
    {
      const typename FEMap::const_iterator it = cache_.find(gt);
      if (it==cache_.end())
      {
        FiniteElementType* fe
            = BernsteinLocalFiniteElementFactory<D,R,dim,k>::create(gt);
        if (fe==0)
          DUNE_THROW(NotImplemented,
              "No Pk/Qk like local finite element available for geometry type "
              << gt << " and order " << k);

        cache_[gt] = fe;
        return *fe;
      }
      return *(it->second);
    }

  protected:
    mutable FEMap cache_;

  };

}

#endif
