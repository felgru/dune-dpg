// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PQKSUBSAMPLED_FACTORY_HH
#define DUNE_PQKSUBSAMPLED_FACTORY_HH

#include <map>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include <dune/localfunctions/lagrange/qksubsampled.hh>

namespace Dune
{

  /** \brief Factory to create any kind of Pk/Qk like element wrapped for the virtual interface
   *
   */
  template<class D, class R, int dim, int s, int k>
  struct PQkSubsampledLocalFiniteElementFactory
  {
    typedef typename QkSubsampledLocalFiniteElement<D,R,dim,s,0>::Traits::LocalBasisType::Traits T;
    typedef LocalFiniteElementVirtualInterface<T> FiniteElementType;
    // typedef PkLocalFiniteElement<D,R,dim,s,k> Pk;
    typedef QkSubsampledLocalFiniteElement<D,R,dim,s,k> Qk;


    //! create finite element for given GeometryType
    static FiniteElementType* create(const GeometryType& gt)
    {
      if (gt.isSimplex())
        // return new LocalFiniteElementVirtualImp<Pk>(Pk());
        DUNE_THROW(Dune::NotImplemented,"No Pk/Qk like local finite element available for geometry type " << gt << ".");

      if (gt.isCube())
        return new LocalFiniteElementVirtualImp<Qk>(Qk());

      DUNE_THROW(Dune::NotImplemented,"No Pk/Qk like local finite element available for geometry type " << gt << ".");
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
   * \tparam s number of subsamples
   * \tparam k Element order
   */
  template<class D, class R, int dim, int s, int k>
  class PQkSubsampledLocalFiniteElementCache
  {
  protected:
    typedef typename QkSubsampledLocalFiniteElement<D,R,dim,s,0>::Traits::
                                                    LocalBasisType::Traits T;
    typedef LocalFiniteElementVirtualInterface<T> FE;
    typedef typename std::map<GeometryType,FE*> FEMap;

  public:
    /** \brief Type of the finite elements stored in this cache */
    typedef FE FiniteElementType;

    /** \brief Default constructor */
    PQkSubsampledLocalFiniteElementCache() {}

    /** \brief Copy constructor */
    PQkSubsampledLocalFiniteElementCache(
            const PQkSubsampledLocalFiniteElementCache& other)
    {
      typename FEMap::iterator it = other.cache_.begin();
      typename FEMap::iterator end = other.cache_.end();
      for(; it!=end; ++it)
        cache_[it->first] = (it->second)->clone();
    }

    ~PQkSubsampledLocalFiniteElementCache()
    {
      typename FEMap::iterator it = cache_.begin();
      typename FEMap::iterator end = cache_.end();
      for(; it!=end; ++it)
        delete it->second;
    }

    //! Get local finite element for given GeometryType
    const FiniteElementType& get(const GeometryType& gt) const
    {
      typename FEMap::const_iterator it = cache_.find(gt);
      if (it==cache_.end())
      {
        FiniteElementType* fe =
            PQkSubsampledLocalFiniteElementFactory<D,R,dim,s,k>::create(gt);
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

#endif
