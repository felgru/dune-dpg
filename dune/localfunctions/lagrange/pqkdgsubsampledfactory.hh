// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PQKDGSUBSAMPLED_FACTORY_HH
#define DUNE_PQKDGSUBSAMPLED_FACTORY_HH

#include <map>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include <dune/localfunctions/lagrange/pkdgsubsampled2d.hh>

namespace Dune
{

  /** \brief Factory that only creates dimension specific local finite elements
   *
   * Empty default implementation
   */
  template<class D, class R, int dim, int s, int k>
  struct DimSpecificPQkDGSubsampledLocalFiniteElementFactory
  { };

  /** \brief Factory that only creates dimension specific local finite elements
   *
   * Specialization for dim=2
   */
  template<class D, class R, int s, int k>
  struct DimSpecificPQkDGSubsampledLocalFiniteElementFactory<D,R,2,s,k>
  {
    typedef typename PkDGSubsampled2DLocalFiniteElement<D,R,s,0>
        ::Traits::LocalBasisType::Traits T;
    typedef PkDGSubsampled2DLocalFiniteElement<D,R,s,k> Pk;
    // typedef QkDGSubsampled2DLocalFiniteElement<D,R,s,k> Qk;

    //! create finite element for given GeometryType
    static LocalFiniteElementVirtualInterface<T>* create(const GeometryType& gt)
    {
      if (gt.isSimplex())
        return new LocalFiniteElementVirtualImp<Pk>(Pk());
      // if (gt.isCube())
      //   return new LocalFiniteElementVirtualImp<Qk>(Qk());
      return nullptr;
    }
  };

  /** \brief Factory to create any kind of Pk/Qk like element wrapped for the virtual interface
   *
   */
  template<class D, class R, int dim, int s, int k>
  struct PQkDGSubsampledLocalFiniteElementFactory
  {
    typedef typename PkDGSubsampled2DLocalFiniteElement<D,R,s,0>
        ::Traits::LocalBasisType::Traits T;
    typedef LocalFiniteElementVirtualInterface<T> FiniteElementType;

    //! create finite element for given GeometryType
    static FiniteElementType* create(const GeometryType& gt)
    {
      return DimSpecificPQkDGSubsampledLocalFiniteElementFactory
               <D,R,dim,s,k>::create(gt);
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
  class PQkDGSubsampledLocalFiniteElementCache
  {
  protected:
    typedef typename PkDGSubsampled2DLocalFiniteElement<D,R,s,0>
        ::Traits::LocalBasisType::Traits T;
    typedef LocalFiniteElementVirtualInterface<T> FE;
    typedef typename std::map<GeometryType,FE*> FEMap;

  public:
    /** \brief Type of the finite elements stored in this cache */
    typedef FE FiniteElementType;

    /** \brief Default constructor */
    PQkDGSubsampledLocalFiniteElementCache() {}

    /** \brief Copy constructor */
    PQkDGSubsampledLocalFiniteElementCache
        (const PQkDGSubsampledLocalFiniteElementCache& other)
    {
      for(const auto& [gt, fe] : other.cache_)
        cache_[gt] = fe->clone();
    }

    ~PQkDGSubsampledLocalFiniteElementCache()
    {
      for(auto&& [gt, fe] : cache_)
        delete fe;
    }

    //! Get local finite element for given GeometryType
    const FiniteElementType& get(const GeometryType& gt) const
    {
      typename FEMap::const_iterator it = cache_.find(gt);
      if (it==cache_.end())
      {
        FiniteElementType* fe =
            PQkDGSubsampledLocalFiniteElementFactory<D,R,dim,s,k>
            ::create(gt);
        if (fe==nullptr)
          DUNE_THROW(Dune::NotImplemented,
            "No Pk/Qk like local finite element available for geometry type "
            << gt << ", " << s << "subsamples and order " << k);

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
