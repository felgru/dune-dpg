// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PQKFACE_FACTORY_HH
#define DUNE_PQKFACE_FACTORY_HH

#include <map>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include <dune/localfunctions/lagrange/p0.hh>
#include <dune/localfunctions/lagrange/pkface2d.hh>
#include <dune/localfunctions/lagrange/qkface2d.hh>

namespace Dune
{

  /** \brief Factory that only creates dimension specific local finite elements
   *
   * Empty default implementation
   */
  template<class D, class R, int d, int k>
  struct DimSpecificPQkFaceLocalFiniteElementFactory
  {
    typedef typename P0LocalFiniteElement<D,R,d>::Traits::LocalBasisType::Traits T;

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
  struct DimSpecificPQkFaceLocalFiniteElementFactory<D,R,3,k>
  {
    typedef typename P0LocalFiniteElement<D,R,3>::Traits::
                                    LocalBasisType::Traits T;

    //! create finite element for given GeometryType
    static LocalFiniteElementVirtualInterface<T>* create(const GeometryType& gt)
    {
      DUNE_THROW(Dune::NotImplemented, "pqkFace not implemented for 3d");
      return 0;
    }
  };

  /** \brief Factory that only creates dimension specific local finite elements
   *
   * Specialization for dim=2
   */
  template<class D, class R, int k>
  struct DimSpecificPQkFaceLocalFiniteElementFactory<D,R,2,k>
  {
    typedef typename P0LocalFiniteElement<D,R,2>::Traits::
                                    LocalBasisType::Traits T;
    typedef PkFace2DLocalFiniteElement<D,R,k> Pk;
    typedef QkFace2DLocalFiniteElement<D,R,k> Qk;

    //! create finite element for given GeometryType
    static LocalFiniteElementVirtualInterface<T>* create(const GeometryType& gt)
    {
      if (gt.isSimplex())
        return new LocalFiniteElementVirtualImp<Pk>(Pk());
      if (gt.isCube())
        return new LocalFiniteElementVirtualImp<Qk>(Qk());
      return 0;
    }
  };

  /** \brief Factory to create any kind of Pk/Qk like element wrapped for the virtual interface
   *
   */
  template<class D, class R, int dim, int k>
  struct PQkFaceLocalFiniteElementFactory
  {
    typedef typename P0LocalFiniteElement<D,R,dim>::Traits::
                                    LocalBasisType::Traits T;
    typedef LocalFiniteElementVirtualInterface<T> FiniteElementType;
    typedef P0LocalFiniteElement<D,R,dim> P0;

    //! create finite element for given GeometryType
    static FiniteElementType* create(const GeometryType& gt)
    {
      if (k==0)
      {
        DUNE_THROW(Dune::NotImplemented, "k=0 does not make sense for faces");
        return new LocalFiniteElementVirtualImp<P0>(P0(gt));
      }

      return DimSpecificPQkFaceLocalFiniteElementFactory<D,R,dim,k>::create(gt);
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
  class PQkFaceLocalFiniteElementCache
  {
  protected:
    typedef typename P0LocalFiniteElement<D,R,dim>::Traits::
                                    LocalBasisType::Traits T;
    typedef LocalFiniteElementVirtualInterface<T> FE;
    typedef typename std::map<GeometryType,FE*> FEMap;

  public:
    /** \brief Type of the finite elements stored in this cache */
    typedef FE FiniteElementType;

    /** \brief Default constructor */
    PQkFaceLocalFiniteElementCache() {}

    /** \brief Copy constructor */
    PQkFaceLocalFiniteElementCache(const PQkFaceLocalFiniteElementCache& other)
    {
      typename FEMap::iterator it = other.cache_.begin();
      typename FEMap::iterator end = other.cache_.end();
      for(; it!=end; ++it)
        cache_[it->first] = (it->second)->clone();
    }

    ~PQkFaceLocalFiniteElementCache()
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
        FiniteElementType* fe = PQkFaceLocalFiniteElementFactory<D,R,dim,k>::create(gt);
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
