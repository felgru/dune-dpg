// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_REFERENCE_REFINEMENTFACTORY_HH
#define DUNE_REFERENCE_REFINEMENTFACTORY_HH

#include <map>
#include <utility>

#include <dune/dpg/common/memory.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/grid/uggrid/uggridfactory.hh>

namespace Dune {
namespace Functions {

  /**
   * \brief Factory to create any kind of refined reference cell
   */
  template<class D, int dim, int level>
  struct ReferenceRefinementFactory
  {
    typedef UGGrid<dim> GridType;


    //! create finite element for given GeometryType
    static std::unique_ptr<GridType> create(const GeometryType& gt)
    {
      const auto& referenceElement
        = ReferenceElements<D, dim>::general(gt);
      auto referenceGeometry = referenceElement.template geometry<0>(0);

      int numVertices = referenceGeometry.corners();
      std::vector<unsigned int> vertices(numVertices);
      GridFactory<GridType> gridFactory;
      for(int i=0; i<numVertices; ++i)
      {
        gridFactory.insertVertex(referenceGeometry.corner(i));
        vertices[i]=i;
      }
      gridFactory.insertElement(gt, vertices);

      std::unique_ptr<GridType> referenceGrid{gridFactory.createGrid()};
      referenceGrid->globalRefine(level);

      return referenceGrid;
    }
  };



  /** \brief A cache that stores all available refinements of
   *         reference cells for the given dimension and refinement level
   *
   * \tparam D Type used for domain coordinates
   * \tparam dim Element dimension
   * \tparam level Refinement level
   */
  template<class D, int dim, int level>
  class ReferenceRefinementCache
  {
  protected:
    typedef UGGrid<dim> Grid;
    typedef typename std::map<GeometryType,std::unique_ptr<Grid>> RefinementMap;

  public:
    /** \brief Type of the grid stored in this cache */
    typedef Grid GridType;

    /** \brief Default constructor */
    ReferenceRefinementCache() {}

    /** \brief Copy constructor does nothing as, we cannot copy the grids
     *         used for the refinements */
    ReferenceRefinementCache(const ReferenceRefinementCache& other)
      : cache_{}
    {};

    //! Get refined reference cell for given GeometryType
    const GridType& get(const GeometryType& gt) const
    {
      typename RefinementMap::const_iterator it = cache_.find(gt);
      if (it==cache_.end())
      {
        std::unique_ptr<GridType> refinement
          = ReferenceRefinementFactory<D,dim,level>::create(gt);
        if (refinement==nullptr)
          DUNE_THROW(Dune::NotImplemented,"No refined reference element available for geometry type " << gt << " and refinement level " << level);

        auto res = cache_.emplace(gt, std::move(refinement));
        return *(res.first->second);
      }
      return *(it->second);
    }

  protected:
    mutable RefinementMap cache_;

  };

}} // end namespace Dune::Functions

#endif
