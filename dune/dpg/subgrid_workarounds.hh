// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_SUBGRID_WORKAROUNDS_HH
#define DUNE_DPG_SUBGRID_WORKAROUNDS_HH

#include <type_traits>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/referenceelements.hh>

namespace Dune {

// Forward declarations for dune-subgrid
template<int dim, class HostGrid, bool MapIndexStorage>
class SubGrid;

template<class GridImp>
class SubGridLeafIntersection;

template<class GridImp, class IntersectionImp>
class Intersection;

namespace detail {
  template<class Intersection>
  struct CenterUnitOuterNormal
  {
    static typename Intersection::GlobalCoordinate
    from(const Intersection& intersection) {
      return intersection.centerUnitOuterNormal();
    }
  };

  template<int dim, class HostGrid, bool MapIndexStorage,
           class IntersectionImp>
  struct CenterUnitOuterNormal
    <Intersection<const SubGrid<dim, HostGrid, MapIndexStorage>,
                  IntersectionImp>>
  {
    using SubGridIntersection
        = Intersection<const SubGrid<dim, HostGrid, MapIndexStorage>,
                       IntersectionImp>;
    using ctype = typename SubGridIntersection::ctype;

    static typename SubGridIntersection::GlobalCoordinate
    from(const SubGridIntersection& intersection) {
      return intersection.unitOuterNormal(referenceElement<ctype, dim-1>
                    (intersection.type()).position(0,0));
    }
  };

  template<class Intersection>
  struct GeometryInInside
  {
    static typename Intersection::LocalGeometry
    from(const Intersection& intersection) {
      return intersection.geometryInInside();
    }
  };

  template<int dim, class HostGrid, bool MapIndexStorage,
           class IntersectionImp>
  struct GeometryInInside
    <Intersection<const SubGrid<dim, HostGrid, MapIndexStorage>,
                  IntersectionImp>>
  {
    using SubGridIntersection
        = Intersection<const SubGrid<dim, HostGrid, MapIndexStorage>,
                       IntersectionImp>;
    using ctype = typename SubGridIntersection::ctype;
    using LocalGeometry = AffineGeometry<ctype, dim-1, dim>;

    static LocalGeometry
    from(const SubGridIntersection& intersection) {
      const auto geometry = intersection.geometry();
      const auto innerGeometry = intersection.inside().geometry();
      const auto referenceElement = Dune::referenceElement(intersection);
      const size_t numVertices = referenceElement.size(dim-1);
      std::vector<FieldVector<ctype, dim>> vertices(numVertices);
      for(size_t i = 0; i < numVertices; i++) {
        vertices[i] = innerGeometry.local(
                        geometry.global(
                          referenceElement.position(i, dim-1)));
      }
      return LocalGeometry(referenceElement, vertices);
    }
  };

  template<class Intersection>
  struct Conforming
  {
    static bool from(const Intersection& intersection) {
      return intersection.conforming();
    }
  };

  template<int dim, class HostGrid, bool MapIndexStorage,
           class IntersectionImp>
  struct Conforming
    <Intersection<const SubGrid<dim, HostGrid, MapIndexStorage>,
                  IntersectionImp>>
  {
    using SubGridIntersection
        = Intersection<const SubGrid<dim, HostGrid, MapIndexStorage>,
                       IntersectionImp>;

    static bool
    from(const SubGridIntersection& intersection) {
      // This assumes that neighboring elements on the same level
      // are always conforming which might not always hold.
      return !intersection.neighbor()
        || intersection.inside().level() == intersection.outside().level();
    }
  };

} // end namespace detail

template<class Intersection>
auto centerUnitOuterNormal(const Intersection& intersection) {
  return detail::CenterUnitOuterNormal<std::decay_t<Intersection>>::
      from(intersection);
}

template<class Intersection>
auto geometryInInside(const Intersection& intersection) {
  return detail::GeometryInInside<std::decay_t<Intersection>>::
      from(intersection);
}

template<class Intersection>
bool conforming(const Intersection& intersection) {
  return detail::Conforming<std::decay_t<Intersection>>::from(intersection);
}

} // end namespace Dune

#endif // DUNE_DPG_SUBGRID_WORKAROUNDS_HH
