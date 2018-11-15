// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_FACES_HH
#define DUNE_DPG_FACES_HH

#include <dune/dpg/assemble_types.hh>
#include <dune/dpg/type_traits.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/referenceelements.hh>

namespace Dune {

template<class Element>
struct FaceComputations {
  using ElementGeometry = typename Element::Geometry;
  using Face = typename Element::template Codim<1>::Entity;
  using ctype = typename ElementGeometry::ctype;
  using GlobalCoordinate = typename ElementGeometry::GlobalCoordinate;
  using ElementCoordinate = typename ElementGeometry::LocalCoordinate;
  using FaceCoordinate = typename Face::Geometry::LocalCoordinate;

  constexpr static size_t cdim = tuple_size<GlobalCoordinate>::value;
  constexpr static size_t facedim = tuple_size<FaceCoordinate>::value;
  constexpr static size_t elemdim = tuple_size<ElementCoordinate>::value;

  using GeometryInElement = AffineGeometry<ctype, facedim, elemdim>;

  FaceComputations(const Face& face, const Element& element)
    : FaceComputations(face, element,
        face.geometry().global({0}), face.geometry().global({1}))
  {
    static_assert(cdim==2, "Computation of unit outer normal for face"
                           " only implemented in 2d!");
  }

  ctype integrationElement() const {
    return integrationElement_;
  }

  const GlobalCoordinate& unitOuterNormal() const {
    return unitOuterNormal_;
  }

  GlobalCoordinate integrationOuterNormal() const {
    GlobalCoordinate res = unitOuterNormal_;
    res *= integrationElement_;
    return res;
  }

  int unitOuterNormalSign() const {
    for (const auto& component : unitOuterNormal_)
    {
      if (component < -1e-10)
      {
        return -1;
      }
      else if (component > 1e-10)
      {
        return 1;
      }
    }
    return 1;
  }

  ElementCoordinate faceToElementPosition(const FaceCoordinate& fc) const {
    return geometryInElement_.global(fc);
  }

  const GeometryInElement& geometryInElement() const {
    return geometryInElement_;
  }

  template<IntegrationType type>
  bool skipFace(const FieldVector<double, cdim>& direction) const
  {
    if(type == IntegrationType::travelDistanceWeighted)
      /* Only integrate over inflow boundaries. */
      return direction * unitOuterNormal() >= 0;
    else return false;
  }

private:
  static GeometryInElement
  getGeometryInElement(const Face& face, const Element& element) {
    const auto geometry = face.geometry();
    const auto innerGeometry = element.geometry();
    const auto referenceElement
        = Dune::referenceElement<ctype, facedim>(face.type());
    const size_t numVertices = referenceElement.size(facedim);
    std::vector<GlobalCoordinate> vertices(numVertices);
    for(size_t i = 0; i < numVertices; i++) {
      vertices[i] = innerGeometry.local(
                      geometry.global(
                        referenceElement.position(i, facedim)));
    }
    return GeometryInElement(referenceElement, vertices);
  }

  FaceComputations(const Face& face, const Element& element,
        GlobalCoordinate globalCorner0, GlobalCoordinate globalCorner1)
    : integrationElement_((globalCorner1 - globalCorner0).two_norm())
    , geometryInElement_{getGeometryInElement(face, element)}
  {
    static_assert(cdim==2, "Computation of unit outer normal for face"
                           " only implemented in 2d!");
    /* This won't work for curvilinear elements, but they don't seem
     * to be supported by UG anyway. */
    unitOuterNormal_
      = { (globalCorner1[1] - globalCorner0[1])
        , (globalCorner0[0] - globalCorner1[0]) };
    unitOuterNormal_ /= unitOuterNormal_.two_norm();

    if((globalCorner0-element.geometry().center())* unitOuterNormal_ < 0) {
      unitOuterNormal_ *= -1;
    }
  }

  const ctype integrationElement_;
  GlobalCoordinate unitOuterNormal_;
  const GeometryInElement geometryInElement_;
};

} // end namespace Dune

#endif // DUNE_DPG_FACES_HH
