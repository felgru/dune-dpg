// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_REFINEDFACES_HH
#define DUNE_DPG_REFINEDFACES_HH

#include <dune/dpg/type_traits.hh>
#include <dune/common/version.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/referenceelements.hh>

namespace Dune {

template<class Element>
struct RefinedFaceComputations {
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

  template<class MacroElement>
  RefinedFaceComputations(const Face& face, const Element& subElement,
      const MacroElement& element)
    : RefinedFaceComputations(face, subElement, element,
        element.geometry().global(face.geometry().global({0})),
        element.geometry().global(face.geometry().global({1})))
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

  ElementCoordinate faceToElementPosition(const FaceCoordinate& fc) const {
    return geometryInElement_.global(fc);
  }

  const GeometryInElement& geometryInElement() const {
    return geometryInElement_;
  }

private:
  static GeometryInElement
  getGeometryInElement(const Face& face, const Element& element) {
    const auto geometry = face.geometry();
    const auto innerGeometry = element.geometry();
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
    const auto referenceElement
        = Dune::referenceElement<ctype, facedim>(face.type());
#else
    const auto& referenceElement
        = ReferenceElements<ctype, facedim>::general(face.type());
#endif
    const size_t numVertices = referenceElement.size(facedim);
    std::vector<GlobalCoordinate> vertices(numVertices);
    for(size_t i = 0; i < numVertices; i++) {
      vertices[i] = innerGeometry.local(
                      geometry.global(
                        referenceElement.position(i, facedim)));
    }
    return GeometryInElement(referenceElement, vertices);
  }

  template<class MacroElement>
  RefinedFaceComputations(const Face& face, const Element& subElement,
        const MacroElement& element,
        GlobalCoordinate globalCorner0, GlobalCoordinate globalCorner1)
    : integrationElement_((globalCorner1 - globalCorner0).two_norm())
    /* This won't work for curvilinear elements, but they don't seem
     * to be supported by UG anyway: */
    , unitOuterNormal_{ (globalCorner1[1] - globalCorner0[1])
                      , (globalCorner0[0] - globalCorner1[0]) }
    , geometryInElement_{getGeometryInElement(face, subElement)}
  {
    static_assert(cdim==2, "Computation of unit outer normal for face"
                           " only implemented in 2d!");
    unitOuterNormal_ /= unitOuterNormal_.two_norm();
    if((globalCorner0
        -element.geometry().global(subElement.geometry().center()))
       * unitOuterNormal_ <= 0) {
      unitOuterNormal_ *= -1.;
    }
  }

  const ctype integrationElement_;
  GlobalCoordinate unitOuterNormal_;
  const GeometryInElement geometryInElement_;
};

} // end namespace Dune

#endif // DUNE_DPG_REFINEDFACES_HH
