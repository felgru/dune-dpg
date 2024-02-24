// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_REFINEDFACES_HH
#define DUNE_DPG_REFINEDFACES_HH

#include <dune/dpg/assemble_types.hh>
#include <dune/dpg/traveldistancenorm.hh>
#include <dune/dpg/type_traits.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/rangegenerators.hh>

namespace Dune {

template<IntegrationType type>
struct RefinedFaceIntegrationData {
  constexpr static int dim = 2;

  template <class Geometry, class SubGeometry>
  RefinedFaceIntegrationData(
      const Geometry&,
      const SubGeometry&,
      const FieldVector<double, dim>&) {}

  template<class LhsSpace, class RhsSpace, class Face, class GeometryInElement>
  QuadratureRule<double, 1> quadratureRule(
      const Face& face,
      unsigned int quadratureOrder,
      unsigned int,
      const GeometryInElement&) const
  {
    QuadratureRule<double, 1> quadFace
      = detail::ChooseQuadrature<LhsSpace, RhsSpace, Face>
        ::Quadrature(face, quadratureOrder);
    return quadFace;
  }

  template<class ElementCoordinate>
  double extraIntegrationWeight(const ElementCoordinate&) const noexcept
  {
    return 1.;
  }
};

template<>
struct RefinedFaceIntegrationData<IntegrationType::travelDistanceWeighted>
{
  constexpr static int dim = 2;

  template <class Geometry, class SubGeometry>
  RefinedFaceIntegrationData(
      const Geometry& geometry,
      const SubGeometry& subGeometryInReferenceElement,
      const FieldVector<double, dim>& direction)
    : referenceBeta(detail::referenceSubElementBeta(geometry,
                    subGeometryInReferenceElement, direction))
  {}

  template<class LhsSpace, class RhsSpace, class Face, class GeometryInElement>
  QuadratureRule<double, 1> quadratureRule(
      const Face& face,
      unsigned int quadratureOrder,
      unsigned int nOutflowFaces,
      const GeometryInElement& geometryInElement)
  {
    QuadratureRule<double, 1> quadFace
      = detail::ChooseQuadrature<LhsSpace, RhsSpace, Face>
        ::Quadrature(face, quadratureOrder);
    if (nOutflowFaces > 1) {
      quadFace = SplitQuadratureRule<double>(
          quadFace,
          detail::splitPointOfInflowFaceInTriangle(
              geometryInElement, referenceBeta));
    }
    return quadFace;
  }

  template<class ElementCoordinate>
  double extraIntegrationWeight(const ElementCoordinate& elementQuadPosSubCell)
  noexcept
  {
    // factor r_K(s)/|beta|
    return detail::travelDistance(elementQuadPosSubCell, referenceBeta);
  }

  private:
    FieldVector<double,dim> referenceBeta;
};

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

  ctype integrationElement() const noexcept {
    return integrationElement_;
  }

  const GlobalCoordinate& unitOuterNormal() const noexcept {
    return unitOuterNormal_;
  }

  int unitOuterNormalSign() const noexcept {
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

  ElementCoordinate
  faceToElementPosition(const FaceCoordinate& fc) const noexcept {
    return geometryInElement_.global(fc);
  }

  const GeometryInElement& geometryInElement() const noexcept {
    return geometryInElement_;
  }

  template<IntegrationType type, class LhsSpace, class RhsSpace>
  QuadratureRule<double, 1> quadratureRule(
      const Face& face,
      unsigned int quadratureOrder,
      unsigned int nOutflowFaces,
      RefinedFaceIntegrationData<type>& integrationData) const
  {
    return integrationData.template quadratureRule<LhsSpace, RhsSpace>
        (face, quadratureOrder, nOutflowFaces, geometryInElement());
  }

  template<IntegrationType type>
  bool skipFace(const FieldVector<double, cdim>& direction) const noexcept
  {
    if constexpr (type == IntegrationType::travelDistanceWeighted)
      /* Only integrate over inflow boundaries. */
      return direction * unitOuterNormal() >= 0;
    else return false;
  }

  template<IntegrationType type, class LocalCoefficients>
  double integrationWeight(
      const LocalCoefficients& localCoefficients,
      const ElementCoordinate elementQuadPos,
      const ElementCoordinate elementQuadPosSubCell,
      const FieldVector<double, cdim> direction,
      const double quadWeight,
      RefinedFaceIntegrationData<type>& integrationData) const
  {
    double integrationWeight;
    if constexpr (type == IntegrationType::normalVector ||
                  type == IntegrationType::travelDistanceWeighted) {
      integrationWeight = localCoefficients.localFactor()(elementQuadPos)
                        * quadWeight
                        * integrationElement_;
      // TODO: scale direction to length 1
      if constexpr (type == IntegrationType::travelDistanceWeighted)
        integrationWeight *= std::fabs(direction * unitOuterNormal_);
      else
        integrationWeight *= direction * unitOuterNormal_;
    } else if constexpr (type == IntegrationType::normalSign) {
      const int sign = unitOuterNormalSign();

      integrationWeight = sign
                        * localCoefficients.localFactor()(elementQuadPos)
                        * quadWeight * integrationElement_;
    }
    integrationWeight *=
        integrationData.extraIntegrationWeight(elementQuadPosSubCell);

    return integrationWeight;
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

template<class SubElement, class Element>
unsigned int outflowFacesOfSubElement
    (SubElement subElement,
     Element element,
     typename Element::Geometry::GlobalCoordinate direction)
{
  unsigned int nOutflowFaces = 0;
  for (const auto& face : subEntities(subElement, Codim<1>{}))
  {
    /* This won't work for curvilinear elements, but they don't seem
     * to be supported by UG anyway. */
    const auto unitOuterNormal
        = RefinedFaceComputations<SubElement>(face, subElement, element)
            .unitOuterNormal();

    const double prod = direction * unitOuterNormal;
    if(prod > 0)
      ++nOutflowFaces;
  }
  return nOutflowFaces;
}

} // end namespace Dune

#endif // DUNE_DPG_REFINEDFACES_HH
