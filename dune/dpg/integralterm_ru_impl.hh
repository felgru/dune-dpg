#include "refinedfaces.hh"

namespace Dune {
namespace detail {

template <IntegrationType type,
          class LhsSpace,
          class RhsSpace>
struct GetLocalMatrix<type, LhsSpace, RhsSpace, true, false>
{
using LhsLocalView = typename LhsSpace::LocalView;
using RhsLocalView = typename RhsSpace::LocalView;

template <class MatrixType,
          class Element,
          class LocalCoefficients>
inline static void interiorImpl(LhsLocalView& lhsLocalView,
                                const RhsLocalView& rhsLocalView,
                                MatrixType& elementMatrix,
                                size_t lhsSpaceOffset,
                                size_t rhsSpaceOffset,
                                unsigned int quadratureOrder,
                                const Element& element,
                                const LocalCoefficients& localCoefficients)
{
  constexpr int dim = Element::mydimension;
  const auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& rhsLocalFiniteElement = rhsLocalView.tree().finiteElement();

  const unsigned int nRhs(rhsLocalFiniteElement.localBasis().size());

  typename detail::ChooseQuadrature<LhsSpace, RhsSpace, Element>::type quad
    = detail::ChooseQuadrature<LhsSpace, RhsSpace, Element>
      ::Quadrature(element, quadratureOrder);

  const auto referenceGridView =
      lhsLocalView.tree().refinedReferenceElementGridView();

  unsigned int subElementOffset = 0;
  unsigned int subElementIndex = 0;
  lhsLocalView.resetSubElements();
  for(const auto& subElement : elements(referenceGridView)) {
    lhsLocalView.bindSubElement(subElement);
    const auto& lhsLocalFiniteElement = lhsLocalView.tree().finiteElement();
    const unsigned int nLhs(lhsLocalFiniteElement.localBasis().size());

    const auto subGeometryInReferenceElement = subElement.geometry();
    for (const auto& quadPoint : quad) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quadPoint.position();
      // Global position of the current quadrature point
      const FieldVector<double,dim> elementQuadPos
          = subGeometryInReferenceElement.global(quadPos);

      // The multiplicative factor in the integral transformation formula
      const double integrationWeight
        = geometry.integrationElement(elementQuadPos)
        * subGeometryInReferenceElement.integrationElement(quadPos)
        * quadPoint.weight() * localCoefficients.localFactor()(elementQuadPos);

      ///////////////////////////////////////
      // evaluate finite element functions //
      ///////////////////////////////////////

      using LhsFunctionEvaluator
        = detail::LocalRefinedFunctionEvaluation<dim, type>;

      const std::vector<FieldVector<double,1> > lhsValues =
          LhsFunctionEvaluator::template evaluateLhs
                  <is_ContinuouslyRefinedFiniteElement<LhsSpace>::value>
                        (lhsLocalFiniteElement,
                         subElementIndex,
                         quadPos,
                         geometry,
                         subGeometryInReferenceElement,
                         localCoefficients);

      using RhsFunctionEvaluator = detail::LocalFunctionEvaluation<dim, type>;

      const std::vector<FieldVector<double,1> > rhsValues =
          RhsFunctionEvaluator::evaluateRhs
                        (rhsLocalFiniteElement,
                         elementQuadPos,
                         geometry,
                         localCoefficients);

      // Compute the actual matrix entries
      for (unsigned int i=0; i<nLhs; i++)
      {
        for (unsigned int j=0; j<nRhs; j++)
        {
          elementMatrix[i+lhsSpaceOffset+subElementOffset]
                       [j+rhsSpaceOffset]
                  += (lhsValues[i] * rhsValues[j]) * integrationWeight;
        }
      }
    }
    if(is_DGRefinedFiniteElement<LhsSpace>::value)
      subElementOffset += lhsLocalFiniteElement.size();
    subElementIndex++;
  }
}


template <class MatrixType,
          class Element,
          class LocalCoefficients>
inline static void
faceImpl(LhsLocalView& lhsLocalView,
         const RhsLocalView& rhsLocalView,
         MatrixType& elementMatrix,
         size_t lhsSpaceOffset,
         size_t rhsSpaceOffset,
         unsigned int quadratureOrder,
         const Element& element,
         const LocalCoefficients& localCoefficients)
{
  constexpr int dim = Element::mydimension;
  const auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& rhsLocalFiniteElement = rhsLocalView.tree().finiteElement();

  const unsigned int nRhs(rhsLocalFiniteElement.localBasis().size());

  const auto referenceGridView =
      lhsLocalView.tree().refinedReferenceElementGridView();

  const auto direction = localCoefficients.localDirection()({0.5,0.5});

  unsigned int subElementOffset = 0;
  unsigned int subElementIndex = 0;
  lhsLocalView.resetSubElements();
  for(const auto& subElement : elements(referenceGridView))
  {
    using SubElement = std::decay_t<decltype(subElement)>;

    lhsLocalView.bindSubElement(subElement);
    const auto& lhsLocalFiniteElement = lhsLocalView.tree().finiteElement();
    const unsigned int nLhs(lhsLocalFiniteElement.localBasis().size());

    const auto subGeometryInReferenceElement = subElement.geometry();

    unsigned int nOutflowFaces = 0;
    for (unsigned short f = 0, fMax = subElement.subEntities(1); f < fMax; f++)
    {
      const auto face = subElement.template subEntity<1>(f);
      const FieldVector<double,dim> unitOuterNormal
          = RefinedFaceComputations<SubElement>(face, subElement, element)
              .unitOuterNormal();

      const double prod = direction * unitOuterNormal;
      if(prod > 0)
        ++nOutflowFaces;
    }

    FieldVector<double,dim> referenceBeta
        = detail::referenceBeta(geometry,
            subGeometryInReferenceElement, direction);

    for (unsigned short f = 0, fMax = subElement.subEntities(1); f < fMax; f++)
    {
      const auto face = subElement.template subEntity<1>(f);
      const auto faceComputations
          = RefinedFaceComputations<SubElement>(face, subElement, element);
      if(faceComputations.template skipFace<type>(direction)) continue;

      using Face = std::decay_t<decltype(face)>;

      const double integrationElement = faceComputations.integrationElement();

      const FieldVector<double,dim> unitOuterNormal
          = faceComputations.unitOuterNormal();

      QuadratureRule<double, 1> quadFace
        = detail::ChooseQuadrature<LhsSpace, RhsSpace, Face>
          ::Quadrature(face, quadratureOrder);
      if (type == IntegrationType::travelDistanceWeighted &&
          nOutflowFaces > 1) {
        quadFace = SplitQuadratureRule<double>(
            quadFace,
            detail::splitPointOfInflowFaceInTriangle(
                faceComputations.geometryInElement(), referenceBeta));
      }

      for (const auto& quadPoint : quadFace) {

        // Position of the current quadrature point in the reference element
        // (face!)
        const FieldVector<double,dim-1>& quadFacePos = quadPoint.position();

        // position of the quadrature point within the subelement
        const FieldVector<double,dim> elementQuadPosSubCell =
                faceComputations.faceToElementPosition(quadFacePos);

        // position of the quadrature point within the reference element
        const FieldVector<double,dim> elementQuadPos =
                subGeometryInReferenceElement.global(elementQuadPosSubCell);

        // The multiplicative factor in the integral transformation formula
        double integrationWeight;
        if(type == IntegrationType::normalVector ||
           type == IntegrationType::travelDistanceWeighted) {
          integrationWeight = localCoefficients.localFactor()(elementQuadPos)
                            * quadPoint.weight()
                            * integrationElement;
          // TODO: scale direction to length 1
          if(type == IntegrationType::travelDistanceWeighted)
            integrationWeight *= std::fabs(direction * unitOuterNormal);
          else
            integrationWeight *= direction * unitOuterNormal;
        } else if(type == IntegrationType::normalSign) {
          const int sign = faceComputations.unitOuterNormalSign();

          integrationWeight = sign
                            * localCoefficients.localFactor()(elementQuadPos)
                            * quadPoint.weight() * integrationElement;
        }

        if(type == IntegrationType::travelDistanceWeighted) {
          integrationWeight *= detail::travelDistance(
              elementQuadPosSubCell,
              referenceBeta);
        }

        //////////////////////////////
        // Left Hand Side Functions //
        //////////////////////////////
        const std::vector<FieldVector<double,1> > lhsValues =
          detail::LocalRefinedFunctionEvaluationHelper
            <is_ContinuouslyRefinedFiniteElement<LhsSpace>::value>::
              evaluateValue(lhsLocalFiniteElement, subElementIndex,
                            elementQuadPosSubCell);

        ///////////////////////////////
        // Right Hand Side Functions //
        ///////////////////////////////
        std::vector<FieldVector<double,1> > rhsValues;
        rhsLocalFiniteElement.localBasis()
            .evaluateFunction(elementQuadPos, rhsValues);

        // Compute the actual matrix entries
        for (size_t i=0; i<nLhs; i++)
        {
          for (size_t j=0; j<nRhs; j++)
          {
            elementMatrix[i+lhsSpaceOffset+subElementOffset]
                         [j+rhsSpaceOffset]
                    += (lhsValues[i] * rhsValues[j]) * integrationWeight;
          }
        }
      }
    }
    if(is_DGRefinedFiniteElement<LhsSpace>::value)
      subElementOffset += lhsLocalFiniteElement.size();
    subElementIndex++;
  }
}
};

}} // end namespace Dune::detail
