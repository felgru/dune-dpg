#include "refinedfaces.hh"

namespace Dune {
namespace detail {

template <IntegrationType type,
          class LhsSpace,
          class RhsSpace>
struct GetLocalMatrix<type, LhsSpace, RhsSpace, true, true>
{
using LhsLocalView = typename LhsSpace::LocalView;
using RhsLocalView = typename RhsSpace::LocalView;

template <class MatrixType,
          class Element,
          class LocalCoefficients>
inline static void interiorImpl(LhsLocalView& lhsLocalView,
                                RhsLocalView& rhsLocalView,
                                MatrixType& elementMatrix,
                                size_t lhsSpaceOffset,
                                size_t rhsSpaceOffset,
                                unsigned int quadratureOrder,
                                const Element& element,
                                const LocalCoefficients& localCoefficients)
{
  constexpr int dim = Element::mydimension;
  const auto geometry = element.geometry();

  typename detail::ChooseQuadrature<LhsSpace, RhsSpace, Element>::type quad
    = detail::ChooseQuadrature<LhsSpace, RhsSpace, Element>
      ::Quadrature(element, quadratureOrder);

  const auto referenceGridView =
      lhsLocalView.tree().refinedReferenceElementGridView();

  unsigned int lhsSubElementOffset = 0;
  unsigned int rhsSubElementOffset = 0;
  unsigned int subElementIndex = 0;
  lhsLocalView.resetSubElements();
  rhsLocalView.resetSubElements();
  for(const auto& subElement : elements(referenceGridView)) {
    lhsLocalView.bindSubElement(subElement);
    // When the IntegralTerm belongs to an InnerProduct, it can happen
    // that lhsLocalView and rhsLocalView are references to the same
    // localView of a test search space. In this case, we have to make
    // sure that we bind only once to the subElement, so that we don't
    // increment the subElement offset twice.
    if(std::addressof(lhsLocalView) != std::addressof(rhsLocalView)) {
      rhsLocalView.bindSubElement(subElement);
    }

    // Get set of shape functions for this subElement
    const auto& lhsLocalFiniteElement = lhsLocalView.tree().finiteElement();
    const auto& rhsLocalFiniteElement = rhsLocalView.tree().finiteElement();

    const unsigned int nLhs(lhsLocalFiniteElement.localBasis().size());
    const unsigned int nRhs(rhsLocalFiniteElement.localBasis().size());

    const auto subGeometryInReferenceElement = subElement.geometry();
    for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();
     // Global position of the current quadrature point
     const FieldVector<double,dim> elementQuadPos
          = subGeometryInReferenceElement.global(quadPos);

      // The multiplicative factor in the integral transformation formula
      const double integrationWeight
        = geometry.integrationElement(elementQuadPos)
        * subGeometryInReferenceElement.integrationElement(quadPos)
        * quad[pt].weight() * localCoefficients.localFactor()(elementQuadPos);

      ///////////////////////////////////////
      // evaluate finite element functions //
      ///////////////////////////////////////

      using FunctionEvaluator
        = detail::LocalRefinedFunctionEvaluation<dim, type>;

      const std::vector<FieldVector<double,1> > lhsValues =
          FunctionEvaluator::template evaluateLhs
                  <is_ContinuouslyRefinedFiniteElement<LhsSpace>::value>
                        (lhsLocalFiniteElement,
                         subElementIndex,
                         quadPos,
                         geometry,
                         subGeometryInReferenceElement,
                         localCoefficients);

      const std::vector<FieldVector<double,1> > rhsValues =
          FunctionEvaluator::template evaluateRhs
                  <is_ContinuouslyRefinedFiniteElement<RhsSpace>::value>
                        (rhsLocalFiniteElement,
                         subElementIndex,
                         quadPos,
                         geometry,
                         subGeometryInReferenceElement,
                         localCoefficients);

      // Compute the actual matrix entries
      for (unsigned int i=0; i<nLhs; i++)
      {
        for (unsigned int j=0; j<nRhs; j++)
        {
          elementMatrix[i+lhsSpaceOffset+lhsSubElementOffset]
                       [j+rhsSpaceOffset+rhsSubElementOffset]
                  += (lhsValues[i] * rhsValues[j]) * integrationWeight;
        }
      }
    }

    if(is_DGRefinedFiniteElement<LhsSpace>::value)
      lhsSubElementOffset += lhsLocalFiniteElement.size();
    if(is_DGRefinedFiniteElement<RhsSpace>::value)
      rhsSubElementOffset += rhsLocalFiniteElement.size();
    subElementIndex++;
  }
}


template <class MatrixType,
          class Element,
          class LocalCoefficients>
inline static void
faceImpl(LhsLocalView& lhsLocalView,
         RhsLocalView& rhsLocalView,
         MatrixType& elementMatrix,
         size_t lhsSpaceOffset,
         size_t rhsSpaceOffset,
         unsigned int quadratureOrder,
         const Element& element,
         const LocalCoefficients& localCoefficients)
{
  constexpr int dim = Element::mydimension;
  const auto geometry = element.geometry();

  const auto referenceGridView =
      lhsLocalView.tree().refinedReferenceElementGridView();

  const auto direction = localCoefficients.localDirection()({0.5,0.5});

  unsigned int lhsSubElementOffset = 0;
  unsigned int rhsSubElementOffset = 0;
  unsigned int subElementIndex = 0;
  lhsLocalView.resetSubElements();
  rhsLocalView.resetSubElements();
  for(const auto& subElement : elements(referenceGridView))
  {
    lhsLocalView.bindSubElement(subElement);
    rhsLocalView.bindSubElement(subElement);

    // Get set of shape functions for this subElement
    const auto& lhsLocalFiniteElement = lhsLocalView.tree().finiteElement();
    const auto& rhsLocalFiniteElement = rhsLocalView.tree().finiteElement();

    const unsigned int nLhs(lhsLocalFiniteElement.localBasis().size());
    const unsigned int nRhs(rhsLocalFiniteElement.localBasis().size());

    using SubElement = std::decay_t<decltype(subElement)>;
    const auto subGeometryInReferenceElement = subElement.geometry();

    unsigned int nOutflowFaces = 0;
    for (unsigned short f = 0, fMax = subElement.subEntities(1); f < fMax; f++)
    {
      const auto face = subElement.template subEntity<1>(f);
      /* This won't work for curvilinear elements, but they don't seem
       * to be supported by UG anyway. */
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
      auto faceComputations
          = RefinedFaceComputations<SubElement>(face, subElement, element);

      using Face = std::decay_t<decltype(face)>;

      const double integrationElement = faceComputations.integrationElement();

      const FieldVector<double,dim> unitOuterNormal
          = faceComputations.unitOuterNormal();

      if(type == IntegrationType::travelDistanceWeighted &&
         direction * unitOuterNormal >= 0) {
        /* Only integrate over inflow boundaries. */
        continue;
      }

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

      for (size_t pt=0, qsize=quadFace.size(); pt < qsize; pt++) {

        // Position of the current quadrature point in the reference element
        // (face!)
        const FieldVector<double,dim-1>& quadFacePos = quadFace[pt].position();

        // position of the quadrature point within the subelement
        const FieldVector<double,dim> elementQuadPosSubCell =
                faceComputations.faceToElementPosition(quadFacePos);

        const FieldVector<double,dim> elementQuadPos =
                subGeometryInReferenceElement.global(elementQuadPosSubCell);

        // The multiplicative factor in the integral transformation formula
        double integrationWeight;
        if(type == IntegrationType::normalVector ||
           type == IntegrationType::travelDistanceWeighted) {
          integrationWeight = localCoefficients.localFactor()(elementQuadPos)
                            * quadFace[pt].weight()
                            * integrationElement;
          // TODO: scale direction to length 1
          if(type == IntegrationType::travelDistanceWeighted)
            integrationWeight *= fabs(direction*unitOuterNormal);
          else
            integrationWeight *= (direction*unitOuterNormal);
        } else if(type == IntegrationType::normalSign) {
          int sign = 1;
          bool signfound = false;
          for (unsigned int i=0;
             i < unitOuterNormal.size() and signfound == false;
             i++)
          {
            if (unitOuterNormal[i]<(-1e-10))
            {
              sign = -1;
              signfound = true;
            }
            else if (unitOuterNormal[i]>(1e-10))
            {
              sign = 1;
              signfound = true;
            }
          }

          integrationWeight = sign
                            * localCoefficients.localFactor()(elementQuadPos)
                            * quadFace[pt].weight() * integrationElement;
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
        const std::vector<FieldVector<double,1> > rhsValues =
          detail::LocalRefinedFunctionEvaluationHelper
            <is_ContinuouslyRefinedFiniteElement<RhsSpace>::value>::
              evaluateValue(rhsLocalFiniteElement, subElementIndex,
                            elementQuadPosSubCell);

        // Compute the actual matrix entries
        for (size_t i=0; i<nLhs; i++)
        {
          for (size_t j=0; j<nRhs; j++)
          {
            elementMatrix[i+lhsSpaceOffset+lhsSubElementOffset]
                         [j+rhsSpaceOffset+rhsSubElementOffset]
                    += (lhsValues[i] * rhsValues[j]) * integrationWeight;
          }
        }
      }
    }
    if(is_DGRefinedFiniteElement<LhsSpace>::value)
      lhsSubElementOffset += lhsLocalFiniteElement.size();
    if(is_DGRefinedFiniteElement<RhsSpace>::value)
      rhsSubElementOffset += rhsLocalFiniteElement.size();
    subElementIndex++;
  }
}
};

}} // end namespace Dune::detail
