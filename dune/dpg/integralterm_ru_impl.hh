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

  std::vector<FieldVector<double,1>> lhsValues;
  lhsValues.reserve(lhsLocalView.maxSize());
  std::vector<FieldVector<double,1>> rhsValues;
  rhsValues.reserve(rhsLocalView.maxSize());

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

      LhsFunctionEvaluator::template evaluateLhs
              <is_ContinuouslyRefinedFiniteElement<LhsSpace>::value>
                    (lhsValues,
                     lhsLocalFiniteElement,
                     subElementIndex,
                     quadPos,
                     geometry,
                     subGeometryInReferenceElement,
                     localCoefficients);

      using RhsFunctionEvaluator = detail::LocalFunctionEvaluation<dim, type>;

      RhsFunctionEvaluator::evaluateRhs
                    (rhsValues,
                     rhsLocalFiniteElement,
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

  std::vector<FieldVector<double,1>> lhsValues;
  lhsValues.reserve(lhsLocalView.maxSize());
  std::vector<FieldVector<double,1>> rhsValues;
  rhsValues.reserve(rhsLocalView.maxSize());

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
    lhsLocalView.bindSubElement(subElement);
    const auto& lhsLocalFiniteElement = lhsLocalView.tree().finiteElement();
    const unsigned int nLhs(lhsLocalFiniteElement.localBasis().size());

    const unsigned int nOutflowFaces
        = outflowFacesOfSubElement(subElement, element, direction);

    const auto subGeometryInReferenceElement = subElement.geometry();
    RefinedFaceIntegrationData<type> integrationData(
        geometry, subGeometryInReferenceElement, direction);

    for (const auto& face : subEntities(subElement, Codim<1>{}))
    {
      using SubElement = std::decay_t<decltype(subElement)>;
      const auto faceComputations
          = RefinedFaceComputations<SubElement>(face, subElement, element);
      if(faceComputations.template skipFace<type>(direction)) continue;

      const QuadratureRule<double, 1> quadFace
        = faceComputations.template quadratureRule<type, LhsSpace, RhsSpace>
              (face, quadratureOrder, nOutflowFaces, integrationData);
      for (const auto& quadPoint : quadFace) {
        // position of the quadrature point within the subelement
        const FieldVector<double,dim> elementQuadPosSubCell =
                faceComputations.faceToElementPosition(quadPoint.position());

        // position of the quadrature point within the reference element
        const FieldVector<double,dim> elementQuadPos =
                subGeometryInReferenceElement.global(elementQuadPosSubCell);

        // The multiplicative factor in the integral transformation formula
        const double integrationWeight
            = faceComputations.template integrationWeight<type>(
                  localCoefficients, elementQuadPos, elementQuadPosSubCell,
                  direction, quadPoint.weight(), integrationData);

        // Left Hand Side Functions
        detail::LocalRefinedFunctionEvaluationHelper
          <is_ContinuouslyRefinedFiniteElement<LhsSpace>::value>::
            evaluateValue(lhsValues, lhsLocalFiniteElement, subElementIndex,
                          elementQuadPosSubCell);

        // Right Hand Side Functions
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
