#include <numeric>

#include <dune/geometry/quadraturerules/splitquadraturerule.hh>
#include <dune/grid/common/rangegenerators.hh>

#include "quadratureorder.hh"
#include "refinedfaces.hh"
#include "traveldistancenorm.hh"

namespace Dune {
namespace detail {

template <IntegrationType type,
          class TestSpace,
          class SolutionSpace>
struct ApplyLocalFunctional<type, TestSpace, SolutionSpace, true, false>
{
using TestLocalView = typename TestSpace::LocalView;
using SolutionLocalView = typename SolutionSpace::LocalView;

template <class VectorType,
          class Element,
          class FunctionalVector>
inline static void interiorImpl(
    TestLocalView& testLocalView,
    const SolutionLocalView& solutionLocalView,
    VectorType& elementVector,
    size_t spaceOffset,
    const Element& element,
    const FunctionalVector& functionalVector)
{
  constexpr int dim = Element::mydimension;
  const auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& solutionLocalFiniteElement
      = solutionLocalView.tree().finiteElement();

  BlockVector<FieldVector<double,1>>
      localFunctionalVector(solutionLocalView.size());
  copyToLocalVector(functionalVector, localFunctionalVector,
                    solutionLocalView);

  std::vector<FieldVector<double,1>> testShapeFunctionValues;
  testShapeFunctionValues.reserve(testLocalView.maxSize());
  std::vector<FieldVector<double,1>> shapeFunctionValues;
  shapeFunctionValues.reserve(solutionLocalView.maxSize());

  const auto referenceGridView =
      testLocalView.tree().refinedReferenceElementGridView();

  unsigned int subElementOffset = 0;
  unsigned int subElementIndex = 0;
  testLocalView.resetSubElements();
  for(const auto& subElement : elements(referenceGridView)) {
    testLocalView.bindSubElement(subElement);
    const auto& testLocalFiniteElement = testLocalView.tree().finiteElement();

    const unsigned int quadratureOrder
        = solutionLocalFiniteElement.localBasis().order()
          + testLocalFiniteElement.localBasis().order();

    typename detail::ChooseQuadrature<TestSpace, SolutionSpace, Element>::type quad
      = detail::ChooseQuadrature<TestSpace, SolutionSpace, Element>
        ::Quadrature(element, quadratureOrder);

    const auto subGeometryInReferenceElement = subElement.geometry();
    for (const auto& quadPoint : quad) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quadPoint.position();
      const FieldVector<double,dim>& quadPosInReferenceElement =
          subGeometryInReferenceElement.global(quadPos);
      // The multiplicative factor in the integral transformation formula
      const double integrationWeight
        = geometry.integrationElement(quadPosInReferenceElement)
        * subGeometryInReferenceElement.integrationElement(quadPos)
        * quadPoint.weight();

      // Evaluate all shape function values at this quadrature point
      detail::LocalRefinedFunctionEvaluationHelper
        <is_ContinuouslyRefinedFiniteElement<TestSpace>::value>::
          evaluateValue(testShapeFunctionValues,
                        testLocalFiniteElement, subElementIndex,
                        quadPos);
      solutionLocalFiniteElement.localBasis().
          evaluateFunction(quadPosInReferenceElement, shapeFunctionValues);

      const double functionalValue =
          std::inner_product(
            localFunctionalVector.begin(), localFunctionalVector.end(),
            shapeFunctionValues.cbegin(), 0.)
          * integrationWeight;
      auto entry = elementVector.begin() + spaceOffset + subElementOffset;
      for (const auto& shapeFunctionValue : testShapeFunctionValues) {
        *entry += functionalValue * shapeFunctionValue;
        ++entry;
      }
    }
    if(is_DGRefinedFiniteElement<TestSpace>::value)
      subElementOffset += testLocalFiniteElement.size();
    subElementIndex++;
  }
}


template <class VectorType,
          class Element,
          class FunctionalVector,
          class LocalCoefficients>
inline static void
faceImpl(TestLocalView& testLocalView,
         const SolutionLocalView& solutionLocalView,
         VectorType& elementVector,
         size_t testSpaceOffset,
         const Element& element,
         const FunctionalVector& functionalVector,
         const LocalCoefficients& localCoefficients)
{
  constexpr int dim = Element::mydimension;
  const auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& solutionLocalFiniteElement
      = solutionLocalView.tree().finiteElement();

  BlockVector<FieldVector<double,1>>
      localFunctionalVector(solutionLocalView.size());
  copyToLocalVector(functionalVector, localFunctionalVector,
                    solutionLocalView);

  std::vector<FieldVector<double,1>> testValues;
  testValues.reserve(testLocalView.maxSize());
  std::vector<FieldVector<double,1>> solutionValues;
  solutionValues.reserve(solutionLocalView.maxSize());

  const auto referenceGridView =
      testLocalView.tree().refinedReferenceElementGridView();

  static_assert(requiredQuadratureOrder
                <typename LocalCoefficients::LocalDirection>::value == 0,
                "LocalDirection has to be constant.");
  const auto direction = localCoefficients.localDirection()({0.5,0.5});

  unsigned int subElementOffset = 0;
  unsigned int subElementIndex = 0;
  testLocalView.resetSubElements();
  for(const auto& subElement : elements(referenceGridView))
  {
    testLocalView.bindSubElement(subElement);
    const auto& testLocalFiniteElement = testLocalView.tree().finiteElement();

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

      static_assert(hasRequiredQuadratureOrder
                    <typename LocalCoefficients::LocalFactor>::value,
                    "There is no requiredQuadratureOrder specialization"
                    " for the LocalFactor.");
      const unsigned int quadratureOrder
          = solutionLocalFiniteElement.localBasis().order()
            + testLocalFiniteElement.localBasis().order()
            + requiredQuadratureOrder
                <typename LocalCoefficients::LocalFactor>::value;

      const QuadratureRule<double, 1> quadFace = faceComputations
          .template quadratureRule<type, TestSpace, SolutionSpace>
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

        // Left Hand Side Shape Functions
        detail::LocalRefinedFunctionEvaluationHelper
          <is_ContinuouslyRefinedFiniteElement<TestSpace>::value>::
            evaluateValue(testValues,
                          testLocalFiniteElement, subElementIndex,
                          elementQuadPosSubCell);

        // Right Hand Side Shape Functions
        solutionLocalFiniteElement.localBasis()
            .evaluateFunction(elementQuadPos, solutionValues);

        const double functionalValue =
            std::inner_product(
              localFunctionalVector.begin(), localFunctionalVector.end(),
              solutionValues.cbegin(), 0.)
            * integrationWeight;
        auto entry = elementVector.begin()
                      + testSpaceOffset + subElementOffset;
        for (const auto& testValue : testValues)
        {
          *entry += functionalValue * testValue;
          ++entry;
        }
      }
    }
    if(is_DGRefinedFiniteElement<TestSpace>::value)
      subElementOffset += testLocalFiniteElement.size();
    subElementIndex++;
  }
}
};

}} // end namespace Dune::detail
