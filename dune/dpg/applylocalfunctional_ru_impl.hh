#include <numeric>
#include <dune/geometry/quadraturerules/splitquadraturerule.hh>
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

  iterateOverLocalIndices(
    solutionLocalView,
    [&](size_t j, auto gj) {
      localFunctionalVector[j] = functionalVector[gj[0]];
    },
    [&](size_t j) { localFunctionalVector[j] = 0; },
    [&](size_t j, auto gj, double wj) {
      localFunctionalVector[j] += wj * functionalVector[gj[0]];
    }
  );

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
      const std::vector<FieldVector<double,1>> testShapeFunctionValues =
        detail::LocalRefinedFunctionEvaluationHelper
          <is_ContinuouslyRefinedFiniteElement<TestSpace>::value>::
            evaluateValue(testLocalFiniteElement, subElementIndex,
                          quadPos);
      std::vector<FieldVector<double,1>> shapeFunctionValues;
      solutionLocalFiniteElement.localBasis().
          evaluateFunction(quadPosInReferenceElement, shapeFunctionValues);

      const double functionalValue =
          std::inner_product(
            localFunctionalVector.begin(), localFunctionalVector.end(),
            shapeFunctionValues.begin(), 0.)
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

  iterateOverLocalIndices(
    solutionLocalView,
    [&](size_t j, auto gj) {
      localFunctionalVector[j] = functionalVector[gj[0]];
    },
    [&](size_t j) { localFunctionalVector[j] = 0; },
    [&](size_t j, auto gj, double wj) {
      localFunctionalVector[j] += wj * functionalVector[gj[0]];
    }
  );

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
    using SubElement = std::decay_t<decltype(subElement)>;

    testLocalView.bindSubElement(subElement);
    const auto& testLocalFiniteElement = testLocalView.tree().finiteElement();

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

      static_assert(hasRequiredQuadratureOrder
                    <typename LocalCoefficients::LocalFactor>::value,
                    "There is no requiredQuadratureOrder specialization"
                    " for the LocalFactor.");
      const unsigned int quadratureOrder
          = solutionLocalFiniteElement.localBasis().order()
            + testLocalFiniteElement.localBasis().order()
            + requiredQuadratureOrder
                <typename LocalCoefficients::LocalFactor>::value;

      QuadratureRule<double, 1> quadFace
        = detail::ChooseQuadrature<TestSpace, SolutionSpace, Face>
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

        ////////////////////////////////////
        // Left Hand Side Shape Functions //
        ////////////////////////////////////
        const std::vector<FieldVector<double,1> > testValues =
          detail::LocalRefinedFunctionEvaluationHelper
            <is_ContinuouslyRefinedFiniteElement<TestSpace>::value>::
              evaluateValue(testLocalFiniteElement, subElementIndex,
                            elementQuadPosSubCell);

        /////////////////////////////////////
        // Right Hand Side Shape Functions //
        /////////////////////////////////////
        std::vector<FieldVector<double,1> > solutionValues;
        solutionLocalFiniteElement.localBasis()
            .evaluateFunction(elementQuadPos, solutionValues);

        const double functionalValue =
            std::inner_product(
              localFunctionalVector.begin(), localFunctionalVector.end(),
              solutionValues.begin(), 0.)
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
