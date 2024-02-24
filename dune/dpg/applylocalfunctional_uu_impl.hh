#include <numeric>

#include <dune/geometry/quadraturerules/splitquadraturerule.hh>
#include <dune/grid/common/rangegenerators.hh>

#include "faces.hh"
#include "quadratureorder.hh"
#include "traveldistancenorm.hh"

namespace Dune {
namespace detail {

template <IntegrationType type,
          class TestSpace,
          class SolutionSpace>
struct ApplyLocalFunctional<type, TestSpace, SolutionSpace, false, false>
{
using TestLocalView = typename TestSpace::LocalView;
using SolutionLocalView = typename SolutionSpace::LocalView;

template <class VectorType,
          class Element,
          class FunctionalVector>
inline static void interiorImpl(
    const TestLocalView& testLocalView,
    const SolutionLocalView& solutionLocalView,
    VectorType& elementVector,
    size_t spaceOffset,
    const Element& element,
    const FunctionalVector& functionalVector)
{
  constexpr int dim = Element::mydimension;
  const auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& testLocalFiniteElement = testLocalView.tree().finiteElement();
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

  const unsigned int quadratureOrder
      = solutionLocalFiniteElement.localBasis().order()
        + testLocalFiniteElement.localBasis().order();

  typename detail::ChooseQuadrature<TestSpace, SolutionSpace, Element>::type quad
    = detail::ChooseQuadrature<TestSpace, SolutionSpace, Element>
      ::Quadrature(element, quadratureOrder);

  for (const auto& quadPoint : quad) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = quadPoint.position();

    // The transformed quadrature weight
    const double integrationWeight
        = quadPoint.weight() * geometry.integrationElement(quadPos);

    testLocalFiniteElement.localBasis()
        .evaluateFunction(quadPos, testShapeFunctionValues);
    solutionLocalFiniteElement.localBasis().
        evaluateFunction(quadPos, shapeFunctionValues);

    const double functionalValue =
        std::inner_product(
          localFunctionalVector.begin(), localFunctionalVector.end(),
          shapeFunctionValues.cbegin(), 0.)
        * integrationWeight;

    auto entry = elementVector.begin() + spaceOffset;
    for (const auto& shapeFunctionValue : testShapeFunctionValues) {
      *entry += functionalValue * shapeFunctionValue;
      ++entry;
    }
  }
}


template <class VectorType,
          class Element,
          class FunctionalVector,
          class LocalCoefficients>
inline static void
faceImpl(const TestLocalView& testLocalView,
         const SolutionLocalView& solutionLocalView,
         VectorType& elementVector,
         size_t spaceOffset,
         const Element& element,
         const FunctionalVector& functionalVector,
         const LocalCoefficients& localCoefficients)
{
  constexpr int dim = Element::mydimension;

  // Get set of shape functions for this element
  const auto& testLocalFiniteElement = testLocalView.tree().finiteElement();
  const auto& solutionLocalFiniteElement
      = solutionLocalView.tree().finiteElement();

  static_assert(requiredQuadratureOrder
                <typename LocalCoefficients::LocalDirection>::value == 0,
                "LocalDirection has to be constant.");
  const auto direction = localCoefficients.localDirection()({0.5,0.5});

  const unsigned int nOutflowFaces = outflowFacesOfElement(element, direction);

  const auto geometry = element.geometry();

  FaceIntegrationData<type> integrationData(geometry, direction);

  BlockVector<FieldVector<double,1>>
      localFunctionalVector(solutionLocalView.size());
  copyToLocalVector(functionalVector, localFunctionalVector,
                    solutionLocalView);

  // Left Hand Side Shape Functions
  std::vector<FieldVector<double,1> > testValues;
  testValues.resize(testLocalView.maxSize());
  // Right Hand Side Shape Functions
  std::vector<FieldVector<double,1> > solutionValues;
  solutionValues.resize(solutionLocalView.maxSize());

  for (const auto& face : subEntities(element, Codim<1>{}))
  {
    const auto faceComputations = FaceComputations<Element>(face, element);
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
      // position of the quadrature point within the element
      const FieldVector<double,dim> elementQuadPos =
              faceComputations.faceToElementPosition(quadPoint.position());

      // The multiplicative factor in the integral transformation formula
      const double integrationWeight
          = faceComputations.template integrationWeight<type>(
              localCoefficients, elementQuadPos, direction, quadPoint,
              integrationData, face);

      testLocalFiniteElement.localBasis()
          .evaluateFunction(elementQuadPos, testValues);

      solutionLocalFiniteElement.localBasis()
          .evaluateFunction(elementQuadPos, solutionValues);

      const double functionalValue =
          std::inner_product(
            localFunctionalVector.begin(), localFunctionalVector.end(),
            solutionValues.begin(), 0.)
          * integrationWeight;

      auto entry = elementVector.begin() + spaceOffset;
      for(const auto& testValue : testValues) {
        *entry += functionalValue * testValue;
        ++entry;
      }
    }
  }
}
};

}} // end namespace Dune::detail
