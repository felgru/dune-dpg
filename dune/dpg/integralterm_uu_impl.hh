#include "faces.hh"

namespace Dune {
namespace detail {

template <IntegrationType type,
          class LhsSpace,
          class RhsSpace>
struct GetLocalMatrix<type, LhsSpace, RhsSpace, false, false>
{
using LhsLocalView = typename LhsSpace::LocalView;
using RhsLocalView = typename RhsSpace::LocalView;

template <class MatrixType,
          class Element,
          class LocalCoefficients>
inline static void interiorImpl(const LhsLocalView& lhsLocalView,
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
  const auto& lhsLocalFiniteElement = lhsLocalView.tree().finiteElement();
  const auto& rhsLocalFiniteElement = rhsLocalView.tree().finiteElement();

  const unsigned int nLhs(lhsLocalFiniteElement.size());
  const unsigned int nRhs(rhsLocalFiniteElement.size());

  typename detail::ChooseQuadrature<LhsSpace, RhsSpace, Element>::type quad
    = detail::ChooseQuadrature<LhsSpace, RhsSpace, Element>
      ::Quadrature(element, quadratureOrder);

  for (const auto& quadPoint : quad) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = quadPoint.position();
    // Global position of the current quadrature point

    // The multiplicative factor in the integral transformation formula
    const double integrationWeight = geometry.integrationElement(quadPos)
                                   * quadPoint.weight()
                                   * localCoefficients.localFactor()(quadPos);

    ///////////////////////////////////////
    // evaluate finite element functions //
    ///////////////////////////////////////

    using FunctionEvaluator = detail::LocalFunctionEvaluation<dim, type>;

    FunctionEvaluator::evaluateLhs
                  (lhsValues,
                   lhsLocalFiniteElement,
                   quadPos,
                   geometry,
                   localCoefficients);

    FunctionEvaluator::evaluateRhs
                  (rhsValues,
                   rhsLocalFiniteElement,
                   quadPos,
                   geometry,
                   localCoefficients);

    // Compute the actual matrix entries
    for (unsigned int i=0; i<nLhs; i++)
    {
      for (unsigned int j=0; j<nRhs; j++)
      {
        elementMatrix[i+lhsSpaceOffset][j+rhsSpaceOffset]
                += (lhsValues[i] * rhsValues[j]) * integrationWeight;
      }
    }
  }
}


template <class MatrixType,
          class Element,
          class LocalCoefficients>
inline static void
faceImpl(const LhsLocalView& lhsLocalView,
         const RhsLocalView& rhsLocalView,
         MatrixType& elementMatrix,
         size_t lhsSpaceOffset,
         size_t rhsSpaceOffset,
         unsigned int quadratureOrder,
         const Element& element,
         const LocalCoefficients& localCoefficients)
{
  constexpr int dim = Element::mydimension;

  // Get set of shape functions for this element
  const auto& lhsLocalFiniteElement = lhsLocalView.tree().finiteElement();
  const auto& rhsLocalFiniteElement = rhsLocalView.tree().finiteElement();

  const unsigned int nLhs(lhsLocalFiniteElement.size());
  const unsigned int nRhs(rhsLocalFiniteElement.size());

  const auto direction = localCoefficients.localDirection()({0.5,0.5});

  const unsigned int nOutflowFaces = outflowFacesOfElement(element, direction);

  const auto geometry = element.geometry();

  std::vector<FieldVector<double,1>> lhsValues;
  lhsValues.reserve(lhsLocalView.maxSize());
  std::vector<FieldVector<double,1>> rhsValues;
  rhsValues.reserve(rhsLocalView.maxSize());

  FaceIntegrationData<type> integrationData(geometry, direction);

  for (const auto& face : subEntities(element, Codim<1>{}))
  {
    const auto faceComputations = FaceComputations<Element>(face, element);
    if(faceComputations.template skipFace<type>(direction)) continue;

    const QuadratureRule<double, 1> quadFace
      = faceComputations.template quadratureRule<type, LhsSpace, RhsSpace>
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

      // Left Hand Side Functions
      lhsLocalFiniteElement.localBasis().evaluateFunction(elementQuadPos,
                                                          lhsValues);

      // Right Hand Side Functions
      rhsLocalFiniteElement.localBasis().evaluateFunction(elementQuadPos,
                                                          rhsValues);

      // Compute the actual matrix entries
      for (size_t i=0; i<nLhs; i++)
      {
        for (size_t j=0; j<nRhs; j++)
        {
          elementMatrix[i+lhsSpaceOffset][j+rhsSpaceOffset]
                  += (lhsValues[i] * rhsValues[j]) * integrationWeight;
        }
      }
    }
  }
}
};

}} // end namespace Dune::detail
