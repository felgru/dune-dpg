#include <numeric>
#include <dune/geometry/quadraturerules/splitquadraturerule.hh>
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
using SolutionLocalIndexSet = typename SolutionSpace::LocalIndexSet;

template <class VectorType,
          class Element,
          class FunctionalVector>
inline static void interiorImpl(
    const TestLocalView& testLocalView,
    const SolutionLocalView& solutionLocalView,
    VectorType& elementVector,
    size_t spaceOffset,
    const SolutionLocalIndexSet& solutionLocalIndexSet,
    const Element& element,
    const FunctionalVector& functionalVector)
{
  constexpr int dim = Element::mydimension;
  const auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& testLocalFiniteElement = testLocalView.tree().finiteElement();
  const auto& solutionLocalFiniteElement = solutionLocalView.tree().finiteElement();

  BlockVector<FieldVector<double,1>>
      localFunctionalVector(solutionLocalView.size());

  iterateOverLocalIndexSet(
    solutionLocalIndexSet,
    [&](size_t j, auto gj) {
      localFunctionalVector[j] = functionalVector[gj[0]];
    },
    [&](size_t j) { localFunctionalVector[j] = 0; },
    [&](size_t j, auto gj, double wj) {
      localFunctionalVector[j] += wj * functionalVector[gj[0]];
    }
  );

  const unsigned int quadratureOrder
      = solutionLocalFiniteElement.localBasis().order()
        + testLocalFiniteElement.localBasis().order();

  typename detail::ChooseQuadrature<TestSpace, SolutionSpace, Element>::type quad
    = detail::ChooseQuadrature<TestSpace, SolutionSpace, Element>
      ::Quadrature(element, quadratureOrder);

  const auto referenceGridView =
      testLocalView.tree().refinedReferenceElementGridView();

  const size_t subElementStride =
      (is_DGRefinedFiniteElement<TestSpace>::value) ?
        testLocalFiniteElement.size() : 0;

  unsigned int subElementOffset = 0;
  unsigned int subElementIndex = 0;
  for(const auto& subElement : elements(referenceGridView)) {
    const auto subGeometryInReferenceElement = subElement.geometry();
    for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();
      const FieldVector<double,dim>& quadPosInReferenceElement =
          subGeometryInReferenceElement.global(quadPos);
      // The multiplicative factor in the integral transformation formula
      const double integrationWeight
        = geometry.integrationElement(quadPosInReferenceElement)
        * subGeometryInReferenceElement.integrationElement(quadPos)
        * quad[pt].weight();

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
      for (size_t i=0, i_max=testLocalFiniteElement.size(); i<i_max; i++) {
        elementVector[i+spaceOffset+subElementOffset]
            += functionalValue * testShapeFunctionValues[i];
      }
    }
    if(is_DGRefinedFiniteElement<TestSpace>::value)
      subElementOffset += subElementStride;
    subElementIndex++;
  }
}


template <class VectorType,
          class Element,
          class FunctionalVector,
          class LocalFactor,
          class DirectionType>
inline static void
faceImpl(const TestLocalView& testLocalView,
         const SolutionLocalView& solutionLocalView,
         VectorType& elementVector,
         size_t testSpaceOffset,
         const SolutionLocalIndexSet& solutionLocalIndexSet,
         const Element& element,
         const FunctionalVector& functionalVector,
         const LocalFactor& localFactor,
         const DirectionType& beta)
{
  constexpr int dim = Element::mydimension;
  const auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& testLocalFiniteElement = testLocalView.tree().finiteElement();
  const auto& solutionLocalFiniteElement
      = solutionLocalView.tree().finiteElement();

  BlockVector<FieldVector<double,1>>
      localFunctionalVector(solutionLocalView.size());

  iterateOverLocalIndexSet(
    solutionLocalIndexSet,
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

  const unsigned int subElementStride =
      (is_DGRefinedFiniteElement<TestSpace>::value) ?
        testLocalFiniteElement.size() : 0;

  unsigned int subElementOffset = 0;
  unsigned int subElementIndex = 0;
  for(const auto& subElement : elements(referenceGridView))
  {
    using SubElement = std::decay_t<decltype(subElement)>;
    const auto subGeometryInReferenceElement = subElement.geometry();

    unsigned int nOutflowFaces = 0;
    for (unsigned short f = 0, fMax = subElement.subEntities(1); f < fMax; f++)
    {
      const auto face = subElement.template subEntity<1>(f);
      const FieldVector<double,dim> unitOuterNormal
          = RefinedFaceComputations<SubElement>(face, subElement, element)
              .unitOuterNormal();

      const double prod = beta * unitOuterNormal;
      if(prod > 0)
        ++nOutflowFaces;
    }

    FieldVector<double,dim> referenceBeta
        = detail::referenceBeta(geometry,
            subGeometryInReferenceElement, beta);

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
         beta * unitOuterNormal >= 0) {
        /* Only integrate over inflow boundaries. */
        continue;
      }

      const unsigned int quadratureOrder
          = solutionLocalFiniteElement.localBasis().order()
            + testLocalFiniteElement.localBasis().order();

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

      for (size_t pt=0, qsize=quadFace.size(); pt < qsize; pt++) {

        // Position of the current quadrature point in the reference element
        // (face!)
        const FieldVector<double,dim-1>& quadFacePos = quadFace[pt].position();

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
          integrationWeight = localFactor(elementQuadPos)
                            * quadFace[pt].weight()
                            * integrationElement;
          // TODO: scale beta to length 1
          if(type == IntegrationType::travelDistanceWeighted)
            integrationWeight *= fabs(beta*unitOuterNormal);
          else
            integrationWeight *= (beta*unitOuterNormal);
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
                            * localFactor(elementQuadPos)
                            * quadFace[pt].weight() * integrationElement;
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
        for (size_t i=0, i_max=testValues.size(); i<i_max; i++)
        {
          elementVector[i+testSpaceOffset+subElementOffset]
                  += functionalValue * testValues[i];
        }
      }
    }
    if(is_DGRefinedFiniteElement<TestSpace>::value)
      subElementOffset += subElementStride;
    subElementIndex++;
  }
}
};

}} // end namespace Dune::detail
