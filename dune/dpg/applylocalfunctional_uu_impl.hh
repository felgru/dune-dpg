#include <numeric>
#include <dune/geometry/quadraturerules/splitquadraturerule.hh>
#include "faces.hh"
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
  const int dim = Element::mydimension;
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

  const unsigned int quadratureOrder
      = solutionLocalFiniteElement.localBasis().order()
        + testLocalFiniteElement.localBasis().order();

  typename detail::ChooseQuadrature<TestSpace, SolutionSpace, Element>::type quad
    = detail::ChooseQuadrature<TestSpace, SolutionSpace, Element>
      ::Quadrature(element, quadratureOrder, nullptr);

  // Loop over all quadrature points
  for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = quad[pt].position();

    // The transformed quadrature weight
    const double integrationWeight
        = quad[pt].weight() * geometry.integrationElement(quadPos);

    // Evaluate all shape function values at this quadrature point
    std::vector<FieldVector<double,1>> testShapeFunctionValues;
    testLocalFiniteElement.localBasis()
        .evaluateFunction(quadPos, testShapeFunctionValues);
    std::vector<FieldVector<double,1>> shapeFunctionValues;
    solutionLocalFiniteElement.localBasis().
        evaluateFunction(quadPos, shapeFunctionValues);

    const double functionalValue =
        std::inner_product(
          localFunctionalVector.begin(), localFunctionalVector.end(),
          shapeFunctionValues.begin(), 0.)
        * integrationWeight;
    for (size_t i=0, i_max=testShapeFunctionValues.size(); i<i_max; i++) {
      elementVector[i+spaceOffset] += functionalValue
                                      * testShapeFunctionValues[i];
    }
  }
}


template <class VectorType,
          class Element,
          class FunctionalVector,
          class FactorType,
          class DirectionType>
inline static void
faceImpl(const TestLocalView& testLocalView,
         const SolutionLocalView& solutionLocalView,
         VectorType& elementVector,
         size_t spaceOffset,
         const SolutionLocalIndexSet& solutionLocalIndexSet,
         const Element& element,
         const FunctionalVector& functionalVector,
         const FactorType& factor,
         const DirectionType& beta)
{
  const int dim = Element::mydimension;

  // Get set of shape functions for this element
  const auto& testLocalFiniteElement = testLocalView.tree().finiteElement();
  const auto& solutionLocalFiniteElement
      = solutionLocalView.tree().finiteElement();

  unsigned int nOutflowFaces = 0;
  for (unsigned short f = 0, fMax = element.subEntities(1); f < fMax; f++)
  {
    auto face = element.template subEntity<1>(f);
    const double prod = beta
        * FaceComputations<Element>(face, element).unitOuterNormal();
    if(prod > 0)
      ++nOutflowFaces;
  }

  const auto geometry = element.geometry();

  FieldVector<double,dim> referenceBeta;
  {
    const auto& jacobianInverse = geometry.jacobianInverseTransposed({0., 0.});
    jacobianInverse.mtv(beta, referenceBeta);
  }

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

  for (unsigned short f = 0, fMax = element.subEntities(1); f < fMax; f++)
  {
    auto face = element.template subEntity<1>(f);
    auto faceComputations = FaceComputations<Element>(face, element);
    if(type == IntegrationType::travelDistanceWeighted &&
       beta * faceComputations.unitOuterNormal() >= 0) {
      /* Only integrate over inflow boundaries. */
      continue;
    }

    const unsigned int quadratureOrder
        = solutionLocalFiniteElement.localBasis().order()
          + testLocalFiniteElement.localBasis().order();

    using Face = std::decay_t<decltype(face)>;
    // TODO: Do we really want to have a transport quadrature rule
    //       on the faces, if one of the FE spaces is a transport space?
    QuadratureRule<double, 1> quadFace
      = detail::ChooseQuadrature<TestSpace, SolutionSpace, Face>
        ::Quadrature(face, quadratureOrder, beta);
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

      // position of the quadrature point within the element
      const FieldVector<double,dim> elementQuadPos =
              faceComputations.faceToElementPosition(quadFacePos);

      const FieldVector<double,dim> globalQuadPos =
              geometry.global(elementQuadPos);

      // The multiplicative factor in the integral transformation
      // formula multiplied with outer normal
      const FieldVector<double,dim> integrationOuterNormal =
              faceComputations.integrationOuterNormal();

      // The multiplicative factor in the integral transformation formula -

      double integrationWeight;
      if(type == IntegrationType::normalVector ||
         type == IntegrationType::travelDistanceWeighted) {
        integrationWeight = detail::evaluateFactor(factor, globalQuadPos)
                          * quadFace[pt].weight();
        if(type == IntegrationType::travelDistanceWeighted)
        {
          // |beta * n|*integrationweight
          integrationWeight *= fabs(beta*integrationOuterNormal);
        }
        else
          integrationWeight *= (beta*integrationOuterNormal);
      } else if(type == IntegrationType::normalSign) {
        const double integrationElement =
            face.geometry().integrationElement(quadFacePos);

        const FieldVector<double,dim>& centerOuterNormal =
               faceComputations.unitOuterNormal();

        int sign = 1;
        bool signfound = false;
        for (unsigned int i=0;
           i<centerOuterNormal.size() and signfound == false;
           i++)
        {
          if (centerOuterNormal[i]<(-1e-10))
          {
            sign = -1;
            signfound = true;
          }
          else if (centerOuterNormal[i]>(1e-10))
          {
            sign = 1;
            signfound = true;
          }
        }

        integrationWeight = sign
                          * detail::evaluateFactor(factor, globalQuadPos)
                          * quadFace[pt].weight() * integrationElement;
      }

      if(type == IntegrationType::travelDistanceWeighted) {
        // factor r_K(s)/|beta|
        integrationWeight *= detail::travelDistance(
            elementQuadPos,
            referenceBeta);
      }

      ////////////////////////////////////
      // Left Hand Side Shape Functions //
      ////////////////////////////////////
      std::vector<FieldVector<double,1> > testValues;
      testLocalFiniteElement.localBasis()
          .evaluateFunction(elementQuadPos, testValues);

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
        elementVector[i+spaceOffset] += functionalValue * testValues[i];
      }
    }
  }
}
};

}} // end namespace Dune::detail
